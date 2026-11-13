cv_fold_id <- function(n, nfolds = 5L, seed = 2026L) {
    kind <- RNGkind()
    saved <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        get(".Random.seed", envir = globalenv())
    }
    on.exit({
        RNGkind(kind[1], kind[2], kind[3])
        if (is.null(saved)) {
            rm(".Random.seed", envir = globalenv())
        } else {
            assign(".Random.seed", saved, envir = globalenv())
        }
    })
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
             sample.kind = "Rejection")
    sample(rep(seq_len(nfolds), length.out = n))
}

cv_select_index <- function(fold_loglik) unname(which.max(colMeans(fold_loglik)))

cv_config <- function(eta, driver_code, nfolds = 5L, fold_seed = 2026L) {
    driver_config(list(name = "cv", nfolds = nfolds, fold_seed = fold_seed), eta,
                  driver_code)
}

fit_cv <- function(dataset, Y, eta, fold_id) kgm_cv(Y, eta, fold_id = fold_id)

run_cv_dataset <- function(dataset, row, config, out, cv_fun) {
    file <- file.path(out, "datasets", paste0(dataset$id, ".rds"))
    signature <- dataset_signature(row, dataset$rep, config)
    if (!is.null(reusable_result(file, signature))) return("reused")
    gen <- make_dataset(row, dataset$rep)
    fold_id <- cv_fold_id(row$n, config$method$nfolds, config$method$fold_seed)
    t0 <- proc.time()[["elapsed"]]
    cv <- cv_fun(dataset, gen$Y, config$eta, fold_id)
    t1 <- proc.time()[["elapsed"]]
    if (!identical(cv$fold_id, fold_id)) stop("the cross-validation used other folds")
    k <- cv_select_index(cv$fold_loglik)
    if (!identical(as.integer(cv$index), k)) {
        stop("the selected grid point does not maximise the mean held-out log-likelihood")
    }
    tables <- fit_tables(data.frame(criterion = "cv"), k, config$eta[k], cv$factors,
                         gen$truth)
    save_atomic(list(id = dataset$id, setting = dataset$setting, rep = dataset$rep,
                     seeds = gen$seeds, signature = signature,
                     data_sha256 = digest::digest(gen$Y, algo = "sha256"),
                     truth_sha256 = digest::digest(gen$truth, algo = "sha256"),
                     fold_id = fold_id, fold_loglik = cv$fold_loglik,
                     fold_tau = cv$fold_tau, tau = cv$tau,
                     selected = tables$selected, selected_modes = tables$selected_modes,
                     seconds = c(cv = t1 - t0)),
                file)
    "done"
}

run_cv_experiment <- function(experiment, settings_file, out, driver_code,
                              workers = 1L, reps = NULL, ids = NULL,
                              cv_fun = fit_cv, eta = kgm_eta_grid()) {
    config <- cv_config(eta, driver_code)
    run_datasets(experiment, settings_file, out, config,
                 function(dataset, row) run_cv_dataset(dataset, row, config, out, cv_fun),
                 workers = workers, reps = reps, ids = ids)
}
