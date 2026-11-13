greedy_clique <- function(adj) {
    adj <- adj | t(adj)
    diag(adj) <- FALSE
    deg <- rowSums(adj)
    best <- if (nrow(adj)) 1L else 0L
    for (v in seq_len(nrow(adj))) {
        members <- v
        nb <- which(adj[v, ])
        for (u in nb[order(-deg[nb])]) if (all(adj[u, members])) members <- c(members, u)
        best <- max(best, length(members))
    }
    as.integer(best)
}

first_minimum <- function(ebic, usable) {
    idx <- which(usable & is.finite(ebic))
    if (!length(idx)) return(NA_integer_)
    idx[which.min(ebic[idx])]
}

refit_within <- function(path, Y, k, gamma, seconds, tol = 0) {
    one <- path
    one$eta <- path$eta[k]
    one$lambda <- path$lambda[k]
    one$fits <- path$fits[k]
    job <- parallel::mcparallel(kgm_score(one, Y, "refit", gamma = gamma, tol = tol),
                               silent = TRUE)
    res <- parallel::mccollect(job, wait = FALSE, timeout = seconds)
    if (is.null(res)) {
        tools::pskill(job$pid, tools::SIGKILL)
        suppressWarnings(parallel::mccollect(job, wait = TRUE))
        return(data.frame(status = "timeout", loglik = NA_real_, ebic = NA_real_,
                          message = "", stringsAsFactors = FALSE))
    }
    value <- res[[1]]
    if (inherits(value, "try-error")) {
        return(data.frame(status = "error", loglik = NA_real_, ebic = NA_real_,
                          message = conditionMessage(attr(value, "condition")),
                          stringsAsFactors = FALSE))
    }
    data.frame(status = if (value$status == "refit_failed") "error" else value$status,
               loglik = value$loglik, ebic = value$ebic, message = value$message,
               stringsAsFactors = FALSE)
}

pilot_refits <- function(path, Y, gamma, seconds, tol = 0) {
    dims <- dim(Y)
    M <- length(dims) - 1L
    p <- dims[seq_len(M)]
    rank_cap <- dims[M + 1L] * prod(p) / p
    memo <- list()
    rows <- lapply(seq_along(path$fits), function(k) {
        mask <- lapply(path$fits[[k]], function(A) {
            keep <- abs(A) > tol
            diag(keep) <- FALSE
            keep
        })
        key <- digest::digest(mask, algo = "sha256")
        if (is.null(memo[[key]])) {
            cliques <- vapply(mask, greedy_clique, 0L)
            res <- if (any(cliques > rank_cap)) {
                data.frame(status = "nonexistent", loglik = NA_real_, ebic = NA_real_,
                           message = "", stringsAsFactors = FALSE)
            } else {
                refit_within(path, Y, k, gamma, seconds, tol = tol)
            }
            res$clique <- paste(cliques, collapse = "/")
            memo[[key]] <<- res
        }
        memo[[key]]
    })
    do.call(rbind, rows)
}

pilot_config <- function(seconds, gamma, eta, driver_code) {
    driver_config(list(name = "refit_pilot", seconds = seconds, gamma = gamma), eta,
                  driver_code)
}

selected_or_none <- function(label, k, path, truth) {
    if (!is.na(k)) return(fit_tables(label, k, path$eta[k], path$fits[[k]], truth))
    none <- fit_tables(label, NA_integer_, NA_real_, path$fits[[1L]], truth)
    none$selected[setdiff(names(none$selected), names(label))] <- NA
    none$selected_modes[setdiff(names(none$selected_modes), c(names(label), "mode"))] <- NA
    none
}

run_pilot_dataset <- function(dataset, row, config, out, path_fun) {
    file <- file.path(out, "datasets", paste0(dataset$id, ".rds"))
    signature <- dataset_signature(row, dataset$rep, config)
    if (!is.null(reusable_result(file, signature))) return("reused")
    gen <- make_dataset(row, dataset$rep)
    gamma <- config$method$gamma
    t0 <- proc.time()[["elapsed"]]
    path <- path_fun(dataset, gen$Y, config$eta)
    t1 <- proc.time()[["elapsed"]]
    rx <- kgm_score(path, gen$Y, "relaxation", gamma = gamma, gradient = row$gradient,
                    tol = 0)
    t2 <- proc.time()[["elapsed"]]
    rf <- pilot_refits(path, gen$Y, gamma, config$method$seconds)
    t3 <- proc.time()[["elapsed"]]
    candidates <- data.frame(k = rx$k, eta = rx$eta, edges = rx$edges,
                             relaxation_status = rx$status,
                             relaxation_loglik = rx$loglik, relaxation_ebic = rx$ebic,
                             refit_status = rf$status, refit_loglik = rf$loglik,
                             refit_ebic = rf$ebic, clique = rf$clique,
                             stringsAsFactors = FALSE)
    chosen <- c(relaxation = first_minimum(rx$ebic, rx$status == "ok"),
                refit = first_minimum(rf$ebic, rf$status == "ok"))
    tables <- lapply(names(chosen), function(m) {
        selected_or_none(data.frame(criterion = m), chosen[[m]], path, gen$truth)
    })
    save_atomic(list(id = dataset$id, setting = dataset$setting, rep = dataset$rep,
                     seeds = gen$seeds, signature = signature,
                     data_sha256 = digest::digest(gen$Y, algo = "sha256"),
                     truth_sha256 = digest::digest(gen$truth, algo = "sha256"),
                     candidates = candidates,
                     selected = do.call(rbind, lapply(tables, `[[`, "selected")),
                     selected_modes = do.call(rbind, lapply(tables, `[[`, "selected_modes")),
                     seconds = c(path = t1 - t0, relaxation = t2 - t1, refit = t3 - t2)),
                file)
    "done"
}

run_pilot_experiment <- function(experiment, settings_file, out, driver_code,
                                 workers = 1L, reps = NULL, ids = NULL,
                                 path_fun = fit_path, eta = kgm_eta_grid(),
                                 seconds = 30L, gamma = 0.5) {
    config <- pilot_config(seconds, gamma, eta, driver_code)
    run_datasets(experiment, settings_file, out, config,
                 function(dataset, row) run_pilot_dataset(dataset, row, config, out, path_fun),
                 workers = workers, reps = reps, ids = ids)
}
