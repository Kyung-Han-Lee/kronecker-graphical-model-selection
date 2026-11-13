fit_path_with_inputs <- function(dataset, Y, eta) {
    dims <- dim(Y)
    M <- length(dims) - 1L
    p <- dims[seq_len(M)]
    n <- dims[M + 1L]
    lambda <- lapply(eta, function(value) kgm_lambda(p, n, value))
    covariances <- list()
    fits <- lapply(seq_along(eta), function(k) {
        fit <- kgmselect:::.fit_tggm(Y, lambda[[k]], return_S = TRUE)
        values <- attr(fit, "S_by_sweep")
        for (sweep in seq_along(values)) {
            for (mode in seq_len(M)) {
                covariances[[length(covariances) + 1L]] <<- data.frame(
                    k = k, sweep = sweep, mode = mode,
                    minimum_eigenvalue = min(eigen(values[[sweep]][[mode]],
                        symmetric = TRUE, only.values = TRUE)$values))
            }
        }
        lapply(fit, identity)
    })
    path <- kgmselect:::.new_kgm_path(eta, lambda, fits, 2L, p, n)
    attr(path, "partial_covariances") <- do.call(rbind, covariances)
    path
}

cached_observations <- function(result) {
    if (is.null(result$inputs$path) || is.null(result$inputs$settings)) {
        stop("candidate inputs were not retained; run with --save_inputs=1")
    }
    generated <- make_dataset(result$inputs$settings, result$rep)
    if (!identical(digest::digest(generated$Y, algo = "sha256"), result$data_sha256)) {
        stop("regenerated observations do not match the cached dataset")
    }
    generated
}

certificate_threshold <- function(data) {
    dims <- dim(data)
    stopifnot(length(dims) == 3L)
    p_all <- prod(dims[1:2])
    n_obs <- dims[3]
    sq <- rowSums(matrix(data^2, nrow = p_all)) / n_obs
    if (any(!is.finite(sq)) || any(sq <= 0)) {
        stop("sample covariance diagonal is nonpositive or nonfinite")
    }
    log(p_all) + mean(log(sq))
}

candidate_metadata <- function(result, experiment, cache) {
    row <- result$inputs$settings
    data.frame(id = paste(experiment, result$id, sep = "__"),
        experiment = experiment, setting = result$setting, M = row$M,
        dimensions = row$dims, family = toupper(row$family),
        degree = if (row$family == "er") row$d else row$m_new,
        n = row$n, rep = result$rep, grid_points = length(result$inputs$path$eta),
        cache = cache, stringsAsFactors = FALSE)
}

candidate_profile <- function(Y, path, k, row, recovery, id) {
    fit <- path$fits[[k]]
    delta <- kgmselect:::.relaxation_deltas(path, k, tol = 0)
    info <- Map(function(A, D) kgmselect:::.relaxation_factor_info(A, D), fit, delta)
    profile <- kgmselect:::.profile_relaxation(Y, fit, delta, gradient = row$gradient)
    coefficients <- kgmselect:::.joint_trace_coefficients(Y, info)
    M <- length(fit)
    p <- vapply(fit, nrow, 0L)
    n <- dim(Y)[M + 1L]
    q <- sum(vapply(info, function(value) isTRUE(value$active), TRUE))
    mu <- lapply(info, function(value) value$mu)
    upper <- vapply(info, function(value) value$u_upper, 0)
    objective <- function(u) {
        determinant <- sum(vapply(seq_len(M), function(m) {
            sum(log1p(u[m] * mu[[m]])) / p[m]
        }, 0))
        terms <- vapply(0:(2^M - 1L), function(mask) {
            active <- as.logical(intToBits(as.integer(mask)))[seq_len(M)]
            prod(u[active])
        }, 0)
        determinant - log(sum(terms * coefficients) / coefficients[1L])
    }
    replay <- profile$loglik_u0 + n * prod(p) * objective(profile$u) / 2
    list(profile = profile, q = q,
         components = list(id = id, n = n, k = k, p_vec = p,
             mu = mu, coefficients = coefficients, upper = upper, q = q,
             fp = recovery$fp[k], fn = recovery$fn[k],
             u_replay = profile$u, loglik_zero = profile$loglik_u0,
             loglik_saved = profile$loglik,
             replay_difference = replay - profile$loglik,
             spectral_probe_max_error = abs(replay - profile$loglik)))
}

export_candidate_dataset <- function(result, experiment, cache, destination, profiles) {
    generated <- cached_observations(result)
    Y <- generated$Y
    path <- result$inputs$path
    row <- result$inputs$settings
    metadata <- candidate_metadata(result, experiment, cache)
    id <- metadata$id
    recovery <- result$recovery
    if (is.null(recovery)) recovery <- candidate_recovery(path, generated$truth)
    metrics <- lapply(seq_along(path$eta), function(k) {
        cbind(k = k, as.data.frame(kgm_metrics(generated$truth, path$fits[[k]], tol = 0)[
            c("tp", "fp", "fn", "tn", "mcc", "fdp", "tpr")]))
    })
    metrics <- do.call(rbind, metrics)
    scores <- result$scores
    metric_rows <- metrics[match(scores$k, metrics$k), setdiff(names(metrics), "k")]
    scores <- cbind(id = id, scores, metric_rows, stringsAsFactors = FALSE)
    append_exact_csv(scores, file.path(destination, "scores.csv"))
    append_exact_csv(metadata, file.path(destination, "paths.csv"))
    supports <- support_records(path)
    append_exact_csv(cbind(id = id, supports, stringsAsFactors = FALSE),
                     file.path(destination, "supports.csv"))
    covariances <- attr(path, "partial_covariances")
    if (is.null(covariances)) stop("partial covariance diagnostics are missing")
    append_exact_csv(cbind(id = id, covariances, stringsAsFactors = FALSE),
                     file.path(destination, "covariances.csv"))
    B <- if (row$M == 2L) certificate_threshold(Y) else NA_real_
    truth_scale <- prod(vapply(generated$truth, norm, 0, type = "F"))
    save_profile <- identical(profiles, "all") || (identical(profiles, "paper") &&
        experiment == "matrix" && row$dims == "116x102" && row$n %in% c(100L, 3200L) &&
        ((row$family == "er" && row$d == 2) || row$family == "ba"))
    components <- list()
    for (k in seq_along(path$eta)) {
        value <- candidate_profile(Y, path, k, row, recovery, id)
        score <- scores[scores$criterion == "relaxation" & scores$gamma == 0.5 & scores$k == k, ]
        if (nrow(score) != 1L || score$status != "ok" ||
            abs(value$profile$loglik - score$loglik) > 1e-5) {
            stop("profile replay differs from the saved score: ", id, " candidate ", k)
        }
        append_exact_csv(data.frame(id = id, k = k, M = row$M, q = value$q),
                         file.path(destination, "intervals.csv"))
        scale <- scores[scores$criterion == "scale" & scores$k == k, ]
        if (!nrow(scale)) stop("scale scores are missing")
        append_exact_csv(data.frame(id = id, k = k, tau = scale$tau[1L],
                                    truth_scale = truth_scale, B = B),
                         file.path(destination, "scales.csv"))
        append_exact_csv(data.frame(id = id, k = k, mode = seq_len(row$M),
            cond_base = value$profile$cond_base,
            cond_relaxed = value$profile$cond_joint),
            file.path(destination, "conditions.csv"))
        if (save_profile) {
            value$components$support_hash <- supports$pattern_hash[k]
            value$components$ebic_saved <- score$ebic
            components[[k]] <- value$components
        }
    }
    if (!save_profile) return(NULL)
    selected <- scores[scores$criterion == "relaxation" & scores$gamma == 0.5, ]
    k <- select_index(selected)
    exact <- recovery$fp == 0 & recovery$fn == 0
    outcome <- if (exact[k]) "exact_selected" else if (any(exact)) {
        "available_not_selected"
    } else "no_true_candidate"
    profile_dir <- file.path(destination, "profiles")
    dir.create(profile_dir, showWarnings = FALSE)
    jsonlite::write_json(list(candidates = components), file.path(profile_dir, paste0(id, ".json")),
                         auto_unbox = TRUE, pretty = FALSE, digits = NA)
    list(id = id, n = as.integer(row$n), rep = as.integer(result$rep),
         family = toupper(row$family), dimensions = row$dims,
         degree = metadata$degree, selected_k = k,
         exact_candidate_count = sum(exact), outcome = outcome,
         grid_points = length(path$eta))
}

prepare_candidate_inputs <- function(outputs, destination, profiles = "paper") {
    if (!profiles %in% c("paper", "all", "none")) stop("profiles must be paper, all or none")
    if (file.exists(destination) || dir.exists(destination)) stop("output already exists: ", destination)
    folders <- c("matrix", "tensor", "refit_comparison")
    folders <- folders[file.exists(file.path(outputs, folders, "run.json"))]
    if (!length(folders)) stop("no matrix, tensor or Refit-comparison results")
    parent <- dirname(destination)
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    staging <- tempfile(".inputs-", tmpdir = parent)
    dir.create(staging)
    on.exit(unlink(staging, recursive = TRUE))
    cases <- list()
    counts <- list()
    for (experiment in folders) {
        directory <- file.path(outputs, experiment)
        run <- jsonlite::read_json(file.path(directory, "run.json"), simplifyVector = TRUE)
        files <- sort(list.files(file.path(directory, "datasets"), pattern = "[.]rds$",
                                 full.names = TRUE))
        if (length(files) != run$datasets$total || run$datasets$failed != 0) {
            stop("incomplete experiment outputs: ", experiment)
        }
        for (file in files) {
            result <- readRDS(file)
            if (is.null(result$inputs)) stop("run ", experiment, " with --save_inputs=1")
            profile <- export_candidate_dataset(result, experiment,
                file.path(experiment, "datasets", basename(file)), staging, profiles)
            if (!is.null(profile)) cases[[length(cases) + 1L]] <- profile
        }
        counts[[experiment]] <- length(files)
    }
    if (length(cases)) {
        jsonlite::write_json(cases, file.path(staging, "cases.json"),
                             auto_unbox = TRUE, pretty = TRUE, digits = NA)
    }
    jsonlite::write_json(list(schema_version = 1L, experiments = counts,
                              profile_paths = length(cases)),
                        file.path(staging, "inputs.json"), auto_unbox = TRUE,
                        pretty = TRUE, digits = NA)
    if (file.exists(destination) || dir.exists(destination)) stop("output already exists: ", destination)
    if (!file.rename(staging, destination)) stop("could not save prepared inputs")
    invisible(destination)
}
