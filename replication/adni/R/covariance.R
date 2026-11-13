check_covariance <- function(opts, root) {
    BLAS_ACTUAL <- load_adni_dependencies()
    progress_every <- 6L
    code_paths <- list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)
    names(code_paths) <- basename(code_paths)
    option_arg <- opts$option
    if (!option_arg %in% c("C", "A")) {
        stop(sprintf("--option 은 C 또는 A 다 (받은 것: %s)", option_arg))
    }
    groups <- split_csv(opts$groups)
    arms <- split_csv(opts$arms)
    if (!all(arms %in% names(ARM_SPECS))) {
        stop(sprintf("모르는 팔: %s. 있는 것은 %s",
                     paste(setdiff(arms, names(ARM_SPECS)), collapse = ", "),
                     paste(names(ARM_SPECS), collapse = ", ")))
    }
    n_eta <- suppressWarnings(as.integer(opts$`n-eta`))
    if (is.na(n_eta) || n_eta < 2L) {
        stop(sprintf("--n-eta 는 2 이상의 정수다 (받은 것: %s)", opts$`n-eta`))
    }
    eta_grid <- make_eta_grid(n_eta)

    out_csv <- opts$`out-csv`

    dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

    aggregated_dir <- if (nzchar(opts$aggregated)) opts$aggregated else file.path(root, "data/processed/antspy/aggregated", option_arg)
    np <- reticulate::import("numpy")

    t_start <- Sys.time()
    started <- format(t_start, "%Y-%m-%dT%H:%M:%S%z")

    spectrum_rows <- function(x_val, m_vec, n_here, eta, i, group, arm) {
        base <- data.frame(
            group = group, arm = arm, eta_index = i, eta = eta,
            stringsAsFactors = FALSE
        )
        lambda <- lambda_from_eta(m_vec, n_here, eta)
        fit <- tryCatch(
            fit_tggm(x_val, lambda, n_sweep = N_SWEEP, return_S = TRUE),
            error = function(e) e
        )
        if (inherits(fit, "error")) {
            return(list(cbind(base, data.frame(
                sweep = NA_integer_, mode = NA_integer_, p_m = NA_integer_,
                eigen_min = NA_real_, eigen_max = NA_real_, asym = NA_real_,
                error = conditionMessage(fit), stringsAsFactors = FALSE
            ))))
        }
        S_by_sweep <- attr(fit, "S_by_sweep")
        stopifnot(length(S_by_sweep) == N_SWEEP)

        rows <- list()
        for (sw in seq_along(S_by_sweep)) {
            for (m in seq_along(S_by_sweep[[sw]])) {
                Sm <- S_by_sweep[[sw]][[m]]
                asym <- max(abs(Sm - t(Sm)))

                ev <- eigen((Sm + t(Sm)) / 2, symmetric = TRUE,
                            only.values = TRUE)$values
                rows[[length(rows) + 1L]] <- cbind(base, data.frame(
                    sweep = sw, mode = m, p_m = nrow(Sm),
                    eigen_min = min(ev), eigen_max = max(ev), asym = asym,
                    error = NA_character_, stringsAsFactors = FALSE
                ))
            }
        }
        rows
    }

    rows <- list()
    for (group in groups) {
        for (arm in arms) {
            d <- load_arm(group, arm, aggregated_dir, np)
            x_val <- aperm(residualize(d$Y, d$X), c(2, 3, 1))
            m_vec <- dim(x_val)[1:2]
            n_here <- dim(x_val)[3]
            cat(sprintf("=== %s / %s  p=(%d,%d) n=%d ===\n",
                        group, arm, m_vec[1], m_vec[2], n_here))
            for (i in seq_along(eta_grid)) {
                rows <- c(rows, spectrum_rows(x_val, m_vec, n_here,
                                              eta_grid[i], i, group, arm))
                if (i %% progress_every == 0L) {
                    cat(sprintf("  [%2d/%2d]\n", i, length(eta_grid)))
                    utils::flush.console()
                }
            }
        }
    }
    res <- do.call(rbind, rows)
    utils::write.csv(res, out_csv, row.names = FALSE)

    run_json <- file.path(dirname(out_csv), "RUN_assumption_check.json")
    run_meta <- list(
        stage = "adni-assumption-check",
        started = started,
        elapsed_sec = round(
            as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1
        ),
        blas_threads = BLAS_ACTUAL,
        option = option_arg,
        groups = groups,
        arms = arms,
        arm_p_time = vapply(arms, function(a) ARM_SPECS[[a]]$p_time, integer(1)),
        eta_grid = eta_grid,
        n_sweep = N_SWEEP,
        sqrt_method = SQRT_METHOD,
        n_failed_points = sum(!is.na(res$error)),
        code_sha256 = code_fingerprints(code_paths),
        r_version = R.version.string,
        pkg = list(
            huge = as.character(utils::packageVersion("huge")),
            rTensor = as.character(utils::packageVersion("rTensor")),
            expm = as.character(utils::packageVersion("expm")),
            reticulate = as.character(utils::packageVersion("reticulate")),
            jsonlite = as.character(utils::packageVersion("jsonlite")),
            RhpcBLASctl = as.character(utils::packageVersion("RhpcBLASctl"))
        )
    )
    writeLines(
        jsonlite::toJSON(run_meta, auto_unbox = TRUE, pretty = TRUE, digits = NA),
        run_json
    )

    cat("\n================ 요약 (경로 전체에 대한 최소값) ================\n")
    ok <- res[is.na(res$error), ]
    if (nrow(ok) == 0L) {
        cat("성공한 격자점이 없다.\n")
    } else {
        agg <- aggregate(
            cbind(eigen_min, ratio = eigen_min / eigen_max) ~ group + arm + mode,
            data = ok, FUN = min
        )
        print(agg[order(agg$arm, agg$group, agg$mode), ], row.names = FALSE)
    }
    n_failed <- sum(!is.na(res$error))
    if (n_failed > 0L) cat(sprintf("\n실패한 격자점 %d 개 (CSV 의 error 열)\n",
                                   n_failed))
    cat(sprintf("\nCSV: %s\nRUN: %s\n", out_csv, run_json))
}
