run_arm <- function(group, arm, aggregated_dir, np, eta_grid) {
    spec <- ARM_SPECS[[arm]]
    cat(sprintf("\n=== %s / %s (%s) ===\n", group, arm, spec$label))

    d <- load_arm(group, arm, aggregated_dir, np)
    resid <- residualize(d$Y, d$X)

    x_val <- aperm(resid, c(2, 3, 1))
    m_vec <- dim(x_val)[1:2]
    n_here <- dim(x_val)[3]
    cat(sprintf("  p_ROI=%d  p_T=%d  n=%d\n", m_vec[1], m_vec[2], n_here))

    rows <- vector("list", length(eta_grid))
    prev_fit <- NULL
    prev_lambda <- NULL
    for (i in seq_along(eta_grid)) {
        eta <- eta_grid[i]
        lambda <- lambda_from_eta(m_vec, n_here, eta)
        t0 <- Sys.time()
        fit <- tryCatch(
            lapply(fit_tggm(x_val, lambda, n_sweep = N_SWEEP), function(M) M),
            error = function(e) e
        )
        scored <- tryCatch(
            {
                if (inherits(fit, "error")) stop(fit)
                score_fit(
                    x_val, fit, n_here, prev_fit = prev_fit,
                    lambda = lambda, prev_lambda = prev_lambda
                )
            },
            error = function(e) {
                message(sprintf(
                    "  ! eta=%.4g 실패: %s", eta, conditionMessage(e)
                ))
                na_score_row()
            }
        )

        if (all(is.finite(scored$ebic_joint))) {
            prev_fit <- fit
            prev_lambda <- lambda
        } else {
            prev_fit <- NULL
            prev_lambda <- NULL
        }
        elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

        rows[[i]] <- data.frame(
            group = group, arm = arm, p_roi = m_vec[1], p_time = m_vec[2],
            n = n_here, eta_index = i, eta = eta,
            lambda_roi = lambda[1], lambda_time = lambda[2],
            scored, elapsed_sec = elapsed,
            stringsAsFactors = FALSE
        )
        cat(sprintf(
            "  [%2d/%2d] eta=%-10.4g edges=(%5s,%6s) tau=%9.3g u=(%.3g,%.3g)  %.1fs\n",
            i, length(eta_grid), eta, scored$edges_roi[1],
            scored$edges_time[1], scored$tau_global[1], scored$u_roi[1],
            scored$u_time[1], elapsed
        ))
        utils::flush.console()
    }
    do.call(rbind, rows)
}

run_analysis <- function(opts, root) {
    BLAS_ACTUAL <- load_adni_dependencies()
    code_paths <- list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)
    names(code_paths) <- basename(code_paths)
    option_arg <- opts$option
    groups <- split_csv(opts$groups)
    arms <- split_csv(opts$arms)
    n_eta <- as.integer(opts$`n-eta`)
    out_dir <- opts$out

    stopifnot(
        length(groups) >= 1L, length(arms) >= 1L,
        all(arms %in% names(ARM_SPECS)),
        option_arg %in% c("C", "A"),
        !is.na(n_eta)
    )
    eta_grid <- make_eta_grid(n_eta)

    aggregated_dir <- if (nzchar(opts$aggregated)) opts$aggregated else file.path(root, "data/processed/antspy/aggregated", option_arg)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    np <- reticulate::import("numpy")

    started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z", tz = "UTC")
    t_start <- Sys.time()

    path <- do.call(rbind, lapply(groups, function(g) {
        do.call(rbind, lapply(arms, function(a) {
            run_arm(g, a, aggregated_dir, np, eta_grid)
        }))
    }))

    path_csv <- file.path(out_dir, "probe_path.csv")
    verdict_csv <- file.path(out_dir, "probe_verdict.csv")
    run_json <- file.path(out_dir, "RUN_probe.json")

    utils::write.csv(path, path_csv, row.names = FALSE)
    verdict <- verdict_table(path)
    utils::write.csv(verdict, verdict_csv, row.names = FALSE)

    run_meta <- list(
        stage = "adni-resample-probe",
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
        gamma_grid = GAMMA_GRID,
        n_sweep = N_SWEEP,
        support_tol = SUPPORT_TOL,
        sqrt_method = SQRT_METHOD,
        joint_u_cap = JOINT_U_CAP,
        joint_boundary_margin = JOINT_BOUNDARY_MARGIN,
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

    cat("\n================ 판정 ================\n")
    print(verdict[verdict$gamma == GAMMA_REPORT, ], row.names = FALSE)
    cat(sprintf(
        "\n경로: %s\n판정: %s\n명세: %s\n", path_csv, verdict_csv, run_json
    ))
}
