bootstrap_indices <- function(n, B, seed) {
    stopifnot(length(n) == 1L, n >= 2L, length(B) == 1L, B >= 1L)

    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
             sample.kind = "Rejection")
    lapply(seq_len(B), function(b) sample.int(n, n, replace = TRUE))
}

score_path <- function(x_val, eta_grid) {
    m_vec <- dim(x_val)[1:2]
    n_here <- dim(x_val)[3]
    rows <- vector("list", length(eta_grid))
    prev_fit <- NULL
    prev_lambda <- NULL
    for (i in seq_along(eta_grid)) {
        lambda <- lambda_from_eta(m_vec, n_here, eta_grid[i])
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
            error = function(e) na_score_row()
        )
        if (all(is.finite(scored$ebic_joint))) {
            prev_fit <- fit
            prev_lambda <- lambda
        } else {
            prev_fit <- NULL
            prev_lambda <- NULL
        }
        rows[[i]] <- data.frame(
            group = "boot", arm = "boot", eta_index = i, eta = eta_grid[i],
            scored, stringsAsFactors = FALSE
        )
    }
    do.call(rbind, rows)
}

select_joint_index <- function(path, gamma = GAMMA_REPORT) {
    sub <- path[abs(path$gamma - gamma) < 1e-12, ]
    sub <- sub[order(sub$eta_index), ]
    valid <- which(is.finite(sub$ebic_joint))
    if (!length(valid)) {
        return(NA_integer_)
    }
    as.integer(sub$eta_index[valid[which.min(sub$ebic_joint[valid])]])
}

pair_partial_correlations <- function(Omega, pairs = BOOT_PAIRS) {
    stopifnot(is.matrix(Omega), nrow(Omega) == ncol(Omega), all(diag(Omega) > 0))
    do.call(rbind, lapply(names(pairs), function(name) {
        ij <- pairs[[name]]
        w <- Omega[ij[1], ij[2]]
        data.frame(
            pair = name,
            pc = -w / sqrt(Omega[ij[1], ij[1]] * Omega[ij[2], ij[2]]),
            selected = abs(w) > SUPPORT_TOL,
            stringsAsFactors = FALSE
        )
    }))
}

percentile_interval <- function(x, level = 0.95) {
    alpha <- (1 - level) / 2
    unname(stats::quantile(x, c(alpha, 1 - alpha), type = 7, na.rm = TRUE))
}

summarise_draws <- function(draws, reference = "CN", contrast = "Dementia",
                            level = 0.95) {
    ok <- draws[draws$ok, ]
    estimate_of <- function(group, pair) {
        e <- ok$pc[ok$b == 0L & ok$group == group & ok$pair == pair]
        if (length(e) == 1L) e else NA_real_
    }
    boot <- ok[ok$b >= 1L, ]
    keys <- unique(boot[, c("group", "pair")])
    per_group <- do.call(rbind, lapply(seq_len(nrow(keys)), function(k) {
        d <- boot[boot$group == keys$group[k] & boot$pair == keys$pair[k], ]
        ci <- percentile_interval(d$pc, level)
        data.frame(
            group = keys$group[k], pair = keys$pair[k],
            estimate = estimate_of(keys$group[k], keys$pair[k]),
            lower = ci[1], upper = ci[2],
            selected_freq = mean(d$selected), n_boot = nrow(d),
            stringsAsFactors = FALSE
        )
    }))
    diffs <- do.call(rbind, lapply(unique(boot$pair), function(pair) {
        a <- boot[boot$group == contrast & boot$pair == pair, c("b", "pc")]
        r <- boot[boot$group == reference & boot$pair == pair, c("b", "pc")]
        m <- merge(a, r, by = "b", suffixes = c("_contrast", "_reference"))
        d <- m$pc_contrast - m$pc_reference
        ci <- percentile_interval(d, level)
        data.frame(
            group = paste0(contrast, "-", reference), pair = pair,
            estimate = estimate_of(contrast, pair) - estimate_of(reference, pair),
            lower = ci[1], upper = ci[2],
            selected_freq = NA_real_, n_boot = length(d),
            stringsAsFactors = FALSE
        )
    }))
    out <- rbind(per_group, diffs)
    rownames(out) <- NULL
    out
}

run_resample <- function(Y, X, idx, eta_grid, gamma = GAMMA_REPORT,
                         pairs = BOOT_PAIRS, refit = FALSE) {
    t0 <- Sys.time()
    failed_grid <- ""
    out <- tryCatch(
        {
            resid <- residualize(Y[idx, , , drop = FALSE], X[idx, , drop = FALSE])
            x_val <- aperm(resid, c(2, 3, 1))
            path <- score_path(x_val, eta_grid)

            failed <- sort(unique(path$eta_index[!is.finite(path$ebic_joint)]))
            if (length(failed)) {
                failed_grid <- paste(failed, collapse = ",")
                stop(sprintf("joint score is not finite at grid points %s", failed_grid))
            }
            best <- select_joint_index(path, gamma)
            if (is.na(best)) stop("결합 점수가 모든 격자점에서 결측이다")
            sel <- selected_fit_at(x_val, eta_grid, best, "joint", n_sweep = N_SWEEP)
            scored <- path[path$eta_index == best & abs(path$gamma - gamma) < 1e-12, ]
            if (!candidate_agrees(scored, sel)) {
                stop(sprintf("the refit at grid point %d differs from the scored candidate", best))
            }
            pcs <- pair_partial_correlations(sel$Omega_list[[1]], pairs)
            pcs$pc_refit <- NA_real_
            refit_converged <- NA
            refit_iter <- NA_integer_
            if (refit) {
                support <- lapply(sel$path_fit, function(M) abs(M) > SUPPORT_TOL)
                ref <- refit_support_mle(x_val, support)
                pcs$pc_refit <- pair_partial_correlations(ref[[1]], pairs)$pc
                refit_converged <- isTRUE(attr(ref, "converged"))
                refit_iter <- as.integer(attr(ref, "iter"))
            }
            list(
                ok = TRUE, error = "", best_index = best,
                edges = as.integer(sel$path_edge_counts), q_u = as.integer(sel$q_u),
                pcs = pcs, refit_converged = refit_converged, refit_iter = refit_iter
            )
        },
        error = function(e) {
            list(
                ok = FALSE, error = conditionMessage(e), best_index = NA_integer_,
                edges = c(NA_integer_, NA_integer_), q_u = NA_integer_,
                pcs = data.frame(
                    pair = names(pairs), pc = NA_real_, selected = NA, pc_refit = NA_real_,
                    stringsAsFactors = FALSE
                ),
                refit_converged = NA, refit_iter = NA_integer_
            )
        }
    )
    out$failed_grid <- failed_grid
    out$elapsed_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    out
}

summarise_family <- function(draws, value = "pc",
                             contrasts = list(c("Dementia", "CN"), c("MCI", "CN")),
                             level = 0.95) {
    ok <- draws[draws$ok & is.finite(draws[[value]]), ]
    pairs <- unique(ok$pair)
    estimate_of <- function(group, pair) {
        e <- ok[[value]][ok$b == 0L & ok$group == group & ok$pair == pair]
        if (length(e) == 1L) e else NA_real_
    }
    boot <- ok[ok$b >= 1L, ]
    groups <- unique(boot$group)
    per_group <- do.call(rbind, lapply(groups, function(g) do.call(rbind, lapply(pairs, function(pair) {
        d <- boot[boot$group == g & boot$pair == pair, ]
        ci <- percentile_interval(d[[value]], level)
        data.frame(
            kind = "group", value = value, contrast = g, pair = pair,
            estimate = estimate_of(g, pair), lower = ci[1], upper = ci[2],
            selected_freq = mean(d$selected), n_boot = nrow(d),
            p_boot = NA_real_, q_bh = NA_real_, excludes_zero = NA,
            stringsAsFactors = FALSE
        )
    }))))
    diffs <- do.call(rbind, lapply(contrasts, function(cg) {
        if (!all(cg %in% groups)) {
            return(NULL)
        }
        rows <- do.call(rbind, lapply(pairs, function(pair) {
            a <- boot[boot$group == cg[1] & boot$pair == pair, c("b", value)]
            r <- boot[boot$group == cg[2] & boot$pair == pair, c("b", value)]
            m <- merge(a, r, by = "b", suffixes = c("_a", "_r"))
            d <- m[[paste0(value, "_a")]] - m[[paste0(value, "_r")]]
            ci <- percentile_interval(d, level)
            B <- length(d)
            p <- min(1, 2 * min((1 + sum(d <= 0)) / (B + 1), (1 + sum(d >= 0)) / (B + 1)))
            data.frame(
                kind = "difference", value = value, contrast = paste0(cg[1], "-", cg[2]), pair = pair,
                estimate = estimate_of(cg[1], pair) - estimate_of(cg[2], pair),
                lower = ci[1], upper = ci[2], selected_freq = NA_real_, n_boot = B,
                p_boot = p, q_bh = NA_real_, excludes_zero = ci[1] > 0 || ci[2] < 0,
                stringsAsFactors = FALSE
            )
        }))
        rows$q_bh <- stats::p.adjust(rows$p_boot, method = "BH")
        rows
    }))
    out <- rbind(per_group, diffs)
    rownames(out) <- NULL
    out
}

candidate_agrees <- function(scored, sel, tol = 1e-8) {
    if (nrow(scored) != 1L) {
        return(FALSE)
    }
    identical(as.integer(scored$q_u), as.integer(sel$q_u)) &&
        isTRUE(all.equal(c(scored$u_roi, scored$u_time), as.numeric(sel$u), tolerance = tol)) &&
        isTRUE(all.equal(scored$joint_ll_gain, sel$ll_gain, tolerance = tol))
}

run_bootstrap <- function(opts, root) {
    BLAS_ACTUAL <- load_adni_dependencies()
    code_paths <- list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)
    names(code_paths) <- basename(code_paths)
    groups <- split_csv(opts$groups)
    arm <- opts$arm
    B <- as.integer(opts$B)
    seed <- as.integer(opts$seed)
    workers <- as.integer(opts$workers)
    with_identity <- identical(opts$identity, "1")
    out_dir <- opts$out
    stopifnot(opts$pairs %in% c("reported", "homotopic"), opts$refit %in% c("0", "1"))
    pairs_list <- if (identical(opts$pairs, "homotopic")) homotopic_pairs(read_aal_labels(opts$`atlas-labels`)) else BOOT_PAIRS
    do_refit <- identical(opts$refit, "1")
    stopifnot(
        all(groups %in% names(GROUP_SEED_OFFSET)), arm %in% names(ARM_SPECS),
        !is.na(B), B >= 0L, !is.na(seed), !is.na(workers), workers >= 1L,
        with_identity || B >= 1L
    )

    code_sha256_at_start <- code_fingerprints(code_paths)

    aggregated_dir <- if (nzchar(opts$aggregated)) opts$aggregated else file.path(root, "data/processed/antspy/aggregated", opts$option)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    log_path <- file.path(out_dir, "progress.log")
    if (file.exists(log_path)) file.remove(log_path)

    np <- reticulate::import("numpy")
    data_by_group <- lapply(groups, function(g) load_arm(g, arm, aggregated_dir, np))
    names(data_by_group) <- groups
    eta_grid <- make_eta_grid(N_ETA_DEFAULT)

    jobs <- list()
    for (g in groups) {
        n_g <- dim(data_by_group[[g]]$Y)[1]
        if (with_identity) {
            jobs[[length(jobs) + 1L]] <- list(group = g, b = 0L, idx = seq_len(n_g))
        }
        if (B >= 1L) {
            draws_g <- bootstrap_indices(n_g, B, seed + GROUP_SEED_OFFSET[[g]])
            for (b in seq_len(B)) {
                jobs[[length(jobs) + 1L]] <- list(group = g, b = b, idx = draws_g[[b]])
            }
        }
    }

    started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z", tz = "UTC")
    t_start <- Sys.time()
    cat(sprintf("작업 %d 개, 워커 %d 개\n", length(jobs), workers))

    results <- parallel::mclapply(jobs, function(job) {
        d <- data_by_group[[job$group]]
        r <- run_resample(d$Y, d$X, job$idx, eta_grid, pairs = pairs_list, refit = do_refit)
        cat(sprintf(
            "%s b=%d ok=%s best=%s %.1fs\n", job$group, job$b, r$ok,
            r$best_index, r$elapsed_sec
        ), file = log_path, append = TRUE)
        data.frame(
            group = job$group, b = job$b, r$pcs,
            ok = r$ok, best_index = r$best_index,
            edges_roi = r$edges[1], edges_second = r$edges[2], q_u = r$q_u,
            refit_converged = r$refit_converged, refit_iter = r$refit_iter,
            failed_grid = r$failed_grid,
            elapsed_sec = r$elapsed_sec, error = r$error,
            stringsAsFactors = FALSE
        )
    }, mc.cores = workers, mc.preschedule = FALSE)

    results <- Map(function(x, job) {
        if (is.data.frame(x)) {
            return(x)
        }
        message <- if (inherits(x, "try-error")) trimws(as.character(x)) else "fork returned no result"
        data.frame(
            group = job$group, b = job$b, pair = names(pairs_list), pc = NA_real_, selected = NA,
            pc_refit = NA_real_, ok = FALSE, best_index = NA_integer_, edges_roi = NA_integer_,
            edges_second = NA_integer_, q_u = NA_integer_, refit_converged = NA, refit_iter = NA_integer_,
            failed_grid = "", elapsed_sec = NA_real_, error = message,
            stringsAsFactors = FALSE
        )
    }, results, jobs)
    draws <- do.call(rbind, results)
    draws_csv <- file.path(out_dir, "bootstrap_draws.csv")
    utils::write.csv(draws, draws_csv, row.names = FALSE)

    summary_csv <- file.path(out_dir, "bootstrap_summary.csv")
    if (B >= 1L && all(c("CN", "Dementia") %in% groups)) {
        utils::write.csv(summarise_draws(draws), summary_csv, row.names = FALSE)
    }

    if (B >= 1L) {
        family <- summarise_family(draws, "pc")
        if (do_refit) family <- rbind(family, summarise_family(draws, "pc_refit"))
        utils::write.csv(family, file.path(out_dir, "bootstrap_family_summary.csv"), row.names = FALSE)
    }

    run_meta <- list(
        stage = "adni-bootstrap-group-pc",
        started = started,
        elapsed_sec = round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1),
        blas_threads = BLAS_ACTUAL,
        option = opts$option,
        groups = groups,
        arm = arm,
        n_by_group = vapply(data_by_group, function(d) dim(d$Y)[1], integer(1)),
        B = B,
        identity = with_identity,
        seed = seed,
        group_seed_offset = as.list(GROUP_SEED_OFFSET[groups]),
        rng_kind = RNGkind(),
        workers = workers,
        eta_grid = eta_grid,
        gamma = GAMMA_REPORT,
        criterion = "joint",
        pairs_mode = opts$pairs,
        pairs = pairs_list,
        refit = do_refit,
        aggregated_dir = aggregated_dir,
        n_refit_not_converged = sum(draws$refit_converged[!duplicated(draws[, c("group", "b")])] %in% FALSE),
        n_sweep = N_SWEEP,
        support_tol = SUPPORT_TOL,
        joint_u_cap = JOINT_U_CAP,
        joint_boundary_margin = JOINT_BOUNDARY_MARGIN,
        n_jobs = length(jobs),
        n_failed = sum(!draws$ok[!duplicated(draws[, c("group", "b")])]),
        n_grid_failed = sum(nzchar(draws$failed_grid[!duplicated(draws[, c("group", "b")])])),
        code_sha256 = code_sha256_at_start,
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
        file.path(out_dir, "RUN_bootstrap.json")
    )
    cat(sprintf("끝: %s, 실패 %d\n", draws_csv, run_meta$n_failed))
}
