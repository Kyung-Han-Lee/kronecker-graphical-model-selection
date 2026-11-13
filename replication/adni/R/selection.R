selected_fit_at <- function(x_val, eta_grid, best_index, criterion,
                            n_sweep = N_SWEEP) {
    stopifnot(
        length(dim(x_val)) == 3L,
        length(eta_grid) >= 1L,
        all(is.finite(eta_grid)), all(eta_grid > 0),
        length(best_index) == 1L, best_index == as.integer(best_index),
        best_index >= 1L, best_index <= length(eta_grid),
        criterion %in% c("plug", "prof", "joint")
    )
    best_index <- as.integer(best_index)
    m_vec <- dim(x_val)[1:2]
    n_here <- dim(x_val)[3]
    fit_one <- function(index) {
        lambda <- lambda_from_eta(m_vec, n_here, eta_grid[index])
        fit <- lapply(
            fit_tggm(x_val, lambda, n_sweep = n_sweep), function(M) M
        )
        list(fit = fit, lambda = lambda)
    }

    current <- fit_one(best_index)
    path_edge_counts <- vapply(
        current$fit,
        function(M) sum(abs(M[upper.tri(M)]) > SUPPORT_TOL),
        integer(1)
    )
    if (criterion != "joint") {
        return(list(
            Omega_list = current$fit,
            path_fit = current$fit,
            path_edge_counts = path_edge_counts,
            u = c(0, 0), q_u = 0L, tau_global = NA_real_, ll_gain = 0
        ))
    }

    if (best_index == 1L) {
        delta <- lapply(
            current$fit, function(M) matrix(0, nrow(M), ncol(M))
        )
    } else {
        previous <- fit_one(best_index - 1L)
        dirs <- Map(
            function(now, before, lambda_now, lambda_before) {
                projected_relaxation_direction(
                    now, before, lambda_now, lambda_before,
                    tol = SUPPORT_TOL
                )
            },
            current$fit, previous$fit,
            as.list(current$lambda), as.list(previous$lambda)
        )
        delta <- lapply(dirs, `[[`, "delta")
    }
    joint <- profile_joint_relaxation(
        x_val, current$fit, delta,
        boundary_margin = JOINT_BOUNDARY_MARGIN, u_cap = JOINT_U_CAP
    )
    list(
        Omega_list = joint$Omega_list,
        path_fit = current$fit,
        path_edge_counts = path_edge_counts,
        u = joint$u,
        q_u = joint$q_u,
        tau_global = joint$tau_global,
        ll_gain = joint$ll_gain
    )
}

verdict_table <- function(path) {
    out <- list()
    for (g in unique(path$group)) {
        for (a in unique(path$arm)) {
            for (gam in unique(path$gamma)) {
                sub <- path[path$group == g & path$arm == a &
                    path$gamma == gam, ]
                if (!nrow(sub)) next
                sub <- sub[order(sub$eta_index), ]
                for (crit in c("plug", "prof", "joint", "joint_uncounted")) {
                    values <- sub[[paste0("ebic_", crit)]]
                    valid <- which(is.finite(values))
                    if (!length(valid)) next
                    k <- valid[which.min(values[valid])]
                    out[[length(out) + 1L]] <- data.frame(
                        group = g, arm = a, gamma = gam, criterion = crit,
                        best_eta = sub$eta[k], best_index = k,
                        n_eta = nrow(sub), n_valid = length(valid),
                        edges_roi = sub$edges_roi[k],
                        edges_time = sub$edges_time[k],
                        is_empty = isTRUE(sub$df[k] == 0L),
                        interior = (k - 1L) %in% valid && (k + 1L) %in% valid,
                        stringsAsFactors = FALSE
                    )
                }
            }
        }
    }
    do.call(rbind, out)
}

save_selected_fits <- function(opts, root) {
    BLAS_ACTUAL <- load_adni_dependencies()
    code_paths <- list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)
    names(code_paths) <- basename(code_paths)
    option_arg <- opts$option
    groups <- split_csv(opts$groups)
    arm <- opts$arm
    gamma_pick <- as.numeric(opts$gamma)
    stopifnot(arm %in% names(ARM_SPECS), gamma_pick %in% GAMMA_GRID)

    aggregated_dir <- if (nzchar(opts$aggregated)) opts$aggregated else file.path(root, "data/processed/antspy/aggregated", option_arg)
    dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)
    np <- reticulate::import("numpy")

    read_selected_eta <- function(verdict_csvs, group, criterion) {
        for (path in verdict_csvs) {
            if (!file.exists(path)) next
            v <- utils::read.csv(path)
            hit <- v[v$group == group & v$arm == arm &
                         abs(v$gamma - gamma_pick) < 1e-12 &
                         v$criterion == criterion, ]
            if (nrow(hit) == 1L) return(hit)
            if (nrow(hit) > 1L) {
                stop(sprintf("%s 에 (%s, %s, %s) 행이 %d 개다",
                             path, group, arm, criterion, nrow(hit)))
            }
        }
        stop(sprintf("판정 CSV 어디에도 (%s, %s, gamma=%g, %s) 가 없다",
                     group, arm, gamma_pick, criterion))
    }

    verdict_csvs <- split_csv(opts$verdict)

    summary_rows <- list()
    for (group in groups) {
        d <- load_arm(group, arm, aggregated_dir, np)
        x_val <- aperm(residualize(d$Y, d$X), c(2, 3, 1))
        m_vec <- dim(x_val)[1:2]
        n_here <- dim(x_val)[3]

        for (criterion in c("plug", "prof", "joint")) {
            hit <- read_selected_eta(verdict_csvs, group, criterion)
            eta <- hit$best_eta
            eta_grid <- make_eta_grid(hit$n_eta)
            if (abs(eta_grid[hit$best_index] - eta) > 1e-12 * max(1, eta)) {
                stop(sprintf("%s/%s 의 eta와 best_index가 일치하지 않는다", group, criterion))
            }
            selected <- selected_fit_at(
                x_val, eta_grid, hit$best_index, criterion, n_sweep = N_SWEEP
            )
            fit <- selected$Omega_list
            matrix_edges <- vapply(
                fit,
                function(M) sum(abs(M[upper.tri(M)]) > SUPPORT_TOL),
                integer(1)
            )

            got <- selected$path_edge_counts
            if (!identical(as.integer(got), as.integer(c(hit$edges_roi, hit$edges_time)))) {
                stop(sprintf(
                    "%s/%s 재적합 엣지 수가 판정과 다르다: (%d, %d) vs (%d, %d)",
                    group, criterion, got[1], got[2],
                    hit$edges_roi, hit$edges_time
                ))
            }

            for (k in seq_along(fit)) {
                utils::write.csv(
                    fit[[k]],
                    file.path(opts$out, sprintf(
                        "Omega%d_%s_%s_%s_g%s.csv", k, group, arm, criterion,
                        sub("[.]", "p", format(gamma_pick))
                    )),
                    row.names = FALSE
                )
            }
            summary_rows[[length(summary_rows) + 1L]] <- data.frame(
                group = group, arm = arm, gamma = gamma_pick, criterion = criterion,
                eta = eta, n = n_here, p_roi = m_vec[1], p_second = m_vec[2],
                edges_roi = got[1], edges_second = got[2],
                matrix_edges_roi = matrix_edges[1],
                matrix_edges_second = matrix_edges[2],
                u_roi = selected$u[1], u_second = selected$u[2],
                tau_global = selected$tau_global,
                stringsAsFactors = FALSE
            )
            cat(sprintf(
                "  %-9s %-5s eta=%-8.4g edges=(%5d,%5d) u=(%.3g,%.3g)\n",
                group, criterion, eta, got[1], got[2],
                selected$u[1], selected$u[2]
            ))
        }
    }
    utils::write.csv(do.call(rbind, summary_rows),
                     file.path(opts$out, "selected_fits.csv"), row.names = FALSE)
    verdict_paths <- normalizePath(verdict_csvs[file.exists(verdict_csvs)])
    run_meta <- make_refit_run_meta(
        option = option_arg,
        groups = groups,
        arm = arm,
        gamma = gamma_pick,
        verdict_paths = verdict_paths,
        verdict_sha256 = code_fingerprints(
            stats::setNames(verdict_paths, verdict_paths)
        ),
        code_sha256 = code_fingerprints(code_paths)
    )
    writeLines(
        jsonlite::toJSON(run_meta, auto_unbox = TRUE, pretty = TRUE, digits = NA),
        file.path(opts$out, "RUN_refit_selected.json")
    )
    cat(sprintf("\n인자 행렬: %s\n", opts$out))
}
