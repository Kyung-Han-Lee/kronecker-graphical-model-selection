score_fit <- function(x_val, fit, n, prev_fit = NULL, lambda,
                      prev_lambda = NULL) {
    stopifnot(
        length(fit) == 2L, length(lambda) == 2L,
        is.null(prev_fit) == is.null(prev_lambda)
    )
    base <- lapply(c(FALSE, TRUE), function(tc) {
        compute_ebic_tlasso(
            x_val, fit,
            gamma = 0, tol = SUPPORT_TOL,
            tau_correct = tc, sqrt_method = SQRT_METHOD
        )
    })
    names(base) <- c("plug", "prof")

    df_fit <- base$plug$df
    stopifnot(
        length(base$prof$df) == 1L, base$prof$df == df_fit + 1L,
        length(base$prof$tau_global) == 1L,
        is.finite(base$prof$tau_global), base$prof$tau_global > 0
    )

    ebic_at <- function(b, gamma) {
        -2 * b$loglik + b$df * log(n) + 4 * gamma * b$penalty
    }

    dirs <- if (is.null(prev_fit)) NULL else Map(
        function(now, previous, lambda_now, lambda_previous) {
            projected_relaxation_direction(
                now, previous, lambda_now, lambda_previous,
                tol = SUPPORT_TOL
            )
        },
        fit, prev_fit, as.list(lambda), as.list(prev_lambda)
    )
    delta <- if (is.null(dirs)) {
        lapply(fit, function(M) matrix(0, nrow(M), ncol(M)))
    } else {
        lapply(dirs, `[[`, "delta")
    }
    n_masked <- if (is.null(dirs)) {
        c(NA_integer_, NA_integer_)
    } else {
        vapply(dirs, `[[`, integer(1), "n_masked")
    }
    t_joint <- system.time(joint <- profile_joint_relaxation(
        x_val, fit, delta,
        boundary_margin = JOINT_BOUNDARY_MARGIN,
        u_cap = JOINT_U_CAP
    ))[["elapsed"]]
    cond_base <- vapply(fit, spd_condition_number, numeric(1))
    cond_joint <- vapply(joint$Omega_list, spd_condition_number, numeric(1))
    df_tau <- df_fit + 1L

    joint_uncounted <- vapply(GAMMA_GRID, function(g) {
        -2 * joint$loglik + df_tau * log(n) +
            4 * g * df_tau * log(prod(dim(x_val)[1:2]))
    }, 0)
    data.frame(
        gamma = GAMMA_GRID,
        edges_roi = base$plug$per_mode_edges[1],
        edges_time = base$plug$per_mode_edges[2],
        df = df_fit,
        loglik_plug = base$plug$loglik,
        loglik_prof = base$prof$loglik,
        loglik_joint = joint$loglik,
        tau_global = base$prof$tau_global,
        tau_joint = joint$tau_global,
        q_u = joint$q_u,
        u_roi = joint$u[1],
        u_time = joint$u[2],
        u_upper_roi = joint$u_upper[1],
        u_upper_time = joint$u_upper[2],
        bound_roi = joint$bound_kind[1],
        bound_time = joint$bound_kind[2],
        cond_base_roi = cond_base[1],
        cond_base_time = cond_base[2],
        cond_joint_roi = cond_joint[1],
        cond_joint_time = cond_joint[2],
        n_masked_roi = n_masked[1],
        n_masked_time = n_masked[2],
        joint_ll_gain = joint$ll_gain,
        joint_convergence = joint$convergence,
        joint_elapsed_sec = t_joint,
        ebic_plug = vapply(GAMMA_GRID, function(g) ebic_at(base$plug, g), 0),
        ebic_prof = vapply(GAMMA_GRID, function(g) ebic_at(base$prof, g), 0),
        ebic_joint = joint_uncounted,
        ebic_joint_uncounted = joint_uncounted,
        stringsAsFactors = FALSE
    )
}

na_score_row <- function() {
    data.frame(
        gamma = GAMMA_GRID, edges_roi = NA_integer_, edges_time = NA_integer_,
        df = NA_integer_, loglik_plug = NA_real_, loglik_prof = NA_real_,
        loglik_joint = NA_real_, tau_global = NA_real_, tau_joint = NA_real_,
        q_u = NA_integer_, u_roi = NA_real_, u_time = NA_real_,
        u_upper_roi = NA_real_, u_upper_time = NA_real_,
        bound_roi = NA_character_, bound_time = NA_character_,
        cond_base_roi = NA_real_, cond_base_time = NA_real_,
        cond_joint_roi = NA_real_, cond_joint_time = NA_real_,
        n_masked_roi = NA_integer_, n_masked_time = NA_integer_,
        joint_ll_gain = NA_real_, joint_convergence = NA_integer_,
        joint_elapsed_sec = NA_real_,
        ebic_plug = NA_real_, ebic_prof = NA_real_,
        ebic_joint = NA_real_, ebic_joint_uncounted = NA_real_,
        stringsAsFactors = FALSE
    )
}
