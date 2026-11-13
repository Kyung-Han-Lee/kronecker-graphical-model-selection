.compute_loglik <- function(x_val, Omega_list, tau_correct = FALSE,
                                  sqrt_method = c("eigen", "expm", "none"),
                                  check_pd = TRUE) {
    shape <- .validate_kronecker_inputs(x_val, Omega_list)
    sqrt_method <- match.arg(sqrt_method)
    if (isTRUE(check_pd)) {
        for (m in seq_along(Omega_list)) {
            if (inherits(try(chol(Omega_list[[m]]), silent = TRUE),
                         "try-error")) {
                stop(sprintf("mode %d factor is not positive definite", m))
            }
        }
    }

    Omega_sqrt <- if (sqrt_method == "none") NULL else
        lapply(Omega_list, function(Om) {
            if (sqrt_method == "expm") {
                s <- expm::sqrtm(Om)
                if (max(abs(Im(s))) > 1e-8) {
                    warning(sprintf(
                        "sqrtm imag part %.2e exceeds tol; Omega may not be SPD",
                        max(abs(Im(s)))))
                }
                return(Re(s))
            }
            e <- eigen((Om + t(Om)) / 2, symmetric = TRUE)
            if (min(e$values) < 0) {
                warning(sprintf(
                    "eigen min value %.2e is negative; Omega may not be SPD",
                    min(e$values)))
            }
            e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
        })

    logdet <- vapply(Omega_list, function(Om) {
        dd <- determinant(Om, logarithm = TRUE)
        if (as.numeric(dd$sign) <= 0) {
            warning(sprintf("determinant sign=%s (expected +1 for SPD)",
                            as.character(as.numeric(dd$sign))))
        }
        as.numeric(dd$modulus)
    }, numeric(1))
    weights <- prod(shape$p_vec) / shape$p_vec

    if (sqrt_method == "none") {
        for (m in seq_along(Omega_list)) {
            if (inherits(try(chol(Omega_list[[m]]), silent = TRUE),
                         "try-error")) {
                stop(sprintf(paste("sqrt_method = 'none' requires positive",
                                   "definite factors; mode %d is not"), m))
            }
        }
        trace_term <- .kronecker_trace(x_val, Omega_list)
    } else {
        trace_term <- .kronecker_squared_norm_trace(x_val, Omega_sqrt)
    }

    if (tau_correct) {
        if (!(trace_term > 0)) {
            stop(sprintf(
                "tau_correct requires tr(S Theta) > 0, got %.3e", trace_term
            ))
        }
        p_total <- prod(shape$p_vec)
        tau_global <- p_total / trace_term
        loglik <- (shape$n / 2) *
            (p_total * log(tau_global) + sum(logdet * weights) - p_total)
        attr(loglik, "tau_global") <- tau_global
        return(loglik)
    }

    (shape$n / 2) * (sum(logdet * weights) - trace_term)
}
