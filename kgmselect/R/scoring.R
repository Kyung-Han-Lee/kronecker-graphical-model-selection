.compute_df <- function(Omega_list, tol = 1e-6) {
    sum(vapply(Omega_list, function(Om) {
        off <- Om
        diag(off) <- 0
        as.integer(sum(abs(off[upper.tri(off)]) > tol))
    }, integer(1)))
}

kgm_score <- function(path, Y,
                      criterion = c("relaxation", "scale", "ebic", "refit"),
                      gamma = 0.5, tol = 1e-6, boundary_margin = 1e-8,
                      u_cap = 10, gradient = c("numeric", "analytic"),
                      optim_control = list(maxit = 100L, factr = 1e10),
                      refit_control = list(max_iter = 200L, tol = 1e-5)) {
    .check_data(Y)
    .check_path(path, Y)
    criterion <- match.arg(criterion)
    gradient <- match.arg(gradient)
    if (!is.numeric(gamma) || !length(gamma) || any(!is.finite(gamma)) ||
            any(gamma < 0)) {
        stop("gamma must contain finite non-negative values")
    }
    dims <- dim(Y)
    M <- length(dims) - 1L
    n <- dims[M + 1L]
    log_p <- log(prod(dims[seq_len(M)]))
    rows <- lapply(seq_along(path$eta), function(k) {
        edges <- .compute_df(path$fits[[k]], tol = tol)

        res <- if (criterion != "refit" &&
                   !.is_numerically_pd(path$fits[[k]])) {
            .unscored("not_pd", "a fitted factor is not numerically positive definite")
        } else {
            tryCatch(
                .score_candidate(path, Y, k, criterion, tol, boundary_margin,
                                 u_cap, gradient, optim_control, refit_control),
                error = function(e) {
                    stop(sprintf("candidate %d: %s", k, conditionMessage(e)),
                         call. = FALSE)
                })
        }

        ebic <- -2 * res$loglik + res$df * log(n) + 4 * gamma * res$df * log_p
        data.frame(k = k, eta = path$eta[k], gamma = gamma,
                   criterion = criterion, edges = edges, df = res$df,
                   loglik = res$loglik, tau = res$tau,
                   ebic = ebic, status = res$status, message = res$message,
                   stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    out <- out[order(out$gamma, out$k), , drop = FALSE]
    rownames(out) <- NULL
    out
}

.score_candidate <- function(path, Y, k, criterion, tol, boundary_margin,
                             u_cap, gradient, optim_control, refit_control) {
    fit <- path$fits[[k]]
    edges <- .compute_df(fit, tol = tol)
    if (criterion == "ebic") {
        ll <- .compute_loglik(Y, fit, tau_correct = FALSE, sqrt_method = "none")
        return(list(loglik = as.numeric(ll), df = edges, tau = 1, status = "ok",
                    message = ""))
    }
    if (criterion == "scale") {
        ll <- .compute_loglik(Y, fit, tau_correct = TRUE, sqrt_method = "none")
        return(list(loglik = as.numeric(ll), df = edges + 1L,
                    tau = attr(ll, "tau_global"), status = "ok", message = ""))
    }
    if (criterion == "relaxation") {
        prof <- .profile_relaxation(Y, fit, .relaxation_deltas(path, k, tol),
                                    boundary_margin = boundary_margin,
                                    u_cap = u_cap, optim_control = optim_control,
                                    gradient = gradient)
        return(list(loglik = prof$loglik, df = edges + 1L,
                    tau = prof$tau_global, status = "ok", message = ""))
    }
    rf <- tryCatch(
        .compute_refit_ebic(Y, fit, gamma = 0, tol = tol,
                            max_iter = refit_control$max_iter,
                            refit_tol = refit_control$tol),
        error = function(e) if (.is_refit_failure(e)) NULL else stop(e))
    if (is.null(rf)) {
        return(.unscored("refit_failed",
                         "the refit reached a factor that is not positive definite"))
    }
    if (!isTRUE(rf$refit_conv)) {
        return(.unscored("nonconverged",
                         sprintf("refit did not converge in %d sweeps",
                                 as.integer(rf$refit_iter))))
    }
    list(loglik = rf$loglik, df = rf$df, tau = NA_real_, status = "ok", message = "")
}

.unscored <- function(status, message) {
    list(loglik = NA_real_, df = NA_integer_, tau = NA_real_, status = status,
         message = message)
}

.is_numerically_pd <- function(fit) {
    all(vapply(fit, function(A) {
        !inherits(try(chol(A), silent = TRUE), "try-error")
    }, logical(1)))
}

.is_refit_failure <- function(e) {
    startsWith(conditionMessage(e), "refit produced a non-PD factor")
}
