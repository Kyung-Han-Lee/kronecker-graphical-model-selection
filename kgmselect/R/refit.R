.partial_whiten_gram <- function(x_val, Omega_list, k) {
    dimen <- dim(x_val)
    K <- length(Omega_list)
    m_vec <- dimen[1:K]
    n <- dimen[K + 1L]

    if (K == 2L) {
        other_Omega <- Omega_list[[if (k == 1L) 2L else 1L]]
        S_mat <- matrix(0, m_vec[k], m_vec[k])
        for (i in seq_len(n)) {
            obs <- matrix(x_val[, , i], m_vec[1L], m_vec[2L])
            S_mat <- S_mat + if (k == 1L) {
                tcrossprod(obs %*% other_Omega, obs)
            } else {
                crossprod(obs, other_Omega %*% obs)
            }
        }
        S_mat <- (S_mat + t(S_mat)) / 2
        return((S_mat / n) * (m_vec[k] / prod(m_vec)))
    }

    S_mat <- .tensor_partial_gram_sum(x_val, Omega_list, k)
    (S_mat / n) * (m_vec[k] / prod(m_vec))
}

.zero_pairs <- function(support_k) {
    p <- nrow(support_k)
    keep <- support_k
    diag(keep) <- TRUE
    idx <- which(!keep & upper.tri(keep), arr.ind = TRUE)
    if (nrow(idx) == 0L) return(NULL)
    rbind(idx, idx[, c(2L, 1L), drop = FALSE])
}

.refit_support_mle <- function(x_val, support_list, max_iter = 50L,
                              tol = 1e-5) {
    dimen <- dim(x_val)
    K <- length(support_list)
    m_vec <- dimen[1:K]
    stopifnot(length(dimen) == K + 1L)

    Omega <- lapply(m_vec, function(mm) diag(mm))
    zeros <- lapply(support_list, .zero_pairs)

    canon_diff <- function(A, B) {
        d <- 0
        for (k in seq_along(A)) {
            a <- A[[k]] / max(norm(A[[k]], type = "F"), .Machine$double.eps)
            b <- B[[k]] / max(norm(B[[k]], type = "F"), .Machine$double.eps)
            d <- d + norm(a - b, type = "F")
        }
        d
    }

    converged <- FALSE
    iter_used <- 0L
    for (iter in seq_len(max_iter)) {
        Omega_old <- Omega
        iter_used <- iter
        for (k in seq_len(K)) {
            S_k <- .partial_whiten_gram(x_val, Omega, k)
            fit <- suppressWarnings(
                glasso::glasso(S_k, rho = 0, zero = zeros[[k]])
            )
            Om_k <- (fit$wi + t(fit$wi)) / 2

            if (!is.null(zeros[[k]])) Om_k[zeros[[k]]] <- 0
            eig_min <- min(eigen(Om_k, symmetric = TRUE,
                                 only.values = TRUE)$values)
            if (!(eig_min > 0)) {
                stop(sprintf(
                    "refit produced a non-PD factor at mode %d (min eig %.3e)",
                    k, eig_min
                ))
            }
            Omega[[k]] <- Om_k
        }

        if (K >= 2L) {
            s <- norm(Omega[[1]], type = "F")
            if (s > 0) {
                Omega[[1]] <- Omega[[1]] / s
                Omega[[2]] <- Omega[[2]] * s
            }
        }

        if (canon_diff(Omega_old, Omega) < tol) {
            converged <- TRUE
            break
        }
    }

    attr(Omega, "iter") <- iter_used
    attr(Omega, "converged") <- converged
    Omega
}

.compute_refit_ebic <- function(x_val, Omega_list, gamma = 0.2,
                                      tol = 1e-6, max_iter = 50L,
                                      refit_tol = 1e-5) {
    K <- length(Omega_list)
    n <- dim(x_val)[length(dim(x_val))]

    support_list <- lapply(Omega_list, function(M) {
        S <- abs(M) > tol
        diag(S) <- FALSE
        S
    })
    df <- sum(vapply(
        support_list, function(S) sum(S & upper.tri(S)), integer(1)
    ))

    Omega_refit <- .refit_support_mle(
        x_val, support_list, max_iter = max_iter, tol = refit_tol
    )
    ll <- as.numeric(.compute_loglik(
        x_val, Omega_refit, sqrt_method = if (K == 2L) "none" else "eigen"
    ))

    p_total <- prod(vapply(Omega_list, nrow, integer(1)))
    ebic_val <- -2 * ll + df * log(n) + 4 * gamma * df * log(p_total)

    list(
        ebic = ebic_val,
        loglik = ll,
        df = df,
        p_total = p_total,
        refit_iter = attr(Omega_refit, "iter"),
        refit_conv = attr(Omega_refit, "converged"),
        Omega_refit = Omega_refit
    )
}

kgm_refit <- function(Y, supports, max_iter = 200L, tol = 1e-5) {
    .check_data(Y)
    .refit_support_mle(Y, supports, max_iter = max_iter, tol = tol)
}
