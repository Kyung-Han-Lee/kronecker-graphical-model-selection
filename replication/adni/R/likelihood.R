.validate_kronecker_inputs <- function(x_val, Omega_list) {
    dims <- dim(x_val)
    M <- length(Omega_list)
    if (M < 1L || is.null(dims) || length(dims) != M + 1L) {
        stop("data must have one sample dimension after the tensor modes")
    }
    p_vec <- dims[seq_len(M)]
    for (m in seq_len(M)) {
        Om <- Omega_list[[m]]
        if (!is.matrix(Om) || nrow(Om) != ncol(Om)) {
            stop(sprintf("mode %d factor must be a square matrix", m))
        }
        if (nrow(Om) != p_vec[m]) {
            stop(sprintf("mode %d factor dimension does not match the data", m))
        }
    }
    list(M = M, p_vec = p_vec, n = dims[M + 1L])
}

.mode_multiply <- function(X, A, mode) {
    dims <- dim(X)
    if (is.null(dims) || length(mode) != 1L || mode < 1L ||
            mode > length(dims) || ncol(A) != dims[mode]) {
        stop("invalid tensor mode multiplication")
    }
    perm <- c(mode, setdiff(seq_along(dims), mode))
    unfolded <- matrix(aperm(X, perm), nrow = dims[mode])
    folded <- array(A %*% unfolded,
                    dim = c(nrow(A), dims[setdiff(seq_along(dims), mode)]))
    aperm(folded, order(perm))
}

.sum_tensor_blocks <- function(x_val, evaluate) {
    dims <- dim(x_val)
    n <- tail(dims, 1L)
    p_vec <- head(dims, -1L)
    stopifnot(length(dims) >= 2L, n >= 1L)
    chunk <- max(1L, min(n, floor(65536 / prod(p_vec))))
    X <- matrix(x_val, nrow = prod(p_vec), ncol = n)
    total <- NULL
    for (first in seq.int(1L, n, by = chunk)) {
        indices <- seq.int(first, min(n, first + chunk - 1L))
        Y <- array(X[, indices, drop = FALSE], c(p_vec, length(indices)))
        value <- evaluate(Y)
        total <- if (is.null(total)) value else total + value
    }
    total
}

.tensor_partial_gram_sum <- function(x_val, factors, mode) {
    M <- length(factors)
    others <- setdiff(seq_len(M), mode)
    .sum_tensor_blocks(x_val, function(Y) {
        Z <- Y
        for (m in others) Z <- .mode_multiply(Z, factors[[m]], m)
        perm <- c(mode, others, M + 1L)
        Yk <- matrix(aperm(Y, perm), nrow = dim(Y)[mode])
        Zk <- matrix(aperm(Z, perm), nrow = dim(Y)[mode])
        G <- tcrossprod(Yk, Zk)
        (G + t(G)) / 2
    })
}

.kronecker_trace_m2 <- function(x_val, factor_list) {
    dims <- dim(x_val)
    H <- matrix(0, dims[2L], dims[2L])
    for (i in seq_len(dims[3L])) {
        Y <- matrix(x_val[, , i], dims[1L], dims[2L])
        H <- H + crossprod(Y, factor_list[[1L]] %*% Y)
    }
    sum(factor_list[[2L]] * H) / dims[3L]
}

.kronecker_trace <- function(x_val, factor_list) {
    shape <- .validate_kronecker_inputs(x_val, factor_list)
    if (shape$M == 2L) return(.kronecker_trace_m2(x_val, factor_list))
    .sum_tensor_blocks(x_val, function(Y) {
        Z <- Y
        for (m in seq_len(shape$M)) {
            Z <- .mode_multiply(Z, factor_list[[m]], m)
        }
        sum(Y * Z)
    }) / shape$n
}

.kronecker_squared_norm_trace <- function(x_val, factor_sqrt_list) {
    shape <- .validate_kronecker_inputs(x_val, factor_sqrt_list)
    if (shape$M > 2L) {
        return(.sum_tensor_blocks(x_val, function(Z) {
            for (m in seq_len(shape$M)) {
                Z <- .mode_multiply(Z, factor_sqrt_list[[m]], m)
            }
            sum(Z^2)
        }) / shape$n)
    }
    X <- matrix(x_val, nrow = prod(shape$p_vec), ncol = shape$n)
    acc <- 0
    for (i in seq_len(shape$n)) {
        Zi <- array(X[, i], dim = shape$p_vec)
        for (m in seq_len(shape$M)) {
            Zi <- .mode_multiply(Zi, factor_sqrt_list[[m]], m)
        }
        acc <- acc + sum(Zi^2)
    }
    acc / shape$n
}

compute_loglik_tlasso <- function(x_val, Omega_list, tau_correct = FALSE,
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

compute_df_tlasso <- function(Omega_list, tol = 1e-6) {
    sum(vapply(Omega_list, function(Om) {
        off <- Om
        diag(off) <- 0
        as.integer(sum(abs(off[upper.tri(off)]) > tol))
    }, integer(1)))
}

compute_ebic_tlasso <- function(x_val, Omega_list, gamma = 0.2, tol = 1e-6,
                                tau_correct = FALSE,
                                sqrt_method = c("eigen", "expm", "none"),
                                check_pd = TRUE) {
    shape <- .validate_kronecker_inputs(x_val, Omega_list)
    sqrt_method <- match.arg(sqrt_method)
    ll <- compute_loglik_tlasso(x_val, Omega_list, tau_correct = tau_correct,
                                sqrt_method = sqrt_method, check_pd = check_pd)
    tau_global <- if (tau_correct) attr(ll, "tau_global") else NULL
    ll <- as.numeric(ll)
    df <- compute_df_tlasso(Omega_list, tol = tol)
    if (tau_correct) df <- df + 1L

    per_mode_s <- vapply(Omega_list, function(Om) {
        as.integer(sum(abs(Om[upper.tri(Om)]) > tol))
    }, integer(1))
    p_total <- prod(shape$p_vec)
    penalty <- df * log(p_total)

    list(
        ebic = -2 * ll + df * log(shape$n) + 4 * gamma * penalty,
        loglik = ll,
        df = df,
        penalty = penalty,
        per_mode_edges = per_mode_s,
        p_total = p_total,
        tau_global = tau_global
    )
}
