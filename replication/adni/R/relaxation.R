.relaxation_support_mask <- function(M, tol) {
    S <- abs(M) > tol
    diag(S) <- FALSE
    S
}

.relaxation_direction_varies <- function(A, delta) {
    f_delta <- norm(delta, type = "F")
    if (!(f_delta > 0)) return(FALSE)
    A_unit <- A / norm(A, type = "F")
    D_unit <- delta / f_delta
    tangent <- D_unit - sum(A_unit * D_unit) * A_unit

    norm(tangent, type = "F") > 8 * .Machine$double.eps
}

projected_relaxation_direction <- function(Omega_now, Omega_prev, lambda_now, lambda_prev,
                              tol = 1e-6, dir_tol = 1e-12) {
    stopifnot(is.matrix(Omega_now), is.matrix(Omega_prev),
              identical(dim(Omega_now), dim(Omega_prev)),
              nrow(Omega_now) == ncol(Omega_now),
              length(lambda_now) == 1L, length(lambda_prev) == 1L,
              is.finite(lambda_now), is.finite(lambda_prev),
              lambda_now > lambda_prev,
              length(tol) == 1L, is.finite(tol), tol >= 0)
    secant <- (Omega_now - Omega_prev) / (lambda_now - lambda_prev)
    keep <- abs(Omega_now) > tol
    diag(keep) <- TRUE
    dropped <- (secant != 0) & !keep & upper.tri(secant)
    delta <- -lambda_now * secant * keep
    delta <- (delta + t(delta)) / 2
    scale_ref <- max(norm(Omega_now, type = "F"), .Machine$double.eps)
    list(delta = delta,
         active = norm(delta, type = "F") > dir_tol * scale_ref &&
             .relaxation_direction_varies(Omega_now, delta),
         n_masked = as.integer(sum(dropped)))
}

spd_condition_number <- function(A) {
    stopifnot(is.matrix(A), nrow(A) == ncol(A))
    vals <- eigen((A + t(A)) / 2, symmetric = TRUE,
                  only.values = TRUE)$values
    if (min(vals) <= 0) stop("A must be positive definite")
    as.numeric(max(vals) / min(vals))
}

relaxation_factor_info <- function(Omega, delta, boundary_margin = 1e-8,
                                u_cap = 10) {
    stopifnot(is.matrix(Omega), is.matrix(delta),
              identical(dim(Omega), dim(delta)),
              nrow(Omega) == ncol(Omega),
              length(boundary_margin) == 1L, is.finite(boundary_margin),
              boundary_margin > 0, boundary_margin < 1,
              length(u_cap) == 1L, is.finite(u_cap), u_cap > 0)
    A <- (Omega + t(Omega)) / 2
    A_chol <- try(chol(A), silent = TRUE)
    if (inherits(A_chol, "try-error")) {
        stop("Omega must be positive definite")
    }
    logdet_A <- 2 * sum(log(diag(A_chol)))
    Dl <- (delta + t(delta)) / 2
    if (!.relaxation_direction_varies(A, Dl)) {
        return(list(A = A, delta = Dl, active = FALSE,
                    logdet_A = logdet_A, mu = numeric(0), u_pd = Inf,
                    u_support = Inf, u_upper = 0, bound_kind = "none"))
    }

    eig <- eigen(A, symmetric = TRUE)
    if (min(eig$values) <= 0) stop("Omega has a non-positive eigenvalue")
    A_inv_sqrt <- eig$vectors %*% (1 / sqrt(eig$values) * t(eig$vectors))
    C <- A_inv_sqrt %*% Dl %*% A_inv_sqrt
    C <- (C + t(C)) / 2
    mu <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
    mu_min <- min(mu)
    u_pd <- if (mu_min >= 0) Inf else -1 / mu_min

    off <- upper.tri(A)
    crossing <- off & A != 0 & Dl != 0 & sign(A) != sign(Dl)
    u_each <- abs(A[crossing]) / abs(Dl[crossing])
    u_support <- if (length(u_each)) min(u_each) else Inf
    if (!(u_pd > 0) || !(u_support > 0)) {
        stop("input already sits on the PD or support boundary")
    }

    cand <- c(pd = u_pd, support = u_support, cap = u_cap)
    list(A = A, delta = Dl, active = TRUE, logdet_A = logdet_A,
         mu = mu, u_pd = u_pd,
         u_support = u_support,
         u_upper = (1 - boundary_margin) * min(cand),
         bound_kind = names(cand)[which.min(cand)])
}

relaxation_transform_factor <- function(info, u) {
    stopifnot(is.list(info), length(u) == 1L, is.finite(u), u >= 0)
    if (!isTRUE(info$active)) {
        if (!identical(as.numeric(u), 0)) {
            stop("u is fixed at zero for a constant normalized direction")
        }
    } else if (u > info$u_upper) {
        stop(sprintf("u %.8g exceeds safe bound %.8g", u, info$u_upper))
    }
    M <- info$A + u * info$delta
    M <- (M + t(M)) / 2
    if (inherits(try(chol(M), silent = TRUE), "try-error")) {
        stop("transformed factor is not positive definite")
    }
    frob <- norm(M, type = "F")
    if (!(is.finite(frob) && frob > 0)) {
        stop("transformed factor has invalid Frobenius norm")
    }
    out <- M / frob
    off <- upper.tri(M)
    base <- info$A[off]
    candidate <- M[off]
    if (!identical(candidate != 0, base != 0)) {
        stop("transformed factor changed the off-diagonal support")
    }
    active <- base != 0
    if (any(sign(candidate[active]) != sign(base[active]))) {
        stop("transformed factor changed an active off-diagonal sign")
    }
    out
}

.relaxation_masks <- function(M) {
    if (length(M) != 1L || !is.finite(M) || M < 1L || M != as.integer(M)) {
        stop("M must be a positive integer")
    }
    out <- as.matrix(expand.grid(
        rep(list(0:1), as.integer(M)), KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE))
    storage.mode(out) <- "integer"
    out
}

.joint_trace_coefficients_m2 <- function(x_val, info) {
    dims <- dim(x_val)
    p1 <- dims[1L]; p2 <- dims[2L]; n <- dims[3L]
    Q10 <- info[[1L]]$A; Q11 <- info[[1L]]$delta
    H0 <- matrix(0, p2, p2); H1 <- matrix(0, p2, p2)
    for (i in seq_len(n)) {
        Yi <- matrix(x_val[, , i], p1, p2)
        H0 <- H0 + crossprod(Yi, Q10 %*% Yi)
        H1 <- H1 + crossprod(Yi, Q11 %*% Yi)
    }

    c(sum(info[[2L]]$A * H0), sum(info[[2L]]$A * H1),
      sum(info[[2L]]$delta * H0), sum(info[[2L]]$delta * H1)) / n
}

.tensor_transform_options <- function(Y, pairs, modes) {
    terms <- list(Y)
    for (i in seq_along(modes)) {
        next_terms <- vector("list", 2L * length(terms))
        for (bit in 0:1) {
            for (j in seq_along(terms)) {
                next_terms[[bit * length(terms) + j]] <-
                    .mode_multiply(terms[[j]], pairs[[i]][[bit + 1L]], modes[i])
            }
        }
        terms <- next_terms
    }
    do.call(cbind, lapply(terms, as.vector))
}

.joint_trace_coefficients <- function(x_val, info) {
    M <- length(info)
    dims <- dim(x_val)
    if (M < 1L || is.null(dims) || length(dims) != M + 1L) {
        stop("data and relaxation information have different mode counts")
    }
    p_vec <- dims[seq_len(M)]
    for (m in seq_len(M)) {
        if (!is.list(info[[m]]) || !is.matrix(info[[m]]$A) ||
                !is.matrix(info[[m]]$delta) ||
                nrow(info[[m]]$A) != p_vec[m] ||
                !identical(dim(info[[m]]$A), dim(info[[m]]$delta))) {
            stop(sprintf("mode %d relaxation information does not match data", m))
        }
    }
    n <- dims[M + 1L]
    if (M == 2L) {
        return(.joint_trace_coefficients_m2(x_val, info))
    }
    left_modes <- seq_len(M %/% 2L)
    right_modes <- seq.int(M %/% 2L + 1L, M)
    left_pairs <- lapply(info[left_modes], function(inf) {
        list(t(inf$A), t(inf$delta))
    })
    right_pairs <- lapply(info[right_modes], function(inf) {
        list(inf$A, inf$delta)
    })
    .sum_tensor_blocks(x_val, function(Y) {
        left <- .tensor_transform_options(Y, left_pairs, left_modes)
        right <- .tensor_transform_options(Y, right_pairs, right_modes)
        as.vector(crossprod(left, right))
    }) / n
}

.joint_trace_value <- function(u, coefficients, masks = NULL) {
    if (is.null(masks)) masks <- .relaxation_masks(length(u))
    if (!is.matrix(masks) || ncol(masks) != length(u) ||
            nrow(masks) != length(coefficients)) {
        stop("trace coefficients and relaxation masks do not match")
    }
    monomials <- apply(masks, 1L, function(mask) {
        active <- which(mask == 1L)
        if (!length(active)) 1 else prod(u[active])
    })
    sum(coefficients * monomials)
}

.joint_loglik <- function(u, info, trace_coef, n, p_vec, masks) {
    if (length(u) != length(info) || length(p_vec) != length(info)) {
        stop("u, factors, and tensor dimensions have different mode counts")
    }
    U <- Map(function(inf, uu) {
        M <- inf$A + uu * inf$delta
        (M + t(M)) / 2
    }, info, u)
    if (any(vapply(U, function(M)
            inherits(try(chol(M), silent = TRUE), "try-error"),
            logical(1)))) {
        return(list(loglik = -Inf, tau_raw = NA_real_, U = U))
    }
    logdet <- vapply(U, function(M) {
        d <- determinant(M, logarithm = TRUE)
        if (as.numeric(d$sign) <= 0) return(-Inf)
        as.numeric(d$modulus)
    }, numeric(1))
    trace_term <- .joint_trace_value(u, trace_coef, masks)
    if (!(is.finite(trace_term) && trace_term > 0) ||
            any(!is.finite(logdet))) {
        return(list(loglik = -Inf, tau_raw = NA_real_, U = U))
    }
    p <- prod(p_vec)
    tau_raw <- p / trace_term
    weights <- p / p_vec
    loglik <- (n / 2) * (p * log(tau_raw) + sum(weights * logdet) - p)
    list(loglik = loglik, tau_raw = tau_raw, U = U)
}

.joint_start_points <- function(upper) {
    q <- length(upper)
    if (q == 0L) return(matrix(numeric(0), nrow = 1L, ncol = 0L))
    vals <- lapply(upper, function(u) unique(c(0, 0.5 * u, 0.9 * u)))
    as.matrix(expand.grid(vals, KEEP.OUT.ATTRS = FALSE,
                          stringsAsFactors = FALSE))
}

.joint_spectral_evaluator <- function(info, trace_coef, n, p_vec) {
    M <- length(info)
    masks <- .relaxation_masks(M)
    stopifnot(length(p_vec) == M, length(trace_coef) == nrow(masks), n >= 1L)
    active <- vapply(info, `[[`, logical(1), "active")
    mu <- lapply(info, `[[`, "mu")
    base_logdet <- vapply(info, `[[`, numeric(1), "logdet_A")
    indices <- lapply(seq_len(M), function(m) which(masks[, m] == 1L))
    p <- prod(p_vec)
    last_u <- last_result <- NULL
    function(u, need_gradient = TRUE) {
        stopifnot(length(u) == M, all(is.finite(u)), all(u[!active] == 0))
        if (!is.null(last_u) && identical(u, last_u)) {
            if (!need_gradient) return(list(loglik = last_result$loglik, gradient = NULL))
            if (!is.null(last_result$gradient)) return(last_result)
        }
        factors <- matrix(1, nrow(masks), M)
        prefix <- matrix(1, nrow(masks), M + 1L)
        for (m in seq_len(M)) {
            factors[indices[[m]], m] <- u[m]
            prefix[, m + 1L] <- prefix[, m] * factors[, m]
        }
        trace_value <- sum(trace_coef * prefix[, M + 1L])
        denominators <- Map(function(v, a) 1 + a * v, mu, u)
        if (!is.finite(trace_value) || trace_value <= 0 ||
                any(vapply(denominators, function(v) any(v <= 0), logical(1)))) {
            return(list(loglik = -Inf, gradient = numeric(M)))
        }
        logdet <- base_logdet + vapply(seq_len(M), function(m) {
            sum(log1p(u[m] * mu[[m]]))
        }, numeric(1))
        loglik <- n / 2 * (p * log(p / trace_value) +
                            sum((p / p_vec) * logdet) - p)
        if (!need_gradient) {
            last_u <<- u
            last_result <<- list(loglik = loglik, gradient = NULL)
            return(last_result)
        }
        suffix <- matrix(1, nrow(masks), M + 1L)
        for (m in rev(seq_len(M))) {
            suffix[, m] <- suffix[, m + 1L] * factors[, m]
        }
        trace_gradient <- vapply(seq_len(M), function(m) {
            rows <- indices[[m]]
            sum(trace_coef[rows] * prefix[rows, m] * suffix[rows, m + 1L])
        }, numeric(1))
        logdet_gradient <- vapply(seq_len(M), function(m) {
            sum(mu[[m]] / denominators[[m]])
        }, numeric(1))
        gradient <- n * p / 2 * (-trace_gradient / trace_value +
                                  logdet_gradient / p_vec)
        gradient[!active] <- 0
        last_u <<- u
        last_result <<- list(loglik = loglik, gradient = gradient)
        last_result
    }
}

.choose_joint_optimum <- function(optima, q_u) {
    valid <- vapply(optima, function(o) {
        !is.null(o) && is.finite(o$value) && length(o$par) == q_u &&
            identical(as.integer(o$convergence), 0L)
    }, logical(1))
    if (!any(valid)) stop("all joint optimizations failed to converge")
    converged <- optima[valid]
    converged[[which.min(vapply(converged, `[[`, numeric(1), "value"))]]
}

profile_joint_relaxation <- function(x_val, Omega_list, delta_list,
                            boundary_margin = 1e-8, u_cap = 10,
                            optim_control = list(maxit = 100L,
                                                 factr = 1e10)) {
    shape <- .validate_kronecker_inputs(x_val, Omega_list)
    if (length(delta_list) != shape$M) {
        stop("factor and relaxation-direction lists have different lengths")
    }
    for (m in seq_len(shape$M)) {
        if (!is.matrix(delta_list[[m]]) ||
                !identical(dim(delta_list[[m]]), dim(Omega_list[[m]]))) {
            stop(sprintf("mode %d relaxation direction has the wrong size", m))
        }
    }
    info <- Map(function(Om, De) relaxation_factor_info(
        Om, De, boundary_margin = boundary_margin, u_cap = u_cap),
        Omega_list, delta_list)
    active <- which(vapply(info, `[[`, logical(1), "active"))
    q_u <- as.integer(length(active))
    trace_coef <- .joint_trace_coefficients(x_val, info)
    masks <- .relaxation_masks(shape$M)
    u_zero <- numeric(shape$M)
    baseline <- .joint_loglik(
        u_zero, info, trace_coef, shape$n, shape$p_vec, masks)
    if (!is.finite(baseline$loglik)) {
        stop("u = 0 profiled likelihood is not finite")
    }

    best_u <- u_zero
    best <- baseline
    convergence <- 0L
    counts <- c("function" = 0L, "gradient" = 0L)

    if (q_u > 0L) {
        upper <- vapply(info[active], `[[`, numeric(1), "u_upper")
        evaluate <- .joint_spectral_evaluator(
            info, trace_coef, shape$n, shape$p_vec)
        objective <- function(z) {
            uu <- u_zero
            uu[active] <- z
            ll <- evaluate(uu, need_gradient = FALSE)$loglik
            if (is.finite(ll)) -ll else .Machine$double.xmax / 100
        }
        starts <- .joint_start_points(upper)
        optima <- vector("list", nrow(starts))
        for (i in seq_len(nrow(starts))) {
            optima[[i]] <- tryCatch(
                optim(starts[i, ], objective, method = "L-BFGS-B",
                      lower = rep(0, q_u), upper = upper,
                      control = optim_control),
                error = function(e) NULL)
        }
        chosen <- .choose_joint_optimum(optima, q_u)
        candidate_u <- u_zero
        candidate_u[active] <- pmin(pmax(chosen$par, 0), upper)
        candidate <- .joint_loglik(
            candidate_u, info, trace_coef, shape$n, shape$p_vec, masks)
        if (candidate$loglik > best$loglik) {
            best_u <- candidate_u
            best <- candidate
            convergence <- as.integer(chosen$convergence)
        }
        good <- optima[vapply(optima, function(o)
            !is.null(o) && identical(as.integer(o$convergence), 0L), logical(1))]
        counts <- colSums(do.call(rbind, lapply(good, function(o) o$counts)))
    }

    normalized <- Map(relaxation_transform_factor, info, best_u)
    norm_product <- prod(vapply(best$U, norm, numeric(1), type = "F"))
    list(
        loglik = as.numeric(best$loglik),
        loglik_u0 = as.numeric(baseline$loglik),
        ll_gain = as.numeric(best$loglik - baseline$loglik),
        u = as.numeric(best_u),
        q_u = q_u,
        tau_global = as.numeric(best$tau_raw * norm_product),
        Omega_list = normalized,
        convergence = convergence,
        counts = counts,
        u_pd = vapply(info, `[[`, numeric(1), "u_pd"),
        u_support = vapply(info, `[[`, numeric(1), "u_support"),
        u_upper = vapply(info, `[[`, numeric(1), "u_upper"),
        bound_kind = vapply(info, `[[`, character(1), "bound_kind"),
        cond_base = vapply(info, function(inf) {
            spd_condition_number(inf$A)
        }, numeric(1)),
        cond_joint = vapply(normalized, spd_condition_number, numeric(1))
    )
}

compute_joint_profile_ebic <- function(
        x_val, Omega_list, delta_list, gamma = 0.2, tol = 1e-6,
        boundary_margin = 1e-8, u_cap = 10,
        optim_control = list(maxit = 100L, factr = 1e10)) {
    stopifnot(length(gamma) == 1L, is.finite(gamma), gamma >= 0)
    prof <- profile_joint_relaxation(x_val, Omega_list, delta_list,
                            boundary_margin = boundary_margin, u_cap = u_cap,
                            optim_control = optim_control)
    M <- length(Omega_list)
    n <- dim(x_val)[M + 1L]
    p_total <- prod(dim(x_val)[seq_len(M)])
    df <- compute_df_tlasso(Omega_list, tol = tol)
    ebic <- -2 * prof$loglik + df * log(n) +
        4 * gamma * df * log(p_total)
    c(prof, list(
        ebic_primary = as.numeric(ebic),
        df = as.integer(df),
        p_total = p_total,
        gamma = gamma
    ))
}
