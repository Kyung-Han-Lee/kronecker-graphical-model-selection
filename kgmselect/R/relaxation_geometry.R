.relaxation_direction_varies <- function(A, delta) {
    f_delta <- norm(delta, type = "F")
    if (!(f_delta > 0)) return(FALSE)
    A_unit <- A / norm(A, type = "F")
    D_unit <- delta / f_delta
    tangent <- D_unit - sum(A_unit * D_unit) * A_unit

    norm(tangent, type = "F") > 8 * .Machine$double.eps
}

.projected_relaxation_direction <- function(Omega_now, Omega_prev, lambda_now, lambda_prev,
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

.spd_condition_number <- function(A) {
    stopifnot(is.matrix(A), nrow(A) == ncol(A))
    vals <- eigen((A + t(A)) / 2, symmetric = TRUE,
                  only.values = TRUE)$values
    if (min(vals) <= 0) stop("A must be positive definite")
    as.numeric(max(vals) / min(vals))
}

.relaxation_factor_info <- function(Omega, delta, boundary_margin = 1e-8,
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

.relaxation_transform_factor <- function(info, u) {
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

kgm_relaxation_direction <- function(Omega_now, Omega_prev, lambda_now,
                                     lambda_prev, tol = 1e-6, dir_tol = 1e-12) {
    .projected_relaxation_direction(Omega_now, Omega_prev, lambda_now,
                                    lambda_prev, tol = tol, dir_tol = dir_tol)
}

kgm_relaxation_interval <- function(Omega, delta, boundary_margin = 1e-8,
                                    u_cap = 10) {
    .relaxation_factor_info(Omega, delta, boundary_margin = boundary_margin,
                            u_cap = u_cap)
}

.relaxation_deltas <- function(path, k, tol) {
    fit <- path$fits[[k]]
    if (k == 1L) return(lapply(fit, function(A) A * 0))
    directions <- Map(.projected_relaxation_direction, fit,
                      path$fits[[k - 1L]], path$lambda[[k]],
                      path$lambda[[k - 1L]], MoreArgs = list(tol = tol))
    lapply(directions, `[[`, "delta")
}
