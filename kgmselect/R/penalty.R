.lambda_from_eta <- function(m_vec, n, eta) {
    p_minus_m <- prod(m_vec) / m_vec
    eta * sqrt(log(m_vec) / (n * p_minus_m))
}

kgm_eta_grid <- function(J = 40L, range = c(1e-4, 1e2)) {
    if (length(J) != 1L || !is.finite(J) || J != as.integer(J) || J < 2L) {
        stop("J must be an integer of at least 2")
    }
    if (length(range) != 2L || any(!is.finite(range)) || range[1] <= 0 ||
            range[2] <= range[1]) {
        stop("range must contain two increasing positive values")
    }
    exp(seq(log(range[1]), log(range[2]), length.out = J))
}

kgm_lambda <- function(p, n, eta) {
    .lambda_from_eta(p, n, eta)
}
