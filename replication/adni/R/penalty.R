lambda_from_eta <- function(m_vec, n, eta) {
    p_minus_m <- prod(m_vec) / m_vec
    eta * sqrt(log(m_vec) / (n * p_minus_m))
}
