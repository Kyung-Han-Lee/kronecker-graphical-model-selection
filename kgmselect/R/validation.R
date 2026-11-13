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

.check_data <- function(Y) {
    dims <- dim(Y)
    if (!is.numeric(Y) || length(dims) < 3L || any(dims < 1L) ||
            any(!is.finite(Y))) {
        stop("Y must be a finite numeric array p_1 x ... x p_M x n with M >= 2")
    }
    invisible(Y)
}

.check_eta <- function(eta) {
    if (!is.numeric(eta) || !length(eta) || any(!is.finite(eta)) ||
            any(eta <= 0) || any(diff(eta) <= 0)) {
        stop("eta must be a positive, strictly increasing numeric vector")
    }
    invisible(eta)
}

.check_path <- function(path, Y) {
    if (!inherits(path, "kgm_path")) {
        stop("path must be a kgm_path object returned by kgm_path()")
    }
    dims <- dim(Y)
    M <- length(dims) - 1L
    if (!identical(as.integer(dims[seq_len(M)]), path$dims) ||
            !identical(as.integer(dims[M + 1L]), path$n)) {
        stop("Y does not have the dimensions of the data used to fit the path")
    }
    invisible(path)
}
