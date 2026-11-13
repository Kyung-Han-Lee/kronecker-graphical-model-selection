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
