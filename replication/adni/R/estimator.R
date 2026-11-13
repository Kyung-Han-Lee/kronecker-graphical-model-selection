.partial_cov <- function(data, sqrt_list, m) {

    dims <- dim(data)
    n_mode <- length(dims) - 1L
    p_vec <- dims[1:n_mode]
    n_obs <- dims[n_mode + 1L]

    sqrt_list[[m]] <- diag(p_vec[m])

    denom <- n_obs * prod(p_vec) / p_vec[m]

    block_len <- prod(p_vec)
    acc <- matrix(0, p_vec[m], p_vec[m])

    if (n_mode == 2L) {
        whiten_gram <- crossprod(sqrt_list[[if (m == 1L) 2L else 1L]])
        for (i in seq_len(n_obs)) {
            obs <- array(data[((i - 1L) * block_len + 1L):(i * block_len)],
                         dim = p_vec)
            acc <- acc + if (m == 1L) {
                obs %*% whiten_gram %*% t(obs)
            } else {
                t(obs) %*% whiten_gram %*% obs
            }
        }
    } else {
        for (i in seq_len(n_obs)) {
            obs <- array(data[((i - 1L) * block_len + 1L):(i * block_len)],
                         dim = p_vec)
            whitened <- rTensor::ttl(rTensor::as.tensor(obs), sqrt_list,
                                     ms = seq_len(n_mode))
            acc <- acc + tcrossprod(rTensor::k_unfold(whitened, m = m)@data)
        }
    }

    acc <- (acc + t(acc)) / 2
    acc / denom
}

fit_tggm <- function(data, lambda, n_sweep = 2L, return_S = FALSE) {

    dims <- dim(data)
    n_mode <- length(dims) - 1L
    p_vec <- dims[1:n_mode]

    stopifnot(n_mode >= 1L, length(lambda) == n_mode,
              all(is.finite(lambda)), all(lambda >= 0),
              length(n_sweep) == 1L, n_sweep >= 1L,
              n_sweep == as.integer(n_sweep),
              length(return_S) == 1L, is.logical(return_S), !is.na(return_S))

    Omega <- lapply(p_vec, function(p) diag(p))
    Omega_sqrt <- lapply(p_vec, function(p) diag(p))

    S_by_sweep <- if (return_S) vector("list", n_sweep) else NULL

    for (sweep_idx in seq_len(n_sweep)) {
        for (m in seq_len(n_mode)) {
            S_mode <- .partial_cov(data, Omega_sqrt, m)
            if (return_S) {
                if (is.null(S_by_sweep[[sweep_idx]])) {
                    S_by_sweep[[sweep_idx]] <- vector("list", n_mode)
                }
                S_by_sweep[[sweep_idx]][[m]] <- S_mode
            }
            fit <- huge::huge(S_mode, lambda = lambda[m], method = "glasso",
                              verbose = FALSE)
            Omega_raw <- as.matrix(fit$icov[[1]])
            frob <- norm(Omega_raw, type = "F")
            if (frob < .Machine$double.eps) stop("적합이 0 행렬이다")
            Omega[[m]] <- Omega_raw / frob

            root <- expm::sqrtm(Omega[[m]])
            if (is.complex(root)) {
                if (max(abs(Im(root))) > 1e-8) {
                    stop("Omega^{1/2} 가 복소다. 모드 ", m)
                }
                root <- Re(root)
            }
            Omega_sqrt[[m]] <- root
        }
    }

    Omega <- lapply(Omega, function(M) {
        M <- (M + t(M)) / 2
        M / norm(M, type = "F")
    })
    if (return_S) attr(Omega, "S_by_sweep") <- S_by_sweep
    Omega
}
