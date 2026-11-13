.stars_subsample_size <- function(n) {
    stopifnot(length(n) == 1L, n >= 2L)
    ratio <- if (n <= 144) 0.8 else 10 * sqrt(n) / n
    as.integer(floor(n * ratio))
}

.stars_subsamples <- function(n, N, seed) {
    stopifnot(length(N) == 1L, N >= 2L)
    old_kind <- RNGkind()
    old_seed <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) .Random.seed else NULL
    on.exit({
        do.call(RNGkind, as.list(old_kind))
        if (is.null(old_seed)) {
            if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
        } else assign(".Random.seed", old_seed, .GlobalEnv)
    }, add = TRUE)
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    b <- .stars_subsample_size(n)
    lapply(seq_len(N), function(s) sample.int(n, b, replace = FALSE))
}

.edge_indicators <- function(fit, tol) {
    unlist(lapply(fit, function(A) abs(A[upper.tri(A)]) > tol), use.names = FALSE)
}

.instability_from_theta <- function(theta) mean(2 * theta * (1 - theta))

.stars_instability <- function(E) .instability_from_theta(rowMeans(E))

.stars_select <- function(D, eta_grid, beta) {
    stopifnot(length(D) == length(eta_grid), all(diff(eta_grid) > 0), all(is.finite(D)))
    J <- length(eta_grid)
    from_sparse <- rev(seq_len(J))
    eligible <- which(cummax(D[from_sparse]) <= beta)
    if (!length(eligible)) {
        return(J)
    }
    from_sparse[max(eligible)]
}

.slice_samples <- function(data, idx) {
    M <- length(dim(data)) - 1L
    do.call(`[`, c(list(data), rep(list(TRUE), M), list(idx), list(drop = FALSE)))
}

.stars_path_instability <- function(data, lambda_grid, N = 20L, seed, n_sweep = 2L, tol = 1e-6) {
    stopifnot(is.list(lambda_grid), length(lambda_grid) >= 1L)
    dims <- dim(data)
    M <- length(dims) - 1L
    n <- dims[M + 1L]
    subs <- .stars_subsamples(n, N, seed)
    J <- length(lambda_grid)
    counts <- NULL
    for (s in seq_len(N)) {
        x <- .slice_samples(data, subs[[s]])
        for (j in seq_len(J)) {
            fit <- .fit_tggm(x, lambda = lambda_grid[[j]], n_sweep = n_sweep)
            e <- .edge_indicators(fit, tol)
            if (is.null(counts)) counts <- matrix(0L, length(e), J)
            counts[, j] <- counts[, j] + e
        }
    }
    apply(counts / N, 2, .instability_from_theta)
}

kgm_stars <- function(path, Y, N = 20L, threshold = 0.05, seed, tol = 1e-6) {
    .check_data(Y)
    .check_path(path, Y)
    D <- .stars_path_instability(Y, path$lambda, N = N, seed = seed,
                                 n_sweep = path$n_sweep, tol = tol)
    index <- vapply(threshold, function(b) .stars_select(D, path$eta, b),
                    integer(1))
    list(index = index, eta = path$eta[index], threshold = threshold,
         instability = D, subsample_size = .stars_subsample_size(path$n))
}
