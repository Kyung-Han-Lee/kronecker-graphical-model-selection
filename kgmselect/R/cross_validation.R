kgm_cv <- function(Y, eta = kgm_eta_grid(), nfolds = 5L, seed = 2026L,
                   n_sweep = 2L, fold_id = NULL) {
    .check_data(Y)
    .check_eta(eta)
    dims <- dim(Y)
    M <- length(dims) - 1L
    p <- dims[seq_len(M)]
    n <- dims[M + 1L]
    if (is.null(fold_id)) {
        if (length(nfolds) != 1L || !is.finite(nfolds) ||
                nfolds != as.integer(nfolds) || nfolds < 3L ||
                nfolds > floor(n / 2)) {
            stop("nfolds must be an integer between 3 and n / 2")
        }
        fold_id <- .with_seed(seed, sample(rep(seq_len(nfolds), length.out = n)))
    } else {
        if (length(fold_id) != n || any(!is.finite(fold_id)) ||
                any(fold_id != as.integer(fold_id)) || any(fold_id < 1L)) {
            stop("fold_id must hold one positive integer label per observation")
        }
        fold_id <- as.integer(fold_id)
        nfolds <- max(fold_id)
        if (nfolds < 3L || !identical(sort(unique(fold_id)), seq_len(nfolds)) ||
                any(tabulate(fold_id, nbins = nfolds) < 2L)) {
            stop("fold_id needs at least three folds of two or more observations")
        }
    }
    J <- length(eta)
    fold_loglik <- matrix(NA_real_, nfolds, J)
    fold_tau <- matrix(NA_real_, nfolds, J)
    for (i in seq_len(nfolds)) {
        held <- fold_id == i
        train <- .slice_samples(Y, which(!held))
        valid <- .slice_samples(Y, which(held))
        n_train <- sum(!held)
        for (j in seq_len(J)) {
            fit <- .fit_tggm(train, lambda = kgm_lambda(p, n_train, eta[j]),
                             n_sweep = n_sweep)
            tau <- attr(.compute_loglik(train, fit, tau_correct = TRUE,
                                        sqrt_method = "none"), "tau_global")
            if (!is.finite(tau) || tau <= 0) stop("invalid training-fold scale")

            scaled <- lapply(fit, function(A) A * tau^(1 / M))
            value <- as.numeric(.compute_loglik(valid, scaled,
                                                tau_correct = FALSE,
                                                sqrt_method = "none"))
            if (!is.finite(value)) stop("non-finite validation log-likelihood")
            fold_loglik[i, j] <- value
            fold_tau[i, j] <- tau
        }
    }
    loglik <- colMeans(fold_loglik)
    index <- which.max(loglik)
    lambda <- kgm_lambda(p, n, eta[index])
    factors <- lapply(.fit_tggm(Y, lambda = lambda, n_sweep = n_sweep), identity)
    tau <- attr(.compute_loglik(Y, factors, tau_correct = TRUE,
                                sqrt_method = "none"), "tau_global")
    structure(list(eta = eta[index], index = index, lambda = lambda,
                   factors = factors, tau = tau, fold_id = fold_id,
                   fold_loglik = fold_loglik, fold_tau = fold_tau,
                   loglik = loglik, eta_grid = eta),
              class = "kgm_cv")
}
