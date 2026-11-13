.finalize_precision <- function(A, margin = 0.5, diag_mode = "eigen",
                                diag_const = NA_real_) {
    stopifnot(is.matrix(A), nrow(A) == ncol(A), isSymmetric(A),
              length(diag_mode) == 1L, is.character(diag_mode),
              diag_mode %in% c("rowsum", "eigen", "fixed"))
    diag(A) <- 0

    if (diag_mode == "fixed") {
        if (!is.na(margin)) {
            stop("margin is not used when diag_mode = \"fixed\"; ",
                 "the margin follows from diag_const and the graph. ",
                 "Pass margin = NA_real_ when calling the generator directly")
        }
        stopifnot(length(diag_const) == 1L, is.finite(diag_const),
                  diag_const > 0)
        lm_A <- abs(min(eigen(A, symmetric = TRUE, only.values = TRUE)$values))
        if (diag_const <= lm_A) {
            stop(sprintf(paste0("diag_const %.6g does not exceed |lambda_min(A)| %.6g, ",
                                "so the factor is not positive definite (lambda_min(Omega) = %.6g). ",
                                "Choose diag_const larger than %.6g"),
                         diag_const, lm_A, diag_const - lm_A, lm_A))
        }
        diag(A) <- diag_const
        return(A)
    }

    if (!is.na(diag_const)) {
        stop("diag_const is used only when diag_mode = \"fixed\"")
    }
    stopifnot(length(margin) == 1L, is.finite(margin), margin > 0)
    diag(A) <- switch(diag_mode,
        rowsum = rowSums(abs(A)) + margin,
        eigen = abs(min(eigen(A, symmetric = TRUE,
                              only.values = TRUE)$values)) + margin)
    A
}

.make_er_precision <- function(p, rho, w = 0.3, seed = 1L, w_spread = 0,
                              margin = NULL, diag_mode = "eigen",
                              diag_const = NA_real_) {
    if (is.null(margin)) {
        margin <- if (diag_mode == "fixed") NA_real_ else 0.5
    }
    stopifnot(p >= 2L, rho > 0, rho <= 1, w > 0,
              !is.null(seed), length(seed) == 1L, is.finite(seed),
              length(w_spread) == 1L, is.finite(w_spread),
              w_spread >= 0, w_spread < w)
    set.seed(seed)
    iu <- upper.tri(diag(p))
    A <- matrix(0, p, p)
    edges <- rbinom(sum(iu), 1L, rho)
    signs <- sample(c(-1, 1), length(edges), replace = TRUE)
    mags <- if (w_spread > 0) {
        runif(length(edges), w - w_spread, w + w_spread)
    } else {
        w
    }
    A[iu] <- edges * signs * mags
    A <- A + t(A)
    .finalize_precision(A, margin = margin, diag_mode = diag_mode,
                        diag_const = diag_const)
}

.make_ba_precision <- function(p, m_new = 1L, w = 0.3, seed = 1L, w_spread = 0,
                              margin = NULL, diag_mode = "eigen",
                              diag_const = NA_real_) {
    if (is.null(margin)) {
        margin <- if (diag_mode == "fixed") NA_real_ else 0.5
    }
    stopifnot(p >= 2L, m_new >= 1L, w > 0,
              !is.null(seed), length(seed) == 1L, is.finite(seed),
              length(w_spread) == 1L, is.finite(w_spread),
              w_spread >= 0, w_spread < w)
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("the igraph package is needed to generate BA graphs")
    }
    set.seed(seed)
    g <- igraph::sample_pa(p, m = m_new, directed = FALSE)
    A <- as.matrix(igraph::as_adjacency_matrix(g))
    A[A != 0] <- 1
    s <- matrix(sample(c(-1, 1), p * p, replace = TRUE), p, p)
    s[lower.tri(s)] <- t(s)[lower.tri(s)]
    mg <- if (w_spread > 0) {
        m <- matrix(runif(p * p, w - w_spread, w + w_spread), p, p)
        m[lower.tri(m)] <- t(m)[lower.tri(m)]
        m
    } else {
        w
    }
    A <- A * s * mg
    .finalize_precision(A, margin = margin, diag_mode = diag_mode,
                        diag_const = diag_const)
}

.sym_sqrt <- function(M) {
    e <- eigen(M, symmetric = TRUE)
    if (min(e$values) < -1e-10) stop("M is not positive semidefinite")
    d <- sqrt(pmax(e$values, 0))
    e$vectors %*% (d * t(e$vectors))
}

.mode_products <- function(Z, Sh, p) {
    K <- length(p)
    Y <- Z
    dim(Y) <- c(p[1], prod(p) / p[1])
    Y <- Sh[[1]] %*% Y
    if (K >= 3L) {
        for (m in 2:(K - 1L)) {
            dim(Y) <- p
            perm <- c(m, seq_len(K)[-m])
            Y <- aperm(Y, perm)
            dim(Y) <- c(p[m], prod(p) / p[m])
            Y <- Sh[[m]] %*% Y
            dim(Y) <- p[perm]
            Y <- aperm(Y, order(perm))
        }
    }
    dim(Y) <- c(prod(p) / p[K], p[K])
    Y <- Y %*% Sh[[K]]
    dim(Y) <- p
    Y
}

.sample_tensor_normal <- function(n, ..., seed = NULL) {
    sigmas <- list(...)
    if (length(sigmas) == 1L && is.list(sigmas[[1]])) sigmas <- sigmas[[1]]
    K <- length(sigmas)
    if (K < 2L) stop("at least two modes are required (got ", K, ")")
    for (k in seq_len(K)) {
        Sk <- sigmas[[k]]
        if (!is.matrix(Sk) || nrow(Sk) != ncol(Sk)) {
            stop("Sigma ", k, " is not a square matrix")
        }
    }
    p <- vapply(sigmas, nrow, integer(1))
    Sh <- lapply(sigmas, .sym_sqrt)

    if (!is.null(seed)) set.seed(seed)
    p_cell <- prod(p)
    out <- array(0, dim = c(p, n))
    for (i in seq_len(n)) {
        Z <- array(rnorm(p_cell), dim = p)
        out[(i - 1L) * p_cell + seq_len(p_cell)] <- .mode_products(Z, Sh, p)
    }
    out
}

kgm_simulate_precision <- function(p, family = c("er", "ba"), degree = 2,
                                   m_new = 1L, w = 0.3, w_spread = 0.25,
                                   min_eigen = 0.5, seed) {
    family <- match.arg(family)
    .with_seed(seed, if (family == "er") {
        .make_er_precision(p, rho = degree / (p - 1), w = w, seed = seed,
                           w_spread = w_spread, margin = min_eigen)
    } else {
        .make_ba_precision(p, m_new = m_new, w = w, seed = seed,
                           w_spread = w_spread, margin = min_eigen)
    })
}

kgm_sample <- function(n, Sigma, seed = NULL) {
    if (!is.list(Sigma)) stop("Sigma must be a list of covariance factors")
    if (is.null(seed)) return(.sample_tensor_normal(n, Sigma))
    .with_seed(seed, .sample_tensor_normal(n, Sigma, seed = seed))
}
