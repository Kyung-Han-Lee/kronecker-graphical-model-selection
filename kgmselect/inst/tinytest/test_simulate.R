golden <- readRDS("fixtures/golden.rds")
s <- golden$sim

er <- kgm_simulate_precision(12L, "er", degree = 2, seed = 1001L)
expect_equal(er, s$er, tolerance = 1e-12)
ev <- eigen(er, symmetric = TRUE, only.values = TRUE)$values
expect_equal(min(ev), 0.5, tolerance = 1e-10)
off <- er[upper.tri(er)]
expect_true(all(abs(off[off != 0]) >= 0.05 & abs(off[off != 0]) <= 0.55))

if (requireNamespace("igraph", quietly = TRUE)) {
    ba <- kgm_simulate_precision(12L, "ba", m_new = 1L, seed = 2001L)
    expect_equal(sum(ba[upper.tri(ba)] != 0), 11L)

    if (identical(as.character(packageVersion("igraph")), s$igraph)) {
        expect_equal(ba, s$ba, tolerance = 1e-12)
    }
}

Y <- kgm_sample(7L, list(solve(s$er), solve(s$ba)), seed = 5011L)
expect_equal(Y, s$sample, tolerance = 1e-12)
expect_equal(dim(Y), c(12L, 12L, 7L))
expect_error(kgm_sample(5L, solve(s$er)))

set.seed(3L)
before <- get(".Random.seed", envir = globalenv())
invisible(kgm_simulate_precision(8L, seed = 4L))
invisible(kgm_sample(3L, list(diag(2), diag(3)), seed = 5L))
expect_identical(get(".Random.seed", envir = globalenv()), before)
