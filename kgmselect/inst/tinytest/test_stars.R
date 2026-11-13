golden <- readRDS("fixtures/golden.rds")
m2 <- golden$m2
st <- golden$stars
path <- kgm_path(m2$Y, m2$eta)
res <- kgm_stars(path, m2$Y, N = st$N, threshold = c(0.05, 0.1), seed = st$seed)
expect_equal(res$instability, st$D, tolerance = 1e-12)
expect_equal(res$index, c(st$index_005, st$index_010))
expect_equal(res$eta, m2$eta[res$index])
expect_equal(res$subsample_size, st$subsample_size)

expect_equal(kgmselect:::.stars_subsample_size(100L), 80L)
expect_equal(kgmselect:::.stars_subsample_size(200L), 141L)
expect_equal(kgmselect:::.stars_select(c(0.2, 0.08, 0.04, 0.06, 0.01), 1:5, 0.05), 5L)
expect_equal(kgmselect:::.stars_select(c(0.2, 0.08, 0.04, 0.03, 0.01), 1:5, 0.05), 3L)
