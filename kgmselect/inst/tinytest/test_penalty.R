lam <- kgm_lambda(c(60L, 50L), n = 500, eta = 1)
expect_equal(lam, c(sqrt(log(60) / (500 * 50)), sqrt(log(50) / (500 * 60))),
             tolerance = 1e-14)
lam3 <- kgm_lambda(c(4L, 5L, 6L), n = 20, eta = 2)
expect_equal(lam3, 2 * sqrt(log(c(4, 5, 6)) / (20 * c(30, 24, 20))),
             tolerance = 1e-14)

g <- kgm_eta_grid()
expect_equal(length(g), 40L)
expect_equal(range(g), c(1e-4, 1e2), tolerance = 1e-12)
expect_identical(g, exp(seq(log(1e-4), log(1e2), length.out = 40L)))
expect_error(kgm_eta_grid(1L))
expect_error(kgm_eta_grid(5L, c(1, 0.5)))
