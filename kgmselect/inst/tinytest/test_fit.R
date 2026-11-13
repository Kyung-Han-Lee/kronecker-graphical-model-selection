source("helper_reference.R")
golden <- readRDS("fixtures/golden.rds")
ref <- reference_for_installed_huge()
m2 <- golden$m2
m3 <- golden$m3

fit4 <- kgm_fit(m2$Y, m2$lambda[[4]])
expect_equal(vapply(fit4, function(A) sqrt(sum(A^2)), 0), c(1, 1),
             tolerance = 1e-12)

path <- kgm_path(m2$Y, m2$eta)
expect_inherits(path, "kgm_path")
expect_identical(path$dims, c(6L, 5L))
expect_identical(path$n, 30L)
expect_equal(path$lambda, m2$lambda, tolerance = 1e-15)
expect_identical(lapply(path$fits, support_of), lapply(m2$fits, support_of))
path3 <- kgm_path(m3$Y, m3$eta)
expect_identical(lapply(path3$fits, support_of), lapply(m3$fits, support_of))

expect_error(kgm_path(m2$Y, rev(m2$eta)))
expect_error(kgm_fit(matrix(1, 3, 3), c(0.1, 0.1)))

if (is.null(ref)) {
    exit_file(sprintf("no reference values for huge %s", packageVersion("huge")))
}
expect_equal(fit4, ref$m2$fits[[4]], tolerance = 1e-10)
expect_equal(path$fits, ref$m2$fits, tolerance = 1e-10)
expect_equal(path3$fits, ref$m3$fits, tolerance = 1e-10)
