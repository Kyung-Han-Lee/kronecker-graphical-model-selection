expect_equal(kgmselect:::.mcc(10, 0, 0, 10), 1)
expect_equal(kgmselect:::.mcc(0, 10, 10, 0), -1)
expect_equal(kgmselect:::.mcc(0, 0, 30, 100), 0)

truth <- list(matrix(c(1, 0.3, 0, 0.3, 1, 0, 0, 0, 1), 3), diag(2))
fit <- list(matrix(c(1, 0.2, 0.1, 0.2, 1, 0, 0.1, 0, 1), 3), diag(2))
r <- kgm_metrics(truth, fit)
expect_equal(c(r$tp, r$fp, r$fn, r$tn), c(1, 1, 0, 2))
expect_equal(r$fdp, 0.5)
expect_equal(r$tpr, 1)
expect_equal(r$mcc, 2 / sqrt(12), tolerance = 1e-12)
expect_equal(r$per_mode$tn, c(1, 1))

empty <- kgm_metrics(truth, list(diag(3), diag(2)))
expect_equal(c(empty$fdp, empty$tpr, empty$mcc), c(0, 0, 0))

golden <- readRDS("fixtures/golden.rds")
g <- golden$metrics
r2 <- kgm_metrics(g$truth, g$fit)
expect_equal(r2$per_mode$tp, c(g$counts$tp1, g$counts$tp2))
expect_equal(r2$per_mode$fp, c(g$counts$fp1, g$counts$fp2))
