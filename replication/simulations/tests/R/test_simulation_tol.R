source("helper_lib.R")

A <- diag(2)
A[1, 2] <- A[2, 1] <- 5e-7
A <- A / sqrt(sum(A^2))
B <- diag(2) / sqrt(2)
truth <- list(matrix(c(1, 0.2, 0.2, 1), 2), diag(2))
set.seed(31)
Y <- array(rnorm(32), c(2, 2, 8))
path <- structure(list(fits = list(list(A, B)), eta = 0.1,
                       lambda = list(c(0.1, 0.1)), dims = c(2L, 2L),
                       n = 8L, n_sweep = 2L), class = "kgm_path")

scores <- score_path(path, Y, c("ebic", "scale", "relaxation"), 0.5, "numeric")
expect_identical(scores$edges, rep(1L, 3))
tables <- selection_tables(scores, path, truth)
expect_true(all(tables$selected$tp == 1))
expect_true(all(tables$selected$fn == 0))
expect_true(all(tables$selected_modes$tp[tables$selected_modes$mode == 1] == 1))
chosen <- fit_tables(data.frame(criterion = "cv"), 1L, 0.1, path$fits[[1]], truth)
expect_equal(chosen$selected$tp, 1)
expect_equal(candidate_recovery(path, truth)$tp, 1)
supports <- support_records(path)
expect_equal(supports$edges, 1)
expect_equal(supports$tiny_edges, 0)
expect_identical(supports$pattern_hash, supports$exact_hash)

expected <- kgm_score(path, Y, "refit", gamma = 0.5, tol = 0)
refitted <- refit_within(path, Y, 1L, 0.5, 30)
expect_identical(refitted$status, expected$status)
expect_equal(refitted$loglik, expected$loglik, tolerance = 1e-12)
