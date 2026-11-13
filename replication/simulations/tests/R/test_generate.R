source("helper_lib.R")
row <- data.frame(setting = "x", M = 2L, family = "er", dims = "12x10", d = 2,
                  m_new = NA, n = 7L, n_idx = 3L, reps = 2L, w = 0.3,
                  w_spread = 0.25, margin = 0.5, gradient = "numeric",
                  stringsAsFactors = FALSE)
seeds <- dataset_seeds(row, 4L)
expect_identical(seeds$graph, c(1004L, 2004L))
expect_identical(seeds$noise, 5043L)

a <- make_dataset(row, 4L)
expect_identical(dim(a$Y), c(12L, 10L, 7L))
expect_identical(a$truth[[1]],
                 kgm_simulate_precision(12L, "er", degree = 2, w = 0.3,
                                        w_spread = 0.25, min_eigen = 0.5,
                                        seed = 1004L))
expect_identical(a$Y, kgm_sample(7L, lapply(a$truth, solve), seed = 5043L))

b <- make_dataset(transform(row, n = 9L, n_idx = 4L), 4L)
expect_identical(a$truth, b$truth)
expect_identical(dim(b$Y), c(12L, 10L, 9L))

ba <- make_dataset(transform(row, family = "ba", d = NA, m_new = 1L, M = 3L,
                             dims = "8x6x5"), 1L)
expect_equal(vapply(ba$truth, function(A) sum(A[upper.tri(A)] != 0), 0),
             c(7, 5, 4))
