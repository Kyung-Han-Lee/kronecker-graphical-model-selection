source("helper_lib.R")
s <- read_settings(file.path("..", "..", "settings", "matrix.csv"))
expect_equal(nrow(s), 140L)
expect_equal(sort(unique(s$n)), c(10L, 30L, 50L, 100L, 200L, 800L, 3200L))
expect_true(all(tapply(s$n_idx, s$n, function(x) length(unique(x))) == 1L))
expect_true(all(s$gradient == "numeric"))
t <- read_settings(file.path("..", "..", "settings", "tensor.csv"))
expect_equal(nrow(t), 24L)
expect_equal(as.vector(table(t$M)), c(12L, 12L))
expect_identical(t$gradient == "numeric", t$n == 10L)

d <- expand_datasets(s)
expect_equal(nrow(d), 14000L)
expect_identical(d$id[1:2], c("m2_er_d1_p40x40_n0100_r001", "m2_er_d1_p40x40_n0100_r002"))
expect_equal(nrow(expand_datasets(t, reps = 2L)), 48L)
expect_identical(dims_from("8x8x10x16"), c(8L, 8L, 10L, 16L))

f <- tempfile(fileext = ".csv")
bad <- s[1:2, ]
bad$family[1] <- "sbm"
write.csv(bad, f, row.names = FALSE, na = "")
expect_error(read_settings(f), "invalid rows")
write.csv(s[c(1, 1), ], f, row.names = FALSE, na = "")
expect_error(read_settings(f), "repeats")
write.csv(s[1:2, names(s) != "n_idx"], f, row.names = FALSE, na = "")
expect_error(read_settings(f), "lacks columns: n_idx")

r <- read_settings(file.path("..", "..", "settings", "refit_comparison.csv"))
expect_equal(nrow(r), 84L)
expect_equal(as.vector(table(r$M)), c(60L, 12L, 12L))
expect_identical(sort(unique(r$n[r$M == 2L])), c(100L, 200L, 800L, 3200L))
expect_true(all(r$n_idx[r$M == 4L & r$n == 5L] == 3L))
expect_equal(nrow(expand_datasets(r)), 8400L)
p <- read_settings(file.path("..", "..", "settings", "refit_pilot.csv"))
expect_identical(p$n, c(2L, 5L, 10L, 20L, 40L))
expect_identical(p$n_idx, 1:5)
expect_equal(nrow(expand_datasets(p)), 250L)

for (name in c("matrix", "tensor", "refit_comparison", "refit_pilot")) {
    q <- read_settings(file.path("..", "..", "settings", "quick", paste0(name, ".csv")))
    expect_true(all(q$reps == 2L), info = name)
}
qs <- read_stars_seeds(file.path("..", "..", "settings", "quick", "stars_seeds.csv"))
expect_identical(qs$seed, stars_seed(qs$source))
expect_identical(qs$dataset, expand_datasets(read_settings(
    file.path("..", "..", "settings", "quick", "matrix.csv")))$id)
