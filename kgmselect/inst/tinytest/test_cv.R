source("helper_reference.R")
golden <- readRDS("fixtures/golden.rds")
ref <- reference_for_installed_huge()
g <- golden$cv
cv <- kgm_cv(g$Y, g$eta, nfolds = 5L, seed = 2026L)
expect_identical(cv$fold_id, as.integer(g$fold_id))
expect_equal(cv$eta, g$eta_selected)
expect_equal(cv$index, which.max(colMeans(cv$fold_loglik)))
expect_identical(support_of(cv$factors), support_of(g$factors))

n <- dim(g$Y)[3]
set.seed(2026L, kind = "Mersenne-Twister", normal.kind = "Inversion",
         sample.kind = "Rejection")
expect_identical(cv$fold_id, sample(rep(seq_len(5L), length.out = n)))

set.seed(1L)
before <- get(".Random.seed", envir = globalenv())
invisible(kgm_cv(g$Y, g$eta[1:3], nfolds = 3L, seed = 7L))
expect_identical(get(".Random.seed", envir = globalenv()), before)

expect_error(kgm_cv(g$Y, g$eta, nfolds = 2L))
expect_error(kgm_cv(g$Y, g$eta, fold_id = rep(1:2, length.out = n)))

if (is.null(ref)) {
    exit_file(sprintf("no reference values for huge %s", packageVersion("huge")))
}
expect_equal(cv$fold_loglik, ref$cv$fold_loglik, tolerance = 1e-10)
expect_equal(cv$fold_tau, ref$cv$fold_tau, tolerance = 1e-10)
expect_equal(cv$tau, ref$cv$tau, tolerance = 1e-10)
expect_equal(cv$factors, ref$cv$factors, tolerance = 1e-10)
