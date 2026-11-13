source("helper_lib.R")

set.seed(2026L, kind = "Mersenne-Twister", normal.kind = "Inversion",
         sample.kind = "Rejection")
expect_identical(cv_fold_id(23L), sample(rep(seq_len(5L), length.out = 23L)))
set.seed(1L)
before <- get(".Random.seed", envir = globalenv())
invisible(cv_fold_id(10L))
expect_identical(get(".Random.seed", envir = globalenv()), before)

expect_identical(cv_select_index(rbind(c(1, 3, 3), c(2, 1, 1))), 2L)

tiny <- data.frame(setting = c("t_er", "t_ba"), M = 2L, family = c("er", "ba"),
                   dims = "6x5", d = c(1, NA), m_new = c(NA, 1L), n = 20L,
                   n_idx = 1L, reps = 2L, w = 0.3, w_spread = 0.25, margin = 0.5,
                   gradient = "numeric", stringsAsFactors = FALSE)
f <- tempfile(fileext = ".csv")
write.csv(tiny, f, row.names = FALSE, na = "")
eta <- kgm_eta_grid(4L, c(0.05, 5))
gen <- make_dataset(read_settings(f)[1, ], 1L)
cv <- kgm_cv(gen$Y, eta, nfolds = 5L, seed = 2026L)
expect_identical(cv$fold_id, cv_fold_id(20L))
expect_identical(cv_select_index(cv$fold_loglik), cv$index)

Y2 <- gen$Y
Y2[, , cv$fold_id == 1L] <- 2 * Y2[, , cv$fold_id == 1L]
cv2 <- kgm_cv(Y2, eta, fold_id = cv$fold_id)
expect_identical(cv2$fold_tau[1, ], cv$fold_tau[1, ])
expect_true(all(cv2$fold_loglik[1, ] != cv$fold_loglik[1, ]))

out <- tempfile()
code <- driver_code_hash(file.path("..", "..", "R"))
run_cv_experiment("tiny", f, out, code, workers = 2L, reps = 1L, eta = eta)
r <- readRDS(file.path(out, "datasets", "t_er_r001.rds"))
expect_identical(r$fold_id, cv$fold_id)
expect_identical(r$fold_loglik, cv$fold_loglik)
expect_identical(r$fold_tau, cv$fold_tau)
expect_identical(r$tau, cv$tau)
expect_identical(r$selected$k, cv$index)
expect_identical(r$selected$eta, eta[cv$index])
m <- kgm_metrics(gen$truth, cv$factors)
expect_identical(r$selected$mcc, m$mcc)
expect_identical(r$selected_modes$tp, m$per_mode$tp)
expect_identical(r$data_sha256, digest::digest(gen$Y, algo = "sha256"))
expect_identical(jsonlite::read_json(file.path(out, "run.json"))$method,
                 list(name = "cv", nfolds = 5L, fold_seed = 2026L))
run_cv_experiment("tiny", f, out, code, workers = 2L, reps = 1L, eta = eta)
expect_equal(jsonlite::read_json(file.path(out, "run.json"))$datasets$reused, 2L)

other_folds <- function(dataset, Y, eta, fold_id) {
    x <- kgm_cv(Y, eta, fold_id = fold_id)
    x$fold_id <- rev(x$fold_id)
    x
}
out2 <- tempfile()
expect_error(run_cv_experiment("tiny", f, out2, code, ids = "t_er_r001", eta = eta,
                               cv_fun = other_folds), "1 datasets failed")
expect_true(grepl("other folds", readLines(file.path(out2, "errors", "t_er_r001.txt"))))
not_max <- function(dataset, Y, eta, fold_id) {
    x <- kgm_cv(Y, eta, fold_id = fold_id)
    x$index <- if (x$index == 1L) 2L else 1L
    x
}
out3 <- tempfile()
expect_error(run_cv_experiment("tiny", f, out3, code, ids = "t_er_r001", eta = eta,
                               cv_fun = not_max), "1 datasets failed")
expect_true(grepl("does not maximise", readLines(file.path(out3, "errors", "t_er_r001.txt"))))
