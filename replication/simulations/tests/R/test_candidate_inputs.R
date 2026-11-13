source("helper_lib.R")

row <- data.frame(setting = "input_test", M = 2L, family = "er", dims = "3x3",
                  d = 1, m_new = NA, n = 20L, n_idx = 1L, reps = 1L,
                  w = 0.3, w_spread = 0.25, margin = 0.5, gradient = "numeric",
                  stringsAsFactors = FALSE)
dataset <- data.frame(id = "input_test_r001", setting = "input_test", rep = 1L)
generated <- make_dataset(row, 1L)
eta <- kgm_eta_grid(4L, c(0.05, 2))
regular <- kgm_path(generated$Y, eta)
retained <- fit_path_with_inputs(dataset, generated$Y, eta)
expect_equal(retained$fits, regular$fits, tolerance = 0)
expect_equal(retained$lambda, regular$lambda, tolerance = 0)
expect_equal(nrow(attr(retained, "partial_covariances")), 16L)
expect_true(all(attr(retained, "partial_covariances")$minimum_eigenvalue > 0))
settings <- tempfile(fileext = ".csv")
write.csv(row, settings, row.names = FALSE, na = "")
output <- tempfile()
run_experiment("matrix", settings, output, c("ebic", "scale", "relaxation"), 0.5,
               driver_code_hash(file.path("..", "..", "R")),
               eta = eta, save_inputs = TRUE)
cached <- readRDS(file.path(output, "datasets", "input_test_r001.rds"))
expect_equal(cached$inputs$path$fits, regular$fits, tolerance = 0)
expect_identical(digest::digest(make_dataset(cached$inputs$settings, cached$rep)$Y, algo = "sha256"),
                 cached$data_sha256)
