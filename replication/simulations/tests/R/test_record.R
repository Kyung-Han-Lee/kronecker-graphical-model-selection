source("helper_lib.R")
row <- data.frame(setting = "x", M = 2L, family = "er", dims = "12x10", d = 2,
                  m_new = NA, n = 7L, n_idx = 3L, reps = 2L, w = 0.3,
                  w_spread = 0.25, margin = 0.5, gradient = "numeric",
                  stringsAsFactors = FALSE)
config <- path_config("scale", 0.5, kgm_eta_grid(5L, c(0.1, 10)),
                      driver_code_hash(file.path("..", "..", "R")))
expect_identical(names(config), c("method", "eta", "fingerprint", "versions",
                                  "driver_code", "platform"))
expect_identical(config$method, list(name = "path", criteria = "scale", gammas = 0.5))
sig <- dataset_signature(row, 1L, config)
expect_identical(sig, dataset_signature(row, 1L, config))
expect_false(identical(sig, dataset_signature(row, 2L, config)))
expect_false(identical(sig, dataset_signature(transform(row, n_idx = 4L), 1L, config)))
expect_false(identical(sig, dataset_signature(transform(row, gradient = "analytic"), 1L, config)))
other_gammas <- modifyList(config, list(method = list(gammas = c(0, 0.5))))
expect_false(identical(sig, dataset_signature(row, 1L, other_gammas)))
expect_false(identical(sig, dataset_signature(row, 1L, modifyList(config, list(eta = rev(config$eta))))))
expect_false(identical(sig, dataset_signature(row, 1L, modifyList(config, list(fingerprint = "x")))))
expect_false(identical(sig, dataset_signature(row, 1L, modifyList(config, list(driver_code = "x")))))
expect_false(identical(sig, dataset_signature(row, 1L, modifyList(config, list(platform = "x")))))

cv <- driver_config(list(name = "cv", nfolds = 5L, fold_seed = 2026L), config$eta,
                    config$driver_code)
expect_false(identical(sig, dataset_signature(row, 1L, cv)))
other_seed <- modifyList(cv, list(method = list(fold_seed = 2027L)))
expect_false(identical(dataset_signature(row, 1L, cv), dataset_signature(row, 1L, other_seed)))
expect_equal(nchar(kgmselect_fingerprint()), 64L)
expect_identical(names(dependency_versions()),
                 c("kgmselect", "huge", "glasso", "expm", "rTensor", "igraph"))
expect_identical(names(platform_info()), c("platform", "blas", "lapack"))

lib <- tempfile()
dir.create(lib)
file.copy(list.files(file.path("..", "..", "R"), full.names = TRUE), lib)
expect_identical(driver_code_hash(lib), config$driver_code)
write("\n", file.path(lib, "scoring.R"), append = TRUE)
expect_false(identical(driver_code_hash(lib), config$driver_code))
expect_error(driver_code_hash(tempfile()), "no driver library")

dir <- tempfile()
dir.create(dir)
f <- file.path(dir, "a.rds")
save_atomic(list(signature = sig, x = 1), f)
expect_identical(list.files(dir, all.files = TRUE, no.. = TRUE), "a.rds")
expect_equal(reusable_result(f, sig)$x, 1)
expect_null(reusable_result(file.path(dir, "b.rds"), sig))
expect_error(reusable_result(f, "other"), "different setup")

expect_identical(parse_options(c("--workers=4", "--out=x"), list(out = "o", workers = "1")),
                 list(out = "x", workers = "4"))
expect_error(parse_options("--threads=4", list(workers = "1")), "unknown option")
