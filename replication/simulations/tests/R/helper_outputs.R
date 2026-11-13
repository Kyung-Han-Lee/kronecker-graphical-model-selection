tiny <- data.frame(setting = c("t_er", "t_ba"), M = 2L, family = c("er", "ba"),
                   dims = "8x6", d = c(1, NA), m_new = c(NA, 1L), n = 20L,
                   n_idx = 1L, reps = 2L, w = 0.3, w_spread = 0.25, margin = 0.5,
                   gradient = "numeric", stringsAsFactors = FALSE)
settings_dir <- tempfile()
dir.create(settings_dir)
settings_file <- file.path(settings_dir, "matrix.csv")
write.csv(tiny, settings_file, row.names = FALSE, na = "")
ids <- expand_datasets(read_settings(settings_file))$id
write.csv(data.frame(setting = "t_er", gamma = 0.5), file.path(settings_dir, "figure1.csv"),
          row.names = FALSE)
seeds <- file.path(settings_dir, "stars_seeds.csv")
write.csv(data.frame(dataset = ids, source = ids, seed = seq_along(ids)), seeds,
          row.names = FALSE)
root <- tempfile()
eta <- kgm_eta_grid(5L, c(0.01, 10))
code <- driver_code_hash(file.path("..", "..", "R"))
run_experiment("matrix", settings_file, file.path(root, "matrix"),
               c("ebic", "scale", "relaxation"), c(0.5, 1), code, workers = 2L, eta = eta)
run_cv_experiment("matrix", settings_file, file.path(root, "cv", "matrix"), code,
                  workers = 2L, eta = eta)
run_stars_experiment("stars", settings_file, seeds, file.path(root, "stars"), code,
                     workers = 2L, eta = eta)
run_units_experiment("units", settings_file, file.path(root, "units"), c(0.1, 1, 10),
                     c("ebic", "scale"), 0.5, code, workers = 2L, eta = eta)
run_pilot_experiment("refit_pilot", settings_file, file.path(root, "refit_pilot"), code,
                     workers = 2L, eta = eta, seconds = 10L)
