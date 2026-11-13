source("helper_lib.R")

expect_identical(stars_seed(c("cell_ba_p116x102_m1_n0010_see6069df88dd_rep001",
                              "cell_er_p116x102_d2_n3200_sd19cc38449a8_rep100")),
                 c(138613999L, 252453985L))

tiny <- data.frame(setting = "t_er", M = 2L, family = "er", dims = "8x6", d = 1,
                   m_new = NA, n = 20L, n_idx = 1L, reps = 2L, w = 0.3,
                   w_spread = 0.25, margin = 0.5, gradient = "numeric",
                   stringsAsFactors = FALSE)
sf <- tempfile(fileext = ".csv")
write.csv(tiny, sf, row.names = FALSE, na = "")
seeds <- data.frame(dataset = c("t_er_r001", "t_er_r002"), source = c("a", "b"),
                    seed = stars_seed(c("a", "b")))
tf <- tempfile(fileext = ".csv")
write.csv(seeds, tf, row.names = FALSE)
expect_identical(read_stars_seeds(tf), seeds)
bad <- seeds
bad$dataset[2] <- bad$dataset[1]
bf <- tempfile(fileext = ".csv")
write.csv(bad, bf, row.names = FALSE)
expect_error(read_stars_seeds(bf), "repeats")

eta <- kgm_eta_grid(5L, c(0.01, 10))
code <- driver_code_hash(file.path("..", "..", "R"))
out <- tempfile()
run_stars_experiment("tiny", sf, tf, out, code, workers = 2L, eta = eta, N = 4L)
r <- readRDS(file.path(out, "datasets", "t_er_r002.rds"))
gen <- make_dataset(read_settings(sf)[1, ], 2L)
path <- kgm_path(gen$Y, eta)
st <- kgm_stars(path, gen$Y, N = 4L, threshold = c(0.05, 0.1), seed = seeds$seed[2],
                tol = 0)
expect_identical(r$stars_seed, seeds$seed[2])
expect_identical(r$instability, st$instability)
expect_identical(r$subsample_size, st$subsample_size)
expect_identical(r$selected$threshold, c(0.05, 0.1))
expect_identical(r$selected$k, st$index)
m <- kgm_metrics(gen$truth, path$fits[[st$index[2]]], tol = 0)
expect_identical(r$selected$mcc[2], m$mcc)
expect_identical(r$selected_modes$tp[r$selected_modes$threshold == 0.1], m$per_mode$tp)
expect_identical(jsonlite::read_json(file.path(out, "run.json"))$method$name, "stars")

run_stars_experiment("tiny", sf, tf, out, code, workers = 2L, eta = eta, N = 4L)
expect_equal(jsonlite::read_json(file.path(out, "run.json"))$datasets$reused, 2L)
seeds2 <- seeds
seeds2$seed[1] <- seeds2$seed[1] + 1L
tf2 <- tempfile(fileext = ".csv")
write.csv(seeds2, tf2, row.names = FALSE)
expect_error(run_stars_experiment("tiny", sf, tf2, out, code, eta = eta, N = 4L),
             "different setup")

seeds3 <- rbind(seeds, data.frame(dataset = "t_er_r003", source = "c", seed = 1L))
tf3 <- tempfile(fileext = ".csv")
write.csv(seeds3, tf3, row.names = FALSE)
expect_error(run_stars_experiment("tiny", sf, tf3, tempfile(), code, eta = eta, N = 4L),
             "not in")

settings_dir <- file.path("..", "..", "settings")
shipped <- read_stars_seeds(file.path(settings_dir, "stars_seeds.csv"))
expect_equal(nrow(shipped), 1400L)
per_setting <- table(sub("_r[0-9]+$", "", shipped$dataset))
expect_identical(as.vector(per_setting), rep(100L, 14L))
expect_true(all(shipped$dataset %in%
                    expand_datasets(read_settings(file.path(settings_dir, "matrix.csv")))$id))
expect_identical(stars_seed(shipped$source), shipped$seed)
