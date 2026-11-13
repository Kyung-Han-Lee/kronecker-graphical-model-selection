source("helper_lib.R")
tiny <- data.frame(setting = c("t_er", "t_ba"), M = 2L, family = c("er", "ba"),
                   dims = "8x6", d = c(1, NA), m_new = c(NA, 1L), n = 20L,
                   n_idx = 1L, reps = 3L, w = 0.3, w_spread = 0.25, margin = 0.5,
                   gradient = "numeric", stringsAsFactors = FALSE)
f <- tempfile(fileext = ".csv")
write.csv(tiny, f, row.names = FALSE, na = "")
out <- tempfile()
eta <- kgm_eta_grid(5L, c(0.01, 10))
code <- driver_code_hash(file.path("..", "..", "R"))
run_experiment("tiny", f, out, c("ebic", "scale"), 0.5, code, workers = 2L, reps = 2L, eta = eta)
expect_identical(sort(list.files(file.path(out, "datasets"))),
                 sort(paste0(c("t_er_r001", "t_er_r002", "t_ba_r001", "t_ba_r002"), ".rds")))
r <- readRDS(file.path(out, "datasets", "t_er_r001.rds"))
expect_equal(nrow(r$scores), 2L * 5L)
expect_equal(nrow(r$selected), 2L)
gen <- make_dataset(read_settings(f)[1, ], 1L)
expect_identical(r$data_sha256, digest::digest(gen$Y, algo = "sha256"))
expect_equal(jsonlite::read_json(file.path(out, "run.json"))$datasets$done, 4L)
expect_identical(jsonlite::read_json(file.path(out, "run.json"))$method,
                 list(name = "path", criteria = list("ebic", "scale"), gammas = 0.5))

mtime <- file.mtime(file.path(out, "datasets", "t_er_r001.rds"))
run_experiment("tiny", f, out, c("ebic", "scale"), 0.5, code, workers = 2L, reps = 2L, eta = eta)
expect_equal(jsonlite::read_json(file.path(out, "run.json"))$datasets$reused, 4L)
expect_identical(file.mtime(file.path(out, "datasets", "t_er_r001.rds")), mtime)
expect_error(run_experiment("tiny", f, out, c("ebic", "scale"), c(0, 0.5), code,
                            workers = 1L, reps = 1L, eta = eta), "different setup")
expect_identical(list.files(file.path(out, "errors")), character())
expect_equal(length(list.files(file.path(out, "datasets"))), 4L)

r_other <- r
r_other$signature <- "other"
saveRDS(r_other, file.path(out, "datasets", "t_er_r001.rds"))
expect_error(run_experiment("tiny", f, out, c("ebic", "scale"), 0.5, code, workers = 1L,
                            reps = 1L, eta = eta), "1 datasets failed")
expect_true(grepl("different setup", readLines(file.path(out, "errors", "t_er_r001.txt"))))
saveRDS(r, file.path(out, "datasets", "t_er_r001.rds"))
unlink(file.path(out, "errors", "t_er_r001.txt"))

lib2 <- tempfile()
dir.create(lib2)
file.copy(list.files(file.path("..", "..", "R"), full.names = TRUE), lib2)
score <- readLines(file.path(lib2, "scoring.R"))
writeLines(sub("which.min(scores$ebic[usable])", "which.max(scores$ebic[usable])",
               score, fixed = TRUE), file.path(lib2, "scoring.R"))
code2 <- driver_code_hash(lib2)
expect_false(identical(code2, code))
expect_error(run_experiment("tiny", f, out, c("ebic", "scale"), 0.5, code2, workers = 1L,
                            reps = 2L, eta = eta), "different setup")
unlink(file.path(out, "run.json"))
expect_error(run_experiment("tiny", f, out, c("ebic", "scale"), 0.5, code2, workers = 1L,
                            reps = 2L, eta = eta), "4 datasets failed")

out2 <- tempfile()
stored <- kgm_path(gen$Y, eta)
run_experiment("tiny", f, out2, "scale", 0.5, code, workers = 1L, ids = "t_er_r001",
               eta = eta, path_fun = function(dataset, Y, eta) stored)
r2 <- readRDS(file.path(out2, "datasets", "t_er_r001.rds"))
expect_identical(r2$selected$k, r$selected$k[r$selected$criterion == "scale"])
expect_equal(length(list.files(file.path(out2, "datasets"))), 1L)

out3 <- tempfile()
run_experiment("tiny", f, out3, "refit", 0.5, code, ids = "t_er_r001", eta = eta,
               support = TRUE)
r3 <- readRDS(file.path(out3, "datasets", "t_er_r001.rds"))
expect_identical(r3$candidates, support_records(kgm_path(gen$Y, eta)))
expect_true(jsonlite::read_json(file.path(out3, "run.json"))$method$support)
expect_null(r$candidates)

expect_identical(r$recovery, candidate_recovery(kgm_path(gen$Y, eta), gen$truth))
