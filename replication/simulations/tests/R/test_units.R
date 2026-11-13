source("helper_lib.R")
tiny <- data.frame(setting = "t_er", M = 2L, family = "er", dims = "8x6", d = 1,
                   m_new = NA, n = 30L, n_idx = 1L, reps = 2L, w = 0.3,
                   w_spread = 0.25, margin = 0.5, gradient = "numeric",
                   stringsAsFactors = FALSE)
f <- tempfile(fileext = ".csv")
write.csv(tiny, f, row.names = FALSE, na = "")
eta <- kgm_eta_grid(6L, c(0.01, 10))
code <- driver_code_hash(file.path("..", "..", "R"))
constants <- c(0.1, 1, 10)
out <- tempfile()
run_units_experiment("tiny", f, out, constants, c("ebic", "scale", "relaxation"), 0.5,
                     code, workers = 2L, eta = eta)
r <- readRDS(file.path(out, "datasets", "t_er_r001.rds"))
expect_identical(sort(unique(r$selected$constant)), constants)
expect_identical(jsonlite::read_json(file.path(out, "run.json"))$method$name, "units")

gen <- make_dataset(read_settings(f)[1, ], 1L)
for (cc in constants) {
    path <- kgm_path(cc * gen$Y, cc^2 * eta)
    sc <- score_path(path, cc * gen$Y, c("ebic", "scale", "relaxation"), 0.5, "numeric")
    got <- r$scores[r$scores$constant == cc, names(sc)]
    rownames(got) <- NULL
    expect_identical(got, sc, info = paste("constant", cc))
}

pick <- function(crit, cc) r$selected$k[r$selected$criterion == crit & r$selected$constant == cc]
expect_identical(vapply(constants, function(cc) pick("scale", cc), 0L), rep(pick("scale", 1), 3L))
expect_identical(vapply(constants, function(cc) pick("relaxation", cc), 0L),
                 rep(pick("relaxation", 1), 3L))
one <- r$scores[r$scores$constant == 1, ]
np <- 30 * 48
for (cc in c(0.1, 10)) {
    predicted <- one$ebic[one$criterion == "ebic"] +
        np * (cc^2 - 1) / one$tau[one$criterion == "scale"]
    got <- r$scores$ebic[r$scores$constant == cc & r$scores$criterion == "ebic"]
    expect_equal(got, predicted, tolerance = 1e-10, info = paste("constant", cc))
}

run_units_experiment("tiny", f, out, constants, c("ebic", "scale", "relaxation"), 0.5,
                     code, workers = 2L, eta = eta)
expect_equal(jsonlite::read_json(file.path(out, "run.json"))$datasets$reused, 2L)
expect_error(run_units_experiment("tiny", f, out, c(0.1, 1), c("ebic", "scale", "relaxation"),
                                  0.5, code, eta = eta), "different setup")
