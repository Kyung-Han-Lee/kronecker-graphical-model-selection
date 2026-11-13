source("helper_lib.R")
for (file in list.files(file.path("..", "..", "summaries"), pattern = "[.]R$", full.names = TRUE)) {
    source(file)
}
source("helper_outputs.R")
x <- collect_outputs(root, settings_dir)
sel <- x$tables$selected

count <- table(sel$experiment)
expect_equal(as.vector(count[c("matrix", "cv", "stars", "units", "refit_pilot")]),
             c(4 * 3 * 2, 4, 4 * 2, 4 * 3 * 2, 4 * 2))
expect_identical(vapply(x$sources, `[[`, "", "folder"),
                 c("matrix", file.path("cv", "matrix"), "stars", "units", "refit_pilot"))
expect_true(all(sel$grid_points == 5L))
expect_true(all(sel$gamma[sel$experiment == "refit_pilot"] == 0.5))
expect_true(all(is.na(sel$gamma[sel$experiment %in% c("cv", "stars")])))
expect_identical(sort(unique(sel$threshold[sel$experiment == "stars"])), c(0.05, 0.1))
expect_identical(sort(unique(sel$constant[sel$experiment == "units"])), c(0.1, 1, 10))
expect_true(all(is.na(sel$constant[sel$experiment != "units"])))
ba <- sel[sel$id == "t_ba_r002", ]
expect_true(all(ba$M == 2L & ba$family == "ba" & ba$dims == "8x6" & ba$degree == 1 &
                    ba$n == 20L & ba$rep == 2L))

r <- readRDS(file.path(root, "matrix", "datasets", "t_ba_r002.rds"))
row <- sel[sel$experiment == "matrix" & sel$id == "t_ba_r002" & sel$criterion == "scale" &
               sel$gamma == 1, ]
want <- r$selected[r$selected$criterion == "scale" & r$selected$gamma == 1, ]
columns <- c("k", "eta", "edges", "empty", "tp", "fp", "fn", "tn", "mcc", "fdp", "tpr")
expect_equal(unlist(row[columns]), unlist(want[columns]))
modes <- r$selected_modes[r$selected_modes$criterion == "scale" & r$selected_modes$gamma == 1, ]
modes <- modes[order(modes$mode), ]
expect_equal(unname(unlist(row[c("tp_1", "fp_1", "fn_1", "tn_1", "tp_2", "fp_2", "fn_2", "tn_2")])),
             unname(c(unlist(modes[1, c("tp", "fp", "fn", "tn")]),
                      unlist(modes[2, c("tp", "fp", "fn", "tn")]))))
for (count_col in c("tp", "fp", "fn", "tn")) {
    expect_equal(sel[[count_col]], sel[[paste0(count_col, "_1")]] + sel[[paste0(count_col, "_2")]])
}

paths <- x$tables$paths
expect_identical(paths$id, sort(ids))
expect_true(all(paths$candidates == 5L))
exact <- r$recovery$fp == 0 & r$recovery$fn == 0
expect_equal(paths$exact_candidates[paths$id == "t_ba_r002"], sum(exact))
expect_identical(paths$adjacent_exact[paths$id == "t_ba_r002"], adjacent_pair(exact))
expect_true(adjacent_pair(c(FALSE, TRUE, TRUE, FALSE)))
expect_false(adjacent_pair(c(TRUE, FALSE, TRUE)))
expect_false(adjacent_pair(TRUE))

f1 <- x$tables$figure1_scores
expect_identical(unique(f1$id), c("t_er_r001", "t_er_r002"))
expect_true(all(f1$experiment == "matrix" & f1$gamma == 0.5))
expect_equal(nrow(f1), 2L * 3L * 5L)
r1 <- readRDS(file.path(root, "matrix", "datasets", "t_er_r002.rds"))
want <- r1$scores[r1$scores$gamma == 0.5, ]
got <- f1[f1$id == "t_er_r002", ]
m <- match(paste(want$criterion, want$k), paste(got$criterion, got$k))
expect_equal(got$ebic[m], want$ebic, tolerance = 0)
expect_equal(got$edges[m], want$edges)
expect_identical(got$status[m], want$status)

pilot <- x$tables$pilot_candidates
expect_equal(nrow(pilot), 4L * 5L)
p1 <- readRDS(file.path(root, "refit_pilot", "datasets", "t_er_r001.rds"))
expect_identical(pilot$refit_status[pilot$id == "t_er_r001"], p1$candidates$refit_status)
expect_identical(pilot$relaxation_status[pilot$id == "t_er_r001"],
                 p1$candidates$relaxation_status)

three <- data.frame(criterion = "scale", gamma = 0.5, k = 2L)
three_modes <- data.frame(criterion = "scale", gamma = 0.5, k = 2L, mode = c(2L, 1L, 3L),
                          tp = c(2, 1, 3), fp = 0, fn = 0, tn = 5)
wide <- widen_modes(three, three_modes, 4L)
expect_equal(unlist(wide[c("tp_1", "tp_2", "tp_3")]), c(tp_1 = 1, tp_2 = 2, tp_3 = 3))
expect_true(all(is.na(wide[c("tp_4", "fp_4", "fn_4", "tn_4")])))

tricky <- c(0.1, 1 / 3, pi * 1e-20, 2^-1074, 1e300, -0.75, NA)
expect_identical(as.numeric(exact_text(tricky)), tricky)
expect_identical(exact_text(c(0.1, 40, 1e-04, 1e-05)), c("0.1", "40", "0.0001", "1e-05"))

d <- tempfile()
write_per_dataset(x, d)
y <- read_per_dataset(d)
expect_identical(sort(names(y)), c("figure1_scores", "paths", "pilot_candidates", "selected"))
for (name in names(y)) {
    expect_identical(names(y[[name]]), names(x$tables[[name]]))
    for (col in names(y[[name]])) {
        expect_equal(y[[name]][[col]], x$tables[[name]][[col]], tolerance = 0,
                     info = paste(name, col))
    }
}
expect_identical(jsonlite::read_json(file.path(d, "sources.json"))[[1]]$experiment, "matrix")
x2 <- x
x2$tables$paths <- NULL
write_per_dataset(x2, d)
expect_false(file.exists(file.path(d, "paths.csv.gz")))

write.csv(transform(tiny, w = 0.31), settings_file, row.names = FALSE, na = "")
expect_error(collect_outputs(root, settings_dir), "made with another matrix.csv")
