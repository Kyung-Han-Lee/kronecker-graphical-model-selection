source("helper_lib.R")
for (file in list.files(file.path("..", "..", "summaries"), pattern = "[.]R$", full.names = TRUE)) {
    source(file)
}

design <- data.frame(setting = c("er10", "er50", "ba10", "ba50"), M = 2L,
                     family = c("er", "er", "ba", "ba"), dims = "8x6",
                     degree = c(2, 2, 1, 1), n = c(10L, 50L, 10L, 50L),
                     stringsAsFactors = FALSE)
datasets <- merge(design, data.frame(rep = 1:3))
datasets$id <- sprintf("%s_r%03d", datasets$setting, datasets$rep)
rows <- function(experiment, criterion, gamma = NA_real_, constant = NA_real_,
                 threshold = NA_real_, settings = design$setting) {
    d <- datasets[datasets$setting %in% settings, ]
    data.frame(experiment = experiment, d, criterion = criterion, gamma = gamma,
               constant = constant, threshold = threshold, grid_points = 5L, k = 3L, eta = 1,
               edges = 2, empty = FALSE, tp = 1, fp = 1, fn = 1, tn = 10, mcc = d$rep / 10,
               fdp = d$rep / 10, tpr = 1 - d$rep / 10, stringsAsFactors = FALSE)
}

matrix_rows <- do.call(rbind, lapply(c(0.5, 1), function(g) {
    ebic <- transform(rows("matrix", "ebic", g), k = 5L, edges = 0, empty = TRUE, mcc = 0,
                      fdp = 0, tpr = 0)
    scale <- transform(rows("matrix", "scale", g), mcc = mcc + 0.01)
    if (g == 1) {
        first <- scale$id == "er10_r001"
        scale$k[first] <- 1L
        scale$edges[first] <- 0
        scale$empty[first] <- TRUE
    }
    relaxation <- transform(rows("matrix", "relaxation", g), fp = ifelse(rep == 3, 1, 0),
                            fn = ifelse(rep == 1, 0, 1))
    rbind(ebic, scale, relaxation)
}))
units_rows <- rbind(
    transform(rows("units", "ebic", 0.5, constant = 0.1), k = 5L, empty = TRUE, mcc = 0),
    transform(rows("units", "ebic", 0.5, constant = 1), k = ifelse(rep == 3, 4L, 3L)),
    transform(rows("units", "ebic", 0.5, constant = 10), k = 1L),
    rows("units", "scale", 0.5, constant = 1),
    rows("units", "relaxation", 0.5, constant = 1))
sel <- rbind(matrix_rows,
             transform(rows("cv", "cv"), fdp = 0.5, tpr = 0.95, mcc = 0.4),
             rows("stars", "stars", threshold = 0.05, settings = c("er10", "er50")),
             transform(rows("stars", "stars", threshold = 0.1, settings = c("er10", "er50")),
                       mcc = 0.3),
             units_rows)

paths <- data.frame(experiment = "matrix", datasets[c("id", "setting", "rep", "M", "family",
                                                      "dims", "degree", "n")],
                    candidates = 5L, exact_candidates = c(1, 0, 2)[datasets$rep],
                    adjacent_exact = datasets$rep == 3)
p <- list(selected = sel, paths = paths)

t1 <- table_b4(p)
expect_identical(t1$setting, c("er10", "er50"))
expect_identical(names(t1)[1:7], c("setting", "M", "family", "degree", "dims", "n",
                                   "replications"))
expect_identical(names(t1)[8:11], c("fdp_mean_scale", "fdp_sd_scale", "fdp_mean_relaxation",
                                    "fdp_sd_relaxation"))
expect_identical(names(t1)[26:27], c("tpr_mean_stars_01", "tpr_sd_stars_01"))
expect_true(all(c("fdp_mean_stars_01", "fdp_sd_stars_01") %in% names(t1)))
expect_equal(t1$replications, c(3L, 3L))
expect_equal(t1$fdp_mean_scale, c(0.2, 0.2))
expect_equal(t1$fdp_sd_scale, c(0.1, 0.1))
expect_equal(t1$fdp_sd_cv, c(0, 0))
expect_equal(t1$tpr_mean_relaxation, c(0.8, 0.8))
expect_equal(t1$fdp_mean_stars, c(0.2, 0.2))

three <- data.frame(experiment = "tensor", id = "t3_r001", setting = "t3", rep = 1L, M = 3L,
                    family = "er", dims = "5x4x3", degree = 1, n = 20L,
                    stringsAsFactors = FALSE)
p2 <- list(selected = rbind(sel, transform(sel[1, ], experiment = "tensor", id = "t3_r001",
                                           setting = "t3", rep = 1L, M = 3L, dims = "5x4x3",
                                           degree = 1, n = 20L, criterion = "relaxation",
                                           gamma = 0.5, fp = 0, fn = 0)),
           paths = rbind(paths, cbind(three, candidates = 5L, exact_candidates = 2,
                                      adjacent_exact = FALSE)))
t2 <- table_b3(p2)
expect_identical(paste(t2$group, t2$M),
                 c("any 2", "any 3", "adjacent 2", "adjacent 3", "single 2", "single 3",
                   "separated 2", "separated 3"))
expect_equal(t2$paths, c(8, 1, 4, 0, 4, 0, 0, 1))
expect_equal(t2$selected, c(4, 1, 0, 0, 4, 0, 0, 1))
expect_equal(t2$selected_percent, c(50, 100, 0, NA, 100, NA, NA, 100))
expect_error(table_b3(list(selected = sel)), "needs the candidate recovery")
expect_error(table_b3(list(selected = sel[sel$criterion != "relaxation", ], paths = paths)),
             "without a Relaxation selection")

t3 <- table_b1(p)
expect_identical(paste(t3$criterion, t3$gamma),
                 c("ebic 0.5", "ebic 1", "scale 0.5", "scale 1", "relaxation 0.5",
                   "relaxation 1"))
expect_equal(t3$interior_percent, c(0, 0, 100, 100 * 11 / 12, 100, 100))
expect_equal(t3$nonempty_n10, c(0, 0, 6, 5, 6, 6))
expect_equal(t3$nonempty_n10_percent, c(0, 0, 100, 100 * 5 / 6, 100, 100))
expect_equal(t3$datasets_n30plus, rep(6L, 6))
expect_equal(t3$nonempty_n30plus_percent, c(0, 0, 100, 100, 100, 100))

t4 <- table_b2(p)
expect_identical(paste(t4$criterion, t4$constant),
                 c("ebic 0.1", "ebic 1", "ebic 10", "scale 1", "relaxation 1"))
expect_equal(t4$largest_percent, c(100, 0, 0, 0, 0))
expect_equal(t4$interior_percent, c(0, 100, 0, 100, 100))
expect_equal(t4$smallest_percent, c(0, 0, 100, 0, 0))
expect_equal(t4$empty_percent, c(100, 0, 0, 0, 0))
expect_equal(t4$same_as_scale, c(0, 8, 0, NA, NA))
expect_equal(t4$same_as_scale_percent, c(0, 100 * 8 / 12, 0, NA, NA))
expect_equal(t4$mcc_mean, c(0, 0.2, 0.2, 0.2, 0.2))
expect_equal(t4$mcc_sd, c(0, rep(sd(rep(c(0.1, 0.2, 0.3), 4)), 4)))

t7 <- table_b6(p)
expect_identical(names(t7)[8:11], c("mcc_mean_scale", "mcc_sd_scale", "mcc_mean_relaxation",
                                    "mcc_sd_relaxation"))
expect_equal(t7$mcc_mean_scale, rep(0.21, 4))
expect_equal(t7$mcc_sd_scale, rep(0.1, 4))
expect_equal(t7$mcc_mean_cv, rep(0.4, 4))
expect_equal(t7$mcc_sd_cv, rep(0, 4))
expect_equal(t7$tpr_mean_relaxation, rep(0.8, 4))
t9 <- table_b7(p)
expect_identical(names(t9)[8:11], c("mcc_mean_scale_gamma0.5", "mcc_sd_scale_gamma0.5",
                                    "mcc_mean_relaxation_gamma0.5",
                                    "mcc_sd_relaxation_gamma0.5"))
expect_equal(t9$mcc_mean_relaxation_gamma1, rep(0.2, 4))
expect_equal(t9$mcc_sd_relaxation_gamma1, rep(0.1, 4))

short <- p
short$selected <- sel[!(sel$experiment == "cv" & sel$id == "er10_r003"), ]
expect_error(table_b6(short), "cv covers other replications of er10")
twice <- p
twice$selected <- rbind(sel, sel[sel$experiment == "cv" & sel$id == "er10_r003", ])
expect_error(table_b6(twice), "setting er10 repeats a replication")

pilot <- data.frame(id = rep(c("p2_r001", "p2_r002"), each = 3), setting = "p2",
                    rep = rep(1:2, each = 3), M = 2L, family = "er", dims = "12x3", degree = 1,
                    n = 2L, k = rep(1:3, 2), eta = 1, edges = 3,
                    relaxation_status = c("not_pd", "ok", "ok", "ok", "ok", "ok"),
                    refit_status = c("nonexistent", "ok", "timeout", "ok", "error",
                                     "nonconverged"), stringsAsFactors = FALSE)
pilot_sel <- data.frame(experiment = "refit_pilot", id = rep(c("p2_r001", "p2_r002"), 2),
                        setting = "p2", rep = rep(1:2, 2), M = 2L, family = "er", dims = "12x3",
                        degree = 1, n = 2L, criterion = rep(c("relaxation", "refit"), each = 2),
                        gamma = 0.5, constant = NA, threshold = NA, grid_points = 3L,
                        k = c(2L, 3L, 2L, NA), eta = 1, edges = 1, empty = FALSE, tp = 1, fp = 0,
                        fn = 0, tn = 1, mcc = c(0.5, 0.7, 0.6, NA), fdp = 0, tpr = 1,
                        stringsAsFactors = FALSE)
t6 <- table_b5(list(selected = pilot_sel, pilot_candidates = pilot))
expect_equal(unlist(t6[c("replications", "candidates", "refit_ok", "refit_timeout",
                         "refit_error", "refit_nonconverged", "refit_nonexistent",
                         "relaxation_not_scored", "both_selected", "same_candidate")]),
             c(replications = 2, candidates = 6, refit_ok = 2, refit_timeout = 1,
               refit_error = 1, refit_nonconverged = 1, refit_nonexistent = 1,
               relaxation_not_scored = 1, both_selected = 1, same_candidate = 1))
expect_equal(c(t6$mcc_mean_relaxation, t6$mcc_mean_refit), c(0.5, 0.6))

expect_true(is.na(t6$mcc_sd_relaxation) && is.na(t6$mcc_sd_refit))
t6b <- table_b5(list(selected = transform(pilot_sel, k = c(2L, 3L, 2L, 3L),
                                        mcc = c(0.5, 0.7, 0.6, 0.8)),
                   pilot_candidates = pilot))
expect_equal(c(t6b$mcc_sd_relaxation, t6b$mcc_sd_refit), rep(sd(c(0.5, 0.7)), 2))
odd <- transform(pilot, refit_status = replace(refit_status, 1, "lost"))
expect_error(table_b5(list(selected = pilot_sel, pilot_candidates = odd)),
             "unknown Refit status lost")

made <- make_tables(p)
expect_identical(names(made$tables), c("TableB4", "TableB1", "TableB2", "TableB6", "TableB7"))
expect_identical(made$skipped, list(TableB3 = "no results of tensor",
                                    TableB5 = "no results of refit_pilot",
                                    Table1 = "no results of tensor"))
