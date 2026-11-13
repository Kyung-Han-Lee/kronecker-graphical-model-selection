source("helper_lib.R")
for (file in list.files(file.path("..", "..", "summaries"), pattern = "[.]R$", full.names = TRUE)) {
    source(file)
}

base <- list(ebic = c(4, 3, 2, 1), scale = c(3, 1, 2, 4), relaxation = c(5, 4, 1, 2))
grid <- expand.grid(k = 1:4, rep = 1:2, criterion = c("relaxation", "ebic", "scale"),
                    stringsAsFactors = FALSE)
scores <- data.frame(experiment = "matrix", id = sprintf("s_r%03d", grid$rep), setting = "s",
                     rep = grid$rep, M = 2L, family = "er", dims = "8x6", degree = 2, n = 50L,
                     criterion = grid$criterion, gamma = 0.5, k = grid$k, eta = 10^(grid$k - 3),
                     edges = ifelse(grid$k > grid$rep, 0, 5),
                     ebic = unlist(Map(function(cr, k, r) base[[cr]][k] * r, grid$criterion,
                                       grid$k, grid$rep)),
                     status = "ok", stringsAsFactors = FALSE)
f1 <- figure1(list(figure1_scores = scores))
curves <- f1$figure1_curves
expect_identical(unique(curves$criterion), c("ebic", "scale", "relaxation"))
expect_equal(curves$centred[curves$criterion == "ebic" & curves$rep == 2], c(3, 1, -1, -3))
medians <- f1$figure1_medians
expect_equal(medians$median[medians$criterion == "scale"], c(0.75, -2.25, -0.75, 2.25))
expect_identical(medians$k[medians$median_minimum], c(4L, 2L, 3L))
expect_identical(medians$k[medians$shading_start], c(2L, 2L, 2L))
expect_equal(medians$eta[medians$criterion == "relaxation"], 10^(-2:1))

expect_error(figure1(list(figure1_scores = transform(scores, edges = 5))),
             "without empty joint support")
expect_error(figure1(list(figure1_scores = transform(scores, status = replace(status, 1,
                                                                              "not_pd")))),
             "could not be scored")
expect_error(figure1(list(figure1_scores = rbind(scores, transform(scores, setting = "t")))),
             "one setting")

f2 <- figure2()$figure2_path
expect_equal(as.vector(table(f2$panel)), c(100, 120))
expect_true(all(f2$lambda_max == 0.9))
a <- f2[f2$panel == "a", ]
expect_equal(range(a$lambda), c(0, 60))
expect_equal(range(f2$lambda[f2$panel == "b"]), c(0, 6))
expect_false(any(f2$lambda > 0.9 * 0.98 & f2$lambda < 0.9 * 1.02))
inv <- solve(figure2_covariance)
inv <- inv / sqrt(sum(inv^2))
expect_equal(unlist(a[1, c("theta", "z", "log_det")]),
             c(theta = atan2(inv[2, 2], inv[1, 1]) * 180 / pi, z = sqrt(2) * inv[1, 2],
               log_det = log(det(inv))))
expect_true(all(f2$z[f2$lambda > 0.9] == 0))
expect_true(all(f2$log_det <= -log(2)))
expect_true(abs(a$theta[100] - 45) < 0.5 && abs(a$log_det[100] + log(2)) < 1e-3)

for (lambda in c(0.1, 0.5, 0.89, 0.9, 2)) {
    omega <- glasso_2x2(figure2_covariance, lambda)
    w <- solve(omega)
    expect_equal(diag(w), diag(figure2_covariance) + lambda)
    if (omega[1, 2] != 0) {
        expect_equal(w[1, 2], figure2_covariance[1, 2] + lambda * sign(omega[1, 2]))
    } else {
        expect_true(abs(w[1, 2] - figure2_covariance[1, 2]) <= lambda)
    }
}

design <- data.frame(setting = c("ba50", "er50", "er10"), M = 2L, family = c("ba", "er", "er"),
                     dims = "8x6", degree = c(1, 2, 2), n = c(50L, 50L, 10L),
                     stringsAsFactors = FALSE)
datasets <- merge(design, data.frame(rep = 1:2))
datasets$id <- sprintf("%s_r%03d", datasets$setting, datasets$rep)
rows <- function(experiment, criterion, mcc, gamma = NA_real_, threshold = NA_real_,
                 settings = design$setting) {
    d <- datasets[datasets$setting %in% settings, ]
    data.frame(experiment = experiment, d, criterion = criterion, gamma = gamma,
               constant = NA_real_, threshold = threshold, grid_points = 5L, k = 3L,
               eta = 1, edges = 2, empty = FALSE, tp = 1, fp = 1, fn = 1, tn = 10,
               mcc = mcc + d$rep / 100, fdp = 0.5, tpr = 0.5, stringsAsFactors = FALSE)
}
sel <- rbind(rows("matrix", "scale", 0.1, gamma = 0.5), rows("matrix", "scale", 0.9, gamma = 1),
             rows("matrix", "relaxation", 0.2, gamma = 0.5), rows("cv", "cv", 0.3),
             rows("stars", "stars", 0.4, threshold = 0.05, settings = c("ba50", "er50")),
             rows("stars", "stars", 0.8, threshold = 0.1, settings = c("ba50", "er50")))
f3 <- figure3(list(selected = sel))$figure3_mcc
expect_identical(unique(f3$setting), c("er50", "ba50"))
expect_identical(f3$method[f3$setting == "er50"],
                 rep(c("scale", "relaxation", "cv", "stars", "stars_01"), each = 2))
expect_equal(f3$mcc[f3$setting == "er50"],
             rep(c(0.1, 0.2, 0.3, 0.4, 0.8), each = 2) + c(0.01, 0.02))
short <- sel[!(sel$experiment == "cv" & sel$id == "er50_r002"), ]
expect_error(figure3(list(selected = short)), "cover different replications")

made <- make_figures(list(selected = sel[sel$experiment != "stars", ]))
expect_identical(names(made$files), "figure2_path")
expect_identical(made$skipped,
                 list(figure1 = "no candidate scores of the setting of settings/figure1.csv",
                      figure3 = "no results of stars"))
made <- make_figures(list(selected = sel, figure1_scores = scores))
expect_identical(names(made$files), figure_value_files)
expect_identical(made$skipped, list())
