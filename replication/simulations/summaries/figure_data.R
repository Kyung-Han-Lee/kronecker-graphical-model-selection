figure_criteria <- c("ebic", "scale", "relaxation")
figure3_methods <- c("scale", "relaxation", "cv", "stars", "stars_01")
figure_value_files <- c("figure1_curves", "figure1_medians", "figure2_path", "figure3_mcc")

figure1 <- function(p) {
    s <- p$figure1_scores
    if (length(unique(s$setting)) != 1L) stop("Figure 1 draws the scores of one setting")
    if (any(s$status != "ok")) stop("Figure 1: a candidate of the setting could not be scored")
    s <- s[order(match(s$criterion, figure_criteria), s$rep, s$k), ]
    s$centred <- s$ebic - stats::ave(s$ebic, s$criterion, s$rep, FUN = mean)

    one <- s[s$criterion == s$criterion[1], ]
    first_empty <- vapply(split(one, one$rep), function(x) {
        k <- x$k[x$edges == 0]
        if (length(k)) min(k) else NA_real_
    }, 0)
    if (anyNA(first_empty)) stop("Figure 1: a path without empty joint support")

    start <- round(stats::median(first_empty))
    medians <- do.call(rbind, lapply(intersect(figure_criteria, s$criterion), function(cr) {
        x <- s[s$criterion == cr, ]
        m <- tapply(x$centred, x$k, stats::median)
        k <- as.integer(names(m))
        data.frame(setting = x$setting[1], criterion = cr, k = k, eta = x$eta[match(k, x$k)],
                   median = as.vector(m), median_minimum = k == k[which.min(m)],
                   shading_start = k == start, stringsAsFactors = FALSE)
    }))
    rownames(s) <- rownames(medians) <- NULL
    list(figure1_curves = s[c("setting", "family", "degree", "dims", "n", "criterion", "gamma",
                              "rep", "k", "eta", "edges", "ebic", "centred")],
         figure1_medians = medians)
}

figure2_covariance <- matrix(c(1, 0.9, 0.9, 2), 2L)
glasso_2x2 <- function(s, lambda) {
    off <- sign(s[1, 2]) * max(abs(s[1, 2]) - lambda, 0)
    solve(matrix(c(s[1, 1] + lambda, off, off, s[2, 2] + lambda), 2L))
}

figure2 <- function() {
    s <- figure2_covariance
    lambda_max <- abs(s[1, 2])
    margin <- 0.02
    before <- seq(0, lambda_max * (1 - margin), length.out = 40L)
    after <- lambda_max * (1 + margin)
    grids <- list(a = c(before, exp(seq(log(after), log(60), length.out = 60L))),
                  b = c(before, seq(after, 6, length.out = 80L)))
    path <- do.call(rbind, lapply(names(grids), function(panel) {
        do.call(rbind, lapply(seq_along(grids[[panel]]), function(i) {
            lambda <- grids[[panel]][i]
            omega <- glasso_2x2(s, lambda)
            omega <- omega / sqrt(sum(omega^2))
            data.frame(panel = panel, point = i, lambda = lambda,
                       theta = atan2(omega[2, 2], omega[1, 1]) * (180 / pi),
                       z = sqrt(2) * omega[1, 2],
                       log_det = log(omega[1, 1] * omega[2, 2] - omega[1, 2]^2),
                       lambda_max = lambda_max, stringsAsFactors = FALSE)
        }))
    }))
    list(figure2_path = path)
}

figure3 <- function(p) {
    sel <- p$selected
    stars <- keep(sel, sel$experiment == "stars" & sel$threshold == 0.05 & sel$M == 2)
    parts <- list(scale = path_rows(sel, "scale", 0.5),
                  relaxation = path_rows(sel, "relaxation", 0.5), cv = cv_rows(sel),
                  stars = stars,
                  stars_01 = keep(sel, sel$experiment == "stars" &
                                      sel$threshold == 0.1 & sel$M == 2))
    x <- do.call(rbind, lapply(names(parts), function(m) {
        r <- keep(parts[[m]], parts[[m]]$setting %in% stars$setting)
        data.frame(r[c("setting", "family", "degree", "dims", "n", "rep", "id")], method = m,
                   mcc = r$mcc, stringsAsFactors = FALSE)
    }))
    reps <- tapply(x$rep, paste(x$setting, x$method), function(r) paste(sort(r), collapse = ","))
    by_setting <- split(reps, sub(" [^ ]+$", "", names(reps)))
    if (any(vapply(by_setting, function(r) length(r) != 5L || length(unique(r)) != 1L, TRUE))) {
        stop("Figure 3: the selection rules cover different replications of a setting")
    }
    x <- x[order(match(x$family, family_order), x$degree, x$n, match(x$method, figure3_methods),
                 x$rep), ]
    rownames(x) <- NULL
    list(figure3_mcc = x)
}

make_figures <- function(per_dataset) {
    files <- list()
    skipped <- list()
    if (is.null(per_dataset$figure1_scores)) {
        skipped$figure1 <- "no candidate scores of the setting of settings/figure1.csv"
    } else {
        files <- c(files, figure1(per_dataset))
    }
    files <- c(files, figure2())
    if (!"stars" %in% per_dataset$selected$experiment) {
        skipped$figure3 <- "no results of stars"
    } else {
        files <- c(files, figure3(per_dataset))
    }
    list(files = files, skipped = skipped)
}
