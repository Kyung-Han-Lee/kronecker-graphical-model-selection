family_order <- c("er", "ba")
criterion_order <- c("ebic", "scale", "relaxation")
refit_statuses <- c("ok", "timeout", "error", "nonconverged", "nonexistent")

table_notes <- data.frame(
    table = c("TableB4", "TableB3", "TableB1", "TableB2", "TableB5", "TableB6", "Table1", "TableB7"),
    builder = c("table_b4", "table_b3", "table_b1", "table_b2", "table_b5", "table_b6", "table_1", "table_b7"),
    label = c("tab:sim-fdr-tpr", "tab:sim-true-model", "tab:app-m2-selection", "tab:app-m2-units", "tab:app-refit-pilot", "tab:app-m2-support", "tab:app-tensor-support", "tab:app-m2-gamma"),
    caption = c("Mean FDP and TPR for (116,102) matrices", "Selection of the true model by Relaxation on the fitted paths that contain it", "Selection frequencies (%) in the matrix simulations", "Selection frequencies (%) at gamma = 0.5 after multiplying the data by a constant and the grid by its square", "Refit and Relaxation at small sample sizes for (p1,p2) = (100,5), ER graphs with d = 2 and gamma = 0.5", "Support recovery for matrix data at gamma = 0.5", "Support recovery for three- and four-way tensors", "MCC of Scale and Relaxation for matrix data at four values of gamma"),
    printed = c("means with SDs in parentheses, three decimals", "percentage selected with one decimal and the counts selected/paths with thousands separated by commas; rows any, adjacent and single by M (the separated row is not printed)", "percentages, whole numbers without decimals and others with two; rows of several gammas with equal values are printed once", "percentages as in Table B.1, and mean MCC with three decimals (-0.000 as 0.000) and its SD in parentheses; constants 0.1, 1, 10 and 100; Scale and Relaxation once, as 'all', since their rows are equal for every constant", "counts, and mean MCC with three decimals (-0.000 as 0.000) and its SD in parentheses", "means with SDs in smaller parentheses, three decimals; n = 10, 100 and 3200", "means with SDs in smaller parentheses, three decimals", "means with SDs in smaller parentheses, three decimals; n = 10, 100 and 3200"),
    needs = c("matrix,cv,stars", "matrix,tensor", "matrix", "units", "refit_pilot", "matrix,cv", "tensor,cv", "matrix"),
    stringsAsFactors = FALSE)

keep <- function(x, condition) x[!is.na(condition) & condition, , drop = FALSE]
percent <- function(count, total) if (total > 0) 100 * count / total else NA_real_

path_rows <- function(sel, criterion, gamma) {
    keep(sel, sel$experiment %in% c("matrix", "tensor") & sel$criterion == criterion &
             sel$gamma == gamma)
}
cv_rows <- function(sel) keep(sel, sel$experiment == "cv")

setting_summary <- function(x, metrics, sd) {
    do.call(rbind, lapply(split(x, x$setting), function(g) {
        if (anyDuplicated(g$rep)) stop("setting ", g$setting[1], " repeats a replication")
        out <- data.frame(setting = g$setting[1], replications = nrow(g),
                          reps = paste(sort(g$rep), collapse = ","), stringsAsFactors = FALSE)
        for (m in metrics) {
            out[[paste0(m, "_mean")]] <- mean(g[[m]])
            if (sd) out[[paste0(m, "_sd")]] <- stats::sd(g[[m]])
        }
        out
    }))
}

side_by_side <- function(parts, settings, metrics, sd, table) {
    out <- NULL
    for (name in names(parts)) {
        rows <- keep(parts[[name]], parts[[name]]$setting %in% settings)
        s <- setting_summary(rows, metrics, sd)
        missing <- setdiff(settings, s$setting)
        if (length(missing)) stop(table, " has no ", name, " results for ", missing[1])
        names(s)[-(1:3)] <- paste0(names(s)[-(1:3)], "_", name)
        if (is.null(out)) {
            out <- s
            next
        }
        both <- merge(out, s, by = "setting", suffixes = c("", ".new"))
        differ <- both$reps != both$reps.new
        if (any(differ)) {
            stop(table, ": ", name, " covers other replications of ", both$setting[differ][1])
        }
        out <- both[setdiff(names(both), c("replications.new", "reps.new"))]
    }
    out
}

value_columns <- function(metrics, rules, sd, by_metric) {
    stats <- if (sd) c("mean", "sd") else "mean"
    grid <- if (by_metric) {
        expand.grid(stat = stats, rule = rules, metric = metrics, stringsAsFactors = FALSE)
    } else {
        expand.grid(stat = stats, metric = metrics, rule = rules, stringsAsFactors = FALSE)
    }
    paste0(grid$metric, "_", grid$stat, "_", grid$rule)
}

with_design <- function(x, rows) {
    design <- unique(rows[c("setting", "M", "family", "degree", "dims", "n")])
    x <- merge(design, x, by = "setting")
    size <- vapply(strsplit(x$dims, "x", fixed = TRUE), function(p) prod(as.numeric(p)), 0)
    x <- x[order(match(x$family, family_order), x$degree, x$M, size, x$dims, x$n), ]
    rownames(x) <- NULL
    x[setdiff(names(x), "reps")]
}

table_b4 <- function(p) {
    sel <- p$selected
    stars <- keep(sel, sel$experiment == "stars" & sel$threshold == 0.05 & sel$M == 2)
    parts <- list(scale = path_rows(sel, "scale", 0.5),
                  relaxation = path_rows(sel, "relaxation", 0.5), cv = cv_rows(sel),
                  stars = stars,
                  stars_01 = keep(sel, sel$experiment == "stars" &
                                      sel$threshold == 0.1 & sel$M == 2))
    x <- with_design(side_by_side(parts, unique(stars$setting), c("fdp", "tpr"), TRUE,
                                  "Table B.4"), sel)
    x[c(names(x)[1:7], value_columns(c("fdp", "tpr"), names(parts), TRUE, TRUE))]
}

table_b3 <- function(p) {
    paths <- p$paths
    if (is.null(paths)) stop("Table B.3 needs the candidate recovery of the paths")
    relaxation <- path_rows(p$selected, "relaxation", 0.5)
    r <- relaxation[match(paste(paths$experiment, paths$id),
                          paste(relaxation$experiment, relaxation$id)), ]
    if (anyNA(r$k)) stop("Table B.3: a path without a Relaxation selection at gamma 0.5")
    selected <- r$fp == 0 & r$fn == 0
    groups <- list(any = paths$exact_candidates > 0, adjacent = paths$adjacent_exact,
                   single = paths$exact_candidates == 1,
                   separated = paths$exact_candidates > 1 & !paths$adjacent_exact)
    do.call(rbind, lapply(names(groups), function(g) {
        do.call(rbind, lapply(sort(unique(paths$M)), function(m) {
            in_group <- groups[[g]] & paths$M == m
            data.frame(group = g, M = m, paths = sum(in_group),
                       selected = sum(selected[in_group]),
                       selected_percent = percent(sum(selected[in_group]), sum(in_group)),
                       stringsAsFactors = FALSE)
        }))
    }))
}

table_b1 <- function(p) {
    x <- keep(p$selected, p$selected$experiment == "matrix")
    keys <- unique(x[c("criterion", "gamma")])
    keys <- keys[order(match(keys$criterion, criterion_order), keys$gamma), ]
    do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
        d <- keep(x, x$criterion == keys$criterion[i] & x$gamma == keys$gamma[i])
        small <- d$n == 10
        large <- d$n >= 30
        interior <- sum(d$k > 1 & d$k < d$grid_points)
        data.frame(criterion = keys$criterion[i], gamma = keys$gamma[i], datasets = nrow(d),
                   interior = interior, interior_percent = percent(interior, nrow(d)),
                   datasets_n10 = sum(small), nonempty_n10 = sum(!d$empty[small]),
                   nonempty_n10_percent = percent(sum(!d$empty[small]), sum(small)),
                   datasets_n30plus = sum(large), nonempty_n30plus = sum(!d$empty[large]),
                   nonempty_n30plus_percent = percent(sum(!d$empty[large]), sum(large)),
                   stringsAsFactors = FALSE)
    }))
}

table_b2 <- function(p) {
    x <- keep(p$selected, p$selected$experiment == "units" & p$selected$gamma == 0.5)
    scale <- keep(x, x$criterion == "scale" & x$constant == 1)
    keys <- unique(x[c("criterion", "constant")])
    keys <- keys[order(match(keys$criterion, criterion_order), keys$constant), ]
    do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
        d <- keep(x, x$criterion == keys$criterion[i] & x$constant == keys$constant[i])
        largest <- sum(d$k == d$grid_points)
        smallest <- sum(d$k == 1)
        interior <- nrow(d) - largest - smallest
        same <- if (keys$criterion[i] == "ebic") {
            sum(d$k == scale$k[match(d$id, scale$id)])
        } else {
            NA_real_
        }
        data.frame(criterion = keys$criterion[i], constant = keys$constant[i],
                   datasets = nrow(d), largest = largest,
                   largest_percent = percent(largest, nrow(d)), interior = interior,
                   interior_percent = percent(interior, nrow(d)), smallest = smallest,
                   smallest_percent = percent(smallest, nrow(d)), empty = sum(d$empty),
                   empty_percent = percent(sum(d$empty), nrow(d)), same_as_scale = same,
                   same_as_scale_percent = if (is.na(same)) NA_real_ else percent(same, nrow(d)),
                   mcc_mean = mean(d$mcc), mcc_sd = stats::sd(d$mcc),
                   stringsAsFactors = FALSE)
    }))
}

table_b5 <- function(p) {
    pilot <- p$pilot_candidates
    if (is.null(pilot)) stop("Table B.5 needs the candidates of the Refit pilot")
    unknown <- setdiff(pilot$refit_status, refit_statuses)
    if (length(unknown)) stop("Table B.5: unknown Refit status ", unknown[1])
    sel <- keep(p$selected, p$selected$experiment == "refit_pilot")
    rows <- lapply(split(pilot, pilot$setting), function(g) {
        ids <- unique(g$id)
        relaxation <- keep(sel, sel$setting == g$setting[1] & sel$criterion == "relaxation")
        refit <- keep(sel, sel$setting == g$setting[1] & sel$criterion == "refit")
        if (!setequal(relaxation$id, ids) || !setequal(refit$id, ids) ||
                nrow(relaxation) != length(ids) || nrow(refit) != length(ids)) {
            stop("Table B.5: the candidates and selections of ", g$setting[1],
                 " cover different datasets")
        }
        refit <- refit[match(relaxation$id, refit$id), ]
        both <- !is.na(relaxation$k) & !is.na(refit$k)
        status <- vapply(refit_statuses, function(s) sum(g$refit_status == s), 0)
        data.frame(setting = g$setting[1], replications = length(ids), candidates = nrow(g),
                   refit_ok = status[["ok"]], refit_timeout = status[["timeout"]],
                   refit_error = status[["error"]],
                   refit_nonconverged = status[["nonconverged"]],
                   refit_nonexistent = status[["nonexistent"]],
                   relaxation_not_scored = nrow(g) - sum(g$relaxation_status == "ok"),
                   both_selected = sum(both),
                   same_candidate = sum(relaxation$k[both] == refit$k[both]),
                   mcc_mean_relaxation = mean(relaxation$mcc[both]),
                   mcc_sd_relaxation = stats::sd(relaxation$mcc[both]),
                   mcc_mean_refit = mean(refit$mcc[both]),
                   mcc_sd_refit = stats::sd(refit$mcc[both]), stringsAsFactors = FALSE)
    })
    with_design(do.call(rbind, rows), pilot)
}

support_table <- function(p, tensor, table) {
    sel <- p$selected
    settings <- unique(sel$setting[sel$experiment == if (tensor) "tensor" else "matrix"])
    parts <- list(scale = path_rows(sel, "scale", 0.5),
                  relaxation = path_rows(sel, "relaxation", 0.5), cv = cv_rows(sel))
    metrics <- c("mcc", "fdp", "tpr")
    x <- with_design(side_by_side(parts, settings, metrics, TRUE, table), sel)
    x[c(names(x)[1:7], value_columns(metrics, names(parts), TRUE, TRUE))]
}
table_b6 <- function(p) support_table(p, FALSE, "Table B.6")
table_1 <- function(p) support_table(p, TRUE, "Table 1")

table_b7 <- function(p) {
    x <- keep(p$selected, p$selected$experiment == "matrix")
    gammas <- sort(unique(x$gamma))
    parts <- list()
    for (g in gammas) {
        for (criterion in c("scale", "relaxation")) {
            parts[[paste0(criterion, "_gamma", g)]] <- path_rows(x, criterion, g)
        }
    }
    y <- with_design(side_by_side(parts, unique(x$setting), "mcc", TRUE, "Table B.7"), x)
    y[c(names(y)[1:7], value_columns("mcc", names(parts), TRUE, FALSE))]
}

make_tables <- function(per_dataset) {
    present <- unique(per_dataset$selected$experiment)
    tables <- list()
    skipped <- list()
    for (i in seq_len(nrow(table_notes))) {
        name <- table_notes$table[i]
        missing <- setdiff(strsplit(table_notes$needs[i], ",", fixed = TRUE)[[1]], present)
        if (length(missing)) {
            skipped[[name]] <- paste("no results of", paste(missing, collapse = ", "))
            next
        }
        tables[[name]] <- get(table_notes$builder[i])(per_dataset)
    }
    list(tables = tables, skipped = skipped)
}
