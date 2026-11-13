source("helper_lib.R")
g <- readRDS(system.file("tinytest", "fixtures", "golden.rds", package = "kgmselect"))$m2
path <- kgm_path(g$Y, g$eta)
scores <- score_path(path, g$Y, c("ebic", "scale"), c(0, 0.5), "numeric")
expect_identical(unique(scores$criterion), c("ebic", "scale"))
for (crit in c("ebic", "scale")) {
    for (gm in c(0, 0.5)) {
        s <- scores[scores$criterion == crit & scores$gamma == gm, ]
        expect_identical(select_index(s), kgm_select(path, g$Y, crit, gamma = gm,
                                                    tol = 0)$index,
                         info = paste(crit, gm))
    }
}

tab <- selection_tables(scores, path, g$truth)
expect_equal(nrow(tab$selected), 4L)
expect_equal(nrow(tab$selected_modes), 8L)
first <- tab$selected[1, ]
m <- kgm_metrics(g$truth, path$fits[[first$k]], tol = 0)
expect_equal(first$mcc, m$mcc)
expect_equal(first$tp, m$tp)
expect_equal(first$edges, sum(vapply(path$fits[[first$k]], function(A) {
    sum(A[upper.tri(A)] != 0)
}, 0)))
expect_identical(tab$selected$empty, tab$selected$edges == 0L)
expect_equal(tab$selected_modes$tp[tab$selected_modes$criterion == first$criterion &
                                      tab$selected_modes$gamma == first$gamma],
             m$per_mode$tp)

s <- scores[scores$criterion == "scale" & scores$gamma == 0.5, ]
k <- select_index(s)
s$status[s$k == k] <- "not_pd"
expect_true(select_index(s) != k)
s$status <- "not_pd"
expect_error(select_index(s), "no candidate")

ft <- fit_tables(data.frame(criterion = first$criterion, gamma = first$gamma),
                 first$k, first$eta, path$fits[[first$k]], g$truth)
row1 <- tab$selected[1, ]
rownames(row1) <- NULL
expect_equal(ft$selected, row1)
modes1 <- tab$selected_modes[tab$selected_modes$criterion == first$criterion &
                                 tab$selected_modes$gamma == first$gamma, ]
rownames(modes1) <- NULL
expect_equal(ft$selected_modes, modes1)

sr <- support_records(path)
expect_identical(sr$k, seq_along(path$fits))
expect_equal(sr$edges, s$edges)
A5 <- path$fits[[5]][[1]]
zero <- which(A5 == 0 & upper.tri(A5), arr.ind = TRUE)[1, ]
tweak <- path
tweak$fits[[5]][[1]][zero[1], zero[2]] <- tweak$fits[[5]][[1]][zero[2], zero[1]] <- 5e-7
positive_tol <- support_records(path, tol = 1e-6)
sr2 <- support_records(tweak, tol = 1e-6)
expect_identical(sr2$tiny_edges[5], positive_tol$tiny_edges[5] + 1L)
expect_identical(sr2$max_tiny[5], 5e-7)
expect_identical(sr2$pattern_hash, positive_tol$pattern_hash)
expect_false(identical(sr2$exact_hash[5], positive_tol$exact_hash[5]))
mask <- lapply(path$fits[[5]], function(A) {
    keep <- abs(A) > 1e-6
    diag(keep) <- FALSE
    keep
})
expect_identical(positive_tol$pattern_hash[5], digest::digest(mask, algo = "sha256"))

rec <- candidate_recovery(path, g$truth)
expect_identical(rec$k, seq_along(path$fits))
expect_equal(unlist(rec[first$k, c("tp", "fp", "fn", "tn")]),
             unlist(first[, c("tp", "fp", "fn", "tn")]))
exact <- path
exact$fits[[1]] <- g$truth
expect_equal(unlist(candidate_recovery(exact, g$truth)[1, c("fp", "fn")]), c(fp = 0, fn = 0))
