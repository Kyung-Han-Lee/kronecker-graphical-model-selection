golden <- readRDS("fixtures/golden.rds")
m2 <- golden$m2
path <- kgm_path(m2$Y, m2$eta)
supports <- lapply(path$fits[[4]], function(A) {
    S <- abs(A) > 1e-6
    diag(S) <- FALSE
    S
})
rf <- kgm_refit(m2$Y, supports)
expect_true(isTRUE(attr(rf, "converged")))
for (m in 1:2) {
    off <- !supports[[m]] & upper.tri(supports[[m]])
    expect_true(all(rf[[m]][off] == 0))
}
expect_equal(sqrt(sum(rf[[1]]^2)), 1, tolerance = 1e-12)
got <- kgm_score(path, m2$Y, "refit", gamma = 0)
want <- m2$scores[m2$scores$criterion == "refit" & m2$scores$gamma == 0, ]
expect_equal(got$loglik, want$loglik[order(want$k)], tolerance = 1e-10)

short <- list(max_iter = 1L, tol = 1e-5)
nc <- kgm_score(path, m2$Y, "refit", gamma = 0, refit_control = short)
expect_true(all(nc$status == "nonconverged"))
expect_true(all(is.na(nc$ebic)))
expect_error(kgm_select(path, m2$Y, "refit", gamma = 0, refit_control = short),
             pattern = "no candidate")

expect_true(any(grepl("refit produced a non-PD factor",
                      deparse(kgmselect:::.refit_support_mle), fixed = TRUE)))
expect_true(kgmselect:::.is_refit_failure(
    simpleError("refit produced a non-PD factor at mode 1 (min eig -1.000e-03)")))
expect_false(kgmselect:::.is_refit_failure(simpleError("subscript out of bounds")))

bad <- path
bad$fits[[3]][[1]] <- diag(c(1, -1, 1, 1, 1, 1)) / sqrt(6)
good <- path
good$fits[[3]][[1]] <- diag(6) / sqrt(6)
rb <- kgm_score(bad, m2$Y, "refit", gamma = 0.5)
rg <- kgm_score(good, m2$Y, "refit", gamma = 0.5)
expect_identical(rb$status[3], "ok")
expect_identical(rb$loglik[3], rg$loglik[3])
