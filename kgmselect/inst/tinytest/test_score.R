source("helper_reference.R")
golden <- readRDS("fixtures/golden.rds")
ref <- reference_for_installed_huge()
m2 <- golden$m2
path <- kgm_path(m2$Y, m2$eta)
n <- dim(m2$Y)[3]
p <- prod(dim(m2$Y)[1:2])

expect_equal(m2$scores$ebic,
             -2 * m2$scores$loglik + m2$scores$df * (log(n) + 4 * m2$scores$gamma * log(p)),
             tolerance = 1e-12)

scores <- list()
for (crit in c("ebic", "scale", "relaxation", "refit")) {
    got <- kgm_score(path, m2$Y, criterion = crit, gamma = c(0, 0.5))
    scores[[crit]] <- got
    want <- m2$scores[m2$scores$criterion == crit, ]
    want <- want[order(want$gamma, want$k), ]
    expect_equal(got$k, want$k, info = crit)
    expect_equal(got$df, want$df, info = crit)
    expect_false("q" %in% names(got), info = crit)
    expect_true(all(got$status == "ok"), info = crit)

    for (g in c(0, 0.5)) {
        a <- got[got$gamma == g, ]
        b <- want[want$gamma == g, ]
        k <- a$k[which.min(a$ebic)]
        expect_equal(k, b$k[which.min(b$ebic)], info = crit)
        expect_identical(support_of(path$fits[[k]]), support_of(m2$fits[[k]]),
                         info = crit)
    }
}

eb <- kgm_select(path, m2$Y, "ebic", gamma = 0.5)
expect_equal(eb$index, length(m2$eta))
expect_equal(sum(vapply(eb$supports, sum, 0)), 0)

sc <- kgm_score(path, m2$Y, "scale", gamma = 0.5)
expect_equal(sc$ebic, -2 * sc$loglik + (sc$edges + 1) * (log(n) + 2 * log(p)),
             tolerance = 1e-10)
rx <- kgm_score(path, m2$Y, "relaxation", gamma = 0.5)
expect_equal(rx$ebic, -2 * rx$loglik + (rx$edges + 1) * (log(n) + 2 * log(p)),
             tolerance = 1e-10)
expect_true(all(rx$loglik >= sc$loglik - 1e-8 * abs(sc$loglik)))

sel <- kgm_select(path, m2$Y, "relaxation", gamma = 0.5)
want <- m2$scores[m2$scores$criterion == "relaxation" & m2$scores$gamma == 0.5, ]
want <- want[order(want$k), ]
expect_equal(sel$index, want$k[which.min(want$ebic)])
expect_equal(sel$eta, m2$eta[sel$index])
expect_equal(vapply(sel$factors, function(A) sqrt(sum(A^2)), 0), c(1, 1),
             tolerance = 1e-12)
expect_true(sel$tau > 0)
fit <- path$fits[[sel$index]]
for (m in 1:2) {
    iu <- upper.tri(fit[[m]])
    expect_identical(sel$factors[[m]][iu] != 0, fit[[m]][iu] != 0)
}

bad <- path
bad$fits[[3]][[1]] <- diag(c(1, -1, 1, 1, 1, 1)) / sqrt(6)
s_bad <- kgm_score(bad, m2$Y, "scale", gamma = 0.5)
expect_equal(s_bad$status[3], "not_pd")
expect_true(is.na(s_bad$ebic[3]))
expect_true(all(s_bad$status[-3] == "ok"))
expect_true(kgm_select(bad, m2$Y, "scale", gamma = 0.5)$index != 3L)

stingy <- list(maxit = 0L, factr = 1e10)
expect_error(kgm_score(path, m2$Y, "relaxation", gamma = 0.5, optim_control = stingy),
             pattern = "all joint optimizations failed")
expect_error(kgm_select(path, m2$Y, "relaxation", gamma = 0.5, optim_control = stingy),
             pattern = "candidate 2")

m3 <- golden$m3
path3 <- kgm_path(m3$Y, m3$eta)
num <- kgm_score(path3, m3$Y, "relaxation", gamma = 0.5, gradient = "numeric")
ana <- kgm_score(path3, m3$Y, "relaxation", gamma = 0.5, gradient = "analytic")
expect_equal(num$loglik, ana$loglik, tolerance = 1e-6)
expect_error(kgm_select(path, m2$Y, "relaxation", gamma = c(0, 0.5)))

if (is.null(ref)) {
    exit_file(sprintf("no reference values for huge %s", packageVersion("huge")))
}
for (crit in names(scores)) {
    want <- ref$m2$scores[ref$m2$scores$criterion == crit, ]
    want <- want[order(want$gamma, want$k), ]
    expect_equal(scores[[crit]]$loglik, want$loglik, tolerance = 1e-10, info = crit)
    expect_equal(scores[[crit]]$ebic, want$ebic, tolerance = 1e-10, info = crit)
}
wn <- ref$m3$scores_numeric[ref$m3$scores_numeric$criterion == "relaxation", ]
wa <- ref$m3$scores_analytic[ref$m3$scores_analytic$criterion == "relaxation", ]
expect_equal(num$loglik, wn$loglik[order(wn$k)], tolerance = 1e-10)
expect_equal(ana$loglik, wa$loglik[order(wa$k)], tolerance = 1e-10)
