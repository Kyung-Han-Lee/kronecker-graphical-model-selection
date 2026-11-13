.n_pass <- 0L
.n_fail <- 0L
check <- function(name, ok) {
    if (isTRUE(ok)) {
        .n_pass <<- .n_pass + 1L
    } else {
        .n_fail <<- .n_fail + 1L
        cat(sprintf("FAIL  %s\n", name))
    }
}

idx <- bootstrap_indices(n = 5L, B = 3L, seed = 1L)
check("재추출 수가 B다", length(idx) == 3L)
check("각 재추출의 크기가 n이다", all(lengths(idx) == 5L))
check("색인이 1..n 안에 있다", all(unlist(idx) >= 1L & unlist(idx) <= 5L))
check("같은 시드는 같은 색인을 준다", identical(idx, bootstrap_indices(5L, 3L, 1L)))
check("다른 시드는 다른 색인을 준다", !identical(idx, bootstrap_indices(5L, 3L, 2L)))

check("뒤대상피질 쌍은 AAL 35·36행이다", identical(BOOT_PAIRS$posterior_cingulate, c(35L, 36L)))
check("쐐기앞소엽 쌍은 AAL 67·68행이다", identical(BOOT_PAIRS$precuneus, c(67L, 68L)))

Omega <- matrix(c(
    2, -0.5, 0, 0,
    -0.5, 2, 0.3, 0,
    0, 0.3, 1.5, 0,
    0, 0, 0, 1
), 4, 4, byrow = TRUE)
pcs <- pair_partial_correlations(Omega, list(a = c(1L, 2L), b = c(3L, 4L), c = c(2L, 3L)))
check("부분상관 이름 순서를 지킨다", identical(pcs$pair, c("a", "b", "c")))
check("부분상관 a = 0.25", isTRUE(all.equal(pcs$pc[1], 0.25)))
check("부분상관 b = 0", isTRUE(all.equal(pcs$pc[2], 0)))
check("부분상관 c = -0.3/sqrt(3)", isTRUE(all.equal(pcs$pc[3], -0.3 / sqrt(3))))
check("0 원소는 선택되지 않은 엣지다", identical(pcs$selected, c(TRUE, FALSE, TRUE)))

synth <- do.call(rbind, lapply(1:6, function(i) {
    data.frame(
        group = "G", arm = "A", eta_index = i, eta = 10^(i - 4),
        gamma = GAMMA_GRID,
        edges_roi = i, edges_time = 2L * i, df = i,
        ebic_plug = 100 - i, ebic_prof = (i - 2)^2,
        ebic_joint = if (i == 4L) NA_real_ else c(9, 5, 7, 0, 1, 3)[i] + GAMMA_GRID,
        ebic_joint_uncounted = 0,
        stringsAsFactors = FALSE
    )
}))
check("결합 선택은 유한한 최소점이다", identical(select_joint_index(synth, 0.5), 5L))
v <- verdict_table(synth)
check("verdict_table과 같은 격자점을 고른다", identical(
    select_joint_index(synth, 0.5),
    as.integer(v$best_index[v$gamma == 0.5 & v$criterion == "joint"])
))
check("모든 점수가 결측이면 NA다", is.na(select_joint_index(
    transform(synth, ebic_joint = NA_real_), 0.5
)))

check("백분위 구간", isTRUE(all.equal(
    percentile_interval(c(1:100, NA)), c(3.475, 97.525)
)))

draws <- data.frame(
    group = rep(c("CN", "Dementia"), each = 6),
    b = rep(0:5, 2),
    pair = "p",
    pc = c(0.60, 0.60, 0.62, 0.58, 0.61, NA, 0.50, 0.50, 0.52, 0.48, 0.49, 0.47),
    selected = c(TRUE, TRUE, TRUE, FALSE, TRUE, NA, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
    ok = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, rep(TRUE, 6)),
    stringsAsFactors = FALSE
)
s <- summarise_draws(draws)
cn <- s[s$group == "CN", ]
dm <- s[s$group == "Dementia", ]
df <- s[s$group == "Dementia-CN", ]
check("요약은 군 두 개와 차이 하나다", nrow(s) == 3L)
check("CN 추정치는 b = 0 값이다", isTRUE(all.equal(cn$estimate, 0.60)))
check("CN 구간은 성공한 b >= 1 의 백분위다", isTRUE(all.equal(
    c(cn$lower, cn$upper), percentile_interval(c(0.60, 0.62, 0.58, 0.61))
)))
check("실패한 재추출은 개수에서 빠진다", cn$n_boot == 4L)
check("선택 비율", isTRUE(all.equal(cn$selected_freq, 0.75)))
check("치매군 추정치", isTRUE(all.equal(dm$estimate, 0.50)))
check("차이 추정치는 b = 0 차이다", isTRUE(all.equal(df$estimate, -0.10)))
check("차이 구간은 같은 b 끼리의 차이로 만든다", isTRUE(all.equal(
    c(df$lower, df$upper),
    percentile_interval(c(0.50, 0.52, 0.48, 0.49) - c(0.60, 0.62, 0.58, 0.61))
)))
check("차이는 두 군이 모두 성공한 b 만 쓴다", df$n_boot == 4L)

set.seed(11)
x_small <- array(stats::rnorm(5 * 4 * 30), dim = c(5, 4, 30))
grid_small <- make_eta_grid(6L)
path_small <- score_path(x_small, grid_small)
check("경로 행 수는 격자점 × gamma 수다", nrow(path_small) == 6L * length(GAMMA_GRID))
check("격자 색인은 1..6이다", identical(sort(unique(path_small$eta_index)), 1:6))
check("결합 점수 열이 있다", "ebic_joint" %in% names(path_small))
check("작은 예에서 결합 점수가 유한하다", all(is.finite(path_small$ebic_joint)))
check("선택 함수가 경로를 받는다", select_joint_index(path_small, 0.5) %in% 1:6)

set.seed(12)
n_syn <- 30L
Y_syn <- array(stats::rnorm(n_syn * 5 * 4), dim = c(n_syn, 5, 4))
X_syn <- cbind(Age = stats::rnorm(n_syn, 70, 5), Sex = rep(0:1, length.out = n_syn), APOE4 = rep(0:2, length.out = n_syn))
one <- run_resample(Y_syn, X_syn, seq_len(n_syn), grid_small, pairs = list(a = c(1L, 2L)))
check("재추출 하나가 성공한다", isTRUE(one$ok))
check("선택 격자점이 격자 안에 있다", one$best_index %in% 1:6)
check("부분상관 행이 쌍 수만큼이다", nrow(one$pcs) == 1L)
check("모드별 엣지 수가 둘이다", length(one$edges) == 2L)
check("재추출 없이 돌린 선택은 같은 경로의 선택과 같다", identical(
    one$best_index,
    select_joint_index(score_path(aperm(residualize(Y_syn, X_syn), c(2, 3, 1)), grid_small), 0.5)
))
bad_X <- X_syn
bad_X[, "Sex"] <- 1
bad <- run_resample(Y_syn, bad_X, seq_len(n_syn), grid_small, pairs = list(a = c(1L, 2L)))
check("공변량이 상수가 되면 실패를 돌려준다", identical(bad$ok, FALSE))
check("실패 사유를 남긴다", is.character(bad$error) && nzchar(bad$error))

check("정상 재추출은 실패 격자가 없다", identical(one$failed_grid, ""))
check("정상 재추출은 반환 후보의 q 를 남긴다", is.integer(one$q_u) && length(one$q_u) == 1L && !is.na(one$q_u))

score_fit_orig <- score_fit
fail_lambda <- lambda_from_eta(c(5L, 4L), n_syn, grid_small[3])
score_fit <- function(x_val, fit, n, prev_fit = NULL, lambda, prev_lambda = NULL) {
    if (isTRUE(all.equal(lambda, fail_lambda))) stop("injected scoring failure")
    score_fit_orig(x_val, fit, n, prev_fit = prev_fit, lambda = lambda, prev_lambda = prev_lambda)
}
broken <- run_resample(Y_syn, X_syn, seq_len(n_syn), grid_small, pairs = list(a = c(1L, 2L)))
score_fit <- score_fit_orig
check("격자점 하나가 실패하면 재추출을 실패로 기록한다", identical(broken$ok, FALSE))
check("실패한 격자 색인을 남긴다", identical(broken$failed_grid, "3"))
check("실패 사유에 격자 색인이 있다", grepl("3", broken$error, fixed = TRUE))
check("실패한 재추출은 부분상관을 내지 않는다", all(is.na(broken$pcs$pc)))

selected_fit_at_orig <- selected_fit_at
selected_fit_at <- function(...) {
    s <- selected_fit_at_orig(...)
    s$q_u <- s$q_u + 1L
    s
}
mismatch <- run_resample(Y_syn, X_syn, seq_len(n_syn), grid_small, pairs = list(a = c(1L, 2L)))
selected_fit_at <- selected_fit_at_orig
check("다시 계산한 후보가 채점한 후보와 다르면 실패로 기록한다", identical(mismatch$ok, FALSE))
check("불일치 사유를 남긴다", grepl("scored", mismatch$error, fixed = TRUE))

hp <- homotopic_pairs(c("A_L", "A_R", "B_L", "Vermis", "B_R", "C_R"))
check("좌우 대응 쌍은 _L·_R 짝만 묶는다", identical(hp, list(A = c(1L, 2L), B = c(3L, 5L))))
labels <- rep("Vermis", 116L)
for (j in seq_len(54L)) labels[2L * j - c(1L, 0L)] <- paste0("Region", j, c("_L", "_R"))
labels[c(35L, 36L)] <- c("Cingulum_Post_L", "Cingulum_Post_R")
labels[c(67L, 68L)] <- c("Precuneus_L", "Precuneus_R")
label_file <- tempfile(fileext = ".txt")
write.table(data.frame(index = seq_along(labels), labels = labels), label_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
aal <- homotopic_pairs(read_aal_labels(label_file))
unlink(label_file)
check("AAL V4 의 좌우 대응 쌍은 54개다", length(aal) == 54L)
check("뒤대상피질 쌍은 35·36행이다", identical(aal$Cingulum_Post, c(35L, 36L)))
check("쐐기앞소엽 쌍은 67·68행이다", identical(aal$Precuneus, c(67L, 68L)))

check("refit 하지 않으면 refit 부분상관은 결측이다", all(is.na(one$pcs$pc_refit)))
two <- run_resample(Y_syn, X_syn, seq_len(n_syn), grid_small, pairs = list(a = c(1L, 2L), b = c(3L, 4L)), refit = TRUE)
check("refit 재추출이 성공한다", isTRUE(two$ok))
check("refit 수렴 여부를 남긴다", is.logical(two$refit_converged) && length(two$refit_converged) == 1L)
x_syn <- aperm(residualize(Y_syn, X_syn), c(2, 3, 1))
sel_syn <- selected_fit_at(x_syn, grid_small, two$best_index, "joint", n_sweep = N_SWEEP)
ref_syn <- refit_support_mle(x_syn, lapply(sel_syn$path_fit, function(M) abs(M) > SUPPORT_TOL))
check("refit 부분상관은 선택 support 의 비벌점 최대우도 인자에서 나온다", isTRUE(all.equal(
    two$pcs$pc_refit, pair_partial_correlations(ref_syn[[1]], list(a = c(1L, 2L), b = c(3L, 4L)))$pc
)))
check("relaxation 부분상관은 refit 여부와 무관하다", isTRUE(all.equal(
    two$pcs$pc, pair_partial_correlations(sel_syn$Omega_list[[1]], list(a = c(1L, 2L), b = c(3L, 4L)))$pc
)))
check("선택되지 않은 쌍의 refit 부분상관은 0이다", all(two$pcs$selected | two$pcs$pc_refit == 0))

fam <- expand.grid(b = 0:4, pair = c("p1", "p2"), group = c("CN", "Dementia"), stringsAsFactors = FALSE)
fam$pc <- 0.6
fam$pc[fam$group == "Dementia" & fam$pair == "p1"] <- 0.5
d2 <- c(0, -0.1, 0.1, -0.2, 0.2)
fam$pc[fam$group == "Dementia" & fam$pair == "p2"] <- 0.6 + d2
fam$selected <- TRUE
fam$ok <- TRUE
fs <- summarise_family(fam, value = "pc", contrasts = list(c("Dementia", "CN"), c("MCI", "CN")))
dif <- fs[fs$kind == "difference", ]
check("요약은 있는 군의 대비만 낸다", identical(unique(dif$contrast), "Dementia-CN"))
check("쌍마다 차이 행이 하나다", identical(dif$pair, c("p1", "p2")))
check("차이 추정치는 b = 0 의 차이다", isTRUE(all.equal(dif$estimate, c(-0.1, 0))))
check("붓스트랩 p 값", isTRUE(all.equal(dif$p_boot, c(0.4, 1))))
check("BH q 값", isTRUE(all.equal(dif$q_bh, stats::p.adjust(c(0.4, 1), "BH"))))
check("구간이 0을 벗어나는지 표시한다", identical(dif$excludes_zero, c(TRUE, FALSE)))
grp <- fs[fs$kind == "group", ]
check("군별 행은 군 × 쌍이다", nrow(grp) == 4L)

cat(sprintf("\n통과 %d / 실패 %d\n", .n_pass, .n_fail))
if (.n_fail > 0L) quit(status = 1L)
