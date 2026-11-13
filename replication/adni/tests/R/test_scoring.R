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

synth_path <- function(n_eta = 8L, prof_argmin = 4L, na_at = integer(0),
                       group = "G", arm = "control") {
    k <- seq_len(n_eta)
    plug <- 100 - k
    prof <- (k - prof_argmin)^2
    joint <- (k - (prof_argmin + 1L))^2
    rows <- lapply(k, function(i) {
        na_here <- i %in% na_at
        data.frame(
            group = group, arm = arm, eta_index = i, eta = 10^(-4 + i),
            gamma = GAMMA_GRID,
            edges_roi = if (na_here) NA_integer_ else as.integer(i),
            edges_time = if (na_here) NA_integer_ else as.integer(10 * i),
            df = if (na_here) NA_integer_ else as.integer(i),
            ebic_plug = if (na_here) NA_real_ else plug[i],
            ebic_prof = if (na_here) NA_real_ else prof[i],
            ebic_joint = if (na_here) NA_real_ else joint[i],
            ebic_joint_uncounted = if (na_here) NA_real_ else joint[i] - 1,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

pick <- function(v, crit) {
    v[v$criterion == crit & v$gamma == GAMMA_REPORT, ]
}

v <- verdict_table(synth_path(n_eta = 8L, prof_argmin = 4L))
check("판정 표가 기준 넷 x gamma 넷 = 16 행", nrow(v) == 16L)
check("plug 는 오른쪽 끝을 고른다", pick(v, "plug")$best_index == 8L)
check("plug 의 끝 선택은 내부가 아니다", isFALSE(pick(v, "plug")$interior))
check("prof 는 내부 최소를 고른다", pick(v, "prof")$best_index == 4L)
check("prof 의 내부 선택이 내부로 기록된다", isTRUE(pick(v, "prof")$interior))
check("joint 는 joint 점수를 고른다", pick(v, "joint")$best_index == 5L)
check(
    "joint uncounted 민감도 판정도 남긴다",
    pick(v, "joint_uncounted")$best_index == 5L
)
check(
    "선택점의 엣지 수가 그 점의 것이다",
    pick(v, "prof")$edges_roi == 4L && pick(v, "prof")$edges_time == 40L
)

v_na <- verdict_table(
    synth_path(n_eta = 8L, prof_argmin = 3L, na_at = c(1L, 2L))
)
check("실패 점을 뺀 유효 점 수를 센다", pick(v_na, "prof")$n_valid == 6L)
check("격자 길이는 그대로 기록된다", pick(v_na, "prof")$n_eta == 8L)
check("prof 가 유효 구간의 왼쪽 끝을 골랐다", pick(v_na, "prof")$best_index == 3L)
check(
    "유효 구간의 끝은 내부가 아니다 (실패 점 옆이라도)",
    isFALSE(pick(v_na, "prof")$interior)
)
v_na2 <- verdict_table(
    synth_path(n_eta = 8L, prof_argmin = 4L, na_at = c(1L, 2L))
)
check("유효 구간 안쪽이면 내부다", isTRUE(pick(v_na2, "prof")$interior))

v_dead <- verdict_table(synth_path(n_eta = 4L, na_at = 1:4))
check("전부 실패하면 판정 행이 없다", is.null(v_dead) || nrow(v_dead) == 0L)

set.seed(42)
p1 <- 8L
p2 <- 6L
n_obs <- 40L
x_val <- array(stats::rnorm(p1 * p2 * n_obs), dim = c(p1, p2, n_obs))
lambda <- lambda_from_eta(c(p1, p2), n_obs, eta = 0.5)
fit <- lapply(fit_tggm(x_val, lambda, n_sweep = N_SWEEP), function(M) M)
scored <- score_fit(x_val, fit, n_obs, lambda = lambda)

check("채점 표가 gamma 격자 길이만큼 나온다", nrow(scored) == length(GAMMA_GRID))
for (gi in seq_along(GAMMA_GRID)) {
    g <- GAMMA_GRID[gi]
    for (tc in c(FALSE, TRUE)) {
        direct <- compute_ebic_tlasso(
            x_val, fit,
            gamma = g, tol = SUPPORT_TOL,
            tau_correct = tc, sqrt_method = SQRT_METHOD
        )
        col <- if (tc) "ebic_prof" else "ebic_plug"
        check(
            sprintf(
                "gamma=%g %s 의 산술이 직접 호출과 일치", g,
                if (tc) "profiled" else "plug-in"
            ),
            abs(scored[[col]][gi] - direct$ebic) < 1e-8 * max(1, abs(direct$ebic))
        )
    }
}
check(
    "기록된 df 에 profiled 의 +1 이 안 들어간다",
    scored$df[1] == compute_ebic_tlasso(
        x_val, fit,
        gamma = 0, tol = SUPPORT_TOL, tau_correct = FALSE,
        sqrt_method = SQRT_METHOD
    )$df
)
check(
    "첫 격자점의 joint 점수는 scale-profiled 점수와 같다",
    max(abs(scored$ebic_joint - scored$ebic_prof)) < 1e-8
)
check("첫 격자점에는 relaxation 좌표가 없다",
      all(scored$q_u == 0L) && all(scored$u_roi == 0) && all(scored$u_time == 0))

lambda_next <- lambda_from_eta(c(p1, p2), n_obs, eta = 0.8)
fit_next <- lapply(
    fit_tggm(x_val, lambda_next, n_sweep = N_SWEEP), function(M) M
)
scored_next <- score_fit(
    x_val, fit_next, n_obs,
    prev_fit = fit, lambda = lambda_next, prev_lambda = lambda
)
check("인접 격자점의 joint 점수가 유한하다",
      all(is.finite(scored_next$ebic_joint)))
check("인접 격자점에는 relaxation 좌표가 있다", all(scored_next$q_u > 0L))
check(
    "joint 점수에 q_u log(n) 이 들어가지 않는다",
    max(abs(scored_next$ebic_joint -
        (-2 * scored_next$loglik_joint + (scored_next$df + 1L) *
            (log(n_obs) + 4 * scored_next$gamma * log(p1 * p2))))) < 1e-8
)
check(
    "ebic_joint_uncounted 열은 joint 점수와 같다",
    identical(scored_next$ebic_joint_uncounted, scored_next$ebic_joint)
)
check("joint profiling은 현재 support의 edge 수를 바꾸지 않는다",
      all(scored_next$edges_roi == scored_next$edges_roi[1]) &&
          all(scored_next$edges_time == scored_next$edges_time[1]))

selected_joint <- selected_fit_at(
    x_val, eta_grid = c(0.5, 0.8), best_index = 2L,
    criterion = "joint", n_sweep = N_SWEEP
)
fit_next_support <- lapply(fit_next, function(M) M != 0)
joint_support <- lapply(selected_joint$Omega_list, function(M) M != 0)
check(
    "선택점 joint 행렬은 현재 경로점의 exact support를 보존한다",
    identical(fit_next_support, joint_support)
)
fit_next_edges <- vapply(
    fit_next, function(M) sum(abs(M[upper.tri(M)]) > SUPPORT_TOL), integer(1)
)
check(
    "선택점 edge 수는 변환 행렬이 아니라 공급된 경로점에서 고정한다",
    identical(selected_joint$path_edge_counts, fit_next_edges)
)
check(
    "선택점 joint 행렬과 relaxation 좌표가 유한하다",
    all(vapply(selected_joint$Omega_list, function(M) all(is.finite(M)), TRUE)) &&
        length(selected_joint$u) == 2L && all(is.finite(selected_joint$u))
)
check(
    "선택점 joint 인자는 단위 Frobenius norm을 갖는다",
    max(abs(vapply(selected_joint$Omega_list, norm, 0, type = "F") - 1)) < 1e-10
)

na_row <- na_score_row()
check("실패 자리표의 열 이름이 채점 표와 같다",
      identical(names(na_row), names(scored)))
check("실패 자리표의 열 형이 채점 표와 같다",
      identical(vapply(na_row, class, ""), vapply(scored, class, "")))

set.seed(7)
Y <- array(stats::rnorm(20 * 4 * 3), dim = c(20L, 4L, 3L))
X <- cbind(Age = stats::rnorm(20), Sex = rep(0:1, 10), APOE4 = stats::rnorm(20))
resid <- residualize(Y, X)
check("잔차 텐서의 차원이 보존된다", identical(dim(resid), dim(Y)))
check(
    "잔차가 설계행렬에 직교한다",
    max(abs(crossprod(cbind(1, X), matrix(resid, nrow = 20)))) < 1e-9
)
check(
    "공변량이 상수면 랭크 결손으로 멈춘다",
    inherits(try(
        residualize(Y, cbind(a = rep(1, 20), b = rep(1, 20), c = rep(1, 20))),
        silent = TRUE
    ), "try-error")
)

write_manifest <- function(dir, records_json) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(
        sprintf('{"groups": [%s]}', paste(records_json, collapse = ",")),
        file.path(dir, "resample_manifest.json")
    )
}
record_json <- function(group, output, problems = character(0)) {
    sprintf(
        '{"group": "%s", "output": %s, "problems": [%s]}',
        group,
        if (is.null(output)) "null" else sprintf('"%s"', output),
        paste(sprintf('"%s"', problems), collapse = ",")
    )
}
tmp_root <- file.path(tempdir(), "probe_manifest_test")
unlink(tmp_root, recursive = TRUE)
sub <- "combat_tr5"

write_manifest(file.path(tmp_root, sub), record_json("Dementia", "x.npy"))
check("무결 기록이 있으면 통과", is.null(
    check_resample_manifest(tmp_root, sub, "Dementia", required = TRUE)
))
check("다른 군의 기록만 있으면 실패 (부분 재실행 잔재)", inherits(try(
    check_resample_manifest(tmp_root, sub, "CN", required = TRUE),
    silent = TRUE
), "try-error"))
write_manifest(file.path(tmp_root, sub), c(
    record_json("Dementia", "x.npy"), record_json("Dementia", "y.npy")
))
check("같은 군의 기록이 둘이면 실패", inherits(try(
    check_resample_manifest(tmp_root, sub, "Dementia", required = TRUE),
    silent = TRUE
), "try-error"))

write_manifest(
    file.path(tmp_root, sub),
    record_json("Dementia", NULL, problems = "빈 방향 과다")
)
check("problems 가 있으면 실패", inherits(try(
    check_resample_manifest(tmp_root, sub, "Dementia", required = TRUE),
    silent = TRUE
), "try-error"))

write_manifest(file.path(tmp_root, sub), record_json("Dementia", NULL))
check("저장되지 않은 기록(output=null)은 실패", inherits(try(
    check_resample_manifest(tmp_root, sub, "Dementia", required = TRUE),
    silent = TRUE
), "try-error"))

check("manifest 가 아예 없으면 재표집 팔은 실패", inherits(try(
    check_resample_manifest(tmp_root, "combat_absent", "Dementia",
                            required = TRUE),
    silent = TRUE
), "try-error"))
check("manifest 가 없어도 control 팔은 통과", is.null(
    check_resample_manifest(tmp_root, "combat", "Dementia", required = FALSE)
))
unlink(tmp_root, recursive = TRUE)

fake_np <- function(arr) list(load = function(path) arr)

arm_root <- file.path(tempdir(), "probe_load_arm_test")
unlink(arm_root, recursive = TRUE)
dir.create(file.path(arm_root, "combat"), recursive = TRUE)
utils::write.csv(
    data.frame(Age = c(70, 75, 80), Sex = c(0, 1, 0), APOE4 = c(0, 1, 2)),
    file.path(arm_root, "combat", "covars_Dementia.csv"),
    row.names = FALSE
)
dir.create(file.path(arm_root, "combat_tr5"), showWarnings = FALSE)
invisible(file.create(file.path(arm_root, "combat_tr5", "Y_Dementia_combat.npy")))
invisible(file.create(file.path(arm_root, "combat", "Y_Dementia_combat.npy")))
write_manifest(file.path(arm_root, "combat_tr5"), record_json("Dementia", "y.npy"))

good <- array(stats::rnorm(3 * 4 * 113), dim = c(3L, 4L, 113L))
loaded <- load_arm("Dementia", "resampled", arm_root, fake_np(good))
check("기대 차원이면 텐서와 공변량이 나온다",
      identical(dim(loaded$Y), c(3L, 4L, 113L)) && identical(dim(loaded$X), c(3L, 3L)))

wrong <- array(stats::rnorm(3 * 4 * 94), dim = c(3L, 4L, 94L))
check("시간 차원이 기대와 다르면 멈춘다", inherits(try(
    load_arm("Dementia", "resampled", arm_root, fake_np(wrong)), silent = TRUE
), "try-error"))

check("control 팔은 187 을 기대한다", inherits(try(
    load_arm("Dementia", "control", arm_root, fake_np(good)), silent = TRUE
), "try-error"))

check("공변량 행수가 피험자 수와 다르면 멈춘다", inherits(try(
    load_arm("Dementia", "resampled", arm_root,
             fake_np(array(0, dim = c(5L, 4L, 113L)))),
    silent = TRUE
), "try-error"))
unlink(arm_root, recursive = TRUE)

empty_path <- synth_path(n_eta = 6L, prof_argmin = 3L)
empty_path$df[empty_path$eta_index == 6L] <- 0L
v_empty <- verdict_table(empty_path)
check("격자 오른쪽 끝의 빈 그래프가 is_empty 로 표시된다",
      isTRUE(pick(v_empty, "plug")$is_empty))
check("내부를 고른 기준은 빈 그래프가 아니다",
      isFALSE(pick(v_empty, "prof")$is_empty))

check("ADNI support tolerance", identical(SUPPORT_TOL, 1e-6))
check("ADNI sweeps and square root", identical(N_SWEEP, 2L) && identical(SQRT_METHOD, "none"))
check("ADNI eta endpoints and grid length", identical(c(ETA_LO, ETA_HI, N_ETA_DEFAULT), c(1e-4, 100, 40)))

check("eta 격자가 양 끝을 덮는다", {
    gr <- make_eta_grid(12L)
    length(gr) == 12L && abs(gr[1] - ETA_LO) < 1e-12 &&
        abs(gr[12] - ETA_HI) < 1e-10
})
check("eta 점 수가 3 미만이면 멈춘다",
      inherits(try(make_eta_grid(2L), silent = TRUE), "try-error"))
check("쉼표 목록의 공백을 턴다",
      identical(split_csv("Dementia, CN ,, MCI"), c("Dementia", "CN", "MCI")))

hash_file <- tempfile()
con <- file(hash_file, open = "wb")
writeChar("abc", con, eos = NULL)
close(con)
expected_hash <-
    "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
check(
    "SHA-256이 알려진 abc 지문과 같다",
    identical(sha256_file(hash_file), expected_hash)
)
fingerprints <- code_fingerprints(c(test_input = hash_file))
check("지문 이름을 보존한다", identical(names(fingerprints), "test_input"))
check(
    "지문 값이 단일 파일 해시와 같다",
    identical(unname(fingerprints[[1]]), sha256_file(hash_file))
)
unlink(hash_file)

refit_meta <- make_refit_run_meta(
    option = "C",
    groups = c("CN", "MCI"), arm = "projected", gamma = 0.5,
    verdict_paths = c("results/a.csv"),
    verdict_sha256 = list(`results/a.csv` = "deadbeef"),
    code_sha256 = list(driver = "c0ffee")
)
check("refit stage를 기록한다", identical(refit_meta$stage, "adni-refit-selected"))
check("refit 군을 보존한다", identical(refit_meta$groups, c("CN", "MCI")))
check("refit 입력 지문을 보존한다", identical(
    refit_meta$verdict_sha256[["results/a.csv"]], "deadbeef"
))
check("refit 코드 지문을 보존한다", identical(
    refit_meta$code_sha256$driver, "c0ffee"
))

cat(sprintf("\n통과 %d / 실패 %d\n", .n_pass, .n_fail))
if (.n_fail > 0L) quit(status = 1L)
