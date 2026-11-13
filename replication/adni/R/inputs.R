split_csv <- function(s) {
    parts <- trimws(strsplit(s, ",")[[1]])
    parts[nzchar(parts)]
}

check_resample_manifest <- function(aggregated_dir, subdir, group,
                                    required = TRUE) {

    manifest_path <- file.path(
        aggregated_dir, subdir, "resample_manifest.json"
    )
    if (!file.exists(manifest_path)) {
        if (!required) {
            return(invisible(NULL))
        }
        stop(sprintf(
            "%s 가 없다. %s/%s 가 검증을 거친 산출물인지 확인할 수 없다.",
            manifest_path, subdir, group
        ))
    }
    manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
    matched <- Filter(
        function(record) identical(record$group, group), manifest$groups
    )
    if (!required && !length(matched)) {
        return(invisible(NULL))
    }
    if (length(matched) != 1L) {
        stop(sprintf(
            "%s 에 %s 의 기록이 %d 개다 (1 개여야 한다). 부분 재실행으로 %s",
            manifest_path, group, length(matched),
            "낡은 텐서가 남아 있을 수 있다."
        ))
    }
    record <- matched[[1]]
    if (length(record$problems)) {
        stop(sprintf(
            "%s/%s 는 재표집 검증에 실패한 산출물이다: %s",
            subdir, group, paste(unlist(record$problems), collapse = "; ")
        ))
    }
    if (!length(record$output) || !nzchar(unlist(record$output)[1])) {
        stop(sprintf(
            "%s/%s 의 manifest 기록에 출력 경로가 없다. 저장되지 않은 실행이다.",
            subdir, group
        ))
    }
    invisible(NULL)
}

load_arm <- function(group, arm, aggregated_dir, np) {
    spec <- ARM_SPECS[[arm]]
    y_path <- file.path(
        aggregated_dir, spec$subdir, sprintf("Y_%s_combat.npy", group)
    )
    c_path <- file.path(
        aggregated_dir, "combat", sprintf("covars_%s.csv", group)
    )
    if (!file.exists(y_path)) stop("텐서가 없다: ", y_path)
    if (!file.exists(c_path)) stop("공변량이 없다: ", c_path)

    Y_loaded <- np$load(y_path)
    Y <- array(as.numeric(Y_loaded), dim = dim(Y_loaded))
    covars <- utils::read.csv(c_path)
    X <- as.matrix(covars[, c("Age", "Sex", "APOE4")])
    stopifnot(length(dim(Y)) == 3L, nrow(X) == dim(Y)[1], !anyNA(X))

    if (dim(Y)[3] != spec$p_time) {
        stop(sprintf(
            "%s 팔의 시간 차원이 %d 인데 %d 를 기대했다 (%s). 재표집 산출물이 다른 --tr-out 으로 만들어졌을 수 있다.",
            arm, dim(Y)[3], spec$p_time, y_path
        ))
    }
    check_resample_manifest(
        aggregated_dir, spec$subdir, group,
        required = isTRUE(spec$requires_manifest)
    )
    list(Y = Y, X = X)
}

residualize <- function(Y, X) {
    n <- dim(Y)[1]
    Y_mat <- matrix(Y, nrow = n)
    design <- cbind(intercept = 1, X)
    qr_design <- qr(design)
    if (qr_design$rank < ncol(design)) {
        stop(sprintf(
            "설계행렬 랭크 결손: rank %d < %d 열. 공변량이 상수이거나 겹친다.",
            qr_design$rank, ncol(design)
        ))
    }
    array(Y_mat - qr.fitted(qr_design, Y_mat), dim = dim(Y))
}

read_aal_labels <- function(path = "deps/atlas/ROI_MNI_V4.txt") {
    as.character(utils::read.table(path, stringsAsFactors = FALSE)[[2]])
}

homotopic_pairs <- function(labels) {
    left <- which(endsWith(labels, "_L"))
    out <- list()
    for (i in left) {
        stem <- substr(labels[i], 1L, nchar(labels[i]) - 2L)
        j <- which(labels == paste0(stem, "_R"))
        if (length(j) == 1L) out[[stem]] <- c(i, j)
    }
    out
}
