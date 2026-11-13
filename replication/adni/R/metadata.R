sha256_file <- function(path) {
    digest::digest(file = path, algo = "sha256")
}

code_fingerprints <- function(paths) {
    stopifnot(
        length(paths) >= 1L,
        !is.null(names(paths)),
        all(nzchar(names(paths)))
    )
    as.list(vapply(paths, sha256_file, character(1)))
}

make_refit_run_meta <- function(option, groups, arm, gamma,
                                verdict_paths, verdict_sha256, code_sha256) {
    list(
        stage = "adni-refit-selected",
        option = option,
        groups = groups,
        arm = arm,
        gamma = gamma,
        verdict_paths = verdict_paths,
        verdict_sha256 = verdict_sha256,
        code_sha256 = code_sha256
    )
}
