parse_adni_arguments <- function(argv, defaults) {
    if ("--help" %in% argv) {
        cat("Options use --name value or --name=value. Defaults:\n")
        for (key in names(defaults)) cat("  --", key, " ", defaults[[key]], "\n", sep = "")
        quit(status = 0L)
    }
    out <- defaults
    i <- 1L
    while (i <= length(argv)) {
        if (!startsWith(argv[i], "--")) stop("expected an option: ", argv[i])
        text <- substring(argv[i], 3L)
        equal <- regexpr("=", text, fixed = TRUE)[1L]
        if (equal > 0L) {
            key <- substr(text, 1L, equal - 1L)
            value <- substring(text, equal + 1L)
        } else {
            key <- text
            if (i == length(argv)) stop("missing value for --", key)
            i <- i + 1L
            value <- argv[i]
        }
        if (!key %in% names(defaults)) stop("unknown option: --", key)
        out[[key]] <- value
        i <- i + 1L
    }
    out
}

resolve_adni_arguments <- function(options, root, caller) {
    resolve <- function(value) {
        parts <- strsplit(value, ",", fixed = TRUE)[[1]]
        paths <- vapply(parts, function(path) {
            if (startsWith(path, "/") || grepl("^[A-Za-z]:", path)) return(path)
            base <- if (startsWith(path, "../outputs/adni/")) root else caller
            file.path(base, path)
        }, "")
        paste(paths, collapse = ",")
    }
    for (key in intersect(c("out", "out-csv", "verdict", "aggregated", "atlas-labels"), names(options))) {
        if (nzchar(options[[key]])) options[[key]] <- resolve(options[[key]])
    }
    options
}

load_adni_dependencies <- function() {
    suppressPackageStartupMessages({
        library(reticulate)
        library(RhpcBLASctl)
        library(jsonlite)
    })
    Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
    if (RhpcBLASctl::blas_get_num_procs() != 1L) stop("BLAS must use one thread")
    RhpcBLASctl::blas_get_num_procs()
}
