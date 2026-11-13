kgmselect_fingerprint <- function() {
    ns <- asNamespace("kgmselect")
    funs <- Filter(function(f) is.function(get(f, envir = ns)),
                   sort(ls(ns, all.names = TRUE)))
    digest::digest(list(version = as.character(utils::packageVersion("kgmselect")),
                        code = lapply(funs, function(f) deparse(get(f, envir = ns)))),
                   algo = "sha256")
}

driver_code_hash <- function(dir) {
    files <- sort(list.files(dir, pattern = "[.]R$", full.names = TRUE))
    if (!length(files)) stop("no driver library in ", dir)
    hashes <- vapply(files, function(f) digest::digest(file = f, algo = "sha256"), "")
    digest::digest(stats::setNames(unname(hashes), basename(files)), algo = "sha256")
}

platform_info <- function() {
    s <- utils::sessionInfo()
    c(platform = R.version$platform, blas = s$BLAS, lapack = s$LAPACK)
}

dependency_versions <- function() {
    pkgs <- c("kgmselect", "huge", "glasso", "expm", "rTensor", "igraph")
    vapply(pkgs, function(p) {
        if (requireNamespace(p, quietly = TRUE)) {
            as.character(utils::packageVersion(p))
        } else {
            NA_character_
        }
    }, "")
}

dataset_signature <- function(row, rep, config) {
    digest::digest(list(setting = as.list(row), rep = as.integer(rep),
                        seeds = dataset_seeds(row, rep), eta = config$eta,
                        method = config$method, tol = 0, n_sweep = 2L,
                        optimisation = list(boundary_margin = 1e-8, u_cap = 10,
                                            maxit = 100L, factr = 1e10),
                        kgmselect = config$fingerprint, driver = config$driver_code,
                        r = R.version.string, packages = config$versions,
                        platform = config$platform),
                   algo = "sha256")
}

reusable_result <- function(file, signature) {
    if (!file.exists(file)) return(NULL)
    old <- readRDS(file)
    if (!identical(old$signature, signature)) {
        stop("a result made with a different setup exists: ", file,
             "; use a new output folder")
    }
    old
}

run_information <- function(experiment, settings_file, config, workers) {
    s <- utils::sessionInfo()
    settings_sha256 <- digest::digest(file = settings_file, algo = "sha256")
    list(experiment = experiment,
         run_signature = digest::digest(list(settings_sha256, config$method,
                                             config$eta, config$fingerprint,
                                             config$driver_code, config$versions,
                                             R.version.string, config$platform),
                                        algo = "sha256"),
         started = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
         settings_file = basename(settings_file),
         settings_sha256 = settings_sha256,
         method = config$method,
         eta = list(length = length(config$eta), range = range(config$eta)),
         workers = workers, r = R.version.string, platform = R.version$platform,
         blas = s$BLAS, lapack = s$LAPACK, packages = as.list(config$versions),
         kgmselect_fingerprint = config$fingerprint, driver_code = config$driver_code)
}
