run_analysis_defaults_test <- function() {
    root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."))
    directory <- tempfile("adni-analysis-defaults-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE))
    staged <- file.path(directory, "replication", "adni")
    dir.create(file.path(staged, "scripts"), recursive = TRUE)
    stopifnot(file.copy(file.path(root, "R"), staged, recursive = TRUE))
    script <- file.path(staged, "scripts", "run_analysis.R")
    stopifnot(file.copy(file.path(root, "scripts", "run_analysis.R"), script))

    input <- file.path(directory, "input")
    dir.create(file.path(input, "combat"), recursive = TRUE)
    dir.create(file.path(input, "combat_pb102"))
    set.seed(37)
    Y <- array(stats::rnorm(40L * 4L * 102L), c(40L, 4L, 102L))
    X <- data.frame(Age = stats::rnorm(40L), Sex = rep(0:1, 20L),
                    APOE4 = rep(0:2, length.out = 40L))
    groups <- c("CN", "MCI", "Dementia")
    np <- reticulate::import("numpy")
    for (group in groups) {
        np$save(file.path(input, "combat_pb102", paste0("Y_", group, "_combat.npy")), Y)
        utils::write.csv(X, file.path(input, "combat", paste0("covars_", group, ".csv")),
                         row.names = FALSE)
    }
    jsonlite::write_json(list(groups = lapply(groups, function(group) {
        list(group = group, output = paste0("Y_", group, "_combat.npy"), problems = list())
    })), file.path(input, "combat_pb102", "resample_manifest.json"), auto_unbox = TRUE)

    invoke <- function(extra, label) {
        log <- file.path(directory, paste0(label, ".txt"))
        arguments <- c("--vanilla", script, "--aggregated", input, "--n-eta", "3", extra)
        status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(arguments),
                           stdout = log, stderr = log)
        if (status != 0L) {
            stop(label, " ADNI invocation failed:\n", paste(tail(readLines(log), 10L), collapse = "\n"))
        }
    }
    invoke(character(), "default")
    output <- file.path(directory, "replication", "outputs", "adni", "joint_profile")
    stopifnot(file.exists(file.path(output, "probe_path.csv")))
    path <- read.csv(file.path(output, "probe_path.csv"))
    stopifnot(identical(unique(path$group), groups),
              all(path$arm == "projected"), all(path$p_time == 102L))

    override <- file.path(directory, "override")
    invoke(c("--groups", "MCI", "--arms", "projected", "--out", override), "override")
    custom <- read.csv(file.path(override, "probe_path.csv"))
    stopifnot(identical(unique(custom$group), "MCI"))
    columns <- c("eta_index", "gamma", "ebic_plug", "ebic_prof", "ebic_joint_uncounted")
    stopifnot(isTRUE(all.equal(path[path$group == "MCI", columns], custom[, columns],
                              tolerance = 0, check.attributes = FALSE)))
    cat("ADNI analysis defaults: 6 checks passed\n")
}
run_analysis_defaults_test()
