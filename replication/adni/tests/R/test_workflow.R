run_workflow_test <- function() {
    root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."))
    directory <- tempfile("adni-workflow-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE))
    input <- file.path(directory, "input")
    dir.create(file.path(input, "combat"), recursive = TRUE)
    dir.create(file.path(input, "combat_pb102"))
    set.seed(37)
    Y <- array(stats::rnorm(40L * 4L * 102L), c(40L, 4L, 102L))
    X <- data.frame(Age = stats::rnorm(40L), Sex = rep(0:1, 20L), APOE4 = rep(0:2, length.out = 40L))
    np <- reticulate::import("numpy")
    np$save(file.path(input, "combat_pb102/Y_CN_combat.npy"), Y)
    utils::write.csv(X, file.path(input, "combat/covars_CN.csv"), row.names = FALSE)
    jsonlite::write_json(list(groups = list(list(group = "CN", output = "Y_CN_combat.npy", problems = list()))),
                        file.path(input, "combat_pb102/resample_manifest.json"), auto_unbox = TRUE)
    caller <- getwd()
    paths <- file.path(directory, "paths")
    run_analysis(list(aggregated = input, option = "C", groups = "CN", arms = "projected",
                      `n-eta` = "3", out = paths), root)
    verdict <- file.path(paths, "probe_verdict.csv")
    stopifnot(file.exists(verdict), identical(getwd(), caller))
    fits <- file.path(directory, "fits")
    save_selected_fits(list(aggregated = input, option = "C", groups = "CN", arm = "projected",
                           gamma = "0.5", verdict = verdict, out = fits), root)
    stopifnot(nrow(read.csv(file.path(fits, "selected_fits.csv"))) == 3L, identical(getwd(), caller))
    labels <- file.path(directory, "labels.txt")
    write.table(data.frame(index = 1:4, labels = c("A_L", "A_R", "B_L", "B_R")),
                labels, row.names = FALSE, col.names = FALSE, quote = FALSE)
    bootstrap <- file.path(directory, "bootstrap")
    run_bootstrap(list(aggregated = input, option = "C", groups = "CN", arm = "projected",
                       B = "0", seed = "2026", workers = "1", identity = "1",
                       out = bootstrap, pairs = "homotopic", refit = "0",
                       `atlas-labels` = labels), root)
    draws <- read.csv(file.path(bootstrap, "bootstrap_draws.csv"))
    metadata <- jsonlite::read_json(file.path(bootstrap, "RUN_bootstrap.json"), simplifyVector = TRUE)
    stopifnot(nrow(draws) == 2L, all(draws$ok), metadata$blas_threads == 1L,
              metadata$n_failed == 0L, identical(getwd(), caller))
    covariance <- file.path(directory, "covariances.csv")
    check_covariance(list(aggregated = input, option = "C", groups = "CN", arms = "projected",
                          `n-eta` = "3", `out-csv` = covariance), root)
    values <- read.csv(covariance)
    stopifnot(nrow(values) == 12L, all(is.finite(values$eigen_min)), identical(getwd(), caller))
    cat("ADNI workflow: 12 checks passed\n")
}
run_workflow_test()
