load_kgmselect <- function() {
    src <- Sys.getenv("KGMSELECT_SRC")
    if (nzchar(src)) {
        suppressMessages(pkgload::load_all(src, quiet = TRUE))
    } else {
        suppressPackageStartupMessages(library(kgmselect))
    }
}

fit_path <- function(dataset, Y, eta) kgm_path(Y, eta)

run_dataset <- function(dataset, row, config, out, path_fun) {
    file <- file.path(out, "datasets", paste0(dataset$id, ".rds"))
    signature <- dataset_signature(row, dataset$rep, config)
    if (!is.null(reusable_result(file, signature))) return("reused")
    gen <- make_dataset(row, dataset$rep)
    t0 <- proc.time()[["elapsed"]]
    path <- if (isTRUE(config$method$save_inputs)) {
        fit_path_with_inputs(dataset, gen$Y, config$eta)
    } else {
        path_fun(dataset, gen$Y, config$eta)
    }
    t1 <- proc.time()[["elapsed"]]
    scores <- score_path(path, gen$Y, config$method$criteria, config$method$gammas,
                         row$gradient)
    t2 <- proc.time()[["elapsed"]]
    tables <- selection_tables(scores, path, gen$truth)
    result <- list(id = dataset$id, setting = dataset$setting, rep = dataset$rep,
                   seeds = gen$seeds, signature = signature,
                   data_sha256 = digest::digest(gen$Y, algo = "sha256"),
                   truth_sha256 = digest::digest(gen$truth, algo = "sha256"),
                   scores = scores, selected = tables$selected,
                   selected_modes = tables$selected_modes,
                   recovery = candidate_recovery(path, gen$truth),
                   seconds = c(path = t1 - t0, scores = t2 - t1))
    if (isTRUE(config$method$support)) result$candidates <- support_records(path)
    if (isTRUE(config$method$save_inputs)) {
        result$inputs <- list(path = path, settings = row)
    }
    save_atomic(result, file)
    "done"
}

run_datasets <- function(experiment, settings_file, out, config, dataset_fun,
                         workers = 1L, reps = NULL, ids = NULL) {
    use_one_thread()
    settings <- read_settings(settings_file)
    datasets <- expand_datasets(settings, reps)
    if (!is.null(ids)) datasets <- datasets[datasets$id %in% ids, , drop = FALSE]
    info <- run_information(experiment, settings_file, config, workers)

    previous <- file.path(out, "run.json")
    if (file.exists(previous) &&
            !identical(jsonlite::read_json(previous)$run_signature, info$run_signature)) {
        stop("results with a different setup exist in ", out, "; use a new output folder")
    }
    dir.create(file.path(out, "datasets"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(out, "errors"), showWarnings = FALSE)
    jsonlite::write_json(info, file.path(out, "run.json"), auto_unbox = TRUE,
                         pretty = TRUE, digits = NA)
    status <- parallel::mclapply(seq_len(nrow(datasets)), function(i) {
        dataset <- datasets[i, ]
        row <- settings[settings$setting == dataset$setting, ]
        tryCatch(dataset_fun(dataset, row), error = function(e) {
            writeLines(conditionMessage(e), file.path(out, "errors", paste0(dataset$id, ".txt")))
            "failed"
        })
    }, mc.cores = workers, mc.preschedule = FALSE)
    status <- vapply(status, function(s) if (is.character(s)) s[1] else "failed", "")
    counts <- table(factor(status, levels = c("done", "reused", "failed")))
    info$finished <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    info$datasets <- as.list(c(total = nrow(datasets), counts))
    jsonlite::write_json(info, file.path(out, "run.json"), auto_unbox = TRUE,
                         pretty = TRUE, digits = NA)
    cat(sprintf("%s: %d datasets, %d done, %d reused, %d failed\n", experiment,
                nrow(datasets), counts[["done"]], counts[["reused"]], counts[["failed"]]))
    if (counts[["failed"]] > 0L) {
        stop(counts[["failed"]], " datasets failed; see ", file.path(out, "errors"))
    }
    invisible(status)
}

run_experiment <- function(experiment, settings_file, out, criteria, gammas,
                           driver_code, workers = 1L, reps = NULL, ids = NULL,
                           path_fun = fit_path, eta = kgm_eta_grid(),
                           support = FALSE, save_inputs = FALSE) {
    config <- path_config(criteria, gammas, eta, driver_code, support, save_inputs)
    run_datasets(experiment, settings_file, out, config,
                 function(dataset, row) run_dataset(dataset, row, config, out, path_fun),
                 workers = workers, reps = reps, ids = ids)
}

use_one_thread <- function() {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
    if (RhpcBLASctl::blas_get_num_procs() != 1L) {
        stop("BLAS uses more than one thread; set OPENBLAS_NUM_THREADS=1")
    }
    invisible(TRUE)
}
