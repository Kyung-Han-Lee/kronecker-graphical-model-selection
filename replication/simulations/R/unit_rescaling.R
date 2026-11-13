units_config <- function(constants, criteria, gammas, eta, driver_code) {
    driver_config(list(name = "units", constants = constants, criteria = criteria,
                       gammas = gammas), eta, driver_code)
}

run_units_dataset <- function(dataset, row, config, out, path_fun) {
    file <- file.path(out, "datasets", paste0(dataset$id, ".rds"))
    signature <- dataset_signature(row, dataset$rep, config)
    if (!is.null(reusable_result(file, signature))) return("reused")
    gen <- make_dataset(row, dataset$rep)
    parts <- lapply(config$method$constants, function(cc) {
        Y <- cc * gen$Y
        t0 <- proc.time()[["elapsed"]]
        path <- path_fun(dataset, Y, cc^2 * config$eta)
        t1 <- proc.time()[["elapsed"]]
        scores <- score_path(path, Y, config$method$criteria, config$method$gammas,
                             row$gradient)
        t2 <- proc.time()[["elapsed"]]
        tables <- selection_tables(scores, path, gen$truth)
        list(scores = cbind(constant = cc, scores),
             selected = cbind(constant = cc, tables$selected),
             selected_modes = cbind(constant = cc, tables$selected_modes),
             seconds = c(constant = cc, path = t1 - t0, scores = t2 - t1))
    })
    stack <- function(name) do.call(rbind, lapply(parts, `[[`, name))
    save_atomic(list(id = dataset$id, setting = dataset$setting, rep = dataset$rep,
                     seeds = gen$seeds, signature = signature,
                     data_sha256 = digest::digest(gen$Y, algo = "sha256"),
                     truth_sha256 = digest::digest(gen$truth, algo = "sha256"),
                     scores = stack("scores"), selected = stack("selected"),
                     selected_modes = stack("selected_modes"), seconds = stack("seconds")),
                file)
    "done"
}

run_units_experiment <- function(experiment, settings_file, out, constants, criteria,
                                 gammas, driver_code, workers = 1L, reps = NULL,
                                 ids = NULL, path_fun = fit_path, eta = kgm_eta_grid()) {
    config <- units_config(constants, criteria, gammas, eta, driver_code)
    run_datasets(experiment, settings_file, out, config,
                 function(dataset, row) run_units_dataset(dataset, row, config, out, path_fun),
                 workers = workers, reps = reps, ids = ids)
}
