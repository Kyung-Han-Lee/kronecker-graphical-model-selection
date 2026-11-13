stars_seed <- function(source, base = 2026L) {
    hex <- vapply(source, function(s) substr(digest::digest(s, algo = "xxhash32"), 1L, 7L),
                  "", USE.NAMES = FALSE)
    (base + strtoi(hex, 16L)) %% .Machine$integer.max
}

read_stars_seeds <- function(file) {
    s <- utils::read.csv(file, stringsAsFactors = FALSE)
    if (!identical(names(s), c("dataset", "source", "seed"))) {
        stop("the seed table needs the columns dataset, source and seed")
    }
    if (anyDuplicated(s$dataset) || anyDuplicated(s$source)) {
        stop("the seed table repeats a dataset or a source name")
    }
    if (!is.integer(s$seed) || anyNA(s$seed)) {
        stop("the seed table holds a seed that is not an integer")
    }
    s
}

stars_config <- function(seeds, eta, driver_code, N = 20L, thresholds = c(0.05, 0.1)) {
    driver_config(list(name = "stars", N = N, thresholds = thresholds,
                       seed_table = digest::digest(seeds, algo = "sha256")),
                  eta, driver_code)
}

run_stars_dataset <- function(dataset, row, config, out, seed, path_fun) {
    file <- file.path(out, "datasets", paste0(dataset$id, ".rds"))
    signature <- dataset_signature(row, dataset$rep, config)
    if (!is.null(reusable_result(file, signature))) return("reused")
    gen <- make_dataset(row, dataset$rep)
    t0 <- proc.time()[["elapsed"]]
    path <- path_fun(dataset, gen$Y, config$eta)
    t1 <- proc.time()[["elapsed"]]
    st <- kgm_stars(path, gen$Y, N = config$method$N,
                    threshold = config$method$thresholds, seed = seed, tol = 0)
    t2 <- proc.time()[["elapsed"]]
    tables <- lapply(seq_along(st$threshold), function(i) {
        k <- st$index[i]
        fit_tables(data.frame(criterion = "stars", threshold = st$threshold[i]), k,
                   path$eta[k], path$fits[[k]], gen$truth)
    })
    save_atomic(list(id = dataset$id, setting = dataset$setting, rep = dataset$rep,
                     seeds = gen$seeds, stars_seed = seed, signature = signature,
                     data_sha256 = digest::digest(gen$Y, algo = "sha256"),
                     truth_sha256 = digest::digest(gen$truth, algo = "sha256"),
                     instability = st$instability, subsample_size = st$subsample_size,
                     selected = do.call(rbind, lapply(tables, `[[`, "selected")),
                     selected_modes = do.call(rbind, lapply(tables, `[[`, "selected_modes")),
                     seconds = c(path = t1 - t0, stars = t2 - t1)),
                file)
    "done"
}

run_stars_experiment <- function(experiment, settings_file, seeds_file, out,
                                 driver_code, workers = 1L, reps = NULL,
                                 path_fun = fit_path, eta = kgm_eta_grid(),
                                 N = 20L, thresholds = c(0.05, 0.1)) {
    seeds <- read_stars_seeds(seeds_file)
    if (!all(seeds$dataset %in% expand_datasets(read_settings(settings_file))$id)) {
        stop("the seed table names datasets that are not in ", basename(settings_file))
    }
    config <- stars_config(seeds, eta, driver_code, N, thresholds)
    seed_of <- stats::setNames(seeds$seed, seeds$dataset)
    run_datasets(experiment, settings_file, out, config,
                 function(dataset, row) {
                     run_stars_dataset(dataset, row, config, out, seed_of[[dataset$id]],
                                       path_fun)
                 },
                 workers = workers, reps = reps, ids = seeds$dataset)
}
