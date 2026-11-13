here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                      value = TRUE)[1])))
here <- dirname(here)
for (f in list.files(file.path(here, "R"), pattern = "[.]R$", full.names = TRUE)) source(f)
opts <- parse_options(commandArgs(trailingOnly = TRUE), list(
    out = file.path(dirname(here), "outputs", "cv"), workers = "1", reps = "",
    experiments = "matrix,tensor", settings_dir = file.path(here, "settings")))
load_kgmselect()
for (experiment in strsplit(opts$experiments, ",", fixed = TRUE)[[1]]) {
    run_cv_experiment(experiment,
                      file.path(opts$settings_dir, paste0(experiment, ".csv")),
                      file.path(opts$out, experiment),
                      driver_code = driver_code_hash(file.path(here, "R")),
                      workers = as.integer(opts$workers),
                      reps = if (nzchar(opts$reps)) as.integer(opts$reps) else NULL)
}
