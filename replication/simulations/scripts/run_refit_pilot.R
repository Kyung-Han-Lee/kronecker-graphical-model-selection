here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                      value = TRUE)[1])))
here <- dirname(here)
for (f in list.files(file.path(here, "R"), pattern = "[.]R$", full.names = TRUE)) source(f)
opts <- parse_options(commandArgs(trailingOnly = TRUE), list(
    out = file.path(dirname(here), "outputs", "refit_pilot"), workers = "1", reps = "",
    settings = file.path(here, "settings", "refit_pilot.csv")))
load_kgmselect()
run_pilot_experiment("refit_pilot", opts$settings, opts$out,
                     driver_code = driver_code_hash(file.path(here, "R")),
                     workers = as.integer(opts$workers),
                     reps = if (nzchar(opts$reps)) as.integer(opts$reps) else NULL,
                     seconds = 30L, gamma = 0.5)
