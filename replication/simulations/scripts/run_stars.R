here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                      value = TRUE)[1])))
here <- dirname(here)
for (f in list.files(file.path(here, "R"), pattern = "[.]R$", full.names = TRUE)) source(f)
opts <- parse_options(commandArgs(trailingOnly = TRUE), list(
    out = file.path(dirname(here), "outputs", "stars"), workers = "1", reps = "",
    settings_dir = file.path(here, "settings")))
load_kgmselect()
run_stars_experiment("stars", file.path(opts$settings_dir, "matrix.csv"),
                     file.path(opts$settings_dir, "stars_seeds.csv"), opts$out,
                     driver_code = driver_code_hash(file.path(here, "R")),
                     workers = as.integer(opts$workers),
                     reps = if (nzchar(opts$reps)) as.integer(opts$reps) else NULL)
