here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                      value = TRUE)[1])))
here <- dirname(here)
for (f in list.files(file.path(here, "R"), pattern = "[.]R$", full.names = TRUE)) source(f)
opts <- parse_options(commandArgs(trailingOnly = TRUE), list(
    out = file.path(dirname(here), "outputs", "units"), workers = "1", reps = "",
    settings = file.path(here, "settings", "matrix.csv")))
load_kgmselect()
run_units_experiment("units", opts$settings, opts$out,
                     constants = c(0.01, 0.1, 1, 10, 100, 1000),
                     criteria = c("ebic", "scale", "relaxation"), gammas = 0.5,
                     driver_code = driver_code_hash(file.path(here, "R")),
                     workers = as.integer(opts$workers),
                     reps = if (nzchar(opts$reps)) as.integer(opts$reps) else NULL)
