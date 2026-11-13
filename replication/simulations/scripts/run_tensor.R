here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                      value = TRUE)[1])))
here <- dirname(here)
for (f in list.files(file.path(here, "R"), pattern = "[.]R$", full.names = TRUE)) source(f)
opts <- parse_options(commandArgs(trailingOnly = TRUE), list(
    out = file.path(dirname(here), "outputs", "tensor"), workers = "1", reps = "", save_inputs = "0",
    settings = file.path(here, "settings", "tensor.csv")))
load_kgmselect()
run_experiment("tensor", opts$settings, opts$out,
               criteria = c("ebic", "scale", "relaxation"), gammas = 0.5,
               driver_code = driver_code_hash(file.path(here, "R")),
               workers = as.integer(opts$workers),
               reps = if (nzchar(opts$reps)) as.integer(opts$reps) else NULL, save_inputs = identical(opts$save_inputs, "1"))
