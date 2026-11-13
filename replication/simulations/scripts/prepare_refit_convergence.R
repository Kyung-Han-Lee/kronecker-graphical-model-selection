here <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=",
    commandArgs(FALSE), value = TRUE)[1]))))
for (file in list.files(file.path(here, "R"), pattern = "[.]R$", full.names = TRUE)) source(file)
options <- parse_options(commandArgs(trailingOnly = TRUE), list(
    outputs = file.path(dirname(here), "outputs"),
    inputs = file.path(dirname(here), "outputs", "inputs"), jobs = "", out = ""))
if (!nzchar(options$jobs) || !nzchar(options$out)) stop("--jobs and --out are required")
load_kgmselect()
prepare_refit_convergence(options$outputs, options$inputs, options$jobs, options$out)
