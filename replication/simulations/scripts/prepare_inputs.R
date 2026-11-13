here <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=",
    commandArgs(FALSE), value = TRUE)[1]))))
for (file in list.files(file.path(here, "R"), pattern = "[.]R$", full.names = TRUE)) source(file)
options <- parse_options(commandArgs(trailingOnly = TRUE), list(
    outputs = file.path(dirname(here), "outputs"),
    out = file.path(dirname(here), "outputs", "inputs"), profiles = "paper"))
load_kgmselect()
prepare_candidate_inputs(options$outputs, options$out, options$profiles)
