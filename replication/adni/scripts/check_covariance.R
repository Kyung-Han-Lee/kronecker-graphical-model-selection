arguments <- commandArgs(trailingOnly = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
for (file in list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)) source(file)
options <- parse_adni_arguments(arguments, list(option = "C", groups = "Dementia,MCI,CN", arms = "control,projected", `n-eta` = as.character(N_ETA_DEFAULT),
    `out-csv` = "../outputs/adni/diagnostics/projected_probe/assumption_check.csv", aggregated = ""))
options <- resolve_adni_arguments(options, root, getwd())
check_covariance(options, root)
