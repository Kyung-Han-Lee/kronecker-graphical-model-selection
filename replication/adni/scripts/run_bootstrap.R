arguments <- commandArgs(trailingOnly = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
for (file in list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)) source(file)
options <- parse_adni_arguments(arguments, list(option = "C", groups = "Dementia,MCI,CN", arm = "projected", B = "500", seed = "2026",
    workers = "1", identity = "1", out = "../outputs/adni/bootstrap", pairs = "reported",
    refit = "0", aggregated = "", `atlas-labels` = file.path(root, "deps/atlas/ROI_MNI_V4.txt")))
options <- resolve_adni_arguments(options, root, getwd())
run_bootstrap(options, root)
