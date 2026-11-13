arguments <- commandArgs(trailingOnly = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
for (file in list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)) source(file)
options <- parse_adni_arguments(arguments, list(option = "C", groups = "Dementia,MCI,CN", arm = "projected", gamma = "0.5", verdict = "../outputs/adni/joint_profile/probe_verdict.csv",
    out = "../outputs/adni/selected_fits", aggregated = ""))
options <- resolve_adni_arguments(options, root, getwd())
save_selected_fits(options, root)
