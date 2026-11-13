root <- dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))))
for (file in list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE)) source(file)
load_adni_dependencies()
for (file in list.files(file.path(root, "tests/R"), pattern = "^test_.*[.]R$", full.names = TRUE)) {
    cat(basename(file), "\n")
    source(file)
}
