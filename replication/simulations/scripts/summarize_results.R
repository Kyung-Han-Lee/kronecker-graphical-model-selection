here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE),
                                                      value = TRUE)[1])))
here <- dirname(here)
for (f in list.files(file.path(here, "R"), pattern = "[.]R$", full.names = TRUE)) source(f)
for (f in list.files(file.path(here, "summaries"), pattern = "[.]R$", full.names = TRUE)) {
    source(f)
}
opts <- parse_options(commandArgs(trailingOnly = TRUE), list(
    outputs = file.path(dirname(here), "outputs"), results = "",
    settings_dir = file.path(here, "settings"), from = "outputs"))
if (!nzchar(opts$results)) opts$results <- file.path(opts$outputs, "results")
if (!opts$from %in% c("outputs", "per_dataset")) stop("--from is outputs or per_dataset")
per_dataset_dir <- file.path(opts$results, "per_dataset")
if (opts$from == "outputs") {
    write_per_dataset(collect_outputs(opts$outputs, opts$settings_dir), per_dataset_dir)
}
per_dataset <- read_per_dataset(per_dataset_dir)
if (is.null(per_dataset$selected)) stop("no per-dataset results in ", per_dataset_dir)
made <- make_tables(per_dataset)
tables_dir <- file.path(opts$results, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(tables_dir, pattern = "^(table[0-9]+|Table(1|B[1-7]))[.]csv$", full.names = TRUE))
for (name in names(made$tables)) {
    write_exact_csv(made$tables[[name]], file.path(tables_dir, paste0(name, ".csv")))
}
figures <- make_figures(per_dataset)
figures_dir <- file.path(opts$results, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(figures_dir, paste0(figure_value_files, ".csv")))
for (name in names(figures$files)) {
    write_exact_csv(figures$files[[name]], file.path(figures_dir, paste0(name, ".csv")))
}
write_manifest(opts$results, per_dataset, made, figures, here)
for (name in names(made$tables)) {
    cat(sprintf("%-8s %4d rows\n", name, nrow(made$tables[[name]])))
}
for (name in names(made$skipped)) cat(sprintf("%-8s skipped: %s\n", name, made$skipped[[name]]))
for (name in names(figures$files)) {
    cat(sprintf("%-15s %5d rows\n", name, nrow(figures$files[[name]])))
}
for (name in names(figures$skipped)) {
    cat(sprintf("%-15s skipped: %s\n", name, figures$skipped[[name]]))
}
