source("helper_lib.R")
source("helper_outputs.R")

script <- normalizePath(file.path("..", "..", "scripts", "summarize_results.R"))
aggregate <- function(...) {
    system2(file.path(R.home("bin"), "Rscript"), c(script, ...), stdout = TRUE, stderr = TRUE)
}
results <- tempfile()
log <- aggregate(paste0("--outputs=", root), paste0("--results=", results),
                 paste0("--settings_dir=", settings_dir))
expect_true(any(grepl("^TableB3 +skipped: no results of tensor$", log)))
expect_true(any(grepl("^figure2_path +220 rows$", log)))
tables <- paste0("TableB", c(1, 2, 4, 5, 6, 7), ".csv")
expect_identical(sort(list.files(file.path(results, "tables"))), sort(tables))
figures <- paste0(c("figure1_curves", "figure1_medians", "figure2_path", "figure3_mcc"), ".csv")
expect_identical(sort(list.files(file.path(results, "figures"))), figures)
m <- jsonlite::read_json(file.path(results, "manifest.json"))
expect_identical(vapply(m$per_dataset$files, `[[`, "", "file"),
                 paste0("per_dataset/", c("figure1_scores", "paths", "pilot_candidates",
                                          "selected"), ".csv.gz"))
expect_identical(m$per_dataset$sources[[1]]$folder, "matrix")
expect_equal(c(m$experiments$matrix$datasets, m$experiments$matrix$selections), c(4, 24))
expect_identical(m$skipped, list(TableB3 = "no results of tensor",
                                 Table1 = "no results of tensor"))
expect_identical(vapply(m$tables, `[[`, "", "label")[1:2],
                 c("tab:sim-fdr-tpr", "tab:app-m2-selection"))
expect_identical(vapply(m$figures, `[[`, "", "file"), paste0("figures/", figures))
expect_equal(vapply(m$figures, `[[`, 0, "rows")[3], 220)
expect_identical(m$figures_skipped, setNames(list(), character()))
for (t in c(m$tables, m$per_dataset$files, m$figures)) {
    expect_identical(t$sha256, digest::digest(file = file.path(results, t$file), algo = "sha256"))
}

files <- list.files(results, recursive = TRUE, full.names = TRUE)
before <- tools::md5sum(files)
aggregate(paste0("--outputs=", root), paste0("--results=", results),
          paste0("--settings_dir=", settings_dir))
expect_identical(tools::md5sum(files), before)
again <- tempfile()
dir.create(again)
file.copy(file.path(results, "per_dataset"), again, recursive = TRUE)
aggregate("--from=per_dataset", paste0("--results=", again))
written <- c(file.path("tables", tables), file.path("figures", figures), "manifest.json")
expect_identical(unname(tools::md5sum(file.path(again, written))),
                 unname(tools::md5sum(file.path(results, written))))

log <- suppressWarnings(aggregate("--from=elsewhere"))
expect_true(any(grepl("--from is outputs or per_dataset", log)))
expect_identical(attr(log, "status"), 1L)
log <- suppressWarnings(aggregate("--from=per_dataset", paste0("--results=", tempfile())))
expect_true(any(grepl("no per-dataset results in", log)))
expect_identical(attr(log, "status"), 1L)

writeLines("obsolete", file.path(results, "tables", "table5.csv"))
writeLines("obsolete", file.path(results, "tables", "TableB3.csv"))
writeLines("keep", file.path(results, "tables", "user_notes.csv"))
aggregate(paste0("--outputs=", root), paste0("--results=", results),
          paste0("--settings_dir=", settings_dir))
expect_false(file.exists(file.path(results, "tables", "table5.csv")))
expect_false(file.exists(file.path(results, "tables", "TableB3.csv")))
expect_identical(readLines(file.path(results, "tables", "user_notes.csv")), "keep")
