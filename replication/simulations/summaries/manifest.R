sha256_of <- function(file) digest::digest(file = file, algo = "sha256")

write_manifest <- function(results, per_dataset, made, figures, here) {
    sel <- per_dataset$selected
    files <- sort(list.files(file.path(results, "per_dataset"), pattern = "[.]csv[.]gz$"))
    counts <- lapply(split(sel, sel$experiment), function(x) {
        rules <- unique(x[c("criterion", "gamma", "constant", "threshold")])
        list(datasets = length(unique(x$id)), selections = nrow(x),
             rules = lapply(seq_len(nrow(rules)), function(i) {
                 r <- as.list(rules[i, ])
                 r[vapply(r, is.na, TRUE)] <- NULL
                 r
             }))
    })
    tables <- lapply(names(made$tables), function(name) {
        note <- table_notes[table_notes$table == name, ]
        file <- file.path(results, "tables", paste0(name, ".csv"))
        list(file = file.path("tables", paste0(name, ".csv")), label = note$label,
             caption = note$caption, printed = note$printed,
             rows = nrow(made$tables[[name]]), sha256 = sha256_of(file))
    })
    manifest <- list(
        aggregate_code = driver_code_hash(file.path(here, "summaries")),
        per_dataset = list(
            files = lapply(files, function(f) {
                list(file = file.path("per_dataset", f),
                     rows = nrow(per_dataset[[sub("[.]csv[.]gz$", "", f)]]),
                     sha256 = sha256_of(file.path(results, "per_dataset", f)))
            }),
            sources = jsonlite::read_json(file.path(results, "per_dataset", "sources.json"))),
        experiments = counts,
        tables = tables,

        skipped = if (length(made$skipped)) made$skipped else stats::setNames(list(), character()),
        figures = lapply(names(figures$files), function(name) {
            file <- file.path("figures", paste0(name, ".csv"))
            list(file = file, rows = nrow(figures$files[[name]]),
                 sha256 = sha256_of(file.path(results, file)))
        }),
        figures_skipped = if (length(figures$skipped)) {
            figures$skipped
        } else {
            stats::setNames(list(), character())
        })
    jsonlite::write_json(manifest, file.path(results, "manifest.json"), auto_unbox = TRUE,
                         pretty = TRUE, digits = NA, null = "null", na = "null")
}
