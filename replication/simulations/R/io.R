save_atomic <- function(object, file) {
    tmp <- tempfile(pattern = ".partial-", tmpdir = dirname(file))
    saveRDS(object, tmp)
    if (!file.rename(tmp, file)) {
        unlink(tmp)
        stop("could not move the result into place: ", file)
    }
    invisible(file)
}

per_dataset_files <- c(selected = "selected", paths = "paths",
                       pilot_candidates = "pilot_candidates", figure1_scores = "figure1_scores")

append_exact_csv <- function(x, path) {
    if (!nrow(x)) return(invisible(path))
    text <- vapply(x, is.character, TRUE)
    for (name in names(x)[vapply(x, is.double, TRUE)]) x[[name]] <- exact_text(x[[name]])
    present <- file.exists(path)
    utils::write.table(x, path, sep = ",", row.names = FALSE, col.names = !present,
                       append = present, na = "", quote = which(text))
    invisible(path)
}

exact_text <- function(x) {
    out <- ifelse(is.na(x), NA_character_, sprintf("%.15g", x))
    for (digits in 16:17) {
        bad <- !is.na(x) & as.numeric(out) != x
        out[bad] <- sprintf("%.*g", digits, x[bad])
    }
    out
}

write_exact_csv <- function(x, file) {
    text <- vapply(x, is.character, TRUE)
    for (col in names(x)[vapply(x, is.double, TRUE)]) x[[col]] <- exact_text(x[[col]])
    con <- if (grepl("[.]gz$", file)) gzfile(file, "w") else file(file, "w")
    on.exit(close(con))
    utils::write.csv(x, con, row.names = FALSE, na = "", quote = which(text))
}

write_per_dataset <- function(collected, dir) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    unlink(c(list.files(dir, pattern = "[.]csv[.]gz$", full.names = TRUE),
             file.path(dir, "sources.json")))
    for (name in names(collected$tables)) {
        write_exact_csv(collected$tables[[name]], file.path(dir, paste0(name, ".csv.gz")))
    }
    jsonlite::write_json(collected$sources, file.path(dir, "sources.json"),
                         auto_unbox = TRUE, pretty = TRUE, digits = NA)
    invisible(dir)
}

read_per_dataset <- function(dir) {
    files <- list.files(dir, pattern = "[.]csv[.]gz$", full.names = TRUE)
    tables <- lapply(files, function(f) {
        utils::read.csv(f, stringsAsFactors = FALSE, na.strings = "")
    })
    names(tables) <- sub("[.]csv[.]gz$", "", basename(files))
    tables
}
