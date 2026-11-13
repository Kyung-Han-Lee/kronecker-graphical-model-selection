collected_experiments <- data.frame(
    experiment = c("matrix", "tensor", "cv", "cv", "stars", "units", "refit_pilot"),
    folder = c("matrix", "tensor", file.path("cv", "matrix"), file.path("cv", "tensor"),
               "stars", "units", "refit_pilot"),
    stringsAsFactors = FALSE)

design_columns <- c("M", "family", "dims", "degree", "n")

choice_columns <- c("criterion", "gamma", "constant", "threshold")

count_columns <- c("tp", "fp", "fn", "tn")

setting_design <- function(settings) {
    data.frame(setting = settings$setting, M = settings$M, family = settings$family,
               dims = settings$dims,
               degree = ifelse(settings$family == "er", settings$d, settings$m_new),
               n = settings$n, stringsAsFactors = FALSE)
}

adjacent_pair <- function(exact) any(exact[-1] & exact[-length(exact)])

widen_modes <- function(selected, modes, max_modes) {
    labels <- intersect(choice_columns, names(selected))
    key <- function(x) do.call(paste, c(unname(as.list(x[labels])), sep = "\r"))
    out <- matrix(NA_real_, nrow(selected), 4L * max_modes, dimnames = list(
        NULL, paste0(rep(count_columns, max_modes), "_", rep(seq_len(max_modes), each = 4L))))
    selected_key <- key(selected)
    modes_key <- key(modes)
    for (i in seq_len(nrow(selected))) {
        m <- modes[modes_key == selected_key[i], , drop = FALSE]
        for (j in seq_len(nrow(m))) {
            out[i, paste0(count_columns, "_", m$mode[j])] <- unlist(m[j, count_columns])
        }
    }
    as.data.frame(out)
}

collect_folder <- function(experiment, outputs, folder, settings_dir, max_modes, figure1) {
    run <- jsonlite::read_json(file.path(outputs, folder, "run.json"), simplifyVector = TRUE)
    settings_file <- file.path(settings_dir, run$settings_file)
    if (!identical(digest::digest(file = settings_file, algo = "sha256"), run$settings_sha256)) {
        stop("the results in ", folder, " were made with another ", run$settings_file)
    }
    design <- setting_design(read_settings(settings_file))
    files <- sort(list.files(file.path(outputs, folder, "datasets"), pattern = "[.]rds$",
                             full.names = TRUE))
    results <- lapply(files, readRDS)
    one <- function(r) {
        id <- data.frame(experiment = experiment, id = r$id, setting = r$setting, rep = r$rep,
                         design[match(r$setting, design$setting), design_columns],
                         stringsAsFactors = FALSE, row.names = NULL)
        sel <- r$selected

        if (!"gamma" %in% names(sel) && length(run$method$gamma) == 1L) {
            sel$gamma <- run$method$gamma
        }
        for (col in setdiff(choice_columns, names(sel))) sel[[col]] <- NA
        sel <- cbind(id[rep(1L, nrow(sel)), ], sel[choice_columns],
                     grid_points = run$eta$length,
                     sel[c("k", "eta", "edges", "empty", count_columns, "mcc", "fdp", "tpr")],
                     widen_modes(r$selected, r$selected_modes, max_modes))
        paths <- if (!is.null(r$recovery)) {
            exact <- r$recovery$fp == 0 & r$recovery$fn == 0
            cbind(id, candidates = nrow(r$recovery), exact_candidates = sum(exact),
                  adjacent_exact = adjacent_pair(exact))
        }
        scores <- if (experiment %in% c("matrix", "tensor") && r$setting %in% figure1$setting) {
            s <- r$scores[r$scores$gamma == figure1$gamma[match(r$setting, figure1$setting)],
                          c("criterion", "gamma", "k", "eta", "edges", "ebic", "status")]
            cbind(id[rep(1L, nrow(s)), ], s)
        }
        pilot <- if (!is.null(r$candidates) && "refit_status" %in% names(r$candidates)) {
            cbind(id[rep(1L, nrow(r$candidates)), names(id) != "experiment"],
                  r$candidates[c("k", "eta", "edges", "relaxation_status", "refit_status")])
        }
        list(selected = sel, paths = paths, pilot_candidates = pilot, figure1_scores = scores)
    }
    bad <- !vapply(results, function(r) r$setting %in% design$setting, TRUE)
    if (any(bad)) stop("results in ", folder, " come from settings outside ", run$settings_file)
    parts <- lapply(results, one)
    tables <- lapply(per_dataset_files, function(name) do.call(rbind, lapply(parts, `[[`, name)))

    record <- list(experiment = experiment, folder = folder, run_signature = run$run_signature,
                   commit = if (length(run$commit)) run$commit else NA_character_,
                   driver_code = run$driver_code,
                   kgmselect_fingerprint = run$kgmselect_fingerprint,
                   settings_file = run$settings_file, settings_sha256 = run$settings_sha256,
                   started = run$started, finished = run$finished,
                   datasets = length(results), failed = run$datasets$failed)
    list(tables = tables, record = record)
}

collect_outputs <- function(outputs, settings_dir) {
    present <- collected_experiments[file.exists(file.path(outputs, collected_experiments$folder,
                                                           "run.json")), ]
    if (!nrow(present)) stop("no driver results in ", outputs)

    settings_files <- unique(vapply(file.path(outputs, present$folder, "run.json"),
                                    function(f) jsonlite::read_json(f)$settings_file, ""))
    max_modes <- max(vapply(file.path(settings_dir, settings_files),
                            function(f) max(read_settings(f)$M), 0))
    figure1_file <- file.path(settings_dir, "figure1.csv")
    figure1 <- if (file.exists(figure1_file)) {
        utils::read.csv(figure1_file, stringsAsFactors = FALSE)
    }
    parts <- lapply(seq_len(nrow(present)), function(i) {
        collect_folder(present$experiment[i], outputs, present$folder[i], settings_dir,
                       max_modes, figure1)
    })
    tables <- lapply(per_dataset_files, function(name) {
        x <- do.call(rbind, lapply(parts, function(p) p$tables[[name]]))
        if (is.null(x)) return(NULL)
        keys <- intersect(c("experiment", "id", choice_columns, "k"), names(x))
        x <- x[do.call(order, c(unname(as.list(x[keys])), na.last = TRUE)), ]
        rownames(x) <- NULL
        x
    })
    list(tables = Filter(Negate(is.null), tables), sources = lapply(parts, `[[`, "record"))
}
