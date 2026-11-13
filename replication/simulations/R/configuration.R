parse_options <- function(args, defaults) {
    opts <- defaults
    for (a in args) {
        m <- regmatches(a, regexec("^--([a-z_]+)=(.*)$", a))[[1]]
        if (length(m) != 3L || !m[2] %in% names(defaults)) stop("unknown option: ", a)
        opts[[m[2]]] <- m[3]
    }
    opts
}

driver_config <- function(method, eta, driver_code) {
    list(method = method, eta = eta, fingerprint = kgmselect_fingerprint(),
         versions = dependency_versions(), driver_code = driver_code,
         platform = platform_info())
}

path_config <- function(criteria, gammas, eta, driver_code, support = FALSE, save_inputs = FALSE) {
    method <- list(name = "path", criteria = criteria, gammas = gammas)
    if (support) method$support <- TRUE
    if (save_inputs) method$save_inputs <- TRUE
    driver_config(method, eta, driver_code)
}

read_settings <- function(file) {
    s <- utils::read.csv(file, stringsAsFactors = FALSE)
    required <- c("setting", "M", "family", "dims", "d", "m_new", "n", "n_idx",
                  "reps", "w", "w_spread", "margin", "gradient")
    missing <- setdiff(required, names(s))
    if (length(missing)) {
        stop("settings file lacks columns: ", paste(missing, collapse = ", "))
    }
    if (anyDuplicated(s$setting)) stop("settings file repeats a setting name")
    bad <- !s$family %in% c("er", "ba") |
        (s$family == "er" & !is.finite(s$d)) |
        (s$family == "ba" & !is.finite(s$m_new)) |
        !s$gradient %in% c("numeric", "analytic") |
        vapply(s$dims, function(x) length(dims_from(x)), 0L) != s$M
    if (any(bad)) {
        stop("settings file has invalid rows: ", paste(s$setting[bad], collapse = ", "))
    }
    s
}

dims_from <- function(dims) as.integer(strsplit(dims, "x", fixed = TRUE)[[1]])

expand_datasets <- function(settings, reps = NULL) {
    rows <- lapply(seq_len(nrow(settings)), function(i) {
        r <- seq_len(if (is.null(reps)) settings$reps[i] else min(reps, settings$reps[i]))
        data.frame(id = sprintf("%s_r%03d", settings$setting[i], r),
                   setting = settings$setting[i], rep = r,
                   stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}
