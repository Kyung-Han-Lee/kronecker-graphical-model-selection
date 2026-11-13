reference_for_installed_huge <- function() {
    installed <- as.character(utils::packageVersion("huge"))
    files <- list.files("fixtures", pattern = "^golden.*[.]rds$", full.names = TRUE)
    for (f in files) {
        g <- readRDS(f)
        if (identical(unname(g$meta$packages[["huge"]]), installed)) return(g)
    }
    NULL
}

support_of <- function(fit, tol = 1e-6) {
    lapply(fit, function(A) abs(A[upper.tri(A)]) > tol)
}
