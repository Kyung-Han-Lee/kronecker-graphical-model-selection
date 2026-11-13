.with_seed <- function(seed, expr) {
    old_kind <- RNGkind()
    has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    old_seed <- if (has_seed) get(".Random.seed", envir = globalenv()) else NULL
    on.exit({
        do.call(RNGkind, as.list(old_kind))
        if (has_seed) {
            assign(".Random.seed", old_seed, envir = globalenv())
        } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
            rm(".Random.seed", envir = globalenv())
        }
    }, add = TRUE)
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
             sample.kind = "Rejection")
    expr
}
