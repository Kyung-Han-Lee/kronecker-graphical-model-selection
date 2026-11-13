score_path <- function(path, Y, criteria, gammas, gradient) {
    do.call(rbind, lapply(criteria, function(criterion) {
        kgm_score(path, Y, criterion = criterion, gamma = gammas,
                  gradient = gradient, tol = 0)
    }))
}

select_index <- function(scores) {
    usable <- scores$status == "ok"
    if (any(!is.finite(scores$ebic[usable]))) {
        stop("a scored candidate has a non-finite criterion value")
    }
    if (!any(usable)) stop("no candidate on the path could be scored")
    scores$k[usable][which.min(scores$ebic[usable])]
}
