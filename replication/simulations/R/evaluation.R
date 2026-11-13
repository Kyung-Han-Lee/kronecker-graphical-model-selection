selection_tables <- function(scores, path, truth, tol = 0) {
    keys <- unique(scores[, c("criterion", "gamma")])
    pooled <- vector("list", nrow(keys))
    modes <- vector("list", nrow(keys))
    for (i in seq_len(nrow(keys))) {
        s <- scores[scores$criterion == keys$criterion[i] &
                        scores$gamma == keys$gamma[i], ]
        k <- select_index(s)
        m <- kgm_metrics(truth, path$fits[[k]], tol = tol)
        edges <- s$edges[s$k == k]
        pooled[[i]] <- data.frame(
            criterion = keys$criterion[i], gamma = keys$gamma[i], k = k,
            eta = path$eta[k], edges = edges, empty = edges == 0L,
            tp = m$tp, fp = m$fp, fn = m$fn, tn = m$tn, mcc = m$mcc,
            fdp = m$fdp, tpr = m$tpr, stringsAsFactors = FALSE)
        modes[[i]] <- cbind(data.frame(criterion = keys$criterion[i],
                                       gamma = keys$gamma[i], k = k,
                                       stringsAsFactors = FALSE),
                            m$per_mode)
    }
    list(selected = do.call(rbind, pooled), selected_modes = do.call(rbind, modes))
}

fit_tables <- function(label, k, eta, fit, truth, tol = 0) {
    m <- kgm_metrics(truth, fit, tol = tol)
    edges <- m$tp + m$fp
    list(selected = cbind(label, data.frame(k = k, eta = eta, edges = edges,
                                            empty = edges == 0L, tp = m$tp,
                                            fp = m$fp, fn = m$fn, tn = m$tn,
                                            mcc = m$mcc, fdp = m$fdp, tpr = m$tpr)),
         selected_modes = cbind(label, data.frame(k = k), m$per_mode))
}

support_records <- function(path, tol = 0) {
    do.call(rbind, lapply(seq_along(path$fits), function(k) {
        fit <- path$fits[[k]]
        mask <- lapply(fit, function(A) {
            keep <- abs(A) > tol
            diag(keep) <- FALSE
            keep
        })
        exact <- lapply(fit, function(A) {
            keep <- A != 0
            diag(keep) <- FALSE
            keep
        })
        values <- unlist(lapply(fit, function(A) abs(A[upper.tri(A)])))
        tiny <- values[values > 0 & values <= tol]
        data.frame(k = k, edges = sum(values > tol), exact_edges = sum(values > 0),
                   tiny_edges = length(tiny), max_tiny = max(c(0, tiny)),
                   pattern_hash = digest::digest(mask, algo = "sha256"),
                   exact_hash = digest::digest(exact, algo = "sha256"),
                   stringsAsFactors = FALSE)
    }))
}

candidate_recovery <- function(path, truth, tol = 0) {
    do.call(rbind, lapply(seq_along(path$fits), function(k) {
        m <- kgm_metrics(truth, path$fits[[k]], tol = tol)
        data.frame(k = k, tp = m$tp, fp = m$fp, fn = m$fn, tn = m$tn)
    }))
}
