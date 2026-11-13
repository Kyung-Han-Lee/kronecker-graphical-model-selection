prepare_refit_convergence <- function(outputs, inputs, jobs_file, destination) {
    if (file.exists(destination)) stop("output already exists")
    paths <- read.csv(file.path(inputs, "paths.csv"), stringsAsFactors = FALSE)
    jobs <- read.csv(jobs_file, stringsAsFactors = FALSE)
    required <- c("id", "k")
    if (!all(required %in% names(jobs)) || !nrow(jobs)) stop("invalid convergence jobs")
    results <- vector("list", nrow(jobs))
    for (i in seq_len(nrow(jobs))) {
        metadata <- paths[paths$id == jobs$id[i], ]
        if (nrow(metadata) != 1L) stop("unknown convergence dataset")
        result <- readRDS(file.path(outputs, metadata$cache))
        generated <- cached_observations(result)
        k <- as.integer(jobs$k[i])
        if (k < 1L || k > length(result$inputs$path$eta)) stop("candidate index is out of range")
        fit <- result$inputs$path$fits[[k]]
        supports <- lapply(fit, function(A) {
            keep <- A != 0
            diag(keep) <- FALSE
            keep
        })
        strict <- kgm_refit(generated$Y, supports, max_iter = 2000L, tol = 1e-8)
        if (!isTRUE(attr(strict, "converged"))) stop("strict refit did not converge")
        saved <- result$scores[result$scores$criterion == "refit" &
            result$scores$gamma == 0.5 & result$scores$k == k, ]
        if (nrow(saved) != 1L || saved$status != "ok") stop("saved Refit score is unavailable")
        ll <- kgmselect:::.compute_loglik(generated$Y, strict,
            sqrt_method = if (length(strict) == 2L) "none" else "eigen")
        results[[i]] <- data.frame(id = jobs$id[i], k = k,
            saved_refit_ll = saved$loglik, strict_refit_ll = as.numeric(ll),
            iterations = attr(strict, "iter"), stringsAsFactors = FALSE)
    }
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    write_exact_csv(do.call(rbind, results), destination)
}
