.mcc <- function(tp, fp, fn, tn) {
    tp <- as.numeric(tp)
    fp <- as.numeric(fp)
    fn <- as.numeric(fn)
    tn <- as.numeric(tn)
    den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    ifelse(den == 0, 0, (tp * tn - fp * fn) / den)
}

.confusion_columns <- function(truth, est, tol) {
    out <- list()
    for (m in seq_along(truth)) {
        Tm <- truth[[m]]
        Em <- est[[m]]
        iu <- upper.tri(Tm)
        tv <- Tm[iu]
        ev <- Em[iu]
        t_nz <- tv != 0
        e_nz <- abs(ev) > tol
        out[[sprintf("tp%d", m)]] <- sum(t_nz & e_nz)
        out[[sprintf("fp%d", m)]] <- sum(!t_nz & e_nz)
        out[[sprintf("fn%d", m)]] <- sum(t_nz & !e_nz)
        out[[sprintf("tn%d", m)]] <- sum(!t_nz & !e_nz)
        out[[sprintf("tps%d", m)]] <- sum(t_nz & e_nz & sign(tv) == sign(ev))
    }
    out
}

kgm_metrics <- function(truth, fit, tol = 1e-6) {
    if (length(truth) != length(fit)) {
        stop("truth and fit must have the same number of modes")
    }
    cc <- .confusion_columns(truth, fit, tol)
    M <- length(truth)
    counts <- function(key) {
        vapply(seq_len(M), function(m) as.numeric(cc[[paste0(key, m)]]), numeric(1))
    }
    per_mode <- data.frame(mode = seq_len(M), tp = counts("tp"),
                           fp = counts("fp"), fn = counts("fn"),
                           tn = counts("tn"))
    tp <- sum(per_mode$tp)
    fp <- sum(per_mode$fp)
    fn <- sum(per_mode$fn)
    tn <- sum(per_mode$tn)
    list(per_mode = per_mode, tp = tp, fp = fp, fn = fn, tn = tn,
         mcc = .mcc(tp, fp, fn, tn),
         fdp = if (tp + fp > 0) fp / (tp + fp) else 0,
         tpr = if (tp + fn > 0) tp / (tp + fn) else 0)
}
