dataset_seeds <- function(row, rep) {
    list(graph = 1000L * seq_len(row$M) + as.integer(rep),
         noise = 5000L + 10L * as.integer(rep) + as.integer(row$n_idx))
}

make_dataset <- function(row, rep) {
    p <- dims_from(row$dims)
    seeds <- dataset_seeds(row, rep)
    truth <- lapply(seq_len(row$M), function(m) {
        if (row$family == "er") {
            kgm_simulate_precision(p[m], "er", degree = row$d, w = row$w,
                                   w_spread = row$w_spread,
                                   min_eigen = row$margin, seed = seeds$graph[m])
        } else {
            kgm_simulate_precision(p[m], "ba", m_new = as.integer(row$m_new),
                                   w = row$w, w_spread = row$w_spread,
                                   min_eigen = row$margin, seed = seeds$graph[m])
        }
    })
    list(truth = truth,
         Y = kgm_sample(row$n, lapply(truth, solve), seed = seeds$noise),
         seeds = seeds)
}
