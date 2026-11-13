kgm_select <- function(path, Y,
                       criterion = c("relaxation", "scale", "ebic", "refit"),
                       gamma = 0.5, tol = 1e-6, boundary_margin = 1e-8,
                       u_cap = 10, gradient = c("numeric", "analytic"),
                       optim_control = list(maxit = 100L, factr = 1e10),
                       refit_control = list(max_iter = 200L, tol = 1e-5)) {
    criterion <- match.arg(criterion)
    gradient <- match.arg(gradient)
    if (length(gamma) != 1L) stop("kgm_select takes a single value of gamma")
    scores <- kgm_score(path, Y, criterion = criterion, gamma = gamma,
                        tol = tol, boundary_margin = boundary_margin,
                        u_cap = u_cap, gradient = gradient,
                        optim_control = optim_control,
                        refit_control = refit_control)
    usable <- scores$status == "ok"
    if (any(!is.finite(scores$ebic[usable]))) {
        stop("a scored candidate has a non-finite criterion value")
    }
    if (!any(usable)) stop("no candidate on the path could be scored")
    k <- scores$k[usable][which.min(scores$ebic[usable])]
    fit <- path$fits[[k]]
    supports <- lapply(fit, function(A) {
        S <- abs(A) > tol
        diag(S) <- FALSE
        S
    })
    u <- NULL
    if (criterion == "relaxation") {
        prof <- .profile_relaxation(Y, fit, .relaxation_deltas(path, k, tol),
                                    boundary_margin = boundary_margin,
                                    u_cap = u_cap, optim_control = optim_control,
                                    gradient = gradient)
        factors <- prof$Omega_list
        tau <- prof$tau_global
        u <- prof$u
    } else if (criterion == "refit") {
        factors <- lapply(.refit_support_mle(Y, supports,
                                             max_iter = refit_control$max_iter,
                                             tol = refit_control$tol), identity)
        tau <- 1
    } else {
        factors <- fit
        tau <- scores$tau[scores$k == k]
    }
    structure(list(criterion = criterion, gamma = gamma, index = k,
                   eta = path$eta[k], lambda = path$lambda[[k]],
                   supports = supports, factors = factors, tau = tau, u = u,
                   scores = scores),
              class = "kgm_selection")
}
