\name{kgm_select}
\alias{kgm_select}
\title{Select a candidate on a path}
\usage{
kgm_select(
  path,
  Y,
  criterion = c("relaxation", "scale", "ebic", "refit"),
  gamma = 0.5,
  tol = 1e-06,
  boundary_margin = 1e-08,
  u_cap = 10,
  gradient = c("numeric", "analytic"),
  optim_control = list(maxit = 100L, factr = 1e+10),
  refit_control = list(max_iter = 200L, tol = 1e-05)
)
}
\arguments{
\item{path}{A \code{kgm_path} object.}

\item{Y}{The data used to fit \code{path}.}

\item{criterion}{One of \code{"relaxation"}, \code{"scale"}, \code{"ebic"} and \code{"refit"}.}

\item{gamma}{A single EBIC parameter.}

\item{tol}{Entries with absolute value above \code{tol} count as edges.}

\item{boundary_margin, u_cap}{Margin \eqn{\vartheta_b} and cap
\eqn{U_{\max}} of the relaxation intervals.}

\item{gradient}{\code{"numeric"} (finite differences) or \code{"analytic"} gradients
in the relaxation profile.}

\item{optim_control}{Control list for \code{\link[stats:optim]{stats::optim()}} with L-BFGS-B.}

\item{refit_control}{List with \code{max_iter} and \code{tol} for the refit.}
}
\value{
Object of class \code{kgm_selection} with \code{criterion}, \code{gamma}, \code{index},
\code{eta}, \code{lambda}, \code{supports}, \code{factors}, \code{tau}, \code{u} (relaxation
coordinates, relaxation only) and \code{scores}.
}
\description{
Scores the path with \code{\link[=kgm_score]{kgm_score()}} and returns the first minimiser of the
criterion over the candidates with status \code{"ok"}; as in Appendix B.6 of
the paper, each criterion is minimised over the candidates it could score.
For relaxation
profiling the factors are the relaxed, Frobenius-normalised factors at the
selected grid point, and the precision matrix is \code{tau} times their
Kronecker product. The supports are those of the fitted factors at the
selected grid point.
}
\examples{
truth <- list(kgm_simulate_precision(6, seed = 1),
              kgm_simulate_precision(5, seed = 2))
Y <- kgm_sample(30, lapply(truth, solve), seed = 3)
path <- kgm_path(Y, kgm_eta_grid(8, c(0.01, 10)))
sel <- kgm_select(path, Y, criterion = "relaxation", gamma = 0.5)
sel$eta
}
