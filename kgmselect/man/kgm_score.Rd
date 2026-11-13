\name{kgm_score}
\alias{kgm_score}
\title{Score the candidates on a path}
\usage{
kgm_score(
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

\item{gamma}{One or more EBIC parameters \eqn{\gamma \ge 0}.}

\item{tol}{Entries with absolute value above \code{tol} count as edges.}

\item{boundary_margin, u_cap}{Margin \eqn{\vartheta_b} and cap
\eqn{U_{\max}} of the relaxation intervals.}

\item{gradient}{\code{"numeric"} (finite differences) or \code{"analytic"} gradients
in the relaxation profile.}

\item{optim_control}{Control list for \code{\link[stats:optim]{stats::optim()}} with L-BFGS-B.}

\item{refit_control}{List with \code{max_iter} and \code{tol} for the refit.}
}
\value{
Data frame with one row per grid point and value of \code{gamma}:
\code{k}, \code{eta}, \code{gamma}, \code{criterion}, \code{edges}, \code{df}, \code{loglik}, \code{tau}, \code{ebic},
\code{status} and \code{message}. A candidate whose fitted factors are not
numerically positive definite has status \code{"not_pd"} and missing scores,
except for \code{"refit"}, which uses only their support (Appendix B.6).
For \code{"refit"}, a refit that does not converge within \code{max_iter} sweeps has
status \code{"nonconverged"} and one that reaches a factor that is not positive
definite has status \code{"refit_failed"}. Any other failure, including a
relaxation profile whose optimisation fails from every starting point,
stops with an error.
}
\description{
Evaluates one selection criterion at every grid point of a fitted path:
\code{"ebic"} is EBIC at the normalised fit, \code{"scale"} the scale-profiled EBIC,
\code{"relaxation"} the relaxation-profiled EBIC and \code{"refit"} EBIC at the
support-constrained maximum likelihood estimate. The two profiled criteria
count the overall scale as one parameter, which adds the same constant to
every candidate and leaves the selection unchanged. The relaxation
coordinates are not counted as parameters (Section 3.3).
}
