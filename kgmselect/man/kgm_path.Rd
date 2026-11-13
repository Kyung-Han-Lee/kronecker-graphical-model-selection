\name{kgm_path}
\alias{kgm_path}
\title{Normalised regularisation path}
\usage{
kgm_path(Y, eta = kgm_eta_grid(), n_sweep = 2L)
}
\arguments{
\item{Y}{Numeric array of dimension \eqn{p_1 \times \cdots \times p_M
  \times n}; the last index runs over observations and the first index
varies fastest in the vectorisation.}

\item{eta}{Increasing grid of penalty multipliers, see \code{\link[=kgm_eta_grid]{kgm_eta_grid()}}.}

\item{n_sweep}{Number of sweeps \eqn{T}.}
}
\value{
Object of class \code{kgm_path}: a list with \code{eta}, \code{lambda} (list of
penalty vectors), \code{fits} (list of factor lists), \code{n_sweep}, \code{dims} and \code{n}.
}
\description{
Fits the normalised cyclic estimator independently at every value of the
grid, with mode-wise penalties from \code{\link[=kgm_lambda]{kgm_lambda()}}.
}
\examples{
truth <- list(kgm_simulate_precision(6, seed = 1),
              kgm_simulate_precision(5, seed = 2))
Y <- kgm_sample(30, lapply(truth, solve), seed = 3)
path <- kgm_path(Y, kgm_eta_grid(8, c(0.01, 10)))
length(path$fits)
}
