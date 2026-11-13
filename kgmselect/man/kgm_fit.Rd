\name{kgm_fit}
\alias{kgm_fit}
\title{Normalised cyclic estimator at one penalty}
\usage{
kgm_fit(Y, lambda, n_sweep = 2L)
}
\arguments{
\item{Y}{Numeric array of dimension \eqn{p_1 \times \cdots \times p_M
  \times n}; the last index runs over observations and the first index
varies fastest in the vectorisation.}

\item{lambda}{Numeric vector of the \eqn{M} mode-wise penalties.}

\item{n_sweep}{Number of sweeps \eqn{T}.}
}
\value{
List of \eqn{M} unit-Frobenius precision factors.
}
\description{
Runs \code{n_sweep} sweeps of mode-wise graphical lasso, divides each solution by
its Frobenius norm and starts from identity factors (Section 2.2 and
Algorithm 1 of the paper). The penalty includes the diagonal.
}
