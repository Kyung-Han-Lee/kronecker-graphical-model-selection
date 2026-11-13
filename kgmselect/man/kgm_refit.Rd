\name{kgm_refit}
\alias{kgm_refit}
\title{Support-constrained maximum likelihood refit}
\usage{
kgm_refit(Y, supports, max_iter = 200L, tol = 1e-05)
}
\arguments{
\item{Y}{Numeric array of dimension \eqn{p_1 \times \cdots \times p_M
  \times n}; the last index runs over observations and the first index
varies fastest in the vectorisation.}

\item{supports}{List of \eqn{M} logical matrices; \code{TRUE} marks off-diagonal
entries that may be nonzero.}

\item{max_iter}{Maximum number of sweeps.}

\item{tol}{Convergence tolerance on the change of the normalised factors.}
}
\value{
List of \eqn{M} factors with attributes \code{iter} and \code{converged}.
}
\description{
Maximises the Kronecker Gaussian likelihood with the off-diagonal support of
every factor fixed, by cyclic unpenalised graphical lasso with zero
constraints. The first factor is Frobenius-normalised and its scale moved to
the second, so the Kronecker product carries the overall scale.
}
