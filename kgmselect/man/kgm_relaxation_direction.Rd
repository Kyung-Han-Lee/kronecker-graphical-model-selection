\name{kgm_relaxation_direction}
\alias{kgm_relaxation_direction}
\title{Support-preserving relaxation direction}
\usage{
kgm_relaxation_direction(
  Omega_now,
  Omega_prev,
  lambda_now,
  lambda_prev,
  tol = 1e-06,
  dir_tol = 1e-12
)
}
\arguments{
\item{Omega_now, Omega_prev}{Fitted factors at grid points \eqn{j} and
\eqn{j - 1}.}

\item{lambda_now, lambda_prev}{Their penalties, \code{lambda_now > lambda_prev}.}

\item{tol}{Entries with absolute value above \code{tol} form the support.}

\item{dir_tol}{Relative size below which the direction counts as zero.}
}
\value{
List with \code{delta}, \code{active} (whether the direction changes the
normalised factor) and \code{n_masked} (secant entries removed by the
projection).
}
\description{
Computes \eqn{\Delta_{m,j} = -\lambda_{m,j}\,\Pi_{m,j}\{(\hat\Omega_{m,j} -
\hat\Omega_{m,j-1}) / (\lambda_{m,j} - \lambda_{m,j-1})\}} for one mode,
where \eqn{\Pi_{m,j}} keeps the diagonal and the fitted off-diagonal support
(Section 3.3 of the paper).
}
