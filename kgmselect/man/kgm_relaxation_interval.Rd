\name{kgm_relaxation_interval}
\alias{kgm_relaxation_interval}
\title{Admissible relaxation interval for one mode}
\usage{
kgm_relaxation_interval(Omega, delta, boundary_margin = 1e-08, u_cap = 10)
}
\arguments{
\item{Omega}{Fitted factor.}

\item{delta}{Relaxation direction from \code{\link[=kgm_relaxation_direction]{kgm_relaxation_direction()}}.}

\item{boundary_margin, u_cap}{Margin \eqn{\vartheta_b} and cap
\eqn{U_{\max}}.}
}
\value{
List with \code{active}, \code{u_pd}, \code{u_support}, \code{u_upper} (right end of
the interval), \code{bound_kind} and quantities reused by the profile.
}
\description{
Returns the interval on which \eqn{\hat\Omega_{m,j} + u\Delta_{m,j}} stays
positive definite and keeps its support and active signs, shrunk by the
boundary margin and capped (Section 3.3 of the paper).
}
