\name{kgm_eta_grid}
\alias{kgm_eta_grid}
\title{Grid of penalty multipliers}
\usage{
kgm_eta_grid(J = 40L, range = c(1e-04, 100))
}
\arguments{
\item{J}{Number of grid points.}

\item{range}{Smallest and largest value.}
}
\value{
Increasing numeric vector of length \code{J}.
}
\description{
Geometrically spaced values of the common multiplier \eqn{\eta}. The default
is the 40-point grid on \eqn{[10^{-4}, 10^2]} used in the paper.
}
\examples{
kgm_eta_grid(5, c(0.01, 1))
}
