\name{kgm_lambda}
\alias{kgm_lambda}
\title{Mode-wise penalties}
\usage{
kgm_lambda(p, n, eta)
}
\arguments{
\item{p}{Integer vector of mode dimensions \eqn{(p_1, \ldots, p_M)}.}

\item{n}{Sample size.}

\item{eta}{Positive multiplier.}
}
\value{
Numeric vector of the \eqn{M} penalties.
}
\description{
Computes \eqn{\lambda_m = \eta \{\log p_m / (n p_{-m})\}^{1/2}}, where
\eqn{p_{-m}} is the product of the other mode dimensions.
}
\examples{
kgm_lambda(c(116, 102), n = 289, eta = 0.0838)
}
