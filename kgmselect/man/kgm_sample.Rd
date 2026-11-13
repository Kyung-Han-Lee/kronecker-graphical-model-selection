\name{kgm_sample}
\alias{kgm_sample}
\title{Tensor-normal sample}
\usage{
kgm_sample(n, Sigma, seed = NULL)
}
\arguments{
\item{n}{Number of observations.}

\item{Sigma}{List of \eqn{M \ge 2} covariance factors.}

\item{seed}{Optional integer seed.}
}
\value{
Array of dimension \eqn{p_1 \times \cdots \times p_M \times n}.
}
\description{
Draws \code{n} mean-zero observations whose vectorisation has covariance
\eqn{\Sigma_M \otimes \cdots \otimes \Sigma_1}, without forming the full
Kronecker matrix. With a seed, the caller's random number state is restored
afterwards.
}
