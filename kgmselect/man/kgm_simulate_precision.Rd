\name{kgm_simulate_precision}
\alias{kgm_simulate_precision}
\title{Sparse precision factor for simulations}
\usage{
kgm_simulate_precision(
  p,
  family = c("er", "ba"),
  degree = 2,
  m_new = 1L,
  w = 0.3,
  w_spread = 0.25,
  min_eigen = 0.5,
  seed
)
}
\arguments{
\item{p}{Dimension.}

\item{family}{\code{"er"} or \code{"ba"} (the latter needs the igraph package).}

\item{degree}{Expected degree of the ER graph.}

\item{m_new}{Edges added per new node of the BA graph.}

\item{w, w_spread}{Centre and half-width of the edge magnitudes.}

\item{min_eigen}{Smallest eigenvalue of the factor.}

\item{seed}{Integer seed.}
}
\value{
A symmetric positive-definite \code{p} by \code{p} matrix.
}
\description{
Draws an Erdos-Renyi graph with expected degree \code{degree} (\code{"er"}) or a
Barabasi-Albert graph with \code{m_new} edges per new node (\code{"ba"}), gives each
edge a random sign and a magnitude uniform on
\eqn{[w - w_{spread}, w + w_{spread}]}, and sets a common diagonal so that
the smallest eigenvalue equals \code{min_eigen} (Section 4.1 and Appendix B of
the paper). The caller's random number state is restored afterwards.
}
