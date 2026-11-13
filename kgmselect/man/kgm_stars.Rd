\name{kgm_stars}
\alias{kgm_stars}
\title{Stability approach to regularisation selection}
\usage{
kgm_stars(path, Y, N = 20L, threshold = 0.05, seed, tol = 1e-06)
}
\arguments{
\item{path}{A \code{kgm_path} object.}

\item{Y}{The data used to fit \code{path}.}

\item{N}{Number of subsamples.}

\item{threshold}{One or more instability thresholds.}

\item{seed}{Integer seed for the subsamples.}

\item{tol}{Entries with absolute value above \code{tol} count as edges.}
}
\value{
List with \code{index} and \code{eta} (one per threshold), \code{threshold},
\code{instability} (one value per grid point) and \code{subsample_size}.
}
\description{
StARS on the normalised path (Liu, Roeder and Wasserman, 2010). The
estimator is refitted at every grid point on \code{N} subsamples drawn without
replacement, with the penalties of the full-data path. The instability of a
grid point is the mean of \eqn{2\theta(1 - \theta)} over the
upper-triangular entries of all modes, where \eqn{\theta} is the edge
frequency across subsamples. The selected grid point is the smallest
penalty whose instability, and that of every larger penalty, is at most
\code{threshold}. Subsample sizes follow \code{huge.select()} in huge 1.6:
\eqn{\lfloor 0.8 n \rfloor} for \eqn{n \le 144} and
\eqn{\lfloor 10 \sqrt n \rfloor} otherwise.
}
