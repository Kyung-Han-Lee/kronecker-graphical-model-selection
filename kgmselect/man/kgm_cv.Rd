\name{kgm_cv}
\alias{kgm_cv}
\title{Cross-validation on the normalised path}
\usage{
kgm_cv(
  Y,
  eta = kgm_eta_grid(),
  nfolds = 5L,
  seed = 2026L,
  n_sweep = 2L,
  fold_id = NULL
)
}
\arguments{
\item{Y}{Numeric array of dimension \eqn{p_1 \times \cdots \times p_M
  \times n}; the last index runs over observations and the first index
varies fastest in the vectorisation.}

\item{eta}{Increasing grid of penalty multipliers, see \code{\link[=kgm_eta_grid]{kgm_eta_grid()}}.}

\item{nfolds}{Number of folds, from 3 to \eqn{n/2}.}

\item{seed}{Seed for the fold labels
\code{sample(rep(seq_len(nfolds), length.out = n))}.}

\item{n_sweep}{Number of sweeps \eqn{T}.}

\item{fold_id}{Optional fold labels; overrides \code{nfolds} and \code{seed}.}
}
\value{
Object of class \code{kgm_cv} with \code{eta}, \code{index}, \code{lambda}, \code{factors},
\code{tau}, \code{fold_id}, \code{fold_loglik}, \code{fold_tau}, \code{loglik} (fold averages) and
\code{eta_grid}.
}
\description{
K-fold cross-validation as in Section 4.2 of the paper. On every training
set the estimator is refitted with penalties computed from the training
sample size, the overall scale is profiled on the training data only, and
the held-out Gaussian log-likelihood is averaged over folds. The grid value
with the largest average is refitted on all observations.
}
