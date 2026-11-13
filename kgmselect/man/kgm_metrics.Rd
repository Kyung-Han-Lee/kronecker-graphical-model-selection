\name{kgm_metrics}
\alias{kgm_metrics}
\title{Support recovery}
\usage{
kgm_metrics(truth, fit, tol = 1e-06)
}
\arguments{
\item{truth}{List of true precision factors.}

\item{fit}{List of estimated factors.}

\item{tol}{Estimated entries with absolute value above \code{tol} count as edges.}
}
\value{
List with \code{per_mode} (data frame of TP, FP, FN and TN by mode) and
the pooled \code{tp}, \code{fp}, \code{fn}, \code{tn}, \code{mcc}, \code{fdp} and \code{tpr}. The MCC, FDP
and TPR are zero when their denominators vanish.
}
\description{
Compares the off-diagonal supports of estimated factors with those of the
true factors over the upper triangles of all modes.
}
