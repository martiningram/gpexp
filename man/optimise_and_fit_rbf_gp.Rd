% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gp_functions.R
\name{optimise_and_fit_rbf_gp}
\alias{optimise_and_fit_rbf_gp}
\title{Optimises and fits an RBF Gaussian Process.}
\usage{
optimise_and_fit_rbf_gp(x_train, y_train, x_new, start_sigma = 10,
  start_l = 10, start_tau = 10)
}
\arguments{
\item{x_train}{The training covariates.}

\item{y_train}{The training observations.}

\item{x_new}{The covariates to predict for.}

\item{start_sigma}{The initial observation noise to use.}

\item{start_l}{The initial length scale to use.}

\item{start_tau}{The initial process standard deviation to use.}
}
\value{
A list containing entries `predictions` and `hyperparameters`,
containing the `mean` and `cov` of the process in `predictions`, and the
kernel's hyperparameters in `hyperparameters`.
}
\description{
This function optimises an RBF kernel using mAP estimation, trains the GP,
and predicts a set of new points.
}
