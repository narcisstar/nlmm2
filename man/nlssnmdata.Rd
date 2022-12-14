% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fakessnmdata.R
\name{nlssnmdata}
\alias{nlssnmdata}
\alias{fake.ssnmdata}
\title{Generate data from a non linear mixed effects regression (nlmer) model under skew scale mixture of normal (ssmn) distributions.}
\usage{
nlssnmdata(
  n = 10,
  n_i = 3,
  n_x = 2,
  n_z,
  n_w = 0,
  nlfun = function(phi, X) {     exp(phi[1] + phi[2] * X[1] + phi[3] * log(X[2])) },
  Xdist = rep("norm", n_x),
  A_i = diag(3),
  Zdist = rep("unif", n_z),
  Zposition = NULL,
  B_i = diag(3),
  Wdist = rep("unif", n_w),
  Wposition = NULL,
  fixef = c(2, -0.5, 1),
  omegabar2 = 1,
  deltao = 2,
  Omegabar_b = diag(q),
  delta_b = rep(0, q),
  smixvar = "Gamma",
  smixvar_b = smixvar,
  nu = switch(smixvar, Bin = c(0.5, 0.5), Gamma = 5, Beta = 5, NULL),
  nu_b = switch(smixvar_b, Bin = c(0.5, 0.5), Gamma = 5, Beta = 5, NULL)
)
}
\arguments{
\item{n}{number of clusters/subjects.}

\item{n_i}{number of observations/measurements per subject.}

\item{n_x}{number of covariates entering the non linear function along with the mixed effects (non linear parameters).}

\item{n_z}{number of covariates with fixed slopes.}

\item{n_w}{number of covariates with random slopes.}

\item{nlfun}{non linear function. Its first two arguments must be (named) \code{phi} (mixed effects) and \code{X} (vector of covariates)}

\item{Xdist}{character vector of length \code{n_x}. The \code{i}-th element specifies the distribution to used to simulate the \code{i}-th element of the vector of covariates entering the non linear function. Each character "\code{name}" must be such that "\code{rname(n)}" generate \code{n} elements from the specified distribution. For exemple, "\code{norm}" corresponds to the standard normal distribution, "\code{mydist}" where \code{mydist = function(n) rpois (n, lambda = 5)} corresponds to a Poisson distribution with mean \code{5}.}

\item{A_i}{numeric \code{r} by \code{p} matrix, common fixed design matrix for all clusters (\code{p} = number of fixed effects).}

\item{Zdist}{character vector of length \code{n_z}. The \code{i}-th element specifies the distribution to used to simulate the \code{i}-th fixed covariate.}

\item{Zposition}{logical matrix of the same size as \code{A_i}. The elements of \code{Zposition} indicate which elements of \code{A_i} are fixed covariates (so the corresponding random effect is a slope). Accordingly, the equality \code{sum(Zposition) == n_z} must hold.}

\item{B_i}{numeric \code{r} by \code{q} matrix, common random design matrix for all clusters (\code{q} = number of random effects).}

\item{Wdist}{character vector of length \code{n_w}. The \code{i}-th element specifies the distribution to used to simulate the \code{i}-th random covariate.}

\item{Wposition}{logical matrix of the same size as \code{B_i}. The elements of \code{Wposition} indicate which elements of \code{B_i} are random covariates (so the corresponding random effect is a slope). Accordingly, the equality \code{sum(Wposition) == n_w} must hold.}

\item{fixef}{numeric \code{p} vector of fixed effects.}

\item{omegabar2, deltao}{numeric scalars (\code{omegabar2 > 0}), residual scale and shape parameters.}

\item{Omegabar_b, delta_b}{numeric \code{q} by \code{q} positive definite matrix and \code{q} vector, scale and shape parameters of random effects.}

\item{smixvar, smixvar_b}{characters, specify the mixing distributions for residuals and random effetcs (see \code{smixvar } at \link{rmssmn}).}

\item{nu, nu_b}{degrees of freedom parameters for residuals and random effetcs (see \code{df} at \link{rmssmn}).}
}
\description{
The function generate data for testing functions related to the fit of \code{nlmer} models.
}
