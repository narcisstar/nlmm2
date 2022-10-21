#' @name SSMNormal
#' @aliases dssmn dmssmn rssmn rmssmn
#' @title (Multivariate) Skew Scale Mixture of Normal Distributions
#' @description Random deviate generator for (Multivariate) Skew Scale Mixture of Normal Distributions
#'
#' @usage
#'
#' dssmn (x, mu = 0, omegabar = 1, delta = 0, smixvar = "Gamma",
#'        df = switch(smixvar, Bin = c(.5, .5), 5), log = FALSE)
#'
#' rssmn (n, mu = 0, omegabar = 1, delta = 0, smixvar = "Gamma",
#'        df = switch(smixvar, Bin = c(.5, .5), 5))
#'
#' dmssmn (x, mu = rep(0, q), Omegabar = diag(q), delta = rep(0, q),
#'         smixvar = "Gamma", df = switch(smixvar, Bin = c(.5, .5), 5),
#'         log = FALSE)
#'
#' rmssmn (n, mu = rep(0, q), Omegabar = diag(q), delta = rep(0, q),
#'         smixvar = "Gamma", df = switch(smixvar, Bin = c(.5, .5), 5))
#'
#' @param x vector or matrix of quantiles. If \code{x} is a matrix, each row is taken to be a quantile.
#' @param n integer, number of deviates to be generated.
#' @param mu numeric, location vector.
#' @param Omegabar,omegabar numeric, scale matrix. \code{omegabar} is a scalar, for the univariate case.
#' @param delta numeric, shape/skewness vector.
#' @param smixvar character, indicates the distribution of the scale mixing variable of the target skew scale mixture of normal (ssmn) distribution. "Gamma" gives the (skew) t, "Beta" gives (skew) slash, "Bin" gives the (skew) contaminated normal, and any other string gives the (skew) normal.
#' @param df numeric, degrees of freedom parameter. It is ignored for (skew) normal distributions. Otherwise, \code{df} is a scalar for all non (skew) normal distributions, except, \code{df} is a vector of two elements each in the open (0, 1) when \code{smixvar = "Bin"}.
#'  Conditions for the first order moment of the ssmn distribution to be defined: \code{df > 1} if \code{smixvar = "Gamma"}, \code{df > .5} if \code{smixvar = "Beta"}.
#'  Conditions for the second order moment of the ssmn distribution to be defined: \code{df > 2} if \code{smixvar = "Gamma"}, \code{df > 1} if \code{smixvar = "Beta"}.
#'
#' @param log logical; if \code{TRUE}, densities \code{d} are given as \code{log(d)}.
#
#' @details A ssmn variable/vector \code{X} is represented as:
#'
#' \code{X = mu + T_0 * delta / sqrt(kappa) + z / sqrt(kappa)} where \code{T_0} is a half normal variable, \code{kappa} is a scale mixing variable, and \code{z} is a normal variable/vector with null mean and covariance matrix \code{Omegabar}.
#'
# Generate variates from the skew scale mixture of normal distribution
# smixvar = distribution of the mixture variable
rmssmn <- function (n, mu = rep(0, q), Omegabar = diag(q),
                    delta = rep(0, q), smixvar = "Gamma",
                    df = switch(smixvar, Bin = c(.5, .5), 5)) {
  if (identical(smixvar, "Bin") && length(df) < 2)
    stop("the (Skew) contaminated model requires a length 2 'df' argument")
  if (all(c(missing(mu), missing(Omegabar), missing(delta)))) {
    mu <- delta <- 0
    q <- Omegabar <- 1
  }
  else {
    if (!missing(mu))
      q <- length(mu)
    else {
      if (!missing(Omegabar))
        q <- NROW(Omegabar)
      else {
        if (!missing(delta))
          q <- length(delta)
      }
    }
  }
  if (q == 1)
    return(rssmn(n, mu = mu, omegabar = sqrt(Omegabar[1]), delta = delta,
                 smixvar = smixvar, df = df))
  if (identical(smixvar, "Bin")) {
    kappa <- rbinom(n, size = 1, prob = df[1]) * df[2]
    zkappa <- kappa == 0
    if (any(zkappa))
      kappa[zkappa] <- 1
  }
  else {
    kappa <- switch(smixvar,
                    Gamma = rgamma(n, shape = df/2, rate = df/2),
                    Beta = rbeta(n, shape1 = df, shape2 = 1),
                    rep(1, n))
  }
  kappa <- 1/sqrt(kappa)
  z <- rmvnorm (n, Var = Omegabar) * replicate(q, kappa)
  if (any(delta != 0)) {
    loc <- outer(delta, abs(rnorm(n)) * kappa)
    if (any(mu != 0)) {
      if (NROW(mu) == n)
        loc <- loc + t(mu)
      else
        loc <- loc + rep(mu, length.out = q)
    }
    res <-  t(loc) + z
  }
  else {
    if (any(mu != 0)) {
      if (nrow(mu) == n)
        loc <- mu
      else
        loc <- t(replicate(n, rep(mu, length.out = q)))
      res <- loc + z
    }
    else
      res <- z
  }
  res
}

# Generate variates from the univariate skew scale mixture of normal distribution
rssmn <- function (n, mu = 0, omegabar = 1, delta = 0, smixvar = "Gamma",
                   df = switch(smixvar, Bin = c(.5, .5), 5)) {
  if (identical(smixvar, "Bin") && length(df) < 2)
    stop("the (Skew) contaminated model requires a length 2 'df' argument")
  if (identical(smixvar, "Bin")) {
    kappa <- rbinom(n, size = 1, prob = df[1]) * df[2]
    zkappa <- kappa == 0
    if (any(zkappa))
      kappa[zkappa] <- 1
  }
  else {
    kappa <- switch(smixvar,
                    Gamma = rgamma(n, shape = df/2, rate = df/2),
                    Beta = rbeta(n, shape1 = df, shape2 = 1),
                    rep(1, n))
  }
  kappa <- 1/sqrt(kappa)
  z <- rnorm (n) * omegabar * kappa
  if (any(delta != 0)) {
    loc <- rep(delta, length.out = n) * abs(rnorm(n)) * kappa
    if (any(mu != 0))
      loc <- loc + rep(mu, length.out = n)
    res <-  loc + z
  }
  else {
    if (any(mu != 0)) {
      loc <- rep(mu, length.out = n)
      res <- loc + z
    }
    else
      res <- z
  }
  res
}
