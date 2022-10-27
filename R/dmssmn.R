# Density of a multivariate skew scale mixture of normal distribution

# not OK: problem with the skewness parameter delta when it has
# more than one non zero component

#' @name SSMNormal
#' @aliases dssmn dmssmn rssmn rmssmn
#' @title (Multivariate) Skew Scale Mixture of Normal Distributions
#' @description Random deviate generator for (Multivariate) Skew Scale Mixture of Normal Distributions
#'
#' @usage
#'
#' dssmn (x, mu = 0, omegabar = 1, delta = 0, mixvar = "Gamma",
#'        df = switch(mixvar, Bin = c(.5, .5), 5), log = FALSE)
#'
#' rssmn (n, mu = 0, omegabar = 1, delta = 0, mixvar = "Gamma",
#'        df = switch(mixvar, Bin = c(.5, .5), 5))
#'
#' dmssmn (x, mu = rep(0, q), Omegabar = diag(q), delta = rep(0, q),
#'         mixvar = "Gamma", df = switch(mixvar, Bin = c(.5, .5), 5),
#'         log = FALSE)
#'
#' rmssmn (n, mu = rep(0, q), Omegabar = diag(q), delta = rep(0, q),
#'         mixvar = "Gamma", df = switch(mixvar, Bin = c(.5, .5), 5))
#'
#' @param x vector or matrix of quantiles. If \code{x} is a matrix, each row is taken to be a quantile.
#' @param n integer, number of deviates to be generated.
#' @param mu numeric, location vector.
#' @param Omegabar,omegabar numeric, scale matrix. \code{omegabar} is a scalar, for the univariate case.
#' @param delta numeric, shape/skewness vector.
#' @param mixvar character, indicates the distribution of the scale mixing variable of the target skew scale mixture of normal (ssmn) distribution. "Gamma" gives the (skew) t, "Beta" gives (skew) slash, "Bin" gives the (skew) contaminated normal, and any other string gives the (skew) normal.
#' @param df numeric, degrees of freedom parameter. It is ignored for (skew) normal distributions. Otherwise, \code{df} is a scalar for all non (skew) normal distributions, except, \code{df} is a vector of two elements each in the open (0, 1) when \code{mixvar = "Bin"}.
#'  Conditions for the first order moment of the ssmn distribution to be defined: \code{df > 1} if \code{mixvar = "Gamma"}, \code{df > .5} if \code{mixvar = "Beta"}.
#'  Conditions for the second order moment of the ssmn distribution to be defined: \code{df > 2} if \code{mixvar = "Gamma"}, \code{df > 1} if \code{mixvar = "Beta"}.
#'
#' @param log logical; if \code{TRUE}, densities \code{d} are given as \code{log(d)}.
#
#' @details A ssmn variable/vector \code{X} is represented as:
#'
#' \code{X = mu + T_0 * delta / sqrt(kappa) + z / sqrt(kappa)} where \code{T_0} is a half normal variable, \code{kappa} is a scale mixing variable, and \code{z} is a normal variable/vector with null mean and covariance matrix \code{Omegabar}.
#'
# Generate variates from the skew scale mixture of normal distribution
# mixvar = distribution of the mixture variable
#' @export
dmssmn <- function (x, mu = rep(0, q), Omegabar = diag(q), delta = rep(0, q),
                    mixvar = "Gamma", df = switch(mixvar, Bin = c(.5, .5), 5),
                    log = FALSE) {
  if (identical(mixvar, "Bin") && length(df) < 2)
    stop("the (Skew) contaminated model requires a length 2 'df' argument")
  if (all(c(missing(mu), missing(Omegabar), missing(delta)))) {
    q <- ncol(x)
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
  if (q == 1) {
    return(dssmn(x, mu = mu, omegabar = sqrt(Omegabar[1]), delta = delta,
                 mixvar = mixvar, df = df, log = log))
  }
  if (is.vector(x))
    x <- matrix(x, ncol = q)
  if (q != ncol(x))
    stop("'x' and 'mu/delta/Omegabar' have non-conforming sizes")
  if (!isSymmetric(Omegabar, tol = sqrt(.Machine$double.eps),
                   check.attributes = FALSE)) {
    stop("Omegabar must be a symmetric matrix")
  }
  if (identical(mixvar, "Beta")) {
    Omega <- Omegabar + outer(delta, delta)
    dec <- tryCatch(chol(Omega), error = function(e) e)
    if (inherits(dec, "error")) {
      mean <- mu + sqrt(2/pi) * delta
      x.is.mu <- colSums(t(x) != mean) == 0
      logpdf <- rep.int(-Inf, nrow(x))
      logpdf[x.is.mu] <- Inf
      names(logpdf) <- rownames(x)
      return(if (log)
        logpdf
        else
          exp(logpdf))
    }
    z <- backsolve(dec, t(x) - mu, transpose = TRUE)
    lambdaz <- c(backsolve(dec, delta, transpose = TRUE))
    den <- sum(lambdaz * lambdaz)
    lambdaz <- c(t(lambdaz) %*% z) / sqrt(1 - den)
    Q <- 0.5 * colSums(z^2)
    lambdaz <- cbind(df, Q, lambdaz)
    df <- lambdaz[, 1]
    Q <- lambdaz[, 2]
    lambdaz <- lambdaz[, 3]
    infinitex <- rowSums(is.infinite(x)) > 0
    ok_0 <- (is.na(lambdaz) + is.infinite(lambdaz) +
               (df <= 0) + infinitex) == 0
    ok <- (ok_0 + is.finite(df)) == 2
    logpdf <- numeric(length(lambdaz))
    logpdf[!ok] <- NaN
    const <- - sum(log(diag(dec))) - 0.5 * q * log(2 * pi)
    if (any(ok)) {
      fdsn <- function(u) {
        exp (pnorm(lambdaz[ok] * sqrt(u), log.p = TRUE) +
               const + (df[ok]- 1 + .5 * q) * logb(u) - Q[ok] * u)
      }
      logpdf[ok] <- logb(2) + logb(df[ok]) + logb(
        cubature::adaptIntegrate(f = fdsn, lowerLimit = 0,
                                 upperLimit = 1,
                                 fDim = sum(ok))$integral)
    }
    ok <- (ok_0 + is.infinite(df)) == 2
    if (any(ok)) {
      logpdf[ok] <- pnorm(lambdaz[ok], log.p = TRUE) + const - Q[ok] + logb(2)
    }
    logpdf <- replace(logpdf, infinitex, -Inf)
    names(logpdf) <- rownames(x)
    return(if (log)
      logpdf
      else exp(logpdf))
  }
  else if (identical(mixvar, "Gamma")) {
    n <- dim(x)[1]
    p <- dim(x)[2]
    res <- EMMIXskew::ddmst(dat = x, n = n, p = p, mean = mu, cov = Omegabar, del = delta, nu = df)
    return(res)
    }
  else if (identical(mixvar, "Bin")) {
    if (length(df) > 2) {
      if (is.matrix(df)) {
        df1 <- df[,1]
        df2 <- df[,2]
      }
      else {
        ndf <- length(df)
        hndf <- floor(ndf/2)
        df1 <- df[1:hndf]
        df2 <- df[(hndf+1):ndf]
      }
    }
    else {
      df1 <- df[1]
      df2 <- df[2]
    }
    Omega <- Omegabar + outer(delta, delta)
    dec <- tryCatch(chol(Omega), error = function(e) e)
    if (inherits(dec, "error")) {
      mean <- mu + sqrt(2/pi) * delta * (1 + df1 * ((1/df2) - 1))
      x.is.mu <- colSums(t(x) != mean) == 0
      logpdf <- rep.int(-Inf, nrow(x))
      logpdf[x.is.mu] <- Inf
      names(logpdf) <- rownames(x)
      return(if (log)
        logpdf
        else
          exp(logpdf))
    }
    z <- backsolve(dec, t(x) - mu, transpose = TRUE)
    lambdaz <- c(backsolve(dec, delta, transpose = TRUE))
    den <- sum(lambdaz * lambdaz)
    lambdaz <- c(t(lambdaz) %*% z) / sqrt(1 - den)
    lambdaz <- cbind(df1, df2, lambdaz)
    df1 <- lambdaz[, 1]
    df2 <- sqrt(lambdaz[, 2])
    lambdaz <- lambdaz[, 3]
    logpdf1 <- -sum(log(diag(dec))) -
      0.5 * q * log(2 * pi) - 0.5 * colSums(z^2)
    logpdf2 <- -sum(log(diag(dec))) + q * logb(df2) -
      0.5 * q * log(2 * pi) - 0.5 * colSums(z^2) * (df2^2)
    infinitex <- rowSums(is.infinite(x)) > 0
    ok <- (infinitex + is.na(lambdaz) + is.infinite(lambdaz)) == 0
    logpdf1[!ok] <- NaN
    logpdf2[!ok] <- NaN
    if (any(ok)) {
      logpdf1[ok] <- pnorm(lambdaz[ok], log.p = TRUE) + logpdf1[ok] +
        logb(2) + logb(1-df1[ok])
      logpdf2[ok] <- pnorm(lambdaz[ok] * df2[ok], log.p = TRUE) +
        logpdf2[ok] + logb(2) + logb(df1[ok])
    }
    logpdf <- logb (exp(logpdf1) + exp(logpdf2))
    logpdf <- replace(logpdf, infinitex, -Inf)
    names(logpdf) <- rownames(x)
    return(if (log)
      logpdf
      else exp(logpdf))
  }
  else {
    if (all(delta == 0)) {
      return(dmvnorm(x, mean = mu, Var = Omegabar, log = log))
    }
    Omega <- Omegabar + outer(delta, delta)
    dec <- tryCatch(chol(Omega), error = function(e) e)
    if (inherits(dec, "error")) {
      mean <- mu + sqrt(2/pi) * delta
      x.is.mu <- colSums(t(x) != mean) == 0
      logpdf <- rep.int(-Inf, nrow(x))
      logpdf[x.is.mu] <- Inf
      names(logpdf) <- rownames(x)
      return(if (log)
        logpdf
        else
          exp(logpdf))
    }
    z <- backsolve(dec, t(x) - mu, transpose = TRUE)
    lambdaz <- c(backsolve(dec, delta, transpose = TRUE))
    den <- 1 / sqrt(1 - sum(lambdaz^2))
    lambdaz <- c(t(lambdaz) %*% z) * den
    logpdf <- - sum(log(diag(dec))) - 0.5 * (q * log(2 * pi) +
                                               colSums(z^2))
    infinitex <- rowSums(is.infinite(x)) > 0
    ok <- (infinitex + is.na(lambdaz) + is.infinite(lambdaz)) == 0
    logpdf[!ok] <- NaN
    if (any(ok)) {
      logpdf[ok] <- pnorm(lambdaz[ok], log.p = TRUE) + logpdf[ok] + logb(2)
    }
    logpdf <- replace(logpdf, infinitex, -Inf)
    names(logpdf) <- rownames(x)
    return(if (log)
      logpdf
      else exp(logpdf))
  }
}
