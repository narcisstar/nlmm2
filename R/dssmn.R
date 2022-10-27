# Density of a univariate skew scale mixture of normal distribution
#' @export
dssmn <- function (x, mu = 0, omegabar = 1, delta = 0, mixvar = "Gamma",
                   df = switch(mixvar, Bin = c(.5, .5), 5), log = FALSE) {
  if (identical(mixvar, "Bin") && length(df) < 2)
    stop("the (Skew) contaminated model requires a length 2 'df' argument")
  if (identical(mixvar, "Beta")) {
    lambda <- delta / omegabar
    omega <- sqrt((omegabar^2) + (delta^2))
    z <- (x - mu) / omega
    zlambda <- cbind(z, df, lambda)
    z <- zlambda[, 1]
    df <- zlambda[, 2]
    lambda <- zlambda[, 3]
    n <- length(z)
    logpdf <- numeric(n)
    ok_0 <- (is.na(lambda) + is.infinite(lambda) + (omegabar <= 0) +
              (df <= 0) + is.infinite(x)) == 0
    ok <- (ok_0 + is.finite(df)) == 2
    logpdf[!ok] <- NaN
    if (any(ok)) {
      fdsn <- function(u) {
        sqrtu <- sqrt(u)
        (u^(df[ok]-.5)) * pnorm(lambda[ok] * z[ok] * sqrtu) * dnorm(z[ok] * sqrtu)
      }
      logpdf[ok] <- logb(df[ok]) + logb(
        cubature::adaptIntegrate(f = fdsn, lowerLimit = 0,
                                 upperLimit = 1,
                                 fDim = sum(ok))$integral)
    }
    ok <- (ok_0 + is.infinite(df)) == 2
    if (any(ok)) {
      logpdf[ok] <- pnorm(lambda[ok] * z[ok], log.p = TRUE) + dnorm(z[ok], log = TRUE)
    }
    logpdf <- replace(logpdf, is.infinite(x), -Inf)
    return(
      if (log)
        logb(2) + logpdf - logb(omega)
      else 2 * exp(logpdf) / omega
    )
  }
  else if (identical(mixvar, "Gamma")) {
    if (all(delta == 0)) {
      logpdf <- dt((x - mu)/omegabar, df = df, log = TRUE) - logb(omegabar)
      return(if (log)
        logpdf
        else
          exp(logpdf))
    }
    else {
      lambda <- delta / omegabar
      omega <- sqrt((omegabar^2) + (delta^2))
      z <- (x - mu) / omega
      logpdf <- dt(z, df = df, log = TRUE) - logb(omega)
      zlambda <- cbind(z, df, lambda)
      z <- zlambda[, 1]
      df <- zlambda[, 2]
      lambda <- zlambda[, 3]
      ok_0 <- (is.na(lambda) + is.infinite(lambda) + (omegabar <= 0) +
                 (df <= 0) + is.infinite(x)) == 0
      ok <- (ok_0 + is.finite(df)) == 2
      logpdf[!ok] <- NaN
      if (any(ok)) {
        ratio <- sqrt((df + 1) / (df + z^2))
        logpdf[ok] <- pt(lambda[ok] * z[ok] * ratio[ok],
                         df = df + 1, log.p = TRUE) +
          logpdf[ok] + logb(2)
      }
      ok <- (ok_0 + is.infinite(df)) == 2
      if (any(ok)) {
        logpdf[ok] <- pnorm(lambda[ok] * z[ok], log.p = TRUE) +
          dnorm(z[ok], log = TRUE) - logb(omega) + logb(2)
      }
      logpdf <- replace(logpdf, is.infinite(x), -Inf)
      return(
        if (log)
          logpdf
        else exp(logpdf)
      )
    }
  }
  else if (identical(mixvar, "Bin")) {
    lambda <- delta / omegabar
    omega <- sqrt((omegabar^2) + (delta^2))
    z <- (x - mu) / omega
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
    zlambda <- cbind(z, df1, df2, lambda)
    z <- zlambda[, 1]
    df1 <- zlambda[, 2]
    df2 <- sqrt(zlambda[, 3])
    lambda <- zlambda[, 4]
    logpdf1 <- (-log(sqrt(2 * pi)) - logb(omega) - z^2/2)
    logpdf2 <- (-log(sqrt(2 * pi)) - logb(omega) + logb(df2) - (df2 * z)^2/2)
    ok <- (is.na(lambda) + is.infinite(lambda) + (omegabar <= 0) +
             (df1 <= 0) + (df1 >= 1) + (df2 <= 0) + (df2 > 1) +
             is.infinite(x)) == 0
    logpdf1[!ok] <- NaN
    logpdf2[!ok] <- NaN
    if (any(ok)) {
      logpdf1[ok] <- pnorm(lambda[ok] * z[ok], log.p = TRUE) +
        logpdf1[ok] + logb(2) + logb(1 - df1[ok])
      logpdf2[ok] <- pnorm(lambda[ok] * z[ok] * df2[ok], log.p = TRUE) +
        logpdf2[ok] + logb(2) + logb(df1[ok])
    }
    logpdf <- logb (exp(logpdf1) + exp(logpdf2))
    logpdf <- replace(logpdf, is.infinite(x), -Inf)
    return(if (log)
      logpdf
      else exp(logpdf))
  }
  else {
    if (all(delta == 0)) {
      return(dnorm(x, mean = mu, sd = omegabar, log = log))
    }
    else {
      lambda <- delta / omegabar
      omega <- sqrt((omegabar^2) + (delta^2))
      z <- (x - mu) / omega
      logpdf <- (-log(sqrt(2 * pi)) - logb(omega) - z^2/2)
      zlambda <- cbind(z, lambda)
      z <- zlambda[, 1]
      lambda <- zlambda[, 2]
      ok <- (is.infinite(x) + is.na(lambda) +
               is.infinite(lambda) + (omegabar <= 0)) == 0
      logpdf[!ok] <- NaN
      if (any(ok)) {
        logpdf[ok] <- pnorm(lambda[ok] * z[ok], log.p = TRUE) +
          logpdf[ok] + logb(2)
      }
      logpdf <- replace(logpdf, is.infinite(x), -Inf)
      return(if (log)
        logpdf
        else exp(logpdf))
    }
  }
}
