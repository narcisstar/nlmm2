#' @name Normal
#' @aliases dmvnorm dmnorm rmvnorm
#' @title Multivariate Normal Distribution
#' @description Multivariate Normal Density function and Random deviate generator
#'
#' @usage dmvnorm (x, mean = rep(0, p), Var = diag(p), log = FALSE)
#'
#' dmnorm (x, mean = rep(0, p), Var = diag(p), log = FALSE)
#'
#' rmvnorm (n, mean = rep(0, nrow(Var)), Var = diag(length(mean)))
#'
#' @param x vector or matrix of quantiles. If x is a matrix, each row is taken to be a quantile.
#' @param n number of deviates to be generated.
#' @param mean mean vector.
#' @param Var covariance matrix.
#' @param log logical; if \code{TRUE}, densities \code{d} are given as \code{log(d)}.
#'
#' @details See \code{heavy::rmnorm} for a faster alternative to \code{rmvnorm}, and \code{dmnorm} is just an alias for \code{dmvnorm}.
#'

# Multivariate normal density
dmvnorm <- function (x, mean = rep(0, p), Var = diag(p), log = FALSE) {
  if (is.vector(x))
    x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  if (!missing(mean)) {
    if (!is.null(dim(mean)))
      dim(mean) <- NULL
    if (length(mean) != p)
      stop("mean and Var have non-conforming size")
  }
  if (!missing(Var)) {
    if (p != ncol(Var))
      stop("x and Var have non-conforming size")
    if (!isSymmetric(Var, tol = sqrt(.Machine$double.eps),
                     check.attributes = FALSE))
      stop("Var must be a symmetric matrix")
  }
  dec <- tryCatch(chol(Var), error = function(e) e)
  if (inherits(dec, "error")) {
    x.is.mu <- colSums(t(x) != mean) == 0
    logretval <- rep.int(-Inf, nrow(x))
    logretval[x.is.mu] <- Inf
  } else {
    tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
    rss <- colSums(tmp^2)
    logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 *
                                                        pi) - 0.5 * rss
  }
  names(logretval) <- rownames(x)
  if (log)
    logretval
  else exp(logretval)
}
dmnorm <- dmvnorm

# Multivariate normal variates
rmvnorm <- function (n, mean = rep(0, nrow(Var)), Var = diag(length(mean))) {
  if (!isSymmetric(Var, tol = sqrt(.Machine$double.eps),
                   check.attributes = FALSE)) {
    stop("Var must be a symmetric matrix")
  }
  if (length(mean) != nrow(Var))
    stop("mean and Var have non-conforming size")
  ev <- eigen(Var, symmetric = TRUE)
  if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
    warning("Var is numerically not positive semidefinite")
  }
  R <- t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
  retval <- matrix(rnorm(n * ncol(Var)), nrow = n, byrow = FALSE) %*% R
  retval <- sweep(retval, 2, mean, "+")
  colnames(retval) <- names(mean)
  retval
}
