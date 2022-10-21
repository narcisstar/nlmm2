#'
#' @name nlssnmdata
#' @aliases fake.ssnmdata
#' @title Generate data from a non linear mixed effects regression (nlmer) model under skew scale mixture of normal (ssmn) distributions.
#'
#' @description The function generate data for testing functions related to the fit of \code{nlmer} models.
#'
#' @param n number of clusters/subjects.
#' @param n_i number of observations/measurements per subject.
#' @param n_x number of covariates entering the non linear function along with the mixed effects (non linear parameters).
#' @param n_z number of covariates with fixed slopes.
#' @param n_w number of covariates with random slopes.
#' @param nlfun non linear function. Its first two arguments must be (named) \code{phi} (mixed effects) and \code{X} (vector of covariates)
#' @param Xdist character vector of length \code{n_x}. The \code{i}-th element specifies the distribution to used to simulate the \code{i}-th element of the vector of covariates entering the non linear function. Each character "\code{name}" must be such that "\code{rname(n)}" generate \code{n} elements from the specified distribution. For exemple, "\code{norm}" corresponds to the standard normal distribution, "\code{mydist}" where \code{mydist = function(n) rpois (n, lambda = 5)} corresponds to a Poisson distribution with mean \code{5}.
#' @param A_i numeric \code{r} by \code{p} matrix, common fixed design matrix for all clusters (\code{p} = number of fixed effects).
#' @param Zdist character vector of length \code{n_z}. The \code{i}-th element specifies the distribution to used to simulate the \code{i}-th fixed covariate.
#' @param Zposition logical matrix of the same size as \code{A_i}. The elements of \code{Zposition} indicate which elements of \code{A_i} are fixed covariates (so the corresponding random effect is a slope). Accordingly, the equality \code{sum(Zposition) == n_z} must hold.
#' @param B_i numeric \code{r} by \code{q} matrix, common random design matrix for all clusters (\code{q} = number of random effects).
#' @param Wdist character vector of length \code{n_w}. The \code{i}-th element specifies the distribution to used to simulate the \code{i}-th random covariate.
#' @param Wposition logical matrix of the same size as \code{B_i}. The elements of \code{Wposition} indicate which elements of \code{B_i} are random covariates (so the corresponding random effect is a slope). Accordingly, the equality \code{sum(Wposition) == n_w} must hold.
#' @param fixef numeric \code{p} vector of fixed effects.
#' @param omegabar2,deltao numeric scalars (\code{omegabar2 > 0}), residual scale and shape parameters.
#' @param Omegabar_b,delta_b numeric \code{q} by \code{q} positive definite matrix and \code{q} vector, scale and shape parameters of random effects.
#' @param smixvar,smixvar_b characters, specify the mixing distributions for residuals and random effetcs (see \code{smixvar } at \link{rmssmn}).
#' @param nu,nu_b degrees of freedom parameters for residuals and random effetcs (see \code{df} at \link{rmssmn}).
#'
nlssnmdata <-
  function (n = 10, n_i = 3, n_x = 2, n_z, n_w = 0,
            nlfun = function(phi, X) {
              exp(phi[1] + phi[2] * X[1] + phi[3] * log(X[2]))
            },
            Xdist = rep("norm", n_x), # Fixed covariates in X_ij
            A_i = diag(3), Zdist = rep("unif", n_z), Zposition = NULL, # Fixed covariates in A_i
            B_i = diag(3), Wdist = rep("unif", n_w), Wposition = NULL, # Random covariates in B_i
            fixef = c(2, -.5, 1), omegabar2 = 1, deltao = 2,
            Omegabar_b = diag(q), delta_b = rep(0, q),
            smixvar = "Gamma", smixvar_b = smixvar,
            nu = switch(smixvar, Bin = c(.5, .5), Gamma = 5, Beta = 5, NULL),
            nu_b = switch(smixvar_b, Bin = c(.5, .5), Gamma = 5, Beta = 5, NULL)) {
    n <- max(2, n)
    n_i <- max(2, n_i)
    N <- n * n_i
    q <- NCOL(B_i)
    if (!missing(Xdist))
      n_x <- length(Xdist)
    X <- sapply(1:n_x, FUN = function(k) {
      eval(str2expression(paste0("r", Xdist[k], "(n = ", N, ")")))
    })
    if (missing(Xdist))
      X[,2] <- abs(X[,2])
    cluster <- rep(1:n, each = n_i)
    subject <- lapply(1:n, FUN = function(i) {
      ((i - 1) * n_i + 1):(i * n_i)
    })
    ckappa_e1 <- switch(smixvar,
                        Bin = sqrt(2/pi) * (1 + nu[1] * (1/sqrt (nu[2]) - 1)),
                        Gamma = sqrt(nu[1]/pi) * gamma(.5*(nu[1]-1)) / gamma(nu[1]/2),
                        Beta = sqrt(2/pi) * nu[1]/(nu[1]-.5),
                        sqrt(2/pi))
    ctau_e1 <- switch(smixvar_b,
                      Bin = sqrt(2/pi) * (1 + nu_b[1] * (1/sqrt (nu_b[2]) - 1)),
                      Gamma = sqrt(nu_b[1]/pi) * gamma(.5*(nu_b[1]-1)) / gamma(nu_b[1]/2),
                      Beta = sqrt(2/pi) * nu_b[1]/(nu_b[1]-.5),
                      sqrt(2/pi))
    b <- rmvnorm(n = n, Var = Omegabar_b)
    U0 <- switch(smixvar_b,
                 Gamma = rgamma(n = n, shape = nu_b/2, rate = nu_b/2),
                 Beta = rbeta(n = n, shape1 = nu_b, shape2 = 1),
                 rep(1, n))
    if (identical(smixvar_b, "Bin")) {
      binu <- rbinom(n, size = 1, prob = nu_b[1])
      if (any(binu))
        U0[binu] <- nu_b[2]
    }
    U <- 1/sqrt(U0)
    V <- abs(rnorm(n=n))
    add <- V * U - ctau_e1
    b <- outer(add, delta_b) + b * U
    if (n_w == 0 || is.null(Wposition)) {
      if (n_z == 0 || is.null(Zposition)) {
        phi <- t(sapply(1:n, FUN = function(i) {
          A_i %*% fixef + B_i %*% b[i,]
        }))
      }
      else {
        if (!missing(Zdist))
          n_z <- length(Zdist)
        Z <- sapply(1:n_z, FUN = function(k) {
          eval(str2expression(paste0("r", Zdist[k], "(n = ", n, ")")))
        })
        if (missing(Xdist))
          Z <- -1 + 2 * Z
        phi <- t(sapply(1:n, FUN = function(i) {
          A_ii <- A_i
          A_ii[Zposition] <- Z[i,]
          A_ii %*% fixef + B_i %*% b[i,]
        }))
      }
    }
    else {
      if (!missing(Wdist))
        n_w <- length(Wdist)
      W <- sapply(1:n_w, FUN = function(k) {
        eval(str2expression(paste0("r", Wdist[k], "(n = ", n, ")")))
      })
      if (missing(Xdist))
        W <- -1 + 2 * W
      if (n_z == 0 || is.null(Zposition)) {
        phi <- t(sapply(1:n, FUN = function(i) {
          B_ii <- B_i
          B_ii[Wposition] <- W[i,]
          A_i %*% fixef + B_ii %*% b[i,]
        }))
      } else {
        if (!missing(Zdist))
          n_z <- length(Zdist)
        Z <- sapply(1:n_z, FUN = function(k) {
          eval(str2expression(paste0("r", Zdist[k], "(n = ", n, ")")))
        })
        if (missing(Xdist))
          Z <- -1 + 2 * Z
        phi <- t(sapply(1:n, FUN = function(i) {
          A_ii <- A_i
          A_ii[Zposition] <- Z[i,]
          B_ii <- B_i
          B_ii[Wposition] <- W[i,]
          A_ii %*% fixef + B_ii %*% b[i,]
        }))
      }
    }
    if (NCOL(X) == 1)
      mu <- c(sapply(1:n, FUN = function(i) {
        X_i <- X[subject[[i]]]
        sapply(X_i, FUN = function(x) {
          nlfun (phi = phi[i,], X = x)
        })
      }))
    else
      mu <- c(sapply(1:n, FUN = function(i) {
        X_i <- X[subject[[i]],]
        c(apply(X_i, MARGIN = 1, FUN = function(x_i) {
          nlfun (phi = phi[i,], X = x_i)
        }))
      }))
    Y <- rssmn (n = N, omegabar = sqrt(omegabar2),
                delta = deltao, smixvar  = smixvar, df = nu) +
      mu - ckappa_e1 * deltao
    yx <- as.data.frame(cbind(Y, X, cluster))
    colnames(yx) <- c("y", paste0("x", 1:n_x), "cluster")
    if (n_w > 0) {
      yx <- cbind(yx, W[cluster,])
      colnames(yx) <- c(colnames(yx), paste0("w", 1:n_w))
    }
    if (n_z > 0) {
      yx <- cbind(yx, Z[cluster,])
      colnames(yx) <- c(colnames(yx), paste0("z", 1:n_z))
    }
    list(data = yx, A = A_i, B = B_i, b = b, U = U0, V = V, resid = Y - mu,
         param = list(fixef = fixef, omegabar2 = omegabar2,
                      deltao = deltao, Omegabar_b = Omegabar_b,
                      delta_b = delta_b, nu = nu, nu_b = nu_b,
                      ckappa_e1 = ckappa_e1, ctau_e1 = ctau_e1),
         smixvar = smixvar, smixvar_b = smixvar_b)
  }

fake.ssnmdata <-
  function (n = 10, n_i = 3, n_x = 2,
            nlfun = function(phi, X) {
              exp(phi[1] + phi[2] * X[1] + phi[3] * log(X[2]))
            }, Xdist = rep("norm", n_x), A_i = diag(3), B_i = diag(3),
            fixef = c(2, -0.5, 1), omegabar2 = 1, deltao = 2, Omegabar_b = diag(q),
            delta_b = rep(0, q), smixvar = "Gamma", smixvar_b = "Gamma",
            nu = switch(smixvar, Bin = c(0.5, 0.5), Gamma=5, Beta=5, NULL),
            nu_b = switch(smixvar_b,Bin = c(0.5, 0.5),  Gamma=5, Beta=5, NULL)) {
    n <- max(2, n)
    n_i <- max(2, n_i)
    N <- n * n_i
    q <- NCOL(B_i)
    if (!missing(Xdist))
      n_x <- length(Xdist)
    X <- sapply(1:n_x, FUN = function(k) {
      eval(str2expression(paste0("r", Xdist[k], "(n = ", N, ")")))
    })
    if (missing(Xdist))
      X[, 2] <- abs(X[, 2])
    cluster <- rep(1:n, each = n_i)
    subject <- lapply(1:n, FUN = function(i) {
      ((i - 1) * n_i + 1):(i * n_i)
    })


    ckappa_e1 <- switch(smixvar, Bin = sqrt(2/pi) * (1 + nu[1] * (1/sqrt(nu[2]) - 1)),
                        Gamma = sqrt(nu[1]/pi) * gamma(0.5*(nu[1] - 1))/gamma(nu[1]/2),
                        Beta = sqrt(2/pi) * nu[1]/(nu[1] - 0.5), sqrt(2/pi))

    ctau_e1 <- switch(smixvar_b, Bin = sqrt(2/pi) * (1 + nu_b[1]*(1/sqrt(nu_b[2]) - 1)),
                      Gamma = exp(0.5 * (log(nu_b) - log(pi)) +
                                    lgamma(0.5 * (nu_b - 1)) - lgamma(0.5 * nu_b)),
                      Beta = sqrt(2/pi) * nu_b[1]/(nu_b[1] - 0.5), sqrt(2/pi))

    b <- rmvnorm(n = n, Var = Omegabar_b)
    U0 <- switch(smixvar_b, Gamma = rgamma(n = n, shape = nu_b/2, rate = nu_b/2),
                 Beta = rbeta(n = n, shape1 = nu_b, shape2 = 1), rep(1, n))
    if (identical(smixvar_b, "Bin")) {
      binu <- rbinom(n, size = n, prob = nu_b[1])
      if (any(binu))
        U0[binu] <- nu_b[2]
    }
    U <- 1/sqrt(U0)
    V <- abs(rnorm(n = n))
    add <- V * U - ctau_e1
    b <- outer(add, delta_b) + b * U
    phi <- t(sapply(1:n, FUN = function(i) {
      A_i %*% fixef + B_i %*% b[i, ]
    }))
    if (NCOL(X) == 1)
      mu <- c(sapply(1:n, FUN = function(i) {
        X_i <- X[subject[[i]]]
        sapply(X_i, FUN = function(x) {
          nlfun(phi = phi[i, ], X = x)
        })
      }))
    else mu <- c(sapply(1:n, FUN = function(i) {
      X_i <- X[subject[[i]], ]
      c(apply(X_i, MARGIN = 1, FUN = function(x_i) {
        nlfun(phi = phi[i, ], X = x_i)
      }))
    }))
    Y <- rssmn(n = N, omegabar = sqrt(omegabar2), delta = deltao,
               smixvar = smixvar, df = nu) + mu - ckappa_e1 * deltao
    yx <- as.data.frame(cbind(Y, X, cluster))
    colnames(yx) <- c("y", paste0("x", 1:n_x), "cluster")
    list(data = yx, A = A_i, B = B_i, b = b, U = U0, V = V, resid = Y -
           mu, param = list(fixef = fixef, omegabar2 = omegabar2,
                            deltao = deltao, Omegabar_b = Omegabar_b, delta_b = delta_b,
                            nu = nu, nu_b = nu_b, ckappa_e1 = ckappa_e1, ctau_e1 = ctau_e1),
         smixvar = smixvar, smixvar_b = smixvar_b)
  }
