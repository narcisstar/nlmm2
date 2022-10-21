.erowMeans <- function (x, na.rm = FALSE, dims = 1)
{
  if (is.null(dim(x)))
    mean(x, na.rm = na.rm)
  else rowMeans(x = x, na.rm = na.rm, dims = dims)
}

#' Make lambda
#'
#' @export
Mlambda <- function(omegabar, delta){
  delta <- c(delta)
  lan <- (solve(sqrtm(omegabar)) %*% delta)/c(sqrt(1 - t(delta) %*% solve(omega) %*% delta))
  return(lan)
}

#' Make dispersion
#' @export

MDispersion <- function(omegabar, delta){
  omegabar + delta %*% t(delta)
}

#' Matrix and numeric square root
#' @description Computes Matrix and numeric number square root using eigen analysis.
#' @param x A matrix or numeric value.
#' @return Square root of x.
#' @examples sqrtm(4) # returns numeric
#' sqrtm(matrix(64)) # returns matrix
#' A <- matrix(data = c(1, 2, 3, 2, 20, 26, 3, 26, 70),
#' ncol = 3); A
#' a <- sqrtm(A); a
#' a %*% a ## ==A; matrix square
#' a * a ## not equals A; simple square
#' @seealso \code{\link{sqrt}}
#' @export
sqrtm <- function(x)
{
  if(is.matrix(x) && (dim(x)[1] > 1))
  {
    x.eig <- eigen(x)
    val <- x.eig$values
    vec <- x.eig$vectors
    s <- vec %*% diag(sqrt(val)) %*% solve(vec)
  }
  else s <- sqrt(x)
  return(s)
}

#' Trace of matrix
#'
#' @export
tr <- function(x)
{
  if(!is.matrix(x)) stop("'x' must be a matrix")
  sum(diag(x))
}

#' Make scale matrix or scale scalar
#' @param omega dispersion parameter, stand for sigma square
#' @param delta numeric, shape/skewness vector coming from \code{\link{Mdelta}}
#' @param lambda numeric, shape/skewness parameter
#' @export
MOmegabar <- function(omega, delta, lambda){
  if(missing(delta) && missing(lambda))
    stop("Specify either 'delta' or 'lambda'.")
  if(missing(delta))
    delta <- Mdelta(omega = omega, lambda = lambda)
  if (NCOL(omega) > 1) {
    if (length(delta) < 2)
      cat("Inconfortable 'omega' and 'delta'.")
    else
      omegabar <- omega - delta %*% t(delta)
  }
  else
    omegabar <- c(omega) - c(delta)^2
  return(omegabar)
}

#' Make delta
#' @param omega dispersion parameter, stand for sigma square
#' @param lambda numeric, shape/skewness parameter
#' @export

Mdelta <- function(omega, lambda){
  if (NCOL(omega) > 1) {
    if (length(lambda) < 2)
      cat("Inconfortable 'omega' and 'lambda'.")
    else
      del <- c(solve(sqrtm(1 + t(lambda) %*% lambda))) * sqrtm(omega) %*% lambda
  }
  else
    del <- c(omega) * c(lambda)/ sqrt(1 + c(lambda)^2)
  return(c(del))
}

#lambda = 1:3
#omega = diag(3)

#' Frame for nlmm
#'
#' @export
nlmer.frame <- function (formula, X.formula, fixef.name,
                         cluster, family = "Gamma", data, subset,
                         Zi = NULL, fixedStart,
                         na.action = options("na.action")$na.action
                         )
{
  mcall <- match.call()
  if(!missing(fixedStart))
    phi = fixedStart
  if(missing(fixef.name) && missing(fixedStart))
    stop("Specify either 'fixef.name' or 'fixedStart'.")

  if (missing(data)) {
    data <- merge.environment(environment(formula), environment(X.formula))
  }
  Tformula <- formula
  Tformula[[3]] <- X.formula[[2]]
  if (missing(subset)) {
    X <- model.frame(Tformula, data = data, na.action = na.action)
  }
  else {
    X <- model.frame(Tformula, data = data, subset = subset,
                     na.action = na.action)
  }
  Y <- X[, 1]
  N <- length(Y)
  yname <- rownames(X)
  X <- X[, -1]
  if (missing(subset)) {
    cluster <- model.frame(cluster, data = data, na.action = na.action)[,
                                                                        1]
  }
  else {
    cluster <- model.frame(cluster, data = data, subset = subset,
                           na.action = na.action)[, 1]
  }
  subject <- unique(cluster)
  if (is.numeric(cluster)) {
    subject.nm <- names(cluster)
    if (is.null(subject.nm)) {
      subject.nm <- paste0("cl.", subject)
    }
    else {
      subject.nm <- unique(subject.nm)
    }
  }
  else {
    subject0 <- subject
    cl0 <- cluster
    subject <- 1:length(subject)
    subject.nm <- names(subject) <- subject0
    cluster <- sapply(cluster, FUN = function(cl) {
      subject[subject0 == cl]
    })
    names(cluster) <- cl0
  }
  n <- length(subject)
  r <- 1:N
  subject <- lapply(subject, FUN = function(cl) {
    r[cluster == cl]
  })
  reorder <- NULL
  for (i in 1:n) {
    reorder <- c(reorder, subject[[i]])
  }
  if(!missing(fixef.name))
    p <- length(fixef.name)
  else
    p <- length(fixedStart)

  if(is.null(Zi))
    q <- p

  #####
  if (NCOL(X) == 1) {
    datalist <- lapply(1:n, FUN = function(i) {
      if(is.null(Zi)){
        list(Y_i = Y[subject[[i]]], X_i = X[subject[[i]]],
             Zi = matrix(1, nrow = length(Y[subject[[i]]]), ncol = q))
      }
      else{
        list(Y_i = Y[subject[[i]]], X_i = X[subject[[i]]],
             Zi = as.matrix(eval(Zi, envir = list(x = X[subject[[i]]]))))
        }
    })

    nlfun <- function(phi, X) {
      c(sapply(X, FUN = function(x) {
        eval(formula[[3]], envir = list(phi = phi, X = x))
      }))
    }
  }
  else {
    if(is.null(Zi)){
      datalist <- lapply(1:n, FUN = function(i) {
        list(Y_i = Y[subject[[i]]], X_i = X[subject[[i]],
        ], Zi = matrix(1, nrow = length(Y[subject[[i]]]), ncol = q))
      })
    }
    else
      datalist <- lapply(1:n, FUN = function(i) {
      list(Y_i = Y[subject[[i]]], X_i = X[subject[[i]],
      ], Zi = as.matrix(eval(Zi, envir = list(phi = phi, x = X[subject[[i]], ]))))
    })
    nlfun <- function(phi, X) {
      c(apply(X, MARGIN = 1, FUN = function(x) {
        eval(formula[[3]], envir = list(phi = phi, X = x))
      }))
    }
  }
  q <- NCOL(datalist[[1]]$Zi)
  names(datalist) <- subject.nm

  if(missing(fixef.name))
    fixef.name <- paste0("alpha.", 1:p)

  if (is.null(yname))
    yname <- paste0("y", 1:N)
  frame <- list(dim = list(N = N, n = n, p = p, q = q),
                nlfun = nlfun, fixef.name = fixef.name, data = datalist)

  frame$family <- family
  frame$call <- mcall
  structure(frame, class = "nlmm.frame")
}



#' Pereira EM algorithm (fit)
#'
#' @export

pereEM <- function(formula, X.formula, fixef.name = fixef.name,
                   Zi = NULL, cluster, family = "Gamma", data,
                   subset, na.action = options("na.action")$na.action,
                   start = list(alpha = c(100, 700, 349),
                                D = diag(c(.5, 1, 1.2)), lambda = rep(4, 3),
                                sigma2 = 3, nu = NULL),
                   pverbose = TRUE, maxit = 40, tol = 0.001, maxit.opt = 16,
                   skew = TRUE){# start = pstart,
  begin <- Sys.time()

  if(is.null(start$nu))
    start$nu <- switch(family, Gamma = 1, Bin = rep(.5, 2), Beta = 1)

  param <- start[c("alpha", "D", "lambda", "sigma2", "nu")]

  frame = nlmer.frame(formula = formula, X.formula = X.formula,
                      fixef.name = fixef.name, cluster = cluster,
                      family = family, data = data, subset = subset, Zi = Zi,
                      fixedStart = param$alpha, na.action = na.action)

  param$Delta <- Mdelta(omega = param$D, lambda = param$lambda)
  param$Gamma <- MOmegabar(omega = param$D, delta = param$Delta)
  #print(param$Delta)

  k= 0
  dif <- dif.logLike <- Inf
  old = old.logLike = 0
  while ((any(max(dif) >= tol) >= tol) && (maxit > k)) {#(abs(dif.logLike)
    k = k + 1
    if(pverbose) cat("iter ", k, "\n")
    mdata <- stepE(param = param, frame = frame)
    param$alpha <- optim(par = param$alpha, fn = ALPHA, mdata = mdata, frame = frame, method = "BFGS", control = list(maxit = maxit.opt))$par #
    names(param$alpha) <- fixef.name

    param$Delta <- colSums(t(sapply(1:frame$dim$n, function(i){
      mdata[[i]]$utbi
    })))/sum(sapply(1:frame$dim$n, function(i){
      mdata[[i]]$ut2i
    }))
    #print(param)

    param$Gamma <- matrix(.erowMeans(sapply(1:frame$dim$n, function(i){
      #mdata[[i]]$ubbti - mdata[[i]]$utbi %*% t(param$Delta) - param$Delta %*% t(mdata[[i]]$utbi) + mdata[[i]]$ut2i * param$Delta %*% t(param$Delta)

      mdata[[i]]$ubbti - mdata[[i]]$utbi %*% t(mdata[[i]]$deltai) - mdata[[i]]$deltai %*% t(mdata[[i]]$utbi) + mdata[[i]]$ut2i * mdata[[i]]$deltai %*% t(mdata[[i]]$deltai)
    }), dims = 1), ncol = frame$dim$q)

    param$sigma2 <- (1/frame$dim$N) * sum(sapply(1:frame$dim$n, function(i){
      c(mdata[[i]]$ui * t(frame$data[[i]]$Y_i - frame$nlfun(param$alpha, frame$data[[i]]$X_i)) %*% (frame$data[[i]]$Y_i - frame$nlfun(param$alpha, frame$data[[i]]$X_i)) - t(frame$data[[i]]$Y_i - frame$nlfun(param$alpha, frame$data[[i]]$X_i)) %*% frame$data[[i]]$Zi %*% mdata[[i]]$ubi - t(mdata[[i]]$ubi) %*% t(frame$data[[i]]$Zi) %*% (frame$data[[i]]$Y_i - frame$nlfun(param$alpha, frame$data[[i]]$X_i)) + tr(frame$data[[i]]$Zi %*% mdata[[i]]$ubbti %*% t(frame$data[[i]]$Zi)))
    }))
    if (any(frame$family == c("Bin", "Gamma", "Beta")))
      param$nu <- optim(par = param$nu, fn = gyi, param = param, frame = frame, mdata = mdata, method = "BFGS", control = list(maxit = maxit.opt))$par

    #print(optim(par = param$nu, fn = gyi, param = param, frame = frame, mdata = mdata, method = "BFGS", control = list(maxit = maxit.opt))$par)

    param$D <- param$Gamma + param$Delta %*% t(param$Delta)

    #print(1 - t(c(param$Delta)) %*% solve(param$D) %*% c(param$Delta))
    #print(solve(sqrtm(param$D)) %*% param$Delta)

    if(skew)
      param$lambda <- (solve(sqrtm(param$D)) %*% param$Delta)/c(sqrt(1 - t(c(param$Delta)) %*% solve(param$D) %*% c(param$Delta)))
    else
      param$lambda <- rep(0, frame$dim$p)

    logLike <- sum(sapply(X = seq_len(frame$dim$n), FUN = function(x) mdata[[x]]$Q1i)) + sum(sapply(X = seq_len(frame$dim$n), FUN = function(x) mdata[[x]]$Q2i))
    if(pverbose)
      cat("Loglik =", logLike, "\n")

    if (!skew){
      dif <- c(param$alpha, param$sigma2, param$D) - old
      old <- c(param$alpha, param$sigma2, param$D)
    }
    else {
      dif <- c(param$alpha, param$sigma2, param$lambda, param$D) - old
      old <- c(param$alpha, param$sigma2, param$lambda, param$D)
    }
    #print(param)


    dif.logLike <- logLike - old.logLike
    old.logLike <- logLike
    #print((dif.logLike >= tol) && (maxit > k) && (k < 2))
  }
  param$Delta <- param$Gamma <- NULL
  endtime <- Sys.time()
  totaltime <- (difftime(endtime, begin, units = "mins"))

  return(list(param = param, logLike = logLike, totaltime= totaltime, iter = k))
}

stepE <- function(param, frame)
{
  lapply(1:frame$dim$n, function(i){
    nu <- param$nu
    Zi <-  frame$data[[i]]$Zi
    mi <- length(frame$data[[i]]$Y_i)
    ind <- i

    #Page 9
    Ri = param$sigma2 * diag(nrow = mi)
    #print(Zi)
    #print(param$Gamma)

    Omegai = Ri + Zi %*% param$Gamma %*% t(Zi)
    Mi = c(solve(sqrtm(1 + t(param$Delta) %*% t(Zi) %*% solve(Omegai) %*% Zi %*% param$Delta)))
    Bi = solve(solve(param$Gamma) + t(Zi) %*% solve(Ri) %*% Zi)
    si = (diag(frame$dim$q) - Bi %*% t(Zi) %*% solve(Ri) %*% Zi) %*% param$Delta
    ri <- frame$data[[i]]$Y_i - frame$nlfun(param$alpha, frame$data[[i]]$X_i)
    mui <- c(Mi^2 * t(param$Delta) %*% t(Zi) %*% solve(Omegai) %*% ri) # devrait etre le vecteur de mean de y_i non
    Gi <- Bi %*% t(Zi) %*% solve(Ri) %*% ri

    # Page 7 eq 11
    ksi <- solve(sqrtm(param$D)) %*% param$lambda# ksi= \xi
    Sigmai <- Zi %*% param$D %*% t(Zi) + Ri
    Lambdai <- solve(solve(param$D) + t(Zi) %*% solve(Ri) %*% Zi)
    lambda_yi <- c((solve(sqrtm(Sigmai)) %*% Zi %*% param$D %*% ksi)/c(sqrtm(1 + t(ksi) %*% Lambdai %*% ksi)))

    Ai <- c(t(lambda_yi) %*% solve(sqrtm(Sigmai)) %*% ri) # Page 9

    #environment(tauir) <- environment()
    #environment(uir) <- environment()
    ui <- c(uir(frame, param, Sigmai, lambda_yi, ri, Ai, r = 1, i=ind))
    taui <- c(tauir(frame, param, Sigmai, lambda_yi, ri, Ai, r = 1, i=ind)); taui

    #Page 9
    uti <- ui * mui + Mi * taui; uti
    ut2i <- ui * mui^2 + Mi^2 + Mi * mui * taui; ut2i
    ubi <- ui * Gi + si %*% uti; ubi
    utbi <- c(Gi * uti + si * ut2i); utbi
    ubbti <- Bi + ui * Gi %*% t(Gi) + Gi %*% t(si) * uti + si %*% t(Gi) * uti + si %*% t(si) * ut2i; ubbti

    # additional qutities page 9
    Q1i <- -0.5*log(det(Ri)) - (ui /2)*t(frame$data[[i]]$Y_i - frame$nlfun(param$alpha, frame$data[[i]]$X_i))%*%solve(Ri)%*%(frame$data[[i]]$Y_i- frame$nlfun(param$alpha, frame$data[[i]]$X_i)) + t(frame$data[[i]]$Y_i- frame$nlfun(param$alpha, frame$data[[i]]$X_i))%*%solve(Ri)%*%Zi%*%ubi - 0.5*tr(solve(Ri)%*%Zi%*%ubbti %*%t(Zi)); Q1i

    Q2i <- -0.5*log(det(param$Gamma)) - 0.5*tr(solve(param$Gamma)%*%ubbti) + utbi%*%solve(param$Gamma)%*%param$Delta - (ut2i /2)%*%t(param$Delta)%*%solve(param$Gamma)%*%param$Delta; Q2i

    # Page 11
    mu_bi <- param$D %*% t(Zi) %*% solve(Sigmai) %*% (frame$data[[i]]$Y_i- frame$nlfun(param$alpha, frame$data[[i]]$X_i))

    bi <- mu_bi + c(tauir(frame, param, Sigmai, lambda_yi, ri, Ai, r = -1, i=ind)/sqrt(1 + t(ksi) %*% Lambdai %*% ksi)) * Lambdai %*% ksi
    yi_hat <- frame$nlfun(param$alpha, frame$data[[i]]$X_i) + Zi %*% bi

    #deltai <- Mdelta(omega = Sigmai, lambda = lambda_yi)
    deltai <- utbi/ut2i

    if (!any(frame$family == "Bin", frame$family == "Gamma", frame$family == "Beta"))
      return(list(ui= ui, ut2i= ut2i, ubi= ubi, utbi= utbi, ubbti= ubbti, Q1i= Q1i, Q2i= Q2i, ri= ri, Zi= Zi, yi_hat = yi_hat, Sigmai = Sigmai, lambda_yi = lambda_yi, deltai = deltai))
    else
      return(list(ui= ui, ut2i= ut2i, ubi= ubi, utbi= utbi, ubbti= ubbti, Q1i= Q1i, Q2i= Q2i, ri= ri, Zi= Zi, yi_hat = yi_hat, nu = nu, Sigmai = Sigmai, lambda_yi = lambda_yi, deltai = deltai))
  })

}


#' uir
#' @description conditional moments. On pages 4 and 5
#'
#' @export
uir <- function(frame, param, Sigmai, lambda_yi, ri, Ai, i = ind, r = 1)
{
  yi <- frame$data[[i]]$Y_i
  mi <- length(yi)

  if(frame$family == "Gamma"){
    #print(Mdelta(omega = Sigmai, lambda = lambda_yi))
    #print(MOmegabar(omega = Sigmai, lambda = lambda_yi))
    #print(param$nu)
    #print(frame$nlfun(param$alpha, frame$data[[i]]$X_i))
    #print(yi)
    #ui <- (EMMIXskew::ddmst(dat = yi, n = 1, p = mi, mean = rep(mui, mi), cov = Sigmai, nu = param$nu, del = rep(0, mi))/ddmst(dat = yi, n =1, p = mi, mean = rep(mui, mi), cov = Sigmai, nu = param$nu, del = lambda_yi)) * ((2^2*gamma((param$nu + mi + 2*1)/2) * (param$nu + t(ri)%*%solve(Sigmai)%*%ri)^-1 )/ gamma((param$nu+mi)/2)) * pt(sqrt((param$nu+mi+ 2)/(param$nu + t(ri)%*%solve(Sigmai)%*%ri))*Ai, param$nu+mi+2) #rep(mui, mi)
    ui <- (dmssmn(x = yi, mu = frame$nlfun(param$alpha, frame$data[[i]]$X_i), Omegabar = MOmegabar(omega = Sigmai, lambda = lambda_yi), df = param$nu, delta = rep(0, mi), mixvar = "Gamma")/dmssmn(x = yi, mu = frame$nlfun(param$alpha, frame$data[[i]]$X_i), Omegabar = MOmegabar(omega = Sigmai, lambda = lambda_yi), df = param$nu, delta = Mdelta(omega = Sigmai, lambda = lambda_yi), mixvar = "Gamma")
    ) *
      ((2^(r + 1) * gamma((param$nu + frame$dim$p + 2* r)/2) * (param$nu + t(ri) %*% solve(Sigmai) %*% ri)^-r )/ gamma((param$nu + frame$dim$p)/2)) * pt(sqrt((param$nu + frame$dim$p + 2 * r)/(param$nu + t(ri) %*% solve(Sigmai) %*% ri)) * Ai, param$nu + frame$dim$p +2 * r) # mean = rep(mui, mi)
  }
  #frame$nlfun(param$alpha, frame$data[[i]]$X_i)

  else if (frame$family == "Bin")
    ui <- (2/dmssmn(x = yi, mu = frame$nlfun(param$alpha, frame$data[[i]]$X_i), Omegabar = MOmegabar(omega = Sigmai, lambda = lambda_yi), df = param$nu, delta = Mdelta(omega = Sigmai, lambda = lambda_yi), mixvar = "Bin")) *
    (param$nu[1] * param$nu[2]^r *dmnorm(x = yi, mean = frame$nlfun(param$alpha, frame$data[[i]]$X_i), Var = param$nu[2]^-1 * Sigmai) * pnorm(sqrt(param$nu[2]) * Ai) + (1 - param$nu[1]) * dmnorm(x = yi, mean = frame$nlfun(param$alpha, frame$data[[i]]$X_i), Var = Sigmai)* pnorm(Ai))


  else if (frame$family == "Beta"){
    ui <- (dmssmn(x = yi, mu = frame$nlfun(param$alpha, frame$data[[i]]$X_i), Omegabar = MOmegabar(omega = Sigmai, lambda = lambda_yi), df = param$nu, delta = rep(0, mi), mixvar = "Beta")/dmssmn(x = yi, mu = frame$nlfun(param$alpha, frame$data[[i]]$X_i), Omegabar = MOmegabar(omega = Sigmai, lambda = lambda_yi), df = param$nu, delta = Mdelta(omega = Sigmai, lambda = lambda_yi), mixvar = "Beta"))
    # * ((2 * gamma((2 * param$nu + frame$dim$p + 2 * r)/2)/ gamma((2 * param$nu + frame$dim$p)/2))) *
    #(2/(t(ri) %*% solve(Sigmai) %*% ri))^r *
     # pgamma()
  }
  else
    ui <- 1

  if (!is.finite(ui)) warning("'ui' return infinite")
    #ui <- 1

  return(ui)
}


#' tauir
#' @description page 4 et 5
#' @export

tauir <- function(frame, param, Sigmai, lambda_yi, ri, Ai, r = 1, i = ind)
{
  yi <- frame$data[[i]]$Y_i
  mi <- length(yi)

  if(frame$family == "Gamma"){
    #taui <- (ddmst(dat = yi, n = 1, p = mi, mean = rep(mui, mi), cov = Sigmai, nu = param$nu)/ddmst(yi, 1, mi, rep(mui, mi), Sigmai, param$nu, del = lambda_yi)) * (2^((r+1)/2) * gamma((param$nu + mi + r)/2) / (sqrt(pi) * gamma((param$nu + mi)/2))) * (((param$nu + t(ri) %*% solve(Sigmai) %*% ri)^((param$nu + mi)/2))/ ((param$nu + t(ri) %*% solve(Sigmai) %*% ri + Ai^2)^((param$nu + mi + r)/2)))
    taui <- (dmssmn(x = yi, mu = frame$nlfun(param$alpha, frame$data[[i]]$X_i), Omegabar = MOmegabar(omega = Sigmai, lambda = lambda_yi), df = param$nu, delta = rep(0, mi), mixvar = "Gamma")/dmssmn(x = yi, mu = frame$nlfun(param$alpha, frame$data[[i]]$X_i), Omegabar = MOmegabar(omega = Sigmai, lambda = lambda_yi), df = param$nu, delta = Mdelta(omega = Sigmai, lambda = lambda_yi), mixvar = "Gamma")) * (2^((r+1)/2) * gamma((param$nu + mi + r)/2) / (sqrt(pi) * gamma((param$nu + mi)/2))) * (((param$nu + t(ri) %*% solve(Sigmai) %*% ri)^((param$nu + mi)/2))/ ((param$nu + t(ri) %*% solve(Sigmai) %*% ri + Ai^2)^((param$nu + mi + r)/2)))

  }
  else if(frame$family == "Bin")
    taui <- (2/dmssmn(x = yi, mu = frame$nlfun(param$alpha, frame$data[[i]]$X_i), Omegabar = MOmegabar(omega = Sigmai, lambda = lambda_yi), df = param$nu, delta = Mdelta(omega = Sigmai, lambda = lambda_yi), mixvar = "Bin")) * (param$nu[1] * param$nu[2]^r *dmnorm(yi, frame$nlfun(param$alpha, frame$data[[i]]$X_i), Var = param$nu[2]^-1 * Sigmai)* dnorm(sqrt(param$nu[2]) * Ai) + (1 - param$nu[1])* dmnorm(yi, frame$nlfun(param$alpha, frame$data[[i]]$X_i), Sigmai)* dnorm(Ai))

  #else if(frame$family == "Beta")
    #taui <- ll
  else
    taui <- dnorm(Ai)/pnorm(Ai)

  if (!is.finite(taui)) warning("'taui returns infinite.")
    #taui <- 1e-10

  return(taui)
}

#' Alpha, fixed effect function
#' @description eq 17
#' @export
ALPHA <- function(alpha = param$alpha, mdata, frame){
  sum(sapply(1:frame$dim$n , function(i){
    c(mdata[[i]]$ui/2 * t(frame$data[[i]]$Y_i - frame$nlfun(alpha, frame$data[[i]]$X_i)) %*% (frame$data[[i]]$Y_i - frame$nlfun(alpha, frame$data[[i]]$X_i)) - t(frame$data[[i]]$Y_i - frame$nlfun(alpha, frame$data[[i]]$X_i)) %*% frame$data[[i]]$Zi %*% mdata[[i]]$ubi)
  }))
}

gyi <- function(nu, param, mdata, frame){
  if (any(frame$family == c("Bin", "Gamma", "Beta")))
    res <- -sum(sapply(1:frame$dim$n , function(i){
      Sigmai = mdata[[i]]$Sigmai
      lambda_yi = mdata[[i]]$lambda_yi

      c(fnlmer::dmssmn(x = frame$data[[i]]$Y_i, mu = frame$nlfun(phi = param$alpha, X = frame$data[[i]]$X_i), Omegabar = MOmegabar(omega = Sigmai, lambda = lambda_yi), df = nu, delta = Mdelta(omega = Sigmai, lambda = lambda_yi), mixvar = frame$family))
    }))

  else if (frame$family == "Bin00") # Not to use
    res <- dmnorm(x = frame$data[[i]]$Y_i, mean = frame$nlfun(phi = param$alpha, X = frame$data[[i]]$X_i), Var = (1/mdata[[i]]$ui) * mdata[[i]]$Sigmai) * c(pnorm(sqrt(mdata[[i]]$ui) * t(mdata[[i]]$lambda_yi) %*% solve(sqrtm(mdata[[i]]$Sigmai)) %*% (frame$data[[i]]$Y_i - frame$nlfun(phi = param$alpha, X = frame$data[[i]]$X_i)))) * dbinom(x, 1, 1 - nu[1])

  #print(res)
  return(res)
}


