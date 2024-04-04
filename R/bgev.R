## Generalised extreme value negative log-likelihood functions

.bgev.d0 <- function(pars, likdata) {
  out <- bgevd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$other)
  out
}

.bgev.d12 <- function(pars, likdata) {
  out <- bgevd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$other)
  out
}

.bgev.d34 <- function(pars, likdata) {
  out <- bgevd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate, likdata$other)
  out
}

.bgevfns <- list(d0 = .bgev.d0, d120 = .bgev.d12, d340 = .bgev.d34)

.bgev_unlink <- list(NULL, function(x) exp(x), function(x) 1.5 / (1 + exp(-x)) - .5)
attr(.bgev_unlink[[2]], "deriv") <- .bgev_unlink[[2]]
attr(.bgev_unlink[[3]], "deriv") <- function(x) 1.5 * exp(-x)/(1 + exp(-x))^2

.qbgev <- function(p, qa, sbeta, xi, psuba, psubb, alpha, beta) {
  hbeta <- .5 * beta
  numer11 <- sbeta / (.ell1(1 - hbeta, xi) - .ell1(hbeta, xi))
  gev_pars <- bgev2gev(qa, sbeta, xi, psuba, psubb, alpha, beta)
  inits <- .qgev(p, gev_pars[ ,1], gev_pars[, 2], gev_pars[, 3])
  .root_NR(inits, .pbgev, p, eps = 1e-5, maxit = 5, tol = 1e-3, 
               qa = qa, sbeta = sbeta, xi = xi, psuba = psuba, psubb = psubb,
               alpha = alpha, beta = beta, log = TRUE)
}

.min_d0 <- function(x, f, p, ...) {
  (f(x, ...) - log(p))^2
}

.root_d0 <- function(x, f, p, ...) {
  f(x, ...) - log(p)
}

.root_d01 <- function(x, f, p, eps = 1e-5, ...) {
  f0 <- .root_d0(x, f, p, ...)
  fh <- .root_d0(x + eps, f, p, ...)
  fl <- .root_d0(x - eps, f, p, ...)
  cbind(f0, .5 * (fh - fl) / eps)
}

.min_d12 <- function(x, f, p, eps = 1e-5, ...) {
  f0 <- .min_d0(x, f, p, ...)
  fh <- .min_d0(x + eps, f, p, ...)
  fl <- .min_d0(x - eps, f, p, ...)
  g <- .5 * (fh - fl) / eps
  h <- (fh + fl - 2 * f0) / (eps * eps)
  out <- -g / h
  attr(out, 'g') <- g
  attr(out, 'h') <- h
  out
}

.root_NR <- function(x, f, p, eps, tol = 1e-3, maxit = 10, ...) {
  cond <- TRUE
  it <- 0
  while(cond) {
    fg <- .root_d01(x, f, p, abs(x) * eps, ...)
    x <- x - fg[, 1] / fg[, 2]
    cond <- it < maxit & any(abs(fg[, 1]) > tol)
  }
  x
}

.min_newton <- function(x, f, p, eps, maxit = 10, ...) {
  cond <- TRUE
  it <- 0
  while(cond) {
    stp <- .min_d12(x, f, p, eps, ...)
    x <- x + stp
    cond <- it < maxit & any(abs(attr(stp, 'g')) > 1e-5)
  }
  x
}
  
  
.bgevfns$q <- .qbgev
.bgevfns$unlink <- .bgev_unlink

.ell1 <- function(a, xi) (-log(a))^(-xi)
.ell2 <- function(a) log(-log(a))

.pbgev <- function(x, qa, sbeta, xi, psuba, psubb, alpha, beta, NAOK = FALSE, log = FALSE) {

  hbeta <- .5 * beta
  lh1 <- (x - qa) * (.ell1(1 - hbeta, xi) - .ell1(hbeta, xi)) / sbeta + .ell1(alpha, xi)
  lh1 <- pmax(lh1, 0)^(-1/xi)
  
  numer11 <- sbeta / (.ell1(1 - hbeta, xi) - .ell1(hbeta, xi))
  iFa <- qa + (.ell1(psuba, xi) - .ell1(alpha, xi)) * numer11
  iFb <- qa + (.ell1(psubb, xi) - .ell1(alpha, xi)) * numer11
  
  denom21 <- .ell2(psuba) - .ell2(psubb)
  numer21 <- iFb - iFa
  numer22 <- .ell2(hbeta) - .ell2(1 - hbeta)
  tqa <- iFa - numer21 * (.ell2(alpha) - .ell2(psuba)) / denom21
  tsbeta <- numer21 * numer22 / denom21
  
  lh2 <- exp(-((x - tqa) * numer22 / tsbeta - .ell2(alpha)))
  
  px <- pbeta(x / (iFb - iFa), 5, 5)
  sx <- pbeta(x / (iFb - iFa), 5, 5, lower.tail = FALSE)
  
  out <- - px * lh1 - sx * lh2
  
  if (!log)
    out <- exp(out)
  
  out

}

bgev2gev <- function(qa, sbeta, xi, psuba, psubb, alpha, beta) {

  hbeta <- .5 * beta
  numer11 <- sbeta / (.ell1(1 - hbeta, xi) - .ell1(hbeta, xi))
  mu <- qa - numer11 * (.ell1(alpha, xi) - 1)
  sigma <- xi * numer11
  cbind(mu, sigma, xi)
  
}

