#' The blended generalised extreme value distribution
#'
#' Density, distribution function, quantile function and random generation for 
#' the blended generalised extreme value distribution with parameters 
#' \code{location}, \code{scale}, \code{shape}, \code{pa}, \code{pb}, 
#' \code{alpha} and \code{beta}.
#' 
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param location,scale,shape location, scale and shape parameters.
#' @param pa,pb,alpha,beta Gumbel to GEV mixing parameters; see Details.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are 
#' \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.
#' 
#' @details
#' 
#' The blended generalised extreme value distribution with \code{location} 
#' parameter \eqn{q_\alpha}, \code{scale} parameter \eqn{s_\beta} and 
#' \code{shape} parameter \eqn{\xi} has cumulative distribution function
#' given by
#' \deqn{H(x \mid q_\alpha, s_\beta, \xi) = F(x \mid q_\alpha, s_\beta, \xi)^{p(x')} G(x \mid \tilde q_\alpha \tilde s_\beta)^{1 - p(x')}}
#' where 
#' \deqn{F(x \mid q_\alpha, s_\beta, \xi) = \exp \left\{ - \left[ \dfrac{x - q_\alpha}{s_\beta(\ell_{1 - \beta / 2, \xi} - \ell_{\beta / 2, \xi})^{-1}} + \ell_{\alpha, \xi} \right]^{-1/\xi}_+  \right\}}
#' with \eqn{\ell_{a, \xi} = (-\log a)^{-\xi}}, \eqn{[x]_+ = \max(0, x)},
#' \code{alpha} and \code{beta} parameters \eqn{\alpha} and \eqn{\beta}, 
#' respectively, 
#' \deqn{G(x \mid q_{\tilde \alpha}, s_{\tilde \beta}) = \exp\left\{ -\exp\left[- \left\{\dfrac{x - \tilde q_\alpha}{\tilde s_\beta(\ell_{1 - \beta / 2} - \ell_{\beta / 2})^{-1}} + \ell_{\alpha} \right\}\right]\right\}}
#' with \eqn{\ell_a = \log(-\log a)},
#' \deqn{\tilde q_\alpha = a - \dfrac{(b - a)(\ell_\alpha - \ell_{p_a}))}{\ell_{p_a} - \ell_{p_b}}~~\text{and}~~\tilde s_\beta = \dfrac{(b - a)(\ell_{\beta / 2} - \ell_{1 - \beta/2})}{\ell_{p_a} - \ell_{p_b}},}
#' with \eqn{a = F^{-1}(p_a \mid q_\alpha, s_\beta, \xi)}, 
#' \eqn{b = F^{-1}(p_b \mid q_\alpha, s_\beta, \xi)}, with \code{pa} and \code{pb}
#' parameters \eqn{p_a} and \eqn{p_b}, respectively, and where 
#' \eqn{x' = (x - a) / (b - a)} and \eqn{p(x)} denotes the cumulative distribution 
#' of the Beta(5, 5) distribution.
#' 
#' Default values for \code{pa}, \code{pb}, \code{alpha} and \code{beta}
#' are taken from Castro-Camilo et al. (2022).
#' 
#' @references 
#' 
#' Castro-Camilo, D., Huser, R., & Rue, H. (2022). Practical strategies for 
#' generalized extreme value-based regression models for extremes. 
#' Environmetrics, 33(6), e2742. \doi{10.1002/env.2742}
#'
#' @examples
#'
#' dbgev(3, 2, 1, .1)
#'
#' @seealso \link{predict.evgam}
#'
#' @return 
#' 
#' \code{dbgev} gives the density,
#' \code{pbgev} gives the distribution function,
#' \code{qbgev} gives the quantile function, and
#' \code{rbgev} generates random deviates.
#' 
#' @export
#' 
dbgev <- function(x, location, scale, shape, pa = 0.05, pb = 0.2, alpha = 0.5, 
                  beta = 0.5, log = FALSE) {

  # hbeta <- .5 * beta
  # ee1 <- .ell1(1 - hbeta, shape)
  # ee2 <- .ell1(hbeta, shape)
  # z1d <- (ee1 - ee2) / scale
  # t4 <- .ell1(alpha, shape)
  # z1 <- (x - location) * z1d + t4
  # t1 <- pmax(z1, 0)^(-1/shape)
  # lF <- - t1
  # lf <- - log(shape) + lF - log(1 + 1 / shape) * z1 + log(t1 - t1) - log(scale)  
  # 
  # numer11 <- scale / (.ell1(1 - hbeta, shape) - .ell1(hbeta, shape))
  # t5 <- .ell1(pa, shape)
  # t6 <- .ell1(pb, shape)
  # iFa <- location + (t5 - t4) / z1d
  # iFb <- location + (t6 - t4) / z1d
  # 
  # t7 <- .ell2(pa)
  # t8 <- .ell2(pb)
  # denom21 <- t7 - t8
  # numer21 <- iFb - iFa
  # numer22 <- .ell2(hbeta) - .ell2(1 - hbeta)
  # t9 <- .ell2(alpha)
  # tqa <- iFa - numer21 * (t9 - t7) / denom21
  # tsbeta <- numer21 * numer22 / denom21
  # z2d <- numer22 / tsbeta
  # 
  # z2 <- (x - tqa) * z2d - t9
  # t2 <- exp(-z2)
  # lG <- - t2
  # # lg <- lG - z2x + log(.ell2(hbeta) - .ell2(1 - hbeta)) + log(denom21) - log(numer21)
  # 
  # t10 <- (iFb - iFa)
  # x2 <- x / t10
  # px <- pbeta(x2, 5, 5)
  # sx <- pbeta(x2, 5, 5, lower.tail = FALSE)
  # 
  # out1 <- px * lF + sx * lG
  # 
  # dx <- dbeta(x2, 5, 5)
  # out21 <- lF * dx / t10
  # out21[!is.finite(out21)] <- 0
  # fiF <- z1^(-(1 + 1 / shape)) * z1d / shape
  # out22 <- px * fiF
  # out23 <- - dx * lG
  # giG <- t2 * z2d
  # out24 <- sx * giG
  # 
  # out2 <- log(out21 + out22 + out23 + out24)
  # 
  # out <- out1 + out2
  
  hbeta <- .5 * beta
  
  a <- .qgev2(pa, location, scale, shape, alpha, beta)
  b <- .qgev2(pb, location, scale, shape, alpha, beta)
  x2 <- (x - a) / (b - a)
  
  px <- pbeta(x2, 5, 5)
  sx <- pbeta(x2, 5, 5, lower.tail = FALSE)
  
  lF <- .pgev2(x, location, scale, shape, alpha, beta, log = TRUE)
  
  tlocation <- a - (b - a) * (.ell2(alpha) - .ell2(pa)) / (.ell2(pa) - .ell2(pb))
  tscale <- (b - a) * (.ell2(hbeta) - .ell2(1 - hbeta)) / (.ell2(pa) - .ell2(pb))
  
  lG <- .pgumb2(x, tlocation, tscale, alpha, beta, log = TRUE)
  
  out <- px * lF + sx * lG
  
  dx <- dbeta(x2, 5, 5) / (b - a)
  
  out2 <- dx * lF
  
  f <- .dgev2(x, location, scale, shape, alpha, beta)
  F <- .pgev2(x, location, scale, shape, alpha, beta)
  
  out2 <- out2 + px * f / F
  
  out2 <- out2 - dx * lG
  
  g <- .dgumb2(x, tlocation, tscale, alpha, beta)
  G <- .pgumb2(x, tlocation, tscale, alpha, beta)
  
  out2 <- out2 + sx * g / G
  
  out <- out + log(out2)

  if (!log)
    out <- exp(out)
  
  out
  
}

.pgev2 <- function(x, location, scale, shape, alpha, beta, log = FALSE) {
  hbeta <- .5 * beta
  x <- (x - location) * (.ell1(1 - hbeta, shape) - .ell1(hbeta, shape)) / scale
  z1 <- pmax(x + .ell1(alpha, shape), 0)
  t1 <- z1^(-1/shape)
  out <- - t1
  out[out == -Inf] <- -1e20
  if (!log)
    out <- exp(out)
  out
}

.qgev2 <- function(p, location, scale, shape, alpha, beta, log = FALSE) {
  hbeta <- .5 * beta
  z1d <- (.ell1(1 - hbeta, shape) - .ell1(hbeta, shape)) / scale
  location + (.ell1(p, shape) - .ell1(alpha, shape)) / z1d
}    

.dgev2 <- function(x, location, scale, shape, alpha, beta, log = FALSE) {
  hbeta <- .5 * beta
  lz1d <- log(.ell1(1 - hbeta, shape) - .ell1(hbeta, shape)) - log(scale)
  x <- (x - location) * exp(lz1d)
  z1 <- pmax(x + .ell1(alpha, shape), 0)
  t1 <- z1^(-1/shape)
  out <- - t1
  out <- out - log(shape)
  out <- out - (1 + 1 / shape) * log(z1)
  out <- out + lz1d
  if (!log)
    out <- exp(out)
  out
}

.dgumb2 <- function(x, location, scale, alpha, beta, log = FALSE) {
  hbeta <- .5 * beta
  lz2d <- log(.ell2(hbeta) - .ell2(1 - hbeta)) - log(scale)
  z2 <- (x - location) * exp(lz2d) - .ell2(alpha)
  t2 <- exp(-z2)
  out <- - t2 - z2 + lz2d
  if (!log)
    out <- exp(out)
  out
}

.pgumb2 <- function(x, location, scale, alpha, beta, log = FALSE) {
  hbeta <- .5 * beta
  lz2d <- log(.ell2(hbeta) - .ell2(1 - hbeta)) - log(scale)
  z2 <- (x - location) * exp(lz2d) - .ell2(alpha)
  t2 <- exp(-z2)
  out <- - t2
  if (!log)
    out <- exp(out)
  out
}

# .pbgev2 <- function(x, location, scale, shape, pa = 0.05, pb = 0.2, alpha = 0.5, 
#                   beta = 0.5, log = FALSE) {
#   
#   hbeta <- .5 * beta
#   
#   a <- .qgev2(pa, location, scale, shape, alpha, beta)
#   b <- .qgev2(pb, location, scale, shape, alpha, beta)
#   x2 <- (x - a) / (b - a)
#   
#   px <- pbeta(x2, 5, 5)
#   sx <- pbeta(x2, 5, 5, lower.tail = FALSE)
#   
#   lF <- .pgev2(x, location, scale, shape, alpha, beta, log = TRUE)
#   
#   tlocation <- a - (b - a) * (.ell2(alpha) - .ell2(pa)) / (.ell2(pa) - .ell2(pb))
#   tscale <- (b - a) * (.ell2(hbeta) - .ell2(1 - hbeta)) / (.ell2(pa) - .ell2(pb))
#     
#   lG <- .pgumb2(x, tlocation, tscale, alpha, beta, log = TRUE)
#   
#   out <- px * lF + sx * lG
#   
#   if (!log)
#     out <- exp(out)
#   
#   out
# }  
# 
# .dbgev2 <- function(x, location, scale, shape, pa = 0.05, pb = 0.2, alpha = 0.5, 
#                     beta = 0.5, log = FALSE) {
#   
#   hbeta <- .5 * beta
#   
#   a <- .qgev2(pa, location, scale, shape, alpha, beta)
#   b <- .qgev2(pb, location, scale, shape, alpha, beta)
#   x2 <- (x - a) / (b - a)
#   
#   px <- pbeta(x2, 5, 5)
#   sx <- pbeta(x2, 5, 5, lower.tail = FALSE)
#   
#   lF <- .pgev2(x, location, scale, shape, alpha, beta, log = TRUE)
#   
#   tlocation <- a - (b - a) * (.ell2(alpha) - .ell2(pa)) / (.ell2(pa) - .ell2(pb))
#   tscale <- (b - a) * (.ell2(hbeta) - .ell2(1 - hbeta)) / (.ell2(pa) - .ell2(pb))
#   
#   lG <- .pgumb2(x, tlocation, tscale, alpha, beta, log = TRUE)
#   
#   out <- px * lF + sx * lG
#   
#   dx <- dbeta(x2, 5, 5) / (b - a)
#   
#   out2 <- dx * lF
# 
#   f <- .dgev2(x, location, scale, shape, alpha, beta)
#   F <- .pgev2(x, location, scale, shape, alpha, beta)
#   
#   out2 <- out2 + px * f / F
#   
#   out2 <- out2 - dx * lG
#   
#   g <- .dgumb2(x, tlocation, tscale, alpha, beta)
#   G <- .pgumb2(x, tlocation, tscale, alpha, beta)
#   
#   out2 <- out2 + sx * g / G
#   
#   out <- out + log(out2)
#   
#   if (!log)
#     out <- exp(out)
#   
#   out
# }  

#' @rdname dbgev
#' 
#' @export
#' 
pbgev <- function(q, location, scale, shape, pa = 0.05, pb = 0.2, alpha = 0.5, 
               beta = 0.5, lower.tail = TRUE, log.p = FALSE) {

  # hbeta <- .5 * beta
  # t1 <- .ell1(1 - hbeta, shape)
  # t2 <- .ell1(hbeta, shape)
  # t3 <- (t1 - t2) / scale
  # t4 <- .ell1(alpha, shape)
  # lh1 <- (q - location) * t3 + t4
  # lh1 <- pmax(lh1, 0)^(-1/shape)
  # 
  # numer11 <- scale / (.ell1(1 - hbeta, shape) - .ell1(hbeta, shape))
  # t5 <- .ell1(pa, shape)
  # t6 <- .ell1(pb, shape)
  # iFa <- location + (t5 - t4) / t3
  # iFb <- location + (t6 - t4) / t3
  # 
  # t7 <- .ell2(pa)
  # t8 <- .ell2(pb)
  # denom21 <- t7 - t8
  # numer21 <- iFb - iFa
  # numer22 <- .ell2(hbeta) - .ell2(1 - hbeta)
  # t9 <- .ell2(alpha)
  # tqa <- iFa - numer21 * (t9 - t7) / denom21
  # tsbeta <- numer21 * numer22 / denom21
  # 
  # lh2 <- exp(-((q - tqa) * numer22 / tsbeta - t9))
  # 
  # q2 <- q / (iFb - iFa)
  # px <- pbeta(q2, 5, 5)
  # sx <- pbeta(q2, 5, 5, lower.tail = FALSE)
  # 
  # out <- - px * lh1 - sx * lh2
  
  hbeta <- .5 * beta
  
  a <- .qgev2(pa, location, scale, shape, alpha, beta)
  b <- .qgev2(pb, location, scale, shape, alpha, beta)
  x2 <- (q - a) / (b - a)
  
  px <- pbeta(x2, 5, 5)
  sx <- pbeta(x2, 5, 5, lower.tail = FALSE)
  
  lF <- .pgev2(q, location, scale, shape, alpha, beta, log = TRUE)
  
  tlocation <- a - (b - a) * (.ell2(alpha) - .ell2(pa)) / (.ell2(pa) - .ell2(pb))
  tscale <- (b - a) * (.ell2(hbeta) - .ell2(1 - hbeta)) / (.ell2(pa) - .ell2(pb))
  
  lG <- .pgumb2(q, tlocation, tscale, alpha, beta, log = TRUE)
  
  out <- px * lF + sx * lG
  
  if (log.p) {
    if (!lower.tail) {
      out <- log(1 - exp(out))
    }
  } else {
    out <- exp(out)
    if (!lower.tail)
      out <- 1 - out
  }

  out
  
}

#' @rdname dbgev
#' 
#' @export
#' 
qbgev <- function(p, location, scale, shape, pa = 0.05, pb = 0.2, alpha = 0.5, 
      beta = 0.5, lower.tail = TRUE, log.p = FALSE) {
  if (log.p & lower.tail)
    p <- exp(p)
  if (!lower.tail) {
    if (log.p) {
      p <- 1 - exp(p)
    } else {
      p <- 1 - p
    }
  }
  hbeta <- .5 * beta
  numer11 <- scale / (.ell1(1 - hbeta, shape) - .ell1(hbeta, shape))
  gev_pars <- bgev2gev(location, scale, shape, pa, pb, alpha, beta)
  inits <- .qgev(p, gev_pars[[1]], gev_pars[[2]], gev_pars[[3]])
  # out <- .root_NR(inits, .pbgev, p, eps = 1e-5, maxit = 5, tol = 1e-3, 
  #          qa = location, sbeta = scale, xi = shape, psuba = pa, psubb = pb,
  #          alpha = alpha, beta = beta, log = TRUE)
  out <- .root_NR(inits, pbgev, p, eps = 1e-5, maxit = 5, tol = 1e-3, 
                  location = location, scale = scale, shape = shape, pa = pa, 
                  pb = pb, alpha = alpha, beta = beta, log.p = TRUE)
  out
}

#' @rdname dbgev
#' 
#' @export
#' 
rbgev <- function(n, location, scale, shape, pa = 0.05, pb = 0.2, alpha = 0.5, 
      beta = 0.5) {
p <- runif(n)
qbgev(p, location, scale, shape, pa, pb, alpha, beta)
}

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

# .qbgev <- function(p, qa, sbeta, xi, psuba, psubb, alpha, beta) {
#   hbeta <- .5 * beta
#   numer11 <- sbeta / (.ell1(1 - hbeta, xi) - .ell1(hbeta, xi))
#   gev_pars <- bgev2gev(qa, sbeta, xi, psuba, psubb, alpha, beta)
#   inits <- .qgev(p, gev_pars[[1]], gev_pars[[2]], gev_pars[[3]])
#   .root_NR(inits, .pbgev, p, eps = 1e-5, maxit = 5, tol = 1e-3, 
#                qa = qa, sbeta = sbeta, xi = xi, psuba = psuba, psubb = psubb,
#                alpha = alpha, beta = beta, log = TRUE)
# }

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
  list(f0, .5 * (fh - fl) / eps)
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

.root_NR <- function(x, f, p, eps, tol = 1e-5, maxit = 10, ...) {
  cond <- TRUE
  it <- 0
  while(cond) {
    fg <- .root_d01(x, f, p, abs(x) * eps, ...)
    x <- x - fg[[1]] / fg[[2]]
    cond <- it < maxit & any(abs(fg[[1]]) > tol)
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
  
  
.bgevfns$q <- qbgev
.bgevfns$unlink <- .bgev_unlink

.ell1 <- function(a, xi) (-log(a))^(-xi)
.ell2 <- function(a) log(-log(a))

# .pbgev <- function(x, qa, sbeta, xi, psuba, psubb, alpha, beta, NAOK = FALSE, log = FALSE) {
# 
#   hbeta <- .5 * beta
#   lh1 <- (x - qa) * (.ell1(1 - hbeta, xi) - .ell1(hbeta, xi)) / sbeta + .ell1(alpha, xi)
#   lh1 <- pmax(lh1, 0)^(-1/xi)
#   
#   numer11 <- sbeta / (.ell1(1 - hbeta, xi) - .ell1(hbeta, xi))
#   iFa <- qa + (.ell1(psuba, xi) - .ell1(alpha, xi)) * numer11
#   iFb <- qa + (.ell1(psubb, xi) - .ell1(alpha, xi)) * numer11
#   
#   denom21 <- .ell2(psuba) - .ell2(psubb)
#   numer21 <- iFb - iFa
#   numer22 <- .ell2(hbeta) - .ell2(1 - hbeta)
#   tqa <- iFa - numer21 * (.ell2(alpha) - .ell2(psuba)) / denom21
#   tsbeta <- numer21 * numer22 / denom21
#   
#   lh2 <- exp(-((x - tqa) * numer22 / tsbeta - .ell2(alpha)))
#   
#   px <- pbeta(x / (iFb - iFa), 5, 5)
#   sx <- pbeta(x / (iFb - iFa), 5, 5, lower.tail = FALSE)
#   
#   out <- - px * lh1 - sx * lh2
#   
#   if (!log)
#     out <- exp(out)
#   
#   out
# 
# }

#' Conversion between blended generalised extreme value distribution and 
#' generalised extreme value distribution parameters
#'
#' Conversion between location, scale and shape parameters of the generalised
#' extreme value distribution and corresponding parameters of blended 
#' generalised extreme value distribution in both directions.
#' 
#' @param location,scale,shape location, scale and shape parameters.
#' @param pa,pb,alpha,beta Gumbel to GEV mixing parameters; see Details.
#' @param simplify logical; should \code{simplify2array()} be called at the end?
#' 
#' @details
#' 
#' The blended generalised extreme value distribution has \code{location} 
#' parameter \eqn{q_\alpha}, \code{scale} \eqn{s_\beta} and 
#' \code{shape} \eqn{\xi} parameters; see, e.g., \code{dbgev}. The generalised 
#' extreme value distribution's location and scale parameters in its typical 
#' form are given by \eqn{\mu = q_\alpha - s_\beta (\ell_{\alpha, \xi} - 1) / (\ell_{1 - \beta/2, \xi} - \ell_{\beta/2, \xi})}
#' and \eqn{\sigma = \xi s_\beta / (\ell_{1 - \beta/2, \xi} - \ell_{\beta/2, \xi})}. 
#' So \code{bgev2gev} gives the mapping \eqn{(q_\alpha, s_\beta, \xi) \mapsto (\mu, \sigma, \xi))}.
#' The reverse mapping, \eqn{(\mu, \sigma, \xi) \mapsto (q_\alpha, s_\beta, \xi))},
#' is given by \code{gev2bgev}.
#' 
#' @references 
#' 
#' Castro-Camilo, D., Huser, R., & Rue, H. (2022). Practical strategies for 
#' generalized extreme value-based regression models for extremes. 
#' Environmetrics, 33(6), e2742. \doi{10.1002/env.2742}
#'
#' @examples
#'
#' bgev2gev(2, 1, .1)
#' gev2bgev(2, 1, .1)
#'
#' @seealso \link{dbgev}
#'
#' @return 
#' 
#' A \code{list} or \code{array} if \code{simplify = TRUE}.
#' 
#' @export
#' 
bgev2gev <- function(location, scale, shape, pa = 0.05, pb = 0.2,
                     alpha = 0.5, beta = 0.5, simplify = FALSE) {

  hbeta <- .5 * beta
  numer11 <- scale / (.ell1(1 - hbeta, shape) - .ell1(hbeta, shape))
  mu <- location - numer11 * (.ell1(alpha, shape) - 1)
  sigma <- shape * numer11
  out <- list(location = mu, scale = sigma, shape = shape)
  if (simplify)
    out <- simplify2array(out)
  out
}

#' @rdname bgev2gev
#' 
#' @export
#' 
gev2bgev <- function(location, scale, shape, pa = 0.05, pb = 0.2,
                     alpha = 0.5, beta = 0.5, simplify = FALSE) {
  
  hbeta <- .5 * beta
  numer <- .ell1(1 - hbeta, shape) - .ell1(hbeta, shape)
  sb <- scale * numer / shape
  qa <- location + sb * (.ell1(alpha, shape) - 1) / numer
  out <- list(location = qa, scale = sb, shape = shape)
  if (simplify)
    out <- simplify2array(out)
  out
}
