## Weibull negative log-likelihood functions

.weib.d0 <- function(pars, likdata) {
likdata$y <- as.matrix(likdata$y)
nhere <- rowSums(is.finite(likdata$y))  
out <- weibd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
out
}

.weib.d12 <- function(pars, likdata, sandwich = FALSE) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- weibd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.weib.d34 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- weibd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.weibfns <- list(d0=.weib.d0, d120=.weib.d12, d340=.weib.d34)

.qweibull <- function(x, scale, shape, logscale=NULL, logshape=NULL) {
if (!missing(logscale)) scale <- exp(logscale)
if (!missing(logshape)) shape <- exp(logshape)
scale * ((-log(1 - x)) ^ (1/shape))
}

.pweibull <- function(x, scale, shape, logscale=NULL, logshape=NULL) {
if (!missing(logscale)) scale <- exp(logscale)
if (!missing(logshape)) shape <- exp(logshape)
1 - exp(-(x / scale)^shape)
}

.dqweibull <- function(x, lscale, lshape) {
.e1 <- -log(1 - x)
.e2 <- exp(lshape)
.e4 <- .e1^(1/.e2) * exp(lscale)
d1 <- .e4
d2 <--(.e4 * log(.e1)/.e2)
cbind(d1, d2)
}
