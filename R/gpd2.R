## Transformed shape GPD negative log-likelihood functions

.gpd2.d0 <- function(pars, likdata) {
likdata$y <- as.matrix(likdata$y)
nhere <- rowSums(is.finite(likdata$y))  
out <- gpd2d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
out
}

.gpd2.d12 <- function(pars, likdata, sandwich = FALSE) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- gpd2d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.gpd2.d34 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- gpd2d34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.gpd2fns <- list(d0=.gpd2.d0, d120=.gpd2.d12, d340=.gpd2.d34)

.qgpd2 <- function(x, scale, shape, logscale = NULL, tshape = NULL) {
if (!missing(logscale)) 
  scale <- exp(logscale)
if (!missing(tshape)) 
  shape <- 1.5 / (1 + exp(-tshape)) - 1
scale * ((1 - x)^(-shape) - 1)/shape
}

.pgpd2 <- function(x, scale, shape, logscale = NULL, tshape=NULL) {
if (!missing(logscale)) 
  scale <- exp(logscale)
if (!missing(tshape)) 
  shape <- 1.5 / (1 + exp(-tshape)) - 1
1 - (1 + shape * (x / scale))^(-1/shape)
}

.dqgpd2 <- function(x, lscale, tshape) {
  .e2 <- exp(-tshape)
  .e3 <- 1 + .e2
  .e5 <- 1.5/.e3 - 1
  .e6 <- 1 - x
  .e7 <- .e6^.e5
  .e9 <- 1/.e7 - 1
  .e10 <- exp(logscale)
  cbind(logscale = .e9 * .e10/.e5, 
        tshape = -(1.5 * ((.e9/.e5 + log(.e6)/.e7) * .e2 * .e10/(.e3^2 * .e5))))
}

.gpd2_unlink <- list(function(x) exp(x), function(x) 1.5 / (1 + exp(-x)) - 1.0)
attr(.gpd2_unlink[[1]], "deriv") <- .gpd2_unlink[[1]]
attr(.gpd2_unlink[[2]], "deriv") <- function(x) 1.5 * exp(-x)/(1 + exp(-x))^2

.gpd2fns$q <- .qgpd2
.gpd2fns$unlink <- .gpd2_unlink
