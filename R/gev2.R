## Generalised extreme value negative log-likelihood functions

.gev2.d0 <- function(pars, likdata) {
likdata$y <- as.matrix(likdata$y)
nhere <- rowSums(is.finite(likdata$y))
if (!likdata$sparse) {
  out <- gev2d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
} else {
  out <- gev2spd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
}
out
}

.gev2.d12 <- function(pars, likdata) {
likdata$y <- as.matrix(likdata$y)
nhere <- rowSums(is.finite(likdata$y))  
if (!likdata$sparse) {
  out <- gev2d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
} else {  
  out <- gev2spd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
}
out
}

.gev2.d34 <- function(pars, likdata) {
likdata$y <- as.matrix(likdata$y)
nhere <- rowSums(is.finite(likdata$y))  
if (!likdata$sparse) {
  out <- gev2d34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
} else {
  out <- gev2spd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
}
out
}

.gev2fns <- list(d0 = .gev2.d0, d120 = .gev2.d12, d340 = .gev2.d34)

.gev2_unlink <- list(NULL, function(x) exp(x), function(x) 1.5 / (1 + exp(-x)) - 1.0)
attr(.gev2_unlink[[2]], "deriv") <- .gev2_unlink[[2]]
attr(.gev2_unlink[[3]], "deriv") <- function(x) 1.5 * exp(-x)/(1 + exp(-x))^2

.gev2fns$q <- .qgev
.gev2fns$unlink <- .gev2_unlink

