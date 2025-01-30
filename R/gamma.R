## Gamma negative log-likelihood functions

.gamma.d0 <- function(pars, likdata) {
likdata$y <- as.matrix(likdata$y)
nhere <- rowSums(is.finite(likdata$y))  
out <- gammad0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
out
}

.gamma.d12 <- function(pars, likdata, sandwich = FALSE) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- gammad12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.gamma.d34 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- gammad34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.gammafns <- list(d0=.gamma.d0, d120=.gamma.d12, d340=.gamma.d34)

.qgamma <- function(x, scale, shape) {
qgamma(x, shape = shape, scale = scale) 
}

.pgamma <- function(x, scale, shape) {
pgamma(x, shape = shape, scale = scale) 
}

.gamma_unlink <- list(function(x) exp(x), function(x) exp(x))
attr(.gamma_unlink[[1]], "deriv") <- .gamma_unlink[[1]]
attr(.gamma_unlink[[2]], "deriv") <- .gamma_unlink[[2]]

.gammafns$q <- .qgamma
.gammafns$unlink <- .gamma_unlink
