## gamma3 negative log-likelihood functions

.gamma3.d0 <- function(pars, likdata) {
likdata$y <- as.matrix(likdata$y)
nhere <- rowSums(is.finite(likdata$y))  
out <- gamma3d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
out
}

.gamma3.d12 <- function(pars, likdata, sandwich = FALSE) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- gamma3d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.gamma3.d34 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- gamma3d34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.gamma3fns <- list(d0=.gamma3.d0, d120=.gamma3.d12, d340=.gamma3.d34)

.qgamma3 <- function(x, location, scale, shape) {
location + qgamma(x, shape = shape, scale = scale) 
}

.pgamma3 <- function(x, location, scale, shape) {
pgamma3(x - location, shape = shape, scale = scale) 
}

.gamma3_unlink <- list(NULL, function(x) exp(x), function(x) exp(x))
attr(.gamma3_unlink[[2]], "deriv") <- .gamma3_unlink[[2]]
attr(.gamma3_unlink[[3]], "deriv") <- .gamma3_unlink[[3]]

.gamma3fns$q <- .qgamma3
.gamma3fns$unlink <- .gamma3_unlink
