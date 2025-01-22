## Three-parameter Weibull negative log-likelihood functions

.weib3.d0 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))
  out <- weib3d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.weib3.d12 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- weib3d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.weib3.d34 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- weib3d34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.weib3fns <- list(d0 = .weib3.d0, d120 = .weib3.d12, d340 = .weib3.d34)

.weib3_unlink <- list(NULL, function(x) exp(x), function(x) exp(x))
attr(.weib3_unlink[[2]], "deriv") <- .weib3_unlink[[2]]
attr(.weib3_unlink[[3]], "deriv") <- .weib3_unlink[[3]]

.weib3fns$q <- NULL
.weib3fns$unlink <- .weib3_unlink

