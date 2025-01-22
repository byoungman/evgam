## Generalised extreme value negative log-likelihood functions

.gevr.d0 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))
  out <- gevrd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.gevr.d12 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- gevrd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.gevr.d34 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- gevrd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.gevrfns <- list(d0 = .gevr.d0, d120 = .gevr.d12, d340 = .gevr.d34)

.gevr_unlink <- list(NULL, function(x) exp(x), function(x) 1.5 / (1 + exp(-x)) - 1.0)
attr(.gevr_unlink[[2]], "deriv") <- .gevr_unlink[[2]]
attr(.gevr_unlink[[3]], "deriv") <- function(x) 1.5 * exp(-x)/(1 + exp(-x))^2

.gevrfns$q <- .qgev
.gevrfns$unlink <- .gevr_unlink

