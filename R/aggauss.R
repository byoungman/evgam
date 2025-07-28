# The conditional extreme value model

.aggauss.d0 <- function(pars, likdata) {
  if (is.null(likdata$args$C))
    likdata$args$C <- .1
  if (is.null(likdata$args$epsilon))
    likdata$args$epsilon <- .1
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- aggaussd0(split(pars, likdata$idpars), 
                  likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], 
                  likdata$y, likdata$args$weights, 
                  likdata$dupid, likdata$duplicate, nhere, 
                  likdata$args$C, likdata$args$epsilon)
  out
}

.aggauss.d12 <- function(pars, likdata) {
  if (is.null(likdata$args$C))
    likdata$args$C <- .1
  if (is.null(likdata$args$epsilon))
    likdata$args$epsilon <- .1
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- aggaussd12(split(pars, likdata$idpars), 
                   likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], 
                   likdata$y, likdata$args$weights, 
                   likdata$dupid, likdata$duplicate, nhere, 
                   likdata$args$C, likdata$args$epsilon)
out
}

.aggauss.d34 <- function(pars, likdata) {
  if (is.null(likdata$args$C))
    likdata$args$C <- .1
  if (is.null(likdata$args$epsilon))
    likdata$args$epsilon <- .1
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- aggaussd34(split(pars, likdata$idpars), 
                   likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], 
                   likdata$y, likdata$args$weights, 
                   likdata$dupid, likdata$duplicate, nhere, 
                   likdata$args$C, likdata$args$epsilon)
out
}

.aggaussfns <- list(d0 = .aggauss.d0, d120 = .aggauss.d12, d340 = .aggauss.d34)

.aggauss_unlink <- list(function(x) x, function(x) exp(x), function(x) exp(x), function(x) exp(x))
attr(.aggauss_unlink[[1]], "deriv") <- function(x) 1 * x
attr(.aggauss_unlink[[2]], "deriv") <- function(x) exp(x)
attr(.aggauss_unlink[[2]], "deriv") <- function(x) exp(x)
attr(.aggauss_unlink[[4]], "deriv") <- function(x) exp(x)

.aggaussfns$q <- NULL
.aggaussfns$unlink <- .aggauss_unlink

