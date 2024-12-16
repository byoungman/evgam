## logit-Gaussian negative log-likelihood functions

.logitgauss.d0 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  if (!likdata$sparse) {
    out <- logitgaussd0(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], ymat=likdata$y, dupid=likdata$dupid, dcate=likdata$duplicate, nhere = nhere)
  } else {
    out <- logitgaussspd0(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], ymat=likdata$y, dupid=likdata$dupid, dcate=likdata$duplicate, nhere = nhere)
  }
  if (!is.finite(out)) out <- 1e20
  out
}

.logitgauss.d12 <- function(pars, likdata, sandwich = FALSE) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  if (!likdata$sparse) {
    out <- logitgaussd12(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], ymat=likdata$y, dupid=likdata$dupid, dcate=likdata$duplicate, nhere = nhere)
  } else {
    out <- logitgaussspd12(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], ymat=likdata$y, dupid=likdata$dupid, dcate=likdata$duplicate, nhere = nhere)
  }
  out
}

.logitgauss.d34 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  if (!likdata$sparse) {
    out <- logitgaussd34(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], ymat=likdata$y, dupid=likdata$dupid, dcate=likdata$duplicate, nhere = nhere)
  } else {
    out <- logitgaussspd34(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], ymat=likdata$y, dupid=likdata$dupid, dcate=likdata$duplicate, nhere = nhere)
  }
  out
}

.logitgaussfns <- list(d0=.logitgauss.d0, d120=.logitgauss.d12, d340=.logitgauss.d34)

.logitgauss_unlink <- list(function(x) x, function(x) exp(x))
attr(.logitgauss_unlink[[1]], "deriv") <- function(x) 0 * x + 1
attr(.logitgauss_unlink[[2]], "deriv") <- .logitgauss_unlink[[2]]

.logitgaussfns$q <- qnorm
.logitgaussfns$unlink <- .logitgauss_unlink

