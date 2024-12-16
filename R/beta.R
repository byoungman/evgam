## Beta negative log-likelihood functions

.beta.d0 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  if (!likdata$sparse) {
    out <- betad0(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], ymat=likdata$y, dupid=likdata$dupid, dcate=likdata$duplicate, nhere = nhere)
  } else {
    out <- betaspd0(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], ymat=likdata$y, dupid=likdata$dupid, dcate=likdata$duplicate, nhere = nhere)
  }
  if (!is.finite(out)) out <- 1e20
  out
}

.beta.d12 <- function(pars, likdata, sandwich = FALSE) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  if (!likdata$sparse) {
    out <- betad12(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], ymat=likdata$y, dupid=likdata$dupid, dcate=likdata$duplicate, nhere = nhere)
  } else {
    out <- betaspd12(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], ymat=likdata$y, dupid=likdata$dupid, dcate=likdata$duplicate, nhere = nhere)
  }
  out
}

.beta.d34 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  if (!likdata$sparse) {
    out <- betad34(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], ymat=likdata$y, dupid=likdata$dupid, dcate=likdata$duplicate, nhere = nhere)
  } else {
    out <- betaspd34(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], ymat=likdata$y, dupid=likdata$dupid, dcate=likdata$duplicate, nhere = nhere)
  }
  out
}

.betafns <- list(d0=.beta.d0, d120=.beta.d12, d340=.beta.d34)

.beta_unlink <- list(function(x) exp(x), function(x) exp(x))
attr(.beta_unlink[[1]], "deriv") <- .beta_unlink[[1]]
attr(.beta_unlink[[2]], "deriv") <- .beta_unlink[[2]]

.betafns$q <- qbeta
.betafns$unlink <- .beta_unlink
