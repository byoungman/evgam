## Poisson negative log-likelihood functions

.pois.d0 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))
  if (!likdata$sparse) {
    out <- poisd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  } else {
    out <- poisspd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  }
  out
}

.pois.d12 <- function(pars, likdata, sandwich = FALSE) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))
  if (!likdata$sparse) {
    out <- poisd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  } else {
    out <- poisspd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  }
  out
}

.pois.d34 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))
  if (!likdata$sparse) {
    out <- poisd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  } else {
    out <- poisspd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  }
  out
}

.poisfns <- list(d0 = .pois.d0, d120 = .pois.d12, d340 = .pois.d34)

.pois_unlink <- list(function(x) exp(x))
attr(.pois_unlink[[1]], "deriv") <- .pois_unlink[[1]]

.poisfns$q <- qpois
.poisfns$unlink <- .pois_unlink

