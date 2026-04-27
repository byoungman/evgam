## Negative binomial negative log-likelihood functions

.negbin.d0 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  w <- likdata$args$weights
  if (is.null(w))
    w <- 1
  if (length(w) < length(likdata$y))
    w <- 0 * likdata$y + w
  off <- likdata$args$offset
  if (is.null(off))
    off <- 1
  if (length(off) < length(likdata$y))
    off <- 0 * likdata$y + off
  nhere <- rowSums(is.finite(likdata$y))
  if (is.null(likdata$args$model))
    likdata$args$model <- 'overdispersion'
  if (likdata$args$model == 'variance') {
    if (!likdata$sparse) {
      out <- negbind0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, as.matrix(w), as.matrix(off))
    } else {
      out <- negbinspd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, as.matrix(w), as.matrix(off))
    }
  } else {
    if (!likdata$sparse) {
      out <- negbinodd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, as.matrix(w), as.matrix(off))
    } else {
      out <- negbinodspd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, as.matrix(w), as.matrix(off))
    }
  }  
  if (!is.finite(out))
    out <- 1e20
  out
}

.negbin.d12 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  w <- likdata$args$weights
  if (is.null(w))
    w <- 1
  if (length(w) < length(likdata$y))
    w <- 0 * likdata$y + w
  off <- likdata$args$offset
  if (is.null(off))
    off <- 1
  if (length(off) < length(likdata$y))
    off <- 0 * likdata$y + off
  nhere <- rowSums(is.finite(likdata$y))  
  if (is.null(likdata$args$model))
    likdata$args$model <- 'overdispersion'
  if (likdata$args$model == 'variance') {
    if (!likdata$sparse) {
      out <- negbind12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, as.matrix(w), as.matrix(off))
    } else {
      out <- negbinspd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, as.matrix(w), as.matrix(off))
    }
  } else {
    if (!likdata$sparse) {
      out <- negbinodd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, as.matrix(w), as.matrix(off))
    } else {
      out <- negbinodspd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, as.matrix(w), as.matrix(off))
    }
  }  
  out
}

.negbin.d34 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  w <- likdata$args$weights
  if (is.null(w))
    w <- 1
  if (length(w) < length(likdata$y))
    w <- 0 * likdata$y + w
  off <- likdata$args$offset
  if (is.null(off))
    off <- 1
  if (length(off) < length(likdata$y))
    off <- 0 * likdata$y + off
  nhere <- rowSums(is.finite(likdata$y))  
  if (is.null(likdata$args$model))
    likdata$args$model <- 'overdispersion'
  if (likdata$args$model == 'variance') {
    if (!likdata$sparse) {
      out <- negbind34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, as.matrix(w), as.matrix(off))
    } else {
      out <- negbinspd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, as.matrix(w), as.matrix(off))
    }
  } else {
    if (!likdata$sparse) {
      out <- negbinodd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, as.matrix(w), as.matrix(off))
    } else {
      out <- negbinodspd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, as.matrix(w), as.matrix(off))
    }
  }  
  out
}

# .negbinfns <- list(d0 = .negbin.d0, d120 = .negbin.d12, d340 = NULL)
.negbinfns <- list(d0 = .negbin.d0, d120 = .negbin.d12, d340 = .negbin.d34)

.negbin_unlink <- list(function(x) exp(x), function(x) exp(x))
attr(.negbin_unlink[[1]], "deriv") <- .negbin_unlink[[1]]
attr(.negbin_unlink[[2]], "deriv") <- .negbin_unlink[[2]]

.negbinfns$q <- qpois
.negbinfns$unlink <- .negbin_unlink

