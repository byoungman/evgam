## Left truncated gamma negative log-likelihood functions

.ltgammab.d0 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))
  left <- likdata$args$left
  if (length(left) == 1)
    left <- array(left, dim(likdata$y))
  alpha <- likdata$args$alpha
  if (length(alpha) == 1)
    alpha <- array(alpha, dim(likdata$y))
  out <- ltgammabd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$y, likdata$dupid, likdata$duplicate, nhere, left, alpha)
  if (!is.finite(out))
    out <- 1e20
  out
}

.ltgammab.d12 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))
  left <- likdata$args$left
  if (length(left) == 1)
    left <- array(left, dim(likdata$y))
  alpha <- likdata$args$alpha
  if (length(alpha) == 1)
    alpha <- array(alpha, dim(likdata$y))
  out <- ltgammabd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$y, likdata$dupid, likdata$duplicate, nhere, left, alpha)
  out
}

.ltgammab.d34 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  left <- likdata$args$left
  if (length(left) == 1)
    left <- array(left, dim(likdata$y))
  alpha <- likdata$args$alpha
  if (length(alpha) == 1)
    alpha <- array(alpha, dim(likdata$y))
  out <- ltgammabd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$y, likdata$dupid, likdata$duplicate, nhere, left, alpha)
  out
}

.ltgammabfns <- list(d0 = .ltgammab.d0, d120 = .ltgammab.d12, d340 = .ltgammab.d34)

.ltgammab_unlink <- list(function(x) exp(x))
attr(.ltgammab_unlink[[1]], "deriv") <- .ltgammab_unlink[[1]]

.ltgammafns$q <- NULL
.ltgammabfns$unlink <- .ltgammab_unlink

