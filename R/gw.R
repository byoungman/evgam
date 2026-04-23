## GW negative log-likelihood functions

.gw.d0 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  w <- likdata$args$weights
  if (is.null(w))
    w <- 1
  if (length(w) < length(likdata$y))
    w <- 0 * likdata$y + w
  p <- likdata$args$p
  if (length(p) == 1)
    p <- 0 * likdata$y + p
  p <- as.matrix(p)
  nhere <- rowSums(is.finite(likdata$y))
  out <- gwd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, w, p)
  # if (!likdata$sparse) {
  #   out <- gwd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, likdata$w)
  # } else {
  #   out <- gwspd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, likdata$w)
  # }
  if (!is.finite(out))
    out <- 1e20
  out
}

.gw.d12 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  w <- likdata$args$weights
  if (is.null(w))
    w <- 1
  if (length(w) < length(likdata$y))
    w <- 0 * likdata$y + w
  p <- likdata$args$p
  if (length(p) == 1)
    p <- 0 * likdata$y + p
  p <- as.matrix(p)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- gwd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, w, p)
  # if (!likdata$sparse) {
  #   out <- gwd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, likdata$w)
  # } else {
  #   out <- gwspd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, likdata$w)
  # }
  out
}

.gw.d34 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  w <- likdata$args$weights
  if (is.null(w))
    w <- 1
  if (length(w) < length(likdata$y))
    w <- 0 * likdata$y + w
  p <- likdata$args$p
  if (length(p) == 1)
    p <- 0 * likdata$y + p
  nhere <- rowSums(is.finite(likdata$y))
  p <- as.matrix(p)
  out <- gwd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, w, p)
  # if (!likdata$sparse) {
  #   out <- gwd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, likdata$w)
  # } else {
  #   out <- gwspd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, likdata$w)
  # }
  out
}

.gwfns <- list(d0 = .gw.d0, d120 = .gw.d12, d340 = .gw.d34)

.gw_unlink <- list(function(x) exp(x), function(x) x)
attr(.gw_unlink[[1]], "deriv") <- .gw_unlink[[1]]
attr(.gw_unlink[[2]], "deriv") <- function(x) 0 * x + 1

.gwfns$q <- qpois
.gwfns$unlink <- .gw_unlink
