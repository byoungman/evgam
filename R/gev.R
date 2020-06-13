## Generalised extreme value negative log-likelihood functions

.gev.d0 <- function(pars, likdata) {
if (!likdata$censored) {
  out <- gevd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
} else {
  id <- likdata$cens.id
  out1 <- gevd0(split(pars, likdata$idpars), likdata$X[[1]][!id, , drop=FALSE], likdata$X[[2]][!id, , drop=FALSE], likdata$X[[3]][!id, , drop=FALSE], likdata$y[!id, 1], likdata$dupid[!id], likdata$duplicate)
  out2 <- gevcd0(split(pars, likdata$idpars), likdata$X[[1]][id, , drop=FALSE], likdata$X[[2]][id, , drop=FALSE], likdata$X[[3]][id, , drop=FALSE], likdata$y[id, , drop=FALSE], likdata$dupid[id], likdata$duplicate)
  out <- out1 + out2
}
out
}

.gev.d12 <- function(pars, likdata) {
if (!likdata$censored) {
  out <- gevd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
} else {
  id <- likdata$cens.id
  out <- matrix(0, likdata$nobs, 9)
  out[!id, ] <- gevd12(split(pars, likdata$idpars), likdata$X[[1]][!id, , drop=FALSE], likdata$X[[2]][!id, , drop=FALSE], likdata$X[[3]][!id, , drop=FALSE], likdata$y[!id, 1], likdata$dupid[!id], likdata$duplicate)
  out[id, ] <- gevcd12(split(pars, likdata$idpars), likdata$X[[1]][id, , drop=FALSE], likdata$X[[2]][id, , drop=FALSE], likdata$X[[3]][id, , drop=FALSE], likdata$y[id, , drop=FALSE], likdata$dupid[id], likdata$duplicate)
}
out
}

.gev.d34 <- function(pars, likdata) {
if (!likdata$censored) {
  out <- gevd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
} else {
  id <- likdata$cens.id
  out <- matrix(0, likdata$nobs, 25)
  out[!id, ] <- gevd34(split(pars, likdata$idpars), likdata$X[[1]][!id, , drop=FALSE], likdata$X[[2]][!id, , drop=FALSE], likdata$X[[3]][!id, , drop=FALSE], likdata$y[!id, 1], likdata$dupid[!id], likdata$duplicate)
  out[id, ] <- gevcd34(split(pars, likdata$idpars), likdata$X[[1]][id, , drop=FALSE], likdata$X[[2]][id, , drop=FALSE], likdata$X[[3]][id, , drop=FALSE], likdata$y[id, , drop=FALSE], likdata$dupid[id], likdata$duplicate)
}
out
}

.gevfns <- list(d0=.gev.d0, d120=.gev.d12, d340=.gev.d34)

.pgev <- function(x, loc, scale, shape, NAOK=FALSE, log=FALSE) {
# function to evaluate daily cdf, e.g. Coles (2001, pp.138)
temp <- 1 + shape * (x - loc) / scale
if (!NAOK) temp <- pmax(temp, 0)
out <- -(temp ^ (-1/shape))
if (!log) out <- exp(out)
out
}

.qgev <- function(p, loc, scale, shape) {
shape[shape == 0] <- 1e-6
shape <- sign(shape) * pmax(abs(shape), 1e-6)
yp <- -log(p)
loc - scale * (1 - yp^(-shape)) / shape
}

.dqgev <- function(p, loc, lscale, shape) {
.e1 <- -log(p)
.e2 <- .e1^shape
.e3 <- 1 - 1/.e2
.e4 <- exp(lscale)
d1 <- 1
d2 <- -(.e3 * .e4/shape)
d3 <- -(.e4 * (log(.e1)/.e2 - .e3/shape)/shape)
cbind(d1, d2, d3)
}
