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

.gevfns <- list(d0 = .gev.d0, d120 = .gev.d12, d340 = .gev.d34)

.pgev <- function(x, loc, scale, shape, NAOK=FALSE, log=FALSE) {
# function to evaluate daily cdf, e.g. Coles (2001, pp.138)
temp <- 1 + shape * (x - loc) / scale
if (!NAOK) temp <- pmax(temp, 0)
out <- -(temp ^ (-1/shape))
if (!log) out <- exp(out)
out
}

.gev_unlink <- list(NULL, function(x) exp(x), NULL)
attr(.gev_unlink[[2]], "deriv") <- .gev_unlink[[2]]

.q_gev <- function(p, pars1, pars2, pars3) {
  loc <- pars1
  scale <- exp(pars2)
  shape <- pars3
  # shape <- sign(shape) * pmax(abs(shape), 1e-6)
  yp <- -log(p)
  loc - scale * (1 - yp^(-shape)) / shape
}

# Deriv::Deriv(.q_gev, paste('pars', 1:3, sep = ''), combine = 'cbind')

.dq_gev <- function (p, pars1, pars2, pars3) {
  .e1 <- -log(p)
  .e2 <- .e1^pars3
  .e3 <- 1 - 1/.e2
  .e4 <- exp(pars2)
  cbind(pars1 = 1, 
        pars2 = -(.e3 * .e4/pars3), 
        pars3 = -(.e4 * (log(.e1)/.e2 - .e3/pars3)/pars3)
        )
}

.gevfns$q <- .q_gev
.gevfns$dq <- .dq_gev
.gevfns$unlink <- .gev_unlink

.gevfns$initfn <- function(lst) {
  inits <- sqrt(6) * sd(as.vector(lst$y), na.rm = TRUE) / pi
  c(mean(as.matrix(lst$y), na.rm = TRUE) - .5772 * inits[1], log(inits[1]), .05)
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

.p_gev <- function(x, pars1, pars2, pars3, log = FALSE) {
  loc <- pars1
  scale <- exp(pars2)
  shape <- pars3
  temp <- 1 + shape * (x - loc) / scale
  temp <- pmax(temp, 0)
  out <- -(temp ^ (-1/shape))
  if (!log) 
    out <- exp(out)
  out
}

.gevfns$p <- .p_gev