## extended generalised Pareto negative log-likelihood functions

## model 1 ##

.egpd1.d0 <- function(pars, likdata) {
if (likdata$censored)
  stop("Censored likelihoods not currently available for extended GPDs.")
egpd1d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.egpd1.d12 <- function(pars, likdata) {
egpd1d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.egpd1.d34 <- function(pars, likdata) {
egpd1d34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.iG1 <- function(v, kappa) v^(1/kappa)

.egpd1fns <- list(d0=.egpd1.d0, d120=.egpd1.d12, d340=.egpd1.d34, m=1, iG=.iG1)

## model 2 ##

.G2 <- function(v, kappa1, kappa2, p) {
p * v^kappa1 + (1 - p) * v^kappa2
}

.iG2i <- function(v, kappa1, kappa2, p) {
vv <- range(c(v^1/kappa1, v^1/kappa2))
d <- diff(vv)
lo <- vv[1]
while(.G2(lo, kappa1, kappa2, p) - v > 0) lo <- max(0, lo - d)
hi <- vv[2]
while(.G2(hi, kappa1, kappa2, p) - v > 0) hi <- min(1, hi + d)
uniroot(function(x) .G2(x, kappa1, kappa2, p) - v, c(lo, hi))$root
}

.iG2 <- function(v, kappa1, kappa2, p) {
n <- length(v)
vapply(seq_len(n), function(i) .iG2i(v[i], kappa1[i], kappa2[i], p[i]), double(1))
}

.egpd2.d0 <- function(pars, likdata) {
if (likdata$censored)
  stop("Censored likelihoods not currently available for extended GPDs.")
egpd2d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.egpd2.d12 <- function(pars, likdata) {
egpd2d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.egpd2.d34 <- function(pars, likdata) {
egpd2d34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.egpd2fns <- list(d0=.egpd2.d0, d120=.egpd2.d12, d340=NULL, m=2, iG=.iG2)

## model 3 ##

.egpd3.d0 <- function(pars, likdata) {
if (likdata$censored)
  stop("Censored likelihoods not currently available for extended GPDs.")
egpd3d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.egpd3.d12 <- function(pars, likdata) {
egpd3d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.egpd3.d34 <- function(pars, likdata) {
egpd3d34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.iG3 <- function(v, delta) 1 - qbeta(1 -v, 1/delta, 2)^(1/delta)

.egpd3fns <- list(d0=.egpd3.d0, d120=.egpd3.d12, d340=.egpd3.d34, m=3, iG=.iG3)

## model 4 ##

.egpd4.d0 <- function(pars, likdata) {
if (likdata$censored)
  stop("Censored likelihoods not currently available for extended GPDs.")
egpd4d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.egpd4.d12 <- function(pars, likdata) {
egpd4d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.egpd4.d34 <- function(pars, likdata) {
egpd4d34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$y[,1], likdata$dupid, likdata$duplicate)
}

.iG4 <- function(v, delta, kappa) 1 - qbeta(1 - v^(2/kappa), 1/delta, 2)^(1/delta)

.egpd4fns <- list(d0=.egpd4.d0, d120=.egpd4.d12, d340=NULL, m=4, iG=.iG4)

