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

.q_egpd1 <- function(p, pars1, pars2, pars3) {
  kappa <- exp(pars3)
  p <- p^(1 / kappa)
  .q_gpd(p, pars1, pars2)
}

.p_egpd1 <- function(x, pars1, pars2, pars3, log = FALSE) {
  out <- .p_gpd(x, pars1, pars2)
  kappa <- exp(pars3)
  out <- out^kappa
  if (log)
    out <- log(out)
  out
}

.egpd1_unlink <- list(function(x) exp(x), NULL, function(x) exp(x))
attr(.egpd1_unlink[[1]], "deriv") <- .egpd1_unlink[[1]]
attr(.egpd1_unlink[[3]], "deriv") <- .egpd1_unlink[[3]]

.egpd1fns$unlink <- .egpd1_unlink

.egpd1fns$initfn <- function(lst) {
  inits <- numeric(3)
  inits[1:2] <- c(log(mean(lst$y, na.rm = TRUE)), .05)
  inits
}

.egpd1fns$q <- .q_egpd1
.egpd1fns$p <- .p_egpd1

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

.egpd2_unlink <- list(function(x) exp(x), NULL, function(x) exp(x), function(x) exp(x), function(x) 1 / (1 + exp(-x)))
attr(.egpd2_unlink[[1]], "deriv") <- .egpd2_unlink[[1]]
attr(.egpd2_unlink[[3]], "deriv") <- .egpd2_unlink[[3]]
attr(.egpd2_unlink[[4]], "deriv") <- .egpd2_unlink[[4]]
attr(.egpd2_unlink[[5]], "deriv") <- function(x) exp(-x)/(1 + exp(-x))^2

.egpd2fns$unlink <- .egpd2_unlink

.egpd2fns$initfn <- function(lst) {
  inits <- numeric(5)
  inits[1:2] <- c(log(mean(lst$y, na.rm = TRUE)), .05)
  inits <- c(inits[1:2], -1, 1, .25)
  inits
}

.q_egpd2 <- function(p, pars1, pars2, pars3, pars4, pars5) {
  kappa1 <- exp(pars3)
  kappa2 <- exp(pars4)
  tau <- 1 / (1 + exp(-pars5))
  p <- .iG2(p, kappa1, kappa2, tau)
  .q_gpd(p, pars1, pars2)
}

.p_egpd2 <- function(x, pars1, pars2, pars3, pars4, pars5, log = FALSE) {
  out <- .p_gpd(x, pars1, pars2)
  kappa1 <- exp(pars3)
  kappa2 <- exp(pars4)
  tau <- 1 / (1 + exp(-pars5))
  out <- tau * out^kappa1 + (1 - tau) * out^kappa2
  if (log)
    out <- log(out)
  out
}

.egpd2fns$q <- .q_egpd2
.egpd2fns$p <- .p_egpd2

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

.iG3 <- function(v, delta) 1 - qbeta(1 - v, 1/delta, 2)^(1/delta)

.egpd3fns <- list(d0=.egpd3.d0, d120=.egpd3.d12, d340=.egpd3.d34, m=3, iG=.iG3)

.egpd3_unlink <- list(function(x) exp(x), NULL, function(x) exp(x))
attr(.egpd3_unlink[[1]], "deriv") <- .egpd3_unlink[[1]]
attr(.egpd3_unlink[[3]], "deriv") <- .egpd3_unlink[[3]]

.egpd3fns$unlink <- .egpd3_unlink

.egpd3fns$initfn <- function(lst) {
  inits <- numeric(3)
  inits[1:2] <- c(log(mean(lst$y, na.rm = TRUE)), .05)
  inits
}

.q_egpd3 <- function(p, pars1, pars2, pars3) {
  delta <- exp(pars3)
  p <- .iG3(p, delta)
  .q_gpd(p, pars1, pars2)
}

.p_egpd3 <- function(x, pars1, pars2, pars3, log = FALSE) {
  out <- .p_gpd(x, pars1, pars2)
  delta <- exp(pars3)
  out <- (1 - out)^delta
  out <- pbeta(out, 1/delta, 2, lower.tail = FALSE, log.p = TRUE)
  if (!log)
    out <- exp(out)
  out
}

.egpd3fns$q <- .q_egpd3
.egpd3fns$p <- .p_egpd3

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

.egpd4fns$initfn <- function(lst) {
  inits <- numeric(4)
  inits[1:2] <- c(log(mean(lst$y, na.rm = TRUE)), .05)
  inits
}

.egpd4_unlink <- list(function(x) exp(x), NULL, function(x) exp(x), function(x) exp(x))
attr(.egpd4_unlink[[1]], "deriv") <- .egpd4_unlink[[1]]
attr(.egpd4_unlink[[3]], "deriv") <- .egpd4_unlink[[3]]
attr(.egpd4_unlink[[4]], "deriv") <- .egpd4_unlink[[4]]

.egpd4fns$unlink <- .egpd4_unlink

.q_egpd4 <- function(p, pars1, pars2, pars3, pars4) {
  delta <- exp(pars3)
  kappa <- exp(pars4)
  p <- .iG4(p, delta, kappa)
  .q_gpd(p, pars1, pars2)
}

.p_egpd4 <- function(x, pars1, pars2, pars3, pars4, log = FALSE) {
  out <- .p_gpd(x, pars1, pars2)
  delta <- exp(pars3)
  kappa <- exp(pars4)
  out <- (1 - out)^delta
  out <- .5 * kappa * pbeta(out, 1 / delta, 2, lower.tail = FALSE, log.p = TRUE)
  if (!log)
    out <- exp(out)
  out
}

.egpd4fns$q <- .q_egpd4
.egpd4fns$p <- .p_egpd4
