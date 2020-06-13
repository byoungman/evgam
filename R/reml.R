## REML functions

.reml0 <- function(pars, likfns, likdata, Sdata, beta=NULL, skipfit=FALSE) {
if (is.null(beta)) beta <- attr(pars, "beta")
sp <- exp(pars)
likdata$S <- .makeS(Sdata, sp)
if (!skipfit) {
fitbeta <- .newton_step(beta, .nllh.pen, .search.pen, likdata=likdata, likfns=likfns, control=likdata$control$inner)
if (any(abs(fitbeta$gradient) > 1)) {
it0 <- likdata$control$inner$itlim
likdata$control$inner$itlim <- 10
fitbeta <- .newton_step(fitbeta$par, .nllh.pen, .search.pen, likdata=likdata, likfns=likfns, control=likdata$control$inner, newton=FALSE, alpha0=.05)
likdata$control$inner$itlim <- it0
fitbeta <- .newton_step(fitbeta$par, .nllh.pen, .search.pen, likdata=likdata, likfns=likfns, control=likdata$control$inner)
}
if (inherits(fitbeta, "try-error")) return(1e20)
} else {
fitbeta <- list(objective=.nllh.pen(beta, likdata=likdata, likfns=likfns))
fitbeta$convergence <- 0
fitbeta$gH <- .gH.pen(beta, likdata=likdata, likfns=likfns)
fitbeta$Hessian <- fitbeta$gH[[2]]
fitbeta$par <- beta
}
logdetSdata <- .logdetS(Sdata, pars)
logdetHdata <- .d0logdetH(fitbeta)
out <- fitbeta$objective + as.numeric(fitbeta$convergence != 0) * 1e20
out <- out + .5 * logdetHdata$d0 - .5 * logdetSdata$d0# - sum(dnorm(pars, 0, 1, log=TRUE))
out <- as.vector(out + likdata$const)
if (!is.finite(out)) return(1e20)
attr(out, "beta") <- fitbeta$par
attr(out, "gradient") <- fitbeta$gradient
attr(out, "Hessian") <- fitbeta$Hessian
return(out)
}

.reml1 <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL) {
if (is.null(beta)) {
    beta <- attr(pars, "beta")
} else {
    attr(pars, "beta") <- beta
}
sp <- exp(pars)
spSl <- lapply(seq_along(sp), function(i) sp[i] * attr(Sdata, "Sl")[[i]])
S <- likdata$S <- Reduce("+", spSl)
dS <- .logdetS(Sdata, pars, deriv=1)
if (is.null(H)) H <- .Hdata(.hess.pen(beta, likdata, likfns))
dH <- .d1logdetH(pars, likdata=likdata, likfns=likfns, Sdata=Sdata, H=H)
d1 <- .5 * sapply(spSl, function(x) crossprod(beta, x %*% beta))
d1 <- d1 - .5 * dS$d1
d1 <- d1 + .5 * dH$d1
d1
}

.reml1.fd <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL, kept=NULL) {
beta <- attr(pars, "beta")
eps <- 1e-4
f0 <- .reml0(pars, likfns, likdata, Sdata, beta=beta, skipfit=TRUE)
f1 <- 0 * pars
for (i in seq_along(pars)) {
  parsi <- pars
  parsi[i] <- parsi[i] + eps
  f1[i] <- .reml0(parsi, likfns, likdata, Sdata, beta=beta)
}
(f1 - f0) / eps
}

.reml12 <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL) {
if (is.null(beta)) {
    beta <- attr(pars, "beta")
} else {
    attr(pars, "beta") <- beta
}
sp <- exp(pars)
spSl <- lapply(seq_along(sp), function(i) sp[i] * attr(Sdata, "Sl")[[i]])
S <- likdata$S <- Reduce("+", spSl)
dS <- .logdetS(Sdata, pars, deriv=2)
if (is.null(H)) H <- .Hdata(.hess.pen(beta, likdata, likfns))
dH <- .d12logdetH(pars, likdata=likdata, likfns=likfns, Sdata=Sdata, H=H)
d1 <- .5 * sapply(spSl, function(x) crossprod(beta, x %*% beta))
d2 <- matrix(0, length(pars), length(pars))
diag(d2) <- d1
d2 <- d2 + crossprod(dH$dbeta$d1, S %*% dH$dbeta$d1)
d1 <- d1 - .5 * dS$d1
d2 <- d2 - .5 * dS$d2
d1 <- d1 + .5 * dH$d1
d2 <- d2 + .5 * dH$d2
d2 <- .5 * (d2 + t(d2))
list(d1=d1, d2=d2)
}

.reml2.fd <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL, kept=NULL) {
beta <- attr(pars, "beta")
eps <- 1e-4
f0 <- .reml1(pars, likfns, likdata, Sdata, beta=beta, H=H)
f1 <- matrix(0, length(pars), length(pars))
for (i in seq_along(pars)) {
  parsi <- pars
  parsi[i] <- parsi[i] + eps
  f1[i,] <- .reml1(parsi, likfns, likdata, Sdata, beta=beta, H=H)
}
f1 <- (f1 - f0) / eps
.5 * (f1 + t(f1))
}

.reml12.1 <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL) {
if (is.null(beta)) {
    beta <- attr(pars, "beta")
} else {
    attr(pars, "beta") <- beta
}
sp <- exp(pars)
spSl <- lapply(seq_along(sp), function(i) sp[i] * attr(Sdata, "Sl")[[i]])
S <- likdata$S <- Reduce("+", spSl)
dS <- .logdetS(Sdata, pars, deriv=2)
if (is.null(H)) H <- .Hdata(.hess.pen(beta, likdata, likfns))
dH <- .d12logdetH(pars, likdata=likdata, likfns=likfns, Sdata=Sdata, H=H)
d1 <- .5 * sapply(spSl, function(x) crossprod(beta, x %*% beta))
d1 <- d1 - .5 * dS$d1
d1 <- d1 + .5 * dH$d1
d1
}

.search.reml <- function(pars, likfns, likdata, Sdata, H=NULL, kept=NULL) {
gH <- .reml12(pars, likfns, likdata, Sdata, H=H)
.search.dir(gH[[1]], gH[[2]], !logical(length(gH[[1]])))
}
