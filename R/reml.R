## REML functions

.reml0 <- function(pars, likfns, likdata, Sdata, beta=NULL) {
if (is.null(beta)) 
  beta <- attr(pars, "beta")
sp <- exp(pars)
likdata$S <- .makeS(Sdata, sp)
fitbeta <- .newton_step(beta, .nllh.pen, .search.pen, likdata=likdata, likfns=likfns, control=likdata$control$inner)
if (!fitbeta$gradconv) {
  it0 <- likdata$control$inner$itlim
  likdata$control$inner$itlim <- 10
  fitbeta <- .newton_step(fitbeta$par, .nllh.pen, .search.pen, likdata=likdata, likfns=likfns, control=likdata$control$inner, newton=FALSE, alpha0=.05)
  likdata$control$inner$itlim <- it0
  fitbeta <- .newton_step(fitbeta$par, .nllh.pen, .search.pen, likdata=likdata, likfns=likfns, control=likdata$control$inner)
}
if (inherits(fitbeta, "try-error")) 
  return(1e20)
# if (!fitbeta$gradconv)
#   return(1e20)
logdetSdata <- .logdetS(Sdata, pars)
logdetHdata <- .d0logdetH(fitbeta)
out <- fitbeta$objective + as.numeric(fitbeta$convergence != 0) * 1e20
out <- out + .5 * logdetHdata$d0 - .5 * logdetSdata$d0
out <- as.vector(out + likdata$const)
if (!is.finite(out)) 
  return(1e20)
attr(fitbeta$par, "dropped") <- !fitbeta$kept
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
spSl <- Map("*", attr(Sdata, "Sl"), exp(pars))
likdata$S <- Reduce("+", spSl)
if (is.null(H)) 
  H <- .Hdata(.hess.pen(beta, likdata, likfns))
d1beta <- .d1beta(pars, beta, spSl, H)
dS <- .logdetS(Sdata, pars, deriv=1)
dH <- .d1logdetH(d1beta, likdata, likfns, spSl, H)
dbSb <- sapply(spSl, function(x) crossprod(beta, x %*% beta))
d1 <- .5 * dbSb
d1 <- d1 - .5 * dS$d1
d1 <- d1 + .5 * dH$d1
d1
}

.reml1.fd <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL, kept=NULL) {
beta <- attr(pars, "beta")
tol <- .Machine$double.eps^(1/4)
eps <- pmax(tol * abs(pars), tol)
f0 <- .reml0(pars, likfns, likdata, Sdata, beta=beta)
f1 <- 0 * pars
for (i in seq_along(pars)) {
  parsi <- pars
  parsi[i] <- parsi[i] + eps[i]
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
nsp <- length(pars)
spSl <- Map("*", attr(Sdata, "Sl"), exp(pars))
likdata$S <- Reduce("+", spSl)
if (is.null(H)) 
  H <- .Hdata(.hess.pen(beta, likdata, likfns))
dbeta <- .d1beta(pars, beta, spSl, H)
dS <- .logdetS(Sdata, pars, deriv=2)
d1H <- .d1H0(dbeta, likdata, likfns)
dbeta <- .d2beta(dbeta, d1H$d1, spSl, H)
d2H <- .d2H0_diag(dbeta, likdata, d1H$GH, H)
# first derivatives
dbSb <- crossprod(beta, dbeta$spSlb)[1,]
d1V <- -dbSb
d1V <- d1V + dS$d1
d1H$d1 <- Map("+", d1H$d1, spSl)
d1H$d1 <- lapply(d1H$d1, function(x) .precond_solve(H$cH, x))
d1V <- d1V - sapply(d1H$d1, function(x) sum(diag(x)))
# second derivatives
d2V <- diag(-dbSb, length(dbSb))
d2V <- d2V + 2 * crossprod(dbeta$d1, H$H0 %*% dbeta$d1)
d2V <- d2V + dS$d2
d2V <- d2V - d2H$d2
for (j in 1:nsp) {
  for (k in 1:j) {
    d2V[j, k] <- d2V[j, k] - sum(diag(d1H$d1[[k]] %*% d1H$d1[[j]]))
    if (j != k)
      d2V[k, j] <- d2V[j, k]
  }
}
list(d1 = -.5 * d1V, d2 = -.5 * d2V)
}

.reml12 <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL) {
if (is.null(beta)) {
  beta <- attr(pars, "beta")
} else {
  attr(pars, "beta") <- beta
}
d1 <- .reml1(pars, likfns, likdata, Sdata, H, beta)
tol <- .Machine$double.eps^(1/4)
eps <- pmax(tol * abs(pars), tol)
d2 <- matrix(NA, length(pars), length(pars))
for (i in seq_along(pars)) {
  parsi <- replace(pars, i, pars[i] + eps[i])
  d2[, i] <- (.reml1(parsi, likfns, likdata, Sdata, NULL, beta) - d1) / eps[i]
}
d2 <- .5 * (d2 + t(d2))
list(d1 = d1, d2 = d2)
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

.reml2.fdfd <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL, kept=NULL) {
beta <- attr(pars, "beta")
eps <- 1e-4
f0 <- .reml1.fd(pars, likfns, likdata, Sdata, beta=beta, H=H)
f1 <- matrix(0, length(pars), length(pars))
for (i in seq_along(pars)) {
  parsi <- pars
  parsi[i] <- parsi[i] + eps
  f1[i,] <- .reml1.fd(parsi, likfns, likdata, Sdata, beta=beta, H=H)
}
f1 <- (f1 - f0) / eps
.5 * (f1 + t(f1))
}

.search.reml <- function(pars, likfns, likdata, Sdata, H=NULL, kept=NULL) {
gH <- .reml12(pars, likfns, likdata, Sdata, H=H)
.search.dir(gH[[1]], gH[[2]], !logical(length(gH[[1]])))
}

## fixed smoothing parameters

.reml0_fixed <- function(pars, likfns, likdata, Sdata, beta=NULL) {
if (is.null(beta)) 
  beta <- attr(pars, "beta")
sp <- exp(pars)
likdata$S <- .makeS(Sdata, sp)
fitbeta <- .newton_step(beta, .nllh.pen, .search.pen, likdata=likdata, likfns=likfns, control=likdata$control$inner)
if (!fitbeta$gradconv) {
  it0 <- likdata$control$inner$itlim
  likdata$control$inner$itlim <- 10
  fitbeta <- .newton_step(fitbeta$par, .nllh.pen, .search.pen, likdata=likdata, likfns=likfns, control=likdata$control$inner, newton=FALSE, alpha0=.05)
  likdata$control$inner$itlim <- it0
  fitbeta <- .newton_step(fitbeta$par, .nllh.pen, .search.pen, likdata=likdata, likfns=likfns, control=likdata$control$inner)
}
list(beta = fitbeta$par)
}
