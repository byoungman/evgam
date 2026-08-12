## REML functions

.pars2rho <- function(pars, likdata, Sdata, type = 'rho') {
  pars <- as.vector(pars)
  if (!is.null(likdata$replacer)) {
    rho0 <- likdata$replacer$start
    rho0[likdata$replacer$rho_change] <- pars
  } else {
    rho0 <- pars
  }
  rho0 <- as.vector(rho0)
  if (type == 'rho0')
    return(rho0)
  rho <- as.vector(rho0 %*% attr(Sdata, 'A'))
  if (!is.null(likdata$replacer))
    rho[likdata$replacer$fixed] <- likdata$replacer$rho_fixed
  rho
}

.reml0 <- function(pars, likfns, likdata, Sdata, beta=NULL) {
  if (is.null(beta)) 
    beta <- attr(pars, "beta")
  rho <- .pars2rho(pars, likdata, Sdata)
  sp <- exp(rho)
  likdata$S <- .makeS(Sdata, sp)
  fitbeta <- .newton_step(beta, .nllh.pen, .search.pen, likdata=likdata, likfns=likfns, control=likdata$control$inner)
  if (!fitbeta$gradconv) {
    if (!likdata$sparse)
      fitbeta <- nlminb(fitbeta$par, .nllh.pen, .grad.pen, .hess.pen, likdata = likdata, likfns = likfns)
    fitbeta <- .newton_step(fitbeta$par, .nllh.pen, .search.pen, likdata=likdata, likfns=likfns, control=likdata$control$inner)
  }
  if (inherits(fitbeta, "try-error")) 
    return(1e20)
  logdetSdata <- try(.logdetS(Sdata, rho), silent = TRUE)
  if (inherits(logdetSdata, 'try-error'))
    logdetSdata <- list(d0 = 1e20)
  logdetHdata <- .d0logdetH(fitbeta, likdata$sparse)
  out <- fitbeta$objective + as.numeric(fitbeta$convergence != 0) * 1e20
  out <- out + .5 * logdetHdata$d0 - .5 * logdetSdata$d0
  out <- as.vector(out + likdata$const)
  if (!is.finite(out)) 
    return(1e20)
  attr(fitbeta$par, "dropped") <- !fitbeta$kept
  attr(out, "beta") <- fitbeta$par
  attr(out, "gradient") <- fitbeta$gradient
  attr(out, "Hessian") <- fitbeta$Hessian
  attr(out, "objective") <- fitbeta$objective
  attr(out, "pars") <- as.vector(pars)
  return(out)
}

.reml1 <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL) {
  if (is.null(beta)) {
    beta <- attr(pars, "beta")
  } else {
    attr(pars, "beta") <- beta
  }
  rho <- .pars2rho(pars, likdata, Sdata)
  sp <- exp(rho)
  spSl <- Map("*", attr(Sdata, "Sl"), sp)
  likdata$S <- Reduce("+", spSl)
  if (is.null(H)) 
    H <- .Hdata(.hess.pen(beta, likdata, likfns))
  d1beta <- .d1beta(rho, beta, spSl, H)
  dS <- .logdetS(Sdata, rho, deriv=1)
  dH <- .d1logdetH(d1beta, likdata, likfns, spSl, H)
  if (!likdata$sparse) {
    dbSb <- sapply(spSl, function(x) base::crossprod(beta, x %*% beta)[1, 1])
  } else {
    dbSb <- sapply(spSl, function(x) Matrix::crossprod(beta, x %*% beta)[1, 1])
  }
  d1 <- .5 * dbSb
  d1 <- d1 - .5 * dS$d1
  d1 <- d1 + .5 * dH$d1
  d1 <- attr(Sdata, 'A') %*% d1
  if (likdata$some_sp_fixed)
    d1 <- d1[likdata$replacer$rho_change]
  d1
}

# .reml1.fd <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL, kept=NULL) {
# beta <- attr(pars, "beta")
# tol <- .Machine$double.eps^(1/4)
# eps <- pmax(tol * abs(pars), tol)
# f0 <- .reml0(pars, likfns, likdata, Sdata, beta=beta)
# f1 <- 0 * pars
# for (i in seq_along(pars)) {
#   parsi <- pars
#   parsi[i] <- parsi[i] + eps[i]
#   f1[i] <- .reml0(parsi, likfns, likdata, Sdata, beta=beta)
# }
# (f1 - f0) / eps
# }

.reml1.fd <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL, kept=NULL) {
  beta <- attr(pars, "beta")
  tol <- .Machine$double.eps^(1/4)
  eps <- pmax(tol * abs(pars), tol)
  fl <- fh <- 0 * pars
  for (i in seq_along(pars)) {
    parsi <- replace(pars, i, pars[i] + eps[i])
    fh[i] <- .reml0(parsi, likfns, likdata, Sdata, beta = beta)
  }
  for (i in seq_along(pars)) {
    parsi <- replace(pars, i, pars[i] - eps[i])
    fl[i] <- .reml0(parsi, likfns, likdata, Sdata, beta = beta)
  }
  .5 * (fh - fl) / eps
}

# .reml12 <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL) {
# if (is.null(beta)) {
#   beta <- attr(pars, "beta")
# } else {
#   attr(pars, "beta") <- beta
# }
# nsp <- length(pars)
# spSl <- Map("*", attr(Sdata, "Sl"), exp(pars))
# likdata$S <- Reduce("+", spSl)
# if (is.null(H)) 
#   H <- .Hdata(.hess.pen(beta, likdata, likfns))
# dbeta <- .d1beta(pars, beta, spSl, H)
# dS <- .logdetS(Sdata, pars, deriv=2)
# d1H <- .d1H0(dbeta, likdata, likfns)
# dbeta <- .d2beta(dbeta, d1H$d1, spSl, H)
# d2H <- .d2H0_diag(dbeta, likdata, d1H$GH, H)
# # first derivatives
# dbSb <- crossprod(beta, dbeta$spSlb)[1,]
# d1V <- -dbSb
# d1V <- d1V + dS$d1
# d1H$d1 <- Map("+", d1H$d1, spSl)
# d1H$d1 <- lapply(d1H$d1, function(x) .precond_solve(H$cH, x))
# d1V <- d1V - sapply(d1H$d1, function(x) sum(diag(x)))
# # second derivatives
# d2V <- diag(-dbSb, length(dbSb))
# d2V <- d2V + 2 * crossprod(dbeta$d1, H$H0 %*% dbeta$d1)
# d2V <- d2V + dS$d2
# d2V <- d2V - d2H$d2
# for (j in 1:nsp) {
#   for (k in 1:j) {
#     d2V[j, k] <- d2V[j, k] - sum(diag(d1H$d1[[k]] %*% d1H$d1[[j]]))
#     if (j != k)
#       d2V[k, j] <- d2V[j, k]
#   }
# }
# list(d1 = -.5 * d1V, d2 = -.5 * d2V)
# }

.reml2 <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL) {
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
dS <- .logdetS(Sdata, pars, deriv = 2)
d1H <- .d1H0(dbeta, likdata, likfns)
dbeta <- .d2beta(dbeta, d1H$d1, spSl, H)
d2H <- .d2logdetH(dbeta, likdata, likfns, spSl, H)
# d2H <- .d2H0_diag(dbeta, likdata, d1H$GH, H)
# first derivatives
dbSb <- crossprod(beta, dbeta$spSlb)[1,]
# d1V <- -dbSb
# d1V <- d1V + dS$d1
# d1H$d1 <- Map("+", d1H$d1, spSl)
# d1H$d1 <- lapply(d1H$d1, function(x) .precond_solve(H$cH, x))
# d1V <- d1V - sapply(d1H$d1, function(x) sum(diag(x)))
# second derivatives

d2V <- crossprod(dbeta$d1, H$H0 %*% dbeta$d1)
diag(d2V) <- diag(d2V) - dbSb
d2V <- d2V + .5 * dS$d2
d2V <- d2V + .5 * d2H$d2
# for (j in 1:nsp) {
#   for (k in 1:j) {
#     d2V[j, k] <- d2V[j, k] - sum(diag(d1H$d1[[k]] %*% d1H$d1[[j]]))
#     if (j != k)
#       d2V[k, j] <- d2V[j, k]
#   }
# }
# list(d1 = -.5 * d1V, d2 = -.5 * d2V)
list(d2 = -d2V)
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
# for (i in seq_along(pars)) {
#   parsi <- replace(pars, i, pars[i] + eps[i])
#   d2[, i] <- (.reml1(parsi, likfns, likdata, Sdata, NULL, beta) - d1) / eps[i]
# }
for (i in seq_along(pars)) {
  phi <- replace(pars, i, pars[i] + eps[i])
  dhi <- .reml1(phi, likfns, likdata, Sdata, NULL, beta)
  plo <- replace(pars, i, pars[i] - eps[i])
  dlo <- .reml1(plo, likfns, likdata, Sdata, NULL, beta)
  d2[, i] <- .5 * (dhi - dlo) / eps[i]
}
d2 <- .5 * (d2 + t(d2))
list(d1 = d1, d2 = d2)
}

.reml12 <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL) {
  list(d1 = .reml1(pars, likfns, likdata, Sdata, H, beta),
       d2 = .reml2(pars, likfns, likdata, Sdata, H, beta)$d2)
}

.reml12 <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL) {
  if (is.null(beta)) {
    beta <- attr(pars, "beta")
  } else {
    attr(pars, "beta") <- beta
  }
  rho <- .pars2rho(pars, likdata, Sdata)
  sp <- exp(rho)
  spSl <- Map("*", attr(Sdata, "Sl"), sp)
  likdata$S <- Reduce("+", spSl)
  if (is.null(H)) 
    H <- .Hdata(.hess.pen(beta, likdata, likfns))
  d1beta <- .d1beta(rho, beta, spSl, H)
  dS <- .logdetS(Sdata, rho, deriv = 2)
  d1H <- .d1H0(d1beta, likdata, likfns)
  d12beta <- .d2beta(d1beta, d1H$d1, spSl, H)
  d12H <- .d12logdetH(d12beta, likdata, likfns, spSl, H)
  if (!likdata$sparse) {
    dbSb <- sapply(spSl, function(x) base::crossprod(beta, x %*% beta)[1, 1])
  } else {
    dbSb <- sapply(spSl, function(x) Matrix::crossprod(beta, x %*% beta)[1, 1])
  }
  dV1 <- -.5 * dbSb
  dV1 <- dV1 + .5 * dS$d1
  dV1 <- dV1 - .5 * d12H$d1
  dV1 <- attr(Sdata, 'A') %*% dV1
  if (likdata$some_sp_fixed)
    dV1 <- dV1[likdata$replacer$rho_change]
  dV2 <- crossprod(d12beta$d1, H$H0 %*% d12beta$d1)
  diag(dV2) <- diag(dV2) - dbSb
  dV2 <- dV2 + .5 * dS$d2
  dV2 <- dV2 - .5 * d12H$d2
  list(d1 = -dV1, d2 = -dV2)
}

.reml12_diag <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL) {
  if (is.null(beta)) {
    beta <- attr(pars, "beta")
  } else {
    attr(pars, "beta") <- beta
  }
  rho <- .pars2rho(pars, likdata, Sdata)
  sp <- exp(rho)
  spSl <- Map("*", attr(Sdata, "Sl"), sp)
  likdata$S <- Reduce("+", spSl)
  if (is.null(H)) 
    H <- .Hdata(.hess.pen(beta, likdata, likfns))
  d1beta <- .d1beta(rho, beta, spSl, H)
  dS <- .logdetS(Sdata, rho, deriv = 2)
  d1H <- .d1H0(d1beta, likdata, likfns)
  d12beta <- .d2beta(d1beta, d1H$d1, spSl, H)
  d12H <- .d12logdetH_diag(d12beta, likdata, likfns, spSl, H)
  if (!likdata$sparse) {
    dbSb <- sapply(spSl, function(x) base::crossprod(beta, x %*% beta)[1, 1])
  } else {
    dbSb <- sapply(spSl, function(x) Matrix::crossprod(beta, x %*% beta)[1, 1])
  }
  dV1 <- -.5 * dbSb
  dV1 <- dV1 + .5 * dS$d1
  dV1 <- dV1 - .5 * d12H$d1
  dV1 <- attr(Sdata, 'A') %*% dV1
  if (likdata$some_sp_fixed)
    dV1 <- dV1[likdata$replacer$rho_change]
  dV2 <- colSums(d12beta$d1 * (H$H0 %*% d12beta$d1))
  dV2 <- dV2 - dbSb
  dV2 <- dV2 + .5 * diag(dS$d2)
  dV2 <- dV2 - .5 * d12H$d2
  list(d1 = -dV1, d2 = -dV2)
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

.search.reml_diag <- function(pars, likfns, likdata, Sdata, H=NULL, kept=NULL) {
  gH <- .reml12_diag(pars, likfns, likdata, Sdata, H=H)
  out <- gH[[1]] / pmax(gH[[2]], 1e-1)
  attr(out, 'gradient') <- gH[[1]]
  out
}
