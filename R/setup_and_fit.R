## Setup and fitting functions
#
# .setup.family
# --- identifies family, its # parameters, their names and its functions
#
# .setup.formula
# --- takes evgam formula and returns something useful for fitting
#
# .predictable.gam
# --- takes a gam(..., fit=FALSE) object and adds stuff so that mgcv::predict.gam works
#
# .X.evgam
# --- gives design matrics for a evgam object
#
# .setup.data
# --- takes data and evgam argument and return something useful for fitting
#
# .setup.pp.data
# --- takes data and converts to something useful for point process model with quadrature points
#
# .sandwich.C
# --- forms C matrix from curvature adjustment as in Chandler & Bate (2007)
#
# .setup.inner.inits
# --- gives initial basis coefficients, having first fitted constant model
#
# .guess
# --- guesses initial log smoothing parameters
#
# .sandwich
# --- updates fitting data if sandwich correction is to be used
#
# .outer
# --- performs outer iterations, i.e. smoothin parameter estimation
#
# .VpVc
# --- computes variance-covariance matrices for basis coefficients
#
# .edf
# --- computes effective degrees of freedom
#
# .swap
# --- swaps elements in early GAM objects, created for fitting, with final estimates
#
# .finalise
# --- adds other useful things to final objects for user functions
#

############ .setup.family ##########################
  
.setup.family <- function(family, pp, egpd, formula, likfns) {
if (family == "gev") {
  lik.fns <- .gevfns
  npar <- 3
  nms <- c("mu", "lpsi", "xi")
} else {
if (family == "gpd") {
  lik.fns <- .gpdfns
  npar <- 2
  nms <- c("lpsi", "xi")
} else {
if (family == "modgpd") {
stop("'family='modgpd'' will return; in the mean time use `family='gpd''")
# gone, but not forgotten
  lik.fns <- NULL
  npar <- 2
  nms <- c("lmodpsi", "xi")
} else {
if (family == "pp") {
  lik.fns <- .ppfns
  npar <- 3
  nms <- c("mu", "lpsi", "xi")
} else {
if (family == "weibull") {
  lik.fns <- .weibfns
  npar <- 2
  nms <- c("llambda", "lk")
} else {
if (family == "exi") {
  lik.fns <- .exifns
  npar <- 1
  nms <- c("location")
} else { 
if (family == "ald") {
  lik.fns <- .aldfns
  npar <- 2
  nms <- c("mu", "lsigma")
} else {
if (family == "gamma") {
  lik.fns <- NULL#gammafns
  npar <- 2
  nms <- c("ltheta", "lk")
} else {
if (family == "orthoggpd") {
stop("'family='orthoggpd'' may not return return")
# gone, and possibly forgotten
  lik.fns <- NULL#ogpdfns
  npar <- 2
  nms <- c("lnu", "xi")
} else {
if (family == "transxigpd") {
stop("'family='transxigpd'' may not return")
# gone, and possibly forgotten
  lik.fns <- NULL#txigpdfns
  npar <- 2
  nms <- c("lpsi", "xi")
} else {
if (family == "transgev") {
stop("'family='transgev'' may not return")
# gone, and possibly forgotten
  lik.fns <- NULL#transgevfns
  npar <- 6
  nms <- c("mu", "lpsi", "xi", "A", "lB", "C")
} else {
if (family == "exponential") {
  lik.fns <- .expfns
  npar <- 1
  nms <- c("llambda")
} else {
if (family == "gauss") {
  lik.fns <- .gaussfns
  npar <- 2
  nms <- c("mu", "logsigma")
} else {
if (family == "egpd") {
  if (is.null(egpd$m))
    egpd$m <- 1
  if (egpd$m == 1) {
    lik.fns <- .egpdfns1
    npar <- 3
    nms <- c("lpsi", "xi", "lkappa")
  } else {
    if (egpd$m == 2) {
      lik.fns <- .egpdfns2
      npar <- 5
      nms <- c("lpsi", "xi", "lkappa1", "lkappa2", "logitp")
      stop("Extended GPD with 'egpd.arg$m == 2' currently not available.")
    } else {
      if (egpd$m == 3) {
        lik.fns <- .egpdfns3
        npar <- 3
        nms <- c("lpsi", "xi", "ldelta")
        stop("Extended GPD with 'egpd.arg$m == 3' currently not available.")
      } else {
        lik.fns <- .egpdfns4
        npar <- 4
        nms <- c("lpsi", "xi", "lkappa", "ldelta")
        stop("Extended GPD with 'egpd.arg$m == 4' currently not available.")
      }
    }
  }
} else {
  if (length(likfns)) {
    lik.fns <- likfns
    family <- "custom"
    npar <- length(formula)
  }
}
}
}
}
}
}
}
}
}
}
}
}
}
}
out <- list(npar=npar, npar2=npar, lik.fns=lik.fns, nms=nms)
}

############ .setup.formula ##########################

.setup.formulae <- function(formula, npar, npar2, data, trace) {
# turn formula into list, which will be repeated.
if (inherits(formula, "formula"))
  formula <- list(formula)
# get variable names
if (npar == 1) {
  if (!(length(formula) %in% c(npar, 1)))
    stop("length(formula) for this family should be 1")
} else {
  if (!(length(formula) %in% c(npar, 1)))
    stop(paste("length(formula) for this family should be", npar, "(or 1 if all parameters are to have the same formula)"))
}
pred.vars <- unique(unlist(lapply(formula, all.vars)))
# check they're all in data
if (!all(pred.vars %in% names(data))) {
  missing.vars <- pred.vars[!(pred.vars %in% names(data))]
  stop(paste("Variable(s) '", paste(missing.vars, collapse=", "), "' not supplied to `data'.", sep=""))
}
# check if first element of list has response; otherwise get response name
terms.list <- lapply(formula, terms.formula, specials=c("s", "te", "ti"))
got.specials <- sapply(lapply(terms.list, function(x) unlist(attr(x, "specials"))), any)
termlabels.list <- lapply(terms.list, attr, "term.labels")
got.intercept <- sapply(terms.list, attr, "intercept") == 1
for (i in seq_along(termlabels.list)) {
  if (length(termlabels.list[[i]]) == 0) {
    if (got.intercept[i]) {
      termlabels.list[[i]] <- "1"
    } else {
      stop(paste("formula element", i, "incorrectly specified"))
    }
  }
}
got.response <- sapply(terms.list, attr, "response") == 1
if (!got.response[1]) {
  stop("formula has no response")
} else {
  response.name <- as.character(formula[[1]])[2]
}
if (any(!got.response)) {
  for (i in which(!got.response)) {
    formula[[i]] <- reformulate(termlabels=termlabels.list[[i]], response=response.name)
  }
}
stripped.formula <- lapply(termlabels.list, function(x) reformulate(termlabels=x))
censored <- FALSE
# now allow for the possibility of a [a, b] response
if (substr(response.name, 1, 5) == "cens(") {
  response.name <- substr(response.name, 6, nchar(response.name) - 1)
  response.name <- gsub(" ", "", response.name)
  response.name <- strsplit(response.name, ",")[[1]]
  if (length(response.name) > 2)
    stop("Censored response can only contain two variables.")
  rr <- response.name[2]
  formula <- lapply(termlabels.list, function(x) reformulate(termlabels=x, response=rr))
  censored <- TRUE
}
#
attr(formula, "response.name") <- response.name
pred.vars <- pred.vars[!(pred.vars %in% response.name)]
attr(formula, "predictor.names") <- pred.vars
attr(formula, "stripped") <- stripped.formula
attr(formula, "censored") <- censored
attr(formula, "smooths") <- got.specials
for (i in seq_along(formula)) {
  attr(formula[[i]], "intercept") <- got.intercept[i]
  attr(formula[[i]], "smooth") <- got.specials[i]
}
formula
}

############ .predictable.gam ##########################

.predictable.gam <- function(G, formula) {
keep <- c("dev.extra", "pterms", "nsdf", "X", "terms", "mf", "smooth", "sp")
G <- G[keep]
G$nb <- ncol(G$X)
G$coefficients <- numeric(G$nb)
old <- c("mf", "pP", "cl")
new <- c("model", "paraPen", "call")
is.in <- !is.na(match(old, names(G)))
if (any(is.in)) names(G)[match(old[is.in], names(G))] <- new[is.in]
G$formula <- formula
class(G) <- "gamlist"
G
}

############ .X.evgam ##########################

.X.evgam <- function(object, newdata) {
object <- object[sapply(object, inherits, what="gamlist")]
if (missing(newdata)) {
  X <- lapply(object, function(x) x$X)
} else {
  for (i in seq_along(object)) 
    class(object[[i]]) <- "gam"
  X <- lapply(object, mgcv::predict.gam, newdata=newdata, type="lpmatrix")
}
names(X) <- names(object)
X
}

############ .setup.data ##########################

.setup.data <- function(data, responsename, formula, family, nms, removeData, 
exiargs, aldargs, pp, knots, maxdata, maxspline, compact, sargs, 
outer, trace) {

## data
for (i in seq_along(responsename)) data <- data[!is.na(data[,responsename[i]]),]

if (nrow(data) > maxdata) {
    id <- sort(sample(nrow(data), maxdata))
    data <- data[id,]
    if (trace >= 0)
      message("`data' truncated to `maxdata' rows. Re-supply `data' to, e.g., `predict.evgam'")
}

if  (compact) {
data.undup <- as.list(data[,unique(unlist(lapply(formula, function(y) unlist(lapply(mgcv::interpret.gam(y)$smooth.spec, function(x) x$term))))), drop=FALSE])
data.undup <- lapply(data.undup, function(x) as.integer(as.factor(x)))
if (length(data.undup) > 1) for (i in 2:length(data.undup)) data.undup[[1]] <- paste(data.undup[[1]], data.undup[[i]], sep=":")
data.undup <- data.undup[[1]]
gc()
unq.id <- which(!duplicated(data.undup))
data.unq <- data.undup[unq.id]
dup.id <- match(data.undup, data.unq)
}

subsampling <- FALSE

## gams
gams <- list()

if (family %in% c("pp", "ppexi")) 
  data <- .setup.pp.data(data, responsename, pp)

for (i in seq_along(formula)) {
  if (nrow(data) > maxspline) {
    id <- sample(nrow(data), maxspline)
    gams[[i]] <- mgcv::gam(formula[[i]], data=data[id,], fit=FALSE, knots=knots, method="REML")
  } else {
    gams[[i]] <- mgcv::gam(formula[[i]], data=data, fit=FALSE, knots=knots, method="REML")
  }
  gams[[i]] <- .predictable.gam(gams[[i]], formula[[i]])
}

gc()

## likelihood
lik.data <- list()
lik.data$control <- list()
lik.data$outer <- outer
lik.data$control$outer <- list(steptol=1e-12, itlim=1e2, fntol=1e-8, gradtol=1e-2, stepmax=3)
lik.data$control$inner <- list(steptol=1e-12, itlim=1e2, fntol=1e-8, gradtol=1e-4, stepmax=1e2)
lik.data$y <- as.matrix(data[,responsename, drop=FALSE])
lik.data$Mp <- sum(unlist(sapply(gams, function(y) c(1, sapply(y$smooth, function(x) x$null.space.dim)))))
lik.data$const <- .5 * lik.data$Mp * log(2 * pi)
lik.data$nobs <- nrow(lik.data$y)
if (attr(formula, "censored")) {
  lik.data$censored <- TRUE
  # check right-censored values not below left censored values
  if (any(lik.data$y[,2] < lik.data$y[,1]))
    stop("For censored response need right >= left in `cens(left, right)'")
  lik.data$cens.id <- lik.data$y[,2] > lik.data$y[,1]
  if (trace >= 0 & sum(lik.data$cens.id) == 0) {
    message("No response data appear to be censored. Switching to uncensored likelihood.")
    lik.data$censored <- FALSE
  }
} else {
  lik.data$censored <- FALSE
}
if (family == "weibull") {
  if (min(lik.data$y) <= 0) 
    stop(expression("Weibull distribution has support (0, \U221E) in evgam."))
}
if (family == "gpd") {
  if (min(lik.data$y) <= 0) 
    stop(expression("GPD has support (0, \U221E) in evgam."))
}
if (family == "exi") {
if (is.null(exiargs$id)) stop("no `id' in `exi.args'.")
if (is.null(exiargs$nexi)) {
  if (trace >= 0)
    message("`exiargs$nexi' assumed to be 2.")
  exiargs$nexi <- 2
}
if (is.null(exiargs$link)) {
  if (trace >= 0)
    message("`exiargs$link' assumed to be `logistic'.")
  exiargs$link <- "logistic"
}
lik.data$exiname <- exiargs$id
lik.data$y <- list(lik.data$y, data[,exiargs$id])
lik.data$nexi <- exiargs$nexi
if (exiargs$link == "cloglog") {
  lik.data$exilink <- 2
  lik.data$linkfn <- function(x) 1 - exp(-exp(x))
  attr(lik.data$linkfn, "deriv") <- function(x) exp(-exp(x)) * exp(x)
}
if (exiargs$link == "logistic") {
  lik.data$exilink <- 1
  lik.data$linkfn <- function(x) 1 / (1 + exp(-x))
  attr(lik.data$linkfn, "deriv") <- function(x) exp(-x)/(1 + exp(-x))^2
}
if (exiargs$link == "probit") {
  lik.data$exilink <- 0
  lik.data$linkfn <- function(x) pnorm(x)
  attr(lik.data$linkfn, "deriv") <- function(x) dnorm(x)
}
attr(lik.data$linkfn, "name") <- exiargs$link
}
if (family %in% c("pp", "ppexi")) {
  lik.data$ppw <- attr(data, "weights") # point process quadrature weights
  lik.data$y <- as.matrix(rbind(as.matrix(attr(data, "quad")[,responsename]), lik.data$y))
  lik.data$ppq <- rep(as.logical(1:0), c(nrow(attr(data, "quad")), nrow(data))) # identify quadrature points
  lik.data$ppcens <- attr(data, "cens")
  lik.data$weights <- attr(data, "cweights")
  lik.data$exi <- attr(data, "exi")
}
if (family == "ald") {
  if (is.null(aldargs$tau)) 
    aldargs$tau <- .5
  if (is.null(aldargs$C)) 
    aldargs$C <- .5
  lik.data$tau <- aldargs$tau
  lik.data$C <- aldargs$C
}
lik.data$sandwich <- !is.null(sargs$id)
if (lik.data$sandwich) 
  lik.data$sandwich.split <- data[,sargs$id]
  
if (!compact) {
  if (nrow(data) > maxspline) {
    lik.data$X <- .X.evgam(gams, data)
  } else {
    lik.data$X <- .X.evgam(gams)
  }
  if (family %in% c("pp", "ppexi")) {
    ppX <- .X.evgam(gams, attr(data, "quad"))
    lik.data$X <- lapply(seq_along(ppX), function(i) rbind(ppX[[i]], lik.data$X[[i]]))
  }
  lik.data$dupid <- 0
  lik.data$duplicate <- 0
} else {
  if (family %in% c("pp", "ppexi"))
    stop("Option compact = TRUE not currently possible for pp model.")
  lik.data$X <- .X.evgam(gams, data[unq.id,])#lapply(gams, .X.evgam, newdata=data[unq.id,])
  lik.data$dupid <- dup.id - 1
  lik.data$duplicate <- 1
}
for (i in seq_along(gams)) {
  if (removeData) 
    gams[[i]]$y <- NULL
}
if (length(lik.data$X) == 1 & length(nms) > 1) {
  for (i in 2:length(nms)) {
    lik.data$X[[i]] <- lik.data$X[[1]]
    gams[[i]] <- gams[[1]]
  }
}
nbk <- sapply(lik.data$X, ncol)
lik.data$nb <- sum(nbk)
lik.data$idpars <- rep(seq_along(lik.data$X), nbk)
lik.data$LAid <- lik.data$idpars > 0
lik.data$subsampling <- subsampling
gotsmooth <- which(sapply(gams, function(x) length(x$sp)) > 0)
lik.data$k <- 1
if (is.null(sargs$id)) {
  lik.data$adjust <- 0
} else {
  if (is.null(sargs$method)) 
    sargs$method <- "magnitude"
  if (sargs$method == "curvature") {
    if (trace > 0) 
      message(paste("Sandwich adjustment method: curvature"))
    lik.data$adjust <- 2
  } else {
    if (trace > 0) 
      message(paste("Sandwich adjustment method: magnitude"))
    lik.data$adjust <- 1
  }
}
if (is.null(sargs$force)) 
  sargs$force <- FALSE
lik.data$force <- sargs$force
list(lik.data=lik.data, gotsmooth=gotsmooth, data=data, gams=gams, sandwich=lik.data$adjust > 0)
}

############ .setup.pp.data ##########################

.setup.pp.data <- function(data, responsename, pp) {

nodes <- pp$nodes
ny <- pp$ny
if (is.null(ny))
  stop("Cannot have NULL pp.args$ny.")
threshold <- pp$threshold
r <- pp$r
if (is.null(threshold) & is.null(r)) 
  stop("Both pp$threshold and pp$r cannot be NULL")

## simple constant partial point process
data$row <- seq_len(nrow(data))
ds <- split(data, data[,pp$id])
wts <- pp$ny
if (length(wts) == 1) {
  wts <- rep(wts, length(ds))
} else {
  wts <- wts[match(names(ds), names(wts))]
}
nobs2 <- sapply(ds, nrow)
## start of original r-largest order statistic stuff
data.quad <- do.call(rbind, lapply(ds, function(x) x[1,]))
if (!is.null(pp$r)) {
enough <- nobs2 >= pp$r
if (any(!enough)) warning(paste(sum(!enough), "unique pp.args$id removed for having fewer than r observations."))
ds <- ds[enough]
wts <- wts[enough]
nid <- sum(enough)
data.quad <- data.quad[enough,]
if (pp$r != -1) {
    du <- sapply(ds, function(x) x[order(x[,responsename], decreasing=TRUE)[pp$r], responsename])
} else {
    du <- sapply(ds, function(x) min(x[, responsename]))
}
## end of ...
} else {
## start of specified threshold stuff
du <- sapply(ds, function(x) x[1, pp$threshold])
nid <- length(du)
## end of specified threshold stuff
}
data.quad[,responsename] <- du
ds <- lapply(seq_len(nid), function(i) subset(ds[[i]], ds[[i]][,responsename] >= du[i]))
out <- dfbind(ds)
attr(out, "weights") <- wts
attr(out, "quad") <- data.quad
if (is.null(pp$cens)) {
  attr(out, "cens") <- NULL
} else {
  attr(out, "cens") <- data[out$row, pp$cens]
}
if (is.null(pp$weights)) {
  attr(out, "cweights") <-rep(1, length(out$row))
} else {
  attr(out, "cweights") <- pp$weights[out$row]
}
out
}

############ .sandwich.C ##########################

.sandwich.C <- function(H, J) {
iJ <- pinv(J)
HA <- crossprod(H, crossprod(iJ, H))
sH <- svd(H)
M <- sqrt(sH$d) * t(sH$v)
sHA <- svd(HA)
MA <- sqrt(sHA$d) * t(sHA$v)
solve(M, MA)
}

############ .setup.inner.inits ##########################

.setup.inner.inits <- function(inits, likdata, likfns, npar, family) {

likdata0 <- likdata
likdata0$X <- lapply(seq_along(likdata$X), function(i) matrix(1, nrow=nrow(likdata$X[[i]]), ncol=1))
# likdata0$pp$X <- lapply(seq_along(likdata$pp$X), function(x) matrix(1, nrow=nrow(likdata$pp$X[[x]]), ncol=1))
likdata0$S <- diag(0, npar)
likdata0$idpars <- seq_len(npar)

if (is.null(inits)) {
  if (npar == 1) 
    inits <- 2
  if (npar == 2) {
    if (family == "ald") {
      inits <- c(quantile(likdata0$y[,1], likdata0$tau), log(sd(likdata0$y[,1])))
    } else {
      inits <- c(log(mean(likdata$y[,1])), .05)
      if (family == "transxigpd") 
        inits[2] <- .9
    }
  }
  if (npar %in% 3:4) {
    inits <- c(sqrt(6) * sd(likdata0$y[,1]) / pi, .05)
    inits <- c(mean(likdata0$y[,1]) - .5772 * inits[1], log(inits[1]), inits[2])
    if (npar == 4) 
      inits <- c(inits, 1)
  }
  if (npar == 6) {
    inits <- c(sqrt(6) * sd(likdata0$y[,1]) / pi, .05)
    inits <- c(mean(likdata0$y[,1]) - .5772 * inits[1], log(inits[1]), inits[2])
    inits <- c(inits, 0, 0, 1)
  }
  likdata0$CH <- diag(length(inits))
  likdata0$compmode <- numeric(length(inits))
  beta0 <- .newton_step_inner(inits, .nllh.nopen, .search.nopen, likdata=likdata0, likfns=likfns, control=likdata$control$inner)$par
} else {
  if (is.list(inits)) {
    betamat <- expand.grid(inits)
    betanllh <- numeric(nrow(betamat))
    for (i in seq_len(nrow(betamat))) {
      beta0 <- unlist(betamat[i,])
      betanllh[i] <- likfns$nllh(beta0, likdata0)
    }
    beta0 <- betamat[which.min(betanllh),]
    print(beta0)
  } else {
    beta0 <- inits
  }
}
beta0 <- unlist(lapply(seq_len(npar), function(i) c(beta0[i], rep(0, ncol(likdata$X[[i]]) - 1))))
compmode <- 0 * beta0
CH <- diag(compmode + 1)
k <- 1
likdata[c("k", "CH", "compmode")] <- list(k, CH, compmode)
diagH <- diag(.gH.nopen(beta0, likdata=likdata, likfns=likfns)[[2]])
if (likdata$sandwich) {
  beta0 <- .newton_step(beta0, .nllh.nopen, .search.nopen, likdata=likdata, likfns=likfns, control=likdata$control$inner)$par
  H <- .gH.nopen(beta0, likdata=likdata, likfns=likfns, sandwich=TRUE)
  if (family == "pp") {
    J0 <- H[[1]]
    J <- J0[,!likdata$ppq]
    J0 <- rowSums(J0[,likdata$ppq])
    J <- split(as.data.frame(t(J)), likdata$sandwich.split)
    wts <- sapply(J, nrow)
    wts <- wts / sum(wts)
    J <- sapply(J, colSums)
    J <- J + J0 %o% wts
    J <- tcrossprod(J)
  } else {
    J <- split(as.data.frame(t(H[[1]])), likdata$sandwich.split)
    J <- sapply(J, colSums)
    J <- tcrossprod(J)
  }
  H <- H[[2]]
  diagH <- diag(H)
  cholH <- try(chol(H), silent=TRUE)
  if (inherits(cholH, "try-error")) {
    if (!likdata$force) {
      stop("Hessian of unpenalised MLE not positive definite.\n  Supply `force=TRUE' to `sandwich.args' to perturb it to be positive definite.")
    } else {
      if (trace >= 0)
        message("Hessian perturbed to be positive definite for sandwich adjustment.")
      iH <- pinv(H)
    }
  } else {
    iH <- chol2inv(cholH)
  }
  if (likdata$adjust == 2) {
    cholJ <- try(chol(J), silent=TRUE)
    if (inherits(cholJ, "try-error") & likdata$adjust == 2) {
      HA <- crossprod(backsolve(cholJ, H, transpose=TRUE))
    } else {
      iHA <- tcrossprod(crossprod(iH, J), iH)
      choliHA <- try(chol(iHA), silent=TRUE)
      if (inherits(choliHA, "try-error")) {
        if (!likdata$force) {
          stop("Sandwich variance not positive definite.\n  Supply `force=TRUE' to `sandwich.args' to perturb it to be positive definite.")
        } else {
          if (trace >= 0)
            message("Sandwich variance perturbed to be positive definite.")
          HA <- pinv(iHA)
        }
      } else {
        HA <- chol2inv(choliHA)
      }
    }
    sH <- svd(H)
    M <- sqrt(sH$d) * t(sH$v)
    sHA <- svd(HA)
    MA <- sqrt(sHA$d) * t(sHA$v)
    CH <- solve(M, MA)
    compmode <- beta0
  } else {
    k <- 1 / mean(diag(crossprod(iH, J)))
  }
}
attr(beta0, "k") <- k
attr(beta0, "CH") <- CH
attr(beta0, "compmode") <- compmode
attr(beta0, "diagH") <- diagH
beta0
}

############ .guess ##########################

.guess <- function(x, d, s) {
okay <- s != 0
val <- d / (d + exp(x) * s)
mean(val[okay]) - .4
}

############ .sandwich ##########################

.sandwich <- function(likdata, beta) {
likdata$k <- attr(beta, "k")
likdata$CH <- attr(beta, "CH")
likdata$compmode <- attr(beta, "compmode")
bigX <- do.call(cbind, likdata$X)
CHX <- bigX %*% likdata$CH
CHX <- lapply(unique(likdata$idpars), function(i) CHX[,likdata$idpars == i])
likdata$CHX <- CHX
likdata
}

############ .outer ##########################

.outer <- function(rho0, beta, likfns, likdata, Sdata, control, correctV, outer, trace) {

attr(rho0, "beta") <- beta

if (outer == "newton") {
  fit.reml <- .newton_step_inner(rho0, .reml0, .search.reml, likfns=likfns, likdata=likdata, Sdata=Sdata, control=likdata$control$outer, trace=trace > 1)
} else {
  if (outer == "fd") {
    fit.reml <- .BFGS(rho0, .reml0, .reml1.fd, likfns=likfns, likdata=likdata, Sdata=Sdata, control=likdata$control$outer, trace=trace > 1)
  } else {
    fit.reml <- .BFGS(rho0, .reml0, .reml1, likfns=likfns, likdata=likdata, Sdata=Sdata, control=likdata$control$outer, trace=trace > 1)
  }
  rho1 <- fit.reml$par
  attr(rho1, "beta") <- fit.reml$beta
  fit.reml$Hessian <- try(.reml12(rho1, likfns=likfns, likdata=likdata, Sdata=Sdata)[[2]], silent=TRUE)
  if (inherits(fit.reml$Hessian, "try-error")) 
    fit.reml$Hessian <- .reml2.fd(rho1, likfns=likfns, likdata=likdata, Sdata=Sdata)
}

fit.reml$invHessian <- .solve_evgam(fit.reml$Hessian)

fit.reml$trace <- trace

if (trace == 1) {
  report <- "\n Final max(|grad|))"
  likdata$S <- .makeS(Sdata, exp(fit.reml$par))
  report <- c(report, paste("   Inner:", signif(max(abs(.gH.pen(fit.reml$beta, likdata, likfns)[[1]])), 3)))
  report <- c(report, paste("   Outer:", signif(max(abs(fit.reml$gradient)), 3)))
  report <- c(report, "", "")
  cat(paste(report, collapse="\n"))
}

fit.reml

}

############ .outer.nosmooth ##########################

.outer.nosmooth <- function(beta, likfns, likdata, control, trace) {

fit.inner <- .newton_step(beta, .nllh.nopen, .search.nopen, likdata=likdata, likfns=likfns, control=likdata$control$inner)

list(beta=fit.inner$par)

}

############ .VpVc ##########################

.VpVc <- function(fitreml, likfns, likdata, Sdata, correctV, sandwich, smooths, trace) {
lsp <- fitreml$par
H0 <- .gH.nopen(fitreml$beta, likdata, likfns)[[2]]
if (smooths) {
  sp <- exp(lsp)
  H <- H0 + likdata$S
} else {
  H <- H0
}
cholH <- try(chol(H), silent=TRUE)
if (inherits(cholH, "try-error") & trace >= 0)
  message("Final Hessian of negative penalized log-likelihood not numerically positive definite.")
Vc <- Vp <- pinv(H)
if (smooths) {
if (correctV) {
cholVp <- try(chol(Vp), silent=TRUE)
if (inherits(cholVp, "try-error")) {
    cholVp <- attr(.perturb(Vp), "chol")
}
attr(lsp, "beta") <- fitreml$beta
spSl <- Map("*", attr(Sdata, "Sl"), exp(lsp))
dbeta <- .d1beta(lsp, fitreml$beta, spSl, .Hdata(H))$d1
Vrho <- fitreml$invHessian
Vbetarho <- tcrossprod(dbeta %*% Vrho, dbeta)
# eps <- 1e-4
# R0 <- .grad.R(fitreml$par, Sdata=Sdata, R=0, eps=1, likfns=likfns, likdata=likdata, H0=H0)
# dR <- lapply(seq_along(sp), function(i) grad.R(replace(fitreml$par, i, fitreml$par[i] + eps), Sdata=Sdata, R=R0, eps=eps, likfns=likfns, likdata=likdata, H0=H0))
VR <- matrix(0, nrow=likdata$nb, ncol=likdata$nb)
# for (k in seq_along(sp)) for (l in seq_along(sp)) VR <- VR + crossprod(dR[[k]] * Vrho[k, l], dR[[l]])
# VR <- .5 * (VR + t(VR))
Vc <- .perturb(Vp + Vbetarho + VR)
}
} else {
  Vrho <- 0
}
list(Vp=Vp, Vc=Vc, Vlsp=Vrho, H0=H0, H=H)
}

############ .edf ##########################

.edf <- function(beta, likfns, likdata, VpVc, sandwich) {
diag(crossprod(VpVc$Vp, VpVc$H0))
}

############ .swap ##########################

.swap <- function(fitreml, gams, likdata, VpVc, gotsmooth, edf, smooths) {
Vp <- VpVc$Vp
Vc <- VpVc$Vc
if (smooths) {
  # bug fixed in evgam_0.1.2 if any(diff(gotsmooth)) != 1
  spl <- split(exp(fitreml$par), unlist(sapply(seq_along(gams), function(x) rep(x, length(gams[[x]]$sp))))) 
  sp <- replace(lapply(seq_along(gams), function(x) NULL), gotsmooth, spl)
}
for (i in seq_along(gams)) {
  idi <- likdata$idpars == i
  gams[[i]]$coefficients <- fitreml$beta[idi]
  names(gams[[i]]$coefficients) <- gams[[i]]$term.names
  gams[[i]]$Vp <- Vp[idi, idi, drop = FALSE]
  gams[[i]]$Vc <- Vc[idi, idi, drop = FALSE]
  if (i %in% gotsmooth) gams[[i]]$sp <- sp[[i]]
  gams[[i]]$edf <- edf[idi]
}
gams
}

############ .finalise ##########################

.finalise <- function(gams, data, likfns, likdata, Sdata, fitreml, VpVc, family, gotsmooth,
formula, responsenm, removeData, edf) {
nms <- c("location", "logscale", "shape")
if (length(gams) == 2) {
  if (family %in% c("ald", "gauss")) {
    nms <- nms[1:2]
  } else {
    nms <- nms[-1]
}}
if (length(gams) == 4) {
  if (is.null(likdata$agg)) {
    nms <- c(nms, "logitdep")
  } else {
    nms <- c(nms, "logdep")
  }
}
if (family == "exponential") nms <- "lograte"
if (family == "weibull") nms[2] <- "logshape"
if (family == "exi") nms <- paste(attr(likdata$linkfn, "name"), "exi", sep="")
names(gams) <- nms
smooths <- length(gotsmooth) > 0
Vp <- VpVc$Vp
Vc <- VpVc$Vc
if (smooths) gams$sp <- exp(fitreml$par)
gams$nobs <- likdata$nobs
gams$logLik <- -1e20
fit.lik <- list(convergence=0)
if (fit.lik$convergence == 0) {
gams$logLik <- -.nllh.nopen(fitreml$beta, likdata, likfns)
gams$logLik <- gams$logLik - likdata$const
}
if (fit.lik$convergence != 0) gams$AIC <- gams$BIC <- 1e20
attr(gams, "df") <- sum(edf)
gams$simulate <- list(mu=fitreml$beta, Sigma=Vp)
gams$family <- family
gams$idpars <- likdata$idpars
# tidy up print names a bit
nms <- names(gams)[seq_along(formula)]
logits <- substr(nms, 1, 5) == "logit"
if (any(logits))
  nms[logits] <- gsub("logit", "", nms[logits])
logs <- substr(nms, 1, 3) == "log"
if (any(logs))
  nms[logs] <- gsub("log", "", nms[logs])
probits <- substr(nms, 1, 6) == "probit"
if (any(probits))
  nms[probits] <- gsub("probit", "", nms[probits])
# end name tidying
gams$predictor.names <- attr(formula, "predictor.names")
formula <- attr(formula, "stripped")
names(formula) <- nms
gams$call <- formula
gams$response.name <- responsenm
gams$gotsmooth <- gotsmooth
if (!removeData) {
    if (family == "pp") {
        gams$data <- attr(data, "quad")
    } else {
        gams$data <- data
    }
}
gams$Vc <- Vc
gams$Vp <- Vp
gams$Vlsp <- VpVc$Vlsp
gams$negREML <- fitreml$objective
gams$coefficients <- fitreml$beta
if (family == "ald") gams$tau <- likdata$tau
if (family == "exi") {
  gams$linkfn <- likdata$linkfn
  gams$exi.name <- likdata$exiname
}
for (i in seq_along(likdata$X)) {
gams[[i]]$X <- likdata$X[[i]]
if (likdata$duplicate == 1) gams[[i]]$X <- gams[[i]]$X[likdata$dupid + 1,]
gams[[i]]$fitted <- as.vector(likdata$X[[i]] %*% gams[[i]]$coefficients)
}
gams$likdata <- likdata
gams$likfns <- likfns
if (smooths) gams$Sdata <- Sdata
gams$formula <- formula
gams$compacted <- likdata$duplicate == 1
if (gams$compacted) gams$compactid <- likdata$dupid + 1
smooth.terms <- unique(lapply(lapply(gams[gotsmooth], function(x) x$smooth), function(y) lapply(y, function(z) z$term)))
smooth.terms <- unique(unlist(smooth.terms, recursive=FALSE))
gams$plotdata <- lapply(smooth.terms, function(x) unique(data[,x, drop=FALSE]))
if (family == "weibull") names(gams)[2] <- "logshape"
if (family == "exponential") names(gams)[1] <- "lograte"
gams$ngam <- length(formula)
for (i in seq_along(gams[nms])[-gotsmooth])
  gams[[i]]$smooth <- NULL
class(gams) <- "evgam"
return(gams)
}

