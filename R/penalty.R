## worker functions

add_in <- function(x, i, v) {
out <- numeric(length(x) + length(i))
out[-i] <- x
out[i] <- v
out
}

## negative log-likelihoods

.nllh.nopen <- function(pars, likdata, likfns, newton=TRUE, id=NULL, vals=NULL) {
if (!is.null(id)) 
  pars <- add_in(pars, id, vals)
pars <- as.vector(likdata$compmode + likdata$CH %*% (pars - likdata$compmode))
out <- likdata$k * likfns$d0(pars, likdata)
if (!is.finite(out)) 
  out <- 1e20
out
}

.nllh.pen <- function(pars, likdata, likfns, newton=TRUE) {
out <- .nllh.nopen(pars, likdata, likfns, newton)
out <- out + .5 * crossprod(pars, as.vector(likdata$S %*% pars))[1, 1]
out
}

## gradients of negative log-likelihoods

.grad.nopen <- function(pars, likdata, likfns, id=NULL, vals=NULL) {
if (!is.null(id)) 
  pars <- add_in(pars, id, vals)
.gH.nopen(pars, likdata, likfns, sandwich=FALSE, deriv=1)[[1]][-id]
}

.grad.pen <- function(pars, likdata, likfns) {
out <- crossprod(pars, likdata$S)[1, ]
pars <- as.vector(likdata$compmode + likdata$CH %*% (pars - likdata$compmode))
temp <- likfns$d120(pars, likdata)
temp <- as.vector(.gH(temp, likdata, deriv=1)[[1]])
temp <- likdata$k * temp
temp <- as.vector(likdata$CH %*% temp)
as.vector(out + temp)
}

## Hessians of negative log-likelihoods

.hess.nopen <- function(pars, likdata, likfns) {
pars <- as.vector(likdata$compmode + likdata$CH %*% (as.vector(pars) - likdata$compmode))
temp <- likfns$d120(pars, likdata)
temp <- .gH(temp, likdata)[[2]]
temp <- likdata$k * temp
temp <- crossprod(likdata$CH, temp) %*% likdata$CH
temp
}

.hess.pen <- function(pars, likdata, likfns) {
.hess.nopen(pars, likdata, likfns) + likdata$S
}

.gH.nopen <- function(pars, likdata, likfns, sandwich=FALSE, deriv=2) {
pars <- as.vector(likdata$compmode + likdata$CH %*% (as.vector(pars) - likdata$compmode))
if ('sandwich' %in% names(formals(likfns$d120))) {
  temp <- likfns$d120(pars, likdata, sandwich)
} else {
  temp <- likfns$d120(pars, likdata)
}
if (is.null(likdata$agg))
  temp <- .gH(temp, likdata, sandwich, deriv)
temp[[1]] <- likdata$k * temp[[1]]
temp[[1]] <- t(temp[[1]] %*% likdata$CH)
if (deriv > 1) {
  temp[[2]] <- likdata$k * temp[[2]]
  temp[[2]] <- crossprod(likdata$CH, temp[[2]]) %*% likdata$CH
  attr(temp, "PP") <- temp[[2]] / norm(temp[[2]], "F")
}
temp
}

.search.nopen <- function(pars, likfns, likdata, kept, newton=TRUE, id=NULL, vals=NULL) {
if (!is.null(id)) 
  pars <- add_in(pars, id, vals)
gH <- .gH.nopen(pars, likdata, likfns)
if (!is.null(id)) {
  gH[[1]] <- gH[[1]][-id, drop=FALSE]
  gH[[2]] <- gH[[2]][-id, -id, drop=FALSE]
}
if (newton) {
  out <- .search.dir(gH[[1]], gH[[2]], kept)
} else {
  out <- gH[[1]]
  attr(out, "gradient") <- gH[[1]]
}
attr(out, "PP") <- attr(gH, "PP")
out
}

.search.pen <- function(pars, likfns, likdata, kept, newton=TRUE) {
gH <- .gH.pen(pars, likdata, likfns)
if (newton) {
  out <- .search.dir(gH[[1]], gH[[2]], kept)
} else {
  out <- gH[[1]]
  attr(out, "gradient") <- gH[[1]]
}
attr(out, "PP") <- attr(gH, "PP")
out
}
