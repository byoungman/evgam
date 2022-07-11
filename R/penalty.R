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
out <- out + .5 * crossprod(pars, likdata$S %*% pars)[1, 1]
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

.Hdata <- function(H) {
out <- list(H0=H)
H2 <- .precondition(H)
H2 <- .perturb(H2)
out$H <- H2
out$dH <- attr(H2, "d")
out$cH <- attr(H2, "chol")
out$iH <- crossprod(backsolve(out$cH, diag(out$dH), transpose=TRUE))
out$kept <- !logical(nrow(H))
out
}

## gradient and Hessians of negative log-likelihoods

.gH <- function(x, likdata, sandwich=FALSE, deriv=2) {
nX <- length(likdata$X)
if (nX == 1) {
  out <- .gH1(x, likdata$X[[1]], likdata$dupid, likdata$duplicate, as.integer(sandwich), deriv)
} else {
  if (nX == 2) {
    out <- .gH2(x, likdata$X[[1]], likdata$X[[2]], likdata$dupid, likdata$duplicate, as.integer(sandwich), deriv)
  } else {
    if (nX == 3) {
      out <- .gH3(x, likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$dupid, likdata$duplicate, as.integer(sandwich), deriv)
    } else {
      if (nX == 4) { # added with evgam_0.1.2 (05/04/2020)
        out <- .gH4(x, likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$dupid, likdata$duplicate, as.integer(sandwich), deriv)
      } else {
        if (nX == 5) { # added with evgam_0.1.5 (29/06/2021)
          out <- .gH5(x, likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$dupid, likdata$duplicate, as.integer(sandwich), deriv)
        } else {
          stop("Number of model parameters not in {1, 2, 3, 4, 5}")
        }
      }
    }
  }
}
out
}

.gH.nopen <- function(pars, likdata, likfns, sandwich=FALSE, deriv=2) {
pars <- as.vector(likdata$compmode + likdata$CH %*% (as.vector(pars) - likdata$compmode))
temp <- likfns$d120(pars, likdata)
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

.gH.pen <- function(pars, likdata, likfns, deriv=2) {
temp <- .gH.nopen(pars, likdata, likfns, deriv=deriv)
temp[[1]] <- temp[[1]] + crossprod(pars, likdata$S)[1, ]
if (deriv > 1) {
  attr(temp, "PP") <- attr(temp, "PP") + likdata$S / norm(likdata$S, "F")
  temp[[2]] <- temp[[2]] + likdata$S
}
temp
}

## Newton search directions

.search.dir <- function(g, H, kept=NULL) {
if (is.null(kept)) 
  kept <- !logical(length(g))
H0 <- H
g0 <- g
g <- g[kept]
H <- H[kept, kept, drop=FALSE]
if (any(!is.finite(g))) 
  stop("Some gradient non-finite")
if (any(!is.finite(H))) 
  stop("Some Hessian non-finite")
H2 <- .precondition(H)
# okay <- is.finite(attr(H2, "d"))
# bad <- !all(okay)
# if (bad) {
#   print("bad")
#   d <- attr(H2, "d")
#   H2 <- H2[okay, okay, drop=FALSE]
#   attr(H2, "d") <- d[okay]
# }
H2 <- .perturb(H2)
R <- attr(H2, "chol")
d <- attr(H2, "d")
out <- numeric(length(kept))
piv <- ipiv <- attr(R, "pivot")
ipiv[piv] <- seq_len(length(piv))
out[kept] <- d * backsolve(R, forwardsolve(R, (g * d)[piv], upper.tri=TRUE, transpose=TRUE))[ipiv]
g0[!kept] <- 0
attr(out, "gradient") <- g0
attr(out, "Hessian") <- H0
attr(out, "cholH") <- R
attr(out, "rank") <- qr(H0)$rank
out
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
attr(out, "Hessian") <- attr(gH, "Hessian")
out
}
