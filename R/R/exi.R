## Exponential negative log-likelihood functions

.exi.d0 <- function(pars, likdata) {
exid0(likdata$y[[1]], likdata$y[[2]], pars, likdata$nexi, likdata$X[[1]], likdata$dupid, likdata$duplicate, likdata$exilink)
}

.exi.d12 <- function(pars, likdata) {
exid12(likdata$y[[1]], likdata$y[[2]], pars, likdata$nexi, likdata$X[[1]], likdata$dupid, likdata$duplicate, likdata$exilink)
}

.exi.d34 <- function(pars, likdata) {
exid34(likdata$y[[1]], likdata$y[[2]], pars, likdata$nexi, likdata$X[[1]], likdata$dupid, likdata$duplicate, likdata$exilink)
}

.exifns <- list(d0=.exi.d0, d120=.exi.d12, d340=.exi.d34)

.runmaxgrp <- function(df, ynm, n) {
# function to replace observations in data frame
# with running n-observation maximum
ndf <- nrow(df)
if (ndf >= n) {
id <- seq_len(ndf - n + 1) + round(.5 * n) - 1
df[id, ynm] <- runmax(df[, ynm], n)
return(df[id,])
} else {
return(NULL)
}
}

#' @param data a data frame
#' @param cons a character string for the variable in \code{data} that identifies consecutive observations
#' @param ynm a character string for the variable in \code{data} that is the observations
#'
#' @return \code{dfrunmax} returns a data frame with observations swapped for \eqn{n}-observation running maximum
#'
#' @rdname runmax
#'
#' @export
#'
dfrunmax <- function(data, cons, ynm, n=2) {
# function to split data into consecutive parts
# and then swap observations with running 
# 2-observation maximum
t1 <- data[,cons]
ind <- c(1, 1 + cumsum(diff(t1) > 1))
t2 <- split(data, ind)
out <- lapply(seq_along(t2), function(x) NA)
some <- sapply(t2, length) > 1
for (i in seq_along(t2)) {
  if (some[i]) {
    out[[i]] <- .runmaxgrp(t2[[i]], ynm=ynm, n=n)
  } else {
    out[[i]] <- t2[[i]]
    out[[i]][,ynm] <- NA
  }
}
do.call(rbind, out)
}

.cons2split <- function(x) {
dx <- diff(x)
ends <- which(dx != 1)
starts <- c(1, ends + 1)
ends <- c(ends, length(x))
rep(seq_along(starts), ends - starts + 1)
}

#' Estimate extremal index using `intervals' method
#'
#' @param x a logical vector or list of logical vectors
#' @param y an integer vector the same length as \code{x}; see Details
#'
#' @return A scalar estimate of the extremal index
#'
#' @examples
#'
#' n <- 1e2
#' x <- runif(n)
#' extremal(x > .9)
#' 
#' y <- sort(sample(n, n - 5))
#' x2 <- x[y]
#' extremal(x2 > .9, y)
#'
#' @details
#' 
#' Intervals estimator of extremal index based on Ferro and Segers (2003)'s moment-based estimator.
#'
#' If \code{x} is supplied and \code{y} is not, \code{x} is assumed to identify consecutive threshold exceedances.
#' If \code{x} is supplied as a list, each list element is assumed to comprise identifiers of consecutive exceedances.
#' If \code{y} is supplied, \code{x} must be a logical vector, and \code{y} gives positions of \code{x} in
#' its original with-missing-values vector: so \code{y} identifies consecutive \code{x}.
#'
#' @references 
#' 
#' Ferro, C. A., & Segers, J. (2003). Inference for clusters of extreme values. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 65(2), 545-556.
#'
#' @rdname extremal
#'
#' @export
#'
extremal <- function(x, y=NULL) {
if (!is.null(y)) 
  x <- split(x, .cons2split(y))
if (!inherits(x, "list")) 
  x <- list(x)
x <- x[sapply(x, sum) > 1]
Tu <- unlist(lapply(x, function(z) diff(seq_along(z)[z])))
Nu <- sum(unlist(x))
if (!any(Tu > 2)) {
  hold <- colSums(cbind(Tu, Tu^2))
} else {
  hold <- colSums(cbind(Tu - 1, (Tu - 1) * (Tu - 2)))
}
pmin(2 * (hold[1])^2/((Nu - 1) * hold[2]), 1)
}
