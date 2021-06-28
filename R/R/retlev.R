#' Quantile estimation of a composite extreme value distribution
#'
#' @param p a scalar giving the quantile of the distribution sought
#' @param loc a scalar, vector or matrix giving the location parameter
#' @param scale as above, but scale parameter
#' @param shape as above, but shape parameter
#' @param m a scalar giving the number of values per return period unit, e.g. 365 for daily data giving annual return levels
#' @param alpha a scalar, vector or matrix of weights if within-block variables not identically distributed and of different frequencies
#' @param theta a scalar, vector or matrix of extremal index values
#' @param family a character string giving the family for which return levels sought
#' @param tau a scalar, vector or matrix of values giving the threshold quantile for the GPD (i.e. 1 - probability of exceedance)
#' @param start a 2-vector giving starting values that bound the return level
#'
#' @details
#'
#' If \eqn{F} is the generalised extreme value or generalised Pareto
#' distribution, \code{qev} solves 
#' \deqn{\prod_{j=1}^n \big\{F(z)\}^{m \alpha_j \theta_j} = p.}
#'
#' For both distributions, location, scale and shape parameters 
#' are given by \code{loc}, \code{scale} and \code{shape}. The 
#' generalised Pareto distribution, for \eqn{\xi \neq 0} and \eqn{z > u}, 
#' is parameterised as \eqn{1 - (1 - \tau) [1 + \xi (z - u) / \psi_u]^{-1/\xi}},
#' where \eqn{u}, \eqn{\psi_u} and \eqn{\xi} are its location, scale and shape
#' parameters, respectively, and \eqn{\tau} corresponds to argument \code{tau}.
#'
#' @examples
#'
#' qev(0.9, c(1, 2), c(1, 1.1), .1, family="gev")
#' qev(0.99, c(1, 2), c(1, 1.1), .1, family="gpd", tau=0.9)
#'
#' @return A scalar or vector of estimates of \code{p}
#'
#' @export
#'
qev <- function(p, loc, scale, shape, m=1, alpha=1, theta=1, family, tau=0, start=NULL) {
if (!(family %in% c("gev", "gpd"))) stop("Invalid family")
loc <- as.matrix(loc)
scale <- as.matrix(scale)
shape <- as.matrix(shape)
nr <- max(nrow(loc), nrow(scale), nrow(shape))
if (length(theta) == 1) theta <- rep(theta, nr)
if (length(alpha) == 1) alpha <- rep(alpha, nr)
alpha <- alpha / sum(alpha)
weights <- m * alpha * theta
nc <- max(ncol(loc), ncol(scale), ncol(shape))
loc <- matrix(loc, nr, nc)
scale <- matrix(scale, nr, nc)
shape <- matrix(shape, nr, nc)
tau <- matrix(tau, nr, nc)
theta <- matrix(theta, nr, nc)
weights <- matrix(weights, nr, nc)
out <- numeric(nc)
for (i in seq_len(nc)) {
  out[i] <- .rlvec(p, loc[,i], scale[,i], shape[,i], m, nr, weights[,i], family, tau[,i], theta[,i], start)
}
out
}

.rlvec <- function(p, loc, scale, shape, m, n, weights, family, tau, theta, start=NULL) {
if (is.null(start)) {
if (family == "gev") {
    start <- range(.qgev(p^(1/n), loc, scale, shape))
} else {
    start <- range(.qgpd(p, loc, scale, shape, 1 - tau, theta, m))
}
dstart <- diff(start)
nullstart <- FALSE
} else {
nullstart <- TRUE
}
while(.rlroot(start[1], loc, scale, shape, tau, weights, p, TRUE, family) > 0) {
    if (nullstart) stop("Invalid `start' values given")
    start[1] <- start[1] - .2 * dstart
}
while(.rlroot(start[2], loc, scale, shape, tau, weights, p, TRUE, family) < 0) {
    if (nullstart) stop("Invalid `start' values given")
    start[2] <- start[2] + .2 * dstart
}
opt <- uniroot(.rlroot, start, loc=loc, scale=scale, shape=shape, tau=tau, weights=weights, p=p, TRUE, family)
opt$root
}

.pmax <- function(x, loc, scale, shape, tau, weights, log, family) {
if (family == "gpd") {
    out <- .pgpd(x, loc, scale, shape, tau, FALSE, TRUE)
} else {
    out <- .pgev(x, loc, scale, shape, FALSE, TRUE)
}
out <- sum(weights * out)
if (!log) out <- exp(out)
out
}

.rlroot <- function(x, loc, scale, shape, tau, weights, p, log, family) {
if (log) p <- log(p)
.pmax(x, loc, scale, shape, tau, weights, log, family) - p
}
