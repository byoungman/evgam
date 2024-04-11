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
#' @param bgev.args a list specifying parameters of the blended GEV distribution; see Details
#' @param start a 2-vector giving starting values that bound the return level
#' @param method a character string giving the numerical estimation procedure; defaults to \code{"uniroot"}
#'
#' @details
#'
#' If \eqn{F} is the generalised extreme value, generalised Pareto or blended
#' GEV distribution, \code{qev} solves 
#' \deqn{\prod_{j=1}^n \big\{F_i(z)\}^{m \alpha_j \theta_j} = p.}
#' for \eqn{i = 1, \ldots, k}. So vectors are supplied as $n$-vectors and matrices 
#' are supplied as \eqn{n \times k} matrices.
#'
#' For all distributions, location, scale and shape parameters are given by 
#' \code{loc}, \code{scale} and \code{shape}. The generalised Pareto 
#' distribution, for \eqn{\xi \neq 0} and \eqn{z > u}, is parameterised as 
#' \eqn{1 - (1 - \tau) [1 + \xi (z - u) / \psi_u]^{-1/\xi}},
#' where \eqn{u}, \eqn{\psi_u} and \eqn{\xi} are its location, scale and shape
#' parameters, respectively, and \eqn{\tau} corresponds to argument \code{tau}. 
#' For the blended GEV distribution \code{pa}, \code{pb}, \code{alpha} and 
#' \code{beta} specify additional parameters of the blended GEV distribution; see
#' \code{family.evgam} for details.
#' 
#' Estimates either use function \code{uniroot} or \code{method = "newton"} uses
#' the Newton-Rhaphson method. The latter is often much quicker if matrices
#' are supplied, i.e. for \eqn{k > 1}.
#'
#' @examples
#'
#' qev(0.9, c(1, 2), c(1, 1.1), .1, family = "gev")
#' qev(0.99, c(1, 2), c(1, 1.1), .1, family = "gpd", tau = 0.9)
#'
#' # an example representative on monthly estimates at two locations
#' qev(0.9, matrix(c(1:12, 2:13), 12, 2), 1.1, .1, family = "gev")
#'
#' # a blended GEV example with default blended GEV specification
#' qev(0.9, matrix(c(1:12, 2:13), 12, 2), 1.1, .1, family = "bgev")

#' @return A scalar or vector of estimates of \code{p}
#'
#' @export
#'
qev <- function(p, loc, scale, shape, m = 1, alpha = 1, theta = 1, family, tau = 0, 
                bgev.args = list(pa = 0.05, pb = 0.2, alpha = 0.5, beta = 0.5), 
                start = NULL, method = "uniroot") {
# if (!(family %in% c("gev", "gpd", "bgev"))) stop("Invalid family")
# loc <- as.matrix(loc)
# scale <- as.matrix(scale)
# shape <- as.matrix(shape)
# nr <- max(nrow(loc), nrow(scale), nrow(shape))
# if (length(theta) == 1) theta <- rep(theta, nr)
# if (length(alpha) == 1) alpha <- rep(alpha, nr)
# alpha <- alpha / sum(alpha)
# weights <- m * alpha * theta
# nc <- max(ncol(loc), ncol(scale), ncol(shape))
# loc <- matrix(loc, nr, nc)
# scale <- matrix(scale, nr, nc)
# shape <- matrix(shape, nr, nc)
# tau <- matrix(tau, nr, nc)
# theta <- matrix(theta, nr, nc)
# weights <- matrix(weights, nr, nc)
if (!(family %in% c("gev", "gpd", "bgev"))) stop("Invalid family")
loc <- as.matrix(loc)
scale <- as.matrix(scale)
shape <- as.matrix(shape)
nr <- max(nrow(loc), nrow(scale), nrow(shape))
nc <- max(ncol(loc), ncol(scale), ncol(shape))
loc <- matrix(loc, nr, nc)
scale <- matrix(scale, nr, nc)
shape <- matrix(shape, nr, nc)
alpha <- matrix(alpha, nr, nc)
alpha <- t(t(alpha)  / colSums(alpha))
theta <- matrix(theta, nr, nc)
weights <- m * alpha * theta
tau <- matrix(tau, nr, nc)
if (method == "uniroot") {
  out <- numeric(nc)
  for (i in seq_len(nc)) {
    out[i] <- .rlvec(p, loc[,i], scale[,i], shape[,i], m, nr, weights[,i], family, 
                     tau[,i], theta[,i], bgev.args, start)
  }
} else {
  out <- .rlvec2(p, loc, scale, shape, m, nr, weights, family, tau, theta, 
                 bgev.args, start)
}
out
}

.rlvec <- function(p, loc, scale, shape, m, n, weights, family, tau, theta, 
                   bgev.args, start = NULL) {
pa <- pb <- alpha <- beta <- 1
if (is.null(start)) {
  if (family == "gpd") {
    start <- range(.qgpd(p, loc, scale, shape, 1 - tau, theta, m))
  } else {
    if (family == "gev") {
      start <- range(.qgev(p^(1/n), loc, scale, shape))
    } else {
      pa <- bgev.args$pa
      pb <- bgev.args$pb
      alpha <- bgev.args$alpha
      beta <- bgev.args$beta
      start <- range(qbgev(p^(1/n), loc, scale, shape, pa, pb, alpha, beta))
    }
  }
  dstart <- diff(start)
  nullstart <- FALSE
} else {
  nullstart <- TRUE
}
while(.rlroot(start[1], loc, scale, shape, tau, weights, p, TRUE, family,
              pa, pb, alpha, beta) > 0) {
    if (nullstart) stop("Invalid `start' values given")
    start[1] <- start[1] - .2 * dstart
}
while(.rlroot(start[2], loc, scale, shape, tau, weights, p, TRUE, family,
              pa, pb, alpha, beta) < 0) {
    if (nullstart) stop("Invalid `start' values given")
    start[2] <- start[2] + .2 * dstart
}
opt <- uniroot(.rlroot, start, loc = loc, scale = scale, shape = shape, tau = tau, 
               weights = weights, p = p, TRUE, family, pa = pa, pb = pb, 
               alpha = alpha, beta = beta)
opt$root
}

.pmax <- function(x, loc, scale, shape, tau, weights, log, family, 
                  pa, pb, alpha, beta) {
if (family == "gpd") {
    out <- .pgpd(x, loc, scale, shape, tau, FALSE, TRUE)
} else {
  if (family == "gev") {
    out <- .pgev(x, loc, scale, shape, FALSE, TRUE)
  } else {
    out <- pbgev(x, loc, scale, shape, pa, pb, alpha, beta, TRUE, TRUE)
  }
}
out <- sum(weights * out)
if (!log) 
  out <- exp(out)
out
}

.rlroot <- function(x, loc, scale, shape, tau, weights, p, log, family,
                    pa, pb, alpha, beta) {
if (log) 
  p <- log(p)
.pmax(x, loc, scale, shape, tau, weights, log, family, pa, pb, alpha, beta) - p
}

.rlvec2 <- function(p, loc, scale, shape, m, n, weights, family, tau, theta, 
                    bgev.args, 
                    start = NULL, tol = 1e-5, eps = 1e-4, maxit = 10) {
  pa <- pb <- alpha <- beta <- 1
  if (is.null(start)) {
    if (family == "gpd") {
      start <- .qgpd(p, loc, scale, shape, 1 - tau, theta, m)
    } else {
      if (family == "gev") {
        start <- .qgev(p^(1/n), loc, scale, shape)
      } else {
        pa <- bgev.args$pa
        pb <- bgev.args$pb
        alpha <- bgev.args$alpha
        beta <- bgev.args$beta
        start <- qbgev(p^(1/n), loc, scale, shape, pa, pb, alpha, beta)
      }
    }
    dstart <- diff(start)
    nullstart <- FALSE
  } else {
    nullstart <- TRUE
  }
  start_obj <- t(apply(start, 1, function(x) .rlroot2(x, loc, scale, shape, tau, 
                                                      weights, p, TRUE, family,
                                                      pa, pb, alpha, beta)))
  rl0 <- start[cbind(apply(abs(start_obj), 2, which.min), seq_len(ncol(loc)))]
  f0 <- .rlroot2(rl0, loc, scale, shape, tau, weights, p, TRUE, family,
                 pa, pb, alpha, beta)
  it <- 1
  while(any(abs(f0) > tol) & it < maxit) {
    f1 <- .rlroot2_d1(rl0, loc, scale, shape, tau, weights, p, TRUE, family, 
                      pa, pb, alpha, beta, eps)
    rl0 <- rl0 - f0 / f1
    f0 <- .rlroot2(rl0, loc, scale, shape, tau, weights, p, TRUE, family,
                   pa, pb, alpha, beta)
    it <- it + 1
  }
  rl0
}

.pmax2 <- function(x, loc, scale, shape, tau, weights, log, family, 
                   pa, pb, alpha, beta) {
  x <- matrix(x, nrow(loc), ncol(loc), byrow = TRUE)
  if (family == "gpd") {
    out <- .pgpd(x, loc, scale, shape, tau, FALSE, TRUE)
  } else {
    if (family == "gev") {
      out <- .pgev(x, loc, scale, shape, FALSE, TRUE)
    } else {
      out <- pbgev(x, loc, scale, shape, pa, pb, alpha, beta, TRUE, TRUE)
    }
  }
  out <- colSums(as.matrix(weights * out))
  if (!log) 
    out <- exp(out)
  out
}

.rlroot2 <- function(x, loc, scale, shape, tau, weights, p, log, family,
                     pa, pb, alpha, beta) {
  if (log) 
    p <- log(p)
  .pmax2(x, loc, scale, shape, tau, weights, log, family, pa, pb, alpha, beta) - p
}

.rlroot2_d1 <- function(x, loc, scale, shape, tau, weights, p, log, family, 
                        pa, pb, alpha, beta, eps) {
  eps <- pmax(1, abs(x)) * eps
  fh <- .rlroot2(x + eps, loc, scale, shape, tau, weights, p, log, family,
                 pa, pb, alpha, beta)
  fl <- .rlroot2(x - eps, loc, scale, shape, tau, weights, p, log, family,
                 pa, pb, alpha, beta)
  .5 * (fh - fl) / eps
}
