#' The Laplace Distribution
#'
#' @description 
#' Density, distribution function, quantile function and random generation for 
#' the Laplace distribution with location \code{location} and scale \code{scale}.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param location location parameter (mean).
#' @param scale scale parameter (must be positive).
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], 
#'   otherwise, P[X > x].
#'
#' @name Laplace
NULL

#' @rdname Laplace
#' @export
dlaplace <- function(x, location = 0, scale = 1, log = FALSE) {
  if (any(scale <= 0)) stop("scale must be positive")
  d <- -abs(x - location) / scale - log(2 * scale)
  if (log) return(d) else return(exp(d))
}

#' @rdname Laplace
#' @export
plaplace <- function(q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (any(scale <= 0)) stop("scale must be positive")
  z <- (q - location) / scale
  p <- ifelse(z <= 0, exp(z) / 2, 1 - exp(-z) / 2)
  if (!lower.tail) p <- 1 - p
  if (log.p) return(log(p)) else return(p)
}

#' @rdname Laplace
#' @export
qlaplace <- function(p, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE) {
  if (any(scale <= 0)) stop("scale must be positive")
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  out <- ifelse(p <= .5, log(2 * p), -log(2 * (1 - p)))
  location + scale * out
}

#' @rdname Laplace
#' @export
rlaplace <- function(n, location = 0, scale = 1) {
  if (any(scale <= 0)) stop("scale must be positive")
  u <- runif(n, -0.5, 0.5)
  std_laplace <- sign(u) * -log(1 - 2 * abs(u))
  location + scale * std_laplace
}
