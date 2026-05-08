#' Convert Cartesian coordinates to Polar coordinates
#'
#' @param x A matrix with exactly 2 rows representing the x and y coordinates.
#' @param norm A numeric value specifying the L-p norm for calculating the radius. Defaults to 2 (Euclidean).
#'
#' @return A list containing:
#' \item{r}{The calculated radii (distances from origin).}
#' \item{phi}{The angles in radians, adjusted to the range [0, 2*pi].}
#' @export
polarise <- function(x, norm = 2) {
  if (nrow(x) != 2)
    stop('This function just works in two dimensions.')
  radius <- colSums(abs(x)^norm)^(1 / norm)
  phi <- atan2(x[2, ], x[1, ])
  phi <- matrix(ifelse(phi < 0, phi + 2 * pi, phi), nrow = 1)
  list(r = radius, phi = phi)
}

#' Fit Geometric Extreme Value Generalized Additive Models
#'
#' Fits a two-stage model consisting of a threshold model (using Asymmetric Laplace Distribution) 
#' and an exceedance model (using a Gamma-type distribution).
#'
#' @param formula A list of two formulas named \code{threshold} and \code{excess}. 
#'   If missing, defaults to smooth angular terms.
#' @param data A data frame containing the variables in the model, including \code{radius} and \code{angle}.
#' @param args A list of additional arguments. \code{tau} (quantile) defaults to 0.8.
#' @param trace Integer vector of length 1 or 2. Control verbose output for both fitting stages.
#' @param knots A list containing knot specifications for \code{threshold} and \code{excess}.
#' @param norm The L-p norm used in the geometric calculations. Defaults to 2.
#'
#' @examples
#' 
#' n <- 1e3
#' gauss_sim_data <- rmvnorm(n, numeric(2), matrix(c(1, .8, .8, 1), 2))
#' laplace_sim_data <- qlaplace(pnorm(gauss_sim_data))
#' polar1 <- polarise(t(laplace_sim_data), norm = 2)
#' data1 <- data.frame(radius = polar1$r, angle = polar1$phi[1, ])
#' fit1 <- geoevgam(data = data1)
#' 
#' \donttest{
#' # example using L_1 norm
#' polar2 <- polarise(t(laplace_sim_data), norm = 1)
#' data2 <- data.frame(radius = polar2$r, angle = polar2$phi[1, ])
#' fit2 <- geoevgam(data = data2) 
#' 
#' # example using different threshold and user-specified splines and knots
#' angles <- seq(0, 2 * pi, by = pi / 4)
#' angles <- sort(c(angles, outer(c(pi / 4, pi + pi / 4), c(-pi / 16, pi / 16), '+')))
#' fmla <- knts <- list()
#' knts$threshold <- knts$excess <- list(angle = angles)
#' fmla$threshold <- fmla$excess <- radius ~ s(angle, bs = 'cp', k = length(angles) - 1)
#' fit3 <- geoevgam(data = data1, formula = fmla, knots = knts, args = list(tau = .9))
#' }
#'
#' @return An object of class \code{geoevgam} containing the fitted models and the norm used.
#' @export
geoevgam <- function(formula, data, args = list(), trace = 0, 
                     knots = list(threshold = NULL, excess = NULL), norm = 2) {
  if (missing(formula)) {
    formula <- list()
    formula$threshold <- radius ~ s(angle, bs = 'cc', k = 12)
    formula$excess <-  radius ~ s(angle, bs = 'cc', k = 12)
    if (is.null(knots$threshold))
      knots$threshold <- list(angle = c(0, 2 * pi))
    if (is.null(knots$excess))
      knots$excess <- list(angle = c(0, 2 * pi))
  }
  if (is.null(args$tau))
    args$tau <- 0.8
  out <- list()
  if (length(trace) == 1)
    trace <- rep(trace, 2)
  if (trace[1] > 0)
    cat('Fitting threshold model...')
  out$threshold <- evgam(formula$threshold, data = data, family = 'ald', 
                         args = list(tau = args$tau), 
                         knots = knots$threshold, trace = trace[1])
  if (trace[1] > 0) {
    cat('done.\n')
    if (trace[2] > 0)
      cat('\n')
  }
  data$fitted_threshold <- fitted(out$threshold)$location
  data_exc <- data[data$radius > data$fitted_threshold, ]
  if (trace[2] > 0)
    cat('Fitting exceedance model...')
  out$excess <- evgam(formula$excess, data = data_exc, family = 'ltgammab', 
                      args = list(left = data_exc$fitted_threshold, alpha = 2), 
                      knots = knots$excess, trace = trace[2])
  out$norm <- norm
  if (trace[1] > 0)
    cat('done.\n')
  class(out) <- 'geoevgam'
  out
}

#' Plot results of a geoevgam fit
#'
#' @param x An object of class \code{geoevgam}.
#' @param type Character vector specifying plot types: "threshold", "excess", "gauge", or "qqplot".
#' @param view Character specifying the coordinate system for gauge plots: "polar" or "cartesian".
#' @param norm Numeric; the norm to use for the gauge plot. Defaults to the norm used in \code{object}.
#' @param margins Character; margin type for QQ-plots. Defaults to "original".
#' @param ... Additional arguments passed to the underlying plot functions.
#' 
#' @examples
#' 
#' n <- 1e3
#' gauss_sim_data <- rmvnorm(n, numeric(2), matrix(c(1, .8, .8, 1), 2))
#' laplace_sim_data <- qlaplace(pnorm(gauss_sim_data))
#' polar1 <- polarise(t(laplace_sim_data), norm = 2)
#' data1 <- data.frame(radius = polar1$r, angle = polar1$phi[1, ])
#' fit1 <- geoevgam(data = data1)
#' plot(fit1)
#' plot(fit1, type = 'gauge', col = 2, lwd = 2, view = 'cartesian')
#' par(mfrow = c(2, 2))
#' plot(fit1, type = 'qqplot')
#' plot(fit1, type = 'qqplot2')
#' plot(fit1, type = 'qqplot', margins = 'uniform')
#' plot(fit1, type = 'qqplot2', margins = 'uniform')
#' 
#' \donttest{
#' # example using L_1 norm
#' polar2 <- polarise(t(laplace_sim_data), norm = 1)
#' data2 <- data.frame(radius = polar2$r, angle = polar2$phi[1, ])
#' fit2 <- geoevgam(data = data2, norm = 1) 
#' # passing norm ensures plotted gauge is correct
#' plot(fit2, type = 'gauge', col = 2, lwd = 2, view = 'cartesian')
#' }
#'
#' @method plot geoevgam
#' @export
plot.geoevgam <- function(x, type = c('threshold', 'excess'), 
                          view = c('polar', 'cartesian'), norm, 
                          margins = 'original', ...) {
  if ('threshold' %in% type) {
    plot(x$threshold, onepage = length(type) == 1, ...)
    title(main = 'Threshold', line = 3)
  }
  if ('excess' %in% type) {
    plot(x$excess, onepage = length(type) == 1, ...)
    title(main = 'Excess', line = 3)
  }
  if ('gauge' %in% type) {
    angles <- seq(0, 2 * pi, length = 100)
    newdata <- data.frame(angle = angles)
    rate <- predict(x$excess, newdata = newdata, type = 'response', std.err = TRUE)$rate
    if ('polar' %in% view)
      plot(angles, rate, type = 'l')
    if ('cartesian' %in% view) {
      unit_circle <- cbind(cos(angles), sin(angles))
      if (missing(norm)) {
        norm <- x$norm
      }
      if (norm != 2) {
        radii <- (abs(unit_circle[, 1])^norm + abs(unit_circle[, 2])^norm)^(1 / norm)
        unit_circle <- unit_circle / radii
      }
      limit_set <- as.data.frame(unit_circle / rate)
      names(limit_set) <- paste('x', 1:2, sep = '_')
      plot(limit_set, type = 'l', xlim = range(c(-1, 1, limit_set[, 1])), 
           ylim = range(c(-1, 1, limit_set[, 2])), ...)
      rect(-1, -1, 1, 1, lty = 2)
    }
  }
  type2 <- tolower(substr(type, 1, 2))
  if ('qqplot' %in% type | 'qqplot2' %in% type)
    predict(x$excess, type = type, margins = margins)
}
