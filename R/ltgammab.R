## Left truncated gamma negative log-likelihood functions

.ltgammab.d0 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  if (is.null(likdata$args$check)) {
    likdata$args$check <- rep(FALSE, nrow(likdata$X[[1]]))
    with_boundary <- FALSE
  } else {
    with_boundary <- TRUE
  }
  X_data <- lapply(likdata$X, function(x) x[!likdata$args$check, , drop = FALSE])
  y_data <- likdata$y[!likdata$args$check, , drop = FALSE]
  nhere <- rowSums(is.finite(y_data))
  left <- likdata$args$left
  if (length(left) == 1)
    left <- array(left, dim(y_data))
  alpha <- likdata$args$alpha
  if (length(alpha) == 1)
    alpha <- array(alpha, dim(y_data))
  out <- ltgammabd0(split(pars, likdata$idpars), X_data[[1]], y_data, likdata$dupid, likdata$duplicate, nhere, as.matrix(left), as.matrix(alpha))
  if (with_boundary) {
    ps <- split(pars, likdata$idpars)
    lograte <- likdata$X[[1]][likdata$args$check, , drop = FALSE] %*% ps[[1]]
    out <- out + sum(.ltgammab_bpen.d0(lograte, likdata))
  }
  if (!is.finite(out))
    out <- 1e20
  out
}

.ltgammab_bpen.d0 <- function(lograte, likdata) {
  rate <- exp(lograte)
  bdry <- likdata$args$boundary / rate[, 1]
  if (is.null(likdata$args$bpen))
    likdata$args$bpen <- 1
  pen <- matrix(0, length(rate), 4)
  up_diff <- (bdry - 1)^2
  idup <- apply(up_diff, 2, which.min)
  down_diff <- (bdry + 1)^2
  iddown <- apply(down_diff, 2, which.min)
  pen[cbind(c(idup, iddown), 1:4)] <- c(up_diff[cbind(idup, 1:2)], down_diff[cbind(iddown, 1:2)])
  likdata$args$bpen * rowSums(pen)
}  

.ltgammab.d12 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  if (is.null(likdata$args$check)) {
    likdata$args$check <- rep(FALSE, nrow(likdata$X[[1]]))
    with_boundary <- FALSE
  } else {
    with_boundary <- TRUE
  }
  X_data <- lapply(likdata$X, function(x) x[!likdata$args$check, , drop = FALSE])
  y_data <- likdata$y[!likdata$args$check, , drop = FALSE]
  nhere <- rowSums(is.finite(y_data))
  left <- likdata$args$left
  if (length(left) == 1)
    left <- array(left, dim(y_data))
  alpha <- likdata$args$alpha
  if (length(alpha) == 1)
    alpha <- array(alpha, dim(y_data))
  out <- ltgammabd12(split(pars, likdata$idpars), X_data[[1]], y_data, likdata$dupid, likdata$duplicate, nhere, as.matrix(left), as.matrix(alpha))
  if (with_boundary) {
    ps <- split(pars, likdata$idpars)
    lograte <- likdata$X[[1]][likdata$args$check, , drop = FALSE] %*% ps[[1]]
    h <- 1e-4
    f0 <- .ltgammab_bpen.d0(lograte, likdata)
    fhi <- .ltgammab_bpen.d0(lograte + h, likdata)
    flo <- .ltgammab_bpen.d0(lograte - h, likdata)
    out <- rbind(out, cbind((fhi - f0) / h, (fhi + flo - 2  * f0) / h^2))
  }
  out
}

.ltgammab.d34 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  if (is.null(likdata$args$check)) {
    likdata$args$check <- rep(FALSE, nrow(likdata$X[[1]]))
    with_boundary <- FALSE
  } else {
    with_boundary <- TRUE
  }
  X_data <- lapply(likdata$X, function(x) x[!likdata$args$check, , drop = FALSE])
  y_data <- likdata$y[!likdata$args$check, , drop = FALSE]
  nhere <- rowSums(is.finite(y_data))
  left <- likdata$args$left
  if (length(left) == 1)
    left <- array(left, dim(y_data))
  alpha <- likdata$args$alpha
  if (length(alpha) == 1)
    alpha <- array(alpha, dim(y_data))
  out <- ltgammabd34(split(pars, likdata$idpars), X_data[[1]], y_data, likdata$dupid, likdata$duplicate, nhere, as.matrix(left), as.matrix(alpha))
  if (with_boundary) {
    ps <- split(pars, likdata$idpars)
    lograte <- likdata$X[[1]][likdata$args$check, , drop = FALSE] %*% ps[[1]]
    h <- 1e-2
    fhi <- .ltgammab_bpen.d0(lograte + h, likdata)
    flo <- .ltgammab_bpen.d0(lograte - h, likdata)
    f2hi <- .ltgammab_bpen.d0(lograte + 2 * h, likdata)
    f2lo <- .ltgammab_bpen.d0(lograte - 2 * h, likdata)
    d3 <- (f2lo - 2 * flo + 2 * fhi - f2lo) / 2 / h^3
    out <- rbind(out, cbind(d3, 0))
  }
  out
}

# .ltgammabfns <- list(d0 = .ltgammab.d0, d120 = .ltgammab.d12, d340 = NULL)
.ltgammabfns <- list(d0 = .ltgammab.d0, d120 = .ltgammab.d12, d340 = .ltgammab.d34)

.ltgammab_unlink <- list(function(x) exp(x))
attr(.ltgammab_unlink[[1]], "deriv") <- .ltgammab_unlink[[1]]

.ltgammafns$q <- NULL
.ltgammabfns$unlink <- .ltgammab_unlink

