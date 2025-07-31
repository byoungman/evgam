# The conditional extreme value model

.condexagg.d0 <- function(pars, likdata) {
  if (is.null(likdata$args$C))
    likdata$args$C <- .01
  if (is.null(likdata$args$epsilon))
    likdata$args$epsilon <- .01
  likdata$y <- as.matrix(likdata$y)
  likdata$args$x <- as.matrix(likdata$args$x)
  # nhere <- rowSums(is.finite(likdata$y))  
  # out <- condexaggd0(split(pars, likdata$idpars), 
  #                     likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], 
  #                     likdata$X[[4]], likdata$X[[5]], likdata$X[[6]], 
  #                     likdata$y, likdata$args$x, likdata$args$weights, 
  #                     likdata$dupid, likdata$duplicate, nhere,
  #                    likdata$args$C, likdata$args$epsilon)
  if (is.null(likdata$args$check)) {
    likdata$args$check <- rep(FALSE, nrow(likdata$X[[1]]))
    constrained <- FALSE
  } else {
    constrained <- TRUE
  }
  X_data <- lapply(likdata$X, function(x) x[!likdata$args$check, , drop = FALSE])
  y_data <- likdata$y[!likdata$args$check, , drop = FALSE]
  likdata$args$weights <- as.matrix(likdata$args$weights[!likdata$args$check, , drop = FALSE])
  nhere <- rowSums(is.finite(y_data))
  likdata$args$x <- as.matrix(likdata$args$x[!likdata$args$check, , drop = FALSE])
  out <- condexaggd0(split(pars, likdata$idpars),
                      X_data[[1]], X_data[[2]], X_data[[3]],
                      X_data[[4]], X_data[[5]], X_data[[6]],
                      y_data, likdata$args$x, likdata$args$weights,
                      likdata$dupid, likdata$duplicate, nhere,
                     likdata$args$C, likdata$args$epsilon)
  if (constrained) {
    X_constr <- lapply(likdata$X, function(x) x[likdata$args$check, , drop = FALSE])
    current <- mapply('%*%', X_constr, split(pars, likdata$idpars))
    out <- out + .parpen.d0(current, likdata)
  }
  out
}

.condexagg.d12 <- function(pars, likdata) {
  if (is.null(likdata$args$C))
    likdata$args$C <- .01
  if (is.null(likdata$args$epsilon))
    likdata$args$epsilon <- .01
  likdata$y <- as.matrix(likdata$y)
  likdata$args$x <- as.matrix(likdata$args$x)
  nhere <- rowSums(is.finite(likdata$y))  
  # out <- condexaggd12(split(pars, likdata$idpars), 
  #                     likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], 
  #                     likdata$X[[4]], likdata$X[[5]], likdata$X[[6]], 
  #                     likdata$y, likdata$args$x, likdata$args$weights, 
  #                     likdata$dupid, likdata$duplicate, nhere,
  #                     likdata$args$C, likdata$args$epsilon)
  if (is.null(likdata$args$check)) {
    likdata$args$check <- rep(FALSE, nrow(likdata$X[[1]]))
    constrained <- FALSE
  } else {
    constrained <- TRUE
  }
  X_data <- lapply(likdata$X, function(x) x[!likdata$args$check, , drop = FALSE])
  y_data <- likdata$y[!likdata$args$check, , drop = FALSE]
  likdata$args$weights <- as.matrix(likdata$args$weights[!likdata$args$check, , drop = FALSE])
  nhere <- rowSums(is.finite(y_data))
  likdata$args$x <- as.matrix(likdata$args$x[!likdata$args$check, , drop = FALSE])
  out <- condexaggd12(split(pars, likdata$idpars), 
                      X_data[[1]], X_data[[2]], X_data[[3]], 
                      X_data[[4]], X_data[[5]], X_data[[6]], 
                      y_data, likdata$args$x, likdata$args$weights, 
                      likdata$dupid, likdata$duplicate, nhere,
                      likdata$args$C, likdata$args$epsilon)
  if (constrained) {
    X_constr <- lapply(likdata$X, function(x) x[likdata$args$check, , drop = FALSE])
    current <- mapply('%*%', X_constr, split(pars, likdata$idpars))
    out <- rbind(out, .parpen.d12(current, likdata))
  }
  out
}

.condexagg.d34 <- function(pars, likdata) {
  if (is.null(likdata$args$C))
    likdata$args$C <- .01
  if (is.null(likdata$args$epsilon))
    likdata$args$epsilon <- .01
  likdata$y <- as.matrix(likdata$y)
  likdata$args$x <- as.matrix(likdata$args$x)
  # nhere <- rowSums(is.finite(likdata$y))  
  # out <- condexaggd34(split(pars, likdata$idpars), 
  #                  likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], 
  #                  likdata$X[[4]], likdata$X[[5]], likdata$X[[6]], 
  #                  likdata$y, likdata$args$x, likdata$args$weights, 
  #                  likdata$dupid, likdata$duplicate, nhere,
  #                  likdata$args$C, likdata$args$epsilon)
  if (is.null(likdata$args$check)) {
    likdata$args$check <- rep(FALSE, nrow(likdata$X[[1]]))
    constrained <- FALSE
  } else {
    constrained <- TRUE
  }
  X_data <- lapply(likdata$X, function(x) x[!likdata$args$check, , drop = FALSE])
  y_data <- likdata$y[!likdata$args$check, , drop = FALSE]
  likdata$args$weights <- as.matrix(likdata$args$weights[!likdata$args$check, , drop = FALSE])
  nhere <- rowSums(is.finite(y_data))
  likdata$args$x <- as.matrix(likdata$args$x[!likdata$args$check, , drop = FALSE])
  out <- condexaggd34(split(pars, likdata$idpars), 
                      X_data[[1]], X_data[[2]], X_data[[3]], 
                      X_data[[4]], X_data[[5]], X_data[[6]], 
                      y_data, likdata$args$x, likdata$args$weights, 
                      likdata$dupid, likdata$duplicate, nhere,
                      likdata$args$C, likdata$args$epsilon)
  if (constrained) {
    out <- rbind(out, 0)
  }
  out
}

.condexaggfns <- list(d0 = .condexagg.d0, d120 = .condexagg.d12, d340 = .condexagg.d34)

.condexagg_unlink <- list(function(x) 2 / (1 + exp(-x)) - 1, 
                       function(x) 1 / (1 + exp(-x)), 
                       function(x) x, 
                       function(x) exp(x),
                       function(x) exp(x),
                       function(x) exp(x))
attr(.condexagg_unlink[[1]], "deriv") <- function(x) 2 * exp(-x)/(1 + exp(-x))^2
attr(.condexagg_unlink[[2]], "deriv") <- function(x) exp(-x)/(1 + exp(-x))^2
attr(.condexagg_unlink[[3]], "deriv") <- function(x) x
attr(.condexagg_unlink[[4]], "deriv") <- function(x) exp(x)
attr(.condexagg_unlink[[5]], "deriv") <- function(x) exp(x)
attr(.condexagg_unlink[[6]], "deriv") <- function(x) exp(x)

.condexaggfns$q <- NULL
.condexaggfns$unlink <- .condexagg_unlink

