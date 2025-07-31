# The conditional extreme value model

.condex.d0 <- function(pars, likdata) {
  # likdata$y <- as.matrix(likdata$y)
  # likdata$args$x <- as.matrix(likdata$args$x)
  # nhere <- rowSums(is.finite(likdata$y))  
  # if (!likdata$sparse) {
  #   out <- condexd0(split(pars, likdata$idpars), 
  #                 likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], 
  #                 likdata$y, likdata$args$x, likdata$args$weights, 
  #                 likdata$dupid, likdata$duplicate, nhere)
  # } else {
  #   out <- condexspd0(split(pars, likdata$idpars), 
  #                   likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], 
  #                   likdata$y, likdata$args$x, likdata$args$weights, 
  #                   likdata$dupid, likdata$duplicate, nhere)
  # }
  if (is.null(likdata$args$check)) {
    likdata$args$check <- rep(FALSE, nrow(likdata$X[[1]]))
    constrained <- FALSE
  } else {
    constrained <- TRUE
  }
  X_data <- lapply(likdata$X, function(x) x[!likdata$args$check, , drop = FALSE])
  y_data <- likdata$y[!likdata$args$check, , drop = FALSE]
  nhere <- rowSums(is.finite(y_data))
  likdata$args$x <- as.matrix(likdata$args$x[!likdata$args$check, , drop = FALSE])
  likdata$args$weights <- as.matrix(likdata$args$weights[!likdata$args$check, , drop = FALSE])
  if (!likdata$sparse) {
    out <- condexd0(split(pars, likdata$idpars), 
                    X_data[[1]], X_data[[2]], X_data[[3]], X_data[[4]], 
                    y_data, likdata$args$x, likdata$args$weights, 
                    likdata$dupid, likdata$duplicate, nhere)
  } else {
    out <- condexspd0(split(pars, likdata$idpars), 
                      X_data[[1]], X_data[[2]], X_data[[3]], X_data[[4]], 
                      y_data, likdata$args$x, likdata$args$weights, 
                      likdata$dupid, likdata$duplicate, nhere)
  }
  if (constrained) {
    X_constr <- lapply(likdata$X, function(x) x[likdata$args$check, , drop = FALSE])
    current <- mapply('%*%', X_constr, split(pars, likdata$idpars))
    out <- out + .parpen.d0(current, likdata)
  }
  out
}

.parpen.d0 <- function(current, likdata) {
  if (is.null(likdata$args$pen))
    likdata$args$pen <- 1e2
  pen <- likdata$args$pen * sum((current - likdata$args$target)^2, na.rm = TRUE)
  pen
}

.condex.d12 <- function(pars, likdata) {
  # likdata$y <- as.matrix(likdata$y)
  # likdata$args$x <- as.matrix(likdata$args$x)
  # nhere <- rowSums(is.finite(likdata$y))  
  # if (!likdata$sparse) {
  #   out <- condexd12(split(pars, likdata$idpars), 
  #                  likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], 
  #                  likdata$y, likdata$args$x, likdata$args$weights, 
  #                  likdata$dupid, likdata$duplicate, nhere)
  # } else {
  #   out <- condexspd12(split(pars, likdata$idpars), 
  #                    likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], 
  #                    likdata$y, likdata$args$x, likdata$args$weights, 
  #                    likdata$dupid, likdata$duplicate, nhere)
  #   
  # }
  if (is.null(likdata$args$check)) {
    likdata$args$check <- rep(FALSE, nrow(likdata$X[[1]]))
    constrained <- FALSE
  } else {
    constrained <- TRUE
  }
  X_data <- lapply(likdata$X, function(x) x[!likdata$args$check, , drop = FALSE])
  y_data <- likdata$y[!likdata$args$check, , drop = FALSE]
  nhere <- rowSums(is.finite(y_data))
  likdata$args$x <- as.matrix(likdata$args$x[!likdata$args$check, , drop = FALSE])
  likdata$args$weights <- as.matrix(likdata$args$weights[!likdata$args$check, , drop = FALSE])
  if (!likdata$sparse) {
    out <- condexd12(split(pars, likdata$idpars), 
                     X_data[[1]], X_data[[2]], X_data[[3]], X_data[[4]], 
                     y_data, likdata$args$x, likdata$args$weights, 
                     likdata$dupid, likdata$duplicate, nhere)
  } else {
    out <- condexspd12(split(pars, likdata$idpars), 
                       X_data[[1]], X_data[[2]], X_data[[3]], X_data[[4]], 
                       y_data, likdata$args$x, likdata$args$weights, 
                       likdata$dupid, likdata$duplicate, nhere)
    
  }
  if (constrained) {
    X_constr <- lapply(likdata$X, function(x) x[likdata$args$check, , drop = FALSE])
    current <- mapply('%*%', X_constr, split(pars, likdata$idpars))
    out <- rbind(out, .parpen.d12(current, likdata))
  }
out
}

.parpen.d12 <- function(current, likdata, h = 1e-4) {
  np <- length(current)
  i_ok <- which(is.finite(likdata$args$target))
  # f0 <- .parpen.d0(current, likdata)
  # g0 <- fh <- fl <- numeric(np)
  # h0 <- numeric(.5 * np * (np + 1))
  # it <- 1
  # for (i in i_ok) {
  #   ph <- replace(current, i, current[i] + h)
  #   fh[i] <- .parpen.d0(ph, likdata)
  #   pl <- replace(current, i, current[i] - h)
  #   fl[i] <- .parpen.d0(pl, likdata)
  #   for (j in i_ok[i_ok <= i]) {
  #     phh <- replace(ph, j, current[j] + h)
  #     phl <- replace(ph, j, current[j] - h)
  #     plh <- replace(pl, j, current[j] + h)
  #     pll <- replace(pl, j, current[j] - h)
  #     fhh <- .parpen.d0(phh, likdata)
  #     fhl <- .parpen.d0(phl, likdata)
  #     flh <- .parpen.d0(plh, likdata)
  #     fll <- .parpen.d0(pll, likdata)
  #     h0[it] <- (fhh - fh[i] - fh[j] + 2 * f0 - fl[i] - fl[j] + fll)
  #     it <- it + 1
  #   }
  # }
  # g0 <- (fh - fl) / 2
  # h0 <- h0 / h
  # c(g0, h0) / h
  if (is.null(likdata$args$pen))
    likdata$args$pen <- 1e2
  g0 <- 2 * (current - likdata$args$target)
  H0 <- matrix(0, np, np)
  H0[i_ok, i_ok] <- 2
  g0[-i_ok] <- 0
  likdata$args$pen * c(g0, H0[!upper.tri(H0)])
}

.condex.d34 <- function(pars, likdata) {
#   likdata$y <- as.matrix(likdata$y)
#   likdata$args$x <- as.matrix(likdata$args$x)
#   nhere <- rowSums(is.finite(likdata$y))  
#   if (!likdata$sparse) {
#     out <- condexd34(split(pars, likdata$idpars), 
#                    likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], 
#                    likdata$y, likdata$args$x, likdata$args$weights, 
#                    likdata$dupid, likdata$duplicate, nhere)
# } else {
#   out <- condexspd34(split(pars, likdata$idpars), 
#                    likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], 
#                    likdata$y, likdata$args$x, likdata$args$weights, 
#                    likdata$dupid, likdata$duplicate, nhere)
#   
# }
  if (is.null(likdata$args$check)) {
    likdata$args$check <- rep(FALSE, nrow(likdata$X[[1]]))
    constrained <- FALSE
  } else {
    constrained <- TRUE
  }
  X_data <- lapply(likdata$X, function(x) x[!likdata$args$check, , drop = FALSE])
  y_data <- likdata$y[!likdata$args$check, , drop = FALSE]
  nhere <- rowSums(is.finite(y_data))
  likdata$args$x <- as.matrix(likdata$args$x[!likdata$args$check, , drop = FALSE])
  if (!likdata$sparse) {
    out <- condexd34(split(pars, likdata$idpars), 
                     X_data[[1]], X_data[[2]], X_data[[3]], X_data[[4]], 
                     y_data, likdata$args$x, likdata$args$weights, 
                     likdata$dupid, likdata$duplicate, nhere)
  } else {
    out <- condexspd34(split(pars, likdata$idpars), 
                       X_data[[1]], X_data[[2]], X_data[[3]], X_data[[4]], 
                       y_data, likdata$args$x, likdata$args$weights, 
                       likdata$dupid, likdata$duplicate, nhere)
    
  }
  if (constrained) {
    out <- rbind(out, 0)
  }
  out
}

.condexfns <- list(d0 = .condex.d0, d120 = .condex.d12, d340 = .condex.d34)

.condex_unlink <- list(function(x) 2 / (1 + exp(-x)) - 1, function(x) 1 / (1 + exp(-x)), function(x) x, function(x) exp(x))
attr(.condex_unlink[[1]], "deriv") <- function(x) 2 * exp(-x)/(1 + exp(-x))^2
attr(.condex_unlink[[2]], "deriv") <- function(x) exp(-x)/(1 + exp(-x))^2
attr(.condex_unlink[[3]], "deriv") <- function(x) x
attr(.condex_unlink[[4]], "deriv") <- function(x) exp(x)

.condexfns$q <- NULL
.condexfns$unlink <- .condex_unlink

