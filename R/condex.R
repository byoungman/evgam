# The conditional extreme value model

.condex.d0 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  likdata$args$x <- as.matrix(likdata$args$x)
  nhere <- rowSums(is.finite(likdata$y))  
  if (!likdata$sparse) {
    out <- condexd0(split(pars, likdata$idpars), 
                  likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], 
                  likdata$y, likdata$args$x, likdata$args$weights, 
                  likdata$dupid, likdata$duplicate, nhere)
  } else {
    out <- condexspd0(split(pars, likdata$idpars), 
                    likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], 
                    likdata$y, likdata$args$x, likdata$args$weights, 
                    likdata$dupid, likdata$duplicate, nhere)
  }
  out
}

.condex.d12 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  likdata$args$x <- as.matrix(likdata$args$x)
  nhere <- rowSums(is.finite(likdata$y))  
  if (!likdata$sparse) {
    out <- condexd12(split(pars, likdata$idpars), 
                   likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], 
                   likdata$y, likdata$args$x, likdata$args$weights, 
                   likdata$dupid, likdata$duplicate, nhere)
  } else {
    out <- condexspd12(split(pars, likdata$idpars), 
                     likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], 
                     likdata$y, likdata$args$x, likdata$args$weights, 
                     likdata$dupid, likdata$duplicate, nhere)
    
}
out
}

.condex.d34 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  likdata$args$x <- as.matrix(likdata$args$x)
  nhere <- rowSums(is.finite(likdata$y))  
  if (!likdata$sparse) {
    out <- condexd34(split(pars, likdata$idpars), 
                   likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], 
                   likdata$y, likdata$args$x, likdata$args$weights, 
                   likdata$dupid, likdata$duplicate, nhere)
} else {
  out <- condexspd34(split(pars, likdata$idpars), 
                   likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], 
                   likdata$y, likdata$args$x, likdata$args$weights, 
                   likdata$dupid, likdata$duplicate, nhere)
  
}
out
}

.condexfns <- list(d0 = .condex.d0, d120 = .condex.d12, d340 = .condex.d34)

.condex_unlink <- list(function(x) 2 / (1 + exp(-x)) - 1, function(x) 1 / (1 + exp(-x)), function(x) x, function(x) exp(x))
attr(.condex_unlink[[2]], "deriv") <- function(x) exp(-x)/(1 + exp(-x))^2
attr(.condex_unlink[[3]], "deriv") <- function(x) x
attr(.condex_unlink[[4]], "deriv") <- function(x) exp(x)

.condexfns$q <- NULL
.condexfns$unlink <- .condex_unlink

