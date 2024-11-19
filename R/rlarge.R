## r-largest generalised extreme value negative log-likelihood functions

.rlarge.d0 <- function(pars, likdata) {
likdata$y <- as.matrix(likdata$y)
out <- rlarged0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate)
out
}

.rlarge.d12 <- function(pars, likdata) {
likdata$y <- as.matrix(likdata$y)
out <- rlarged12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate)
out
}

.rlarge.d34 <- function(pars, likdata) {
likdata$y <- as.matrix(likdata$y)
out <- rlarged34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate)
out
}

.rlargefns <- list(d0 = .rlarge.d0, d120 = .rlarge.d12, d340 = .rlarge.d34)

.rlarge_unlink <- list(NULL, function(x) exp(x), function(x) 1.5 / (1 + exp(-x)) - 1.0)
attr(.rlarge_unlink[[2]], "deriv") <- .rlarge_unlink[[2]]
attr(.rlarge_unlink[[3]], "deriv") <- function(x) 1.5 * exp(-x)/(1 + exp(-x))^2

.rlargefns$q <- .qgev
.rlargefns$unlink <- .rlarge_unlink

