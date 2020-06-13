## Weibull negative log-likelihood functions

.weib.d0 <- function(pars, likdata) {
ny <- ncol(likdata$y)
out <- lapply(seq_len(ny), function(i) weibd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y[,i], likdata$dupid, likdata$duplicate))
out <- Reduce("+", out) / ny
out
}

.weib.d12 <- function(pars, likdata) {
ny <- ncol(likdata$y)
out <- lapply(seq_len(ny), function(i) weibd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y[,i], likdata$dupid, likdata$duplicate))
out <- Reduce("+", out) / ny
out
}

.weib.d34 <- function(pars, likdata) {
ny <- ncol(likdata$y)
out <- lapply(seq_len(ny), function(i) weibd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y[,i], likdata$dupid, likdata$duplicate))
out <- Reduce("+", out) / ny
out
}

.weibfns <- list(d0=.weib.d0, d120=.weib.d12, d340=.weib.d34)

.qweibull <- function(x, scale, shape, logscale=NULL, logshape=NULL) {
if (!missing(logscale)) scale <- exp(logscale)
if (!missing(logshape)) shape <- exp(logshape)
scale * ((-log(1 - x)) ^ (1/shape))
}

.pweibull <- function(x, scale, shape, logscale=NULL, logshape=NULL) {
if (!missing(logscale)) scale <- exp(logscale)
if (!missing(logshape)) shape <- exp(logshape)
1 - exp(-(x / scale)^shape)
}

.dqweibull <- function(x, lscale, lshape) {
.e1 <- -log(1 - x)
.e2 <- exp(lshape)
.e4 <- .e1^(1/.e2) * exp(lscale)
d1 <- .e4
d2 <--(.e4 * log(.e1)/.e2)
cbind(d1, d2)
}
