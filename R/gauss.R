## Gaussian negative log-likelihood functions

.gauss.d0 <- function(pars, likdata) {
ny <- ncol(likdata$y)
out <- lapply(seq_len(ny), function(i) gaussd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y[,i], likdata$dupid, likdata$duplicate))
out <- Reduce("+", out) / ny
out
}

.gauss.d12 <- function(pars, likdata) {
ny <- ncol(likdata$y)
out <- lapply(seq_len(ny), function(i) gaussd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y[,i], likdata$dupid, likdata$duplicate))
out <- Reduce("+", out) / ny
out
}

.gauss.d34 <- function(pars, likdata) {
ny <- ncol(likdata$y)
out <- lapply(seq_len(ny), function(i) gaussd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y[,i], likdata$dupid, likdata$duplicate))
out <- Reduce("+", out) / ny
out
}

.gaussfns <- list(d0=.gauss.d0, d120=.gauss.d12, d340=.gauss.d34)
