## Exponential negative log-likelihood functions

.exp.d0 <- function(pars, likdata) {
ny <- ncol(likdata$y)
out <- lapply(seq_len(ny), function(i) expd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$y[,i], likdata$dupid, likdata$duplicate))
out <- Reduce("+", out) / ny
out
}

.exp.d12 <- function(pars, likdata, sandwich = FALSE) {
ny <- ncol(likdata$y)
out <- lapply(seq_len(ny), function(i) expd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$y[,i], likdata$dupid, likdata$duplicate))
out <- Reduce("+", out) / ny
out
}

.exp.d34 <- function(pars, likdata) {
ny <- ncol(likdata$y)
out <- lapply(seq_len(ny), function(i) expd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$y[,i], likdata$dupid, likdata$duplicate))
out <- Reduce("+", out) / ny
out
}

.expfns <- list(d0=.exp.d0, d120=.exp.d12, d340=.exp.d34)
