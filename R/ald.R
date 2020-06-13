## Asymmetric Laplace negative log-likelihood functions

.ald.d0 <- function(pars, likdata) {
ny <- ncol(likdata$y)
out <- lapply(seq_len(ny), function(i) aldd0(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], yvec=likdata$y[,i], tau=likdata$tau, C=likdata$C, dupid=likdata$dupid, dcate=likdata$duplicate))
out <- Reduce("+", out) / ny
if (!is.finite(out)) out <- 1e20
out
}

.ald.d12 <- function(pars, likdata) {
ny <- ncol(likdata$y)
out <- lapply(seq_len(ny), function(i) aldd12(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], yvec=likdata$y[,i], tau=likdata$tau, C=likdata$C, dupid=likdata$dupid, dcate=likdata$duplicate))
out <- Reduce("+", out) / ny
out
}

.ald.d34 <- function(pars, likdata) {
ny <- ncol(likdata$y)
out <- lapply(seq_len(ny), function(i) aldd34(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], yvec=likdata$y[,i], tau=likdata$tau, C=likdata$C, dupid=likdata$dupid, dcate=likdata$duplicate))
out <- Reduce("+", out) / ny
out
}

.aldfns <- list(d0=.ald.d0, d120=.ald.d12, d340=.ald.d34)

