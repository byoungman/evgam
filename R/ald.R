## Asymmetric Laplace negative log-likelihood functions

.ald.d0 <- function(pars, likdata) {
# ny <- ncol(likdata$y)
# out <- lapply(seq_len(ny), function(i) aldd0(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], yvec=likdata$y[,i], tau=likdata$tau, C=likdata$C, dupid=likdata$dupid, dcate=likdata$duplicate))
# out <- Reduce("+", out) / ny
likdata$y <- as.matrix(likdata$y)
nhere <- rowSums(is.finite(likdata$y))  
out <- aldd0(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], ymat=likdata$y, tau=likdata$tau, C=likdata$C, dupid=likdata$dupid, dcate=likdata$duplicate, nhere = nhere)
if (!is.finite(out)) out <- 1e20
out
}

.ald.d12 <- function(pars, likdata, sandwich = FALSE) {
# ny <- ncol(likdata$y)
# out <- lapply(seq_len(ny), function(i) aldd12(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], yvec=likdata$y[,i], tau=likdata$tau, C=likdata$C, dupid=likdata$dupid, dcate=likdata$duplicate))
# out <- Reduce("+", out) / ny
likdata$y <- as.matrix(likdata$y)
nhere <- rowSums(is.finite(likdata$y))  
out <- aldd12(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], ymat=likdata$y, tau=likdata$tau, C=likdata$C, dupid=likdata$dupid, dcate=likdata$duplicate, nhere = nhere)
out
}

.ald.d34 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- aldd34(split(pars, likdata$idpars), X1=likdata$X[[1]], X2=likdata$X[[2]], ymat=likdata$y, tau=likdata$tau, C=likdata$C, dupid=likdata$dupid, dcate=likdata$duplicate, nhere = nhere)
  out
}

.aldfns <- list(d0=.ald.d0, d120=.ald.d12, d340=.ald.d34)

