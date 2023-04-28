## Aggregated GEV functions

.gevagg.d0 <- function(pars, likdata) {
pars <- split(pars, likdata$idpars)
pmat <- mapply("%*%", likdata$X, pars)
# identify which are aggregated
agg <- likdata$agg$agg
# point-referenced data
out <- ldgev(likdata$y[!agg, 1], pmat[!agg, 1], pmat[!agg, 2], pmat[!agg, 3])
if (!is.finite(out))
  return(1e20)
# # aggregated data
pdfr0 <- mapply("%*%", likdata$agg$X[1:2], pars[1:2])
pdfr0[,2] <- exp(pdfr0[,2])
pdfr <- ragged_mean_mat(pdfr0, likdata$agg$nxy)
out <- out + evgam:::ldgevagg(likdata$y[agg, 1], pdfr[,1], log(pdfr[,2]), pmat[agg, 3], pmat[agg, 4])
if (!is.finite(out))
  return(1e20)
pmin(out, 1e20)
}

.gevagg.d12 <- function(pars, likdata, sandwich = FALSE) {
pars <- split(pars, likdata$idpars)
pmat <- mapply("%*%", likdata$X, pars)
# identify which are aggregated
agg <- likdata$agg$agg
# point-referenced data
d12 <- ldgev12(likdata$y[!agg, 1], pmat[!agg, 1], pmat[!agg, 2], pmat[!agg, 3])
d1 <- do.call(cbind, lapply(1:4, function(i) likdata$X[[i]][!agg, , drop=FALSE] * d12[,i]))
d2 <- diag(0 * likdata$idpars, length(likdata$idpars))
id <- evgam:::.indices(4)[[2]]
for (i in 1:4) {
  for (j in 1:i) {
    temp <- crossprod(likdata$X[[i]][!agg, , drop=FALSE], likdata$X[[j]][!agg, , drop=FALSE] * d12[,id[i, j]])
    d2[likdata$idpars == i, likdata$idpars == j] <- temp
  }
}
# aggregated data
pdfr0 <- mapply("%*%", likdata$agg$X[1:2], pars[1:2])
pdfr0[,2] <- exp(pdfr0[,2])
pdfr <- ragged_mean_mat(pdfr0, likdata$agg$nxy)
d12 <- evgam:::ldgevagg12(likdata$y[agg, 1], pdfr[,1], log(pdfr[,2]), pmat[agg, 3], pmat[agg, 4])
dmb_dm <- ragged_mean_mat(likdata$agg$X[[1]], likdata$agg$nxy)
d1_1 <- d12[,1] * dmb_dm
dpb_dp <- ragged_mean_mat(pdfr0[,2] * likdata$agg$X[[2]], likdata$agg$nxy)
d1_2 <- (d12[,2] / pdfr[,2]) * dpb_dp
d1_34 <- do.call(cbind, lapply(3:4, function(i) likdata$X[[i]][agg, , drop=FALSE] * d12[,i]))
d1 <- rbind(d1, cbind(d1_1, d1_2, d1_34))
d21 <- evgam:::ldgevagg12(likdata$y[agg, 1], pdfr[,1], log(pdfr[,2]), pmat[agg, 3], pmat[agg, 4])
likdata$X[[1]] <- dmb_dm
likdata$X[[2]] <- dpb_dp / pdfr[,2]
for (i in 3:4)
  likdata$X[[i]] <- likdata$X[[i]][agg, , drop=FALSE]
for (i in 1:4) {
  idi <- likdata$idpars == i
  for (j in 1:i) {
    idj <- likdata$idpars == j
    d2[idi, idj] <- d2[idi, idj] + crossprod(likdata$X[[i]], d12[, id[i, j]] * likdata$X[[j]])
    if (j < i)
      d2[idj, idi] <- t(d2[idi, idj, drop=FALSE])
  }
}
t1 <- crossprod(likdata$agg$X[[2]], (pdfr0[,2] * rep(d12[,2] / pdfr[,2], likdata$agg$nxy)) * likdata$agg$X[[2]] * likdata$agg$weights)
t2 <- crossprod(likdata$X[[2]], d12[,2] * likdata$X[[2]])
d2[likdata$idpars == 2, likdata$idpars == 2] <- d2[likdata$idpars == 2, likdata$idpars == 2] + t1 - t2
if (!sandwich)
  d1 <- evgam:::.sum_col(d1)
list(d1, d2)
}

.gevaggfns <- list(d0=.gevagg.d0, d120=.gevagg.d12, d340=NULL)

.sum_col <- function(x) crossprod(rep(1, nrow(x)), x)[1,]
