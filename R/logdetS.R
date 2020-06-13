## functions for calculating log(|S|)

.makeBlockMatrix <- function(lst) {
lst <- lapply(lst, as.matrix)
nrows <- sapply(lst, nrow)
starts <- cumsum(c(0, nrows))
n <- starts[length(starts)]
out <- matrix(0, n, n)
for (i in seq_along(lst)) out[starts[i] + seq_len(nrow(lst[[i]])), starts[i] + seq_len(nrow(lst[[i]]))] <- lst[[i]]
out
}

.SimilarityTransform <- function(S) {
eigS <- eigen(S, symmetric=TRUE)
U <- eigS$vectors
D <- diag(eigS$values)
k <- sum(D > .Machine$double.eps^.8 * max(D))
U[,1:k]
}

.ldetS <- function(lsp, Sl, useHessian, S=NULL) {
sp <- exp(lsp)
for (i in seq_along(Sl)) Sl[[i]] <- sp[i] * Sl[[i]]
if (is.null(S)) S <- Reduce("+", Sl)
out <- list()
U <- try(.SimilarityTransform(S))
if (inherits(U, "try-error")) {
out$ldetS <- 1e20
out$gradldetS <- rep(1e20, length(lsp))
out$hessldetS <- diag(out$gradldetS)
} else {
S1 <- crossprod(U, S %*% U)
out$ldetS <- determinant(S1)[[1]]
Sl <- lapply(Sl, function(x) crossprod(U, x %*% U))
Sl <- lapply(Sl, function(x) .5 * (t(x) + x))
S1 <- pinv(S1)
S1Sl <- lapply(Sl, function(x) S1 %*% x)
out$gradldetS <- sapply(S1Sl, function(x) sum(diag(x)))
if (useHessian) {
H <- diag(out$gradldetS)
for (j in seq_along(Sl)) for (k in 1:j) {
temp <- sp[j] * sp[k] * sum(t(S1Sl[[k]]) * S1Sl[[j]])
if (i == j) {
H[j, k] <- H[j, k] - temp
} else {
H[j, k] <- H[k, j] <- -temp
}
}
out$hessldetS <- .5 * (H + t(H))
}
}
return(out)
}

.logdetSmultiple <- function(smthlst, rho, join=FALSE, useHessian=FALSE) {
sp <- exp(rho)
if (length(smthlst$S) == 1) {
out <- list(ldetS=rho * smthlst$rank, gradldetS=smthlst$rank, hessldetS=0)
} else {
out <- .ldetS(rho, smthlst$S, useHessian=useHessian)
}
return(out)
}

.logdetS <- function(smth, rho, deriv=0) {
smth <- subset(smth, sapply(smth, length) == 5)
id <- rep(seq_along(smth), sapply(smth, function(x) length(x$S)))
rho <- split(rho, id)
temp <- lapply(seq_along(rho), function(i) .logdetSmultiple(smth[[i]], rho[[i]], useHessian=deriv > 1))
ldetS <- sum(unlist(sapply(temp, function(x) x$ldetS)))
gradldetS <- unlist(sapply(temp, function(x) x$gradldetS))
out <- list(d0=ldetS, d1=gradldetS)
if (deriv > 1) {
hessldetS <- lapply(temp, function(x) x$hessldetS)
out$d2 <- .makeBlockMatrix(hessldetS)
}
return(out)
}
