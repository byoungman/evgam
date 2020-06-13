## functions for calculating log(|H|)

.makeIndices <- function(K) {

i1 <- array(NA, rep(K, 1))
i2 <- array(NA, rep(K, 2))
i3 <- array(NA, rep(K, 3))
i4 <- array(NA, rep(K, 4))

ii <- ij <- ijk <- ijkl <- 1

for (i in 1:K) {
i1[i] <- ii
ii <- ii + 1
for (j in i:K) {
i2[i, j] <- i2[j, i] <- ij
ij <- ij + 1
for (k in j:K) {
i3[i, j, k] <- i3[i, k, j] <- i3[j, i, k] <- i3[j, k, i] <- i3[k, i, j] <- i3[k, j, i] <- ijk
ijk <- ijk + 1
for (l in k:K) {
i4[i, j, k, l] <- i4[i, j, l, k] <- i4[i, k, j, l] <- i4[i, k, l, j] <- ijkl
i4[i, l, j, k] <- i4[i, l, k, j] <- ijkl
i4[j, i, k, l] <- i4[j, i, l, k] <- i4[j, k, i, l] <- i4[j, k, l, i] <- ijkl
i4[j, l, i, k] <- i4[j, l, k, i] <- ijkl
i4[k, i, j, l] <- i4[k, i, l, j] <- i4[k, j, i, l] <- i4[k, j, l, i] <- ijkl
i4[k, l, i, j] <- i4[k, l, j, i] <- ijkl
i4[l, i, j, k] <- i4[l, i, k, j] <- i4[l, j, i, k] <- i4[l, j, k, i] <- ijkl
i4[l, k, i, j] <- i4[l, k, j, i] <- ijkl
ijkl <- ijkl + 1
}
}
}
}

i2 <- i2 + max(i1)
i4 <- i4 + max(i3, na.rm=TRUE)

list(i1=i1, i2=i2, i3=i3, i4=i4)

}

.d1H0_diag <- function(GH, CH, X, idpars, dbeta, spSl, H) {

nb <- nrow(dbeta$d1)
nsp <- ncol(dbeta$d1)
nX <- length(X)
n <- nrow(X[[1]])

iHtX <- t(.dbind(X))
iHtX <- t(.precond_solve(H$cH, iHtX))
VX <- matrix(0, nrow(iHtX), ncol(iHtX))

X <- lapply(seq_len(nX), function(i) X[[i]] %*% attr(CH, "list")[[i]])
d1eta <- lapply(seq_len(nX), function(i) X[[i]] %*% dbeta$d1[idpars == i, , drop=FALSE])

ind <- .makeIndices(nX)

trd1H <- rep(NA, nsp)

for (k in 1:nsp) {
for (i in 1:nX) for (j in 1:nX) {
v <- numeric(n)
for (r in 1:nX) {
v <- v + d1eta[[r]][,k] * GH[,ind$i3[i, j, r]]
}
if (i == j) {
    rind <- 1:n + (i-1)*n
    VX[rind, idpars == i] <- v * X[[i]]
} else {
    rind1 <- 1:n + (i-1)*n
    rind2 <- 1:n + (j-1)*n
    VX[rind2, idpars == i] <- v * X[[i]]
    VX[rind1, idpars == j] <- v * X[[j]]
}
}
trd1H[k] <- sum(iHtX * VX)
}

list(d1=trd1H)

}


.d1H0 <- function(GH, CH, X, idpars, dbeta, spSl) {

nb <- nrow(dbeta$d1)
nsp <- ncol(dbeta$d1)
nX <- length(X)
n <- nrow(X[[1]])

dbeta$d1 <- CH %*% dbeta$d1
d1H <- array(NA, c(nb, nb, nsp))

ind <- .makeIndices(nX)

for (i in 1:nX) for (j in 1:nX) {
bX <- matrix(NA, n, nb)
for (k in 1:nX) {
bX[,idpars == k] <- X[[k]] * GH[,ind$i3[i, j, k]]
}
v <- bX %*% dbeta$d1
for (l in 1:nsp) {
d1H[idpars == i, idpars == j, l] <- crossprod(X[[i]], v[,l] * X[[j]])
}
}
for (l in 1:nsp) {
    d1H[,,l] <- crossprod(CH, as.matrix(d1H[,,l]) %*% CH)
}

list(d1=d1H)

}

.dbind <- function(x) {
nx <- length(x)
nc <- sapply(x, ncol)
nr <- sapply(x, nrow)
cends <- cumsum(nc)
cstarts <- c(1, cends[seq_len(nx - 1)] + 1)
rends <- cumsum(nr)
rstarts <- c(1, rends[seq_len(nx - 1)] + 1)
out <- matrix(0, rends[nx], cends[nx])
for (i in seq_len(nx)) out[rstarts[i]:rends[i], cstarts[i]:cends[i]] <- x[[i]]
out
}

.d2H0_diag <- function(GH, CH, X, idpars, dbeta, spSl, H) {

nb <- nrow(dbeta$d1)
nsp <- ncol(dbeta$d1)
nX <- length(X)
n <- nrow(X[[1]])

iHtX <- tcrossprod(t(CH), .dbind(X))
iHtX <- t(.precond_solve(H$cH, iHtX))
VX <- matrix(0, nrow(iHtX), ncol(iHtX))

dbeta$d1 <- CH %*% dbeta$d1
dbeta$d2 <- apply(dbeta$d2, 2:3, function(x) CH %*% x)
d1eta <- lapply(seq_len(nX), function(i) X[[i]] %*% dbeta$d1[idpars == i, , drop=FALSE])
d2eta <- lapply(seq_len(nX), function(i) apply(dbeta$d2[idpars == i, , , drop=FALSE], 2:3, function(x) X[[i]] %*% x))

ind <- .makeIndices(nX)

trd2H <- matrix(NA, nsp, nsp)

for (k in 1:nsp) for (l in 1:nsp) {
for (i in 1:nX) for (j in 1:nX) {
v <- numeric(nrow(GH))
temp <- matrix(0, nX, nX)
for (r in 1:nX) {
v <- v + d2eta[[r]][,k,l] * GH[,ind$i3[i, j, r]]
for (t in 1:nX) {
v <- v + d1eta[[r]][,k] * d1eta[[t]][,l] * GH[,ind$i4[i, j, r, t]]
}
}
if (i == j) {
    rind <- 1:n + (i-1)*n
    VX[rind, idpars == i] <- v * X[[i]]
} else {
    rind1 <- 1:n + (i-1)*n
    rind2 <- 1:n + (j-1)*n
    VX[rind2, idpars == i] <- v * X[[i]]
    VX[rind1, idpars == j] <- v * X[[j]]
}
VX <- VX %*% CH
}
trd2H[k,l] <- sum(iHtX * VX)
}

list(d2=trd2H)

}

.d2H0 <- function(GH, CH, X, idpars, dbeta, spSl) {

nb <- nrow(dbeta$d1)
nsp <- ncol(dbeta$d1)
nX <- length(X)

d2H <- array(NA, c(nb, nb, nsp, nsp))
X <- lapply(seq_len(nX), function(i) X[[i]] %*% attr(CH, "list")[[i]])
d1eta <- lapply(seq_len(nX), function(i) X[[i]] %*% dbeta$d1[idpars == i, , drop=FALSE])
d2eta <- lapply(seq_len(nX), function(i) apply(dbeta$d2[idpars == i, , , drop=FALSE], 2:3, function(x) X[[i]] %*% x))

ind <- .makeIndices(nX)

for (k in 1:nsp) for (l in 1:nsp) {
for (i in 1:nX) for (j in 1:nX) {
v <- numeric(nrow(GH))
for (r in 1:nX) {
v <- v + d2eta[[r]][,k,l] * GH[,ind$i3[i, j, r]]
for (t in 1:nX) {
v <- v + d1eta[[r]][,k] * d1eta[[t]][,l] * GH[,ind$i4[i, j, r, t]]
}
}
d2H[idpars == i, idpars == j, k, l] <- crossprod(X[[i]], v * X[[j]])
}
}

list(d2=d2H)

}

.dbeta <- function(lsp, spSl, H, beta, likdata, likfns, deriv=1) {
sp <- exp(lsp)
nsp <- length(sp)
beta <- likdata$compmode + likdata$CH %*% (beta - likdata$compmode)
dbeta <- list()
dbeta$d1 <- -.precond_solve(H$cH, sapply(spSl, function(x) x %*% beta))
dbeta$GH <- likdata$k * likfns$d340(beta, likdata)
if (deriv > 1) {
  gradH <- .d1H0(dbeta$GH, likdata$CH, likdata$X, likdata$idpars, dbeta, spSl)$d1
  d2beta <- array(NA, c(length(beta), nsp, nsp))
  for (k in seq_len(nsp)) for (l in seq_len(k)) {
    temp0 <- gradH[,,l] %*% dbeta$d1[,k] + spSl[[l]] %*% dbeta$d1[,k] + spSl[[k]] %*% dbeta$d1[,l]
    temp1 <- -.precond_solve(H$cH, temp0)
    if (k == l) temp1 <- temp1 + dbeta$d1[,k]
    d2beta[, k, l] <- temp1
    if (k != l) d2beta[, l, k] <- d2beta[, k, l]
  }
  dbeta$d2 <- d2beta
  dbeta$gradH <- gradH
}
dbeta
}

.d0logdetH <- function(x) {
if (is.null(x$cholHessian)) {
  out <- as.vector(determinant(x$Hessian)[[1]])
} else {
  out <- 2 * sum(log(diag(x$cholHessian)))
  out <- out - 2 * sum(log(attr(x$cholHessian, "d")))
}
list(d0 = out) 
}

.d1logdetH <- function(lsp, likdata, likfns, Sdata, H) {
beta <- attr(lsp, "beta")
sp <- exp(lsp)
spSl <- lapply(seq_along(sp), function(i) sp[i] * attr(Sdata, "Sl")[[i]])
if (!all(H$kept)) browser()
dbeta <- .dbeta(lsp, spSl, H, beta, likdata, likfns, deriv=1)
d1 <- .d1H0_diag(dbeta$GH, likdata$CH, likdata$X, likdata$idpars, dbeta, spSl, H)$d1
d1 <- d1 + sapply(spSl, function(x) sum(diag(.precond_solve(H$cH, x))))
list(d1=d1, dbeta=dbeta)
}

.d1logdetH <- function(lsp, likdata, likfns, Sdata, H) {
beta <- attr(lsp, "beta")
sp <- exp(lsp)
spSl <- lapply(seq_along(sp), function(i) sp[i] * attr(Sdata, "Sl")[[i]])
if (!all(H$kept)) browser()
dbeta <- .dbeta(lsp, spSl, H, beta, likdata, likfns, deriv=2)
dbeta$gradH <- dbeta$gradH + unlist(spSl)
dbeta$gradH <- array(apply(dbeta$gradH, 3, function(x) .precond_solve(H$cH, x)), dim(dbeta$gradH))
list(d1=apply(dbeta$gradH, 3, function(x) sum(diag(x))))
}

.d12logdetH <- function(lsp, likdata, likfns, Sdata, H) {
beta <- attr(lsp, "beta")
sp <- exp(lsp)
spSl <- lapply(seq_along(sp), function(i) sp[i] * attr(Sdata, "Sl")[[i]])
dbeta <- .dbeta(lsp, spSl, H, beta, likdata, likfns, deriv=2)
d12 <- .d2H0_diag(dbeta$GH, likdata$CH, likdata$X, likdata$idpars, dbeta, spSl, H)$d2
diag(d12) <- diag(d12) + sapply(spSl, function(x) sum(diag(.precond_solve(H$cH, x))))
dbeta$gradH <- dbeta$gradH + unlist(spSl)
dbeta$gradH <- array(apply(dbeta$gradH, 3, function(x) .precond_solve(H$cH, x)), dim(dbeta$gradH))
d22 <- matrix(0, length(lsp), length(lsp))
for (k in seq_along(lsp)) for (l in seq_len(k)) {
  d22[k, l] <- sum(t(dbeta$gradH[,,l]) * dbeta$gradH[,,k])
  if (l != k) d22[l, k] <- d22[k, l]
}
out <- list(d1=apply(dbeta$gradH, 3, function(x) sum(diag(x))))
out$d2 <- d12 - d22
out$dbeta <- dbeta
out
}
