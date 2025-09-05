## functions for calculating log(|H|)

# .indices <- function(i) {
# out <- list(list(i1 = structure(1, .Dim = 1L), i2 = structure(2, .Dim = c(1L, 
# 1L)), i3 = structure(1, .Dim = c(1L, 1L, 1L)), i4 = structure(2, .Dim = c(1L, 
# 1L, 1L, 1L))), list(i1 = structure(c(1, 2), .Dim = 2L), i2 = structure(c(3, 
# 4, 4, 5), .Dim = c(2L, 2L)), i3 = structure(c(1, 2, 2, 3, 2, 
# 3, 3, 4), .Dim = c(2L, 2L, 2L)), i4 = structure(c(5, 6, 6, 7, 
# 6, 7, 7, 8, 6, 7, 7, 8, 7, 8, 8, 9), .Dim = c(2L, 2L, 2L, 2L))), 
#     list(i1 = structure(c(1, 2, 3), .Dim = 3L), i2 = structure(c(4, 
#     5, 6, 5, 7, 8, 6, 8, 9), .Dim = c(3L, 3L)), i3 = structure(c(1, 
#     2, 3, 2, 4, 5, 3, 5, 6, 2, 4, 5, 4, 7, 8, 5, 8, 9, 3, 5, 
#     6, 5, 8, 9, 6, 9, 10), .Dim = c(3L, 3L, 3L)), i4 = structure(c(11, 
#     12, 13, 12, 14, 15, 13, 15, 16, 12, 14, 15, 14, 17, 18, 15, 
#     18, 19, 13, 15, 16, 15, 18, 19, 16, 19, 20, 12, 14, 15, 14, 
#     17, 18, 15, 18, 19, 14, 17, 18, 17, 21, 22, 18, 22, 23, 15, 
#     18, 19, 18, 22, 23, 19, 23, 24, 13, 15, 16, 15, 18, 19, 16, 
#     19, 20, 15, 18, 19, 18, 22, 23, 19, 23, 24, 16, 19, 20, 19, 
#     23, 24, 20, 24, 25), .Dim = c(3L, 3L, 3L, 3L))), list(i1 = structure(c(1, 
#     2, 3, 4), .Dim = 4L), i2 = structure(c(5, 6, 7, 8, 6, 9, 
#     10, 11, 7, 10, 12, 13, 8, 11, 13, 14), .Dim = c(4L, 4L)), 
#         i3 = structure(c(1, 2, 3, 4, 2, 5, 6, 7, 3, 6, 8, 9, 
#         4, 7, 9, 10, 2, 5, 6, 7, 5, 11, 12, 13, 6, 12, 14, 15, 
#         7, 13, 15, 16, 3, 6, 8, 9, 6, 12, 14, 15, 8, 14, 17, 
#         18, 9, 15, 18, 19, 4, 7, 9, 10, 7, 13, 15, 16, 9, 15, 
#         18, 19, 10, 16, 19, 20), .Dim = c(4L, 4L, 4L)), i4 = structure(c(21, 
#         22, 23, 24, 22, 25, 26, 27, 23, 26, 28, 29, 24, 27, 29, 
#         30, 22, 25, 26, 27, 25, 31, 32, 33, 26, 32, 34, 35, 27, 
#         33, 35, 36, 23, 26, 28, 29, 26, 32, 34, 35, 28, 34, 37, 
#         38, 29, 35, 38, 39, 24, 27, 29, 30, 27, 33, 35, 36, 29, 
#         35, 38, 39, 30, 36, 39, 40, 22, 25, 26, 27, 25, 31, 32, 
#         33, 26, 32, 34, 35, 27, 33, 35, 36, 25, 31, 32, 33, 31, 
#         41, 42, 43, 32, 42, 44, 45, 33, 43, 45, 46, 26, 32, 34, 
#         35, 32, 42, 44, 45, 34, 44, 47, 48, 35, 45, 48, 49, 27, 
#         33, 35, 36, 33, 43, 45, 46, 35, 45, 48, 49, 36, 46, 49, 
#         50, 23, 26, 28, 29, 26, 32, 34, 35, 28, 34, 37, 38, 29, 
#         35, 38, 39, 26, 32, 34, 35, 32, 42, 44, 45, 34, 44, 47, 
#         48, 35, 45, 48, 49, 28, 34, 37, 38, 34, 44, 47, 48, 37, 
#         47, 51, 52, 38, 48, 52, 53, 29, 35, 38, 39, 35, 45, 48, 
#         49, 38, 48, 52, 53, 39, 49, 53, 54, 24, 27, 29, 30, 27, 
#         33, 35, 36, 29, 35, 38, 39, 30, 36, 39, 40, 27, 33, 35, 
#         36, 33, 43, 45, 46, 35, 45, 48, 49, 36, 46, 49, 50, 29, 
#         35, 38, 39, 35, 45, 48, 49, 38, 48, 52, 53, 39, 49, 53, 
#         54, 30, 36, 39, 40, 36, 46, 49, 50, 39, 49, 53, 54, 40, 
#         50, 54, 55), .Dim = c(4L, 4L, 4L, 4L))))
# out[[i]]
# }

.index_arrays <- function(derivs, i) {
  
  if (i == 1) {
    
    ind <- lapply(seq_along(derivs), function(i) array(i, rep(1, derivs[i])))
  
  } else {
  
    out <- list()
    
    for (j in derivs) {
      t1 <- apply(expand.grid(lapply(1:j, function(x) 1:i)), 1, sort, simplify = FALSE)
      t2 <- do.call(rbind, lapply(t1, sort))
      t3 <- t2[!duplicated(t2), , drop = FALSE]
      out[[j]] <- list(t2, t3)
    }
    
    give_index <- function(x) apply(x, 1, paste0, collapse = '')
    i1 <- lapply(lapply(out[derivs], '[[', 1), give_index)
    i0 <- unique(unlist(i1))
    i2 <- lapply(i1, match, i0)
    # reps <- log(sapply(i2, length)) / log(i)
    reps <- logb(sapply(i2, length), base = i)
    ind <- mapply(array, i2, dim = lapply(reps, function(z) rep(i, z)))
  
  }
  
  return(ind)
    
}

.indices <- function(i) {
  d12 <- .index_arrays(1:2, i)
  d34 <- .index_arrays(3:4, i)
  out <- c(d12, d34)
  names(out) <- paste('i', seq_along(out), sep = '')
  out
}

.d1H0 <- function(dbeta, likdata, likfns) {

X <- likdata$X
idpars <- likdata$idpars
CH <- likdata$CH

nb <- nrow(dbeta$d1)
nsp <- ncol(dbeta$d1)
nX <- length(X)
n <- nrow(X[[1]])

ind <- .indices(nX)

beta <- likdata$compmode + likdata$CH %*% (dbeta$d0 - likdata$compmode)

GH <- likdata$k * likfns$d340(beta, likdata)

dbeta$d1 <- CH %*% dbeta$d1

d1H <- array(NA, c(nb, nb, nsp))

for (i in 1:nX) {
  for (j in 1:nX) {
    bX <- matrix(NA, n, nb)
    for (k in 1:nX)
      bX[,idpars == k] <- X[[k]] * GH[,ind$i3[i, j, k]]
    v <- bX %*% dbeta$d1
    for (l in 1:nsp)
      d1H[idpars == i, idpars == j, l] <- crossprod(X[[i]], v[,l] * X[[j]])
  }
}

d1H <- lapply(1:nsp, function(i) crossprod(CH, d1H[, , i] %*% CH))

list(d1=d1H, GH=GH)

}

.d2H0_diag <- function(dbeta, likdata, GH, H) {

X <- likdata$X
idpars <- likdata$idpars
CH <- likdata$CH

nb <- nrow(dbeta$d1)
nsp <- ncol(dbeta$d1)
nX <- length(X)
n <- nrow(X[[1]])

ind <- .indices(nX)

dbeta$d1 <- CH %*% dbeta$d1

d1eta <- lapply(seq_len(nX), function(i) X[[i]] %*% dbeta$d1[idpars == i, , drop=FALSE])
d2eta <- lapply(seq_len(nX), function(i) apply(dbeta$d2[idpars == i, , , drop=FALSE], 2:3, function(x) X[[i]] %*% x))

eH <- eigen(H$H)
V <- CH %*% eH$vectors
XV <- lapply(seq_len(nX), function(i) X[[i]] %*% (H$dH[idpars == i] * V[idpars == i, , drop=FALSE]))

d2H <- array(0, c(nb, nsp, nsp))

for (k in 1:nsp) {
  for (l in 1:k) {

for (i in 1:nX) {
  for (j in 1:nX) {
    v <- numeric(n)
    for (r in 1:nX) {
      v <- v + GH[,ind$i3[i, j, r]] * d2eta[[r]][, k, l]
      for (t in 1:nX)
        v <- v + d1eta[[r]][,k] * d1eta[[t]][,l] * GH[,ind$i4[i, j, r, t]]    
    }
    d2H[, k, l] <- d2H[, k, l] + colSums(XV[[i]] * XV[[j]] * v)
  }
}
if (l != k)
  d2H[, l, k] <- d2H[, k, l]
}
}

trd2H <- apply(d2H, 2:3, function(x) sum(x / eH$values))

list(d2=trd2H)

}

.d2beta <- function(dbeta, gradH, spSl, H) {
nsp <- ncol(dbeta$d1)
nb <- nrow(dbeta$d1)
d2beta <- array(NA, c(nb, nsp, nsp))
for (k in 1:nsp) {
  for (l in 1:k) {
    temp <- gradH[[k]] %*% dbeta$d1[,l]
    temp <- temp - spSl[[k]] %*% dbeta$d1[,l]
    temp <- temp - spSl[[l]] %*% dbeta$d1[,k]
    if (k == l) 
      temp <- temp - dbeta$spSlb[[k]]
    temp <- .precond_solve(H$cH, temp)
    d2beta[, k, l] <- temp
    if (l != k)
      d2beta[, l, k] <- d2beta[, k, l]
  }
}
dbeta$d2 <- d2beta
dbeta
}
