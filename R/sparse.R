.joinSmooth <- function(lst, sparse = FALSE) {
  nms <- c("S", "first.para", "last.para", "rank", "null.space.dim")
  nbi <- sapply(lst, function(x) x$nb)
  starts <- cumsum(c(0, nbi))
  if (length(lst) == 1) {
    lst <- lst[[1]]$smooth
  } else {
    lst <- lapply(lst, function(x) x$smooth)
    smooths <- sapply(lst, length) > 0
    for (i in which(smooths)) {
      for (j in seq_along(lst[[i]])) {
        lst[[i]][[j]]$first.para <- lst[[i]][[j]]$first.para + starts[i]
        lst[[i]][[j]]$last.para <- lst[[i]][[j]]$last.para + starts[i]
      }
    }
    lst <- unlist(lst, recursive=FALSE)
  }
  out <- lapply(lst, function(x) subset(x, names(x) %in% nms))
  lst2 <- list()
  nb <- tail(starts, 1)
  for (i in seq_along(lst)) {
    temp <- list()
    ind <- lst[[i]]$first.para:lst[[i]]$last.para
    for (j in seq_along(lst[[i]]$S)) {
      temp2 <- matrix(0, nb, nb)
      temp2[ind, ind] <- lst[[i]]$S[[j]]
      #     }
      temp[[j]] <- temp2
    }
    lst2[[i]] <- temp
  }
  lst2 <- unlist(lst2, recursive=FALSE)
  for (i in seq_along(out)) {
    if (length(out[[i]]$null.space.dim) == 0)
      out[[i]]$null.space.dim <- 0
  }
  if (sparse)
    lst2 <- lapply(lst2, as, 'dgCMatrix')
  attr(out, "Sl") <- lst2
  attr(out, "nb") <- nb
  out
}

.gH_dense <- function(x, likdata, sandwich=FALSE, deriv=2) {
  nX <- length(likdata$X)
  if (nX == 1) {
    out <- .gH1(x, likdata$X[[1]], likdata$dupid, likdata$duplicate, as.integer(sandwich), deriv)
  } else {
    if (nX == 2) {
      out <- .gH2(x, likdata$X[[1]], likdata$X[[2]], likdata$dupid, likdata$duplicate, as.integer(sandwich), deriv)
    } else {
      if (nX == 3) {
        out <- .gH3(x, likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$dupid, likdata$duplicate, as.integer(sandwich), deriv)
      } else {
        if (nX == 4) { # added with evgam_0.1.2 (05/04/2020)
          out <- .gH4(x, likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$dupid, likdata$duplicate, as.integer(sandwich), deriv)
        } else {
          if (nX == 5) { # added with evgam_0.1.5 (29/06/2021)
            out <- .gH5(x, likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$dupid, likdata$duplicate, as.integer(sandwich), deriv)
          } else {
            if (nX == 6) { # added with evgam_1.0.1 (31/08/2022)
              out <- .gH6(x, likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$X[[6]], likdata$dupid, likdata$duplicate, as.integer(sandwich), deriv)
            } else {
              stop("Number of model parameters not in {1, 2, 3, 4, 5, 6}")
            }
          }
        }
      }
    }
  }
  out
}

.gH_sparse <- function(x, likdata, sandwich=FALSE, deriv=2) {
  nX <- length(likdata$X)
  if (nX == 1) {
    out <- .gHsp1(x, likdata$X[[1]], likdata$dupid, likdata$duplicate, as.integer(sandwich), deriv)
  } else {
    if (nX == 2) {
      out <- .gHsp2(x, likdata$X[[1]], likdata$X[[2]], likdata$dupid, likdata$duplicate, as.integer(sandwich), deriv)
    } else {
      if (nX == 3) {
        out <- .gHsp3(x, likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$dupid, likdata$duplicate, as.integer(sandwich), deriv)
      } else {
        if (nX == 4) { # added with evgam_0.1.2 (05/04/2020)
          out <- .gHsp4(x, likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$dupid, likdata$duplicate, as.integer(sandwich), deriv)
        } else {
          if (nX == 5) { # added with evgam_0.1.5 (29/06/2021)
            out <- .gHsp5(x, likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$dupid, likdata$duplicate, as.integer(sandwich), deriv)
          } else {
            if (nX == 6) { # added with evgam_1.0.1 (31/08/2022)
              out <- .gHsp6(x, likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$X[[4]], likdata$X[[5]], likdata$X[[6]], likdata$dupid, likdata$duplicate, as.integer(sandwich), deriv)
            } else {
              stop("Number of model parameters not in {1, 2, 3, 4, 5, 6}")
            }
          }
        }
      }
    }
  }
  out
}

.gH <- function(x, likdata, sandwich=FALSE, deriv=2) {
  if (!likdata$sparse) {
    out <- .gH_dense(x, likdata, sandwich, deriv)
  } else {
    out <- .gH_sparse(x, likdata, sandwich, deriv)
  }
  out
}

.gH.pen.dense <- function(pars, likdata, likfns, deriv=2) {
  temp <- .gH.nopen(pars, likdata, likfns, deriv=deriv)
  temp[[1]] <- temp[[1]] + crossprod(pars, likdata$S)[1, ]
  if (deriv > 1) {
    attr(temp, "PP") <- attr(temp, "PP") + likdata$S / norm(likdata$S, "F")
    temp[[2]] <- temp[[2]] + likdata$S
  }
  temp
}

.gH.pen.sparse <- function(pars, likdata, likfns, deriv=2) {
  temp <- .gH.nopen(pars, likdata, likfns, deriv=deriv)
  if (!likdata$sparse) {
    temp[[1]] <- temp[[1]] + base::crossprod(pars, likdata$S)[1, ]
  } else {
    temp[[1]] <- temp[[1]] + Matrix::crossprod(pars, likdata$S)[1, ]
  }
  if (deriv > 1) {
    attr(temp, "PP") <- attr(temp, "PP") + likdata$S / Matrix::norm(likdata$S, "F")
    temp[[2]] <- temp[[2]] + likdata$S
  }
  temp
}

.gH.pen <- function(pars, likdata, likfns, deriv=2) {
  if (likdata$sparse) {
    out <- .gH.pen.sparse(pars, likdata, likfns, deriv)
  } else {
    out <- .gH.pen.dense(pars, likdata, likfns, deriv)
  }
  out
}

.search.dir.dense <- function(g, H, kept=NULL) {
  if (is.null(kept))
    kept <- !logical(length(g))
  H0 <- H
  g0 <- g
  g <- g[kept]
  H <- H[kept, kept, drop=FALSE]
  if (any(!is.finite(g)))
    stop("Some gradient non-finite")
  if (any(!is.finite(H)))
    stop("Some Hessian non-finite")
  H2 <- .precondition(H)
  H2 <- .perturb(H2)
  R <- attr(H2, "chol")
  d <- attr(H2, "d")
  out <- numeric(length(kept))
  piv <- ipiv <- attr(R, "pivot")
  ipiv[piv] <- seq_len(length(piv))
  out[kept] <- .precond_solve_dense(attr(H2, "chol"), g)
  g0[!kept] <- 0
  attr(out, "gradient") <- g0
  attr(out, "Hessian") <- H0
  attr(out, "cholH") <- R
  attr(out, "rank") <- attr(H2, 'rank0')
  out
}

.search.dir.sparse <- function(g, H, kept=NULL) {
  if (is.null(kept))
    kept <- !logical(length(g))
  H0 <- H
  g0 <- g
  g <- g[kept]
  H <- H[kept, kept, drop=FALSE]
  if (any(!is.finite(g)))
    stop("Some gradient non-finite")
  if (any(!is.finite(H)))
    stop("Some Hessian non-finite")
  out <- numeric(length(kept))
  H2 <- .precondition(H)
  H2 <- .perturb(H2)
  attr(attr(H2, 'chol'), 'H2') <- H2
  out[kept] <- .precond_solve_sparse(attr(H2, 'chol'), g)
  g0[!kept] <- 0
  attr(out, "gradient") <- g0
  attr(out, "Hessian") <- H2
  attr(out, "cholH") <- attr(H2, 'chol')
  attr(out, "rank") <- attr(H2, 'rank0')
  out
}

.search.dir <- function(g, H, kept=NULL) {
  if (inherits(H, 'matrix')) {
    out <- .search.dir.dense(g, H, kept)
  } else {
    out <- .search.dir.sparse(g, H, kept)
  }
  out
}

.Hdata_dense <- function(H) {
  out <- list(H0=H)
  H2 <- .precondition(H)
  H2 <- .perturb(H2)
  out$H <- H2
  out$dH <- attr(H2, "d")
  out$cH <- attr(H2, "chol")
  out$iH <- out$dH * t(chol2inv(out$cH) * out$dH)#crossprod(backsolve(out$cH, diag(out$dH), transpose=TRUE))
  out$kept <- !logical(nrow(H))
  out
}

.Hdata_sparse <- function(H) {
  out <- list(H0=H)
  H2 <- .precondition(H)
  H2 <- .perturb(H2)
  out$H <- H2
  out$dH <- attr(H2, "d")
  out$cH <- attr(H2, "chol")
  out$iH <- out$dH * Matrix::t(Matrix::chol2inv(out$cH) * out$dH)# out$iH <- Matrix::solve(out$cH, Matrix::Diagonal(out$dH), transpose=TRUE))
  out$kept <- !logical(nrow(H))
  out
}

.Hdata <- function(H) {
  if (inherits(H, 'matrix')) {
    out <- .Hdata_dense(H)
  } else {
    out <- .Hdata_sparse(H)
  }
  out
}

.perturb_dense <- function(A) {
  d <- attr(A, "d")
  eps <- 1e-16
  cholA <- suppressWarnings(chol(A, pivot=TRUE))
  rk <- rk0 <- attr(cholA, "rank")
  while(rk < nrow(A)) {
    diag(A) <- diag(A) + eps
    cholA <- suppressWarnings(chol(A, pivot=TRUE))
    rk <- attr(cholA, "rank")
    eps <- 1e2 * eps
  }
  attr(cholA, 'd') <- d
  attr(A, "chol") <- cholA
  attr(A, "rank") <- rk
  attr(A, "rank0") <- rk0
  return(A)
}

.perturb_sparse <- function(A) {
  d <- attr(A, 'd')
  eps <- 1e-12
  rk0 <- Matrix::rankMatrix(A, method = 'qr')
  # cholA <- Matrix::Cholesky(A)
  cholA <- suppressWarnings(try(Matrix::chol(A), silent = TRUE))
  # while(any(Matrix::diag(cholA) <= 0)) {
  while(inherits(cholA, 'try-error')) {
    Matrix::diag(A) <- Matrix::diag(A) + eps
    # cholA <- Matrix::Cholesky(A)
    cholA <- suppressWarnings(try(Matrix::chol(A), silent = TRUE))
    #rk <- Matrix::rankMatrix(A, method = 'qr')
    eps <- 1e2 * eps
  }
  # cholA <- Matrix::Cholesky(A)
  attr(cholA, 'd') <- d
  attr(A, "chol") <- cholA
  # attr(A, "chol2") <- Matrix::Cholesky(A)
  attr(A, "rank") <- Matrix::rankMatrix(A, method = 'qr')
  attr(A, "rank0") <- rk0
  return(A)
}

.perturb <- function(A) {
  if (inherits(A, 'matrix')) {
    out <- .perturb_dense(A)
  } else {
    out <- .perturb_sparse(A)
  }
  out
}

.precondition_dense <- function(H) {
  d <- 1/sqrt(abs(pmax(diag(H), 1e-8)))
  out <- d * t(t(H) * d)
  attr(out, "d") <- d
  out
}

.precondition_sparse <- function(H) {
  d <- 1/sqrt(abs(pmax(Matrix::diag(H), 1e-8)))
  out <- d * Matrix::t(Matrix::t(H) * d)
  attr(out, "d") <- d
  out
}

.precondition <- function(H) {
  if (inherits(H, 'matrix')) {
    out <- .precondition_dense(H)
  } else {
    out <- .precondition_sparse(H)
  }
  out
}

.precond_solve_dense <- function(L, x) {
  d <- attr(L, "d")
  x <- d * as.matrix(x)
  piv <- ipiv <- attr(L, "pivot")
  ipiv[piv] <- seq_len(nrow(L))
  d * backsolve(L, backsolve(L, x[piv, , drop=FALSE], upper.tri=TRUE, transpose=TRUE))[ipiv, , drop=FALSE]
}

.precond_solve_sparse <- function(L, x) {
  d <- attr(L, 'd')
  # d * Matrix::solve(L, x * d)
  d * Matrix::solve(L, Matrix::solve(Matrix::t(L), x * d))
}

.precond_solve <- function(L, x) {
  if (is.null(attr(L, 'd')))
    attr(L, 'd') <- rep(1, nrow(L))
  if (inherits(L, 'matrix')) {
    out <- .precond_solve_dense(L, x)
  } else {
    out <- .precond_solve_sparse(L, x)
  }
  out
}

.d0logdetH <- function(x) {
  if (is.null(x$cholHessian)) {
    out <- as.vector(determinant(x$Hessian)[[1]])
  } else {
    cH <- x$cholHessian
    if (inherits(cH, 'Cholesky')) {
      out <- 2 * as.vector(Matrix::determinant(cH)$modulus)
    } else {
      if (inherits(cH, c('dtCMatrix', 'dCHMsimpl', 'dtrMatrix'))) {
        out <- 2 * sum(log(Matrix::diag(cH)))
        if (inherits(cH, c('dtCMatrix', 'dCHMsimpl', 'dtrMatrix'))) {
          out <- out - 2 * sum(log(attr(cH, "d")))
        }
      } else {
        out <- 2 * sum(log(diag(cH)))
        out <- out - 2 * sum(log(attr(cH, "d")))
      }
    }
  }
  list(d0 = out)
}

.d1beta <- function(lsp, beta, spSl, H) {
  out <- list(d0 = beta)
  spSlb <- sapply(spSl, function(x) x %*% beta)
  if (inherits(spSlb, 'list'))
    spSlb <- sapply(spSlb, as.vector)
  out$d1 <- -.precond_solve(H$cH, spSlb)
  out$spSlb <- spSlb
  out
}

.VpVc <- function(fitreml, likfns, likdata, Sdata, correctV, sandwich, smooths, trace) {
  lsp <- fitreml$par
  H0 <- .gH.nopen(fitreml$beta, likdata, likfns)[[2]]
  if (smooths) {
    sp <- exp(lsp)
    H <- H0 + likdata$S
  } else {
    H <- H0
  }
  Hd <- .Hdata(H)
  if (!likdata$sparse) {
    cholH <- try(chol(H), silent = TRUE)
  } else {
    cholH <- suppressWarnings(try(Matrix::chol(H), silent = TRUE))
  }
  if (inherits(cholH, "try-error") & trace >= 0)
    message("Final Hessian of negative penalized log-likelihood not numerically positive definite.")
  Vc <- Vp <- Hd$iH
  if (smooths) {
    if (correctV) {
      cholVp <- try(chol(Vp), silent=TRUE)
      if (inherits(cholVp, "try-error")) {
        cholVp <- attr(.perturb(Vp), "chol")
      }
      attr(lsp, "beta") <- fitreml$beta
      spSl <- Map("*", attr(Sdata, "Sl"), exp(lsp))
      dbeta <- .d1beta(lsp, fitreml$beta, spSl, Hd)$d1
      Vrho <- fitreml$invHessian
      if (!likdata$sparse) {
        Vbetarho <- base::tcrossprod(dbeta %*% Vrho, dbeta)
      } else {
        Vbetarho <- Matrix::tcrossprod(dbeta %*% Vrho, dbeta)
      }
      VR <- matrix(0, nrow=likdata$nb, ncol=likdata$nb)
      Vc <- .perturb(Vp + Vbetarho + VR)
    } else {
      Vrho <- 0
    }
  } else {
    Vrho <- 0
  }
  list(Vp=Vp, Vc=Vc, Vlsp=Vrho, H0=H0, H=H)
}

.choltr <- function(x, b) {
  if (inherits(x, 'matrix')) {
    sum(diag(.precond_solve(x, b)))
  } else {
    sum(Matrix::diag(.precond_solve(x, b)))
  }
}

.d1H0_diag <- function(dbeta, likdata, likfns, H) {

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

  d1eta <- lapply(seq_len(nX), function(i) X[[i]] %*% dbeta$d1[idpars == i, , drop=FALSE])

  trd1H <- numeric(nsp)
  
  for (l in 1:nsp) {

    if (!likdata$sparse) {
      d1Hl <- matrix(0, nb, nb)
    } else {
      d1Hl <- Matrix::Matrix(0, nrow = nb, ncol = nb, sparse = TRUE)
    }

    for (i in 1:nX) {
      for (j in 1:nX) {
        v <- numeric(n)
        for (k in 1:nX)
          v <- v + GH[,ind$i3[i, j, k]] * d1eta[[k]][, l]
        if (!likdata$sparse) {
          d1Hl[idpars == i, idpars == j] <- crossprod(X[[i]], X[[j]] * v)
        } else {
          d1Hl[idpars == i, idpars == j] <- Matrix::crossprod(X[[i]], X[[j]] * v)
        }
      }  
    }

    trd1H[l] <- .choltr(H$cH, d1Hl)
  
  }

  list(d1 = trd1H)

}

.d1logdetH <- function(dbeta, likdata, likfns, spSl, H) {
  d1 <- .d1H0_diag(dbeta, likdata, likfns, H)$d1
  d1 <- d1 + sapply(spSl, function(x) .choltr(H$cH, x))
  list(d1 = d1, dbeta = dbeta)
}
