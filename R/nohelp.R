# ## generic functions ##
# 
.pivchol_rmvn <- function(n, mu, Sig) {
R <- suppressWarnings(chol(Sig, pivot = TRUE))
piv <- order(attr(R, "pivot"))  ## reverse pivoting index
r <- attr(R, "rank")  ## numerical rank
V <- R[1:r, piv]
Y <- crossprod(V, matrix(rnorm(n * r), r))
Y + as.vector(mu)
}

.precondition <- function(H) {
d <- 1/sqrt(abs(pmax(diag(H), 1e-8)))
out <- d * t(t(H) * d)
attr(out, "d") <- d
out
}

.perturb <- function(A) {
d <- attr(A, "d")
eps <- 1e-16
cholA <- suppressWarnings(chol(A, pivot=TRUE))
while(attr(cholA, "rank") < nrow(A)) {
diag(A) <- diag(A) + eps
cholA <- suppressWarnings(chol(A, pivot=TRUE))
eps <- 1e2 * eps
}
attr(A, "chol") <- cholA
return(A)
}

.solve_evgam <- function(H) {
H2 <- .precondition(H)
H2 <- .perturb(H2)
D <- attr(H2, "d")
D <- diag(D, length(D))
R <- attr(H2, "chol")
crossprod(backsolve(R, D, transpose=TRUE))
}

.precond_solve <- function(L, x) {
d <- attr(L, "d")
x <- d * as.matrix(x)
piv <- ipiv <- attr(L, "pivot")
ipiv[piv] <- seq_len(nrow(L))
d * backsolve(L, backsolve(L, x[piv, , drop=FALSE], upper.tri=TRUE, transpose=TRUE))[ipiv, , drop=FALSE]
}

.fdHess <- function(pars, gr, ..., eps=1e-4) {
np <- length(pars)
g0 <- gr(pars, ...)
out <- vapply(seq_len(np), function(i) gr(replace(pars, i, pars[i] + eps), ...) - g0, double(np))
out <- out + t(out)
out <- .5 * out / eps
out
}

.bdiag <- function(x) {
dims <- sapply(x, dim)
x <- x[dims[1,] > 0 | dims[2,] > 0]
dims <- sapply(x, dim)
row_ends <- cumsum(dims[1,])
row_starts <- c(1, row_ends + 1)
col_ends <- cumsum(dims[2,])
col_starts <- c(1, col_ends + 1)
out <- matrix(0, tail(row_ends, 1), tail(col_ends, 1))
for (i in seq_along(x))
  out[row_starts[i]:row_ends[i], col_starts[i]:col_ends[i]] <- x[[i]]
out
}

## smoothing matrix manipulation functions

.makeS <- function(lst, sp) {
Reduce("+", Map("*", attr(lst, "Sl"), sp))
}

.joinSmooth <- function(lst) {
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
attr(out, "Sl") <- lst2
attr(out, "nb") <- nb
out
}

.grad.R <- function(pars, Sdata, R, eps, likfns, likdata, H0=NULL) {
S <- .makeS(Sdata, exp(pars))
if (is.null(H0)) H0 <- .hess.nopen(likdata$beta, likdata, likfns)
H <- H0 + S
H <- H[likdata$LAid, likdata$LAid]
cholH <- try(chol(H), silent=TRUE)
if (inherits(cholH, "try-error")) {
iH <- pinv(H)
} else {
iH <- chol2inv(cholH)
}
choliH <- try(chol(iH), silent=TRUE)
if (inherits(choliH, "try-error")) choliH <- attr(.perturb(iH), "chol")
(choliH - R) / eps
}

## other stuff

.padWithVec <- function(df, vec, nm) {
df <- df[1,]
id <- cbind(1, seq_along(vec))
df <- df[id[,1],]
df[,nm] <- vec
df
}

.addNodes <- function(df1, df2, id) {
df1 <- lapply(df1, rep, nrow(df2))
for (i in seq_along(id)) df1[[id[i]]] <- df2[,i]
as.data.frame(df1)
}

.outer.list <- function(lst) {
out <- lst[[1]]
if (length(lst) > 1) {
for (i in 2:length(lst)) out <- out %o% lst[[i]]
}
c(out)
}

.updateControl <- function(lst0, lst) {
for (i in seq_along(lst0)) {
nmi <- names(lst0)[i]
i2 <- which(names(lst) == nmi)
lsti <- lst[[i2]]
for (j in seq_along(lst0[[i]])) {
nmj <- names(lst0[[i]])[j]
j2 <- which(names(lsti) == nmj)
lsti[[j2]] <- lst0[[i]][[j]]
}
lst[[i2]] <- lsti
}
lst
}

## p-values for smooths
##
## To avoid use of mgcv::: this code is directly taken from mgcv 1.8-14 
## by Professor Simon Wood
##
## Thank you, Simon.
##
.smoothTest <- function(b,X,V,eps=.Machine$double.eps^.5) {
## Forms Cox, Koh, etc type test statistic, and
## obtains null distribution by simulation...
## if b are coefs f=Xb, cov(b) = V. z is a vector of 
## i.i.d. N(0,1) deviates

  qrx <- qr(X)
  R <- qr.R(qrx)
  V <- R%*%V[qrx$pivot,qrx$pivot]%*%t(R)
  V <- (V + t(V))/2
  ed <- eigen(V,symmetric=TRUE)
  k <- length(ed$values)
  ## could truncate, but it doesn't improve power in correlated case!
  f <- t(ed$vectors[,1:k])%*%R%*%b
  t <- sum(f^2)
  k <- ncol(X)
  lambda <- as.numeric(ed$values[1:k])
  pval <- .liu2(t,lambda) ## should really use Davies
  list(stat=t,pval=pval)  
} 

.liu2 <- function(x, lambda, h = rep(1,length(lambda)),lower.tail=FALSE) {
# Evaluate Pr[sum_i \lambda_i \chi^2_h_i < x] approximately.
# Code adapted from CompQuadForm package of Pierre Lafaye de Micheaux 
# and directly from....
# H. Liu, Y. Tang, H.H. Zhang, A new chi-square approximation to the 
# distribution of non-negative definite quadratic forms in non-central 
# normal variables, Computational Statistics and Data Analysis, Volume 53, 
# (2009), 853-856. Actually, this is just Pearson (1959) given that
# the chi^2 variables are central. 
# Note that this can be rubbish in lower tail (e.g. lambda=c(1.2,.3), x = .15)
  
#  if (FALSE) { ## use Davies exact method in place of Liu et al/ Pearson approx.
#    require(CompQuadForm)
#    r <- x
#    for (i in 1:length(x)) r[i] <- davies(x[i],lambda,h)$Qq
#    return(pmin(r,1))
#  }

  if (length(h) != length(lambda)) stop("lambda and h should have the same length!")
 
  lh <- lambda*h
  muQ <- sum(lh)
  
  lh <- lh*lambda
  c2 <- sum(lh)
  
  lh <- lh*lambda
  c3 <- sum(lh)
  
  s1 <- c3/c2^1.5
  s2 <- sum(lh*lambda)/c2^2

  sigQ <- sqrt(2*c2)

  t <- (x-muQ)/sigQ

  if (s1^2>s2) {
    a <- 1/(s1-sqrt(s1^2-s2))
    delta <- s1*a^3-a^2
    l <- a^2-2*delta
  } else {
    a <- 1/s1
    delta <- 0
    l <- c2^3/c3^2
  }

  muX <- l+delta
  sigX <- sqrt(2)*a
  
  return(pchisq(t*sigX+muX,df=l,ncp=delta,lower.tail=lower.tail))

}

.simf <- function(x,a,df,nq=50) {
## suppose T = sum(a_i \chi^2_1)/(chi^2_df/df). We need
## Pr[T>x] = Pr(sum(a_i \chi^2_1) > x *chi^2_df/df). Quadrature 
## used here. So, e.g.
## 1-pf(4/3,3,40);simf(4,rep(1,3),40);1-pchisq(4,3)
  p <- (1:nq-.5)/nq
  q <- qchisq(p,df)
  x <- x*q/df
  pr <- sum(.liu2(x,a)) ## Pearson/Liu approx to chi^2 mixture
  pr/nq 
}

.testStat <- function(p,X,V,rank=NULL,type=0,res.df= -1) {
## Implements Wood (2013) Biometrika 100(1), 221-228
## The type argument specifies the type of truncation to use.
## on entry `rank' should be an edf estimate
## 0. Default using the fractionally truncated pinv.
## 1. Round down to k if k<= rank < k+0.05, otherwise up.
## 2. Naive rounding.
## 3. Round up.
## 4. Numerical rank estimation, tol=1e-3
## res.df is residual dof used to estimate scale. <=0 implies
## fixed scale.

  qrx <- qr(X,tol=0)
  R <- qr.R(qrx)
  V <- R%*%V[qrx$pivot,qrx$pivot,drop=FALSE]%*%t(R)
  V <- (V + t(V))/2
  ed <- eigen(V,symmetric=TRUE)
  ## remove possible ambiguity from statistic...
  siv <- sign(ed$vectors[1,]);siv[siv==0] <- 1
  ed$vectors <- sweep(ed$vectors,2,siv,"*")

  k <- max(0,floor(rank)) 
  nu <- abs(rank - k)     ## fractional part of supplied edf
  if (type < -.5) { ## Crude modification of Cox and Koh
    res <- .smoothTest(p,X,V)
    res$rank <- rank
    return(res)
  } else  if (type==1) { ## round up is more than .05 above lower
    if (rank > k + .05||k==0) k <- k + 1
    nu <- 0;rank <- k
  } else if (type==2) { ## naive round
    nu <- 0;rank <- k <- max(1,round(rank))
    warning("p-values may give low power in some circumstances")
  } else if (type==3) { ## round up
    nu <- 0; rank <- k <- max(1,ceiling(rank))
    warning("p-values un-reliable")
  } else if (type==4) { ## rank estimation
    rank <- k <- max(sum(ed$values>1e-3*max(ed$values)),1) 
    nu <- 0
    warning("p-values may give very low power")
  }

  if (nu>0) k1 <- k+1 else k1 <- k

  ## check that actual rank is not below supplied rank+1
  r.est <- sum(ed$values > max(ed$values)*.Machine$double.eps^.9)
  if (r.est<k1) {k1 <- k <- r.est;nu <- 0;rank <- r.est}

  ## Get the eigenvectors...
  # vec <- qr.qy(qrx,rbind(ed$vectors,matrix(0,nrow(X)-ncol(X),ncol(X))))
  vec <- ed$vectors
  if (k1<ncol(vec)) vec <- vec[,1:k1,drop=FALSE]

  ## deal with the fractional part of the pinv...
  if (nu>0&&k>0) {
     if (k>1) vec[,1:(k-1)] <- t(t(vec[,1:(k-1)])/sqrt(ed$val[1:(k-1)]))
     b12 <- .5*nu*(1-nu)
     if (b12<0) b12 <- 0
     b12 <- sqrt(b12)
     B <- matrix(c(1,b12,b12,nu),2,2)
     ev <- diag(ed$values[k:k1]^-.5,nrow=k1-k+1)
     B <- ev%*%B%*%ev
     eb <- eigen(B,symmetric=TRUE)
     rB <- eb$vectors%*%diag(sqrt(eb$values))%*%t(eb$vectors)
     vec1 <- vec
     vec1[,k:k1] <- t(rB%*%diag(c(-1,1))%*%t(vec[,k:k1]))
     vec[,k:k1] <- t(rB%*%t(vec[,k:k1]))
  } else {
    vec1 <- vec <- if (k==0) t(t(vec)*sqrt(1/ed$val[1])) else
            t(t(vec)/sqrt(ed$val[1:k]))
    if (k==1) rank <- 1
  }
  ## there is an ambiguity in the choise of test statistic, leading to slight
  ## differences in the p-value computation depending on which of 2 alternatives 
  ## is arbitrarily selected. Following allows both to be computed and p-values
  ## averaged (can't average test stat as dist then unknown) 
  d <- t(vec)%*%(R%*%p)
  d <- sum(d^2) 
  d1 <- t(vec1)%*%(R%*%p)
  d1 <- sum(d1^2)
  ##d <- d1 ## uncomment to avoid averaging

  rank1 <- rank ## rank for lower tail pval computation below

  ## note that for <1 edf then d is not weighted by EDF, and instead is 
  ## simply refered to a chi-squared 1

  if (nu>0) { ## mixture of chi^2 ref dist
     if (k1==1) rank1 <- val <- 1 else { 
       val <- rep(1,k1) ##ed$val[1:k1]
       rp <- nu+1
       val[k] <- (rp + sqrt(rp*(2-rp)))/2
       val[k1] <- (rp - val[k])
     }
   
     if (res.df <= 0) pval <- (.liu2(d,val) + .liu2(d1,val))/2 else ##  pval <- davies(d,val)$Qq else
     pval <- (.simf(d,val,res.df) + .simf(d1,val,res.df))/2
  } else { pval <- 2 }
  ## integer case still needs computing, also liu/pearson approx only good in 
  ## upper tail. In lower tail, 2 moment approximation is better (Can check this 
  ## by simply plotting the whole interesting range as a contour plot!)
  if (pval > .5) { 
    if (res.df <= 0) pval <- (pchisq(d,df=rank1,lower.tail=FALSE)+pchisq(d1,df=rank1,lower.tail=FALSE))/2 else
    pval <- (pf(d/rank1,rank1,res.df,lower.tail=FALSE)+pf(d1/rank1,rank1,res.df,lower.tail=FALSE))/2
  }
  list(stat=d,pval=min(1,pval),rank=rank)
} ## end of testStat


.giveSmoothTest <- function(obj, j) {
id <- obj$smooth[[j]]$first.para:obj$smooth[[j]]$last.para
p <- coef(obj)[id]
X <- obj$X[,id]
V <- obj$Vp[id, id]
edf1 <- sum(obj$edf[id])
edf1 <- min(length(id), edf1)
stat <- .testStat(p, X, V, edf1)
out <- data.frame(stat[[3]], length(id), stat[[1]], stat[[2]])
rownames(out) <- obj$smooth[[j]]$label
out
}

.giveParametricTest <- function(obj) {
np <- length(coef(obj))
id <- !logical(np)
for (j in seq_along(obj$smooth)) {
if (!is.null(obj$smooth[[j]]$first.para)) {
idsmooth <- obj$smooth[[j]]$first.para:obj$smooth[[j]]$last.para
id[idsmooth] <- FALSE
}
}
id2 <- which(id)
beta <- coef(obj)[id]
ese <- sqrt(obj$Vp[cbind(id2, id2)])
tval <- beta / ese
pval <- pnorm(-abs(tval))
out <- try(data.frame(beta, ese, tval, pval))
rownames(out) <- colnames(obj$X[,id, drop=FALSE])
out
}

.tidySmoothTable <- function(tab) {
for (i in 1:3) tab[,i] <- round(tab[,i], 2)
tab[,4] <- signif(tab[,4], 3)
tab[,4] <- replace(unlist(tab[,4]), unlist(tab[,4]) < 2e-16, "<2e-16") 
tab
}

.tidyParametricTable <- function(tab) {
for (i in 1:3) tab[,i] <- round(tab[,i], 2)
tab[,4] <- signif(tab[,4], 3)
tab[,4] <- replace(unlist(tab[,4]), unlist(tab[,4]) < 2e-16, "<2e-16") 
tab
}

.smooth.summary.evgam <- function(object) {
idgamlist <- sapply(object, class) == "gamlist" 
npar <- sum(idgamlist)
out <- lapply(object$gotsmooth, function(i) do.call(rbind, lapply(seq_along(object[[i]]$smooth), function(j) .giveSmoothTest(object[[i]], j))))
for (i in seq_along(out)) names(out[[i]]) <- c("edf", "max.df", "Chi.sq", "Pr(>|t|)")
nms <- names(object)[object$gotsmooth]
nms[nms == "mu"] <- "location"
nms[nms == "lpsi"] <- "logscale"
nms[nms == "xi"] <- "shape"
names(out) <- nms
out
}

.parametric.summary.evgam <- function(object) {
idgamlist <- sapply(object, class) == "gamlist" 
npar <- sum(idgamlist)
out <- lapply(which(idgamlist), function(i) .giveParametricTest(object[[i]]))
for (i in seq_along(out)) names(out[[i]]) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
nms <- names(object)[idgamlist]
nms[nms == "mu"] <- "location"
nms[nms == "lpsi"] <- "logscale"
nms[nms == "xi"] <- "shape"
names(out) <- nms
out
}

## Newton functions

.zerosinvec <- function(x, id) {
out <- numeric(length(id))
out[id] <- x
out
}

.zerosinmat <- function(x, id) {
out <- matrix(0, length(id), length(id))
out[id, id] <- x
out
}

.newton_step_inner <- function(pars, fn, sfn, ..., control, trace=0, kept=NULL, alpha0=1) {

steptol <- control$steptol
itlim <- control$itlim
fntol <- control$fntol
gradtol <- control$gradtol
stepmax <- control$stepmax

pars0 <- pars
if (is.null(kept)) 
  kept <- !logical(length(pars))

it <- 1
okay <- TRUE
f0 <- fn(pars, ...)
step1 <- NULL

while (okay) {
  if (it > 1) 
    g0 <- g
  if (!is.null(step1)) {
    step0 <- step1
    g <- attr(step0, "gradient")
  } else {
    attr(pars, "beta") <- attr(f0, "beta")
    step0 <- sfn(pars, ..., kept=kept)
    g <- attr(step0, "gradient")
  }
  if (trace) 
    .itreport(f0, g, it - 1)
  if (mean(abs(g)) < gradtol) {
    report <- c("gradient tolerance reached")
    break
  }

  step0 <- sign(step0) * pmin(abs(step0), stepmax)
  alpha <- alpha0
  report <- NULL
  ls <- TRUE
  while(ls & is.null(report)) {
    step <- alpha * step0
    stepokay <- mean(abs(step)) > steptol
    if (!stepokay) {
      report <- c("step tolerance reached")
    } else {
      theta1 <- pars - step
      f1 <- fn(theta1, ...)
      d <- f1 - f0
      if (!is.finite(d)) 
        d <- 10
      if (d < 0) {
        attr(theta1, "beta") <- attr(f1, "beta")
        step1 <- try(sfn(theta1, ..., kept=kept), silent=TRUE)
      if (inherits(step1, "try-error")) 
        d <- 1
      if (any(!is.finite(attr(step1, "gradient")))) 
        d <- 1
      }  
      if (d < 0) {
        f0 <- f1
        pars <- theta1
        ls <- FALSE
      } else {
        if (d < fntol) 
          report <- c("function tolerance reached")
        alpha <- .5 * alpha
      }
    }
  }
  if (!is.null(report)) 
    break
  it <- it + 1
  if (it == itlim) {
    report <- c("iteration limit reached")
    okay <- FALSE
  }
}

if (trace) 
  cat(paste("\n ", it, "iterations:", report, "\n"))
out <- list(pars=as.vector(pars), objective=f0)
out$gradient <- attr(step0, "gradient")
out$Hessian <- attr(step0, "Hessian")
# if (!is.null(attr(step0, "PP"))) {
#   drop <- .rank2drop(attr(step0, "PP"))
#   if (length(drop) > 0) {
#     if (any(!kept)) browser()
#     kept[which(kept)[drop]] <- FALSE
#   }
# }
if (!is.null(attr(step0, "PP")))
  kept <- .new.kept(attr(step0, "PP"), kept)
out$kept <- kept  
out$cholHessian <- attr(step0, "cholH")
out$diagHessian <- attr(step0, "diagH")
out$rankHessian <- attr(step0, "rank")
out$convergence <- 0
out$report <- report
out$iterations <- it
out$gradconv <- substr(report, 1, 4) == "grad"
if (!is.null(attr(pars, "beta"))) 
  out$beta <- attr(pars, "beta")
out
}

.rank2drop <- function(x) {
nc <- ncol(x)
R <- qr(x)
r <- R$rank
if (r < nc) {
  drop <- R$pivot[-seq_len(r)]
} else {
  drop <- integer(0)
}
drop
}

.new.kept <- function(x, kept) {
x <- x[kept, kept, drop=FALSE]
nc <- ncol(x)
R <- qr(x)
r <- R$rank
if (r < nc) {
  drop <- R$pivot[-seq_len(r)]
} else {
  drop <- integer(0)
}
if (length(drop) > 0)
  kept[which(kept)[drop]] <- FALSE
kept
}

.newton_step <- function(pars, fn, sfn, ..., control, trace=0, alpha0=1) {

nkept <- length(pars)
fit0 <- .newton_step_inner(pars, fn, sfn, ..., control=control, trace=trace, alpha0=alpha0)
nkept1 <- sum(fit0$kept)

while(nkept > nkept1) {
  nkept <- nkept1
  # fit0$par[!fit0$kept] <- 0
  fit0 <- .newton_step_inner(fit0$par, fn, sfn, ..., control=control, trace=trace, kept=fit0$kept, alpha0=alpha0)
  nkept1 <- sum(fit0$kept)
}

fit0

}

.itreport <- function(f, g, it) {
    report <- paste("\n Outer iteration ", it, ":", sep="")
    rep1 <- paste("  Outer max(|grad|):", signif(max(abs(g)), 3))
    rep2 <- paste("  Inner max(|grad|): ", signif(max(abs(attr(f, "gradient"))), 3), ".", sep="")
    report <- c(report, paste(rep1, rep2, sep="; "))
    cat(paste(report, collapse="\n"))
}

.BFGS <- function(pars, fn, gfn, ..., control, trace=0) {

steptol <- control$steptol
itlim <- control$itlim
fntol <- control$fntol
gradtol <- control$gradtol
stepmax <- control$stepmax

it <- 1
okay <- TRUE
f0 <- fn(pars, ...)
g1 <- NULL
I <- iH <- H <- diag(length(pars))

while (okay) {
if (it > 1) g0 <- g
if (!is.null(g1)) {
g <- g1
} else {
attr(pars, "beta") <- attr(f0, "beta")
g <- gfn(pars, ...)
}
if (trace) .itreport(f0, g, it - 1)
if (mean(abs(g)) < gradtol) {
report <- c("gradient tolerance reached")
break
}
step0 <- crossprod(iH, g)
step0 <- sign(step0) * pmin(abs(step0), stepmax)
alpha <- 1
report <- NULL
ls <- TRUE
while(ls & is.null(report)) {
step <- alpha * step0
stepokay <- all(abs(step) > steptol)
if (!stepokay) {
report <- c("step tolerance reached")
} else {
theta1 <- pars - step
f1 <- fn(theta1, ...)
d <- f1 - f0
if (d < 0) {
attr(theta1, "beta") <- attr(f1, "beta")
g1 <- gfn(theta1, ...)
if (any(!is.finite(g1))) d <- 1
yk <- g1 - g
denom <- sum(- yk * step)
t1 <- I - tcrossprod(- step, yk) / denom
t2 <- I - tcrossprod(yk, - step) / denom
t3 <- tcrossprod(- step) / denom
iH <- t1 %*% iH %*% t2 + t3
if (any(!is.finite(iH))) d <- 1
}
if (d < 0) {
f0 <- f1
pars <- theta1
ls <- FALSE
} else {
if (d < fntol) {
report <- c("function tolerance reached")
}
alpha <- .5 * alpha
}
}
}
if (!is.null(report)) break
it <- it + 1
if (it == itlim) {
report <- c("iteration limit reached")
okay <- FALSE
}
}
if (trace) cat(paste("\n ", it, "iterations:", report, "\n"))
out <- list(par=as.vector(pars), objective=f0)
out$gradient <- g
out$convergence <- 0
out$report <- report
out$iterations <- it
if (!is.null(attr(pars, "beta"))) out$beta <- attr(pars, "beta")
out
}
