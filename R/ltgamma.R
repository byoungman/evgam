## Left truncated gamma negative log-likelihood functions

.ltgamma.d0 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))
  left <- likdata$args$left
  if (length(left) == 1)
    left <- array(left, dim(likdata$y))
  out <- ltgammad0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, as.matrix(left))
  if (!is.finite(out))
    out <- 1e20
  out
}

.ltgamma.d12 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  left <- likdata$args$left
  if (length(left) == 1)
    left <- array(left, dim(likdata$y))
  out <- ltgammad12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, as.matrix(left))
  out
}

.ltgamma.d34 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  left <- likdata$args$left
  if (length(left) == 1)
    left <- array(left, dim(likdata$y))
  out <- ltgammad34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y, likdata$dupid, likdata$duplicate, nhere, as.matrix(left))
  out
}

# .ltgammafns <- list(d0 = .ltgamma.d0, d120 = .ltgamma.d12, d340 = NULL)
.ltgammafns <- list(d0 = .ltgamma.d0, d120 = .ltgamma.d12, d340 = .ltgamma.d34)

.ltgamma_unlink <- list(function(x) exp(x), function(x) exp(x))
attr(.ltgamma_unlink[[1]], "deriv") <- .ltgamma_unlink[[1]]
attr(.ltgamma_unlink[[2]], "deriv") <- .ltgamma_unlink[[2]]

# .ltgammafns$q <- qgamma
.ltgammafns$unlink <- .ltgamma_unlink

.q_ltgamma <- function(p, pars1, pars2, left) {
  alpha <- exp(pars1)
  beta <- exp(pars2)
  F_L <- pgamma(left, shape = alpha, rate = beta)
  p_adj <- F_L + p * (1 - F_L)
  out <- qgamma(p_adj, shape = alpha, rate = beta)
  out
}

.ltgammafns$q <- .q_ltgamma

.ltgammafns$initfn <- function(lst) {
  ybar <- mean(lst$y - lst$args$left, na.rm = TRUE)
  vbar <- var(as.vector(lst$y - lst$args$left), na.rm = TRUE)
#  inits <- ybar / vbar
#  inits <- c(log(ybar * inits), log(inits))
  inits <- c(0, -log(ybar))
  inits
}

.p_ltgamma <- function(x, pars1, pars2, left, log = FALSE) {
  alpha <- exp(pars1)
  beta <- exp(pars2)
  F_y <- pgamma(x, shape = alpha, rate = beta)
  F_L <- pgamma(left, shape = alpha, rate = beta)
  out <- (F_y - F_L) / (1 - F_L)
  out[x < left] <- 0
  if (log) 
    out <- log(out)
  out
}

.ltgammafns$p <- .p_ltgamma

