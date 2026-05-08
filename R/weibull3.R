## Three-parameter Weibull negative log-likelihood functions

.weib3.d0 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))
  out <- weib3d0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.weib3.d12 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- weib3d12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.weib3.d34 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- weib3d34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.weib3fns <- list(d0 = .weib3.d0, d120 = .weib3.d12, d340 = .weib3.d34)

.weib3_unlink <- list(NULL, function(x) exp(x), function(x) exp(x))
attr(.weib3_unlink[[2]], "deriv") <- .weib3_unlink[[2]]
attr(.weib3_unlink[[3]], "deriv") <- .weib3_unlink[[3]]

.weib3fns$unlink <- .weib3_unlink

.weib3fns$initfn <- function(lst) {
  yy <- lst$y
  inits <- min(yy, na.rm = TRUE) - .1
  yy <- yy - inits
  inits <- c(inits, log(mean(yy, na.rm = TRUE)), .5)
  inits
}

.q_weib3 <- function(p, pars1, pars2, pars3) {
  scale <- exp(pars2)
  shape <- exp(pars3)
  pars1 + scale * ((-log(1 - p)) ^ (1/shape))
}

# Deriv::Deriv(.q_weib3, paste('pars', 1:3, sep = ''), combine = 'cbind')

.dq_weib3 <- function(p, pars1, pars2, pars3) {
  .e1 <- -log(1 - p)
  .e2 <- exp(pars3)
  .e4 <- .e1^(1/.e2) * exp(pars2)
  cbind(pars1 = 1, 
        pars2 = .e4, 
        pars3 = -(.e4 * log(.e1)/.e2)
        )
}

.weib3fns$q <- .q_weib3
.weib3fns$dq <- .dq_weib3

.p_weib3 <- function(x, pars1, pars2, pars3, log = FALSE) {
  scale <- exp(pars2)
  shape <- exp(pars3)
  out <- 1 - exp(-((x - pars1) / scale)^shape)
  if (log)
    out <- log(out)
  out
}

.weib3fns$p <- .p_weib3

