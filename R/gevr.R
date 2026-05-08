## Generalised extreme value negative log-likelihood functions

.gevr.d0 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))
  out <- gevrd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.gevr.d12 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- gevrd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.gevr.d34 <- function(pars, likdata) {
  likdata$y <- as.matrix(likdata$y)
  nhere <- rowSums(is.finite(likdata$y))  
  out <- gevrd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, nhere)
  out
}

.gevrfns <- list(d0 = .gevr.d0, d120 = .gevr.d12, d340 = .gevr.d34)

.gevr_unlink <- list(NULL, function(x) exp(x), function(x) 1.5 / (1 + exp(-x)) - 1.0)
attr(.gevr_unlink[[2]], "deriv") <- .gevr_unlink[[2]]
attr(.gevr_unlink[[3]], "deriv") <- function(x) 1.5 * exp(-x)/(1 + exp(-x))^2

.gevrfns$q <- .qgev
.gevrfns$unlink <- .gevr_unlink

.gevrfns$initfn <- function(lst) {
  inits <- c(sqrt(6) * sd(as.vector(lst$y), na.rm = TRUE) / pi, .85)
  inits <- c(mean(lst$y, na.rm = TRUE) - .5772 * inits[1], log(inits[1]), inits[2])
  inits
}

.q_gevr <- function(p, pars1, pars2, pars3) {
  loc <- pars1
  scale <- exp(pars2)
  shape <- 1.5 / (1 + exp(-pars3)) - 1.0
  shape[shape == 0] <- 1e-6
  shape <- sign(shape) * pmax(abs(shape), 1e-6)
  yp <- -log(p)
  loc - scale * (1 - yp^(-shape)) / shape
}

# Deriv::Deriv(.q_gevr, paste('pars', 1:3, sep = ''), combine = 'cbind')

.dq_gevr <- function (p, pars1, pars2, pars3) {
  .e2 <- exp(-pars3)
  .e3 <- 1 + .e2
  .e5 <- 1.5/.e3 - 1
  .e6 <- -log(p)
  .e7 <- .e6^.e5
  .e8 <- 1 - 1/.e7
  .e9 <- exp(pars2)
  cbind(pars1 = 1, 
        pars2 = -(.e8 * .e9/.e5), 
        pars3 = -(1.5 * (.e2 * .e9 * (log(.e6)/.e7 - .e8/.e5)/(.e3^2 * .e5)))
  )
}

.gevrfns$q <- .q_gevr
.gevrfns$dq <- .dq_gevr
.gevrfns$unlink <- .gevr_unlink

.p_gevr <- function(x, pars1, pars2, pars3, log = FALSE) {
  loc <- pars1
  scale <- exp(pars2)
  shape <- 1.5 / (1 + exp(-pars3)) - 1.0
  temp <- 1 + shape * (x - loc) / scale
  temp <- pmax(temp, 0)
  out <- -(temp ^ (-1/shape))
  if (!log) 
    out <- exp(out)
  out
}

.gevrfns$p <- .p_gevr