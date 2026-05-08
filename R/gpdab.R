## Generalised Pareto with constrained parameters negative log-likelihood functions

.gpdab.d0 <- function(pars, likdata) {
out <- gpdabd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y[,1], likdata$gpdab, likdata$dupid, likdata$duplicate)
out
}

.gpdab.d12 <- function(pars, likdata, sandwich = FALSE) {
out <- gpdabd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y[,1], likdata$gpdab, likdata$dupid, likdata$duplicate)
out
}

.gpdab.d34 <- function(pars, likdata) {
out <- gpdabd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$y[,1], likdata$gpdab, likdata$dupid, likdata$duplicate)
out
}

.gpdabfns <- list(d0=.gpdab.d0, d120=.gpdab.d12, d340=.gpdab.d34)

.gpdabfns$initfn <- function(lst) {
  inits <- c(log(mean(lst$y, na.rm = TRUE)), .05)
  if (inits[1] < lst$args$lohi[1] | inits[1] > lst$args$lohi[3])
      inits[1] <- mean(lst$args$lohi[c(1, 3)])
  if (inits[2] < lst$args$lohi[2] | inits[2] > lst$args$lohi[4])
      inits[2] <- mean(lst$args$lohi[c(2, 4)])
  inits[1] <- exp(inits[1])
  if (lst$args$lohi[4] < 0) {
    inits[2] <- lst$args$lohi[2] + .95 * (lst$args$lohi[4] - lst$args$lohi[2])
    inits[1] <- -1.1 * inits[2] * max(lst$y, na.rm = TRUE)
  }
  ab <- c(lst$args$lohi[1:2], lst$args$lohi[3:4] - lst$args$lohi[1:2])
  inits <- -log(ab[3:4] / (inits - ab[1:2]) - 1)
  attr(inits, 'ab') <- ab
  inits
}

# .qgpd <- function(p, loc, scale, shape, zeta=1, theta=1, m=1) {
# shape <- sign(shape) * pmax(abs(shape), 1e-6)
# out <- 1 - p ^ (1 / (m * theta))
# loc + scale * ((out / zeta)^(-shape) - 1) / shape
# }
# 
# .pgpd <- function(x, loc, scale, shape, tau, NAOK=FALSE, log=FALSE) {
# # function to evaluate daily cdf, e.g. Coles (2001, pp.138)
# below <- x < loc
# x <- pmax(x, loc)
# temp <- 1 + shape * (x - loc) / scale
# if (!NAOK) temp <- pmax(temp, 0)
# out <- temp ^ (-1/shape)
# out <- 1 - (1 - tau) * out
# if (log) out <- log(out)
# if (NAOK) out[below] <- NA
# out
# }
# 
# .dqgpd <- function(p, lscale, shape) {
# shape <- sign(shape) * pmax(abs(shape), 1e-6)
# .e1 <- 1 - p
# .e2 <- .e1^shape
# .e4 <- 1/.e2 - 1
# .e5 <- exp(lscale)
# d1 <- .e4 * .e5/shape
# d2 <- -((.e4/shape + log(.e1)/.e2) * .e5/shape)
# cbind(d1, d2)
# }

.q_gpdab <- function(p, pars1, pars2, a1, a2, b1, b2) {
  scale <- 1 / (1 + exp(-pars1))
  scale <- a1 + (b1 - a1) * scale
  shape <- 1 / (1 + exp(-pars2))
  shape <- a2 + (b2 - a2) * shape
  scale * ((1 - p)^(-shape) - 1)/shape
}

# Deriv::Deriv(.q_gpdab, paste('pars', 1:2, sep = ''), combine = 'cbind')

.dq_gpdab <- function(p, pars1, pars2, a1, a2, b1, b2) {
  .e2 <- exp(-pars2)
  .e3 <- 1 + .e2
  .e4 <- b2 - a2
  .e6 <- .e4/.e3 + a2
  .e7 <- 1 - p
  .e9 <- .e7^.e6
  .e10 <- exp(-pars1)
  .e11 <- 1 + .e10
  .e13 <- 1/.e9 - 1
  .e14 <- b1 - a1
  cbind(pars1 = .e13 * .e14 * .e10/(.e6 * .e11^2), 
        pars2 = -((.e13/.e6 + log(.e7)/.e9) * (.e14/.e11 + a1) * .e4 * .e2/(.e6 * .e3^2))
        )
}

.gpdab_unlink <- list(function(x, a1, b1) a1 + (b1 - a1) / (1 + exp(-x)), 
                      function(x, a2, b2) a2 + (b2 - a2) / (1 + exp(-x)))
attr(.gpdab_unlink[[1]], "deriv") <- function(x, a1, b1) (b1 - a1) * exp(-x)/(1 + exp(-x))^2
attr(.gpdab_unlink[[2]], "deriv") <- function(x, a2, b2) (b2 - a2) * exp(-x)/(1 + exp(-x))^2

.gpdabfns$unlink <- .gpdab_unlink
.gpdabfns$q <- .q_gpdab
.gpdabfns$dq <- .dq_gpdab

.p_gpdab <- function(x, pars1, pars2, a1, a2, b1, b2, log = FALSE) {
  scale <- 1 / (1 + exp(-pars1))
  scale <- a1 + (b1 - a1) * scale
  shape <- 1 / (1 + exp(-pars2))
  shape <- a2 + (b2 - a2) * shape
  temp <- 1 + shape * x / scale
  out <- temp ^ (-1/shape)
  out <- 1 - out
  if (log) 
    out <- log(out)
  out
}

.gpdabfns$p <- .p_gpdab

