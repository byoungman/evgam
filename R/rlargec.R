## r-largest generalised extreme value negative log-likelihood functions

.rlargec.d0 <- function(pars, likdata) {
likdata$y <- as.matrix(likdata$y)
out <- rlargecd0(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, likdata$args$drop)
out
}

.rlargec.d12 <- function(pars, likdata) {
likdata$y <- as.matrix(likdata$y)
out <- rlargecd12(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, likdata$args$drop)
out
}

.rlargec.d34 <- function(pars, likdata) {
likdata$y <- as.matrix(likdata$y)
out <- rlargecd34(split(pars, likdata$idpars), likdata$X[[1]], likdata$X[[2]], likdata$X[[3]], likdata$y, likdata$dupid, likdata$duplicate, likdata$args$drop)
out
}

.rlargecfns <- list(d0 = .rlargec.d0, d120 = .rlargec.d12, d340 = .rlargec.d34)

.rlargecfns$q <- .qgev
.rlargecfns$unlink <- .rlarge_unlink

