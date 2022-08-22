# # #' 
# # #' @export
# # #' 
# # evgam <- function(formula, data, family="gev", aggregated = NULL, correctV=TRUE, 
# # rho0=0, inits=NULL, outer="bfgs", control=NULL, removeData=FALSE, trace=0, 
# # knots=NULL, maxdata=1e20, maxspline=1e20, compact=FALSE, 
# # ald.args=list(), exi.args=list(), pp.args=list(), sandwich.args=list()) {
# # 
# # ## setup family
# # if (is.null(aggregated)) {
# #   family.info <- .setup.family(family, pp.args)
# # } else {
# # ## aggregation stuff
# #   family.info <- .setup.family.agg(family, pp.args)
# # }
# # 
# # if (is.null(family.info$lik.fns$d340))
# #   outer <- "fd"
# # 
# # ## setup formulae
# # formula <- .setup.formulae(formula, family.info$npar, family.info$npar2, data, trace)
# # response.name <- attr(formula, "response.name")
# # 
# # ## setup mgcv objects and data
# # temp.data <- .setup.data(data, response.name, formula, family, family.info$nms, 
# #   removeData, exi.args, ald.args, pp.args, knots, maxdata, 
# #   maxspline, compact, sandwich.args, tolower(outer), trace)
# # data <- temp.data$data
# # 
# # if (!is.null(aggregated)) {
# #   agg <- list()
# #   agg$agg <- aggregated$id
# #   agg$X <- .X.evgam(temp.data$gams, aggregated$data)
# #   agg$nxy <- aggregated$nagg
# #   temp.data$lik.data$agg <- agg
# # }
# # 
# # ## initialise inner iteration
# # beta <- .setup.inner.inits(inits, temp.data$lik.data, family.info$lik.fns, family.info$npar, family)
# # lik.data <- .sandwich(temp.data$lik.data, beta)
# # if (trace > 0 & lik.data$adjust > 0) 
# #   cat(paste("\n Sandwich correct lambda =", signif(lik.data$k, 3), "\n"))
# # 
# # ## check whether any smoothing parameters need estimating
# # 
# # smooths <- length(temp.data$gotsmooth) > 0
# # 
# # if (smooths) {
# # 
# # ## initialise outer iteration
# # S.data <- .joinSmooth(temp.data$gams)
# # nsp <- length(attr(S.data, "Sl"))
# # if (is.null(rho0)) {
# #     diagSl <- sapply(attr(S.data, "Sl"), diag)
# #     rho0 <- apply(diagSl, 2, function(y) uniroot(.guess, c(-1e2, 1e2), d=attr(beta, "diagH"), s=y)$root)
# # } else {
# #     if (length(rho0) == 1) rho0 <- rep(rho0, nsp)
# # }
# # 
# # lik.data$S <- .makeS(S.data, exp(rho0))
# # 
# # ## perform outer iteration
# # fit.reml <- .outer(rho0, beta, family.info$lik.fns, lik.data, S.data, control, correctV, lik.data$outer, trace)
# # 
# # sp <- exp(fit.reml$par)
# # lik.data$S <- .makeS(S.data, sp)
# # 
# # } else {
# # 
# # S.data <- NULL
# # fit.reml <- .outer.nosmooth(beta, family.info$lik.fns, lik.data, control, trace)
# # 
# # }
# # 
# # ## covariance matrices
# # VpVc <- .VpVc(fit.reml, family.info$lik.fns, lik.data, S.data, correctV=correctV, sandwich=temp.data$sandwich, smooths=smooths, trace=trace)
# # 
# # ## effective degrees of freedom
# # edf <- .edf(fit.reml$beta, family.info$lik.fns, lik.data, VpVc, temp.data$sandwich)
# # 
# # ## update mgcv objects
# # names(temp.data$gams) <- family.info$nms
# # gams <- .swap(fit.reml, temp.data$gams, lik.data, VpVc, temp.data$gotsmooth, edf, smooths)
# # 
# # ## add extra things that make an evgam object
# # ## differ from a list of mgcv objects
# # 
# # gams <- .finalise(gams, data, family.info$lik.fns, lik.data, S.data, fit.reml, VpVc, family, temp.data$gotsmooth, formula, response.name, removeData, edf)
# # 
# # return(gams)
# # }
# 
# .setup.family.agg <- function(family, pp) {
# lik.fns <- .gevaggfns
# npar <- 4
# nms <- c("mu", "lpsi", "xi", "ltheta")
# out <- list(npar=npar, npar2=npar, lik.fns=lik.fns, nms=nms)
# }
# 
# .setup.inner.inits <- function(inits, likdata, likfns, npar, family) {
# 
# likdata0 <- likdata
# likdata0$X <- lapply(seq_along(likdata$X), function(i) matrix(1, nrow=nrow(likdata$X[[i]]), ncol=1))
# if (!is.null(likdata$agg)) 
#   likdata0$agg$X <- lapply(seq_along(likdata$agg$X), function(i) matrix(1, nrow=nrow(likdata$agg$X[[i]]), ncol=1))
# # likdata0$pp$X <- lapply(seq_along(likdata$pp$X), function(x) matrix(1, nrow=nrow(likdata$pp$X[[x]]), ncol=1))
# likdata0$S <- diag(0, npar)
# likdata0$idpars <- seq_len(npar)
# 
# if (is.null(inits)) {
#   if (npar == 1) 
#     inits <- 2
#   if (npar == 2) {
#     if (family == "ald") {
#       inits <- c(quantile(likdata0$y[,1], likdata0$tau), log(sd(likdata0$y[,1])))
#     } else {
#       inits <- c(log(mean(likdata$y[,1])), .05)
#       if (family == "transxigpd") 
#         inits[2] <- .9
#     }
#   }
#   if (npar %in% 3:4) {
#     inits <- c(sqrt(6) * sd(likdata0$y[,1]) / pi, .05)
#     inits <- c(mean(likdata0$y[,1]) - .5772 * inits[1], log(inits[1]), inits[2])
#     if (npar == 4) 
#       inits <- c(inits, 0)
#   }
#   if (npar == 6) {
#     inits <- c(sqrt(6) * sd(likdata0$y[,1]) / pi, .05)
#     inits <- c(mean(likdata0$y[,1]) - .5772 * inits[1], log(inits[1]), inits[2])
#     inits <- c(inits, 0, 0, 1)
#   }
#   likdata0$CH <- diag(length(inits))
#   likdata0$compmode <- numeric(length(inits))
#   if (!is.null(likdata$agg))
#     inits[3] <- -.663
#   test.finite <- .nllh.nopen(inits, likdata0, likfns)
#   if (npar >= 3 & test.finite == 1e20)
#     inits[1] <- min(likdata$y[,1]) - .1 * diff(range(likdata$y[,1]))
#   beta0 <- .newton_step(inits, .nllh.nopen, .search.nopen, likdata=likdata0, likfns=likfns, control=likdata$control$inner)$par
# } else {
#   if (is.list(inits)) {
#     betamat <- expand.grid(inits)
#     betanllh <- numeric(nrow(betamat))
#     for (i in seq_len(nrow(betamat))) {
#       beta0 <- unlist(betamat[i,])
#       betanllh[i] <- likfns$nllh(beta0, likdata0)
#     }
#     beta0 <- betamat[which.min(betanllh),]
#     print(beta0)
#   } else {
#     beta0 <- inits
#   }
# }
# beta0 <- unlist(lapply(seq_len(npar), function(i) c(beta0[i], rep(0, ncol(likdata$X[[i]]) - 1))))
# compmode <- 0 * beta0
# CH <- diag(compmode + 1)
# k <- 1
# likdata[c("k", "CH", "compmode")] <- list(k, CH, compmode)
# diagH <- diag(.gH.nopen(beta0, likdata=likdata, likfns=likfns)[[2]])
# if (likdata$sandwich) {
#   if (is.null(likdata$sandwich.lambda)) {
#     beta0 <- .newton_step(beta0, .nllh.nopen, .search.nopen, likdata=likdata, likfns=likfns, control=likdata$control$inner)
#     kept <- beta0$kept
# #     if (any(max(abs(beta0$gradient)) > 1)) {
# #       beta0 <- nlm(.nllh.nopen.nlm, beta0$par, likdata = likdata, likfns = likfns, check.analyticals = FALSE, print.level = 1, iterlim = 1e3)
# #       beta0 <- nlminb(beta0$estimate, .nllh.nopen, .grad.nopen, .hess.nopen, likdata = likdata, likfns = likfns)
# # #     g <- .gH.nopen(beta0$par, likdata=likdata, likfns=likfns)[[1]]
# # #     if (any(max(abs(g)) > 1)) {
# # #      
# # #       beta0 <- nlminb(beta0$estimate, .nllh.nopen, .grad.nopen, .hess.nopen, likdata = likdata, likfns = likfns)
# # #     }
# #     }
#     beta0 <- beta0$par
#     H <- .gH.nopen(beta0, likdata=likdata, likfns=likfns, sandwich=TRUE)
#     g <- rowSums(H[[1]])
#     if (family == "pp") {
#       J0 <- H[[1]]
#       J <- J0[,!likdata$ppq]
#       J0 <- rowSums(J0[,likdata$ppq])
#       J <- split(as.data.frame(t(J)), likdata$sandwich.split)
#       wts <- sapply(J, nrow)
#       wts <- wts / sum(wts)
#       J <- sapply(J, colSums)
#       J <- J + J0 %o% wts
#       J <- tcrossprod(J)
#     } else {
#       J <- split(as.data.frame(t(H[[1]])), likdata$sandwich.split)
#       J <- sapply(J, colSums)
#       J <- tcrossprod(J - rowSums(J))
#     }
#     H <- H[[2]]
#     diagH <- diag(H)
#     cholH <- try(chol(H), silent=TRUE)
#     if (inherits(cholH, "try-error")) {
#       if (!likdata$force) {
#         stop("Hessian of unpenalised MLE not positive definite.\n  Supply `force=TRUE' to `sandwich.args' to perturb it to be positive definite.")
#       } else {
#         if (trace >= 0)
#           message("Hessian perturbed to be positive definite for sandwich adjustment.")
#         iH <- pinv(H)
#       }
#     } else {
#       iH <- chol2inv(cholH)
#     }
#     if (likdata$adjust == 2) {
#       cholJ <- try(chol(J), silent=TRUE)
#       if (inherits(cholJ, "try-error") & likdata$adjust == 2) {
#         HA <- crossprod(backsolve(cholJ, H, transpose=TRUE))
#       } else {
#         iHA <- tcrossprod(crossprod(iH, J), iH)
#         choliHA <- try(chol(iHA), silent=TRUE)
#         if (inherits(choliHA, "try-error")) {
#           if (!likdata$force) {
#             stop("Sandwich variance not positive definite.\n  Supply `force=TRUE' to `sandwich.args' to perturb it to be positive definite.")
#           } else {
#             if (trace >= 0)
#               message("Sandwich variance perturbed to be positive definite.")
#             HA <- pinv(iHA)
#           }
#         } else {
#           HA <- chol2inv(choliHA)
#         }
#       }
#       sH <- svd(H)
#       M <- sqrt(sH$d) * t(sH$v)
#       sHA <- svd(HA)
#       MA <- sqrt(sHA$d) * t(sHA$v)
#       CH <- solve(M, MA)
#       compmode <- beta0
#     } else {
#       k <- 1 / mean(diag(crossprod(iH, J)))
#     }
#   } else {
#   k <- likdata$sandwich.lambda
#   g <- NA
#   }
# }
# attr(beta0, "k") <- k
# attr(beta0, "CH") <- CH
# attr(beta0, "compmode") <- compmode
# attr(beta0, "diagH") <- diagH
# if (likdata$sandwich)
#   attr(beta0, "inner_grad") <- g
# beta0
# }
# 
# 
# # .outer <- function(rho0, beta, likfns, likdata, Sdata, control, correctV, outer, trace) {
# # 
# # attr(rho0, "beta") <- beta
# # 
# # if (outer == "newton") {
# #   fit.reml <- .newton_step_inner(rho0, .reml0, .search.reml, likfns=likfns, likdata=likdata, Sdata=Sdata, control=likdata$control$outer, trace=trace > 1)
# # } else {
# #   if (outer == "fd") {
# #     fit.reml <- .BFGS(rho0, .reml0, .reml1.fd, likfns=likfns, likdata=likdata, Sdata=Sdata, control=likdata$control$outer, trace=trace > 1)
# #   } else {
# #     fit.reml <- .BFGS(rho0, .reml0, .reml1, likfns=likfns, likdata=likdata, Sdata=Sdata, control=likdata$control$outer, trace=trace > 1)
# #   }
# #   rho1 <- fit.reml$par
# #   attr(rho1, "beta") <- fit.reml$beta
# #   if (correctV) {
# #     fit.reml$Hessian <- try(.reml12(rho1, likfns=likfns, likdata=likdata, Sdata=Sdata)[[2]], silent=TRUE)
# #     if (inherits(fit.reml$Hessian, "try-error")) 
# #       fit.reml$Hessian <- try(.reml2.fd(rho1, likfns=likfns, likdata=likdata, Sdata=Sdata), silent=TRUE)
# #     if (inherits(fit.reml$Hessian, "try-error")) 
# #       fit.reml$Hessian <- .reml2.fdfd(rho1, likfns=likfns, likdata=likdata, Sdata=Sdata)
# #   }
# # }
# # 
# # if (correctV)
# #   fit.reml$invHessian <- .solve_evgam(fit.reml$Hessian)
# # 
# # fit.reml$trace <- trace
# # 
# # if (trace == 1) {
# #   report <- "\n Final max(|grad|))"
# #   likdata$S <- .makeS(Sdata, exp(fit.reml$par))
# #   report <- c(report, paste("   Inner:", signif(max(abs(.gH.pen(fit.reml$beta, likdata, likfns)[[1]])), 3)))
# #   report <- c(report, paste("   Outer:", signif(max(abs(fit.reml$gradient)), 3)))
# #   report <- c(report, "", "")
# #   cat(paste(report, collapse="\n"))
# # }
# # 
# # fit.reml
# # 
# # }
# # 
# # .reml2.fdfd <- function(pars, likfns, likdata, Sdata, H=NULL, beta=NULL, kept=NULL) {
# # beta <- attr(pars, "beta")
# # eps <- 1e-4
# # f0 <- .reml1.fd(pars, likfns, likdata, Sdata, beta=beta, H=H)
# # f1 <- matrix(0, length(pars), length(pars))
# # for (i in seq_along(pars)) {
# #   parsi <- pars
# #   parsi[i] <- parsi[i] + eps
# #   f1[i,] <- .reml1.fd(parsi, likfns, likdata, Sdata, beta=beta, H=H)
# # }
# # f1 <- (f1 - f0) / eps
# # .5 * (f1 + t(f1))
# # }
# # 
# 
# ## Aggregated generalised extreme value negative log-likelihood functions
# 
# 
# # .gevagg.d00 <- function(pars, likdata) {
# # pars <- split(pars, likdata$idpars)
# # pmat <- mapply("%*%", likdata$X, pars)
# # # identify which are aggregated
# # agg <- likdata$agg$agg
# # # point-referenced data
# # out <- evgam:::ldgev(likdata$y[!agg, 1], pmat[!agg, 1], pmat[!agg, 2], pmat[!agg, 3])
# # if (!is.finite(out))
# #   return(1e20)
# # # # aggregated data
# # pdfr0 <- mapply("%*%", likdata$agg$X[1:2], pars[1:2])
# # pdfr0[,2] <- exp(pdfr0[,2])
# # # browser()
# # # pdfr <- matrix(NA, sum(agg), 2)
# # # for (i in 1:2) {
# # #   likdata$agg$rep_index$y <- pdfr0[, i]
# # #   pdfr[, i] <- aggregate(y ~ id, likdata$agg$rep_index, mean)$y
# # # }
# # # if (length(likdata$agg$nxy) > 1) {
# # #   pdfr <- rowsum(likdata$agg$weights * pdfr0, likdata$agg$index)
# # # } else {
# # #   pdfr <- mean_arr(pdfr0, likdata$agg$nxy)
# # # }
# # pdfr <- ragged_mean_mat(pdfr0, likdata$agg$nxy)
# # # # pdfr <- apply(pdfr, 2, ragged_mean, id1 = likdata$agg$rep_index, nid = sum(agg))
# # out <- out + evgam:::ldgevagg(likdata$y[agg, 1], pdfr[,1], log(pdfr[,2]), pmat[agg, 3], pmat[agg, 4])
# # if (!is.finite(out))
# #   return(1e20)
# # pmin(out, 1e20)
# # }
# 
# .gevagg.d0 <- function(pars, likdata) {
# pars <- split(pars, likdata$idpars)
# pmat <- mapply("%*%", likdata$X, pars)
# # identify which are aggregated
# agg <- likdata$agg$agg
# # point-referenced data
# out <- evgam:::ldgev(likdata$y[!agg, 1], pmat[!agg, 1], pmat[!agg, 2], pmat[!agg, 3])
# if (!is.finite(out))
#   return(1e20)
# # # aggregated data
# pdfr0 <- mapply("%*%", likdata$agg$X[1:2], pars[1:2])
# pdfr0[,2] <- exp(pdfr0[,2])
# pdfr <- ragged_mean_mat(pdfr0, likdata$agg$nxy)
# out <- out + evgam:::ldgevagg(likdata$y[agg, 1], pdfr[,1], log(pdfr[,2]), pmat[agg, 3], pmat[agg, 4])
# if (!is.finite(out))
#   return(1e20)
# pmin(out, 1e20)
# }
# 
# # .gevagg.d12 <- function(pars, likdata) {
# # pars <- split(pars, likdata$idpars)
# # pmat <- mapply("%*%", likdata$X, pars)
# # # identify which are aggregated
# # agg <- likdata$agg$agg
# # # point-referenced data
# # d12 <- evgam:::ldgev12(likdata$y[!agg, 1], pmat[!agg, 1], pmat[!agg, 2], pmat[!agg, 3])
# # d1 <- unlist(lapply(1:4, function(i) evgam:::.sum_col(likdata$X[[i]][!agg, , drop=FALSE] * d12[,i])))
# # d2 <- diag(0 * likdata$idpars, length(likdata$idpars))
# # id <- evgam:::.indices(4)[[2]]
# # for (i in 1:4) {
# #   for (j in 1:i) {
# #     temp <- crossprod(likdata$X[[i]][!agg, , drop=FALSE], likdata$X[[j]][!agg, , drop=FALSE] * d12[,id[i, j]])
# #     d2[likdata$idpars == i, likdata$idpars == j] <- temp
# #   }
# # }
# # # aggregated data
# # pdfr0 <- mapply("%*%", likdata$agg$X[1:2], pars[1:2])
# # pdfr0[,2] <- exp(pdfr0[,2])
# # # if (length(likdata$agg$nxy) > 1) {
# # #   pdfr <- rowsum(likdata$agg$weights * pdfr0, likdata$agg$index)
# # # } else {
# # #   pdfr <- mean_arr(pdfr0, likdata$agg$nxy)
# # # }
# # pdfr <- ragged_mean_mat(pdfr0, likdata$agg$nxy)
# # d12 <- evgam:::ldgevagg12(likdata$y[agg, 1], pdfr[,1], log(pdfr[,2]), pmat[agg, 3], pmat[agg, 4])
# # # if (length(likdata$agg$nxy) > 1) {
# # #   dmb_dm <- rowsum(likdata$agg$weights * likdata$agg$X[[1]], likdata$agg$index)
# # # } else {
# # #   dmb_dm <- mean_arr(likdata$agg$X[[1]], likdata$agg$nxy)
# # # }
# # dmb_dm <- ragged_mean_mat(likdata$agg$X[[1]], likdata$agg$nxy)
# # d1_1 <- evgam:::.sum_col(d12[,1] * dmb_dm)
# # # if (length(likdata$agg$nxy) > 1) {
# # #   dpb_dp <- rowsum(pdfr0[,2] * likdata$agg$weights * likdata$agg$X[[2]], likdata$agg$index)
# # # } else {
# # #   dpb_dp <- mean_arr(pdfr0[,2] * likdata$agg$X[[2]], likdata$agg$nxy)
# # # }
# # dpb_dp <- ragged_mean_mat(pdfr0[,2] * likdata$agg$X[[2]], likdata$agg$nxy)
# # d1_2 <- evgam:::.sum_col((d12[,2] / pdfr[,2]) * dpb_dp)
# # d1_34 <- unlist(lapply(3:4, function(i) evgam:::.sum_col(likdata$X[[i]][agg, , drop=FALSE] * d12[,i])))
# # d1 <- d1 + c(d1_1, d1_2, d1_34)
# # d21 <- evgam:::ldgevagg12(likdata$y[agg, 1], pdfr[,1], log(pdfr[,2]), pmat[agg, 3], pmat[agg, 4])
# # likdata$X[[1]] <- dmb_dm
# # likdata$X[[2]] <- dpb_dp / pdfr[,2]
# # for (i in 3:4)
# #   likdata$X[[i]] <- likdata$X[[i]][agg, , drop=FALSE]
# # for (i in 1:4) {
# #   idi <- likdata$idpars == i
# #   for (j in 1:i) {
# #     idj <- likdata$idpars == j
# #     d2[idi, idj] <- d2[idi, idj] + crossprod(likdata$X[[i]], d12[, id[i, j]] * likdata$X[[j]])
# #     if (j < i)
# #       d2[idj, idi] <- t(d2[idi, idj, drop=FALSE])
# #   }
# # }
# # # if (length(likdata$agg$nxy) > 1) {
# # #   t1 <- crossprod(likdata$agg$X[[2]], (pdfr0[,2] * rep(d12[,2] / pdfr[,2], likdata$agg$nxy)) * likdata$agg$X[[2]] * likdata$agg$weights)# / likdata$agg$nxy
# # # } else {
# # #   t1 <- crossprod(likdata$agg$X[[2]], (pdfr0[,2] * rep(d12[,2] / pdfr[,2], each = likdata$agg$nxy)) * likdata$agg$X[[2]]) / likdata$agg$nxy
# # # }
# # t1 <- crossprod(likdata$agg$X[[2]], (pdfr0[,2] * rep(d12[,2] / pdfr[,2], likdata$agg$nxy)) * likdata$agg$X[[2]] * likdata$agg$weights)
# # t2 <- crossprod(likdata$X[[2]], d12[,2] * likdata$X[[2]])
# # d2[likdata$idpars == 2, likdata$idpars == 2] <- d2[likdata$idpars == 2, likdata$idpars == 2] + t1 - t2
# # list(d1, d2)
# # }
# 
# .gevagg.d12 <- function(pars, likdata, sandwich = FALSE) {
# pars <- split(pars, likdata$idpars)
# pmat <- mapply("%*%", likdata$X, pars)
# # identify which are aggregated
# agg <- likdata$agg$agg
# # point-referenced data
# d12 <- evgam:::ldgev12(likdata$y[!agg, 1], pmat[!agg, 1], pmat[!agg, 2], pmat[!agg, 3])
# d1 <- do.call(cbind, lapply(1:4, function(i) likdata$X[[i]][!agg, , drop=FALSE] * d12[,i]))
# d2 <- diag(0 * likdata$idpars, length(likdata$idpars))
# id <- evgam:::.indices(4)[[2]]
# for (i in 1:4) {
#   for (j in 1:i) {
#     temp <- crossprod(likdata$X[[i]][!agg, , drop=FALSE], likdata$X[[j]][!agg, , drop=FALSE] * d12[,id[i, j]])
#     d2[likdata$idpars == i, likdata$idpars == j] <- temp
#   }
# }
# # aggregated data
# pdfr0 <- mapply("%*%", likdata$agg$X[1:2], pars[1:2])
# pdfr0[,2] <- exp(pdfr0[,2])
# pdfr <- ragged_mean_mat(pdfr0, likdata$agg$nxy)
# d12 <- evgam:::ldgevagg12(likdata$y[agg, 1], pdfr[,1], log(pdfr[,2]), pmat[agg, 3], pmat[agg, 4])
# dmb_dm <- ragged_mean_mat(likdata$agg$X[[1]], likdata$agg$nxy)
# d1_1 <- d12[,1] * dmb_dm
# dpb_dp <- ragged_mean_mat(pdfr0[,2] * likdata$agg$X[[2]], likdata$agg$nxy)
# d1_2 <- (d12[,2] / pdfr[,2]) * dpb_dp
# d1_34 <- do.call(cbind, lapply(3:4, function(i) likdata$X[[i]][agg, , drop=FALSE] * d12[,i]))
# d1 <- rbind(d1, cbind(d1_1, d1_2, d1_34))
# d21 <- evgam:::ldgevagg12(likdata$y[agg, 1], pdfr[,1], log(pdfr[,2]), pmat[agg, 3], pmat[agg, 4])
# likdata$X[[1]] <- dmb_dm
# likdata$X[[2]] <- dpb_dp / pdfr[,2]
# for (i in 3:4)
#   likdata$X[[i]] <- likdata$X[[i]][agg, , drop=FALSE]
# for (i in 1:4) {
#   idi <- likdata$idpars == i
#   for (j in 1:i) {
#     idj <- likdata$idpars == j
#     d2[idi, idj] <- d2[idi, idj] + crossprod(likdata$X[[i]], d12[, id[i, j]] * likdata$X[[j]])
#     if (j < i)
#       d2[idj, idi] <- t(d2[idi, idj, drop=FALSE])
#   }
# }
# t1 <- crossprod(likdata$agg$X[[2]], (pdfr0[,2] * rep(d12[,2] / pdfr[,2], likdata$agg$nxy)) * likdata$agg$X[[2]] * likdata$agg$weights)
# t2 <- crossprod(likdata$X[[2]], d12[,2] * likdata$X[[2]])
# d2[likdata$idpars == 2, likdata$idpars == 2] <- d2[likdata$idpars == 2, likdata$idpars == 2] + t1 - t2
# if (!sandwich)
#   d1 <- evgam:::.sum_col(d1)
# list(d1, d2)
# }
# 
# 
# # .gevagg.d0 <- function(pars, likdata) {
# # pars <- split(pars, likdata$idpars)
# # pmat <- mapply("%*%", likdata$X, pars)
# # # identify which are aggregated
# # agg <- likdata$agg$agg
# # # point-referenced data
# # out <- ldgev(likdata$y[!agg, 1], pmat[!agg, 1], pmat[!agg, 2], pmat[!agg, 3])
# # if (!is.finite(out))
# #   return(1e20)
# # # # aggregated data
# # pdfr <- mapply("%*%", likdata$agg$X[1:2], pars[1:2])
# # pdfr[,2] <- exp(pdfr[,2])
# # pdfr <- mean_arr(pdfr, likdata$agg$nxy)
# # out <- out + ldgevagg(likdata$y[agg, 1], pdfr[,1], log(pdfr[,2]), pmat[agg, 3], pmat[agg, 4])
# # if (!is.finite(out))
# #   return(1e20)
# # pmin(out, 1e20)
# # }
# # 
# # .gevagg.d12 <- function(pars, likdata) {
# # pars <- split(pars, likdata$idpars)
# # pmat <- mapply("%*%", likdata$X, pars)
# # # identify which are aggregated
# # agg <- likdata$agg$agg
# # # point-referenced data
# # d12 <- ldgev12(likdata$y[!agg, 1], pmat[!agg, 1], pmat[!agg, 2], pmat[!agg, 3])
# # d1 <- unlist(lapply(1:4, function(i) .sum_col(likdata$X[[i]][!agg, , drop=FALSE] * d12[,i])))
# # d2 <- diag(0 * likdata$idpars, length(likdata$idpars))
# # id <- .indices(4)[[2]]
# # for (i in 1:4) {
# #   for (j in 1:i) {
# #     temp <- crossprod(likdata$X[[i]][!agg, , drop=FALSE], likdata$X[[j]][!agg, , drop=FALSE] * d12[,id[i, j]])
# #     d2[likdata$idpars == i, likdata$idpars == j] <- temp
# #   }
# # }
# # # aggregated data
# # pdfr0 <- mapply("%*%", likdata$agg$X[1:2], pars[1:2])
# # pdfr0[,2] <- exp(pdfr0[,2])
# # pdfr <- mean_arr(pdfr0, likdata$agg$nxy)
# # d12 <- ldgevagg12(likdata$y[agg, 1], pdfr[,1], log(pdfr[,2]), pmat[agg, 3], pmat[agg, 4])
# # dmb_dm <- mean_arr(likdata$agg$X[[1]], likdata$agg$nxy)
# # d1_1 <- .sum_col(d12[,1] * dmb_dm)
# # dpb_dp <- mean_arr(pdfr0[,2] * likdata$agg$X[[2]], likdata$agg$nxy)
# # d1_2 <- .sum_col((d12[,2] / pdfr[,2]) * dpb_dp)
# # d1_34 <- unlist(lapply(3:4, function(i) .sum_col(likdata$X[[i]][agg, , drop=FALSE] * d12[,i])))
# # d1 <- d1 + c(d1_1, d1_2, d1_34)
# # d21 <- ldgevagg12(likdata$y[agg, 1], pdfr[,1], log(pdfr[,2]), pmat[agg, 3], pmat[agg, 4])
# # likdata$X[[1]] <- dmb_dm
# # likdata$X[[2]] <- dpb_dp / pdfr[,2]
# # for (i in 3:4)
# #   likdata$X[[i]] <- likdata$X[[i]][agg, , drop=FALSE]
# # for (i in 1:4) {
# #   idi <- likdata$idpars == i
# #   for (j in 1:i) {
# #     idj <- likdata$idpars == j
# #     d2[idi, idj] <- d2[idi, idj] + crossprod(likdata$X[[i]], d12[, id[i, j]] * likdata$X[[j]])
# #     if (j < i)
# #       d2[idj, idi] <- t(d2[idi, idj, drop=FALSE])
# #   }
# # }
# # t1 <- crossprod(likdata$agg$X[[2]], (pdfr0[,2] * rep(d12[,2] / pdfr[,2], each=likdata$agg$nxy)) * likdata$agg$X[[2]]) / likdata$agg$nxy
# # t2 <- crossprod(likdata$X[[2]], d12[,2] * likdata$X[[2]])
# # d2[likdata$idpars == 2, likdata$idpars == 2] <- d2[likdata$idpars == 2, likdata$idpars == 2] + t1 - t2
# # list(d1, d2)
# # }
# 
# .gevaggfns <- list(d0=.gevagg.d0, d120=.gevagg.d12, d340=NULL)
# 
# # .gH.nopen <- function(pars, likdata, likfns, sandwich=FALSE, deriv=2) {
# # pars <- as.vector(likdata$compmode + likdata$CH %*% (as.vector(pars) - likdata$compmode))
# # temp <- likfns$d120(pars, likdata, sandwich)
# # if (is.null(likdata$agg))
# #   temp <- .gH(temp, likdata, sandwich, deriv)
# # temp[[1]] <- likdata$k * temp[[1]]
# # temp[[1]] <- t(temp[[1]] %*% likdata$CH)
# # if (deriv > 1) {
# #   temp[[2]] <- likdata$k * temp[[2]]
# #   temp[[2]] <- crossprod(likdata$CH, temp[[2]]) %*% likdata$CH
# #   attr(temp, "PP") <- temp[[2]] / norm(temp[[2]], "F")
# # }
# # temp
# # }
# 
# # #' 
# # #' @export
# # #'
# # predict.evgam <- function(object, newdata, type="link", prob=NULL, se.fit=FALSE, marginal=TRUE, 
# # exi = FALSE, trace = 0, ...) {
# # 
# # ## a few checks
# # 
# # family <- object$family
# # 
# # if (type == "qqplot" & !(family %in% c("gev", "gpd", "weibull")))
# #   stop(paste("`type = 'qqplot'' not yet available for family '", family, "'", sep=""))
# # 
# # if (family == "pp") {
# #   family <- "gev"
# #   if (missing(newdata)) {
# #     newdata <- object$data
# #     if (trace >= 0)
# #       message("Predictions for point process model given for quadrature points, not original data frame")
# #   }
# # }
# # 
# # if (!is.null(prob)) 
# #   type <- "quantile"
# # if (type == "quantile" & is.null(prob)) 
# #   stop("non-NULL `prob' required if `type = quantile'")
# # 
# # ## end checks
# # 
# # ## a few things to set aside
# # 
# # if (se.fit) {
# # 
# # if (marginal) {
# #   V.type <- "Vc" 
# # } else {
# #   V.type <- "Vp"
# # }
# # 
# # conf.pars <- list(object$coefficients, object[[V.type]], object$idpars)
# # 
# # }
# # 
# # if (type == "qqplot")
# #   response.name <- object$response.name
# #   
# # if (family == "exi")
# #   linkfn <- object$linkfn
# #   
# # got.newdata <- !missing(newdata)
# # 
# # if (got.newdata) {
# #   pred.vars <- object$predictor.names
# #   missing.vars <- pred.vars[!(pred.vars %in% names(newdata))]
# #   if (length(missing.vars) > 0)
# #     stop(paste("Variable(s) '", paste(missing.vars, collapse=", "), "' not supplied to `newdata'.", sep=""))
# # }
# # 
# # if (got.newdata) {
# #   ndat <- nrow(newdata)
# # } else {
# #   if (is.null(object$data)) 
# #     stop("Supply `evgam' with `removeData = FALSE' if not supplying `newdata'.")
# #   ndat <- nrow(object$data)
# # }
# # 
# # ## end stuff to set aside
# # 
# # ## X creation starts
# # 
# # X <- .X.evgam(object, newdata)
# # nX <- length(X)
# # nms <- names(object)[seq_len(nX)]
# # 
# # ## X creation ends
# # 
# # if (type == "lpmatrix") {
# # 
# # return(X)
# # 
# # } else { ## start links
# # 
# # out <- lapply(seq_len(nX), function(i) X[[i]] %*% object[[i]]$coefficients)
# # names(out) <- names(X)
# # out <- as.data.frame(lapply(out, function(x) x[,1]))
# # 
# # if (se.fit) { ## start working with standard errors
# # 
# # if (type != "quantile") {
# #   std.err <- lapply(seq_len(nX), function(i) sqrt(rowSums(X[[i]] * (X[[i]] %*% object[[i]][[V.type]]))))
# #   std.err <- as.data.frame(std.err)
# #   names(std.err) <- nms
# # }
# # 
# # } ## end of working with standard errors (for now)
# # 
# # if (type %in% c("response", "quantile", "qqplot")) {
# # 
# # if (type == "qqplot") se.fit <- FALSE
# # 
# # if (family != "exi") {
# # 
# # unlink <- which(substr(nms, 1, 3) == "log")
# # for (i in unlink) {
# #   out[,i] <- exp(out[,i])
# #   if (substr(nms[i], 1, 5) == "logit") {
# #     temp <- exp(-out[,i])
# #     out[,i] <- 1 / (1 + temp)
# #     if (se.fit & type == "response")
# #       std.err[,i] <- temp * out[,i] * std.err[,i] / (1 + temp)
# #   }
# #   if (se.fit & type == "response")
# #     std.err[,i] <- out[,i] * std.err[,i]
# # }
# # 
# # if (exi & ncol(out) == 4) {
# #   out[,4] <- out[,4] ^ out[,3]
# #   out[,1] <- out[,1] - out[,2] * (1 - out[,4]) / out[,3]
# #   out[,2] <- out[,2] * out[,4]
# #   out <- out[, 1:3, drop=FALSE]
# #   nms <- nms[1:3]
# # }
# # 
# # } else {
# # 
# #   if (se.fit) 
# #     std.err[,1] <- attr(linkfn, "deriv")(out[,1]) * std.err[,1]
# #   out[,1] <- linkfn(out[,1])
# # 
# # }
# # 
# # nms <- gsub("cloglog", "", nms)
# # nms <- gsub("probit", "", nms)
# # nms <- gsub("logit", "", nms)
# # nms <- gsub("log", "", nms)
# # names(out) <- nms
# # 
# # if (type == "qqplot") { ## start qqplot
# # 
# # pit <- !all(apply(out, 2, function(x) all(diff(x) < 1e-12)))
# # x <- ppoints(ndat)
# # y <- newdata[,response.name]
# # if (is.null(y))
# #   stop("No response data.")
# #   
# # if (!pit) {
# #   if (!(family %in% c("gev", "gpd", "weibull")))
# #     stop("Unsuitable `family' for `type == 'qqplot''")
# #   if (family == "gev")
# #     x <- .qgev(x, out[,1], out[,2], out[,3])
# #   if (family == "gpd")
# #     x <- .qgpd(x, 0, out[,1], out[,2])
# #   if (family == "weibull") 
# #     x <- .qweibull(x, out[,1], out[,2])
# # } else {
# #   if (trace > 0)
# #     message("Margins converted to unit exponential by probability integral transformation.")
# #   x <- qexp(x)
# #   if (family == "gev")
# #     y <- .pgev(y, out[,1], out[,2], out[,3])
# #   if (family == "gpd")
# #     y <- .pgpd(y, 0, out[,1], out[,2])
# #   if (family == "weibull") 
# #     y <- .pweibull(y, out[,1], out[,2])
# #   y <- qexp(y)
# # }
# # 
# # qqplot(x, y)
# # abline(0, 1)
# # 
# # } else { ## end qqplot
# # 
# # if (se.fit & type == "response") 
# #   names(std.err) <- nms
# # 
# # if (type == "quantile") { ## convert response to quantile predictions
# # 
# # pars <- out
# # nprob <- length(prob)
# # out <- matrix(NA, ndat, nprob)
# # 
# # for (j in seq_len(nprob)) {
# #   if (family == "gpd") {
# #     out[, j] <- .qgpd(prob[j], 0, pars[,1], pars[,2])
# #   } else {
# #     if (family == "gev") {
# #       out[, j] <- .qgev(prob[j], pars[,1], pars[,2], pars[,3])
# #     } else {
# #       if (family == "weibull") {
# #         out[, j] <- .qweibull(prob[j], scale=pars[,1], shape=pars[,2])
# #       } else {
# #         stop("invalid family")
# #       } 
# #     }
# #   }
# # }
# # 
# # if (se.fit) { ## standard errors for quantile predictions using Delta method
# # 
# # Sigma <- array(NA, dim=c(ndat, nX, nX))
# # idp <- conf.pars[[3]]
# # for (i in seq_len(ndat)) {
# #   # rather inefficient: wants vectorising or Rcpp-ising
# #   for (j in seq_len(nX)) {
# #     for (k in j:nX) {
# #       xj <- X[[j]][i,]
# #       xk <- X[[k]][i,]
# #       V <- conf.pars[[2]][idp == j, idp == k, drop=FALSE]
# #       Sigma[i, j, k] <- sum(xj * (V %*% xk))
# #       if (k != j) 
# #         Sigma[i, k, j] <- Sigma[i, j, k]
# #     }
# #   }
# # }
# # 
# # std.err <- matrix(NA, ndat, nprob)
# # 
# # for (j in seq_len(nprob)) {
# #   if (family == "gev") {
# #     jac <- .dqgev(prob[j], pars[,1], log(pars[,2]), pars[,3])
# #   }
# #   if (family == "ggpd") {
# #     jac <- .dqgpd(prob[j], log(pars[,2]), pars[,3])
# #   }
# #   if (family == "weibull") {
# #     jac <- .dqweibull(prob[j], log(pars[,2]), pars[,3])
# #   }
# #   for (i in seq_len(ndat)) {
# #     std.err[i, j] <- sum(jac[i,] * (Sigma[i,,] %*% jac[i,]))
# #   }
# # }
# # 
# # std.err <- as.data.frame(sqrt(std.err))
# # names(std.err) <- paste("q", round(prob, 3), sep=":")
# # 
# # } ## end std.err for quantile
# # 
# # out <- as.data.frame(out)
# # names(out) <- paste("q", round(prob, 3), sep=":")
# # 
# # } ## end of response to quantile
# # 
# # if (se.fit) {
# # 
# # out <- list(fitted = out, se.fit = std.err)
# # 
# # }
# # 
# # } ## end not qqplot
# # 
# # } ## end of not link
# # 
# # if (type != "qqplot") return(out)
# # 
# # } ## end of link
# # 
# # }
