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

.newton_step_inner <- function(pars, fn, sfn, ..., control, trace=0, kept = NULL,
                               attr2pass = 'beta') {
#                                f2return = list(
#                                  finalnames = c('gradient', 'Hessian', 'cholHessian', 'diagHessian', 'rankHessian'),
#                                  workingnames = c('gradient', 'Hessian', 'cholH', 'diagH', 'rank')),
#                                s2return = list(finalnames = c('gradient', 'Hessian', 'cholHessian', 'diagHessian', 'rankHessian'),
#                                                workingnames = c('gradient', 'Hessian', 'cholH', 'diagH', 'rank'))
# ) {
#   
  steptol <- control$steptol
  itlim <- control$itlim
  fntol <- control$fntol
  gradtol <- control$gradtol
  dgradtol <- control$dgradtol
  stepmax <- control$stepmax
  alpha0 <- control$alpha0
  rho0 <- control$rho0
  
  if (is.null(dgradtol))
    dgradtol <- 1e-4
  
  nms <- names(pars)
  
  pars0 <- pars
  if (is.null(kept))
    kept <- !logical(length(pars))
  
  it <- 0
  okay <- TRUE
  f0 <- fn(pars, ...)
  step1 <- NULL
  evals <- 1
  attr(pars, 'first') <- FALSE
  
  while (okay) {
    it <- it + 1
    if (!is.null(attr2pass))
      attributes(pars)[attr2pass] <- attributes(f0)[attr2pass]
    if (it == 1) {
      step0 <- sfn(pars, ..., kept = kept)
      gg <- 1e6 * 1:3
    }
    # else {
    #   step0 <- step1
    # }
    g <- attr(step0, "gradient")
    gg <- c(max(abs(g)), gg)
    
    if (trace)
      .itreport(f0, g, it - 1)
    
    if (gg[1] < gradtol) {
      report <- c("gradient tolerance reached")
      break
    }
    
    if (abs(diff(range(gg[1:4]))) < dgradtol) {
      report <- c("gradient difference tolerance reached")
      break
    }
    
    biggest <- max(abs(step0))
    if (biggest > stepmax)
      step0 <- step0 * stepmax / biggest
    
    alpha <- alpha0
    report <- NULL
    ls <- TRUE
    while(ls & is.null(report)) {
      step <- alpha * step0
      stepokay <- sqrt(sum(step^2)) > steptol
      if (!stepokay) {
        report <- c("step tolerance reached")
      } else {
        theta1 <- pars - step
        if (!is.null(attr2pass))
          attributes(theta1)[attr2pass] <- attributes(f0)[attr2pass]
        names(theta1) <- nms
        f1 <- fn(theta1, ...)
        evals <- evals + 1
        d <- f1 - f0
        if (!is.finite(d))
          d <- 10
        if (d < 0) {
          if (!is.null(attr2pass))
            attributes(theta1)[attr2pass] <- attributes(f1)[attr2pass]
          step1 <- try(sfn(theta1, ..., kept=kept), silent=TRUE)
          if (inherits(step1, "try-error"))
            d <- 1
          if (any(!is.finite(attr(step1, "gradient"))) | any(!is.finite(attr(step1, "PP"))))
            d <- 1
        }
        if (d < 0) {
          f0 <- f1
          pars <- theta1
          if (!is.null(attr2pass))
            attributes(pars)[attr2pass] <- attributes(f0)[attr2pass]
          ls <- FALSE
          step0 <- step1
          if (abs(d) < fntol)
            report <- c("function tolerance reached")
        } else {
          alpha <- ifelse(alpha == 1, .1 * alpha, rho0 * alpha)
        }
      }
    }
    if (!is.null(report))
      break
    if (it == itlim) {
      report <- c("iteration limit reached")
      okay <- FALSE
    }
  }
  
  if (trace)
    cat(paste("\n ", it, "iterations:", report, "\n"))
  out <- list(par=as.vector(pars), objective=f0)
  out <- c(out, attributes(pars))
  # out$gradient <- attr(step0, "gradient")
  # out$Hessian <- attr(step0, "Hessian")
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
  # out$cholHessian <- attr(step0, "cholH")
  # out$diagHessian <- attr(step0, "diagH")
  # out$rankHessian <- attr(step0, "rank")
  out$convergence <- 0
  out$report <- report
  out$iterations <- it
  out$gradconv <- substr(report, 1, 4) == "grad"
  # if (!is.null(attr(pars, "beta")))
  #   out$beta <- attr(pars, "beta")
  if (!is.null(attr(f0, 'gradient'))) {
    out <- c(out, attributes(f0))
  } else {
    out <- c(out, attributes(step0))
  }
  # if (!is.null(attr2return))
  #   out[attr2return$finalnames] <- attributes(f0)[attr2return$workingnames]
  out$fevaluations <- evals
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

.newton_step <- function(pars, fn, sfn, ..., control, trace=0, alpha0 = NULL) {
  
  nkept <- length(pars)
  fit0 <- .newton_step_inner(pars, fn, sfn, ..., control=control, trace=trace)
  nkept1 <- sum(fit0$kept)
  
  while(nkept > nkept1) {
    nkept <- nkept1
    # fit0$par[!fit0$kept] <- 0
    fit0 <- .newton_step_inner(fit0$par, fn, sfn, ..., control=control, trace=trace, kept=fit0$kept)
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

# .BFGS <- function(pars, fn, gfn, ..., control, trace=0) {
#   
#   steptol <- control$steptol
#   itlim <- control$itlim
#   fntol <- control$fntol
#   gradtol <- control$gradtol
#   stepmax <- control$stepmax
#   dgradtol <- control$dgradtol
#   
#   it <- 0
#   okay <- TRUE
#   f0 <- fn(pars, ...)
#   g1 <- NULL
#   I <- iH <- diag(length(pars))
#   
#   while (okay) {
#     it <- it + 1
#     attr(pars, "beta") <- attr(f0, "beta")
#     if (it == 1) {
#       g <- gfn(pars, ...)
#       gg <- 1e6 * 1:3
#     } else {
#       g <- g1
#     }
#     gg <- c(max(abs(g)), gg)
#     
#     if (trace)
#       .itreport(f0, g, it - 1)
#     
#     if (gg[1] < gradtol) {
#       report <- c("gradient tolerance reached")
#       break
#     }
#     
#     if (abs(diff(range(gg[1:4]))) < dgradtol) {
#       report <- c("gradient difference tolerance reached")
#       break
#     }
#     
#     # if (it > 1) g0 <- g
#     # if (!is.null(g1)) {
#     #   g <- g1
#     # } else {
#     #   attr(pars, "beta") <- attr(f0, "beta")
#     #   g <- gfn(pars, ...)
#     # }
#     # if (trace) 
#     #   .itreport(f0, g, it - 1)
#     # if (mean(abs(g)) < gradtol) {
#     #   report <- c("gradient tolerance reached")
#     #   break
#     # }
#     step0 <- crossprod(iH, g)
#     biggest <- max(abs(step0))
#     if (biggest > stepmax)
#       step0 <- step0 * stepmax / biggest
#     # step0 <- sign(step0) * pmin(abs(step0), stepmax)
#     alpha <- 1
#     report <- NULL
#     ls <- TRUE
#     while(ls & is.null(report)) {
#       step <- alpha * step0
#       stepokay <- sqrt(sum(step^2)) > steptol
#       if (!stepokay) {
#         report <- c('step tolerance reached')
#       } else {
#         theta1 <- pars - step
#         f1 <- fn(theta1, ...)
#         d <- f1 - f0
#         if (d < 0) {
#           attr(theta1, "beta") <- attr(f1, "beta")
#           g1 <- gfn(theta1, ...)
#           if (any(!is.finite(g1))) {
#             d <- 1
#           } else {
#             yk <- g1 - g
#             denom <- sum(- yk * step)
#             if (denom <= 0 || !is.finite(denom)) {
#               d <- 1
#             } else {
#               t1 <- I - tcrossprod(- step, yk) / denom
#               t2 <- I - tcrossprod(yk, - step) / denom
#               t3 <- tcrossprod(- step) / denom
#               iH <- t1 %*% iH %*% t2 + t3
#               if (any(!is.finite(iH))) 
#                 d <- 1
#             }
#           }
#         }
#         if (d < 0) {
#           f0 <- f1
#           pars <- theta1
#           ls <- FALSE
#           g <- g1
#         } else {
#           if (abs(d) < fntol) {
#             report <- 'function tolerance reached'
#           } else {
#             alpha <- .5 * alpha
#           }
#         }
#       }
#     }
#     if (!is.null(report)) 
#       break
#     if (it == itlim) {
#       report <- 'iteration limit reached'
#       okay <- FALSE
#     }
#   }
#   if (trace) cat(paste("\n ", it, "iterations:", report, "\n"))
#   out <- list(par=as.vector(pars), objective=f0)
#   out$gradient <- g
#   out$convergence <- 0
#   out$report <- report
#   out$iterations <- it
#   if (!is.null(attr(pars, "beta"))) out$beta <- attr(pars, "beta")
#   out
# }

.BFGS <- function(pars, fn, gfn, ..., control, trace = 0, 
                  attr2pass = 'beta') {
                  # attr2return = list(finalnames = 'beta',
                  #                    workingnames = 'beta')) {
                  # 
  steptol <- control$steptol
  itlim <- control$itlim
  fntol <- control$fntol
  gradtol <- control$gradtol
  stepmax <- control$stepmax
  dgradtol <- control$dgradtol
  
  if (is.null(dgradtol))
    dgradtol <- 1e-4
  
  nms <- names(pars)
  
  it <- 0
  okay <- TRUE
  f0 <- fn(pars, ...)
  g1 <- NULL
  I <- iH <- diag(length(pars))
  attr(pars, 'first') <- FALSE
  
  while (okay) {
    it <- it + 1
    if (!is.null(attr2pass))
      attributes(pars)[attr2pass] <- attributes(f0)[attr2pass]
    if (it == 1) {
      g <- gfn(pars, ...)
      gg <- 1e6 * 1:3
    } else {
      g <- g1
    }
    gg <- c(max(abs(g)), gg)
    
    if (trace)
      .itreport(f0, g, it - 1)
    
    if (gg[1] < gradtol) {
      report <- c("gradient tolerance reached")
      break
    }
    
    if (abs(diff(range(gg[1:4]))) < dgradtol) {
      report <- c("gradient difference tolerance reached")
      break
    }
    
    # if (it > 1) g0 <- g
    # if (!is.null(g1)) {
    #   g <- g1
    # } else {
    #   attr(pars, "beta") <- attr(f0, "beta")
    #   g <- gfn(pars, ...)
    # }
    # if (trace) 
    #   .itreport(f0, g, it - 1)
    # if (mean(abs(g)) < gradtol) {
    #   report <- c("gradient tolerance reached")
    #   break
    # }
    step0 <- crossprod(iH, g)
    biggest <- max(abs(step0))
    if (biggest > stepmax)
      step0 <- step0 * stepmax / biggest
    # step0 <- sign(step0) * pmin(abs(step0), stepmax)
    alpha <- 1
    report <- NULL
    ls <- TRUE
    while(ls & is.null(report)) {
      step <- alpha * step0
      stepokay <- sqrt(sum(step^2)) > steptol
      if (!stepokay) {
        report <- c('step tolerance reached')
      } else {
        theta1 <- pars - step
        if (!is.null(attr2pass))
          attributes(theta1)[attr2pass] <- attributes(f0)[attr2pass]
        names(theta1) <- nms
        f1 <- fn(theta1, ...)
        d <- f1 - f0
        if (d < 0) {
          if (!is.null(attr2pass))
            attributes(theta1)[attr2pass] <- attributes(f1)[attr2pass]
          g1 <- gfn(theta1, ...)
          if (any(!is.finite(g1))) {
            d <- 1
          } else {
            yk <- g1 - g
            denom <- sum(- yk * step)
            if (denom <= 0 || !is.finite(denom)) {
              d <- 1
            } else {
              t1 <- I - tcrossprod(- step, yk) / denom
              t2 <- I - tcrossprod(yk, - step) / denom
              t3 <- tcrossprod(- step) / denom
              iH <- t1 %*% iH %*% t2 + t3
              if (any(!is.finite(iH))) 
                d <- 1
            }
          }
        }
        if (d < 0) {
          f0 <- f1
          pars <- theta1
          if (!is.null(attr2pass))
            attributes(pars)[attr2pass] <- attributes(f0)[attr2pass]
          ls <- FALSE
          g <- g1
        } else {
          if (abs(d) < fntol) {
            report <- 'function tolerance reached'
          } else {
            alpha <- .5 * alpha
          }
        }
      }
    }
    if (!is.null(report)) 
      break
    if (it == itlim) {
      report <- 'iteration limit reached'
      okay <- FALSE
    }
  }
  if (trace) cat(paste("\n ", it, "iterations:", report, "\n"))
  out <- list(par=as.vector(pars), objective=f0)
  out <- c(out, attributes(pars))
  out$convergence <- 0
  out$report <- report
  out$iterations <- it
  if (!is.null(attr(f0, 'gradient'))) {
    out <- c(out, attributes(f0))
  }
  out$gradient <- g
  out$Hessian <- NULL
  # if (!is.null(attr2return))
  #   out[attr2return$finalnames] <- attributes(f0)[attr2return$workingnames]
  out
}
