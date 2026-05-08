#' Predictions from a fitted \code{evgam} object
#'
#' @param object a fitted \code{evgam} object
#' @param newdata a data frame
#' @param type a character string giving the type of prediction sought; see Details. Defaults to \code{"link"}
#' @param prob a scalar or vector of probabilities for quantiles to be estimated if \code{type == "quantile"}; defaults to 0.5
#' @param se.fit a logical: should estimated standard errors be returned? Defaults to \code{FALSE}
#' @param marginal a logical: should uncertainty estimates integrate out smoothing parameter uncertainty? Defaults to \code{TRUE}
#' @param exi a logical: if a dependent GEV is fitted should the independent parameters be returned? Defaults to \code{FALSE}
#' @param as.gev a logical: should blended GEV parameters be converted to GEV parameters? Defaults to \code{FALSE}
#' @param margins a character string giving the QQ-plot margins. Defaults to original scale
#' @param trace an integer where higher values give more output. -1 suppresses everything. Defaults to 0
#' @param ... passed to \code{plot} if \code{type == 'qqplot'} or \code{type == 'qqplot2'}
#'
#' @details
#'
#' There are five options for \code{type}: 1) \code{"link"} distribution parameters 
#' transformed to their model fitting scale; 2) \code{"response"} as 1), but on their 
#' original scale; 3) "lpmatrix" a list of design matrices; 4) "quantile"
#' estimates of distribution quantile(s); and 5) "qqplot" a quantile-quantile
#' plot.
#'
#' @references 
#'
#' Youngman, B. D. (2022). evgam: An R Package for Generalized Additive Extreme
#' Value Modules. Journal of Statistical Software. To appear. \doi{10.18637/jss.v103.i03}
#'
#' @return A data frame or list of predictions, or a plot if \code{type == "qqplot"}
#'
#' @examples
#' 
#' data(fremantle)
#' fmla_gev <- list(SeaLevel ~ s(Year, k=5, bs="cr"), ~ 1, ~ 1)
#' m_gev <- evgam(fmla_gev, fremantle, family = "gev")
#' # prediction of link GEV parameter for fremantle data
#' predict(m_gev)
#' # predictions for Year 1989
#' y1989 <- data.frame(Year = 1989)
#' # link GEV parameter predictions
#' predict(m_gev, y1989)
#' # GEV parameter predictions
#' predict(m_gev, y1989, type= "response")
#' # 10-year return level predictions
#' predict(m_gev, y1989, type= "quantile", prob = .9)
#' # 10- and 100-year return level predictions
#' predict(m_gev, y1989, type= "quantile", prob = c(.9, .99))
#' 
#' @export
#'
predict.evgam <- function(object, newdata, type="link", prob=NULL, se.fit=FALSE, 
                          marginal=TRUE, exi = FALSE, as.gev = FALSE, margins = 'original',
                          trace = 0, ...) {
  
  ## a few checks
  
  family <- object$family
  
  if (type == 'qqplot2') {
    residual.qq <- TRUE
    type <- 'qqplot'
  } else {
    residual.qq <- FALSE
  }
  
  if (!(type %in% c('response', 'link', 'quantile', 'qqplot', 'lpmatrix')))
    stop(paste('type', type, 'not supported.'))
  
  if (family == "pp") {
    family <- "gev"
    if (missing(newdata)) {
      newdata <- object$data
      if (trace >= 0)
        message("Predictions for point process model given for quadrature points, not original data frame")
    }
  }
  
  if (!is.null(prob)) 
    type <- "quantile"
  if (type == "quantile" & is.null(prob)) 
    stop("non-NULL `prob' required if `type = quantile'")

  q_fn <- object$likfns$q
  dq_fn <- object$likfns$dq
  if (type == 'quantile') {
    if (is.null(q_fn))
      stop("custom.fns$q needs supplying for type = 'quantile' and family = 'custom'; see example.")
    if (se.fit)
      if (is.null(dq_fn))
        stop("custom.fns$dq needs supplying for type = 'quantile' if se.fit = TRUE and family = 'custom'; see example.")
  }
  p_fn <- object$likfns$p
  unlink_fns <- object$likfns$unlink
  if (type == 'response' & is.null(unlink_fns))
    stop("custom.fns$unlink needs supplying for type = 'response' and family = 'custom'.")

  if (type == "qqplot") {
    if (family != 'custom') {
      if (is.null(p_fn) | is.null(p_fn))
        stop(paste("`type = 'qqplot'' not yet available for family '", family, "'", sep=""))
    } else {
      stop("custom.fns$p and custom.fns$q needs supplying for type = 'qqplot' and family = 'custom'; see example.")
    }
  }
  ## end checks
  
  ## a few things to set aside
  
  if (se.fit) {
    
    if (marginal) {
      V.type <- "Vc" 
    } else {
      V.type <- "Vp"
    }
    
    conf.pars <- list(object$coefficients, object[[V.type]], object$idpars)
    
  }
  
  if (type == "qqplot")
    response.name <- object$response.name
  
  if (family == "exi")
    linkfn <- object$linkfn
  
  if (family == 'gpdab')
    gpdab <- matrix(object$gpdab, 2)
  
  got.newdata <- !missing(newdata)
  
  if (got.newdata) {
    pred.vars <- object$predictor.names
    missing.vars <- pred.vars[!(pred.vars %in% names(newdata))]
    if (length(missing.vars) > 0)
      stop(paste("Variable(s) '", paste(missing.vars, collapse=", "), "' not supplied to `newdata'.", sep=""))
  }
  
  if (got.newdata) {
    ndat <- nrow(newdata)
  } else {
    if (is.null(object$data)) 
      stop("Supply `evgam' with `removeData = FALSE' if not supplying `newdata'.")
    ndat <- nrow(object$data)
  }
  
  ## end stuff to set aside
  
  ## X creation starts
  
  if (!got.newdata)
    newdata <- NULL
  X <- .X.evgam(object, newdata, object$sparse)
  nX <- length(X)
  nms <- names(object)[seq_len(nX)]
  
  ## X creation ends
  
  if (type == "lpmatrix") {
    
    return(X)
    
  } else { ## start links
    
    out <- lapply(seq_len(nX), function(i) X[[i]] %*% object[[i]]$coefficients)
    names(out) <- names(X)
    out <- as.data.frame(lapply(out, function(x) x[,1]))
    out0 <- out
    
    if (se.fit) { ## start working with standard errors
      
      if (type != "quantile") {
        std.err <- lapply(seq_len(nX), function(i) sqrt(rowSums(X[[i]] * (X[[i]] %*% object[[i]][[V.type]]))))
        std.err <- as.data.frame(std.err)
        names(std.err) <- nms
      }
      
    } ## end of working with standard errors (for now)
    
    if (type %in% c("response", "quantile", "qqplot")) {
      
      if (type == "qqplot") 
        se.fit <- FALSE
      
      if (type == 'response') {
        
        for (i in seq_along(nms)) {
          
          lsti <- list(x = out0[[i]])
          if (!is.null(object$likdata$other)) {
            lsti <- c(lsti, object$likdata$other)
            if (!is.null(unlink_fns[[i]])) {
              fmls <- formalArgs(unlink_fns[[i]])
              lsti <- lsti[names(lsti) %in% fmls]
            }
          }

          if (se.fit) {
            if (!is.null(attr(unlink_fns[[i]], "deriv"))) {
              std.err[, i] <- do.call(attr(unlink_fns[[i]], "deriv"), lsti) * std.err[, i]
            }
          }
          
          if (!is.null(unlink_fns[[i]])) {
            out[, i] <- do.call(unlink_fns[[i]], lsti)
          }
          
        }
        
        if (family == 'bgev' & as.gev)
          out <- bgev2gev(out[, 1], out[, 2], out[, 3], object$likdata$other$pa, 
                          object$likdata$other$pb, object$likdata$other$alpha, 
                          object$likdata$other$beta, simplify = TRUE)
        
        names(out) <- object$rnms
        
      }
        
      if (se.fit & type == "response") 
        names(std.err) <- object$rnms
        
      if (type %in% c("quantile", "qqplot")) { ## convert response to quantile predictions
          
        pars <- out
        nprob <- length(prob)
        out <- matrix(NA, ndat, nprob)
        
        qnms <- paste('pars', 1:ncol(out0), sep = '')
        out0 <- as.data.frame(out0)
        out0 <- as.list(out0)
        names(out0) <- qnms
        
        if (!is.null(object$likdata$other)) {
          out0 <- c(out0, object$likdata$other)
        }
        
      }    
        
      if (type == 'quantile') {
        
        fmls <- formalArgs(q_fn)
        out0 <- out0[names(out0) %in% fmls]
        
        for (j in seq_len(nprob)) {
          
          pj <- prob[j]
          out0$p <- pj
          out[, j] <- do.call(q_fn, out0)
          
          if (se.fit) { ## standard errors for quantile predictions using Delta method
            
            Sigma <- array(NA, dim=c(ndat, nX, nX))
            idp <- conf.pars[[3]]
            for (i in seq_len(ndat)) {
              # rather inefficient: wants vectorising or Rcpp-ising
              for (j in seq_len(nX)) {
                for (k in j:nX) {
                  xj <- X[[j]][i,]
                  xk <- X[[k]][i,]
                  V <- conf.pars[[2]][idp == j, idp == k, drop=FALSE]
                  Sigma[i, j, k] <- sum(xj * (V %*% xk))
                  if (k != j) 
                    Sigma[i, k, j] <- Sigma[i, j, k]
                }
              }
            }
            
            std.err <- matrix(NA, ndat, nprob)
            
            for (j in seq_len(nprob)) {
              pj <- prob[j]
              out0$p <- pj
              
              jac <- do.call(dq_fn, out0)
              
              t1 <- apply(Sigma, 3, function(x) rowSums(x * jac))
              std.err[, j] <- rowSums(t1 * jac)
              
            }
            
            std.err <- as.data.frame(sqrt(std.err))
            names(std.err) <- paste("q", round(prob, 3), sep=":")
            
          } ## end std.err for quantile
        } ## end loop over probs for quantile 
        
        out <- as.data.frame(out)
        names(out) <- paste("q", round(prob, 3), sep=":")
        
      } ## end of response is quantile
      
      if (type == 'qqplot') {
        
        fmls <- formalArgs(p_fn)
        out0 <- out0[names(out0) %in% fmls]
        
        if (got.newdata) {
          y <- newdata[, response.name]
        } else {
          y <- object$likdata$y
          if (trace > 0)
            message("Response data taken from evgam object.")
        }
        x <- ppoints(sum(is.finite(y)))
        x <- replace(rep(NA, length(y)), is.finite(y), x)
        
        if (trace > 0)
          message("Margins converted to unit exponential by probability integral transformation.")
        
        if (margins %in% c('uniform', 'exponential')) {
          out0$x <- y
          y <- do.call(p_fn, out0)
          if (margins == 'exponential') {
            y <- qexp(y)
            x <- qexp(x)
          }
        } else {
          out0$p <- x
          x <- do.call(q_fn, out0)
        }
        
        x <- sort(x)
        y <- sort(y)
        
        if (residual.qq) {
          y <- y - x
          ylab <- "Empirical - Theoretical"
          ablist <- list(h = 0)
        } else {
          ylab <- "Empirical"
          ablist <- list(a = 0, b = 1)
        }
        
        if (residual.qq) {
          ylim <-  max(abs(y), na.rm = TRUE) * c(-1, 1)
        } else {
          ylim <- range(y, na.rm = TRUE)
        }
        plot(x, y, ylim = ylim, xlab = "Theoretical", ylab = ylab, 
             main = "Quantile-quantile plot", ...)
        do.call(abline, ablist)
        
      } ## end of qqplot
      
      if (se.fit)
        out <- list(fitted = out, se.fit = std.err)

    } 
      
    if (type != "qqplot") return(out)
    
  } ## end of link
  
}
