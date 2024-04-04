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
#' @param trace an integer where higher values give more output. -1 suppresses everything. Defaults to 0
#' @param ... unused
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
                          marginal=TRUE, exi = FALSE, as.gev = FALSE, trace = 0, ...) {
  
  ## a few checks
  
  family <- object$family
  
  if (!(type %in% c('response', 'link', 'quantile', 'qqplot', 'lpmatrix')))
      stop(paste('type', type, 'not supported.'))
  
  if (type == "qqplot" & !(family %in% c("gev", "gpd", "weibull")))
    stop(paste("`type = 'qqplot'' not yet available for family '", family, "'", sep=""))
  
  if (family == "pp") {
    family <- "gev"
    if (missing(newdata)) {
      newdata <- object$data
      if (trace >= 0)
        message("Predictions for point process model given for quadrature points, not original data frame")
    }
  }
  
  if (family == "egpd") {
    egpd_m <- object$likfns$m
    egpd_iG <- object$likfns$iG
  }
  
  if (family %in% c("custom", "bgev")) {
    q_fn <- object$likfns$q
    unlink_fns <- object$likfns$unlink
  }
  
  if (!is.null(prob)) 
    type <- "quantile"
  if (type == "quantile" & is.null(prob)) 
    stop("non-NULL `prob' required if `type = quantile'")
  if (type %in% "quantile" & family == "custom") {
    if (is.null(q_fn))
      stop("custom.fns$q needs supplying for type = 'quantile' and family = 'custom'.")
    if (is.null(unlink_fns))
      stop("custom.fns$unlink needs supplying for type = 'response' and family = 'custom'.")
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
  
  X <- .X.evgam(object, newdata)
  nX <- length(X)
  nms <- names(object)[seq_len(nX)]
  
  ## X creation ends
  
  if (type == "lpmatrix") {
    
    return(X)
    
  } else { ## start links
    
    out <- lapply(seq_len(nX), function(i) X[[i]] %*% object[[i]]$coefficients)
    names(out) <- names(X)
    out <- as.data.frame(lapply(out, function(x) x[,1]))
    
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
      
      if (family %in% c("custom", "bgev")) {
        
        for (i in seq_along(nms)) {
          
          if (se.fit) {
            if (!is.null(attr(unlink_fns[[i]], "deriv")))
              std.err[, i] <- attr(unlink_fns[[i]], "deriv")(out[, i]) * std.err[, i]
          }
          
          if (!is.null(unlink_fns[[i]]))
            out[, i] <- unlink_fns[[i]](out[, i])
          
        }
        
        if (family == 'bgev' & as.gev)
          out <- bgev2gev(out[, 1], out[, 2], out[, 3], object$likdata$other$pa, object$likdata$other$pb, object$likdata$other$alpha, object$likdata$other$beta)
        
      } else {
        
        if (family != "exi") {
          
          unlink <- which(substr(nms, 1, 3) == "log")
          for (i in unlink) {
            if (substr(nms[i], 1, 5) == "logit") {
              temp <- exp(-out[, i])
              out[, i] <- 1 / (1 + temp)
              if (se.fit & type == "response")
                std.err[, i] <- temp * out[, i] * std.err[, i] / (1 + temp)
              if (family == "gpdab")
                out[, i] <- gpdab[i, 1] + gpdab[i, 2] * out[, i]
            } else {
              out[, i] <- exp(out[, i])
            }
            if (se.fit & type == "response") {
              std.err[, i] <- out[, i] * std.err[, i]
              if (family == "gpdab")
                std.err[, i] <- gpdab[i, 2] * std.err[, i]
            }
          }
          
          if (exi & ncol(out) == 4) {
            out[, 4] <- out[, 4] ^ out[, 3]
            out[, 1] <- out[, 1] - out[, 2] * (1 - out[, 4]) / out[, 3]
            out[, 2] <- out[, 2] * out[, 4]
            out <- out[, 1:3, drop = FALSE]
            nms <- nms[1:3]
          }
          
        } else {
          
          if (se.fit) 
            std.err[, 1] <- attr(linkfn, "deriv")(out[, 1]) * std.err[, 1]
          out[, 1] <- linkfn(out[, 1])
          
        }
        
      }
      
      nms <- gsub("cloglog", "", nms)
      nms <- gsub("probit", "", nms)
      nms <- gsub("logit", "", nms)
      nms <- gsub("log", "", nms)
      nms <- gsub("trans", "", nms)
      names(out) <- nms
      
      if (type == "qqplot") { ## start qqplot
        
        pit <- !all(apply(out, 2, function(x) all(diff(x) < 1e-12)))
        x <- ppoints(ndat)
        y <- newdata[, response.name]
        if (is.null(y))
          stop("No response data.")
        
        if (!pit) {
          if (!(family %in% c("gev", "gpd", "gpdab", "weibull")))
            stop("Unsuitable `family' for `type == 'qqplot''")
          if (family == "gev")
            x <- .qgev(x, out[,1], out[,2], out[,3])
          if (family == "gpd")
            x <- .qgpd(x, 0, out[,1], out[,2])
          if (family == "weibull") 
            x <- .qweibull(x, out[,1], out[,2])
        } else {
          if (trace > 0)
            message("Margins converted to unit exponential by probability integral transformation.")
          x <- qexp(x)
          if (family == "gev")
            y <- .pgev(y, out[,1], out[,2], out[,3])
          if (family %in% c("gpd", "gpdab"))
            y <- .pgpd(y, 0, out[,1], out[,2], 0)
          if (family == "weibull") 
            y <- .pweibull(y, out[,1], out[,2])
          y <- qexp(y)
        }
        
        plot(x, sort(y), xlab = "Theoretical", ylab = "Empirical", main = "Quantile-quantile plot")
        abline(0, 1)
        
      } else { ## end qqplot
        
        if (se.fit & type == "response") 
          names(std.err) <- nms
        
        if (type == "quantile") { ## convert response to quantile predictions
          
          pars <- out
          nprob <- length(prob)
          out <- matrix(NA, ndat, nprob)
          
          for (j in seq_len(nprob)) {
            
            pj <- prob[j]
            
            if (family == "egpd") {
              if (egpd_m %in% c(1, 3)) {
                pj <- egpd_iG(pj, pars[, 3])
              } else {
                if (egpd_m == 2) {
                  pj <- egpd_iG(pj, pars[, 3], pars[, 4], pars[, 5])
                } else {
                  pj <- egpd_iG(pj, pars[, 3], pars[, 4])
                }
              }
            }
            
            if (family %in% c("custom", "bgev")) {
              
              if (family == 'bgev') {
                out[, j] <- q_fn(pj, pars[,1], pars[,2], pars[,3], 
                                 object$likdata$other[1], object$likdata$other[2], object$likdata$other[3], object$likdata$other[4])
              } else {
                
                if (length(nms) > 4)
                  stop ("Currently on predictions with non-NULL prob and family = 'custom' only possible for fewer than five parameters.")
                
                if (length(nms) == 1) {
                  out[, j] <- q_fn(pj, pars[,1])
                } else {
                  if (length(nms) == 2) {
                    out[, j] <- q_fn(pj, pars[,1], pars[,2])
                  } else {
                    if (length(nms) == 3) {
                      out[, j] <- q_fn(pj, pars[,1], pars[,2], pars[,3])
                    } else {
                      out[, j] <- q_fn(pj, pars[,1], pars[,2], pars[,3], pars[,4])
                    }
                  }
                }
              }
              
            } else {
              
              if (family %in% c("gpd", "egpd", "gpdab")) {
                out[, j] <- .qgpd(pj, 0, pars[,1], pars[,2])
              } else {
                if (family == "gev") {
                  out[, j] <- .qgev(pj, pars[,1], pars[,2], pars[,3])
                } else {
                  if (family == "weibull") {
                    out[, j] <- .qweibull(pj, scale=pars[,1], shape=pars[,2])
                  } else {
                    stop("invalid family")
                  } 
                }
              }
            }
            
          }
          
          if (se.fit) { ## standard errors for quantile predictions using Delta method
            
            if (family %in% c("egpd", "custom", "bgev")) 
              stop("Standard errors not yet available for this family.")
            
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
              if (family == "gev") {
                jac <- .dqgev(prob[j], pars[,1], log(pars[,2]), pars[,3])
              }
              if (family == "gpd") {
                jac <- .dqgpd(prob[j], log(pars[,2]), pars[,3])
              }
              if (family == "weibull") {
                jac <- .dqweibull(prob[j], log(pars[,2]), pars[,3])
              }
              for (i in seq_len(ndat)) {
                std.err[i, j] <- sum(jac[i,] * (Sigma[i,,] %*% jac[i,]))
              }
            }
            
            std.err <- as.data.frame(sqrt(std.err))
            names(std.err) <- paste("q", round(prob, 3), sep=":")
            
          } ## end std.err for quantile
          
          out <- as.data.frame(out)
          names(out) <- paste("q", round(prob, 3), sep=":")
          
        } ## end of response to quantile
        
        if (se.fit) {
          
          out <- list(fitted = out, se.fit = std.err)
          
        }
        
      } ## end not qqplot
      
    } ## end of not link
    
    if (type != "qqplot") return(out)
    
  } ## end of link
  
}
