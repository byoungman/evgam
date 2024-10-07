#' Distribution families in \pkg{evgam}
#'
#' Various distributions can be fitted with function \code{evgam} using 
#' \code{family = "..."}, where options for \code{...} are given below. The default is 
#' \code{family = "gev"}.
#'
#' @details
#' 
#' The following families are currently available using \code{evgam(..., family = "...")}.
#'
#' \itemize{
#'   \item \code{"ald"}, the asymmetric Laplace distribution. This is primarily 
#' intended for quantile regression, as in Yu & Moyeed (2001);
#'   \item \code{"gev"} (default), the generalised extreme value (GEV) distribution; 
#'   \item \code{"exp"}, the exponential distribution;
#'   \item \code{"gpd"}, the generalised Pareto distribution; 
#'   \item \code{"gauss"}, the Gaussian distribution.
#'   \item \code{"pp"}, the point process model for extremes. This is implemented 
#' through \eqn{r}-largest order statistics. See Details below. 
#'   \item \code{"weibull"}, the Weibull distribution;
#'   \item \code{"exi"}, estimation if the extremal index. See Schlather & Tawn (2003) and Details below.
#'   \item \code{"egpd"}, the extended generalised Pareto distribution. See Naveau
#' et al. (2016) and Details below;
#'   \item \code{"bgev"}, the blended GEV distribution. See Castro-Camilo et al (2022) and Details
#'   \item \code{"condex"}, the conditional extreme value model. See Details
#'   \item \code{"custom"}, custom distributions. See \code{custom.evgam} for an example of use. 
#' }
#'
#' Arguments for the asymmetric Laplace distribution are given by \code{ald.args}. A 
#' scalar \code{tau} defines the quantile sought, which has no default. The scalar
#' \code{C} specifies the curvature parameter of Oh et al. (2011).
#'
#' Arguments for extremal index estimation are given by \code{exi.args}. A character
#' string \code{id} specifies the variable in \code{data}over which an \code{nexi}
#' (default 2) running max. has been taken. The \code{link} is specified as a character string,
#' which is one of \code{"logistic"}, \code{"probit"}, \code{"cloglog"}; defaults to \code{"logistic"}.
#' 
#' Arguments for the point process model are given by \code{pp.args}. An integer \code{r}
#' specifies the number of order statistics from which the model will be estimated.
#' If \code{r = -1}, all \code{data} will be used. The character string \code{id} specifies the variable 
#' in \code{data} over which the point process isn't integrated; e.g. if a map 
#' of parameter estimates related to extremes over time is sought, integration isn't 
#' over locations. The scalar \code{nper} number of data per period of interest; scalar or
#' integer vector \code{ny} specifies the number of periods; if \code{length(ny) > 1}
#' then \code{names(ny)} must ne supplied and must match to every unique \code{id}. 
#' logical \code{correctny} specifies whether \code{ny} is 
#' corrected to adjust proportionally for data missingness.
#'
#' Arguments for the point process model are given by \code{bgev.args}. Probabilities 
#' \code{pa} and \code{pb} specify the lower and upper probabilities at the which
#' the Gumbel distribution blends into a GEV distribution. Then \code{alpha} and
#' \code{beta} specify the quantile and its range, respectively, used to parameterise
#' the GEV distribution. Defaults are \code{pa = 0.05} and \code{pb = 0.2} and
#' \code{alpha = beta = 0.5}, as used in Castro-Camilo et al (2022). 
#' 
#' Arguments for extended Pareto distribution are given by \code{egpd.args}. An 
#' integer, \code{model}, specifies which model from Naveau et at. (2016) to fit. 
#' The first two parameters of each are the GPD's log scale and shape parameters, 
#' \eqn{(\log \psi, \xi)}. Then, in the notation of Naveau et at. (2016) the 
#' remaining parameters are \eqn{(\log \kappa)}, 
#' \eqn{(\log \kappa_1, \log \kappa_2, \textrm{logit}(p))}, 
#' \eqn{(\log \delta)} and \eqn{(\log \delta, \log \kappa)} for models i, ii, iii
#' and iv, respectively, which are specified with \code{model = 1, 2, 3} or 
#' \code{4}, respectively.
#'
#' See \code{\link{evgam}} for examples.
#'
#' @seealso \link{evgam}
#'
#' @references 
#' 
#' Castro-Camilo, D., Huser, R., & Rue, H. (2022). Practical strategies for 
#' generalized extreme value-based regression models for extremes. 
#' Environmetrics, 33(6), e2742. \doi{10.1002/env.2742}
#' 
#' Naveau, P., Huser, R., Ribereau, P., and Hannart, A. (2016), Modeling 
#' jointly low, moderate, and heavy rainfall intensities without a threshold 
#' selection, Water Resources Research, 52, 2753-2769. \doi{10.1002/2015WR018552}
#'
#' Oh, H. S., Lee, T. C., & Nychka, D. W. (2011). Fast nonparametric 
#' quantile regression with arbitrary smoothing methods. Journal of 
#' Computational and Graphical Statistics, 20(2), 510-526. 
#' \doi{10.1198/jcgs.2010.10063}
#'
#' Schlather, M., & Tawn, J. A. (2003). A dependence measure for multivariate and 
#' spatial extreme values: Properties and inference. Biometrika, 90(1), 139-156.
#' \doi{10.1093/biomet/90.1.139}
#'
#' Youngman, B. D. (2022). evgam: An R Package for Generalized Additive Extreme
#' Value Models. Journal of Statistical Software. \doi{10.18637/jss.v103.i03}
#'
#' Yu, K., & Moyeed, R. A. (2001). Bayesian quantile regression. 
#' Statistics & Probability Letters, 54(4), 437-447.\doi{10.1016/S0167-7152(01)00124-9}
#' 
#' @name family.evgam
#' 
NULL

#' Custom distributions with \pkg{evgam}
#'
#' Users can supply code to fit an arbitrary distribution with \code{evgam} using
#' \code{family = "custom"}. See Details for more information and the example below,
#' which demonstrates fitting of the Gumbel distribution.
#'
#' @details
#' 
#' Users should supply a list to \code{evgam} called \code{custom.fns}, which  
#' which comprises the following functions: \code{d0}, 
#' which evaluates the negative log-likelihood; code{d120}, which evaluates its 
#' first and second derivatives; and, optionally, \code{d340}, which evaluates its 
#' third and fourth derivatives. The list may also contain \code{q}, which evaluates 
#' the custom distribution's quantile function, for use with \code{predict, prob
#' = ...)}. The list may also contain \code{unlink}, which gives functions to 
#' reverse any link functions (unlink) used with linear predictors, so that 
#' \code{predict(..., type = "response")} works meaningfully, and \code{unlink}
#' may also have an attribute \code{deriv}, which gives derivatives of unlink
#' functions, which allows \code{predict(..., type = "response", std.err = TRUE)} 
#' to return meaningful values.
#'
#' Trying to mimic the Gumbel example below for your chosen distribution is almost
#' certainly the easiest was to get \code{family = "custom"} to do what you want.
#' 
#' @examples
#'
#' # Gumbel custom likelihood
#' 
#' # negative log likelihood
#' gumd0 <- function(pars, likdata) {
#' 
#' # this is needed, and should be fine
#' pars <- split(pars, likdata$idpars)
#' 
#' # likdata$X should be set up by evgam
#' mu <- likdata$X[[1]] %*% pars[[1]]
#' lsigma <- likdata$X[[2]] %*% pars[[2]]
#' y <- likdata$y
#' 
#' y <- (y - mu) * exp(-lsigma)
#' nllh <- sum(lsigma - y + exp(y))
#' 
#' return(nllh)
#' 
#' }
#' 
#' # first and second derivatives of neg. log lik. in order
#' # d_mu, d_lsigma, d_{mu, mu}, d_{mu, lsigma}, d_{lsigma, lsigma}
#' gumd12 <- function(pars, likdata) {
#' 
#' # this is needed, and should be fine
#' pars <- split(pars, likdata$idpars)
#' 
#' # likdata$X should be set up by evgam
#' mu <- likdata$X[[1]] %*% pars[[1]]
#' lsigma <- likdata$X[[2]] %*% pars[[2]]
#' y <- likdata$y
#' 
#' out <- matrix(0, length(y), 5)
#' 
#' ee1 <- exp(lsigma)
#' ee2 <- y - mu
#' ee3 <- ee2/ee1
#' ee4 <- exp(ee3)
#' ee5 <- (ee3 + 1) * ee4
#' ee6 <- 1 - ee4
#' 
#' # first derivatives
#' out[, 1] <- ee6/ee1
#' out[, 2] <- ee6 * ee2/ee1 + 1
#' 
#' # second derivatives
#' out[, 3] <- ee4/ee1^2
#' out[, 4] <- -((1 - ee5)/ee1)
#' out[, 5] <- (ee5 - 1) * ee2/ee1
#' 
#' return(out)
#' 
#' }
#' 
#' gum_fns <- list(d0 = gumd0, d120 = gumd12, d340 = NULL)
#' 
#' unlink <- list(NULL, function(x) exp(x))
#' attr(unlink[[2]], "deriv") <- unlink[[2]]
#' qgumbel <- function(p, location, scale) location - scale * log(-log(p))
#' 
#' gum_fns$q <- qgumbel
#' gum_fns$unlink <- unlink
#'
#'
#' \donttest{
#' 
#' data(COprcp)
#' COprcp <- cbind(COprcp, COprcp_meta[COprcp$meta_row,])
#' COprcp$year <- format(COprcp$date, "%Y")
#' COprcp_gev <- aggregate(prcp ~ year + meta_row, COprcp, max)
#' COprcp_gev <- cbind(COprcp_gev, COprcp_meta[COprcp_gev$meta_row,])
#' 
#' inits <- sqrt(6 * var(COprcp_gev$prcp) / pi^2)
#' inits <- c(mean(COprcp_gev$prcp) - inits[1] * 0.5772156649, log(inits[1]))
#' fmla_gum <- list(location = prcp ~ s(lon, lat) + s(elev), logscale = ~ s(lon, lat))
#' m <- evgam(fmla_gum, data = COprcp_gev, family = 'custom', custom.fns = gum_fns, 
#' trace = 2, inits = inits)
#' 
#' predict(m, COprcp_gev[1:10,])
#' predict(m, COprcp_gev[1:10,], type = 'response')
#' predict(m, COprcp_gev[1:10,], prob = .99)
#'
#' }
#' 
#' @seealso \link{family.evgam}
#'
#' @references 
#' 
#' Youngman, B. D. (2022). evgam: An R Package for Generalized Additive Extreme
#' Value Modules. Journal of Statistical Software. To appear. \doi{10.18637/jss.v103.i03}
#' 
#' @name custom.family.evgam
#' 
NULL
