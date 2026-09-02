## evgam 1.0.2

### Changes:

* Smooths can now use `id` argument of `mgcv::s()` to allow multiple smooths to have same smoothing parameter.

* Objects in workspace can now be used in `evgam()`'s `formula` argument (as opposed to all variables needing to be part of `data`).

* Added function and class `geoevgam` with S3 method `plot` for fitting and plotting geometric extreme value models. (Thanks to Callum Murphy-Barltrop, Jenny Wadsworth and Miguel de Carvalho.)

* Added `type = 'qqplot2'` to `predict()` for residual-based QQ-plots. (Thanks to Jenny Wadsworth.)

* Added `family = 'gw'` for fitting generalised Weibull distribution. (Thanks to Cees de Valk.)

* Added `family = 'poisson'` for fitting Poisson distribution.

* Added `args(dist = 'aggauss')` for fitting `family = 'condex'` with asymmetric generalised Gaussian residuals. (Thanks again to Kristina Bratkova and Aiden Farrell.)

* Added `family = 'aggauss'` for fitting asymmetric generalised Gaussian distribution. (Thanks to Kristina Bratkova and Aiden Farrell.)

* GPD model with shape parameter constrained to [-0.5, 1.0] added with `family = "gpd2"`.

* Added `ltgamma` and `ltgammab` families for the left truncated gamma distribution with unknown and know shape, respectively. Use `args = list(lower = )` to give the scalar, vector or matrix of left-truncation points, and `args = list(alpha = )` to give scalar, vector or matrix of gamma distribution shape parameters with family `ltgammab`.

* Added `condex` family, for conditional extreme value model of Heffernan and Tawn (2004, https://doi.org/10.1111/j.1467-9868.2004.02050.x). For an example see https://byoungman.github.io/evgam/condex/.

* Added `beta` and `logitgauss` families.

* Added option `sparse = TRUE`, which coerces matrices to sparse matrices through package `Matrix` where possible. Particularly useful for `bs = 're'` or `bs = 'mrf'`.

* Added `df2matdf()` for turning a vector response to a matrix response if explanatory variable combinations are repeated.

* Added functionality to fit extended generalised Pareto distribution through `evgam(..., family = "egpd")`. See Naveau et al. (W2016, https://doi.org/10.1002/2015WR018552) and `family.evgam`. (Thanks to Touqeer Ahmad.)

* Added functionality to fit blended generalised extreme value (GEV) distribution through `evgam(..., family = "bgev")`. See Castro-Camilo et al. (2022, https://doi.org/10.1002/env.2742) and `family.evgam`. (Thanks to Jordan Richards.)

* Also added `dbgev()`, `pbgev()`, `qbgev()` and `rbgev()` for density, distribution function, quantile function and random generation, respectively, for the blended GEV distribution.

* Added functionality to fit models via custom likelihood functions, i.e. extending those available in `evgam` through `family = ...`. See `custom.family.evgam`.

* Added functionality to constrain both GPD parameters using `gpd.args = list(lower = ..., upper = ...)`. (Thanks to Callum Murphy-Barltrop.)

* GEV model with shape parameter constrained to [-0.5, 1.0] added with `family = "gev2"`.

### Bug fixes:

* Slight change to calculation of log determinant of penalty matrices and their derivatives for tensor smooths. Should implement algorithm in Appendix B of Wood (2011, https://doi.org/10.1111/j.1467-9868.2010.00749.x)

* That only variables are checked as being supplied to `data` is now properly detected. (Thanks, Simon Brown.)

* That values of smoothing parameters supplied to `evgam()` are properly recognised has been fixed.

## evgam 1.0.1

### Changes:

* Makevars and Makevars.win updated.

## evgam 1.0.0

### Changes:

* Version increased to 1.0.0 to reflect publication of Youngman (2022, JSS, \doi{10.18637/jss.v103.i03}).

* References to Youngman (2022, JSS, \doi{10.18637/jss.v103.i03}) added, where appropriate.

### Bug fixes:

* That all variables have been supplied to `data` is now properly detected.

## evgam 0.1.4

### Changes:

* None.

### Bug fixes:

* An error is thrown if there are fewer than r data for any pp.args$id, as opposed to r + 1 incorrectly implemented previously. (Thanks, Yousra El Bachir.)

## evgam 0.1.3

### Changes:

* `plot()` for an evgam object now calls `mgcv::plot.gam()` to plot smooths (with thanks to Debbie Dupuis for triggering this). `plot()` no longer has the `addMap` option, for adding map outlines via `maps::map()`; instead using one-figure devices with `maps::map()` separately is recommended.

* Calculations of log(|S|_+) for penalty matrix S now fully implements Wood (JRSSB, 2011(73)1, Appendix B).

* Calculations of log(|H|) for Hessian H now use diagonality simpifications; see Wood (book: GAMs in R 2nd ed. (2017) pp. 286).

## evgam 0.1.2

### Changes:

* The Fremantle data from package ismev have been added, and are used for examples. Usage is `data(fremantle)`, as in ismev.

* `colplot()` adds the option to add a legend, which defaults to `FALSE`.

* `logLik.evgam()` now returns an object of class `'logLik'`, allowing, e.g., `AIC()` and `BIC()` to be used.

* `extremal0()` has gone, as `extremal()` can now do the same.

* `evgam()`'s trace argument now allows -1, which suppresses any information on the console.

### Bug fixes:

* Negative response data now work okay with `family = "ald"`.

* `evgams()`'s formula argument may have smooths and parametric-only terms in any order. (Previously, smooths had
    to come first, so `formula = list(response ~ s(), ~ 1, ~ s())` broke.)
    
* `predict.evgam(object)` with `missing(newdata)` only gave one set predictions for `object$data`. It now gives predictions
   for all rows of `object$data` (as it should).

## evgam 0.1.1

### Changes:

* `plot.evgam()` now has informative y-axis labels for one-dimensional smooths.

### Bug fixes:

* Compilation flag with clang++ in gradHess.cpp addressed.

* `simulate.evgam()` correctly labels variables for `family = "response"`.

## evgam 0.1.0

* Initial release.
