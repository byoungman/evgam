## evgam 0.1.5

### Changes:

* 'gamma' argument

* custom likelihoods

* extended GPD

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
