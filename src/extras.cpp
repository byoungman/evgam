// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export]]
arma::mat armapinv(arma::mat x, double tol)
{
if (tol < 0) {
  return arma::pinv(x);
} else {
  return arma::pinv(x, tol);
}
}

// [[Rcpp::export]]
arma::mat armaginv(arma::mat x, double tol)
{
return arma::pinv(x, tol, "std");
}
