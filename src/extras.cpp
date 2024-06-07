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
// 
// // [[Rcpp::export]]
// arma::mat armaspsolve(arma::sp_mat a, arma::mat b)
// {
//   return arma::spsolve(a, b, "lapack");
// }
// 
// // [[Rcpp::export]]
// arma::mat armaspLsolve(arma::sp_mat L, arma::mat b)
// {
//   b = arma::spsolve(L.t(), b, "lapack");
//   return arma::spsolve(L, b, "lapack");
// }
