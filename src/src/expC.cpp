// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// //' Exponential distribution negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each Exponential parameter
// //' @param X1 a design matrix for the Exponential log rate parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return expd0 a scalar, the negative log-liklihood
// //' @return expd12 a matrix, first then second derivatives w.r.t. Exponential parameters
// //' @return expd34 a matrix, third then fourth derivatives w.r.t. Exponential parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double expd0(const Rcpp::List& pars, const arma::mat& X1, arma::vec yvec, const arma::uvec dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
int nobs = yvec.size();

if (dcate == 1) {
    lpsivec = lpsivec.elem(dupid);
}

double y, lpsi;
double nllh=0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];
nllh += y * exp(lpsi) - lpsi;
    
}

return(nllh);

}

// //' @rdname expd0
// [[Rcpp::export]]
arma::mat expd12(const Rcpp::List& pars, arma::mat X1, arma::vec yvec, const arma::uvec dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 2);

if (dcate == 1) {
    lpsivec = lpsivec.elem(dupid);
}

double y, lpsi, ee2;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];

ee2 = y * exp(lpsi);

out(j, 0) = ee2 - 1;
out(j, 1) = ee2;

}

return out;

}

// //' @rdname expd0
// [[Rcpp::export]]
arma::mat expd34(const Rcpp::List& pars, arma::mat X1, arma::vec yvec, const arma::uvec dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 2);

if (dcate == 1) {
    lpsivec = lpsivec.elem(dupid);
}

double y, lpsi, ee2;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];

ee2 = y * exp(lpsi);

out(j, 0) = ee2;
out(j, 1) = ee2;

}

return out;

}
