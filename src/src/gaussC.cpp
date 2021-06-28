// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// //' Gaussian negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each Gaussian distribution parameter
// //' @param X1 a design matrix for the GEV log scale parameter
// //' @param X2 a design matrix for the GEV shape parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return gaussd0 a scalar, the negative log-likelihood
// //' @return gaussd12 a matrix, first then second derivatives w.r.t. Gauss. dist. parameters
// //' @return gaussd34 a matrix, third then fourth derivatives w.r.t. Gauss. dist. parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double gaussd0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, arma::vec yvec, const arma::uvec& dupid, int dcate)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lsigvec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = yvec.size();

double nllh = 0.0;

double y, mu, lsig;

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lsigvec = lsigvec.elem(dupid);
}

for (int j=0; j < nobs; j++) {

mu = muvec[j];
lsig = lsigvec[j];
y = yvec[j];

y = y - mu;
y = y * y;

nllh += .5 * y / exp(2 * lsig) + lsig;
    
}

return(nllh);

}

// //' @rdname gaussd0
// [[Rcpp::export]]
arma::mat gaussd12(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, arma::vec yvec, const arma::uvec& dupid, int dcate)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lsigvec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 5);

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lsigvec = lsigvec.elem(dupid);
}

double y, mu, lsig;
double ee2, ee3, ee4, ee5; 

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lsig = lsigvec[j];

ee2 = exp(2 * lsig);
ee3 = y - mu;
ee4 = ee3/ee2;
ee5 = ee3 * ee3;

out(j, 0) = -ee4;
out(j, 1) = 1 - ee5/ee2;

out(j, 2)= 1/ee2;
out(j, 3) = 2 * ee4;
out(j, 4) = 2 * ee5/ee2;

}

return out;

}

// //' @rdname gaussd0
// [[Rcpp::export]]
arma::mat gaussd34(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, arma::vec yvec, const arma::uvec& dupid, int dcate)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lsigvec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 9);

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lsigvec = lsigvec.elem(dupid);
}

double y, mu, lsig;
double ee2, ee3, ee4; 

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lsig = lsigvec[j];

ee2 = exp(2 * lsig);
ee3 = y - mu;
ee4 = ee3 * ee3;

out(j, 0) = 0;
out(j, 1) = -(2/ee2);
out(j, 2) = -(4 * ee3/ee2);
out(j, 3) = -(4 * ee4/ee2);

out(j, 4) = 0;
out(j, 5) = 0;
out(j, 6) = 4/ee2;
out(j, 7) = 8 * ee3/ee2;
out(j, 8) = 8 * ee4/ee2;

}

return out;

}
