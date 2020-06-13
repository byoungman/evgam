// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// //' Weibull distribution negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each Weibull parameter
// //' @param X1 a design matrix for the Weibull log scale parameter
// //' @param X2 a design matrix for the Weibull log shape parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return weibd0 a scalar, the negative log-liklihood
// //' @return weibd12 a matrix, first then second derivatives w.r.t. Weibull parameters
// //' @return weibd34 a matrix, third then fourth derivatives w.r.t. Weibull parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double weibd0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::vec yvec, arma::uvec dupid, int dcate)
{
    
arma::vec llambdavec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lkvec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = yvec.size();

if (dcate == 1) {
    llambdavec = llambdavec.elem(dupid);
    lkvec = lkvec.elem(dupid);
}

double y, ll, lk;
double ee1;
double nllh=0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
ll = llambdavec[j];
lk = lkvec[j];

ee1 = exp(lk);

nllh -= ((ee1 - 1) * (log(y) - ll) + lk - (R_pow(y/exp(ll), ee1) + ll));

}

return(nllh);

}

// //' @rdname weibd0
// [[Rcpp::export]]
arma::mat weibd12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::vec yvec, arma::uvec dupid, int dcate)
{
    
arma::vec llambdavec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lkvec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 5);

if (dcate == 1) {
    llambdavec = llambdavec.elem(dupid);
    lkvec = lkvec.elem(dupid);
}

double y, ll, lk;
double ee1, ee2, ee3, ee4, ee5, ee7, ee8, ee9;

for (int j=0; j < nobs; j++) {

y = yvec[j];
ll = llambdavec[j];
lk = lkvec[j];

ee1 = exp(lk);
ee2 = exp(ll);
ee3 = y/ee2;
ee4 = ee1 - 1;
ee5 = R_pow(ee3, ee4);
ee7 = log(y) - ll;
ee8 = R_pow(ee3, ee1);
ee9 = ee1 * ee7;

out(j, 0) = (1 - y * ee5/ee2) * ee1;
out(j, 1) = -((1 - ee8) * ee1 * ee7 + 1);
out(j, 2) = ee4 * (R_pow(y, 2)/(R_pow(ee2, 2) * R_pow(ee3, 2)) - 1) + y * (ee5 + y * ee4 * R_pow(ee3, ee1 - 2)/ee2) * ee1/ee2;
out(j, 3) = (1 - y * (ee5 + ee9 * ee5)/ee2) * ee1;
out(j, 4) = -((1 - (ee8 + ee9 * ee8)) * ee1 * ee7);

}

return out;

}

// //' @rdname weibd0
// [[Rcpp::export]]
arma::mat weibd34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::vec yvec, arma::uvec dupid, int dcate)
{
    
arma::vec llambdavec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lkvec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 9);

if (dcate == 1) {
    llambdavec = llambdavec.elem(dupid);
    lkvec = lkvec.elem(dupid);
}

double y, ll, lk;
double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee17, ee18, ee19;
double ee20, ee21, ee22, ee23, ee26, ee27, ee28, ee29;
double ee30, ee32, ee33;

for (int j=0; j < nobs; j++) {

y = yvec[j];
ll = llambdavec[j];
lk = lkvec[j];

ee1 = exp(lk);
ee2 = exp(ll);
ee3 = y/ee2;
ee4 = ee1 - 1;
ee5 = ee1 - 2;
ee6 = R_pow(ee3, ee4);
ee7 = R_pow(ee3, ee5);
ee9 = log(y) - ll;
ee10 = R_pow(ee3, ee1);
ee11 = R_pow(ee3, 2);
ee12 = ee1 - 3;
ee13 = R_pow(ee3, ee12);
ee14 = ee1 * ee9;
ee15 = 2 * ee7;
ee17 = R_pow(ee2, 2) * ee11;
ee18 = R_pow(y, 2);
ee19 = ee7 + ee15;
ee20 = 2 * ee6;
ee21 = ee14 * ee6;
ee22 = 2 * ee10;
ee23 = ee18/ee17;
ee26 = (ee6 + ee20 + ee21) * ee1 * ee9 + ee6;
ee27 = ee4 * ee7;
ee28 = ee6 + ee21;
ee29 = ee19 + y * ee5 * ee13/ee2;
ee30 = 1/ee11;
ee32 = 2 * ee23 - 3;
ee33 = ee14 * ee10;

out(j, 0) = (1 + ee18 * ee32/ee17) * ee4 - y * (ee6 + y * ee29 * ee4/ee2) * ee1/ee2;
out(j, 1) = ee1 * (y * (ee28 + y * ((ee4 * ee9 * ee7 + ee7) * ee1 +
   ee27 + ee30)/ee2)/ee2 - 1);
out(j, 2) = (1 - y * ee26/ee2) * ee1;
out(j, 3) = -((1 - ((ee10 + ee22 + ee33) * ee1 * ee9 + ee10)) * ee1 * ee9);
out(j, 4) = ee4 * (ee18 * (7 + ee18 * (8 * ee23 - 14)/ee17)/ee17 -
   1) + y * (ee6 + y * (ee19 + 4 * ee7 + y * (ee13 + ee13 +
   4 * ee13 + y * ee12 * R_pow(ee3, ee1 - 4)/ee2) * ee5/ee2) * ee4/ee2) * ee1/ee2;
out(j, 5) = (1 + y * (y * (ee32/ee11 - ((ee19 * ee4 * ee9 +
   ee7 + ee15 + y * ((ee5 * ee9 * ee13 + ee13) * ee4 + ee5 * ee13)/ee2) * ee1 +
   ee29 * ee4))/ee2 - ee28)/ee2) * ee1;
out(j, 6) = ee1 * (y * (ee26 + y * ((((ee19 + ee14 * ee7) * ee4 +
   2 * (ee1 * ee7)) * ee9 + ee7 + ee7 + ee7) * ee1 + ee27 +
   ee30)/ee2)/ee2 - 1);
out(j, 7) = (1 - y * (((ee20 + 4 * ee6 + ee21) * ee1 * ee9 +
   ee6 + ee6 + ee6 + ee20 + ee20) * ee1 * ee9 + ee6)/ee2) * ee1;
out(j, 8) = -((1 - (((ee22 + 4 * ee10 + ee33) * ee1 * ee9 +
   ee10 + ee10 + ee10 + ee22 + ee22) * ee1 * ee9 + ee10)) * ee1 * ee9);

}

return out;

}
