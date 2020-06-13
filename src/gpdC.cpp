// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double xieps = 0.0;
const double xieps3 = 0.0;

// //' Generalized Pareto distribution (GPD) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each GPD parameter
// //' @param X1 a design matrix for the GEV log scale parameter
// //' @param X2 a design matrix for the GEV shape parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return gpdd0 a scalar, the negative log-likelihood
// //' @return gpdd12 a matrix, first then second derivatives w.r.t. GPD parameters
// //' @return gpdd34 a matrix, third then fourth derivatives w.r.t. GPD parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double gpdd0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, arma::vec yvec, const arma::uvec& dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xivec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = yvec.size();

if (dcate == 1) {
    lpsivec = lpsivec.elem(dupid);
    xivec = xivec.elem(dupid);
}

double y, lpsi, xi;
double ee1;
double nllh=0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];

if (fabs(xi) >= xieps3) {

ee1 = xi * y / exp(lpsi);

if (ee1 <= -1.0) {
    nllh = 1e20;
    break;
} else {

nllh += lpsi + (1.0 / xi + 1.0) * log1p(ee1);
}

} else {

nllh += lpsi + y / exp(lpsi);
    
}

}

return(nllh);

}

// //' @rdname gpdd0
// [[Rcpp::export]]
arma::mat gpdd12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::vec yvec, const arma::uvec dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xivec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 5);

if (dcate == 1) {
    lpsivec = lpsivec.elem(dupid);
    xivec = xivec.elem(dupid);
}

double y, lpsi, xi;
double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9, ee10, ee12; 

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];

if (fabs(xi) >= xieps) {

ee1 = exp(lpsi);
ee2 = xi * y;
ee3 = ee2/ee1;
ee4 = (1 + ee3) * ee1;
ee5 = 1/xi;
ee6 = 1 + ee5;
ee7 = xi * xi;
ee8 = log1p(ee3);
ee9 = ee2 * ee6;
ee10 = ee2/ee4;
ee12 = y * ee6/ee4;

out(j, 0) = 1 - ee9/ee4;
out(j, 1) = ee12 - ee8/ee7;

out(j, 2)= -(ee9 * (ee10 - 1)/ee4);
out(j, 3) = -(y * ((1 - ee10) * ee6 - ee5)/ee4);
out(j, 4) = -((y/ee4 - 2 * (ee8/xi))/ee7 + y * (1/ee7 + ee12)/ee4);

} else {
    
ee1 = exp(lpsi);

out(j, 0) = 1 - y/ee1;
out(j, 0) = 0;

out(j, 0) = 1 * y/ee1;
out(j, 0) = 0;
out(j, 0) = 0;

}

}

return out;

}

// //' @rdname gpdd0
// [[Rcpp::export]]
arma::mat gpdd34(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::vec yvec, const arma::uvec dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xivec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 9);

if (dcate == 1) {
    lpsivec = lpsivec.elem(dupid);
    xivec = xivec.elem(dupid);
}

double y, lpsi, xi;
double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9; 
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee19; 
double ee20, ee22, ee24, ee25, ee27, ee28; 
double ee30, ee31, ee33, ee34, ee35; 

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];

if (fabs(xi) >= xieps) {

ee1 = exp(lpsi);
ee2 = xi * y;
ee3 = ee2/ee1;
ee4 = 1 + ee3;
ee5 = ee4 * ee1;
ee6 = ee2/ee5;
ee7 = xi * xi;
ee8 = y/ee5;
ee9 = 1 + 1/xi;
ee10 = 2 * ee6;
ee11 = 2 * ee8;
ee12 = 2/xi;
ee13 = 2 - ee10;
ee14 = ee11 + ee12;
ee15 = ee11 + 6/xi;
ee16 = 4 * ee4;
ee17 = y * ee9;
ee19 = ee15/xi + y * ee14/ee5;
ee20 = 1 - ee6;
ee22 = 1 + ee2 * (ee10 - 3)/ee5;
ee24 = 2 * (1 + 2 * ee3) + ee16;
ee25 = ee12 + ee8;
ee27 = 6 * ee3;
ee28 = 8 * ee3;
ee30 = log1p(ee3)/xi;
ee31 = ee2 * ee9;
ee33 = ee2 * (3 - ee10)/ee5;
ee34 = ee6 - 1;
ee35 = ee17/ee5;

// third derivatives
// 1=log(scale), 2=shape
// order: 111, 112, 122, 222

out(j, 0) = -(ee31 * ee22/ee5);
out(j, 1) = -(y * (ee9 * (ee33 - 1) - ee34/xi)/ee5);
out(j, 2) = y * (y * (ee9 * ee13 - ee12)/ee5)/ee5;
out(j, 3) = ((ee11 - 6 * 
    ee30)/xi + y * ee25/ee5)/ee7 + y * (ee25/ee7 + y * (1/ee7 + 
    2 * ee35)/ee5)/ee5;

// fourth derivatives
// 1=log(scale), 2=shape
// order: 1111, 1112, 1122, 1222, 2222

out(j, 4) = -(ee31 * (ee2 * 
    (7 + ee2 * ((ee28 - ee24)/ee4 - 6)/ee5)/ee5 - 1)/ee5);
out(j, 5) = -(y * (ee9 * (1 + ee2 * (ee2 * ((ee24 - 
ee28)/ee4 + 6)/ee5 - 7)/ee5) - ee22/xi)/ee5);
out(j, 6) = -(y * 
    (ee17 * (4 - ee2 * ((ee16 - ee27)/ee4 + 6)/ee5)/ee5 - 
    (2 * ee33 - (2 + 2 * ee34))/ee7)/ee5);
out(j, 7) = -(y * 
    (((4 * ee20 - 6)/xi + 2 * (ee20/xi) + 2 * y * ee13/ee5)/ee7 + 
    y * (ee13/ee7 + y * ((2 * ee4 - ee27)/ee4 + 4) * 
    ee9/ee5)/ee5)/ee5);
out(j, 8) = (((24 * ee30 - 
    6 * ee8)/xi - y * ee15/ee5)/xi - y * ee19/ee5)/ee7 - 
    y * (ee19/ee7 + y * (ee14/ee7 + y * (2/ee7 + 6 * 
    ee35)/ee5)/ee5)/ee5;

} else {

ee2 = 1 * y/exp(lpsi);

// third derivatives
// 1=log(scale), 2=shape
// order: 111, 112, 122, 222

out(j, 0) = -ee2;
out(j, 1) = 0;
out(j, 2) = 0;
out(j, 3) = 0; 

// fourth derivatives
// 1=log(scale), 2=shape
// order: 1111, 1112, 1122, 1222, 2222

out(j, 4) = ee2;
out(j, 5) = 0;
out(j, 6) = 0;
out(j, 7) = 0;
out(j, 8) = 0;
        
}

}

return out;

}
