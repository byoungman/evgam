// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// //' Three-parameter Weibull distribution negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each Weibull parameter
// //' @param X1 a design matrix for the Weibull location parameter
// //' @param X2 a design matrix for the Weibull log scale parameter
// //' @param X3 a design matrix for the Weibull log shape parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return weib3d0 a scalar, the negative log-likelihood
// //' @return weib3d12 a matrix, first then second derivatives w.r.t. Weibull parameters
// //' @return weib3d34 a matrix, third then fourth derivatives w.r.t. Weibull parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double weib3d0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
    
arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec p3vec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = nhere.size();

if (dcate == 1) {
  p1vec = p1vec.elem(dupid);
  p2vec = p2vec.elem(dupid);
  p3vec = p3vec.elem(dupid);
}

double y, pars1, pars2, pars3, lambda, kappa, res;
double nllh = 0.0;

for (int j=0; j < nobs; j++) {
  
  pars1 = p1vec[j];
  pars2 = p2vec[j];
  pars3 = p3vec[j];
  lambda = exp(pars2);
  kappa = exp(pars3);

  for (int l=0; l < nhere[j]; l++) {
    
    y = ymat(j, l);
    
    res = y - pars1;
    
    if (res <= 0) {
      nllh = 1e20;
      break;
    }
    
    nllh -= log(kappa) - kappa * log(lambda) + (kappa - 1) * log(res) - R_pow(res / lambda, kappa);
  
  }
  
}

return(nllh);

}

// //' @rdname weib3d0
// [[Rcpp::export]]
arma::mat weib3d12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
    
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec p3vec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
    p3vec = p3vec.elem(dupid);
  }
  
  double y, pars1, pars2, pars3;
  
  double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13;
  
  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];
    pars3 = p3vec[j];

    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      
      ee1 = exp(pars3);
      ee2 = y - pars1;
      ee3 = exp(pars2);
      ee4 = ee2/ee3;
      ee5 = ee1 - 1;
      ee6 = log(ee2);
      ee7 = R_pow(ee4, ee5);
      ee8 = ee6 - pars2;
      ee9 = R_pow(ee4, (ee1 - 2));
      ee10 = 1 + ee1 * ee8;
      ee11 = (ee7 + ee9 * ee5 * ee2/ee3) * ee1;
      ee12 = ee7 * ee10;
      ee13 = R_pow(ee4, ee1);
      
      out(j, 0) += ee5/ee2 - ee7 * ee1/ee3;
      out(j, 1) += (1 - ee7 * ee2/ee3) * ee1;
      out(j, 2) += (ee13 * ee8 + pars2 - ee6) * ee1 - 1;
      out(j, 3) += (ee9 * ee1/R_pow(ee3, 2) + 1/R_pow(ee2, 2)) * ee5;
      out(j, 4) += ee11/ee3;
      out(j, 5) += (1/ee2 - ee12/ee3) * ee1;
      out(j, 6) += ee11 * ee2/ee3;
      out(j, 7) += (1 - ee12 * ee2/ee3) * ee1;
      out(j, 8) += (ee13 * ee10 * ee8 + pars2 - ee6) * ee1;
      
    }

}

return out;

}

// //' @rdname weib3d0
// [[Rcpp::export]]
arma::mat weib3d34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
    
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec p3vec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
    p3vec = p3vec.elem(dupid);
  }
  
  double y, pars1, pars2, pars3;
  
  double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee17;
  double ee20, ee21, ee22, ee23, ee26, ee27, ee29;
  double ee30, ee31, ee32, ee35, ee36, ee37, ee38, ee39;
  double ee40, ee41, ee42, ee43, ee44, ee47, ee48, ee49;
  
  arma::mat out = arma::mat(nobs, 25, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];
    pars3 = p3vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      
      ee1 = exp(pars3);
      ee2 = y - pars1;
      ee3 = exp(pars2);
      ee4 = log(ee2);
      ee5 = ee2/ee3;
      ee6 = ee4 - pars2;
      ee7 = ee1 - 1;
      ee8 = ee1 - 2;
      ee9 = R_pow(ee5, ee8);
      ee10 = ee1 * ee6;
      ee11 = R_pow(ee5, ee7);
      ee12 = ee1 - 3;
      ee13 = R_pow(ee5, ee12);
      ee14 = ee7 * ee6;
      ee15 = 1 + ee10;
      ee17 = (ee14 + 2) * ee1 - 1;
      ee20 = (3 + ee10) * ee1 * ee6 + 1;
      ee21 = R_pow(ee3, 2);
      ee22 = ee17 * ee9;
      ee23 = ee20 * ee11;
      ee26 = (ee8 * ee6 + 1) * ee7 + ee1 - 2;
      ee27 = ee11 * ee15;
      ee29 = 2 * (1 + 2 * ee7);
      ee30 = (((ee15 * ee7 + 4 * ee1 - 2) * ee6 + 4) * ee1 - 1) *  ee9;
      ee31 = (ee26 * ee13 * ee2/ee3 + (ee14 + 1) * ee9) * ee1;
      ee32 = (ee9 + ee13 * ee8 * ee2/ee3) * ee7;
      ee35 = ((6 + ee10) * ee1 * ee6 + 7) * ee1 * ee6 + 1;
      ee36 = 1 + 2 * ee8;
      ee37 = R_pow(ee3, 3);
      ee38 = ((ee31 + ee17 * (2 * ee9) + ee32) * ee2/ee3 + ee27) *  ee1;
      ee39 = (ee30 * ee2/ee3 + ee23) * ee1;
      ee40 = (ee22 * ee2/ee3 + ee27) * ee1;
      ee41 = ee35 * ee11;
      ee42 = ((ee5 + 2) * ee36 * ee3/ee2 + 3) * ee7;
      ee43 = ee13 * ee1;
      ee44 = R_pow(ee5, ee1);
      ee47 = 1/ee2;
      ee48 = 1/R_pow(ee2, 2);
      ee49 = 2/R_pow(ee2, 3);
      
      out(j, 0) += (ee49 - ee43 * ee8/ee37) * ee7;
      out(j, 1) += -(ee9 * ee1 * ee1 * ee7/ee21);
      out(j, 2) += (ee22/ee21 + ee48) * ee1;
      out(j, 3) += -((ee7/ee3 + ee29) * ee9 * ee2/ee3 * ee1);
      out(j, 4) += ee40/ee3;
      out(j, 5) += (ee47 - ee23/ee3) * ee1;
      out(j, 6) += -((ee7 * ee2/ee3 + ee29) * ee9 * ee2/ee3 * ee1);
      out(j, 7) += ee40 * ee2/ee3;
      out(j, 8) += (1 - ee23 * ee2/ee3) * ee1;
      out(j, 9) += (ee20 * ee44 * ee6 + pars2 - ee4) * ee1;
      out(j, 10) += (R_pow(ee5, (ee1 - 4)) * ee1 * ee8 * ee12/R_pow(ee3, 4) +
        6/R_pow(ee2, 4)) * ee7;
      out(j, 11) += ee43 * ee1 * ee7 * ee8/ee37;
      out(j, 12) += (ee49 - (ee26 * ee1 + ee7 * ee8) * ee13/ee37) *
        ee1;
      out(j, 13) += (ee36 * (1/ee3 + 2) * ee3 + 2) * ee9 * ee7/ee21 *
        ee1;
      out(j, 14) += -((ee31 + ee22 + ee32) * ee1/ee21);
      out(j, 15) += (ee30/ee21 + ee48) * ee1;
      out(j, 16) += (ee42/ee3 + ee29) * ee9 * ee1 * ee2/ee3;
      out(j, 17) += -(ee38/ee3);
      out(j, 18) += ee39/ee3;
      out(j, 19) += (ee47 - ee41/ee3) * ee1;
      out(j, 20) += (ee42 * ee2/ee3 + ee29) * ee9 * ee1 * ee2/ee3;
      out(j, 21) += -(ee38 * ee2/ee3);
      out(j, 22) += ee39 * ee2/ee3;
      out(j, 23) += (1 - ee41 * ee2/ee3) * ee1;
      out(j, 24) += (ee35 * ee44 * ee6 + pars2 - ee4) * ee1;
      
    }

}

return out;

}
