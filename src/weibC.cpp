// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// //' Two-parameter Weibull distribution negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each Weibull parameter
// //' @param X1 a design matrix for the Weibull log scale parameter
// //' @param X2 a design matrix for the Weibull log shape parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return weibd0 a scalar, the negative log-likelihood
// //' @return weibd12 a matrix, first then second derivatives w.r.t. Weibull parameters
// //' @return weibd34 a matrix, third then fourth derivatives w.r.t. Weibull parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double weibd0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
    
arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = nhere.size();

if (dcate == 1) {
  p1vec = p1vec.elem(dupid);
  p2vec = p2vec.elem(dupid);
}

double y, pars1, pars2;

double lambda, kappa;
double nllh = 0.0;

for (int j=0; j < nobs; j++) {
  
  pars1 = p1vec[j];
  pars2 = p2vec[j];
  lambda = exp(pars1);
  kappa = exp(pars2);

  for (int l=0; l < nhere[j]; l++) {
    
    y = ymat(j, l);
    
    nllh += kappa * pars1 - pars2 + (1 - kappa) * log(y) + R_pow(y / lambda, kappa);
  
  }
  
}

return(nllh);

}

// //' @rdname weibd0
// [[Rcpp::export]]
arma::mat weibd12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
    
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, pars1, pars2;
  
  double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  
  arma::mat out = arma::mat(nobs, 5, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];

    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      
      ee1 = exp(pars2);
      ee2 = exp(pars1);
      ee3 = log(y);
      ee4 = y/ee2;
      ee5 = ee1 - 1;
      ee6 = ee3 - pars1;
      ee7 = R_pow(ee4, ee5);
      ee8 = R_pow(ee4, ee1);
      ee9 = 1 + ee1 * ee6;
      
      out(j, 0) += (1 - y * ee7/ee2) * ee1;
      out(j, 1) += (ee6 * ee8 + pars1 - ee3) * ee1 - 1;
      out(j, 2) += y * (ee7 + y * ee5 * R_pow(ee4, (ee1 - 2))/ee2) *
        ee1/ee2;
      out(j, 3) += (1 - y * ee9 * ee7/ee2) * ee1;
      out(j, 4) += (ee9 * ee6 * ee8 + pars1 - ee3) * ee1;
      
    }

}

return out;

}

// //' @rdname weibd0
// [[Rcpp::export]]
arma::mat weibd34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
    
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, pars1, pars2;
  
  double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee13, ee14, ee15, ee18;
  double ee20, ee21, ee22, ee23, ee25;
  
  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];

    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);

      ee1 = exp(pars2);
      ee2 = exp(pars1);
      ee3 = log(y);
      ee4 = ee3 - pars1;
      ee5 = y/ee2;
      ee6 = ee1 - 1;
      ee7 = ee1 - 2;
      ee8 = ee1 * ee4;
      ee9 = R_pow(ee5, ee7);
      ee10 = R_pow(ee5, ee6);
      ee13 = (3 + ee8) * ee1 * ee4 + 1;
      ee14 = ee6 * ee4;
      ee15 = 1 + ee8;
      ee18 = ((6 + ee8) * ee1 * ee4 + 7) * ee1 * ee4 + 1;
      ee20 = (ee14 + 2) * ee1 - 1;
      ee21 = ee15 * ee10;
      ee22 = R_pow(ee5, (ee1 - 3));
      ee23 = R_pow(ee5, ee1);
      ee25 = 2 * (1 + 2 * ee6);
      
      out(j, 0) += -(ee1 * (y * ((ee25 + y * ee6/ee2) * ee9)/ee2));
      out(j, 1) += y * (ee21 + y * ee20 * ee9/ee2) * ee1/ee2;
      out(j, 2) += (1 - y * ee13 * ee10/ee2) * ee1;
      out(j, 3) += (ee13 * ee4 * ee23 + pars1 - ee3) * ee1;
      out(j, 4) += y * ((ee25 + y * ((1 + 2 * ee7) * (2 + ee5) * ee2/
        y + 3) * ee6/ee2) * ee9) * ee1/ee2;
      out(j, 5) += -(y * (ee21 + y * (((ee14 + 1) * ee9 + y * ((ee7 *
        ee4 + 1) * ee6 + ee1 - 2) * ee22/ee2) * ee1 + ee20 * (2 *
        ee9) + (ee9 + y * ee7 * ee22/ee2) * ee6)/ee2) * ee1/ee2);
      out(j, 6) += y * (ee13 * ee10 + y * (((ee15 * ee6 + 4 * ee1 -
        2) * ee4 + 4) * ee1 - 1) * ee9/ee2) * ee1/ee2;
      out(j, 7) += (1 - y * ee18 * ee10/ee2) * ee1;
      out(j, 8) += (ee18 * ee4 * ee23 + pars1 - ee3) * ee1;
      
    }

}

return out;

}
