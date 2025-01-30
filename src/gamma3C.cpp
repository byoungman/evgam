// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// //' Three-parameter Gamma distribution negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each Gamma parameter
// //' @param X1 a design matrix for the Gamma location parameter
// //' @param X2 a design matrix for the Gamma log scale parameter
// //' @param X3 a design matrix for the Gamma log shape parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return gamma3d0 a scalar, the negative log-likelihood
// //' @return gamma3d12 a matrix, first then second derivatives w.r.t. Gamma parameters
// //' @return gamma3d34 a matrix, third then fourth derivatives w.r.t. Gamma parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double gamma3d0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
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

double y, pars1, pars2, pars3, theta, alpha, res;
double nllh = 0.0;

for (int j=0; j < nobs; j++) {
  
  pars1 = p1vec[j];
  pars2 = p2vec[j];
  pars3 = p3vec[j];
  theta = exp(pars2);
  alpha = exp(pars3);

  for (int l=0; l < nhere[j]; l++) {
    
    y = ymat(j, l);
    
    res = y - pars1;
    
    if (res <= 0) {
      nllh = 1e20;
      break;
    }
    
    nllh += (1 - alpha) * log(res) + (res / theta) + alpha * pars2 + lgamma(alpha);
  
  }
  
}

return(nllh);

}

// //' @rdname gamma3d0
// [[Rcpp::export]]
arma::mat gamma3d12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
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
  
  double ee1, ee2, ee3, ee4, ee5, ee6, ee7;
  
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
      ee4 = 1 - ee1;
      ee5 = 1/ee3;
      ee6 = R::digamma(ee1);
      ee7 = log(ee2);
      
      out(j, 0) += -(ee4/ee2 + ee5);
      out(j, 1) += ee1 - ee2/ee3;
      out(j, 2) += (ee6 + pars2 - ee7) * ee1;
      out(j, 3) += -(ee4/R_pow(ee2, 2));
      out(j, 4) += ee5;
      out(j, 5) += ee1/ee2;
      out(j, 6) += 1 * ee2/ee3;
      out(j, 7) += ee1;
      out(j, 8) += (ee6 + ee1 * R::trigamma(ee1) + pars2 - ee7) * ee1;
      
    }

}

return out;

}

// //' @rdname gamma3d0
// [[Rcpp::export]]
arma::mat gamma3d34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
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
  
  double ee1, ee2, ee3, ee5, ee6, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15;
  
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
      ee5 = R_pow(ee2, 3);
      ee6 = 1 - ee1;
      ee8 = 1 * ee2/ee3;
      ee9 = 1/ee3;
      ee10 = R::digamma(ee1);
      ee11 = ee1/ee2;
      ee12 = ee1/R_pow(ee2, 2);
      ee13 = log(ee2);
      ee14 = R::psigamma(ee1, 2);
      ee15 = R::trigamma(ee1);
      
      out(j, 0) += -(2 * (ee6/ee5));
      out(j, 1) += 0;
      out(j, 2) += ee12;
      out(j, 3) += -ee9;
      out(j, 4) += 0;
      out(j, 5) += ee11;
      out(j, 6) += -ee8;
      out(j, 7) += 0;
      out(j, 8) += ee1;
      out(j, 9) += ((3 * ee15 + ee1 * ee14) * ee1 + ee10 + pars2 -
        ee13) * ee1;
      out(j, 10) += -(6 * (ee6/R_pow(ee2, 4)));
      out(j, 11) += 0;
      out(j, 12) += 2 * (ee1/ee5);
      out(j, 13) += 0;
      out(j, 14) += 0;
      out(j, 15) += ee12;
      out(j, 16) += ee9;
      out(j, 17) += 0;
      out(j, 18) += 0;
      out(j, 19) += ee11;
      out(j, 20) += ee8;
      out(j, 21) += 0;
      out(j, 22) += 0;
      out(j, 23) += ee1;
      out(j, 24) += (((6 * ee14 + ee1 * R::psigamma(ee1, 3)) * ee1 +
        7 * ee15) * ee1 + ee10 + pars2 - ee13) * ee1;
      
    }

}

return out;

}
