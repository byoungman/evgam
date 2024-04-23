// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// //' Asymmetric Laplace distribution (ALD) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each ALD parameter
// //' @param X1 a design matrix for the ALD location parameter
// //' @param X2 a design matrix for the ALD log scale parameter
// //' @param yvec a vector
// //' @param tau a scalar, the quantile sought
// //' @param C a scalar, for the Ho et al. (2000) correction
// //' @param dupid a scalar or vector, identifying duplicates in X1 and X2
// //' @return aldd0 a scalar, the negative log-liklihood
// //' @return aldd12 a matrix, first then second derivatives w.r.t. ALD parameters
// //' @return aldd34 a matrix, third then fourth derivatives w.r.t. ALD parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double aldd0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::vec tau, arma::vec C, arma::uvec dupid, int dcate, arma::uvec nhere)
{
  
  arma::vec mu = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lsigma = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    mu = mu.elem(dupid);
    lsigma = lsigma.elem(dupid);
  }
  
  double nllh=0.0, res;
  
  for (int j=0; j < nobs; j++) {
    
    for (int l=0; l < nhere[j]; l++) {
      
    res = (ymat(j, l) - mu[j]) / exp(lsigma[j]);
    nllh += lsigma[j];
    
    if (res <= -C[0]) {
      
      nllh += (tau[0] - 1.0) * (2.0 * res + C[0]);
      
    } else {
      
      if (res < 0.0) {
        
        nllh += (1.0 - tau[0]) * res * res / C[0];
        
      } else {
        
        if (res <= C[0]) {
          
          nllh += tau[0] * res * res / C[0];
          
        } else {
          
          nllh += tau[0] * (2.0 * res - C[0]);
          
        }}}
    }
    
  }
  
  return(nllh);
  
}

// //' @rdname aldd0
// [[Rcpp::export]]
arma::mat aldd12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::vec tau, arma::vec C, arma::uvec dupid, int dcate, arma::uvec nhere)
{   
  
  arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  arma::mat out = arma::mat(nobs, 5, arma::fill::zeros);
  
  if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
  }
  
  double ee1, ee2, ee3, ee5, ee6, ee7;
  double y, mu, lpsi, res;
  
  for (int j=0; j < nobs; j++) {
    
    mu = muvec[j];
    lpsi = lpsivec[j];

    for (int l=0; l < nhere[j]; l++) {
      
    y = ymat(j, l);
    res = (y - mu) / exp(lpsi);
    
    if (res <= -C[0]) {
      
      ee1 = exp(lpsi);
      ee2 = tau[0] - 1;
      ee6 = 2 * (ee2 * (y - mu)/ee1);
      ee7 = 2 * (ee2/ee1);
      
      out(j, 0) += -ee7;
      out(j, 1) += 1 - ee6;
      
      out(j, 2) += 0;
      out(j, 3) += ee7;
      out(j, 4) += ee6;
      
    } else {
      
      if (res < 0.0) {
        
        ee1 = 1 - tau[0];
        ee2 = C[0] * exp(2 * lpsi);
        ee3 = y - mu;
        ee5 = ee1 * ee3/ee2;
        ee7 = ee1 * ee3 * ee3/ee2;
        
        out(j, 0) += -(2 * ee5);
        out(j, 1) += 1 - 2 * ee7;
        
        out(j, 2) += 2 * (ee1/ee2);
        out(j, 3) += 4 * ee5;
        out(j, 4) += 4 * ee7;
        
      } else {
        
        if (res <= C[0]) {
          
          ee1 = C[0] * exp(2 * lpsi);
          ee2 = y - mu;
          ee5 = tau[0] * ee2/ee1;
          ee7 = tau[0] * ee2 * ee2/ee1;
          
          out(j, 0) += -(2 * ee5);
          out(j, 1) += 1 - 2 * ee7;
          
          out(j, 2) += 2 * (tau[0]/ee1);
          out(j, 3) += 4 * ee5;
          out(j, 4) += 4 * ee7;
          
        } else {
          
          ee1 = exp(lpsi);
          ee2 = 2 * (tau[0] * (y - mu)/ee1);
          ee3 = 2 * (tau[0]/ee1);
          
          out(j, 0) += -ee3;
          out(j, 1) += 1 - ee2;
          
          out(j, 2) += 0;
          out(j, 3) += ee3;
          out(j, 4) += ee2;
          
        }}}
    }
  }
  
  return out;
  
}

// //' @rdname aldd0
// [[Rcpp::export]]
arma::mat aldd34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::vec tau, arma::vec C, arma::uvec dupid, int dcate, arma::uvec nhere)
{   
  
  arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  
  if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
  }
  
  double ee1, ee2, ee3, ee5, ee6, ee7, ee8;
  double y, mu, lpsi, res;
  
  for (int j=0; j < nobs; j++) {
    
    mu = muvec[j];
    lpsi = lpsivec[j];

    for (int l=0; l < nhere[j]; l++) {
      
    y = ymat(j, l);
    res = (y - mu) / exp(lpsi);
    
    if (res <= -C[0]) {
      
      ee1 = exp(lpsi);
      ee2 = tau[0] - 1;
      ee6 = 2 * (ee2 * (y - mu)/ee1);
      ee7 = 2 * (ee2/ee1);
      
      out(j, 0) += 0;
      out(j, 1) += 0;
      out(j, 2) += -ee7;
      out(j, 3) += -ee6;
      
      out(j, 4) += 0;
      out(j, 5) += 0;
      out(j, 6) += 0;
      out(j, 7) += ee7;
      out(j, 8) += ee6;
      
    } else {
      
      if (res < 0.0) {
        
        ee1 = 1 - tau[0];
        ee2 = C[0] * exp(2 * lpsi);
        ee3 = y - mu;
        ee5 = ee1 * ee3/ee2;
        ee7 = ee1 * ee3 * ee3/ee2;
        ee8 = ee1/ee2;
        
        out(j, 0) += 0;
        out(j, 1) += -(4 * ee8);
        out(j, 2) += -(8 * ee5);
        out(j, 3) += -(8 * ee7);
        
        out(j, 4) += 0;
        out(j, 5) += 0;
        out(j, 6) += 8 * ee8;
        out(j, 7) += 16 * ee5;
        out(j, 8) += 16 * ee7;
        
      } else {
        
        if (res <= C[0]) {
          
          ee1 = C[0] * exp(2 * lpsi);
          ee2 = y - mu;
          ee5 = tau[0] * ee2/ee1;
          ee7 = tau[0] * ee2 * ee2/ee1;
          ee8 = tau[0]/ee1;
          
          out(j, 0) += 0;
          out(j, 1) += -(4 * ee8);
          out(j, 2) += -(8 * ee5);
          out(j, 3) += -(8 * ee7);
          
          out(j, 4) += 0;
          out(j, 5) += 0;
          out(j, 6) += 8 * ee8;
          out(j, 7) += 16 * ee5;
          out(j, 8) += 16 * ee7;
          
        } else {
          
          ee1 = exp(lpsi);
          ee2 = 2 * (tau[0] * (y - mu)/ee1);
          ee3 = 2 * (tau[0]/ee1);
          
          out(j, 0) += 0;
          out(j, 1) += 0;
          out(j, 2) += -ee3;
          out(j, 3) += -ee2;
          
          out(j, 4) += 0;
          out(j, 5) += 0;
          out(j, 6) += 0;
          out(j, 7) += ee3;
          out(j, 8) += ee2;
          
        }}}
    
  }
    
  }
  
  return out;
  
}
