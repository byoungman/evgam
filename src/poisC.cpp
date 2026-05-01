// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// //' Poisson distribution negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for Poisson log mean
// //' @param X1 a design matrix for the Poisson log mean parameter
// //' @param ymat a matrix
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return poisd0 a scalar, the negative log-likelihood
// //' @return poisd12 a matrix, first then second derivatives w.r.t. Poisson parameters
// //' @return poisd34 a matrix, third then fourth derivatives w.r.t. Poisson parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double poisd0(Rcpp::List pars, arma::mat X1, arma::mat ymat, arma::uvec dupid, 
              int dcate, arma::ivec nhere, arma::mat wmat, arma::mat offmat)
{
  
  arma::vec lmuvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    lmuvec = lmuvec.elem(dupid);
  }
  
  double y, lmu, mu, w, offset;
  double nllh = 0.0;
  
  for (int j=0; j < nobs; j++) {
    
    lmu = lmuvec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      w = wmat(j, l);
      offset = offmat(j, l);
      mu = exp(log(offset) + lmu);
      
      nllh += w * (mu - y * lmu);
        
    } 
    
  }
  
  return(nllh);
  
}

// //' @rdname poisd0
// [[Rcpp::export]]
arma::mat poisd12(Rcpp::List pars, arma::mat X1, arma::mat ymat, arma::uvec dupid, 
                  int dcate, arma::ivec nhere, arma::mat wmat, arma::mat offmat)
{
  
  arma::vec lmuvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    lmuvec = lmuvec.elem(dupid);
  }
  
  double y, lmu, mu, w, offset;
  arma::mat out = arma::mat(nobs, 2, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    lmu = lmuvec[j];

    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      w = wmat(j, l);
      offset = offmat(j, l);
      mu = exp(log(offset) + lmu);
      
      out(j, 0) += w * (mu - y);
      out(j, 1) += w * mu;

    }
    
  }
  
  return out;
  
}

// //' @rdname poisd0
// [[Rcpp::export]]
arma::mat poisd34(Rcpp::List pars, arma::mat X1, arma::mat ymat, arma::uvec dupid, 
                  int dcate, arma::ivec nhere, arma::mat wmat, arma::mat offmat)
{
  
  arma::vec lmuvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    lmuvec = lmuvec.elem(dupid);
  }
  
  double lmu, mu, w, offset;
  arma::mat out = arma::mat(nobs, 2, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    lmu = lmuvec[j];

    for (int l=0; l < nhere[j]; l++) {
      
      w = wmat(j, l);
      offset = offmat(j, l);
      mu = exp(log(offset) + lmu);
      
      out(j, 0) += w * mu;
      out(j, 1) += w * mu;
      
    }
    
  }
  
  return out;
  
}

// //' @param pars a list of vectors of coefficients for Poisson log mean
// //' @param X1 a sparse design matrix for the Poisson log mean parameter
// //' @param ymat a matrix
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return poisd0 a scalar, the negative log-likelihood
// //' @return poisd12 a matrix, first then second derivatives w.r.t. Poisson parameters
// //' @return poisd34 a matrix, third then fourth derivatives w.r.t. Poisson parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double poisspd0(Rcpp::List pars, arma::sp_mat X1, arma::mat ymat, arma::uvec dupid, 
                int dcate, arma::ivec nhere, arma::mat wmat, arma::mat offmat)
{
  
  arma::vec lmuvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    lmuvec = lmuvec.elem(dupid);
  }
  
  double y, lmu, mu, w, offset;
  double nllh = 0.0;
  
  for (int j=0; j < nobs; j++) {
    
    lmu = lmuvec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      w = wmat(j, l);
      offset = offmat(j, l);
      mu = exp(log(offset) + lmu);
      
      nllh += w * (mu - y * lmu);
      
    } 
    
  }
  
  return(nllh);
  
}

// //' @rdname poisspd0
// [[Rcpp::export]]
arma::mat poisspd12(Rcpp::List pars, arma::sp_mat X1, arma::mat ymat, arma::uvec dupid, 
                    int dcate, arma::ivec nhere, arma::mat wmat, arma::mat offmat)
{
  
  arma::vec lmuvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    lmuvec = lmuvec.elem(dupid);
  }
  
  double y, lmu, mu, w, offset;
  arma::mat out = arma::mat(nobs, 2, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    lmu = lmuvec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      w = wmat(j, l);
      offset = offmat(j, l);
      mu = exp(log(offset) + lmu);
      
      out(j, 0) += w * (mu - y);
      out(j, 1) += w * mu;
      
    }
    
  }
  
  return out;
  
}

// //' @rdname poisspd0
// [[Rcpp::export]]
arma::mat poisspd34(Rcpp::List pars, arma::sp_mat X1, arma::mat ymat, arma::uvec dupid, 
                    int dcate, arma::ivec nhere, arma::mat wmat, arma::mat offmat)
{
  
  arma::vec lmuvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    lmuvec = lmuvec.elem(dupid);
  }
  
  double lmu, mu, w, offset;
  arma::mat out = arma::mat(nobs, 2, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    lmu = lmuvec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      w = wmat(j, l);
      offset = offmat(j, l);
      mu = exp(log(offset) + lmu);
      
      out(j, 0) += w * mu;
      out(j, 1) += w * mu;
      
    }
    
  }
  
  return out;
  
}
