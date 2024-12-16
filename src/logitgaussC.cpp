// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double logtwopi = log(2.0 * M_PI);

// //' Logit Gaussian distribution negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each Beta parameter
// //' @param X1 a design matrix for the location parameter
// //' @param X2 a design matrix for the log scale parameter
// //' @param ymat a matrix
// //' @param dupid a scalar or vector, identifying duplicates in X1 and X2
// //' @return logitgaussd0 a scalar, the negative log-likelihood
// //' @return logitgaussd12 a matrix, first then second derivatives w.r.t. logit Gaussian parameters
// //' @return logitgaussd34 a matrix, third then fourth derivatives w.r.t. logit Gaussian parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double logitgaussd0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
  
  arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lsigmavec = X2 * Rcpp::as<arma::vec>(pars[1]);
  
  int nobs = nhere.size();
  
  if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lsigmavec = lsigmavec.elem(dupid);
  }
  
  double nllh = 0.0;
  double y, mu, lsigma, sigma, res;
  
  for (int j=0; j < nobs; j++) {
    
    mu = muvec[j];
    lsigma = lsigmavec[j];
    sigma = exp(lsigma);
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      res = (1.0 / (1.0 + exp(-y)) - mu) / sigma;
      nllh += 0.5 * res * res + 0.5 * logtwopi + log(1 - y) + log(y) + lsigma;
      
    }
  }
  
  return(nllh);
  
}

// //' @rdname logitgaussd0
// [[Rcpp::export]]
arma::mat logitgaussd12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{   
  
  arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lsigmavec = X2 * Rcpp::as<arma::vec>(pars[1]);
  
  int nobs = nhere.size();
  
  if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lsigmavec = lsigmavec.elem(dupid);
  }
  
  arma::mat out = arma::mat(nobs, 5, arma::fill::zeros);
  
  double y, mu, lsigma, sigma;
  double ee2, ee6, ee7, ee9;
  
  for (int j=0; j < nobs; j++) {
    
    mu = muvec[j];
    lsigma = lsigmavec[j];
    sigma = exp(lsigma);

    ee2 = R_pow(exp(lsigma), 2);
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);

      ee6 = 1/(1 + exp(-y)) - mu;
      ee7 = ee6/ee2;
      ee9 = R_pow(ee6, 2)/ee2;

      out(j, 0) += -ee7;
      out(j, 1) += 1 - ee9;
      out(j, 2) += 1/ee2;
      out(j, 3) += 2 * ee7;
      out(j, 4) += 2 * ee9;
      
      }
  }
  
  return out;
  
}

// //' @rdname logitgaussd0
// [[Rcpp::export]]
arma::mat logitgaussd34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{   
  
  arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lsigmavec = X2 * Rcpp::as<arma::vec>(pars[1]);
  
  int nobs = nhere.size();
  
  if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lsigmavec = lsigmavec.elem(dupid);
  }
  
  double y, mu, lsigma, sigma;
  double ee2, ee6, ee7, ee9;

  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    mu = muvec[j];
    lsigma = lsigmavec[j];
    sigma = exp(lsigma);
    
    ee2 = R_pow(exp(lsigma), 2);

    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);

      ee6 = 1/(1 + exp(-y)) - mu;
      ee7 = ee6/ee2;
      ee9 = R_pow(ee6, 2)/ee2;
      
      out(j, 0) += 0;
      out(j, 1) += -(2/ee2);
      out(j, 2) += -(4 * ee7);
      out(j, 3) += -(4 * ee9);
      out(j, 4) += 0;
      out(j, 5) += 0;
      out(j, 6) += 4/ee2;
      out(j, 7) += 8 * ee7;
      out(j, 8) += 8 * ee9;
      
    }
      
  }
  
  return out;
  
}

// //' Logit Gaussian distribution negative log-likelihood using sparse matrices
// //'
// //' @param pars a list of vectors of coefficients for each Beta parameter
// //' @param X1 a sparse design matrix for the location parameter
// //' @param X2 a sparse design matrix for the log scale parameter
// //' @param ymat a matrix
// //' @param dupid a scalar or vector, identifying duplicates in X1 and X2
// //' @return logitgaussspd0 a scalar, the negative log-likelihood
// //' @return logitgaussspd12 a matrix, first then second derivatives w.r.t. logit Gaussian parameters
// //' @return logitgaussspd34 a matrix, third then fourth derivatives w.r.t. logit Gaussian parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double logitgaussspd0(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
  
  arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lsigmavec = X2 * Rcpp::as<arma::vec>(pars[1]);
  
  int nobs = nhere.size();
  
  if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lsigmavec = lsigmavec.elem(dupid);
  }
  
  double nllh = 0.0;
  double y, mu, lsigma, sigma, res;
  
  for (int j=0; j < nobs; j++) {
    
    mu = muvec[j];
    lsigma = lsigmavec[j];
    sigma = exp(lsigma);
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      res = (1.0 / (1.0 + exp(-y)) - mu) / sigma;
      nllh += 0.5 * res * res + 0.5 * logtwopi + log(1 - y) + log(y) + lsigma;
      
    }
  }
  
  return(nllh);
  
}

// //' @rdname logitgaussspd0
// [[Rcpp::export]]
arma::mat logitgaussspd12(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{   
  
  arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lsigmavec = X2 * Rcpp::as<arma::vec>(pars[1]);
  
  int nobs = nhere.size();
  
  if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lsigmavec = lsigmavec.elem(dupid);
  }
  
  arma::mat out = arma::mat(nobs, 5, arma::fill::zeros);
  
  double y, mu, lsigma, sigma;
  double ee2, ee6, ee7, ee9;
  
  for (int j=0; j < nobs; j++) {
    
    mu = muvec[j];
    lsigma = lsigmavec[j];
    sigma = exp(lsigma);
    
    ee2 = R_pow(exp(lsigma), 2);
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      
      ee6 = 1/(1 + exp(-y)) - mu;
      ee7 = ee6/ee2;
      ee9 = R_pow(ee6, 2)/ee2;
      
      out(j, 0) += -ee7;
      out(j, 1) += 1 - ee9;
      out(j, 2) += 1/ee2;
      out(j, 3) += 2 * ee7;
      out(j, 4) += 2 * ee9;
      
    }
  }
  
  return out;
  
}

// //' @rdname logitgaussspd0
// [[Rcpp::export]]
arma::mat logitgaussspd34(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{   
  
  arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lsigmavec = X2 * Rcpp::as<arma::vec>(pars[1]);
  
  int nobs = nhere.size();
  
  if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lsigmavec = lsigmavec.elem(dupid);
  }
  
  double y, mu, lsigma, sigma;
  double ee2, ee6, ee7, ee9;
  
  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    mu = muvec[j];
    lsigma = lsigmavec[j];
    sigma = exp(lsigma);
    
    ee2 = R_pow(exp(lsigma), 2);
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      
      ee6 = 1/(1 + exp(-y)) - mu;
      ee7 = ee6/ee2;
      ee9 = R_pow(ee6, 2)/ee2;
      
      out(j, 0) += 0;
      out(j, 1) += -(2/ee2);
      out(j, 2) += -(4 * ee7);
      out(j, 3) += -(4 * ee9);
      out(j, 4) += 0;
      out(j, 5) += 0;
      out(j, 6) += 4/ee2;
      out(j, 7) += 8 * ee7;
      out(j, 8) += 8 * ee9;
      
    }
    
  }
  
  return out;
  
}
