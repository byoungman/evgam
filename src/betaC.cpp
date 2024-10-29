// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

double digamma(double x) {
  return R::digamma(x);
}

double trigamma(double x) {
  return R::trigamma(x);
}

double psigamma(double x, int y) {
  return R::psigamma(x, y);
}

// //' Beta distribution negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each Beta parameter
// //' @param X1 a design matrix for the first Beta shape location parameter
// //' @param X2 a design matrix for the second Beta shape parameter
// //' @param ymat a matrix
// //' @param dupid a scalar or vector, identifying duplicates in X1 and X2
// //' @return betad0 a scalar, the negative log-likelihood
// //' @return betad12 a matrix, first then second derivatives w.r.t. Beta parameters
// //' @return betad34 a matrix, third then fourth derivatives w.r.t. Beta parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double betad0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
  
  arma::vec lshape1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lshape2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  
  int nobs = nhere.size();
  
  if (dcate == 1) {
    lshape1vec = lshape1vec.elem(dupid);
    lshape2vec = lshape2vec.elem(dupid);
  }
  
  double nllh = 0.0;
  double y, lshape1, lshape2, shape1, shape2, lb12;
  
  for (int j=0; j < nobs; j++) {
    
    lshape1 = lshape1vec[j];
    lshape2 = lshape2vec[j];
    shape1 = exp(lshape1);
    shape2 = exp(lshape2);
    lb12 = lgamma(shape1) + lgamma(shape2) - lgamma(shape1 + shape2);

    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      nllh += (1 - shape1) * log(y) + (1 - shape2) * log(1 - y) + lb12;
      
    }
  }
  
  return(nllh);
  
}

// //' @rdname betad0
// [[Rcpp::export]]
arma::mat betad12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{   
  
  arma::vec lshape1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lshape2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  
  int nobs = nhere.size();
  
  if (dcate == 1) {
    lshape1vec = lshape1vec.elem(dupid);
    lshape2vec = lshape2vec.elem(dupid);
  }
  
  arma::mat out = arma::mat(nobs, 5, arma::fill::zeros);
  
  double lshape1, lshape2, y;

  double ee1, ee2, ee3, ee4, ee5, ee6, ee7;
  double ee10, ee11, ee12;

  for (int j=0; j < nobs; j++) {
    
    lshape1 = lshape1vec[j];
    lshape2 = lshape2vec[j];
    
    ee1 = exp(lshape1);
    ee2 = exp(lshape2);
    ee3 = ee1 + ee2;
    ee4 = digamma(ee3);
    ee5 = digamma(ee1);
    ee6 = digamma(ee2);
    ee7 = trigamma(ee3);
    ee10 = R_pow(ee4, 2) + ee7;
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);

      ee11 = log(1 - y);
      ee12 = log(y);

      out(j, 0) += (ee5 - (ee4 + ee12)) * ee1;
      out(j, 1) += (ee6 - (ee4 + ee11)) * ee2;
      out(j, 2) += ((1 - ee5 * ee1) * ee5 + (R_pow(ee5, 2) + trigamma(ee1) -
        ee10) * ee1 - ((1 - ee4 * ee1) * ee4 + ee12)) * ee1;
      out(j, 3) += -(ee1 * ee2 * ee7);
      out(j, 4) += ((1 - ee6 * ee2) * ee6 + (R_pow(ee6, 2) + trigamma(ee2) -
        ee10) * ee2 - ((1 - ee4 * ee2) * ee4 + ee11)) * ee2;
      
      }
  }
  
  return out;
  
}

// //' @rdname betad0
// [[Rcpp::export]]
arma::mat betad34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{   
  
  arma::vec lshape1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lshape2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  
  int nobs = nhere.size();
  
  if (dcate == 1) {
    lshape1vec = lshape1vec.elem(dupid);
    lshape2vec = lshape2vec.elem(dupid);
  }
  
  double lshape1, lshape2, y;

  double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee10;

  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    lshape1 = lshape1vec[j];
    lshape2 = lshape2vec[j];

    ee1 = exp(lshape1);
    ee2 = exp(lshape2);
    ee3 = ee1 + ee2;
    ee4 = psigamma(ee3, 2);
    ee5 = trigamma(ee3);
    ee6 = -(ee1 * (ee1 * ee4 + ee5) * ee2);
    ee7 = -(ee1 * ee2 * (ee2 * ee4 + ee5));
    ee8 = 3 * ee5;
    ee10 = digamma(ee3);
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);

      out(j, 0) += ((3 * trigamma(ee1) + ee1 * (psigamma(ee1, 2) - ee4) - ee8) * ee1 + digamma(ee1) - (ee10 + log(y))) * ee1;
      out(j, 1) += ee6;
      out(j, 2) += ee7;
      out(j, 3) += ((3 * trigamma(ee2) + ee2 * (psigamma(ee2, 2) - ee4) - ee8) * ee2 + digamma(ee2) - (ee10 + log(1 - y))) * ee2;
      
    }
      
  }
  
  return out;
  
}

// //' Beta distribution negative log-likelihood using sparse design matrices
// //'
// //' @param pars a list of vectors of coefficients for each Beta parameter
// //' @param X1 a sparse design matrix for the first Beta shape location parameter
// //' @param X2 a sparse design matrix for the second Beta shape parameter
// //' @param ymat a matrix
// //' @param dupid a scalar or vector, identifying duplicates in X1 and X2
// //' @return betaspd0 a scalar, the negative log-likelihood
// //' @return betaspd12 a matrix, first then second derivatives w.r.t. Beta parameters
// //' @return betaspd34 a matrix, third then fourth derivatives w.r.t. Beta parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double betaspd0(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{

  arma::vec lshape1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lshape2vec = X2 * Rcpp::as<arma::vec>(pars[1]);

  int nobs = nhere.size();

  if (dcate == 1) {
    lshape1vec = lshape1vec.elem(dupid);
    lshape2vec = lshape2vec.elem(dupid);
  }

  double nllh = 0.0;
  double y, lshape1, lshape2, shape1, shape2, lb12;

  for (int j=0; j < nobs; j++) {

    lshape1 = lshape1vec[j];
    lshape2 = lshape2vec[j];
    shape1 = exp(lshape1);
    shape2 = exp(lshape2);
    lb12 = lgamma(shape1) + lgamma(shape2) - lgamma(shape1 + shape2);

    for (int l=0; l < nhere[j]; l++) {

      y = ymat(j, l);
      nllh += (1 - shape1) * log(y) + (1 - shape2) * log(1 - y) + lb12;

    }
  }

  return(nllh);

}

// //' @rdname betaspd0
// [[Rcpp::export]]
arma::mat betaspd12(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{

  arma::vec lshape1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lshape2vec = X2 * Rcpp::as<arma::vec>(pars[1]);

  int nobs = nhere.size();

  if (dcate == 1) {
    lshape1vec = lshape1vec.elem(dupid);
    lshape2vec = lshape2vec.elem(dupid);
  }

  arma::mat out = arma::mat(nobs, 5, arma::fill::zeros);

  double lshape1, lshape2, y;

  double ee1, ee2, ee3, ee4, ee5, ee6, ee7;
  double ee10, ee11, ee12;

  for (int j=0; j < nobs; j++) {

    lshape1 = lshape1vec[j];
    lshape2 = lshape2vec[j];

    ee1 = exp(lshape1);
    ee2 = exp(lshape2);
    ee3 = ee1 + ee2;
    ee4 = digamma(ee3);
    ee5 = digamma(ee1);
    ee6 = digamma(ee2);
    ee7 = trigamma(ee3);
    ee10 = R_pow(ee4, 2) + ee7;

    for (int l=0; l < nhere[j]; l++) {

      y = ymat(j, l);

      ee11 = log(1 - y);
      ee12 = log(y);

      out(j, 0) += (ee5 - (ee4 + ee12)) * ee1;
      out(j, 1) += (ee6 - (ee4 + ee11)) * ee2;
      out(j, 2) += ((1 - ee5 * ee1) * ee5 + (R_pow(ee5, 2) + trigamma(ee1) -
        ee10) * ee1 - ((1 - ee4 * ee1) * ee4 + ee12)) * ee1;
      out(j, 3) += -(ee1 * ee2 * ee7);
      out(j, 4) += ((1 - ee6 * ee2) * ee6 + (R_pow(ee6, 2) + trigamma(ee2) -
        ee10) * ee2 - ((1 - ee4 * ee2) * ee4 + ee11)) * ee2;

    }
  }

  return out;

}

// //' @rdname betaspd0
// [[Rcpp::export]]
arma::mat betaspd34(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{

  arma::vec lshape1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lshape2vec = X2 * Rcpp::as<arma::vec>(pars[1]);

  int nobs = nhere.size();

  if (dcate == 1) {
    lshape1vec = lshape1vec.elem(dupid);
    lshape2vec = lshape2vec.elem(dupid);
  }

  double lshape1, lshape2, y;

  double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee10;

  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);

  for (int j=0; j < nobs; j++) {

    lshape1 = lshape1vec[j];
    lshape2 = lshape2vec[j];

    ee1 = exp(lshape1);
    ee2 = exp(lshape2);
    ee3 = ee1 + ee2;
    ee4 = psigamma(ee3, 2);
    ee5 = trigamma(ee3);
    ee6 = -(ee1 * (ee1 * ee4 + ee5) * ee2);
    ee7 = -(ee1 * ee2 * (ee2 * ee4 + ee5));
    ee8 = 3 * ee5;
    ee10 = digamma(ee3);

    for (int l=0; l < nhere[j]; l++) {

      y = ymat(j, l);

      out(j, 0) += ((3 * trigamma(ee1) + ee1 * (psigamma(ee1, 2) - ee4) - ee8) * ee1 + digamma(ee1) - (ee10 + log(y))) * ee1;
      out(j, 1) += ee6;
      out(j, 2) += ee7;
      out(j, 3) += ((3 * trigamma(ee2) + ee2 * (psigamma(ee2, 2) - ee4) - ee8) * ee2 + digamma(ee2) - (ee10 + log(1 - y))) * ee2;

    }

  }

  return out;

}
