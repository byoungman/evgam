// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// //' Gamma distribution negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each Gamma parameter
// //' @param X1 a design matrix for the Gamma log scale parameter
// //' @param X2 a design matrix for the Gamma log shape parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return gammad0 a scalar, the negative log-likelihood
// //' @return gammad12 a matrix, first then second derivatives w.r.t. gammaull parameters
// //' @return gammad34 a matrix, third then fourth derivatives w.r.t. gammaull parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double gammad0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
    
arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = nhere.size();

if (dcate == 1) {
  p1vec = p1vec.elem(dupid);
  p2vec = p2vec.elem(dupid);
}

double y, pars1, pars2;

double theta, alpha;
double nllh = 0.0;

for (int j=0; j < nobs; j++) {
  
  pars1 = p1vec[j];
  pars2 = p2vec[j];
  theta = exp(pars1);
  alpha = exp(pars2);

  for (int l=0; l < nhere[j]; l++) {
    
    y = ymat(j, l);
    
    nllh += (1 - alpha) * log(y) + (y / theta) + alpha * pars1 + lgamma(alpha);
  
  }
  
}

return(nllh);

}

// //' @rdname gammad0
// [[Rcpp::export]]
arma::mat gammad12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
    
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, pars1, pars2;
  
  double ee1, ee2, ee3, ee4;
  
  arma::mat out = arma::mat(nobs, 5, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];

    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      
      ee1 = exp(pars2);
      ee2 = R::digamma(ee1);
      ee3 = exp(pars1);
      ee4 = log(y);
      
      out(j, 0) += ee1 - y/ee3;
      out(j, 1) += (ee2 + pars1 - ee4) * ee1;
      out(j, 2) += 1 * y/ee3;
      out(j, 3) += ee1;
      out(j, 4) += (ee2 + ee1 * R::trigamma(ee1) + pars1 - ee4) * ee1;
      
    }

}

return out;

}

// //' @rdname gammad0
// [[Rcpp::export]]
arma::mat gammad34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
    
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, pars1, pars2;
  
  double ee1, ee3, ee4, ee5, ee6, ee7;
  
  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];

    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);

      ee1 = exp(pars2);
      ee3 = 1 * y/exp(pars1);
      ee4 = R::digamma(ee1);
      ee5 = log(y);
      ee6 = R::psigamma(ee1, 2);
      ee7 = R::trigamma(ee1);
      
      out(j, 0) += -ee3;
      out(j, 1) += 0;
      out(j, 2) += ee1;
      out(j, 3) += ((3 * ee7 + ee1 * ee6) * ee1 + ee4 + pars1 - ee5) *
        ee1;
      out(j, 4) += ee3;
      out(j, 5) += 0;
      out(j, 6) += 0;
      out(j, 7) += ee1;
      out(j, 8) += (((6 * ee6 + ee1 * R::psigamma(ee1, 3)) * ee1 + 7 *
        ee7) * ee1 + ee4 + pars1 - ee5) * ee1;
      
    }

}

return out;

}
