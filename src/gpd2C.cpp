// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// //' Generalized Pareto distribution (GPD) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each GPD parameter
// //' @param X1 a design matrix for the GPD log scale parameter
// //' @param X2 a design matrix for the GPD transformed shape parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return gpd2d0 a scalar, the negative log-likelihood
// //' @return gpd2d12 a matrix, first then second derivatives w.r.t. gpd2ull parameters
// //' @return gpd2d34 a matrix, third then fourth derivatives w.r.t. gpd2ull parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double gpd2d0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
    
arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = nhere.size();

if (dcate == 1) {
  p1vec = p1vec.elem(dupid);
  p2vec = p2vec.elem(dupid);
}

double nllh = 0.0;
double pars1, pars2, xi, ixi, y, ee1;

for (int j=0; j < nobs; j++) {
  
  pars1 = p1vec[j];
  pars2 = p2vec[j];
  xi = 1.5 / (1.0 + exp(-pars2)) - 1.0;
  ixi = 1.0 / xi;
  
  for (int l=0; l < nhere[j]; l++) {
    
    y = ymat(j, l);
    ee1 = xi * y / exp(pars1);
    
    if (ee1 <= -1.0) {
      nllh = 1e20;
      break;
    } 
    
    nllh += pars1 + (1.0 + ixi) * log1p(ee1);
  
  }
  
}

return(nllh);

}

// //' @rdname gpd2d0
// [[Rcpp::export]]
arma::mat gpd2d12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
    
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, pars1, pars2;
  
  double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee16, ee17, ee18, ee19;
  
  arma::mat out = arma::mat(nobs, 5, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];

    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      
      ee2 = exp(-pars2);
      ee3 = 1 + ee2;
      ee5 = 1.5/ee3 - 1;
      ee6 = exp(pars1);
      ee7 = y * ee5;
      ee8 = ee7/ee6;
      ee9 = 1 + ee8;
      ee10 = ee9 * ee6;
      ee11 = R_pow(ee3, 2);
      ee12 = 1 + 1/ee5;
      ee13 = R_pow(ee5, 2);
      ee14 = y * ee12;
      ee16 = ee11 * ee9 * ee6;
      ee17 = log1p(ee8);
      ee18 = ee14 * ee5;
      ee19 = ee7/ee10;
      
      out(j, 0) += 1 - ee18/ee10;
      out(j, 1) += (1.5 * (ee14/ee10) - 1.5 * (ee17/ee13)) * ee2/
        ee11;
      out(j, 2) += -(ee18 * (ee19 - 1)/ee10);
      out(j, 3) += -(y * (ee12 * (1.5 - 1.5 * ee19) - 1.5/ee5) * ee2/
        ee16);
      out(j, 4) += -(((2.25 * (y * ee2/ee16) - ((4.5/(ee3 * ee5) -
        3) * ee2/ee3 + 1.5) * ee17)/ee13 + y * (((2.25 * (y/(ee3 *
        ee9 * ee6)) - 3) * ee2/ee3 + 1.5) * ee12 + 2.25 * (ee2/(ee11 *
        ee13)))/ee10) * ee2/ee11);
      
    }

}

return out;

}

// //' @rdname gpd2d0
// [[Rcpp::export]]
arma::mat gpd2d34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
    
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, pars1, pars2;
  
  double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee21, ee22, ee25, ee26, ee27, ee28, ee29;
  double ee30, ee33, ee39;
  double ee40, ee41, ee42, ee43, ee44, ee47, ee48, ee49;
  double ee50, ee53, ee54, ee55, ee58, ee59;
  double ee61, ee64, ee65, ee66, ee67, ee69;
  double ee70, ee71, ee72, ee73, ee74, ee77, ee78, ee79;
  double ee80, ee81, ee84, ee85;
  
  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];

    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);

      ee2 = exp(-pars2);
      ee3 = 1 + ee2;
      ee5 = 1.5/ee3 - 1;
      ee6 = exp(pars1);
      ee7 = y * ee5;
      ee8 = ee7/ee6;
      ee9 = 1 + ee8;
      ee10 = R_pow(ee3, 2);
      ee11 = 1.5 - 3 * (ee2/ee3);
      ee12 = ee9 * ee6;
      ee13 = ee2/ee10;
      ee14 = 1 + 2 * ee2;
      ee15 = 3 * ee14;
      ee16 = 3 * ee3;
      ee17 = ee11 * ee5;
      ee18 = 12 * ee2;
      ee19 = 2.25 * ee13;
      ee21 = ee10 * ee9 * ee6;
      ee22 = ee15 + ee16;
      ee25 = y/(ee3 * ee9 * ee6);
      ee26 = ee19 - ee17;
      ee27 = y * ee2;
      ee28 = 1 + 1/ee5;
      ee29 = ee7/ee12;
      ee30 = ee3 * ee5;
      ee33 = (2.25 * ee25 - 3) * ee2/ee3 + 1.5;
      ee39 = (4.5/ee30 - 3) * ee2/ee3 + 1.5;
      ee40 = 1.5 - ((ee22 - ee18)/ee3 + 3) * ee2/ee3;
      ee41 = 3 * ee26;
      ee42 = ee27/(ee10 * ee6);
      ee43 = 3 * ee17;
      ee44 = 4.5 * ee11;
      ee47 = ((ee18 + y * (ee44 + 6.75 * (ee27/ee21))/ee12 - ee22)/ee3 -  3) * ee2/ee3 + 1.5;
      ee48 = ee9 * ee11;
      ee49 = R_pow(ee5, 2);
      ee50 = (ee41 - (27 * ee13 + ee43))/ee5;
      ee53 = (4.5 * ee25 - 3) * ee2/ee3 + 1.5;
      ee54 = 3 * ee11;
      ee55 = 4 * ee2;
      ee58 = (((ee50 - ee54)/ee5 + ee15 + ee16 - ee18)/ee3 + 3) *  ee2/ee3 - 1.5;
      ee59 = ee40 * ee5;
      ee61 = 2.25 * ee42 - ee48;
      ee64 = ((6 * (2 * ee14 + ee55) + 8 * ee22 - 96 * ee2)/ee3 +  12) * ee2;
      ee65 = ee10 * ee5;
      ee66 = ee10 * ee49;
      ee67 = ee33 + y * (ee19 - ee53 * ee5)/ee12;
      ee69 = 1 + ee7 * (2 * ee29 - 3)/ee12;
      ee70 = 1.5 - 1.5 * ee29;
      ee71 = 2 * ee26;
      ee72 = 3 * (ee33 * ee39);
      ee73 = 3 * (1 + ee55);
      ee74 = 3 * ee40;
      ee77 = 6 * ee3;
      ee78 = 9 * ee14;
      ee79 = 9 * ee13;
      ee80 = log1p(ee8);
      ee81 = y * ee28;
      ee84 = ee7 * (4.5 - 3 * ee29)/ee12 - 1.5;
      ee85 = ee29 - 1;
      
      out(j, 0) += -(ee81 * ee69 * ee5/ee12);
      out(j, 1) += -(y * (ee28 * ee84 - 1.5 * (ee85/ee5)) * ee2/ee21);
      out(j, 2) += y * ((((3 * ee70 - 4.5)/ee30 + 3) * ee2/ee3 - 1.5)/
        ee5 + ee67 * ee28) * ee2/ee21;
      out(j, 3) += ((ee58 * ee80 + y * (1.5 * ee33 + 3 * ee39) * ee2/
        ee21)/ee49 + y * (ee47 * ee28 + (1.5 * ee39 + 3 * ee33) *
          ee2/ee66)/ee12) * ee2/ee10;
      out(j, 4) += -(ee81 * ee5 * (ee7 * (7 + y * ((8 * ee8 - (2 *
        (1 + 2 * ee8) + 4 * ee9))/ee9 - 6) * ee5/ee12)/ee12 - 1)/ee12);
      out(j, 5) += -(y * (ee28 * (1.5 + ee7 * (y * ((2 * (1.5 * ee9 +
        1.5 * ee8) + 6 * ee9 - 12 * ee8)/ee9 + 9) * ee5/ee12 - 10.5)/
          ee12) - 1.5 * (ee69/ee5)) * ee2/ee21);
      out(j, 6) += -(y * ((ee33 + y * (ee71 + ee19 - (ee53 + y * ((ee5 *
        (2 * ee61 - 18 * ee42) + 9 * (ee9 * ee2/ee10))/ee9 +
        ee79)/ee12) * ee5)/ee12) * ee28 - (3 * (ee2 * ee84/ee65) - ee39 *
        ee85)/ee5) * ee2/ee21);
      out(j, 7) += -(y * (((((ee50 + 3 * (ee39 * ee70) + 4.5 * ee67 -
        ee54)/ee5 + ee15 + ee16 - ee18)/ee3 + 3) * ee2/ee3 - 1.5)/
          ee5 + (ee47 + y * ((ee44 + y * ((ee5 * (3 * ee61 - 27 * ee42) +
            3 * (ee9 * ee26))/ee9 + ee41)/ee12) * ee2/ee10 - ee59)/
              ee12) * ee28) * ee2/ee21);
      out(j, 8) += ((y * (4.5 * ee58 - (1.5 * ee47 + ee72)) * ee2/
        ee21 - (((((ee11 * (6 * ee26 - 18 * ee13) + (12 * (ee41 - ee43) +
          9 * (ee71 + ee79) - 324 * ee13) * ee2/ee65 + 3 * (4.5 *
          (ee11 * ee2/ee10) - ee59) - 6 * ee59)/ee5 - ee74)/ee5 + ee73 +
          ee77 + ee78 - ee64)/ee3 + 3) * ee2/ee3 - 1.5) * ee80)/
            ee49 + y * ((((ee73 + ee77 + ee78 + y * (y * ((4.5 * ee61 - (40.5 *
              ee42 + 9 * ee48))/ee9 - 9 * ee11) * ee2/ee21 - (ee53 *
              ee11 + 2 * (R_pow(ee11, 2) + 1.5 * ee40) + ee74))/ee12 - ee64)/
                ee3 + 3) * ee2/ee3 - 1.5) * ee28 + 
                  (1.5 * ee58 - (ee72 +
                  4.5 * ee47)) * ee2/ee66)/ee12) * ee2/ee10;
      
    }

}

return out;

}
