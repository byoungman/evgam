// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// //' GW negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each GW parameter
// //' @param X1 a design matrix for the log scale parameter
// //' @param X2 a design matrix for the shape parameter
// //' @param ymat a matrix
// //' @param leftmat a matrix
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return lgw0 a scalar, the negative log-likelihood
// //' @return lgw12 a matrix, first then second derivatives w.r.t. parameters
// //' @return lgw34 a matrix, third then fourth derivatives w.r.t. parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double gwd0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, 
            arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat wmat, arma::mat pmat)
{

  int nobs = nhere.size();

  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);

  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }

  double y, w, pars1, pars2, scale, theta, prob, base, part1, part2, part3;

  double nllh = 0.0;

  for (int j=0; j < nobs; j++) {

    pars1 = p1vec[j];
    pars2 = p2vec[j];
    scale = exp(pars1);
    theta = pars2;
    
    for (int l=0; l < nhere[j]; l++) {

      y = ymat(j, l);
      prob = pmat(j, l);
      w = wmat(j, l);
      
      base = 1 + (theta * y) / scale;
      part1 = log(-log(prob)) - log(scale);
      part2 = (1/theta - 1) * log(base);
      part3 = log(prob) * R_pow(base, 1/theta);

      nllh -= w * (part1 + part2 + part3);

    }

  }

return(nllh);

}

// //' @rdname gwd0
// [[Rcpp::export]]
arma::mat gwd12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, 
               arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat wmat, arma::mat pmat)
{
  
  int nobs = nhere.size();
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, w, pars1, pars2, prob;
  
  arma::mat out = arma::mat(nobs, 5, arma::fill::zeros);
  
  double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee15, ee16;
  double ee20, ee21, ee22, ee24, ee26;
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];

    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      prob = pmat(j, l);
      w = wmat(j, l);
      
      ee1 = exp(pars1);
      ee2 = pars2 * y;
      ee3 = ee2/ee1;
      ee4 = 1 + ee3;
      ee5 = 1/pars2;
      ee6 = ee5 - 1;
      ee7 = ee4 * ee1;
      ee8 = log1p(ee3);
      ee9 = R_pow(ee4, ee6);
      ee10 = log(prob);
      ee11 = ee8/pars2;
      ee12 = R_pow(ee4, ee5);
      ee13 = y * ee9;
      ee15 = y * ee6/ee7;
      ee16 = R_pow(ee4, (ee5 - 2));
      ee20 = pars2 * ee6;
      ee21 = ee2/ee7;
      ee22 = R_pow(pars2, 2);
      ee24 = ee13/ee1 - ee12 * ee8/pars2;
      ee26 = y/ee7 - 2 * ee11;
      
      out(j, 0) += w * (1 + y * (ee9 * ee10 + ee20/ee4)/ee1);
      out(j, 1) += w * (-((ee10 * ee24 - ee11)/pars2 + ee15));
      out(j, 2) += w * (y * (ee20 * (ee21 - 1)/ee4 - (ee9 + ee2 *
        ee16 * ee6/ee1) * ee10)/ee1);
      out(j, 3) += w * (y * (((1 - ee21) * ee6 - ee5)/ee4 + ee10 *
        (y * ee16 * ee6/ee1 - ee9 * ee8/ee22))/ee1);
      out(j, 4) += w * (-((ee10 * (ee13 * (ee15 - (1 + ee11)/pars2)/
        ee1 - (ee12 * ee26 + ee8 * ee24/pars2)/pars2) - ee26/pars2)/
          pars2 - y * (1/ee22 + ee15)/ee7));  
      
    }
    
  }
  
  return out;
  
}

// //' @rdname gwd0
// [[Rcpp::export]]
arma::mat gwd34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, 
                arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat wmat, arma::mat pmat)
{
  
  int nobs = nhere.size();
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, w, pars1, pars2, prob;
  
  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  
  double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee22, ee24, ee25, ee28;
  double ee31, ee32, ee33, ee34, ee35, ee36, ee37, ee38, ee39;
  double ee40, ee41, ee42, ee43, ee48;
  double ee50, ee52, ee54, ee55, ee56, ee57;
  double ee61, ee65;
  double ee73, ee74, ee75, ee76, ee77;
  double ee80, ee83, ee84, ee87, ee89;
  double ee90, ee96, ee98;
  double ee100, ee101, ee102, ee103, ee105, ee106, ee108;
  double ee112;
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      prob = pmat(j, l);
      w = wmat(j, l);
      
      ee1 = exp(pars1);
      ee2 = pars2 * y;
      ee3 = ee2/ee1;
      ee4 = 1 + ee3;
      ee5 = 1/pars2;
      ee6 = ee4 * ee1;
      ee7 = log1p(ee3);
      ee8 = ee5 - 1;
      ee9 = ee5 - 2;
      ee10 = R_pow(pars2, 2);
      ee11 = y/ee6;
      ee12 = R_pow(ee4, ee9);
      ee13 = ee7/pars2;
      ee14 = R_pow(ee4, ee8);
      ee15 = ee5 - 3;
      ee16 = R_pow(ee4, ee15);
      ee17 = 2 * ee13;
      ee18 = ee11 - ee17;
      ee19 = 2 * ee11;
      ee20 = 2/pars2;
      ee22 = ee14 * ee7/ee10;
      ee24 = ee12 * ee7/ee10;
      ee25 = R_pow(ee4, ee5);
      ee28 = y * ee12 * ee8/ee1;
      ee31 = y * ee16 * ee9/ee1;
      ee32 = ee31 - ee24;
      ee33 = ee2/ee6;
      ee34 = ee28 - ee22;
      ee35 = log(prob);
      ee36 = y * ee14;
      ee37 = ee7 * ee34;
      ee38 = y * ee8;
      ee39 = ee12/ee10;
      ee40 = ee8 * ee32;
      ee41 = 2 * ee33;
      ee42 = ee20 + ee11;
      ee43 = ee38/ee6;
      ee48 = (ee19 - 6 * ee13)/pars2 + y * ee42/ee6;
      ee50 = y * (ee40 - ee39)/ee1;
      ee52 = ee36/ee1 - ee25 * ee7/pars2;
      ee54 = ee14 * ee18 + ee37;
      ee55 = ee5 - 4;
      ee56 = ee19 + 6/pars2;
      ee57 = ee54/ee10;
      ee61 = ee19 + ee20;
      ee65 = y * R_pow(ee4, ee55) * ee15/ee1 - ee16 * ee7/ee10;
      ee73 = ee56/pars2 + y * ee61/ee6;
      ee74 = 2 - ee41;
      ee75 = 2 + ee2 * ee9/ee6;
      ee76 = 4 * ee4;
      ee77 = ee50 - ee57;
      ee80 = y * (ee9 * ee65 - ee16/ee10)/ee1 - (ee12 * ee18 +  ee7 * ee32)/ee10;
      ee83 = ee36 * (ee43 - (1 + ee13)/pars2)/ee1 - (ee25 * ee18 +  ee7 * ee52/pars2)/pars2;
      ee84 = y * ee9;
      ee87 = ee48 * ee14;
      ee89 = ((24 * ee13 - 6 * ee11)/pars2 - y * ee56/ee6)/pars2 -  y * ee73/ee6;
      ee90 = (1 + ee2 * (ee8 * ee75 + ee5)/ee1) * ee12;
      ee96 = 1 + ee2 * (ee41 - 3)/ee6;
      ee98 = 2 * (1 + 2 * ee3) + ee76;
      ee100 = 6 * ee3;
      ee101 = 8 * ee3;
      ee102 = ee7 * ee77;
      ee103 = pars2 * ee8;
      ee105 = ee2 * (3 - ee41)/ee6;
      ee106 = ee33 - 1;
      ee108 = y * (ee8 * ee80 - ee12 * (2 * (ee84/ee6) - (2 +  ee17)/pars2)/ee10)/ee1;
      ee112 = y * (ee50 - ((ee14 * (ee11 - 3 * ee13) + ee37)/pars2 +  ee14 * (ee43 - (2 + ee13)/pars2) + ee28)/pars2)/ee1 -  ((2 * (ee52 * ee18) + ee7 * ee83)/pars2 - ee48 * ee25)/pars2;
      
      out(j, 0) += w * (y * (ee90 * ee35 + pars2 * ee96 * ee8/ee4)/ee1);
      out(j, 1) += w * (y * ((ee8 * (ee105 - 1) - ee106/pars2)/ee4 - ee35 *
        (y * (ee12 * ee9 + ee103 * ee32)/ee1 - ee22))/ee1);
      out(j, 2) += w * (y * (ee35 * ee77 - 2 * ((ee8 * ee74 - ee20)/R_pow(ee4, 2)))/ee1);
      out(j, 3) += w * (-((ee48/pars2 + ee35 * ee112)/pars2 + y * (ee42/
        ee10 + y * (1/ee10 + 2 * ee43)/ee6)/ee6));
      out(j, 4) += w * (y * (ee103 * (ee2 * (7 + ee2 * ((ee101 - ee98)/
        ee4 - 6)/ee6)/ee6 - 1)/ee4 - (ee90 + ee2 * (ee12 * ee75 + ee16 *
          (2 + ee2 * (ee9 * (2 + ee2 * ee15/ee6) + ee20 - 2)/ee1)) *
          ee8/ee1) * ee35)/ee1);
      out(j, 5) += w * (y * (((1 + ee2 * (ee2 * ((ee98 - ee101)/ee4 + 6)/
        ee6 - 7)/ee6) * ee8 - ee96/pars2)/ee4 + ee35 * (y * (ee12 *
          ee55 + pars2 * (3 * ee40 + y * (ee8 *   (pars2 * ee9 * ee65 -
          2 * ee16) - ee16 * ee9)/ee1))/ee1 - ee22))/ee1);
      out(j, 6) += w * (y * ((ee38 * (4 - ee2 * ((ee76 - ee100)/ee4 + 6)/
        ee6)/ee6 - (2 * ee105 - (2 + 2 * ee106))/ee10)/ee4 - ee35 *
          (y * (ee16 * (2 * (ee4 * ee7/ee10) - 2 * (ee84/ee1)) + ee8 *
          (pars2 * ee80 + ee31 - ee24) - ee39)/ee1 - ee57))/ee1);
      out(j, 7) += w * (y * ((((6 * (1 - ee33) - 6)/pars2 + 2 * (y * ee74/
        ee6))/ee10 + y * (ee74/ee10 + y * ((2 * ee4 - ee100)/ee4 +
          4) * ee8/ee6)/ee6)/ee4 + ee35 * (ee108 - (2 * (ee34 * ee18) +
          ee102 - ee87)/ee10))/ee1);
      out(j, 8) += w * (-((ee89/pars2 + ee35 * (y * (ee108 - (((2 * ee18 -
        6) * ee34 + ee102 - ((ee14 * (ee11 - (ee17 + 6)) + 2 * ee54 +
        ee37)/pars2 +   ee87))/pars2 + 3 * ee50)/pars2)/ee1 -
        ((3 * (ee83 * ee18) + ee7 * ee112 - 3 * (ee48 * ee52))/pars2 -
        ee89 * ee25)/pars2))/pars2 - y * (ee73/ee10 + y * (ee61/
          ee10 + y * (2/ee10 + 6 * ee43)/ee6)/ee6)/ee6));  
    }
    
  }
  
  return out;
  
}
