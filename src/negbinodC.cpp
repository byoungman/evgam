// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// //' Negative binomial distribution negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for neg bin log mean and log overdispersion
// //' @param X1 a design matrix for the neg bin log mean parameter
// //' @param X2 a design matrix for the neg bin log overdispersion mean parameter
// //' @param ymat a matrix
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return negbinodd0 a scalar, the negative log-likelihood
// //' @return negbinodd12 a matrix, first then second derivatives w.r.t. negbinson parameters
// //' @return negbinodd34 a matrix, third then fourth derivatives w.r.t. negbinson parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double negbinodd0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, 
                arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat wmat)
{
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, w, pars1, pars2;
  double mu, size, p;
  double nllh = 0.0;
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      w = wmat(j, l);

      mu = exp(pars1);
      // size = exp(pars2);
      size = exp(-pars2);
      p = size / (size + mu);
      nllh -= lgamma(y + size) - lgamma(size) - lgamma(y + 1) + size * log(p) + y * log(1 - p);

    } 
    
  }
  
  return(nllh);
  
}

// //' @rdname negbinodd0
// [[Rcpp::export]]
arma::mat negbinodd12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, 
                    arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat wmat)
{
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, w, pars1, pars2;
  // double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  // double ee10, ee11;
  double ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12;

  arma::mat out = arma::mat(nobs, 5, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];

    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      w = wmat(j, l);
      
      // ee1 = exp(pars2);
      // ee2 = exp(pars1);
      // ee3 = ee2 + ee1;
      // ee4 = ee1/ee3;
      // ee5 = 1 - ee4;
      // ee6 = ee1 + y;
      // ee7 = R_pow(ee4, 2);
      // ee8 = ee5 * ee3;
      // ee9 = Rf_digamma(ee6);
      // ee10 = Rf_digamma(ee1);
      // ee11 = log(ee3);
      // 
      // out(j, 0) += w * (-(ee2 * ee1 * (y/ee8 - 1)/ee3));
      // out(j, 1) += w * (-((1 + ee9 + pars2 - (ee6/ee3 + ee10 + ee11)) *
      //   ee1));
      // out(j, 2) += w * (-(((y * (1 - (2 + ee1/ee8) * ee2/ee3)/ee5 -
      //   ee2 * (R_pow(ee1, 2)/(R_pow(ee3, 2) * ee7) - 2))/ee3 - 1) *
      //   ee2 * ee1/ee3));
      // out(j, 3) += w * (-((((ee5 * ee1/(ee3 * ee7) + 2) * ee1 + y)/
      //   ee3 - 2) * ee2 * ee1/ee3));
      // out(j, 4) += w * (-((3 + ee9 + ee1 * (Rf_trigamma(ee6) - (((R_pow(ee5, 2)/
      //   ee7 - 2) * ee1/ee3 + 5)/ee3 + Rf_trigamma(ee1))) + pars2 -
      //     (ee10 + ee11 + y * ((1 - (3 - 2 * ee4) * ee1/ee3)/ee5 +
      //     ee4)/ee3)) * ee1));

ee2 = exp(-pars2);
ee3 = exp(pars1);
ee4 = ee2 + ee3;
ee5 = ee2/ee4;
ee6 = 1 - ee5;
ee7 = ee2 + y;
ee8 = R_pow(ee5, 2);
ee9 = ee6 * ee4;
ee10 = Rf_digamma(ee7);
ee11 = Rf_digamma(ee2);
ee12 = log(ee4);

out(j, 0) += w * (-(ee2 * ee3 * (y/ee9 - 1)/ee4));
out(j, 1) += w * (-((ee7/ee4 + ee11 + ee12 + pars2 - (1 + ee10)) *
  ee2));
out(j, 2) += w * (-(((y * (1 - (2 + ee2/ee9) * ee3/ee4)/ee6 -
  (R_pow(ee2, 2)/(R_pow(ee4, 2) * ee8) - 2) * ee3)/ee4 - 1) *
  ee2 * ee3/ee4));
out(j, 3) += w * ((((ee6 * ee2/(ee4 * ee8) + 2) * ee2 + y)/
  ee4 - 2) * ee2 * ee3/ee4);
out(j, 4) += w * ((ee11 + ee12 + pars2 + y * ((1 - (3 - 2 *
  ee5) * ee2/ee4)/ee6 + ee5)/ee4 - (3 + ee10 + ee2 * (Rf_trigamma(ee7) -
  (((R_pow(ee6, 2)/ee8 - 2) * ee2/ee4 + 5)/ee4 + Rf_trigamma(ee2))))) *
  ee2);

    }
    
  }
  
  return out;
  
}

// //' @rdname negbinodd0
// [[Rcpp::export]]
arma::mat negbinodd34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, 
                    arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat wmat)
{
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, w, pars1, pars2;
  // double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  // double ee10, ee11, ee13, ee14, ee15, ee16, ee17, ee19;
  // double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
  // double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee38;
  // double ee40, ee41, ee43, ee44, ee46, ee48, ee49;
  // double ee51, ee53, ee54, ee55, ee57, ee58;
  // double ee61, ee64, ee65, ee66, ee67, ee68, ee69;
  // double ee71, ee77, ee78, ee79;
  // double ee80, ee81, ee82, ee83, ee84, ee86, ee87, ee88;
  // double ee100, ee103, ee105, ee106, ee107, ee108, ee109;
  // double ee110, ee111, ee112, ee113, ee114, ee115, ee116, ee117, ee119;
  // double ee120, ee121, ee122, ee123, ee124, ee126, ee128, ee129;
  // double ee130, ee131, ee132;

double ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee14, ee15, ee16, ee17, ee19;
double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee38, ee39;
double ee40, ee42, ee44, ee45, ee47, ee48, ee49;
double ee51, ee52, ee54, ee55, ee56, ee59;
double ee62, ee63, ee64, ee65, ee66, ee67, ee68;
double ee70, ee72, ee77;
double ee81, ee82, ee83, ee84, ee85, ee86, ee87, ee89;
double ee90, ee91;
double ee103, ee104, ee106, ee107, ee108, ee109;
double ee110, ee111, ee112, ee113, ee114, ee115, ee116, ee117, ee118, ee119;
double ee121, ee122, ee123, ee124, ee125, ee126, ee129;
double ee131, ee132, ee133, ee134, ee135;

  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      w = wmat(j, l);

      // ee1 = exp(pars2);
      // ee2 = exp(pars1);
      // ee3 = ee2 + ee1;
      // ee4 = ee1/ee3;
      // ee5 = 1 - ee4;
      // ee6 = 2 * ee4;
      // ee7 = 2 * ee1;
      // ee8 = 2 * ee3;
      // ee9 = ee7 + ee2;
      // ee10 = R_pow(ee4, 2);
      // ee11 = 2 * ee9;
      // ee13 = (3 - ee6) * ee1/ee3;
      // ee14 = ee2/ee3;
      // ee15 = 2 * ee14;
      // ee16 = 1 - ee15;
      // ee17 = 1 - ee6;
      // ee19 = 2 * ee2 + ee1;
      // ee20 = 2 * ee19;
      // ee21 = 4 * ee3;
      // ee22 = 1 - ee13;
      // ee23 = 8 * ee1;
      // ee24 = R_pow(ee1, 2);
      // ee25 = R_pow(ee3, 2);
      // ee26 = ee20 + ee8;
      // ee27 = ee11 + ee21;
      // ee28 = 4 * ee1;
      // ee29 = ee3 * ee10;
      // ee30 = ee25 * ee10;
      // ee31 = ee11 + ee8;
      // ee32 = 8 * ee2;
      // ee33 = R_pow(ee5, 2);
      // ee34 = (ee27 - ee23)/ee3;
      // ee35 = 6 * ee1;
      // ee36 = ee1 + y;
      // ee38 = (ee26 - ee32)/ee3 + 2;
      // ee40 = ee8 - ee35;
      // ee41 = 4 * (ee11 + ee28);
      // ee43 = (ee34 + 6) * ee1/ee3;
      // ee44 = ee17 * ee5;
      // ee46 = ee5 * ee24/ee30;
      // ee48 = ee40 * ee2/ee3;
      // ee49 = 3 * ee4;
      // ee51 = ee38 * ee2/ee3;
      // ee53 = ((ee31 - ee23)/ee3 + 2) * ee1/ee3;
      // ee54 = (ee48 + ee7)/ee3;
      // ee55 = ee16 * ee5;
      // ee57 = ee44 + 1 - ee13;
      // ee58 = ee5 * ee3;
      // ee61 = (2 * (ee33/ee10) - 2) * ee1/ee3;
      // ee64 = 1 - (7 - ee43) * ee1/ee3;
      // ee65 = 2 * (ee28 + ee2);
      // ee66 = ee8 + ee7;
      // ee67 = 4 * ee17;
      // ee68 = 4 * ee2;
      // ee69 = 8 * ee46;
      // ee71 = R_pow(ee22, 2) + ee64 * ee5;
      // ee77 = 2 * (ee33 + 1 - ee13);
      // ee78 = 2 * ee22;
      // ee79 = 2 * (1 - (4 - ee49) * ee1/ee3);
      // ee80 = 2 * ee16;
      // ee81 = 2 * ee17;
      // ee82 = 64 * ee1;
      // ee83 = 8 * ee4;
      // ee84 = ee71/ee5;
      // ee86 = (ee61 + 3) * ee1/ee3;
      // ee87 = ((4 * (ee20 + ee68) + 8 * ee26 - 64 * ee2)/ee3 +  8) * ee2;
      // ee88 = (1 - ee54) * ee5;
      // ee100 = ee34 + 12;
      // ee103 = (2 * (ee24/ee30) - 2) * ee2/ee3 + 1;
      // ee105 = (2 + 2 * (ee1/ee58)) * ee2/ee3;
      // ee106 = 1 - ee53;
      // ee107 = 2 * (ee16 * ee17 + 1 - ee54);
      // ee108 = 2 * (R_pow(ee16, 2) + 1 - ee51);
      // ee109 = 2 * ee57;
      // ee110 = 2 * (R_pow(ee17, 2) + 1 - ee53);
      // ee111 = 2 * (ee5 * ee1/ee29);
      // ee112 = 2 * ee38;
      // ee113 = 2 * (ee68 + ee1);
      // ee114 = 4 * ee55;
      // ee115 = 4 * ee16;
      // ee116 = ee41 + 8 * ee31;
      // ee117 = ee41 + 8 * ee27;
      // ee119 = 4 * ee66 - 40 * ee1;
      // ee120 = 6 * ee19;
      // ee121 = 6 * ee3;
      // ee122 = 8 * ee22;
      // ee123 = 8 * ee66;
      // ee124 = Rf_digamma(ee36);
      // ee126 = Rf_digamma(ee1) + log(ee3);
      // ee128 = ee2 * ee1/ee25;
      // ee129 = Rf_psigamma(ee36, 2);
      // ee130 = Rf_psigamma(ee1, 2);
      // ee131 = Rf_trigamma(ee36);
      // ee132 = Rf_trigamma(ee1);
      // 
      // out(j, 0) += -(((y * (1 - (((1 + ee80 - ee105) * ee1/ee5 + ee20 +
      //   ee8 - ee32)/ee3 + 2) * ee2/ee3)/ee5 - (((ee103 + ee80) *
      //   ee24/ee29 + ee32 - ee26)/ee3 - 2) * ee2)/ee3 - 1) * ee2 *
      //   ee1/ee3);
      // out(j, 1) += -((((ee55 * ee1/ee29 + 2) * ee1 + y * (1 - ((((ee81 +
      //   ee6)/ee5 - 6) * ee1 + ee8) * ee2/ee3 + (1 + ee15) * ee1)/
      //     ee3)/ee5 - ((((1 + ee81 - 2 * ee46) * ee1/ee29 + 6) * ee1 -
      //       ee8)/ee3 - 2) * ee2)/ee3 - 2) * ee2 * ee1/ee3);
      // out(j, 2) += -(((((((ee5 * (2 - (ee111 + 4) * ee1/ee3) + 3 -
      //   (5 - ee6) * ee1/ee3)/ee10 - 8) * ee1 + ee11 + ee8)/ee3 + 6) *
      //   ee1 + y * ((ee57/ee5 + (ee23 - ee31)/ee3 - 1) * ee1/ee3 +
      //   1)/ee5)/ee3 - 4) * ee2 * ee1/ee3);
      // out(j, 3) += -(((3 * ee131 + ee1 * (ee129 - ee130) - (((ee5 *
      //   (ee78 + 4 - (ee61 + 6) * ee1/ee3)/ee10 - ee100) * ee1/ee3 +
      //   19)/ee3 + 3 * ee132)) * ee1 + 7 + ee124 + pars2 - (ee126 +
      //   y * (((((ee27 - ee35)/ee3 + 3) * ee1/ee3 - 6) * ee1/ee3 + 1)/
      //     ee5 + (2 * (ee22/ee5) + ee6) * ee1/ee3)/ee3)) * ee1);
      // out(j, 4) += -(((y * (1 - ((((1 - ee105) * ee16 + 2 + ee108 -
      //   (((2 * (ee55 + ee128) + ee114 - 8 * ee128)/ee5 + ee115) *
      //   ee1/ee58 + ee112) * ee2/ee3) * ee1/ee5 + ee113 + ee21 + ee120 -
      //   ee87)/ee3 + 2) * ee2/ee3)/ee5 - (((((2 * (1 - 3 * ee14) +
      //   8 * ee16 + 8 * (ee2 * ee24/(R_pow(ee3, 3) * ee10))) * ee24/
      //     ee30 - ee112) * ee2/ee3 + ee103 * ee16 + 2 + ee108) * ee24/
      //       ee29 + ee87 - (ee113 + ee21 + ee120))/ee3 - 2) * ee2)/ee3 -
      //         1) * ee2 * ee1/ee3);
      // out(j, 5) += -(((((1 - ee51) * ee5 * ee1/ee29 + 2) * ee1 + y *
      //   (1 - ((((((ee115 - ((2 * (1 - ee49) + ee67 + ee83) * ee2/
      //     ee58 + 2)) * ee1 - (ee40/ee3 + 2) * ee2)/ee3 + 2 + ee107)/ee5 -
      //       (2 + 4 * (ee26/ee3))) * ee1 + ee20 + ee21 - ee119 * ee2/
      //         ee3) * ee2/ee3 + (ee51 + 1) * ee1)/ee3)/ee5 - (((((((((2 * (2 -
      //           ee49) + ee67 - ee69) * ee2/ee3 - ee114) * ee1/ee29 - 2) *
      //           ee1 + (((2 * (ee1/ee29) + 6) * ee1 - ee8)/ee3 - 2) * ee2)/
      //             ee3 + ee16 * (3 - (2 + ee111) * ee1/ee3) + 2 + ee107) * ee1/
      //               ee10 + 4 * ee26)/ee3 + 2) * ee1 + (ee119/ee3 + 8) * ee2 - (4 *
      //                 ee19 + ee121))/ee3 - 2) * ee2)/ee3 - 2) * ee2 * ee1/ee3);
      // out(j, 6) += -(((((((ee16 * (3 - (ee61 + 5) * ee1/ee3) + 2 *
      //   ee88)/ee10 - 8) * ee1 + ee11 + ee8)/ee3 + 6) * ee1 + y * ((((ee88 +
      //   ee22 * ee16)/ee5 + ((ee80 + 6) * ee1 - (ee48 + ee11 +
      //   ee8))/ee3 - 1) * ee1 - (((((ee79 + 4 * ee44)/ee5 + ee67 +
      //   ee83) * ee1/ee3 + ee110)/ee5 + 6 - (ee41 + ee123 - ee82)/ee3) *
      //   ee1 + ee8) * ee2/ee3)/ee3 + 1)/ee5 - ((((((1 + ee110 +
      //   ee67 - (ee5 * (4 + 8 * ee17 - ee69) + ee77) * ee24/ee30)/ee10 -
      //   64) * ee1 + ee41 + ee123)/ee3 + 6) * ee1 - ee121)/ee3 -
      //   2) * ee2)/ee3 - 4) * ee2 * ee1/ee3);
      // out(j, 7) += -(((((((ee17 * (1 + ee78 - ee86) + ee5 * (2 * ee106 +
      //   7 - ((((ee5 * (6 - ee69) + ee77 + 4 * ee57)/ee10 - 8) *
      //   ee1 + ee11 + ee8)/ee3 + 14) * ee1/ee3) + 7 - (19 - ee100 *
      //   ee1/ee3) * ee1/ee3)/ee10 - ((ee116 - ee82)/ee3 + 32)) * ee1 +
      //   10 * ee3 + 12 * ee9 + ee65)/ee3 + 14) * ee1 + y * ((((ee43 +
      //   ee109 - 7) * ee1/ee3 + 1 + 2 * (ee106 * ee5) + 3 * (ee22 *
      //   ee17))/ee5 + (((ee109 + ee79)/ee5 + (ee116 - (ee31 + 48 *
      //   ee1))/ee3 + ee81 + 6) * ee1 - (ee65 + ee21 + 6 * ee9))/ee3 -
      //   1) * ee1/ee3 + 1)/ee5)/ee3 - 8) * ee2 * ee1/ee3);
      // out(j, 8) += -((((6 * ee129 + ee1 * (Rf_psigamma(ee36, 3) - Rf_psigamma(ee1, 3)) -
      //   6 * ee130) * ee1 + 7 * ee131 - (((((1 - ee86) *
      //   ee22 + ee5 * (10 + 2 * ee64 + ee122 - (((ee5 * (8 - ee69) +
      //   ee77 + ee122) * ee5/ee10 -   8) * ee1/ee3 + 18) * ee1/
      //     ee3) + 2 * ee71)/ee10 - ((ee65 + 22 * ee9 + 36 * ee3 - ((ee117 -
      //       ee82)/ee3 +   56) * ee1)/ee3 + 50)) * ee1/ee3 + 65)/ee3 +
      //       7 * ee132)) * ee1 + 15 + ee124 + pars2 - (ee126 + y * (((((((ee27 +
      //       56 * ee1 - ee117)/ee3 - 18) * ee1 + 14 * ee9 + ee65 +
      //       20 * ee3)/ee3 + ee78 + 7) * ee1/ee3 + ee84 - 14) * ee1/
      //         ee3 + 1)/ee5 + (((ee79 + 4 * ee22)/ee5 + ee83) * ee1/ee3 + (2 *
      //           ee84 + 4 * (ee22 * ee1/ee3))/ee5) * ee1/ee3)/ee3)) * ee1);      

      ee2 = exp(-pars2);
      ee3 = exp(pars1);
      ee4 = ee2 + ee3;
      ee5 = ee2/ee4;
      ee6 = 1 - ee5;
      ee7 = 2 * ee5;
      ee8 = 2 * ee2;
      ee9 = 2 * ee4;
      ee10 = ee8 + ee3;
      ee11 = R_pow(ee5, 2);
      ee12 = 2 * ee10;
      ee14 = (3 - ee7) * ee2/ee4;
      ee15 = ee3/ee4;
      ee16 = 1 - 2 * ee15;
      ee17 = 1 - ee7;
      ee19 = 2 * ee3 + ee2;
      ee20 = 2 * ee19;
      ee21 = 4 * ee4;
      ee22 = 1 - ee14;
      ee23 = R_pow(ee2, 2);
      ee24 = R_pow(ee4, 2);
      ee25 = 8 * ee2;
      ee26 = ee12 + ee21;
      ee27 = ee20 + ee9;
      ee28 = ee24 * ee11;
      ee29 = 4 * ee2;
      ee30 = ee12 + ee9;
      ee31 = 6 * ee2;
      ee32 = ee4 * ee11;
      ee33 = ee9 - ee31;
      ee34 = 8 * ee3;
      ee35 = R_pow(ee6, 2);
      ee36 = (ee26 - ee25)/ee4;
      ee38 = ee33 * ee3/ee4;
      ee39 = ee2 + y;
      ee40 = (ee27 - ee34)/ee4;
      ee42 = 4 * (ee12 + ee29);
      ee44 = (ee36 + 6) * ee2/ee4;
      ee45 = ee17 * ee6;
      ee47 = ee6 * ee23/ee28;
      ee48 = ee40 + 2;
      ee49 = 3 * ee5;
      ee51 = ((ee30 - ee25)/ee4 + 2) * ee2/ee4;
      ee52 = (ee38 + ee8)/ee4;
      ee54 = ee45 + 1 - ee14;
      ee55 = ee16 * ee6;
      ee56 = ee6 * ee4;
      ee59 = (2 * (ee35/ee11) - 2) * ee2/ee4;
      ee62 = 1 - (7 - ee44) * ee2/ee4;
      ee63 = 2 * ee16;
      ee64 = 2 * (ee29 + ee3);
      ee65 = ee9 + ee8;
      ee66 = 4 * ee17;
      ee67 = 4 * ee3;
      ee68 = 8 * ee47;
      ee70 = ee48 * ee3/ee4;
      ee72 = R_pow(ee22, 2) + ee62 * ee6;
      ee77 = (2 * (ee23/ee28) - 2) * ee3/ee4 + 1;
      ee81 = 2 * (ee35 + 1 - ee14);
      ee82 = 2 * ((4 - ee49) * ee2/ee4 - 1);
      ee83 = 2 * ee22;
      ee84 = 2 * ee17;
      ee85 = 64 * ee2;
      ee86 = 8 * ee5;
      ee87 = ee72/ee6;
      ee89 = (ee59 + 3) * ee2/ee4;
      ee90 = ((4 * (ee20 + ee67) + 8 * ee27 - 64 * ee3)/ee4 +  8) * ee3;
      ee91 = (1 - ee52) * ee6;
      ee103 = ee36 + 12;
      ee104 = ee77 + ee63;
      ee106 = (2 + 2 * (ee2/ee56)) * ee3/ee4;
      ee107 = 1 - ee51;
      ee108 = 2 * (ee17 * ee16 + 1 - ee52);
      ee109 = 2 * ee54;
      ee110 = 2 * (R_pow(ee17, 2) + 1 - ee51);
      ee111 = 2 * (R_pow(ee16, 2) + 1 - ee70);
      ee112 = 2 * (ee6 * ee2/ee32);
      ee113 = 2 * ee48;
      ee114 = 2 * (ee67 + ee2);
      ee115 = 4 * ee55;
      ee116 = ee66 + ee86;
      ee117 = 4 * ee16;
      ee118 = ee42 + 8 * ee30;
      ee119 = ee42 + 8 * ee26;
      ee121 = 4 * ee65 - 40 * ee2;
      ee122 = 6 * ee19;
      ee123 = 6 * ee4;
      ee124 = 8 * ee22;
      ee125 = 8 * ee65;
      ee126 = Rf_digamma(ee39);
      ee129 = Rf_digamma(ee2) + log(ee4) + pars2;
      ee131 = ee2 * ee3/ee24;
      ee132 = Rf_psigamma(ee39, 2);
      ee133 = Rf_psigamma(ee2, 2);
      ee134 = Rf_trigamma(ee39);
      ee135 = Rf_trigamma(ee2);

      out(j, 0) += -(((y * (1 - (((1 + ee63 - ee106) * ee2/ee6 + ee20 +
        ee9 - ee34)/ee4 + 2) * ee3/ee4)/ee6 - ((ee104 * ee23/
          ee32 + ee34 - ee27)/ee4 - 2) * ee3)/ee4 - 1) * ee2 * ee3/ee4);
      out(j, 1) += ((y * (1 - ((((ee84 + ee7)/ee6 + 2) * ee3/ee4 +
        1) * ee2 + ee38)/ee4)/ee6 - ((((ee84 - 2 * ee47) * ee3/ee4 -
        ee55) * ee2/ee32 - 2) * ee2 + (((6 + ee2/ee32) * ee2 - ee9)/
          ee4 - 2) * ee3))/ee4 - 2) * ee2 * ee3/ee4;
      out(j, 2) += -(((((((ee6 * (2 - (ee112 + 4) * ee2/ee4) + 3 -
        (5 - ee7) * ee2/ee4)/ee11 - 8) * ee2 + ee12 + ee9)/ee4 + 6) *
        ee2 + y * ((ee54/ee6 + (ee25 - ee30)/ee4 - 1) * ee2/ee4 +
        1)/ee6)/ee4 - 4) * ee2 * ee3/ee4);
      out(j, 3) += -((ee129 + y * (((((ee26 - ee31)/ee4 + 3) * ee2/
        ee4 - 6) * ee2/ee4 + 1)/ee6 + (2 * (ee22/ee6) + ee7) * ee2/
          ee4)/ee4 - ((3 * ee134 + ee2 * (ee132 - ee133) - (((ee6 * (ee83 +
            4 - (ee59 + 6) * ee2/ee4)/ee11 - ee103) * ee2/ee4 + 19)/
              ee4 + 3 * ee135)) * ee2 + 7 + ee126)) * ee2);
      out(j, 4) += -(((y * (1 - ((((1 - ee106) * ee16 + 2 + ee111 -
        (((2 * (ee55 + ee131) + ee115 - 8 * ee131)/ee6 + ee117) *
        ee2/ee56 + ee113) * ee3/ee4) * ee2/ee6 + ee114 + ee21 + ee122 -
        ee90)/ee4 + 2) * ee3/ee4)/ee6 - (((((2 * (1 - 3 * ee15) +
        8 * ee16 + 8 * (ee23 * ee3/(R_pow(ee4, 3) * ee11))) * ee23/
          ee28 - ee113) * ee3/ee4 + ee77 * ee16 + 2 + ee111) * ee23/
            ee32 + ee90 - (ee114 + ee21 + ee122))/ee4 - 2) * ee3)/ee4 - 1) *
              ee2 * ee3/ee4);
      out(j, 5) += ((y * (1 - (((((((2 * (ee49 - 1) - ee116) * ee3/
        ee56 + ee117 - 2) * ee2 - (ee33/ee4 + 2) * ee3)/ee4 + 2 + ee108)/
          ee6 + ee40 + 2) * ee3/ee4 + 1) * ee2 + (ee20 + ee21 -
            ((2 + 4 * (ee27/ee4)) * ee2 + ee121 * ee3/ee4)) * ee3/ee4)/
              ee4)/ee6 - ((((((((2 * (2 - ee49) + ee66 - ee68) * ee3/ee4 -
                ee115) * ee2/ee32 - 2) * ee2 - ee38)/ee4 + (1 - (2 + ee112) *
                ee2/ee4) * ee16 + 1 + ee108) * ee3/ee4 - (1 - ee70) * ee6) *
                ee2/ee32 - 2) * ee2 + ((((ee104 * ee2/ee11 + 4 * ee27)/ee4 +
                2) * ee2 + (ee121/ee4 + 8) * ee3 - (4 * ee19 + ee123))/
                  ee4 - 2) * ee3))/ee4 - 2) * ee2 * ee3/ee4;
      out(j, 6) += -(((y * (((((((ee82 - 4 * ee45)/ee6 - ee116) *
        ee2/ee4 - ee110) * ee3/ee4 + ee91 + ee22 * ee16)/ee6 + ((ee63 +
        6) * ee2 - (ee38 + ee12 + ee9))/ee4 - 1) * ee2 - ((6 - (ee42 +
        ee125 - ee85)/ee4) * ee2 + ee9) * ee3/ee4)/ee4 + 1)/
          ee6 - ((((((1/ee11 - 64) * ee2 + ee42 + ee125)/ee4 + 6) * ee2 -
            ee123)/ee4 - 2) * ee3 + (((((ee110 + ee66 - (ee6 * (4 + 8 *
            ee17 - ee68) + ee81) * ee23/ee28) * ee3/ee4 - (ee16 * (3 -
            (ee59 + 5) * ee2/ee4) + 2 * ee91))/ee11 + 8) * ee2 - ee30)/
              ee4 - 6) * ee2))/ee4 - 4) * ee2 * ee3/ee4);
      out(j, 7) += -(((y * (((((ee82 - ee109)/ee6 + (ee30 + 48 * ee2 -
        ee118)/ee4 - (ee84 + 6)) * ee2 + ee64 + ee21 + 6 * ee10)/
          ee4 + 1 - ((ee44 + ee109 - 7) * ee2/ee4 + 1 + 2 * (ee107 *
            ee6) + 3 * (ee22 * ee17))/ee6) * ee2/ee4 - 1)/ee6 - ((((ee17 *
            (1 + ee83 - ee89) + ee6 * (2 * ee107 + 7 - ((((ee6 * (6 -
            ee68) + ee81 + 4 * ee54)/ee11 - 8) * ee2 + ee12 + ee9)/ee4 +
            14) * ee2/ee4) + 7 - (19 - ee103 * ee2/ee4) * ee2/ee4)/ee11 -
            ((ee118 - ee85)/ee4 + 32)) * ee2 + 10 * ee4 + 12 * ee10 +
            ee64)/ee4 + 14) * ee2)/ee4 + 8) * ee2 * ee3/ee4);
      out(j, 8) += -((((6 * ee132 + ee2 * (Rf_psigamma(ee39, 3) - Rf_psigamma(ee2, 3)) -
        6 * ee133) * ee2 + 7 * ee134 - (((((1 - ee89) *
        ee22 + ee6 * (10 + 2 * ee62 + ee124 - (((ee6 * (8 - ee68) +
        ee81 + ee124) * ee6/ee11 - 8) * ee2/ee4 + 18) * ee2/ee4) +
        2 * ee72)/ee11 - ((ee64 + 22 * ee10 + 36 * ee4 - ((ee119 -
        ee85)/ee4 + 56) * ee2)/ee4 + 50)) * ee2/ee4 + 65)/ee4 + 7 *
        ee135)) * ee2 + 15 + ee126 + y * ((((ee82 - 4 * ee22)/ee6 -
        ee86) * ee2/ee4 - (2 * ee87 + 4 * (ee22 * ee2/ee4))/ee6) *
        ee2/ee4 - ((((((ee26 + 56 * ee2 - ee119)/ee4 - 18) * ee2 +
        14 * ee10 + ee64 + 20 * ee4)/ee4 + ee83 + 7) * ee2/ee4 + ee87 -
        14) * ee2/ee4 + 1)/ee6)/ee4 - ee129) * ee2);
    
        }

  }

  return out;
  
}

// //' @param pars a list of vectors of coefficients for negbinson log mean
// //' @param X1 a sparse design matrix for the negbinson log mean parameter
// //' @param ymat a matrix
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return negbinodspd0 a scalar, the negative log-likelihood
// //' @return negbinodpdd12 a matrix, first then second derivatives w.r.t. negbinson parameters
// //' @return negbinodpdd34 a matrix, third then fourth derivatives w.r.t. negbinson parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double negbinodspd0(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::mat ymat, 
                  arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat wmat) {

arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);

int nobs = nhere.size();

if (dcate == 1) {
  p1vec = p1vec.elem(dupid);
  p2vec = p2vec.elem(dupid);
}

double y, w, pars1, pars2, sigsq;
// double mu, alpha, ee1;
double mu, size, p;
double nllh = 0.0;

for (int j=0; j < nobs; j++) {
  
  pars1 = p1vec[j];
  pars2 = p2vec[j];
  
  for (int l=0; l < nhere[j]; l++) {
    
    y = ymat(j, l);
    w = wmat(j, l);
    
    mu = exp(pars1);
    size = exp(-pars2);
    p = size / (size + mu);
    nllh -= lgamma(y + size) - lgamma(size) - lgamma(y + 1) + size * log(p) + y * log(1 - p);
    
  } 
  
}

return(nllh);

}

// //' @rdname negbinodspd0
// [[Rcpp::export]]
arma::mat negbinodspd12(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::mat ymat, 
                     arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat wmat)
{
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, w, pars1, pars2;

  double ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12;

  arma::mat out = arma::mat(nobs, 5, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      w = wmat(j, l);
      
      ee2 = exp(-pars2);
      ee3 = exp(pars1);
      ee4 = ee2 + ee3;
      ee5 = ee2/ee4;
      ee6 = 1 - ee5;
      ee7 = ee2 + y;
      ee8 = R_pow(ee5, 2);
      ee9 = ee6 * ee4;
      ee10 = Rf_digamma(ee7);
      ee11 = Rf_digamma(ee2);
      ee12 = log(ee4);

      out(j, 0) += w * (-(ee2 * ee3 * (y/ee9 - 1)/ee4));
      out(j, 1) += w * (-((ee7/ee4 + ee11 + ee12 + pars2 - (1 + ee10)) *
        ee2));
      out(j, 2) += w * (-(((y * (1 - (2 + ee2/ee9) * ee3/ee4)/ee6 -
        (R_pow(ee2, 2)/(R_pow(ee4, 2) * ee8) - 2) * ee3)/ee4 - 1) *
        ee2 * ee3/ee4));
      out(j, 3) += w * ((((ee6 * ee2/(ee4 * ee8) + 2) * ee2 + y)/
        ee4 - 2) * ee2 * ee3/ee4);
      out(j, 4) += w * ((ee11 + ee12 + pars2 + y * ((1 - (3 - 2 *
        ee5) * ee2/ee4)/ee6 + ee5)/ee4 - (3 + ee10 + ee2 * (Rf_trigamma(ee7) -
        (((R_pow(ee6, 2)/ee8 - 2) * ee2/ee4 + 5)/ee4 + Rf_trigamma(ee2))))) *
        ee2);
    }
    
  }
  
  return out;
  
}

// //' @rdname negbinodspd0
// [[Rcpp::export]]
arma::mat negbinodspd34(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::mat ymat, 
                      arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat wmat)
{
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, w, pars1, pars2;

  double ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee14, ee15, ee16, ee17, ee19;
  double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
  double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee38, ee39;
  double ee40, ee42, ee44, ee45, ee47, ee48, ee49;
  double ee51, ee52, ee54, ee55, ee56, ee59;
  double ee62, ee63, ee64, ee65, ee66, ee67, ee68;
  double ee70, ee72, ee77;
  double ee81, ee82, ee83, ee84, ee85, ee86, ee87, ee89;
  double ee90, ee91;
  double ee103, ee104, ee106, ee107, ee108, ee109;
  double ee110, ee111, ee112, ee113, ee114, ee115, ee116, ee117, ee118, ee119;
  double ee121, ee122, ee123, ee124, ee125, ee126, ee129;
  double ee131, ee132, ee133, ee134, ee135;

  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      w = wmat(j, l);
      
      ee2 = exp(-pars2);
      ee3 = exp(pars1);
      ee4 = ee2 + ee3;
      ee5 = ee2/ee4;
      ee6 = 1 - ee5;
      ee7 = 2 * ee5;
      ee8 = 2 * ee2;
      ee9 = 2 * ee4;
      ee10 = ee8 + ee3;
      ee11 = R_pow(ee5, 2);
      ee12 = 2 * ee10;
      ee14 = (3 - ee7) * ee2/ee4;
      ee15 = ee3/ee4;
      ee16 = 1 - 2 * ee15;
      ee17 = 1 - ee7;
      ee19 = 2 * ee3 + ee2;
      ee20 = 2 * ee19;
      ee21 = 4 * ee4;
      ee22 = 1 - ee14;
      ee23 = R_pow(ee2, 2);
      ee24 = R_pow(ee4, 2);
      ee25 = 8 * ee2;
      ee26 = ee12 + ee21;
      ee27 = ee20 + ee9;
      ee28 = ee24 * ee11;
      ee29 = 4 * ee2;
      ee30 = ee12 + ee9;
      ee31 = 6 * ee2;
      ee32 = ee4 * ee11;
      ee33 = ee9 - ee31;
      ee34 = 8 * ee3;
      ee35 = R_pow(ee6, 2);
      ee36 = (ee26 - ee25)/ee4;
      ee38 = ee33 * ee3/ee4;
      ee39 = ee2 + y;
      ee40 = (ee27 - ee34)/ee4;
      ee42 = 4 * (ee12 + ee29);
      ee44 = (ee36 + 6) * ee2/ee4;
      ee45 = ee17 * ee6;
      ee47 = ee6 * ee23/ee28;
      ee48 = ee40 + 2;
      ee49 = 3 * ee5;
      ee51 = ((ee30 - ee25)/ee4 + 2) * ee2/ee4;
      ee52 = (ee38 + ee8)/ee4;
      ee54 = ee45 + 1 - ee14;
      ee55 = ee16 * ee6;
      ee56 = ee6 * ee4;
      ee59 = (2 * (ee35/ee11) - 2) * ee2/ee4;
      ee62 = 1 - (7 - ee44) * ee2/ee4;
      ee63 = 2 * ee16;
      ee64 = 2 * (ee29 + ee3);
      ee65 = ee9 + ee8;
      ee66 = 4 * ee17;
      ee67 = 4 * ee3;
      ee68 = 8 * ee47;
      ee70 = ee48 * ee3/ee4;
      ee72 = R_pow(ee22, 2) + ee62 * ee6;
      ee77 = (2 * (ee23/ee28) - 2) * ee3/ee4 + 1;
      ee81 = 2 * (ee35 + 1 - ee14);
      ee82 = 2 * ((4 - ee49) * ee2/ee4 - 1);
      ee83 = 2 * ee22;
      ee84 = 2 * ee17;
      ee85 = 64 * ee2;
      ee86 = 8 * ee5;
      ee87 = ee72/ee6;
      ee89 = (ee59 + 3) * ee2/ee4;
      ee90 = ((4 * (ee20 + ee67) + 8 * ee27 - 64 * ee3)/ee4 +  8) * ee3;
      ee91 = (1 - ee52) * ee6;
      ee103 = ee36 + 12;
      ee104 = ee77 + ee63;
      ee106 = (2 + 2 * (ee2/ee56)) * ee3/ee4;
      ee107 = 1 - ee51;
      ee108 = 2 * (ee17 * ee16 + 1 - ee52);
      ee109 = 2 * ee54;
      ee110 = 2 * (R_pow(ee17, 2) + 1 - ee51);
      ee111 = 2 * (R_pow(ee16, 2) + 1 - ee70);
      ee112 = 2 * (ee6 * ee2/ee32);
      ee113 = 2 * ee48;
      ee114 = 2 * (ee67 + ee2);
      ee115 = 4 * ee55;
      ee116 = ee66 + ee86;
      ee117 = 4 * ee16;
      ee118 = ee42 + 8 * ee30;
      ee119 = ee42 + 8 * ee26;
      ee121 = 4 * ee65 - 40 * ee2;
      ee122 = 6 * ee19;
      ee123 = 6 * ee4;
      ee124 = 8 * ee22;
      ee125 = 8 * ee65;
      ee126 = Rf_digamma(ee39);
      ee129 = Rf_digamma(ee2) + log(ee4) + pars2;
      ee131 = ee2 * ee3/ee24;
      ee132 = Rf_psigamma(ee39, 2);
      ee133 = Rf_psigamma(ee2, 2);
      ee134 = Rf_trigamma(ee39);
      ee135 = Rf_trigamma(ee2);

      out(j, 0) += -(((y * (1 - (((1 + ee63 - ee106) * ee2/ee6 + ee20 +
        ee9 - ee34)/ee4 + 2) * ee3/ee4)/ee6 - ((ee104 * ee23/
          ee32 + ee34 - ee27)/ee4 - 2) * ee3)/ee4 - 1) * ee2 * ee3/ee4);
      out(j, 1) += ((y * (1 - ((((ee84 + ee7)/ee6 + 2) * ee3/ee4 +
        1) * ee2 + ee38)/ee4)/ee6 - ((((ee84 - 2 * ee47) * ee3/ee4 -
        ee55) * ee2/ee32 - 2) * ee2 + (((6 + ee2/ee32) * ee2 - ee9)/
          ee4 - 2) * ee3))/ee4 - 2) * ee2 * ee3/ee4;
      out(j, 2) += -(((((((ee6 * (2 - (ee112 + 4) * ee2/ee4) + 3 -
        (5 - ee7) * ee2/ee4)/ee11 - 8) * ee2 + ee12 + ee9)/ee4 + 6) *
        ee2 + y * ((ee54/ee6 + (ee25 - ee30)/ee4 - 1) * ee2/ee4 +
        1)/ee6)/ee4 - 4) * ee2 * ee3/ee4);
      out(j, 3) += -((ee129 + y * (((((ee26 - ee31)/ee4 + 3) * ee2/
        ee4 - 6) * ee2/ee4 + 1)/ee6 + (2 * (ee22/ee6) + ee7) * ee2/
          ee4)/ee4 - ((3 * ee134 + ee2 * (ee132 - ee133) - (((ee6 * (ee83 +
            4 - (ee59 + 6) * ee2/ee4)/ee11 - ee103) * ee2/ee4 + 19)/
              ee4 + 3 * ee135)) * ee2 + 7 + ee126)) * ee2);
      out(j, 4) += -(((y * (1 - ((((1 - ee106) * ee16 + 2 + ee111 -
        (((2 * (ee55 + ee131) + ee115 - 8 * ee131)/ee6 + ee117) *
        ee2/ee56 + ee113) * ee3/ee4) * ee2/ee6 + ee114 + ee21 + ee122 -
        ee90)/ee4 + 2) * ee3/ee4)/ee6 - (((((2 * (1 - 3 * ee15) +
        8 * ee16 + 8 * (ee23 * ee3/(R_pow(ee4, 3) * ee11))) * ee23/
          ee28 - ee113) * ee3/ee4 + ee77 * ee16 + 2 + ee111) * ee23/
            ee32 + ee90 - (ee114 + ee21 + ee122))/ee4 - 2) * ee3)/ee4 - 1) *
              ee2 * ee3/ee4);
      out(j, 5) += ((y * (1 - (((((((2 * (ee49 - 1) - ee116) * ee3/
        ee56 + ee117 - 2) * ee2 - (ee33/ee4 + 2) * ee3)/ee4 + 2 + ee108)/
          ee6 + ee40 + 2) * ee3/ee4 + 1) * ee2 + (ee20 + ee21 -
            ((2 + 4 * (ee27/ee4)) * ee2 + ee121 * ee3/ee4)) * ee3/ee4)/
              ee4)/ee6 - ((((((((2 * (2 - ee49) + ee66 - ee68) * ee3/ee4 -
                ee115) * ee2/ee32 - 2) * ee2 - ee38)/ee4 + (1 - (2 + ee112) *
                ee2/ee4) * ee16 + 1 + ee108) * ee3/ee4 - (1 - ee70) * ee6) *
                ee2/ee32 - 2) * ee2 + ((((ee104 * ee2/ee11 + 4 * ee27)/ee4 +
                2) * ee2 + (ee121/ee4 + 8) * ee3 - (4 * ee19 + ee123))/
                  ee4 - 2) * ee3))/ee4 - 2) * ee2 * ee3/ee4;
      out(j, 6) += -(((y * (((((((ee82 - 4 * ee45)/ee6 - ee116) *
        ee2/ee4 - ee110) * ee3/ee4 + ee91 + ee22 * ee16)/ee6 + ((ee63 +
        6) * ee2 - (ee38 + ee12 + ee9))/ee4 - 1) * ee2 - ((6 - (ee42 +
        ee125 - ee85)/ee4) * ee2 + ee9) * ee3/ee4)/ee4 + 1)/
          ee6 - ((((((1/ee11 - 64) * ee2 + ee42 + ee125)/ee4 + 6) * ee2 -
            ee123)/ee4 - 2) * ee3 + (((((ee110 + ee66 - (ee6 * (4 + 8 *
            ee17 - ee68) + ee81) * ee23/ee28) * ee3/ee4 - (ee16 * (3 -
            (ee59 + 5) * ee2/ee4) + 2 * ee91))/ee11 + 8) * ee2 - ee30)/
              ee4 - 6) * ee2))/ee4 - 4) * ee2 * ee3/ee4);
      out(j, 7) += -(((y * (((((ee82 - ee109)/ee6 + (ee30 + 48 * ee2 -
        ee118)/ee4 - (ee84 + 6)) * ee2 + ee64 + ee21 + 6 * ee10)/
          ee4 + 1 - ((ee44 + ee109 - 7) * ee2/ee4 + 1 + 2 * (ee107 *
            ee6) + 3 * (ee22 * ee17))/ee6) * ee2/ee4 - 1)/ee6 - ((((ee17 *
            (1 + ee83 - ee89) + ee6 * (2 * ee107 + 7 - ((((ee6 * (6 -
            ee68) + ee81 + 4 * ee54)/ee11 - 8) * ee2 + ee12 + ee9)/ee4 +
            14) * ee2/ee4) + 7 - (19 - ee103 * ee2/ee4) * ee2/ee4)/ee11 -
            ((ee118 - ee85)/ee4 + 32)) * ee2 + 10 * ee4 + 12 * ee10 +
            ee64)/ee4 + 14) * ee2)/ee4 + 8) * ee2 * ee3/ee4);
      out(j, 8) += -((((6 * ee132 + ee2 * (Rf_psigamma(ee39, 3) - Rf_psigamma(ee2, 3)) -
        6 * ee133) * ee2 + 7 * ee134 - (((((1 - ee89) *
        ee22 + ee6 * (10 + 2 * ee62 + ee124 - (((ee6 * (8 - ee68) +
        ee81 + ee124) * ee6/ee11 - 8) * ee2/ee4 + 18) * ee2/ee4) +
        2 * ee72)/ee11 - ((ee64 + 22 * ee10 + 36 * ee4 - ((ee119 -
        ee85)/ee4 + 56) * ee2)/ee4 + 50)) * ee2/ee4 + 65)/ee4 + 7 *
        ee135)) * ee2 + 15 + ee126 + y * ((((ee82 - 4 * ee22)/ee6 -
        ee86) * ee2/ee4 - (2 * ee87 + 4 * (ee22 * ee2/ee4))/ee6) *
        ee2/ee4 - ((((((ee26 + 56 * ee2 - ee119)/ee4 - 18) * ee2 +
        14 * ee10 + ee64 + 20 * ee4)/ee4 + ee83 + 7) * ee2/ee4 + ee87 -
        14) * ee2/ee4 + 1)/ee6)/ee4 - ee129) * ee2);

    }
    
  }
  
  return out;
  
}
