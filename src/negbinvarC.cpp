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
// //' @return negbind0 a scalar, the negative log-likelihood
// //' @return negbind12 a matrix, first then second derivatives w.r.t. negbinson parameters
// //' @return negbind34 a matrix, third then fourth derivatives w.r.t. negbinson parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double negbind0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, 
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
  double mu, size, p, sigsq;
  double nllh = 0.0;
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      w = wmat(j, l);

      mu = exp(pars1);
      // size = exp(pars2);
      // size = exp(-pars2);
      sigsq = mu + exp(pars2);
      size = mu * mu / (sigsq - mu);
      p = size / (size + mu);
      nllh -= lgamma(y + size) - lgamma(size) - lgamma(y + 1) + size * log(p) + y * log(1 - p);

    } 
    
  }
  
  return(nllh);
  
}

// //' @rdname negbind0
// [[Rcpp::export]]
arma::mat negbind12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, 
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
  // double ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  // double ee10, ee11, ee12;
  double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee19;
  double ee20, ee21, ee22, ee23, ee27, ee29;
  double ee30, ee32, ee33, ee34, ee35, ee38;
  double ee41, ee42, ee43, ee44, ee45, ee46, ee47;

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

      // ee2 = exp(-pars2);
      // ee3 = exp(pars1);
      // ee4 = ee2 + ee3;
      // ee5 = ee2/ee4;
      // ee6 = 1 - ee5;
      // ee7 = ee2 + y;
      // ee8 = R_pow(ee5, 2);
      // ee9 = ee6 * ee4;
      // ee10 = Rf_digamma(ee7);
      // ee11 = Rf_digamma(ee2);
      // ee12 = log(ee4);
      // 
      // out(j, 0) += w * (-(ee2 * ee3 * (y/ee9 - 1)/ee4));
      // out(j, 1) += w * (-((ee7/ee4 + ee11 + ee12 + pars2 - (1 + ee10)) *
      //   ee2));
      // out(j, 2) += w * (-(((y * (1 - (2 + ee2/ee9) * ee3/ee4)/ee6 -
      //   (R_pow(ee2, 2)/(R_pow(ee4, 2) * ee8) - 2) * ee3)/ee4 - 1) *
      //   ee2 * ee3/ee4));
      // out(j, 3) += w * ((((ee6 * ee2/(ee4 * ee8) + 2) * ee2 + y)/
      //   ee4 - 2) * ee2 * ee3/ee4);
      // out(j, 4) += w * ((ee11 + ee12 + pars2 + y * ((1 - (3 - 2 *
      //   ee5) * ee2/ee4)/ee6 + ee5)/ee4 - (3 + ee10 + ee2 * (Rf_trigamma(ee7) -
      //   (((R_pow(ee6, 2)/ee8 - 2) * ee2/ee4 + 5)/ee4 + Rf_trigamma(ee2))))) *
      //   ee2);
      
      ee1 = exp(pars1);
      ee2 = exp(pars2);
      ee3 = ee1/ee2;
      ee4 = 1 + ee3;
      ee5 = R_pow(ee1, 2);
      ee6 = R_pow((ee4 * ee1), 2);
      ee7 = ee5/ee2;
      ee8 = 1 + 2 * ee3;
      ee10 = ee1/(ee4 * ee2);
      ee11 = ee6 * ee2;
      ee12 = R_pow(ee1, 3);
      ee13 = 1 - ee10;
      ee14 = 1/ee4;
      ee15 = 2/ee4;
      ee16 = ee7 + y;
      ee19 = ee14 - ee12/ee11;
      ee20 = ee15 - ee8 * ee5/ee6;
      ee21 = Rf_digamma(ee16);
      ee22 = Rf_digamma(ee7);
      ee23 = log1p(ee3);
      ee27 = 2 * (ee8 * ee4 * ee5/ee6);
      ee29 = R_pow(ee10, 2) * ee2;
      ee30 = 4/ee4;
      ee32 = pars1 - (ee23 + pars2);
      ee33 = Rf_trigamma(ee16);
      ee34 = Rf_trigamma(ee7);
      ee35 = ee13 * ee2;
      ee38 = ee8 * (4 - ee27) + 1 + 4 * ee3;
      ee41 = ee19 * ee20;
      ee42 = R_pow(ee19, 2);
      ee43 = R_pow(ee20, 2);
      ee44 = 2 * (ee4 * ee12/ee11);
      ee45 = 2 * ee32;
      ee46 = 2 * ee21;
      ee47 = 2 * ee22;

      out(j, 0) += w * (-(((ee4 * ee20 + ee45 + ee46 - ee47) * ee1 -
        y * ee20/ee13) * ee1/ee2));
      out(j, 1) += w * (-(ee1 * (y * ee19/ee13 - (ee4 * ee19 + ee21 +
        pars1 - (ee22 + ee23 + pars2)) * ee1)/ee2));
      out(j, 2) += w * (-(((ee4 * (4 * ee20 + ee30 - ee38 * ee5/ee6) +
        (4 * ee33 - (ee43/ee29 + 4 * ee34)) * ee5/ee2 + 4 * ee32 +
        4 * ee21 - 4 * ee22) * ee1 - y * ((ee43/ee35 - ee38 * ee1/
          ee6) * ee1 + ee30)/ee13) * ee1/ee2));
      out(j, 3) += w * (((ee4 * (2 * ee19 + ee30 - ((8 - ee27) * ee1/
        ee2 + 2) * ee5/ee6) + (2 * ee33 - (ee41/ee29 + 2 * ee34)) *
          ee5/ee2 + ee45 + ee46 - ee47) * ee1 - y * ((ee41/ee35 - ((6 -
          ee27) * ee1/ee2 + 1) * ee1/ee6) * ee1 + ee15)/ee13) * ee1/
            ee2);
      out(j, 4) += w * (-(ee1 * (y * (((3 - ee44) * ee5/ee6 - ee42/
        ee13) * ee1/ee2 - ee14)/ee13 - ((ee42/ee29 + ee34 - ee33) *
          ee5/ee2 + ((5 - ee44) * ee12/ee11 - 3/ee4) * ee4 + ee22 + ee23 +
          pars2 - (ee21 + pars1)) * ee1)/ee2));
    
    }
    
  }
  
  return out;
  
}

// //' @rdname negbind0
// [[Rcpp::export]]
arma::mat negbind34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, 
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

  // double ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  // double ee10, ee11, ee12, ee14, ee15, ee16, ee17, ee19;
  // double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
  // double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee38, ee39;
  // double ee40, ee42, ee44, ee45, ee47, ee48, ee49;
  // double ee51, ee52, ee54, ee55, ee56, ee59;
  // double ee62, ee63, ee64, ee65, ee66, ee67, ee68;
  // double ee70, ee72, ee77;
  // double ee81, ee82, ee83, ee84, ee85, ee86, ee87, ee89;
  // double ee90, ee91;
  // double ee103, ee104, ee106, ee107, ee108, ee109;
  // double ee110, ee111, ee112, ee113, ee114, ee115, ee116, ee117, ee118, ee119;
  // double ee121, ee122, ee123, ee124, ee125, ee126, ee129;
  // double ee131, ee132, ee133, ee134, ee135;
  
  double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee17, ee19;
  double ee20, ee21, ee23, ee24, ee25, ee26, ee27, ee29;
  double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee37, ee38;
  double ee40, ee43, ee45, ee48;
  double ee51, ee52, ee55, ee56, ee57, ee58, ee59;
  double ee60, ee61, ee62, ee63, ee65, ee66, ee67, ee68, ee69;
  double ee70, ee71, ee72, ee73, ee74, ee75, ee76, ee77, ee78, ee79;
  double ee80, ee81, ee82, ee83, ee84, ee85, ee87, ee88, ee89;
  double ee90, ee98, ee99;
  double ee100, ee101, ee102, ee104, ee106, ee107;
  double ee110, ee111, ee112, ee113, ee114, ee117;
  double ee122, ee123, ee124, ee125, ee127, ee129;
  double ee133, ee138;
  double ee141, ee143, ee147, ee148, ee149;
  double ee150, ee151, ee153, ee154, ee155, ee156, ee157, ee158;
  double ee160, ee161, ee162, ee163, ee166, ee167, ee168;
  double ee172, ee173, ee174, ee175, ee176, ee177, ee178, ee179;
  double ee180, ee181, ee182, ee183, ee185, ee187, ee189;
  double ee194, ee195, ee197, ee198;
  double ee200, ee201, ee203, ee205, ee207;
  double ee210, ee211, ee212, ee213, ee215, ee217, ee218;
  double ee220, ee222, ee223, ee224, ee225, ee226, ee229;
  double ee234;
  double ee242, ee244, ee246, ee248;
  double ee251, ee253, ee255, ee256, ee257, ee258;
  double ee263, ee266, ee267, ee268;
  double ee271, ee273, ee274, ee275, ee276;
  double ee281, ee282, ee284, ee285, ee286, ee287, ee288, ee289;
  double ee290, ee291, ee292, ee293, ee294, ee295, ee296, ee297, ee298, ee299;
  double ee300, ee302, ee303, ee304, ee305, ee306, ee307, ee308, ee309;
  double ee310, ee311, ee312, ee313, ee314, ee315, ee316, ee317, ee318, ee319;
  double ee320, ee321, ee322;
  
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

      // ee2 = exp(-pars2);
      // ee3 = exp(pars1);
      // ee4 = ee2 + ee3;
      // ee5 = ee2/ee4;
      // ee6 = 1 - ee5;
      // ee7 = 2 * ee5;
      // ee8 = 2 * ee2;
      // ee9 = 2 * ee4;
      // ee10 = ee8 + ee3;
      // ee11 = R_pow(ee5, 2);
      // ee12 = 2 * ee10;
      // ee14 = (3 - ee7) * ee2/ee4;
      // ee15 = ee3/ee4;
      // ee16 = 1 - 2 * ee15;
      // ee17 = 1 - ee7;
      // ee19 = 2 * ee3 + ee2;
      // ee20 = 2 * ee19;
      // ee21 = 4 * ee4;
      // ee22 = 1 - ee14;
      // ee23 = R_pow(ee2, 2);
      // ee24 = R_pow(ee4, 2);
      // ee25 = 8 * ee2;
      // ee26 = ee12 + ee21;
      // ee27 = ee20 + ee9;
      // ee28 = ee24 * ee11;
      // ee29 = 4 * ee2;
      // ee30 = ee12 + ee9;
      // ee31 = 6 * ee2;
      // ee32 = ee4 * ee11;
      // ee33 = ee9 - ee31;
      // ee34 = 8 * ee3;
      // ee35 = R_pow(ee6, 2);
      // ee36 = (ee26 - ee25)/ee4;
      // ee38 = ee33 * ee3/ee4;
      // ee39 = ee2 + y;
      // ee40 = (ee27 - ee34)/ee4;
      // ee42 = 4 * (ee12 + ee29);
      // ee44 = (ee36 + 6) * ee2/ee4;
      // ee45 = ee17 * ee6;
      // ee47 = ee6 * ee23/ee28;
      // ee48 = ee40 + 2;
      // ee49 = 3 * ee5;
      // ee51 = ((ee30 - ee25)/ee4 + 2) * ee2/ee4;
      // ee52 = (ee38 + ee8)/ee4;
      // ee54 = ee45 + 1 - ee14;
      // ee55 = ee16 * ee6;
      // ee56 = ee6 * ee4;
      // ee59 = (2 * (ee35/ee11) - 2) * ee2/ee4;
      // ee62 = 1 - (7 - ee44) * ee2/ee4;
      // ee63 = 2 * ee16;
      // ee64 = 2 * (ee29 + ee3);
      // ee65 = ee9 + ee8;
      // ee66 = 4 * ee17;
      // ee67 = 4 * ee3;
      // ee68 = 8 * ee47;
      // ee70 = ee48 * ee3/ee4;
      // ee72 = R_pow(ee22, 2) + ee62 * ee6;
      // ee77 = (2 * (ee23/ee28) - 2) * ee3/ee4 + 1;
      // ee81 = 2 * (ee35 + 1 - ee14);
      // ee82 = 2 * ((4 - ee49) * ee2/ee4 - 1);
      // ee83 = 2 * ee22;
      // ee84 = 2 * ee17;
      // ee85 = 64 * ee2;
      // ee86 = 8 * ee5;
      // ee87 = ee72/ee6;
      // ee89 = (ee59 + 3) * ee2/ee4;
      // ee90 = ((4 * (ee20 + ee67) + 8 * ee27 - 64 * ee3)/ee4 +  8) * ee3;
      // ee91 = (1 - ee52) * ee6;
      // ee103 = ee36 + 12;
      // ee104 = ee77 + ee63;
      // ee106 = (2 + 2 * (ee2/ee56)) * ee3/ee4;
      // ee107 = 1 - ee51;
      // ee108 = 2 * (ee17 * ee16 + 1 - ee52);
      // ee109 = 2 * ee54;
      // ee110 = 2 * (R_pow(ee17, 2) + 1 - ee51);
      // ee111 = 2 * (R_pow(ee16, 2) + 1 - ee70);
      // ee112 = 2 * (ee6 * ee2/ee32);
      // ee113 = 2 * ee48;
      // ee114 = 2 * (ee67 + ee2);
      // ee115 = 4 * ee55;
      // ee116 = ee66 + ee86;
      // ee117 = 4 * ee16;
      // ee118 = ee42 + 8 * ee30;
      // ee119 = ee42 + 8 * ee26;
      // ee121 = 4 * ee65 - 40 * ee2;
      // ee122 = 6 * ee19;
      // ee123 = 6 * ee4;
      // ee124 = 8 * ee22;
      // ee125 = 8 * ee65;
      // ee126 = Rf_digamma(ee39);
      // ee129 = Rf_digamma(ee2) + log(ee4) + pars2;
      // ee131 = ee2 * ee3/ee24;
      // ee132 = Rf_psigamma(ee39, 2);
      // ee133 = Rf_psigamma(ee2, 2);
      // ee134 = Rf_trigamma(ee39);
      // ee135 = Rf_trigamma(ee2);
      // 
      // out(j, 0) += -(((y * (1 - (((1 + ee63 - ee106) * ee2/ee6 + ee20 +
      //   ee9 - ee34)/ee4 + 2) * ee3/ee4)/ee6 - ((ee104 * ee23/
      //     ee32 + ee34 - ee27)/ee4 - 2) * ee3)/ee4 - 1) * ee2 * ee3/ee4);
      // out(j, 1) += ((y * (1 - ((((ee84 + ee7)/ee6 + 2) * ee3/ee4 +
      //   1) * ee2 + ee38)/ee4)/ee6 - ((((ee84 - 2 * ee47) * ee3/ee4 -
      //   ee55) * ee2/ee32 - 2) * ee2 + (((6 + ee2/ee32) * ee2 - ee9)/
      //     ee4 - 2) * ee3))/ee4 - 2) * ee2 * ee3/ee4;
      // out(j, 2) += -(((((((ee6 * (2 - (ee112 + 4) * ee2/ee4) + 3 -
      //   (5 - ee7) * ee2/ee4)/ee11 - 8) * ee2 + ee12 + ee9)/ee4 + 6) *
      //   ee2 + y * ((ee54/ee6 + (ee25 - ee30)/ee4 - 1) * ee2/ee4 +
      //   1)/ee6)/ee4 - 4) * ee2 * ee3/ee4);
      // out(j, 3) += -((ee129 + y * (((((ee26 - ee31)/ee4 + 3) * ee2/
      //   ee4 - 6) * ee2/ee4 + 1)/ee6 + (2 * (ee22/ee6) + ee7) * ee2/
      //     ee4)/ee4 - ((3 * ee134 + ee2 * (ee132 - ee133) - (((ee6 * (ee83 +
      //       4 - (ee59 + 6) * ee2/ee4)/ee11 - ee103) * ee2/ee4 + 19)/
      //         ee4 + 3 * ee135)) * ee2 + 7 + ee126)) * ee2);
      // out(j, 4) += -(((y * (1 - ((((1 - ee106) * ee16 + 2 + ee111 -
      //   (((2 * (ee55 + ee131) + ee115 - 8 * ee131)/ee6 + ee117) *
      //   ee2/ee56 + ee113) * ee3/ee4) * ee2/ee6 + ee114 + ee21 + ee122 -
      //   ee90)/ee4 + 2) * ee3/ee4)/ee6 - (((((2 * (1 - 3 * ee15) +
      //   8 * ee16 + 8 * (ee23 * ee3/(R_pow(ee4, 3) * ee11))) * ee23/
      //     ee28 - ee113) * ee3/ee4 + ee77 * ee16 + 2 + ee111) * ee23/
      //       ee32 + ee90 - (ee114 + ee21 + ee122))/ee4 - 2) * ee3)/ee4 - 1) *
      //         ee2 * ee3/ee4);
      // out(j, 5) += ((y * (1 - (((((((2 * (ee49 - 1) - ee116) * ee3/
      //   ee56 + ee117 - 2) * ee2 - (ee33/ee4 + 2) * ee3)/ee4 + 2 + ee108)/
      //     ee6 + ee40 + 2) * ee3/ee4 + 1) * ee2 + (ee20 + ee21 -
      //       ((2 + 4 * (ee27/ee4)) * ee2 + ee121 * ee3/ee4)) * ee3/ee4)/
      //         ee4)/ee6 - ((((((((2 * (2 - ee49) + ee66 - ee68) * ee3/ee4 -
      //           ee115) * ee2/ee32 - 2) * ee2 - ee38)/ee4 + (1 - (2 + ee112) *
      //           ee2/ee4) * ee16 + 1 + ee108) * ee3/ee4 - (1 - ee70) * ee6) *
      //           ee2/ee32 - 2) * ee2 + ((((ee104 * ee2/ee11 + 4 * ee27)/ee4 +
      //           2) * ee2 + (ee121/ee4 + 8) * ee3 - (4 * ee19 + ee123))/
      //             ee4 - 2) * ee3))/ee4 - 2) * ee2 * ee3/ee4;
      // out(j, 6) += -(((y * (((((((ee82 - 4 * ee45)/ee6 - ee116) *
      //   ee2/ee4 - ee110) * ee3/ee4 + ee91 + ee22 * ee16)/ee6 + ((ee63 +
      //   6) * ee2 - (ee38 + ee12 + ee9))/ee4 - 1) * ee2 - ((6 - (ee42 +
      //   ee125 - ee85)/ee4) * ee2 + ee9) * ee3/ee4)/ee4 + 1)/
      //     ee6 - ((((((1/ee11 - 64) * ee2 + ee42 + ee125)/ee4 + 6) * ee2 -
      //       ee123)/ee4 - 2) * ee3 + (((((ee110 + ee66 - (ee6 * (4 + 8 *
      //       ee17 - ee68) + ee81) * ee23/ee28) * ee3/ee4 - (ee16 * (3 -
      //       (ee59 + 5) * ee2/ee4) + 2 * ee91))/ee11 + 8) * ee2 - ee30)/
      //         ee4 - 6) * ee2))/ee4 - 4) * ee2 * ee3/ee4);
      // out(j, 7) += -(((y * (((((ee82 - ee109)/ee6 + (ee30 + 48 * ee2 -
      //   ee118)/ee4 - (ee84 + 6)) * ee2 + ee64 + ee21 + 6 * ee10)/
      //     ee4 + 1 - ((ee44 + ee109 - 7) * ee2/ee4 + 1 + 2 * (ee107 *
      //       ee6) + 3 * (ee22 * ee17))/ee6) * ee2/ee4 - 1)/ee6 - ((((ee17 *
      //       (1 + ee83 - ee89) + ee6 * (2 * ee107 + 7 - ((((ee6 * (6 -
      //       ee68) + ee81 + 4 * ee54)/ee11 - 8) * ee2 + ee12 + ee9)/ee4 +
      //       14) * ee2/ee4) + 7 - (19 - ee103 * ee2/ee4) * ee2/ee4)/ee11 -
      //       ((ee118 - ee85)/ee4 + 32)) * ee2 + 10 * ee4 + 12 * ee10 +
      //       ee64)/ee4 + 14) * ee2)/ee4 + 8) * ee2 * ee3/ee4);
      // out(j, 8) += -((((6 * ee132 + ee2 * (Rf_psigamma(ee39, 3) - Rf_psigamma(ee2, 3)) -
      //   6 * ee133) * ee2 + 7 * ee134 - (((((1 - ee89) *
      //   ee22 + ee6 * (10 + 2 * ee62 + ee124 - (((ee6 * (8 - ee68) +
      //   ee81 + ee124) * ee6/ee11 - 8) * ee2/ee4 + 18) * ee2/ee4) +
      //   2 * ee72)/ee11 - ((ee64 + 22 * ee10 + 36 * ee4 - ((ee119 -
      //   ee85)/ee4 + 56) * ee2)/ee4 + 50)) * ee2/ee4 + 65)/ee4 + 7 *
      //   ee135)) * ee2 + 15 + ee126 + y * ((((ee82 - 4 * ee22)/ee6 -
      //   ee86) * ee2/ee4 - (2 * ee87 + 4 * (ee22 * ee2/ee4))/ee6) *
      //   ee2/ee4 - ((((((ee26 + 56 * ee2 - ee119)/ee4 - 18) * ee2 +
      //   14 * ee10 + ee64 + 20 * ee4)/ee4 + ee83 + 7) * ee2/ee4 + ee87 -
      //   14) * ee2/ee4 + 1)/ee6)/ee4 - ee129) * ee2);
    
    ee1 = exp(pars1);
    ee2 = exp(pars2);
    ee3 = ee1/ee2;
    ee4 = 1 + ee3;
    ee5 = R_pow((ee4 * ee1), 2);
    ee6 = R_pow(ee1, 2);
    ee7 = 2 * ee3;
    ee8 = 1 + ee7;
    ee9 = R_pow(ee1, 3);
    ee10 = ee5 * ee2;
    ee11 = 4 * ee3;
    ee12 = 1/ee4;
    ee13 = ee6/ee2;
    ee14 = 2/ee4;
    ee17 = ee8 * ee4 * ee6/ee5;
    ee19 = ee1/(ee4 * ee2);
    ee20 = 1 + ee11;
    ee21 = ee12 - ee9/ee10;
    ee23 = ee8 * ee6/ee5;
    ee24 = 2 * ee17;
    ee25 = 2 * ee8;
    ee26 = ee14 - ee23;
    ee27 = ee13 + y;
    ee29 = ee4 * ee9/ee10;
    ee30 = 1 - ee19;
    ee31 = R_pow(ee19, 2);
    ee32 = 4/ee4;
    ee33 = R_pow(ee4, 2);
    ee34 = 2 * ee29;
    ee35 = R_pow(ee8, 2);
    ee36 = ee20 * ee4;
    ee37 = 3 - ee34;
    ee38 = ee8 * (4 - ee24);
    ee40 = ee38 + 1 + ee11;
    ee43 = (6 - ee24) * ee1/ee2 + 1;
    ee45 = 1 + ee25 + ee11;
    ee48 = ee37 * ee9/ee10 - ee12;
    ee51 = ee32 - ee40 * ee6/ee5;
    ee52 = R_pow(ee2, 2);
    ee55 = ee14 - ee43 * ee6/ee5;
    ee56 = 1 + 2 * ee4;
    ee57 = 8 * ee3;
    ee58 = ee35 + ee36;
    ee59 = ee56 + ee7;
    ee60 = ee33 * ee9;
    ee61 = Rf_trigamma(ee27);
    ee62 = Rf_trigamma(ee13);
    ee63 = 4 - 8 * ee17;
    ee65 = ee59 * ee1/ee2;
    ee66 = Rf_psigamma(ee27, 2);
    ee67 = Rf_psigamma(ee13, 2);
    ee68 = ee31 * ee52;
    ee69 = 4 * ee8;
    ee70 = ee30 * ee2;
    ee71 = R_pow(ee21, 2);
    ee72 = 2 * ee58;
    ee73 = ee45 * ee4;
    ee74 = ee4 * ee31;
    ee75 = 2 * ee65;
    ee76 = 2 * ee20;
    ee77 = 4 * ee20;
    ee78 = ee21 * ee26;
    ee79 = 8/ee4;
    ee80 = ee74 * ee52;
    ee81 = ee60/ee2;
    ee82 = 8 * ee29;
    ee83 = ee60/ee10;
    ee84 = ee8 * ee63;
    ee85 = R_pow(ee26, 2);
    ee87 = ee63 * ee1/ee2;
    ee88 = 2 * ee36;
    ee89 = 4 * ee45;
    ee90 = 6 * ee20;
    ee98 = ee8 * (12 - ((ee84 + ee89) * ee4 + ee72) * ee6/ee5) +  1 + ee90 + ee57;
    ee99 = 2 * (ee73 * ee6/ee5);
    ee100 = 8 * ee83;
    ee101 = Rf_digamma(ee27);
    ee102 = Rf_digamma(ee13);
    ee104 = ee6 * ee66/ee2;
    ee106 = ee6 * ee67/ee2;
    ee107 = log1p(ee3);
    ee110 = (ee4 * (10 - ee82) + ee25) * ee9/ee10 - 7;
    ee111 = ee31 * ee2;
    ee112 = 1 + ee69;
    ee113 = 1 + ee57;
    ee114 = 16/ee4;
    ee117 = ee110 * ee9/ee10 + ee12;
    ee122 = ee8 * (4 - ((ee87 + ee76) * ee4 + ee75) * ee6/ee5) +  (20 - ee99) * ee1/ee2 + 1;
    ee123 = ee4 * (ee77 + ee11);
    ee124 = ee21 * ee51;
    ee125 = ee100 - ee25;
    ee127 = pars1 - (ee107 + pars2);
    ee129 = ee98 * ee6/ee5;
    ee133 = ee48 * ee26;
    ee138 = ee21 * ee55;
    ee141 = (14 - (ee123 - ee8 * ee125) * ee6/ee5) * ee1/ee2 +  1;
    ee143 = ee112 + ee77 + ee57;
    ee147 = 2 * ee73 + ee72;
    ee148 = ee75 + ee88;
    ee149 = 4 * ee4;
    ee150 = 4 * (2 * (ee5 * ee8) + 4 * ee81);
    ee151 = 64 * ee81;
    ee153 = ee71 * ee1/ee2;
    ee154 = 4 * ee61;
    ee155 = 4 * ee62;
    ee156 = ee79 - ee129;
    ee157 = Rf_psigamma(ee27, 3);
    ee158 = Rf_psigamma(ee13, 3);
    ee160 = ee48 * ee30;
    ee161 = ee8 * (ee151 - ee150);
    ee162 = ee78 * ee1;
    ee163 = ee78 * ee6;
    ee166 = ee55 * ee51;
    ee167 = ee85 * ee6;
    ee168 = R_pow(ee51, 2);
    ee172 = ee88 + 2 * ee35;
    ee173 = ee25 + ee149;
    ee174 = ee25 + ee11;
    ee175 = 2 * ee143;
    ee176 = ee76 + ee69;
    ee177 = 2 * ee51;
    ee178 = 2 * ee61;
    ee179 = 2 * ee62;
    ee180 = 4 * (ee5 * ee148);
    ee181 = 4 * (1 + 12 * ee3 + ee25);
    ee182 = 8 - 16 * ee17;
    ee183 = 8 * (ee21 * ee6/ee80);
    ee185 = ee6 * ee157/ee2;
    ee187 = ee6 * ee158/ee2;
    ee189 = ee122 * ee1/ee5;
    ee194 = ee141 * ee6/ee5 - ee14;
    ee195 = ee160 + ee153;
    ee197 = ee48 * ee4 + ee71 * ee6/ee68;
    ee198 = ee133 - ee138;
    ee200 = ee48/ee4 - ee71;
    ee201 = ee30 * ee55;
    ee203 = ee4 * ee55 - ee163/ee68;
    ee205 = ee4 * ee51 - ee167/ee68;
    ee207 = ee162/ee2;
    ee210 = 2 * (ee78/ee80);
    ee211 = 2 * (ee71/(ee74 * ee2));
    ee212 = 2 * ee48;
    ee213 = 2 * ee55;
    ee215 = 2 * ee104 + ee178;
    ee217 = 2 * ee106 + ee179;
    ee218 = 4 * ee65;
    ee220 = 4 * ee104 + ee154;
    ee222 = 4 * ee106 + ee155;
    ee223 = ee32 - ee122 * ee6/ee5;
    ee224 = 6 * ee3;
    ee225 = ee104 + ee61;
    ee226 = ee106 + ee62;
    ee229 = ee117 * ee26;
    ee234 = ((ee4 * (ee181 - (((4 * (2 * (ee5 * ee59) + 4 *  (ee8 * ee33 * ee6)) - 64 * (ee8 * ee33 * ee6)) * ee1/ee2 +  ee180) * ee8/ee5 + 4 * (ee147 * ee1/ee2)) * ee6/ee5) +  2 * ((ee112 + ee149 + ee11) * ee1/ee2)) * ee8 + ee45 *  (4 * ee59 - 8 * (ee8 * ee33 * ee6/ee5)) * ee1/ee2 + 2 *  (ee58 * ee20)) * ee1/ee5;
    ee242 = (((4 * ee148 - (ee161 - ee180)/ee5) * ee9/ee10 -  2 * ee113) * ee8 * ee4 - (2 * ((ee56 + ee224) * ee8) +  4 * (ee59 * ee20)) * ee1/ee2) * ee1/ee5;
    ee244 = ee98 * ee1/ee5;
    ee246 = (ee161 - 4 * (ee5 * ee172))/ee5 - 4 * ee172;
    ee248 = (ee4 * (20 - 16 * ee29) + ee69) * ee1/ee2;
    ee251 = ee141 * ee1/ee5;
    ee253 = (ee37 * ee1/ee5 + ee211) * ee6/ee2;
    ee255 = (ee37 * ee6/ee5 - 2 * (ee71/ee30)) * ee1/ee2;
    ee256 = ee48 * ee51;
    ee257 = ee30 * ee51;
    ee258 = ee8 * (16 * ee83 - ee69);
    ee263 = ee8 * (32 - ((ee8 * ee182 + ee175 + 4 * ee176) *  ee4 + 4 * ee58) * ee6/ee5) + (16/ee2 - ((ee4 * (4 * ee143 -  ((4 * (ee5 * ee147) + 4 * (2 * (ee58 * ee5) + 4 * (ee35 *  ee33 * ee6)) - 64 * (ee35 * ee33 * ee6))/ee5 + 4 *  ee147) * ee8 * ee6/ee5) + 2 * (ee113 * ee4 + 3 *  (ee8 * ee20))) * ee8 + ee45 * (6 * ee58 - 8 * (ee35 *  ee33 * ee6/ee5))) * ee1/ee5) * ee1 + 1 + 24 * ee20 +  8 * ee113;
    ee266 = ee8 * (8 - ((ee87 + ee25 + ee76) * ee4 + ee75) *  ee6/ee5) + (24 - ee99) * ee1/ee2 + 2;
    ee267 = ee8/ee5;
    ee268 = ee20 * (ee100 - 8 * ee8);
    ee271 = ee4 * (ee175 + 2 * ee176) * ee6/ee5;
    ee273 = ee4 * (4 - ee82) + ee25;
    ee274 = ee4 * (ee181 + 4 * ee174 + ee57);
    ee275 = ee21 * ee156;
    ee276 = (ee177 + ee32 - (ee40/ee5 + 2 * (ee85/ee80)) * ee6) *  ee26;
    ee281 = ee85 * ee1/ee2;
    ee282 = (ee151 - (4 * (ee5 * ee173) + ee150))/ee5;
    ee284 = ee182 * ee1/ee2;
    ee285 = 12 * ee20;
    ee286 = 12 * ee26;
    ee287 = 2 * (ee117 * ee21 + R_pow(ee48, 2));
    ee288 = 2 * ee117;
    ee289 = 2 * (ee194 * ee26 - R_pow(ee55, 2));
    ee290 = 2 * ee195;
    ee291 = 2 * ee200;
    ee292 = 2 * ee138;
    ee293 = 2 * (ee166 + ee26 * ee223);
    ee294 = 2 * (ee26 * ee156 + ee168);
    ee295 = 2 * ee174;
    ee296 = 2 * ee220;
    ee297 = 2 * ee222;
    ee298 = 2 * ee127;
    ee299 = 2 * ee101;
    ee300 = 2 * ee102;
    ee302 = 4 * ee201 + 8 * ee207;
    ee303 = 4 * ee21;
    ee304 = 4 * ee173;
    ee305 = 4 * ee55;
    ee306 = 4 * ee26;
    ee307 = 4 * ee127;
    ee308 = 4 * ee101;
    ee309 = 4 * ee102;
    ee310 = 4 * ee66;
    ee311 = 4 * ee67;
    ee312 = 6 * ee205;
    ee313 = 6 * ee113;
    ee314 = 6/ee4;
    ee315 = 8 * ee153;
    ee316 = 8 * ee48;
    ee317 = 8 * ee21;
    ee318 = 8 * ee127;
    ee319 = 8 * ee101;
    ee320 = 8 * ee102;
    ee321 = 8 * ee61;
    ee322 = 8 * ee62;

    out(j, 0) += -(((ee4 * (ee286 + ee79 - ee129) + (16 * ee61 +
      ee296 - (ee276/ee111 + 16 * ee62 + ee297)) * ee6/ee2 + ee312 +
      ee318 + ee319 - ee320) * ee1 - y * ((((2 * (ee85/ee70) -
      ee40 * ee1/ee5) * ee1 + ee177 + ee32) * ee26/ee70 - ee244) *
      ee1 + ee79)/ee30) * ee1/ee2);
    out(j, 1) += ((ee4 * (ee303 + ee306 + ee79 - ee266 * ee6/ee5) +
      (2 * ee215 + ee321 - ((ee124 + (ee213 + ee14 - (ee267 +
      ee210) * ee6) * ee26)/ee111 + 2 * ee217 + ee322)) * ee6/ee2 +
      4 * ee203 + ee307 + ee308 - ee309) * ee1 - y * (((ee124 +
      (2 * (ee162/ee70) + ee213) * ee26)/ee70 - ee189) * ee1 + ee32)/
        ee30) * ee1/ee2;
    out(j, 2) += ((((ee21 * (ee79 - ((2 + 2 * ee43 + ee11)/ee5 +
      ee210) * ee6) - ee133)/ee111 + 2 * ee226 + ee155 - (2 * ee225 +
      ee154)) * ee6/ee2 + (((28 - (ee273 * ee8 + ee123) * ee6/
        ee5) * ee1/ee2 + 4) * ee6/ee5 - (ee303 + ee79)) * ee4 + 2 *
          ee197 + ee300 - (ee298 + ee299)) * ee1 - y * (((ee133 - ((2 *
          (ee78/ee70) - 2 * (ee43 * ee1/ee5)) * ee1 + ee32) * ee21)/
            ee70 + ee251) * ee1 - ee14)/ee30) * ee1/ee2;
    out(j, 3) += -(ee1 * (y * ((ee110 * ee6/ee5 - (ee255 + ee212 -
      ee12) * ee21/ee30) * ee1/ee2 + ee12)/ee30 - ((((((6 - ee34) *
      ee1/ee5 + ee211) * ee6/ee2 + ee212 - ee32) * ee21/ee31 +
      ee6 * (ee66 - ee67))/ee2 + 3 * ee61 - 3 * ee62) * ee6/ee2 +
      (((ee4 * (16 - ee82) + ee25) * ee9/ee10 - 19) * ee9/ee10 +
      7/ee4) * ee4 + ee101 + pars1 - (ee102 + ee107 + pars2)) * ee1)/
        ee2);
    out(j, 4) += -(((ee4 * (ee114 + 32 * ee26 - ee263 * ee6/ee5) +
      (12 * ee220 + 2 * ((16 * ee66 + 2 * (4 * ee185 + ee310)) *
      ee6/ee2 + ee321) + 48 * ee61 - (((ee114 - (((10 * ee51 - 8 *
      (ee167/ee80))/ee4 + 2 * (ee85 + ee51/ee4)) * ee26/ee68 +
      2 * (ee98/ee5)) * ee6) * ee26 + ee168 + ee294)/ee111 + 12 *
      ee222 + 2 * ((16 * ee67 + 2 * (4 * ee187 + ee311)) * ee6/ee2 +
      ee322) + 48 * ee62)) * ee6/ee2 + 16 * ee127 + 16 * ee101 +
      24 * ee205 + 8 * (ee4 * ee156 - ee276 * ee6/ee68) - 16 * ee102) *
      ee1 - y * (((((((2 * (ee257 - ee281) + 4 * ee257 + 8 *
      ee281)/ee30 + 6 * ee51) * ee26/ee70 - 2 * ee244) * ee1 +
      ee114) * ee26 + ee168 + ee294)/ee70 - ee263 * ee1/ee5) * ee1 +
      ee114)/ee30) * ee1/ee2);
    out(j, 5) += ((ee4 * (ee286 + ee114 + ee317 - (((80 - ee271)/
      ee2 - ee234) * ee1 + ee8 * (24 - ((ee84 + ee284 + ee295 + ee89) *
        ee4 + ee72 + ee218) * ee6/ee5) + ee285 + 2) * ee6/ee5) +
        (2 * ((2 * (2 * ee185 + 2 * ee66) + 8 * ee66) * ee6/ee2 +
        ee154) + ee296 + 24 * ee61 + 8 * ee215 - ((ee275 + (ee177 +
        ee79 - (((((2 - ee183) * ee26 + ee305)/ee4 + 2 * (ee78 + ee55/
          ee4)) * ee26 + 6 * (ee124/ee4))/ee68 + ee266/ee5) * ee6) *
            ee26 + ee166 + ee293)/ee111 + 2 * ((2 * (2 * ee187 + 2 *
            ee67) + 8 * ee67) * ee6/ee2 + ee155) + ee297 + 24 * ee62 + 8 *
            ee217)) * ee6/ee2 + 12 * ee203 + 6 * (ee4 * ee223 - (ee124 +
            (ee213 - 2 * (ee163/ee80)) * ee26) * ee6/ee68) + ee312 +
            ee318 + ee319 - ee320) * ee1 + y * (((((ee189 + ((2 * (ee207 -
            ee201) - ee302) * ee26/ee30 - 6 * ee124)/ee70) * ee1 - ee32) *
            ee26 - (ee275 + ee166 + ee293))/ee70 + (((72 - ee271)/
              ee2 - ee234) * ee1 + ee8 * (12 - ((ee284 + ee295) * ee4 + ee218) *
                ee6/ee5) + 1 + ee90) * ee1/ee5) * ee1 - ee79)/ee30) *
                ee1/ee2;
    out(j, 6) += (((((ee242 + (96 - (ee273 * ee45 + ee274 - ee258) *
      ee6/ee5)/ee2) * ee1 + ee8 * (16 - (ee4 * (2 * ee87 + ee25 +
      ee77) + ee218) * ee6/ee5) + 4) * ee6/ee5 - (ee114 + ee306 +
      ee317)) * ee4 + ((ee21 * (ee114 - ((2 + 2 * ee122 + 2 *
      ee38 + ee57)/ee5 + 2 * (ee124/ee80)) * ee6) + ee26 * (ee14 +
      ee305 - (((ee26 * (4 - ee183) + 8 * ee55) * ee21/ee4 - 2 *
      (ee200 * ee26))/ee68 + ee267) * ee6) - (ee256 + ee289))/ee111 +
      12 * ee62 + 2 * ((2 * (ee187 + ee67) + ee311) * ee6/ee2 +
      ee179) + 4 * ee217 + 4 * ee226 - (12 * ee61 + 2 * ((2 * (ee185 +
      ee66) + ee310) * ee6/ee2 + ee178) + 4 * ee215 + 4 *
      ee225)) * ee6/ee2 + 4 * (ee194 * ee4 + (ee21 * (ee32 - (2 *
      (ee43/ee5) + ee210) * ee6) - ee133) * ee6/ee68) + 4 * ee197 +
      ee309 - (ee307 + ee308 + 8 * ee203)) * ee1 - y * ((((ee242 +
      (52 - (ee274 - (ee45 * ee125 + ee258)) * ee6/ee5)/ee2) *
      ee1 + 1 + ee69) * ee1/ee5 + (((2 * (ee195 * ee26) - ee21 * ee302)/
        ee30 - 4 * ee138) * ee26 * ee1/ee70 + ee256 + ee289 -
          ((2 * (ee124/ee70) - 2 * ee189) * ee1 + ee79) * ee21)/ee70) *
          ee1 - ee32)/ee30) * ee1/ee2;
    out(j, 7) += ((((((((ee246 * ee6/ee5 + 12) * ee1/ee2 + ee285 +
      ee313) * ee4 + ee248 + (ee4 * (6 - 24 * ee29) + 6 * ee8) *
      ee8 - ee268) * ee6/ee5 - 92) * ee1/ee2 - 8) * ee6/ee5 + ee114 +
      6 * ee21) * ee4 + ((((((ee21 * (6 - ee183)/ee4 - ee291) *
      ee26 + (ee292 - 4 * ee198)/ee4)/ee68 + (3 + 3 * ee141 +
      6 * ee43 + ee224)/ee5) * ee6 - 24/ee4) * ee21 + ee48 * (3 *
      ee55 + ee314 - 3 * ee23) - ee229)/ee111 + 2 * ((3 * ee66 + ee185) *
      ee6/ee2 + ee61) + 6 * ee225 + 6 * ee61 - (2 * ((3 *
      ee67 + ee187) * ee6/ee2 + ee62) + 6 * ee226 + 6 * ee62)) * ee6/
        ee2 + 2 * (ee117 * ee4 + (ee253 + ee212 - ee12) * ee21 *
          ee6/ee68) + ee298 + ee299 - (ee300 + 6 * ee197)) * ee1 - y *
          ((((((ee246 * ee9/ee10 + ee313) * ee4 + ee248 - ee268) * ee6/
            ee5 - 30) * ee1/ee2 - 1) * ee1/ee5 + (ee229 - ((((((ee290 -
              ee315) * ee26 + 2 * (ee198 * ee30))/ee30 + 2 * ee198 - ee292)/
                ee70 + 3 * ee251) * ee1 - ee314) * ee21 + 3 * (ee48 * ee55)))/
                  ee70) * ee1 + ee14)/ee30) * ee1/ee2;
    out(j, 8) += -(ee1 * (y * ((((((24 * ee4 + ee304 - ee282) *
      ee9/ee10 - 34) * ee4 - (14 * ee8 + ee76)) * ee9/ee10 + 15) *
      ee6/ee5 - ((ee255 - ee12) * ee48 + ee21 * (ee288 - ((ee290 +
      4 * ee160 - ee315)/ee30 + 4 * ee48) * ee21 * ee1/ee70) + ee287)/
        ee30) * ee1/ee2 - ee12)/ee30 - ((((((((ee21 * (ee183 -
          8) + ee316)/ee4 + ee291) * ee21/ee111 - (18 - ee82) * ee1/
            ee5) * ee6/ee2 + 10/ee4 + ee288 - ee316) * ee21 + (ee253 - ee12) *
              ee48 + ee287)/ee31 + (6 * ee67 + ee6 * (ee158 - ee157)/
                ee2 - 6 * ee66) * ee6)/ee2 + 7 * ee62 - 7 * ee61) * ee6/ee2 +
                  (((((ee304 + 56 * ee4 - ee282) * ee9/ee10 - 86) * ee4 - (ee76 +
                  22 * ee8)) * ee9/ee10 + 65) * ee9/ee10 - 15/ee4) * ee4 +
                  ee102 + ee107 + pars2 - (ee101 + pars1)) * ee1)/ee2);

        }

  }

  return out;
  
}

// //' @param pars a list of vectors of coefficients for negbinson log mean
// //' @param X1 a sparse design matrix for the negbinson log mean parameter
// //' @param ymat a matrix
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return negbind0 a scalar, the negative log-likelihood
// //' @return negbind12 a matrix, first then second derivatives w.r.t. negbinson parameters
// //' @return negbind34 a matrix, third then fourth derivatives w.r.t. negbinson parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double negbinspd0(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::mat ymat, 
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
    // size = exp(pars2);
    // size = exp(-pars2);
    sigsq = mu + exp(pars2);
    size = mu * mu / (sigsq - mu);
    p = size / (size + mu);
    nllh -= lgamma(y + size) - lgamma(size) - lgamma(y + 1) + size * log(p) + y * log(1 - p);
    
  } 
  
}

return(nllh);

}

// //' @rdname negbinspd0
// [[Rcpp::export]]
arma::mat negbinspd12(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::mat ymat, 
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

  // double ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  // double ee10, ee11, ee12;
  // 
  double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee19;
  double ee20, ee21, ee22, ee23, ee27, ee29;
  double ee30, ee32, ee33, ee34, ee35, ee38;
  double ee41, ee42, ee43, ee44, ee45, ee46, ee47;

  arma::mat out = arma::mat(nobs, 5, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      w = wmat(j, l);
      
      // ee2 = exp(-pars2);
      // ee3 = exp(pars1);
      // ee4 = ee2 + ee3;
      // ee5 = ee2/ee4;
      // ee6 = 1 - ee5;
      // ee7 = ee2 + y;
      // ee8 = R_pow(ee5, 2);
      // ee9 = ee6 * ee4;
      // ee10 = Rf_digamma(ee7);
      // ee11 = Rf_digamma(ee2);
      // ee12 = log(ee4);
      // 
      // out(j, 0) += w * (-(ee2 * ee3 * (y/ee9 - 1)/ee4));
      // out(j, 1) += w * (-((ee7/ee4 + ee11 + ee12 + pars2 - (1 + ee10)) *
      //   ee2));
      // out(j, 2) += w * (-(((y * (1 - (2 + ee2/ee9) * ee3/ee4)/ee6 -
      //   (R_pow(ee2, 2)/(R_pow(ee4, 2) * ee8) - 2) * ee3)/ee4 - 1) *
      //   ee2 * ee3/ee4));
      // out(j, 3) += w * ((((ee6 * ee2/(ee4 * ee8) + 2) * ee2 + y)/
      //   ee4 - 2) * ee2 * ee3/ee4);
      // out(j, 4) += w * ((ee11 + ee12 + pars2 + y * ((1 - (3 - 2 *
      //   ee5) * ee2/ee4)/ee6 + ee5)/ee4 - (3 + ee10 + ee2 * (Rf_trigamma(ee7) -
      //   (((R_pow(ee6, 2)/ee8 - 2) * ee2/ee4 + 5)/ee4 + Rf_trigamma(ee2))))) *
      //   ee2);

      ee1 = exp(pars1);
      ee2 = exp(pars2);
      ee3 = ee1/ee2;
      ee4 = 1 + ee3;
      ee5 = R_pow(ee1, 2);
      ee6 = R_pow((ee4 * ee1), 2);
      ee7 = ee5/ee2;
      ee8 = 1 + 2 * ee3;
      ee10 = ee1/(ee4 * ee2);
      ee11 = ee6 * ee2;
      ee12 = R_pow(ee1, 3);
      ee13 = 1 - ee10;
      ee14 = 1/ee4;
      ee15 = 2/ee4;
      ee16 = ee7 + y;
      ee19 = ee14 - ee12/ee11;
      ee20 = ee15 - ee8 * ee5/ee6;
      ee21 = Rf_digamma(ee16);
      ee22 = Rf_digamma(ee7);
      ee23 = log1p(ee3);
      ee27 = 2 * (ee8 * ee4 * ee5/ee6);
      ee29 = R_pow(ee10, 2) * ee2;
      ee30 = 4/ee4;
      ee32 = pars1 - (ee23 + pars2);
      ee33 = Rf_trigamma(ee16);
      ee34 = Rf_trigamma(ee7);
      ee35 = ee13 * ee2;
      ee38 = ee8 * (4 - ee27) + 1 + 4 * ee3;
      ee41 = ee19 * ee20;
      ee42 = R_pow(ee19, 2);
      ee43 = R_pow(ee20, 2);
      ee44 = 2 * (ee4 * ee12/ee11);
      ee45 = 2 * ee32;
      ee46 = 2 * ee21;
      ee47 = 2 * ee22;
      
      out(j, 0) += w * (-(((ee4 * ee20 + ee45 + ee46 - ee47) * ee1 -
        y * ee20/ee13) * ee1/ee2));
      out(j, 1) += w * (-(ee1 * (y * ee19/ee13 - (ee4 * ee19 + ee21 +
        pars1 - (ee22 + ee23 + pars2)) * ee1)/ee2));
      out(j, 2) += w * (-(((ee4 * (4 * ee20 + ee30 - ee38 * ee5/ee6) +
        (4 * ee33 - (ee43/ee29 + 4 * ee34)) * ee5/ee2 + 4 * ee32 +
        4 * ee21 - 4 * ee22) * ee1 - y * ((ee43/ee35 - ee38 * ee1/
          ee6) * ee1 + ee30)/ee13) * ee1/ee2));
      out(j, 3) += w * (((ee4 * (2 * ee19 + ee30 - ((8 - ee27) * ee1/
        ee2 + 2) * ee5/ee6) + (2 * ee33 - (ee41/ee29 + 2 * ee34)) *
          ee5/ee2 + ee45 + ee46 - ee47) * ee1 - y * ((ee41/ee35 - ((6 -
          ee27) * ee1/ee2 + 1) * ee1/ee6) * ee1 + ee15)/ee13) * ee1/
            ee2);
      out(j, 4) += w * (-(ee1 * (y * (((3 - ee44) * ee5/ee6 - ee42/
        ee13) * ee1/ee2 - ee14)/ee13 - ((ee42/ee29 + ee34 - ee33) *
          ee5/ee2 + ((5 - ee44) * ee12/ee11 - 3/ee4) * ee4 + ee22 + ee23 +
          pars2 - (ee21 + pars1)) * ee1)/ee2));
    }
    
  }
  
  return out;
  
}

// //' @rdname negbinspd0
// [[Rcpp::export]]
arma::mat negbinspd34(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::mat ymat, 
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

  // double ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  // double ee10, ee11, ee12, ee14, ee15, ee16, ee17, ee19;
  // double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
  // double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee38, ee39;
  // double ee40, ee42, ee44, ee45, ee47, ee48, ee49;
  // double ee51, ee52, ee54, ee55, ee56, ee59;
  // double ee62, ee63, ee64, ee65, ee66, ee67, ee68;
  // double ee70, ee72, ee77;
  // double ee81, ee82, ee83, ee84, ee85, ee86, ee87, ee89;
  // double ee90, ee91;
  // double ee103, ee104, ee106, ee107, ee108, ee109;
  // double ee110, ee111, ee112, ee113, ee114, ee115, ee116, ee117, ee118, ee119;
  // double ee121, ee122, ee123, ee124, ee125, ee126, ee129;
  // double ee131, ee132, ee133, ee134, ee135;
  
  double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee17, ee19;
  double ee20, ee21, ee23, ee24, ee25, ee26, ee27, ee29;
  double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee37, ee38;
  double ee40, ee43, ee45, ee48;
  double ee51, ee52, ee55, ee56, ee57, ee58, ee59;
  double ee60, ee61, ee62, ee63, ee65, ee66, ee67, ee68, ee69;
  double ee70, ee71, ee72, ee73, ee74, ee75, ee76, ee77, ee78, ee79;
  double ee80, ee81, ee82, ee83, ee84, ee85, ee87, ee88, ee89;
  double ee90, ee98, ee99;
  double ee100, ee101, ee102, ee104, ee106, ee107;
  double ee110, ee111, ee112, ee113, ee114, ee117;
  double ee122, ee123, ee124, ee125, ee127, ee129;
  double ee133, ee138;
  double ee141, ee143, ee147, ee148, ee149;
  double ee150, ee151, ee153, ee154, ee155, ee156, ee157, ee158;
  double ee160, ee161, ee162, ee163, ee166, ee167, ee168;
  double ee172, ee173, ee174, ee175, ee176, ee177, ee178, ee179;
  double ee180, ee181, ee182, ee183, ee185, ee187, ee189;
  double ee194, ee195, ee197, ee198;
  double ee200, ee201, ee203, ee205, ee207;
  double ee210, ee211, ee212, ee213, ee215, ee217, ee218;
  double ee220, ee222, ee223, ee224, ee225, ee226, ee229;
  double ee234;
  double ee242, ee244, ee246, ee248;
  double ee251, ee253, ee255, ee256, ee257, ee258;
  double ee263, ee266, ee267, ee268;
  double ee271, ee273, ee274, ee275, ee276;
  double ee281, ee282, ee284, ee285, ee286, ee287, ee288, ee289;
  double ee290, ee291, ee292, ee293, ee294, ee295, ee296, ee297, ee298, ee299;
  double ee300, ee302, ee303, ee304, ee305, ee306, ee307, ee308, ee309;
  double ee310, ee311, ee312, ee313, ee314, ee315, ee316, ee317, ee318, ee319;
  double ee320, ee321, ee322;

  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      w = wmat(j, l);
      
      // ee2 = exp(-pars2);
      // ee3 = exp(pars1);
      // ee4 = ee2 + ee3;
      // ee5 = ee2/ee4;
      // ee6 = 1 - ee5;
      // ee7 = 2 * ee5;
      // ee8 = 2 * ee2;
      // ee9 = 2 * ee4;
      // ee10 = ee8 + ee3;
      // ee11 = R_pow(ee5, 2);
      // ee12 = 2 * ee10;
      // ee14 = (3 - ee7) * ee2/ee4;
      // ee15 = ee3/ee4;
      // ee16 = 1 - 2 * ee15;
      // ee17 = 1 - ee7;
      // ee19 = 2 * ee3 + ee2;
      // ee20 = 2 * ee19;
      // ee21 = 4 * ee4;
      // ee22 = 1 - ee14;
      // ee23 = R_pow(ee2, 2);
      // ee24 = R_pow(ee4, 2);
      // ee25 = 8 * ee2;
      // ee26 = ee12 + ee21;
      // ee27 = ee20 + ee9;
      // ee28 = ee24 * ee11;
      // ee29 = 4 * ee2;
      // ee30 = ee12 + ee9;
      // ee31 = 6 * ee2;
      // ee32 = ee4 * ee11;
      // ee33 = ee9 - ee31;
      // ee34 = 8 * ee3;
      // ee35 = R_pow(ee6, 2);
      // ee36 = (ee26 - ee25)/ee4;
      // ee38 = ee33 * ee3/ee4;
      // ee39 = ee2 + y;
      // ee40 = (ee27 - ee34)/ee4;
      // ee42 = 4 * (ee12 + ee29);
      // ee44 = (ee36 + 6) * ee2/ee4;
      // ee45 = ee17 * ee6;
      // ee47 = ee6 * ee23/ee28;
      // ee48 = ee40 + 2;
      // ee49 = 3 * ee5;
      // ee51 = ((ee30 - ee25)/ee4 + 2) * ee2/ee4;
      // ee52 = (ee38 + ee8)/ee4;
      // ee54 = ee45 + 1 - ee14;
      // ee55 = ee16 * ee6;
      // ee56 = ee6 * ee4;
      // ee59 = (2 * (ee35/ee11) - 2) * ee2/ee4;
      // ee62 = 1 - (7 - ee44) * ee2/ee4;
      // ee63 = 2 * ee16;
      // ee64 = 2 * (ee29 + ee3);
      // ee65 = ee9 + ee8;
      // ee66 = 4 * ee17;
      // ee67 = 4 * ee3;
      // ee68 = 8 * ee47;
      // ee70 = ee48 * ee3/ee4;
      // ee72 = R_pow(ee22, 2) + ee62 * ee6;
      // ee77 = (2 * (ee23/ee28) - 2) * ee3/ee4 + 1;
      // ee81 = 2 * (ee35 + 1 - ee14);
      // ee82 = 2 * ((4 - ee49) * ee2/ee4 - 1);
      // ee83 = 2 * ee22;
      // ee84 = 2 * ee17;
      // ee85 = 64 * ee2;
      // ee86 = 8 * ee5;
      // ee87 = ee72/ee6;
      // ee89 = (ee59 + 3) * ee2/ee4;
      // ee90 = ((4 * (ee20 + ee67) + 8 * ee27 - 64 * ee3)/ee4 +  8) * ee3;
      // ee91 = (1 - ee52) * ee6;
      // ee103 = ee36 + 12;
      // ee104 = ee77 + ee63;
      // ee106 = (2 + 2 * (ee2/ee56)) * ee3/ee4;
      // ee107 = 1 - ee51;
      // ee108 = 2 * (ee17 * ee16 + 1 - ee52);
      // ee109 = 2 * ee54;
      // ee110 = 2 * (R_pow(ee17, 2) + 1 - ee51);
      // ee111 = 2 * (R_pow(ee16, 2) + 1 - ee70);
      // ee112 = 2 * (ee6 * ee2/ee32);
      // ee113 = 2 * ee48;
      // ee114 = 2 * (ee67 + ee2);
      // ee115 = 4 * ee55;
      // ee116 = ee66 + ee86;
      // ee117 = 4 * ee16;
      // ee118 = ee42 + 8 * ee30;
      // ee119 = ee42 + 8 * ee26;
      // ee121 = 4 * ee65 - 40 * ee2;
      // ee122 = 6 * ee19;
      // ee123 = 6 * ee4;
      // ee124 = 8 * ee22;
      // ee125 = 8 * ee65;
      // ee126 = Rf_digamma(ee39);
      // ee129 = Rf_digamma(ee2) + log(ee4) + pars2;
      // ee131 = ee2 * ee3/ee24;
      // ee132 = Rf_psigamma(ee39, 2);
      // ee133 = Rf_psigamma(ee2, 2);
      // ee134 = Rf_trigamma(ee39);
      // ee135 = Rf_trigamma(ee2);
      // 
      // out(j, 0) += -(((y * (1 - (((1 + ee63 - ee106) * ee2/ee6 + ee20 +
      //   ee9 - ee34)/ee4 + 2) * ee3/ee4)/ee6 - ((ee104 * ee23/
      //     ee32 + ee34 - ee27)/ee4 - 2) * ee3)/ee4 - 1) * ee2 * ee3/ee4);
      // out(j, 1) += ((y * (1 - ((((ee84 + ee7)/ee6 + 2) * ee3/ee4 +
      //   1) * ee2 + ee38)/ee4)/ee6 - ((((ee84 - 2 * ee47) * ee3/ee4 -
      //   ee55) * ee2/ee32 - 2) * ee2 + (((6 + ee2/ee32) * ee2 - ee9)/
      //     ee4 - 2) * ee3))/ee4 - 2) * ee2 * ee3/ee4;
      // out(j, 2) += -(((((((ee6 * (2 - (ee112 + 4) * ee2/ee4) + 3 -
      //   (5 - ee7) * ee2/ee4)/ee11 - 8) * ee2 + ee12 + ee9)/ee4 + 6) *
      //   ee2 + y * ((ee54/ee6 + (ee25 - ee30)/ee4 - 1) * ee2/ee4 +
      //   1)/ee6)/ee4 - 4) * ee2 * ee3/ee4);
      // out(j, 3) += -((ee129 + y * (((((ee26 - ee31)/ee4 + 3) * ee2/
      //   ee4 - 6) * ee2/ee4 + 1)/ee6 + (2 * (ee22/ee6) + ee7) * ee2/
      //     ee4)/ee4 - ((3 * ee134 + ee2 * (ee132 - ee133) - (((ee6 * (ee83 +
      //       4 - (ee59 + 6) * ee2/ee4)/ee11 - ee103) * ee2/ee4 + 19)/
      //         ee4 + 3 * ee135)) * ee2 + 7 + ee126)) * ee2);
      // out(j, 4) += -(((y * (1 - ((((1 - ee106) * ee16 + 2 + ee111 -
      //   (((2 * (ee55 + ee131) + ee115 - 8 * ee131)/ee6 + ee117) *
      //   ee2/ee56 + ee113) * ee3/ee4) * ee2/ee6 + ee114 + ee21 + ee122 -
      //   ee90)/ee4 + 2) * ee3/ee4)/ee6 - (((((2 * (1 - 3 * ee15) +
      //   8 * ee16 + 8 * (ee23 * ee3/(R_pow(ee4, 3) * ee11))) * ee23/
      //     ee28 - ee113) * ee3/ee4 + ee77 * ee16 + 2 + ee111) * ee23/
      //       ee32 + ee90 - (ee114 + ee21 + ee122))/ee4 - 2) * ee3)/ee4 - 1) *
      //         ee2 * ee3/ee4);
      // out(j, 5) += ((y * (1 - (((((((2 * (ee49 - 1) - ee116) * ee3/
      //   ee56 + ee117 - 2) * ee2 - (ee33/ee4 + 2) * ee3)/ee4 + 2 + ee108)/
      //     ee6 + ee40 + 2) * ee3/ee4 + 1) * ee2 + (ee20 + ee21 -
      //       ((2 + 4 * (ee27/ee4)) * ee2 + ee121 * ee3/ee4)) * ee3/ee4)/
      //         ee4)/ee6 - ((((((((2 * (2 - ee49) + ee66 - ee68) * ee3/ee4 -
      //           ee115) * ee2/ee32 - 2) * ee2 - ee38)/ee4 + (1 - (2 + ee112) *
      //           ee2/ee4) * ee16 + 1 + ee108) * ee3/ee4 - (1 - ee70) * ee6) *
      //           ee2/ee32 - 2) * ee2 + ((((ee104 * ee2/ee11 + 4 * ee27)/ee4 +
      //           2) * ee2 + (ee121/ee4 + 8) * ee3 - (4 * ee19 + ee123))/
      //             ee4 - 2) * ee3))/ee4 - 2) * ee2 * ee3/ee4;
      // out(j, 6) += -(((y * (((((((ee82 - 4 * ee45)/ee6 - ee116) *
      //   ee2/ee4 - ee110) * ee3/ee4 + ee91 + ee22 * ee16)/ee6 + ((ee63 +
      //   6) * ee2 - (ee38 + ee12 + ee9))/ee4 - 1) * ee2 - ((6 - (ee42 +
      //   ee125 - ee85)/ee4) * ee2 + ee9) * ee3/ee4)/ee4 + 1)/
      //     ee6 - ((((((1/ee11 - 64) * ee2 + ee42 + ee125)/ee4 + 6) * ee2 -
      //       ee123)/ee4 - 2) * ee3 + (((((ee110 + ee66 - (ee6 * (4 + 8 *
      //       ee17 - ee68) + ee81) * ee23/ee28) * ee3/ee4 - (ee16 * (3 -
      //       (ee59 + 5) * ee2/ee4) + 2 * ee91))/ee11 + 8) * ee2 - ee30)/
      //         ee4 - 6) * ee2))/ee4 - 4) * ee2 * ee3/ee4);
      // out(j, 7) += -(((y * (((((ee82 - ee109)/ee6 + (ee30 + 48 * ee2 -
      //   ee118)/ee4 - (ee84 + 6)) * ee2 + ee64 + ee21 + 6 * ee10)/
      //     ee4 + 1 - ((ee44 + ee109 - 7) * ee2/ee4 + 1 + 2 * (ee107 *
      //       ee6) + 3 * (ee22 * ee17))/ee6) * ee2/ee4 - 1)/ee6 - ((((ee17 *
      //       (1 + ee83 - ee89) + ee6 * (2 * ee107 + 7 - ((((ee6 * (6 -
      //       ee68) + ee81 + 4 * ee54)/ee11 - 8) * ee2 + ee12 + ee9)/ee4 +
      //       14) * ee2/ee4) + 7 - (19 - ee103 * ee2/ee4) * ee2/ee4)/ee11 -
      //       ((ee118 - ee85)/ee4 + 32)) * ee2 + 10 * ee4 + 12 * ee10 +
      //       ee64)/ee4 + 14) * ee2)/ee4 + 8) * ee2 * ee3/ee4);
      // out(j, 8) += -((((6 * ee132 + ee2 * (Rf_psigamma(ee39, 3) - Rf_psigamma(ee2, 3)) -
      //   6 * ee133) * ee2 + 7 * ee134 - (((((1 - ee89) *
      //   ee22 + ee6 * (10 + 2 * ee62 + ee124 - (((ee6 * (8 - ee68) +
      //   ee81 + ee124) * ee6/ee11 - 8) * ee2/ee4 + 18) * ee2/ee4) +
      //   2 * ee72)/ee11 - ((ee64 + 22 * ee10 + 36 * ee4 - ((ee119 -
      //   ee85)/ee4 + 56) * ee2)/ee4 + 50)) * ee2/ee4 + 65)/ee4 + 7 *
      //   ee135)) * ee2 + 15 + ee126 + y * ((((ee82 - 4 * ee22)/ee6 -
      //   ee86) * ee2/ee4 - (2 * ee87 + 4 * (ee22 * ee2/ee4))/ee6) *
      //   ee2/ee4 - ((((((ee26 + 56 * ee2 - ee119)/ee4 - 18) * ee2 +
      //   14 * ee10 + ee64 + 20 * ee4)/ee4 + ee83 + 7) * ee2/ee4 + ee87 -
      //   14) * ee2/ee4 + 1)/ee6)/ee4 - ee129) * ee2);
      // 
      ee1 = exp(pars1);
      ee2 = exp(pars2);
      ee3 = ee1/ee2;
      ee4 = 1 + ee3;
      ee5 = R_pow((ee4 * ee1), 2);
      ee6 = R_pow(ee1, 2);
      ee7 = 2 * ee3;
      ee8 = 1 + ee7;
      ee9 = R_pow(ee1, 3);
      ee10 = ee5 * ee2;
      ee11 = 4 * ee3;
      ee12 = 1/ee4;
      ee13 = ee6/ee2;
      ee14 = 2/ee4;
      ee17 = ee8 * ee4 * ee6/ee5;
      ee19 = ee1/(ee4 * ee2);
      ee20 = 1 + ee11;
      ee21 = ee12 - ee9/ee10;
      ee23 = ee8 * ee6/ee5;
      ee24 = 2 * ee17;
      ee25 = 2 * ee8;
      ee26 = ee14 - ee23;
      ee27 = ee13 + y;
      ee29 = ee4 * ee9/ee10;
      ee30 = 1 - ee19;
      ee31 = R_pow(ee19, 2);
      ee32 = 4/ee4;
      ee33 = R_pow(ee4, 2);
      ee34 = 2 * ee29;
      ee35 = R_pow(ee8, 2);
      ee36 = ee20 * ee4;
      ee37 = 3 - ee34;
      ee38 = ee8 * (4 - ee24);
      ee40 = ee38 + 1 + ee11;
      ee43 = (6 - ee24) * ee1/ee2 + 1;
      ee45 = 1 + ee25 + ee11;
      ee48 = ee37 * ee9/ee10 - ee12;
      ee51 = ee32 - ee40 * ee6/ee5;
      ee52 = R_pow(ee2, 2);
      ee55 = ee14 - ee43 * ee6/ee5;
      ee56 = 1 + 2 * ee4;
      ee57 = 8 * ee3;
      ee58 = ee35 + ee36;
      ee59 = ee56 + ee7;
      ee60 = ee33 * ee9;
      ee61 = Rf_trigamma(ee27);
      ee62 = Rf_trigamma(ee13);
      ee63 = 4 - 8 * ee17;
      ee65 = ee59 * ee1/ee2;
      ee66 = Rf_psigamma(ee27, 2);
      ee67 = Rf_psigamma(ee13, 2);
      ee68 = ee31 * ee52;
      ee69 = 4 * ee8;
      ee70 = ee30 * ee2;
      ee71 = R_pow(ee21, 2);
      ee72 = 2 * ee58;
      ee73 = ee45 * ee4;
      ee74 = ee4 * ee31;
      ee75 = 2 * ee65;
      ee76 = 2 * ee20;
      ee77 = 4 * ee20;
      ee78 = ee21 * ee26;
      ee79 = 8/ee4;
      ee80 = ee74 * ee52;
      ee81 = ee60/ee2;
      ee82 = 8 * ee29;
      ee83 = ee60/ee10;
      ee84 = ee8 * ee63;
      ee85 = R_pow(ee26, 2);
      ee87 = ee63 * ee1/ee2;
      ee88 = 2 * ee36;
      ee89 = 4 * ee45;
      ee90 = 6 * ee20;
      ee98 = ee8 * (12 - ((ee84 + ee89) * ee4 + ee72) * ee6/ee5) +  1 + ee90 + ee57;
      ee99 = 2 * (ee73 * ee6/ee5);
      ee100 = 8 * ee83;
      ee101 = Rf_digamma(ee27);
      ee102 = Rf_digamma(ee13);
      ee104 = ee6 * ee66/ee2;
      ee106 = ee6 * ee67/ee2;
      ee107 = log1p(ee3);
      ee110 = (ee4 * (10 - ee82) + ee25) * ee9/ee10 - 7;
      ee111 = ee31 * ee2;
      ee112 = 1 + ee69;
      ee113 = 1 + ee57;
      ee114 = 16/ee4;
      ee117 = ee110 * ee9/ee10 + ee12;
      ee122 = ee8 * (4 - ((ee87 + ee76) * ee4 + ee75) * ee6/ee5) +  (20 - ee99) * ee1/ee2 + 1;
      ee123 = ee4 * (ee77 + ee11);
      ee124 = ee21 * ee51;
      ee125 = ee100 - ee25;
      ee127 = pars1 - (ee107 + pars2);
      ee129 = ee98 * ee6/ee5;
      ee133 = ee48 * ee26;
      ee138 = ee21 * ee55;
      ee141 = (14 - (ee123 - ee8 * ee125) * ee6/ee5) * ee1/ee2 +  1;
      ee143 = ee112 + ee77 + ee57;
      ee147 = 2 * ee73 + ee72;
      ee148 = ee75 + ee88;
      ee149 = 4 * ee4;
      ee150 = 4 * (2 * (ee5 * ee8) + 4 * ee81);
      ee151 = 64 * ee81;
      ee153 = ee71 * ee1/ee2;
      ee154 = 4 * ee61;
      ee155 = 4 * ee62;
      ee156 = ee79 - ee129;
      ee157 = Rf_psigamma(ee27, 3);
      ee158 = Rf_psigamma(ee13, 3);
      ee160 = ee48 * ee30;
      ee161 = ee8 * (ee151 - ee150);
      ee162 = ee78 * ee1;
      ee163 = ee78 * ee6;
      ee166 = ee55 * ee51;
      ee167 = ee85 * ee6;
      ee168 = R_pow(ee51, 2);
      ee172 = ee88 + 2 * ee35;
      ee173 = ee25 + ee149;
      ee174 = ee25 + ee11;
      ee175 = 2 * ee143;
      ee176 = ee76 + ee69;
      ee177 = 2 * ee51;
      ee178 = 2 * ee61;
      ee179 = 2 * ee62;
      ee180 = 4 * (ee5 * ee148);
      ee181 = 4 * (1 + 12 * ee3 + ee25);
      ee182 = 8 - 16 * ee17;
      ee183 = 8 * (ee21 * ee6/ee80);
      ee185 = ee6 * ee157/ee2;
      ee187 = ee6 * ee158/ee2;
      ee189 = ee122 * ee1/ee5;
      ee194 = ee141 * ee6/ee5 - ee14;
      ee195 = ee160 + ee153;
      ee197 = ee48 * ee4 + ee71 * ee6/ee68;
      ee198 = ee133 - ee138;
      ee200 = ee48/ee4 - ee71;
      ee201 = ee30 * ee55;
      ee203 = ee4 * ee55 - ee163/ee68;
      ee205 = ee4 * ee51 - ee167/ee68;
      ee207 = ee162/ee2;
      ee210 = 2 * (ee78/ee80);
      ee211 = 2 * (ee71/(ee74 * ee2));
      ee212 = 2 * ee48;
      ee213 = 2 * ee55;
      ee215 = 2 * ee104 + ee178;
      ee217 = 2 * ee106 + ee179;
      ee218 = 4 * ee65;
      ee220 = 4 * ee104 + ee154;
      ee222 = 4 * ee106 + ee155;
      ee223 = ee32 - ee122 * ee6/ee5;
      ee224 = 6 * ee3;
      ee225 = ee104 + ee61;
      ee226 = ee106 + ee62;
      ee229 = ee117 * ee26;
      ee234 = ((ee4 * (ee181 - (((4 * (2 * (ee5 * ee59) + 4 *  (ee8 * ee33 * ee6)) - 64 * (ee8 * ee33 * ee6)) * ee1/ee2 +  ee180) * ee8/ee5 + 4 * (ee147 * ee1/ee2)) * ee6/ee5) +  2 * ((ee112 + ee149 + ee11) * ee1/ee2)) * ee8 + ee45 *  (4 * ee59 - 8 * (ee8 * ee33 * ee6/ee5)) * ee1/ee2 + 2 *  (ee58 * ee20)) * ee1/ee5;
      ee242 = (((4 * ee148 - (ee161 - ee180)/ee5) * ee9/ee10 -  2 * ee113) * ee8 * ee4 - (2 * ((ee56 + ee224) * ee8) +  4 * (ee59 * ee20)) * ee1/ee2) * ee1/ee5;
      ee244 = ee98 * ee1/ee5;
      ee246 = (ee161 - 4 * (ee5 * ee172))/ee5 - 4 * ee172;
      ee248 = (ee4 * (20 - 16 * ee29) + ee69) * ee1/ee2;
      ee251 = ee141 * ee1/ee5;
      ee253 = (ee37 * ee1/ee5 + ee211) * ee6/ee2;
      ee255 = (ee37 * ee6/ee5 - 2 * (ee71/ee30)) * ee1/ee2;
      ee256 = ee48 * ee51;
      ee257 = ee30 * ee51;
      ee258 = ee8 * (16 * ee83 - ee69);
      ee263 = ee8 * (32 - ((ee8 * ee182 + ee175 + 4 * ee176) *  ee4 + 4 * ee58) * ee6/ee5) + (16/ee2 - ((ee4 * (4 * ee143 -  ((4 * (ee5 * ee147) + 4 * (2 * (ee58 * ee5) + 4 * (ee35 *  ee33 * ee6)) - 64 * (ee35 * ee33 * ee6))/ee5 + 4 *  ee147) * ee8 * ee6/ee5) + 2 * (ee113 * ee4 + 3 *  (ee8 * ee20))) * ee8 + ee45 * (6 * ee58 - 8 * (ee35 *  ee33 * ee6/ee5))) * ee1/ee5) * ee1 + 1 + 24 * ee20 +  8 * ee113;
      ee266 = ee8 * (8 - ((ee87 + ee25 + ee76) * ee4 + ee75) *  ee6/ee5) + (24 - ee99) * ee1/ee2 + 2;
      ee267 = ee8/ee5;
      ee268 = ee20 * (ee100 - 8 * ee8);
      ee271 = ee4 * (ee175 + 2 * ee176) * ee6/ee5;
      ee273 = ee4 * (4 - ee82) + ee25;
      ee274 = ee4 * (ee181 + 4 * ee174 + ee57);
      ee275 = ee21 * ee156;
      ee276 = (ee177 + ee32 - (ee40/ee5 + 2 * (ee85/ee80)) * ee6) *  ee26;
      ee281 = ee85 * ee1/ee2;
      ee282 = (ee151 - (4 * (ee5 * ee173) + ee150))/ee5;
      ee284 = ee182 * ee1/ee2;
      ee285 = 12 * ee20;
      ee286 = 12 * ee26;
      ee287 = 2 * (ee117 * ee21 + R_pow(ee48, 2));
      ee288 = 2 * ee117;
      ee289 = 2 * (ee194 * ee26 - R_pow(ee55, 2));
      ee290 = 2 * ee195;
      ee291 = 2 * ee200;
      ee292 = 2 * ee138;
      ee293 = 2 * (ee166 + ee26 * ee223);
      ee294 = 2 * (ee26 * ee156 + ee168);
      ee295 = 2 * ee174;
      ee296 = 2 * ee220;
      ee297 = 2 * ee222;
      ee298 = 2 * ee127;
      ee299 = 2 * ee101;
      ee300 = 2 * ee102;
      ee302 = 4 * ee201 + 8 * ee207;
      ee303 = 4 * ee21;
      ee304 = 4 * ee173;
      ee305 = 4 * ee55;
      ee306 = 4 * ee26;
      ee307 = 4 * ee127;
      ee308 = 4 * ee101;
      ee309 = 4 * ee102;
      ee310 = 4 * ee66;
      ee311 = 4 * ee67;
      ee312 = 6 * ee205;
      ee313 = 6 * ee113;
      ee314 = 6/ee4;
      ee315 = 8 * ee153;
      ee316 = 8 * ee48;
      ee317 = 8 * ee21;
      ee318 = 8 * ee127;
      ee319 = 8 * ee101;
      ee320 = 8 * ee102;
      ee321 = 8 * ee61;
      ee322 = 8 * ee62;
      
      out(j, 0) += -(((ee4 * (ee286 + ee79 - ee129) + (16 * ee61 +
        ee296 - (ee276/ee111 + 16 * ee62 + ee297)) * ee6/ee2 + ee312 +
        ee318 + ee319 - ee320) * ee1 - y * ((((2 * (ee85/ee70) -
        ee40 * ee1/ee5) * ee1 + ee177 + ee32) * ee26/ee70 - ee244) *
        ee1 + ee79)/ee30) * ee1/ee2);
      out(j, 1) += ((ee4 * (ee303 + ee306 + ee79 - ee266 * ee6/ee5) +
        (2 * ee215 + ee321 - ((ee124 + (ee213 + ee14 - (ee267 +
        ee210) * ee6) * ee26)/ee111 + 2 * ee217 + ee322)) * ee6/ee2 +
        4 * ee203 + ee307 + ee308 - ee309) * ee1 - y * (((ee124 +
        (2 * (ee162/ee70) + ee213) * ee26)/ee70 - ee189) * ee1 + ee32)/
          ee30) * ee1/ee2;
      out(j, 2) += ((((ee21 * (ee79 - ((2 + 2 * ee43 + ee11)/ee5 +
        ee210) * ee6) - ee133)/ee111 + 2 * ee226 + ee155 - (2 * ee225 +
        ee154)) * ee6/ee2 + (((28 - (ee273 * ee8 + ee123) * ee6/
          ee5) * ee1/ee2 + 4) * ee6/ee5 - (ee303 + ee79)) * ee4 + 2 *
            ee197 + ee300 - (ee298 + ee299)) * ee1 - y * (((ee133 - ((2 *
            (ee78/ee70) - 2 * (ee43 * ee1/ee5)) * ee1 + ee32) * ee21)/
              ee70 + ee251) * ee1 - ee14)/ee30) * ee1/ee2;
      out(j, 3) += -(ee1 * (y * ((ee110 * ee6/ee5 - (ee255 + ee212 -
        ee12) * ee21/ee30) * ee1/ee2 + ee12)/ee30 - ((((((6 - ee34) *
        ee1/ee5 + ee211) * ee6/ee2 + ee212 - ee32) * ee21/ee31 +
        ee6 * (ee66 - ee67))/ee2 + 3 * ee61 - 3 * ee62) * ee6/ee2 +
        (((ee4 * (16 - ee82) + ee25) * ee9/ee10 - 19) * ee9/ee10 +
        7/ee4) * ee4 + ee101 + pars1 - (ee102 + ee107 + pars2)) * ee1)/
          ee2);
      out(j, 4) += -(((ee4 * (ee114 + 32 * ee26 - ee263 * ee6/ee5) +
        (12 * ee220 + 2 * ((16 * ee66 + 2 * (4 * ee185 + ee310)) *
        ee6/ee2 + ee321) + 48 * ee61 - (((ee114 - (((10 * ee51 - 8 *
        (ee167/ee80))/ee4 + 2 * (ee85 + ee51/ee4)) * ee26/ee68 +
        2 * (ee98/ee5)) * ee6) * ee26 + ee168 + ee294)/ee111 + 12 *
        ee222 + 2 * ((16 * ee67 + 2 * (4 * ee187 + ee311)) * ee6/ee2 +
        ee322) + 48 * ee62)) * ee6/ee2 + 16 * ee127 + 16 * ee101 +
        24 * ee205 + 8 * (ee4 * ee156 - ee276 * ee6/ee68) - 16 * ee102) *
        ee1 - y * (((((((2 * (ee257 - ee281) + 4 * ee257 + 8 *
        ee281)/ee30 + 6 * ee51) * ee26/ee70 - 2 * ee244) * ee1 +
        ee114) * ee26 + ee168 + ee294)/ee70 - ee263 * ee1/ee5) * ee1 +
        ee114)/ee30) * ee1/ee2);
      out(j, 5) += ((ee4 * (ee286 + ee114 + ee317 - (((80 - ee271)/
        ee2 - ee234) * ee1 + ee8 * (24 - ((ee84 + ee284 + ee295 + ee89) *
          ee4 + ee72 + ee218) * ee6/ee5) + ee285 + 2) * ee6/ee5) +
          (2 * ((2 * (2 * ee185 + 2 * ee66) + 8 * ee66) * ee6/ee2 +
          ee154) + ee296 + 24 * ee61 + 8 * ee215 - ((ee275 + (ee177 +
          ee79 - (((((2 - ee183) * ee26 + ee305)/ee4 + 2 * (ee78 + ee55/
            ee4)) * ee26 + 6 * (ee124/ee4))/ee68 + ee266/ee5) * ee6) *
              ee26 + ee166 + ee293)/ee111 + 2 * ((2 * (2 * ee187 + 2 *
              ee67) + 8 * ee67) * ee6/ee2 + ee155) + ee297 + 24 * ee62 + 8 *
              ee217)) * ee6/ee2 + 12 * ee203 + 6 * (ee4 * ee223 - (ee124 +
              (ee213 - 2 * (ee163/ee80)) * ee26) * ee6/ee68) + ee312 +
              ee318 + ee319 - ee320) * ee1 + y * (((((ee189 + ((2 * (ee207 -
              ee201) - ee302) * ee26/ee30 - 6 * ee124)/ee70) * ee1 - ee32) *
              ee26 - (ee275 + ee166 + ee293))/ee70 + (((72 - ee271)/
                ee2 - ee234) * ee1 + ee8 * (12 - ((ee284 + ee295) * ee4 + ee218) *
                  ee6/ee5) + 1 + ee90) * ee1/ee5) * ee1 - ee79)/ee30) *
                  ee1/ee2;
      out(j, 6) += (((((ee242 + (96 - (ee273 * ee45 + ee274 - ee258) *
        ee6/ee5)/ee2) * ee1 + ee8 * (16 - (ee4 * (2 * ee87 + ee25 +
        ee77) + ee218) * ee6/ee5) + 4) * ee6/ee5 - (ee114 + ee306 +
        ee317)) * ee4 + ((ee21 * (ee114 - ((2 + 2 * ee122 + 2 *
        ee38 + ee57)/ee5 + 2 * (ee124/ee80)) * ee6) + ee26 * (ee14 +
        ee305 - (((ee26 * (4 - ee183) + 8 * ee55) * ee21/ee4 - 2 *
        (ee200 * ee26))/ee68 + ee267) * ee6) - (ee256 + ee289))/ee111 +
        12 * ee62 + 2 * ((2 * (ee187 + ee67) + ee311) * ee6/ee2 +
        ee179) + 4 * ee217 + 4 * ee226 - (12 * ee61 + 2 * ((2 * (ee185 +
        ee66) + ee310) * ee6/ee2 + ee178) + 4 * ee215 + 4 *
        ee225)) * ee6/ee2 + 4 * (ee194 * ee4 + (ee21 * (ee32 - (2 *
        (ee43/ee5) + ee210) * ee6) - ee133) * ee6/ee68) + 4 * ee197 +
        ee309 - (ee307 + ee308 + 8 * ee203)) * ee1 - y * ((((ee242 +
        (52 - (ee274 - (ee45 * ee125 + ee258)) * ee6/ee5)/ee2) *
        ee1 + 1 + ee69) * ee1/ee5 + (((2 * (ee195 * ee26) - ee21 * ee302)/
          ee30 - 4 * ee138) * ee26 * ee1/ee70 + ee256 + ee289 -
            ((2 * (ee124/ee70) - 2 * ee189) * ee1 + ee79) * ee21)/ee70) *
            ee1 - ee32)/ee30) * ee1/ee2;
      out(j, 7) += ((((((((ee246 * ee6/ee5 + 12) * ee1/ee2 + ee285 +
        ee313) * ee4 + ee248 + (ee4 * (6 - 24 * ee29) + 6 * ee8) *
        ee8 - ee268) * ee6/ee5 - 92) * ee1/ee2 - 8) * ee6/ee5 + ee114 +
        6 * ee21) * ee4 + ((((((ee21 * (6 - ee183)/ee4 - ee291) *
        ee26 + (ee292 - 4 * ee198)/ee4)/ee68 + (3 + 3 * ee141 +
        6 * ee43 + ee224)/ee5) * ee6 - 24/ee4) * ee21 + ee48 * (3 *
        ee55 + ee314 - 3 * ee23) - ee229)/ee111 + 2 * ((3 * ee66 + ee185) *
        ee6/ee2 + ee61) + 6 * ee225 + 6 * ee61 - (2 * ((3 *
        ee67 + ee187) * ee6/ee2 + ee62) + 6 * ee226 + 6 * ee62)) * ee6/
          ee2 + 2 * (ee117 * ee4 + (ee253 + ee212 - ee12) * ee21 *
            ee6/ee68) + ee298 + ee299 - (ee300 + 6 * ee197)) * ee1 - y *
            ((((((ee246 * ee9/ee10 + ee313) * ee4 + ee248 - ee268) * ee6/
              ee5 - 30) * ee1/ee2 - 1) * ee1/ee5 + (ee229 - ((((((ee290 -
                ee315) * ee26 + 2 * (ee198 * ee30))/ee30 + 2 * ee198 - ee292)/
                  ee70 + 3 * ee251) * ee1 - ee314) * ee21 + 3 * (ee48 * ee55)))/
                    ee70) * ee1 + ee14)/ee30) * ee1/ee2;
      out(j, 8) += -(ee1 * (y * ((((((24 * ee4 + ee304 - ee282) *
        ee9/ee10 - 34) * ee4 - (14 * ee8 + ee76)) * ee9/ee10 + 15) *
        ee6/ee5 - ((ee255 - ee12) * ee48 + ee21 * (ee288 - ((ee290 +
        4 * ee160 - ee315)/ee30 + 4 * ee48) * ee21 * ee1/ee70) + ee287)/
          ee30) * ee1/ee2 - ee12)/ee30 - ((((((((ee21 * (ee183 -
            8) + ee316)/ee4 + ee291) * ee21/ee111 - (18 - ee82) * ee1/
              ee5) * ee6/ee2 + 10/ee4 + ee288 - ee316) * ee21 + (ee253 - ee12) *
                ee48 + ee287)/ee31 + (6 * ee67 + ee6 * (ee158 - ee157)/
                  ee2 - 6 * ee66) * ee6)/ee2 + 7 * ee62 - 7 * ee61) * ee6/ee2 +
                    (((((ee304 + 56 * ee4 - ee282) * ee9/ee10 - 86) * ee4 - (ee76 +
                    22 * ee8)) * ee9/ee10 + 65) * ee9/ee10 - 15/ee4) * ee4 +
                    ee102 + ee107 + pars2 - (ee101 + pars1)) * ee1)/ee2);
    }
    
  }
  
  return out;
  
}
