// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

double pgam(double y, double a, double b)
{
  return Rf_pgamma(y, a, b, 1, 0);
}

const double eps0 = 1e-8;
const double eps1 = 1e-4;

// double dpgama(double y, double a, double b, double eps)
// {
//   return 0.5 * (pgam(y, a + eps, b) - pgam(y, a - eps, b)) / eps;
// }
// 
// double dpgamb(double y, double a, double b) {
//   return y * R_pow(b * y, a - 1) * exp(-b*y) / gamma(a);
// }
//   
// double dpgama(double y, double a, double b) 
// {
//   return 0.5 * (pgam(y, a + eps0, b) - pgam(y, a - eps0, b)) / eps0;
// }
// 
// double dpgamaa(double y, double a, double b) 
// {
//   return (pgam(y, a + eps1, b) + pgam(y, a - eps1, b) - 2 * pgam(y, a, b)) / eps0;
// }
// 
// double dpgamab(double y, double a, double b) 
// {
//   return .25 * (pgam(y, a + eps1, b + eps1) + pgam(y, a - eps1, b - eps1) - 
//                 pgam(y, a + eps1, b - eps1) + pgam(y, a - eps1, b + eps1)) / eps0;
// }

// [[Rcpp::export]]
double ldenom(double left, double alpha, double beta) {
  double bl = beta * left;
  double out = 1.0 / tgamma(alpha + 1.0);
  for (int l=1; l < 13; l++) {
    out += R_pow(bl, l) / tgamma(alpha + l + 1.0);
  }
  out = out * R_pow(bl, alpha) * exp(-bl);
  return log(1 - out);
}

// double fj(double y, double p1, double p2, double left) {
//   double alpha = exp(p1);
//   double beta = exp(p2);
//   double ldenom = log(1 - pgam(left, alpha, beta));
//   double lnum = (alpha - 1) * log(y) - beta * y + alpha * p2 - lgamma(alpha);
//   return ldenom - lnum;
// }
//     
// arma::rowvec dfj(double y, double p1, double p2, double left) {
//   double f0 = fj(y, p1, p2, left);
//   arma::rowvec out(5);
//   double p1h = p1 + eps1;
//   double p2h = p2 + eps1;
//   double p1l = p1 - eps1;
//   double p2l = p2 - eps1;
//   double f1h = fj(y, p1h, p2, left);
//   double f1l = fj(y, p1l, p2, left);
//   double f2h = fj(y, p1, p2h, left);
//   double f2l = fj(y, p1, p2l, left);
//   out[0] = 0.5 * (f1h - f1l) / eps1;
//   out[1] = 0.5 * (f2h - f2l) / eps1;
//   out[2] = 0.25 * (f1h + f1l - 2 * f0) / eps0;
//   out[4] = 0.25 * (f2h + f2l - 2 * f0) / eps0;
//   double fhh = fj(y, p1h, p2h, left);
//   double fll = fj(y, p1l, p2l, left);
//   out[3] = 0.5 * (fhh - f1h - f2h + 2 * f0 - f1l - f2l + fll) / eps0;
//   return out;
// }
// 

// //' Left-truncated Gamma distribution negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each GEV parameter
// //' @param X1 a design matrix for the shape parameter
// //' @param X2 a design matrix for the rate parameter
// //' @param ymat a matrix
// //' @param leftmat a matrix
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return ltgammad0 a scalar, the negative log-likelihood
// //' @return ltgamma12 a matrix, first then second derivatives w.r.t. parameters
// //' @return ltgamma34 a matrix, third then fourth derivatives w.r.t. parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double ltgammad0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat leftmat)
{
    
  int nobs = nhere.size();
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, left, p1, p2, alpha, beta, lnum;

  double nllh = 0.0;

  for (int j=0; j < nobs; j++) {

    p1 = p1vec[j];
    p2 = p2vec[j];
    alpha = exp(p1);
    beta = exp(p2);

    for (int l=0; l < nhere[j]; l++) {

      y = ymat(j, l);
      left = leftmat(j, l);

      lnum = (alpha - 1) * log(y) - beta * y + alpha * p2 - lgamma(alpha);
      nllh += ldenom(left, alpha, beta) - lnum;

    }

  }

return(nllh);

}

// //' @rdname ltgammad0
// [[Rcpp::export]]
arma::mat ltgammad12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat leftmat)
{
    
int nobs = nhere.size();
arma::mat out = arma::mat(nobs, 5, arma::fill::zeros);
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);

  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, left, p1, p2;
  
  double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
  double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee39;
  double ee41, ee42, ee43, ee44, ee45, ee46, ee47, ee48, ee49;
  double ee50, ee51, ee52, ee53, ee56, ee57, ee58;
  double ee60, ee61, ee62, ee65, ee68;
  double ee70, ee72, ee74, ee75, ee76;
  
  for (int j=0; j < nobs; j++) {
    
    p1 = p1vec[j];
    p2 = p2vec[j];

    for (int l=0; l < nhere[j]; l++) {

      y = ymat(j, l);
      left = leftmat(j, l);
      
      ee1 = exp(p1);
      ee2 = exp(p2);
      ee3 = left * ee2;
      ee4 = R_pow(ee3, 2);
      ee5 = 10 + ee1;
      ee6 = 11 + ee1;
      ee7 = 12 + ee1;
      ee8 = 13 + ee1;
      ee9 = 2 + ee1;
      ee10 = 3 + ee1;
      ee11 = 4 + ee1;
      ee12 = 5 + ee1;
      ee13 = 6 + ee1;
      ee14 = 7 + ee1;
      ee15 = 8 + ee1;
      ee16 = 9 + ee1;
      ee17 = 1 + ee1;
      ee18 = tgamma(ee5);
      ee19 = tgamma(ee6);
      ee20 = tgamma(ee7);
      ee21 = tgamma(ee8);
      ee22 = tgamma(ee9);
      ee23 = tgamma(ee10);
      ee24 = tgamma(ee11);
      ee25 = tgamma(ee12);
      ee26 = tgamma(ee13);
      ee27 = tgamma(ee14);
      ee28 = tgamma(ee15);
      ee29 = tgamma(ee16);
      ee30 = tgamma(ee17);
      ee31 = R_pow(ee3, ee1);
      ee32 = ee4/ee21;
      ee33 = ee3/ee18;
      ee34 = ee3/ee20;
      ee35 = ee3/ee26;
      ee36 = ee3/ee28;
      ee39 = (((((ee32 + 1/ee19 + ee34) * ee4 + 1/ee29 + ee33) *  ee4 + 1/ee27 + ee36) * ee4 + 1/ee25 + ee35) * ee4 + 1/ee23 +  ee3/ee24) * ee4 + 1/ee30 + ee3/ee22;
      ee41 = exp(-ee3);
      ee42 = Rf_digamma(ee5);
      ee43 = Rf_digamma(ee6);
      ee44 = Rf_digamma(ee7);
      ee45 = Rf_digamma(ee8);
      ee46 = Rf_digamma(ee9);
      ee47 = Rf_digamma(ee10);
      ee48 = Rf_digamma(ee11);
      ee49 = Rf_digamma(ee12);
      ee50 = Rf_digamma(ee13);
      ee51 = Rf_digamma(ee14);
      ee52 = Rf_digamma(ee15);
      ee53 = Rf_digamma(ee16);
      ee56 = 1 - ee39 * ee41 * ee31;
      ee57 = Rf_digamma(ee17);
      ee58 = ee1 - 1;
      ee60 = log(left) + p2;
      ee61 = R_pow(ee3, ee58);
      ee62 = 1/ee22;
      ee65 = (((((11/ee20 + 12 * (ee3/ee21)) * ee4 + 10 * (ee3/ee19) +  9/ee18) * ee4 + 7/ee28 + 8 * (ee3/ee29)) * ee4 + 5/ee26 +  6 * (ee3/ee27)) * ee4 + 3/ee24 + 4 * (ee3/ee25)) * ee4 +  ee62 + 2 * (ee3/ee23);
      ee68 = (((((ee43/ee19 + ee45 * ee4/ee21 + left * ee44 *  ee2/ee20) * ee4 + ee53/ee29 + left * ee42 * ee2/ee18) *  ee4 + ee51/ee27 + left * ee52 * ee2/ee28) * ee4 + ee49/ee25 +  left * ee50 * ee2/ee26) * ee4 + ee47/ee23 + left * ee48 *  ee2/ee24) * ee4 + ee57/ee30 + left * ee46 * ee2/ee22;
      ee70 = ee1 * ee61 - ee31;
      ee72 = ee39 * ee70 + ee65 * ee31;
      ee74 = ee39 * ee60 - ee68;
      ee75 = Rf_digamma(ee1);
      ee76 = log(y);
      
      out(j, 0) += -((ee74 * ee41 * ee31/ee56 + ee76 + p2 - ee75) *
        ee1);
      out(j, 1) += -(ee1 + ee2 * (left * ee72 * ee41/ee56 - y));
      out(j, 2) += -(((((ee39 * (1 + ee1 * ee60) - 2 * (ee68 * ee1)) *
        ee60 - ((((((((1 - 2 * (ee43 * ee1)) * ee43 + (R_pow(ee43, 2) +
        Rf_trigamma(ee6)) * ee1)/ee19 + ((1 - 2 * (ee45 * ee1)) *
        ee45 + (R_pow(ee45, 2) + Rf_trigamma(ee8)) * ee1) * ee4/ee21 +
        left * ((1 - 2 * (ee44 * ee1)) * ee44 + (R_pow(ee44, 2) +
        Rf_trigamma(ee7)) * ee1) * ee2/ee20) * ee4 + ((1 - 2 * (ee53 *
        ee1)) * ee53 + (R_pow(ee53, 2) + Rf_trigamma(ee16)) * ee1)/ee29 +
        left * ((1 - 2 * (ee42 * ee1)) * ee42 + (R_pow(ee42, 2) +
        Rf_trigamma(ee5)) * ee1) * ee2/ee18) * 
        ee4 + ((1 - 2 *
        (ee51 * ee1)) * ee51 + (R_pow(ee51, 2) + Rf_trigamma(ee14)) * ee1)/
          ee27 + left * ((1 - 2 * (ee52 * ee1)) * ee52 + (R_pow(ee52, 2) +
            Rf_trigamma(ee15)) * ee1) * ee2/ee28) * ee4 + ((1 - 2 *
            (ee49 * ee1)) * ee49 + (R_pow(ee49, 2) + Rf_trigamma(ee12)) *
            ee1)/ee25 + left * ((1 - 2 * (ee50 * ee1)) * ee50 + (R_pow(ee50, 2) +
            Rf_trigamma(ee13)) * ee1) * ee2/ee26) * ee4 + ((1 -
            2 * (ee47 * ee1)) * ee47 + (R_pow(ee47, 2) + Rf_trigamma(ee10)) *
            ee1)/ee23 + left * ((1 - 2 * (ee48 * ee1)) * ee48 + (R_pow(ee48, 
                                 2) +
                                   Rf_trigamma(ee11)) * ee1) * ee2/ee24) * ee4 +
                                   ((1 - 2 * (ee57 * ee1)) * ee57 + (R_pow(ee57, 2) + Rf_trigamma(ee17)) *
                                   ee1)/ee30 + left * ((1 - 2 * (ee46 * ee1)) * ee46 +
                                   (R_pow(ee46, 2) + Rf_trigamma(ee9)) * ee1) * ee2/ee22)) * ee31 +
                                   R_pow(ee74, 2) * ee41 * ee1 * R_pow(ee3, (2 * ee1))/ee56) *
                                   ee41/ee56 + ee76 + p2 - (ee75 + ee1 * Rf_trigamma(ee1))) *
                                   ee1);
      out(j, 3) += -((((ee39 * ((ee1 - ee3) * ee60 + 1) + left * ee65 *
        ee2 * ee60) * ee31 + left * ((ee72 * ee74 * ee41/ee56 -
        ((((((11 * (ee44/ee20) + 12 * (left * ee45 * ee2/ee21)) *
        ee4 + 10 * (left * ee43 * ee2/ee19) + 9 * (ee42/ee18)) * ee4 +
        7 * (ee52/ee28) + 8 * (left * ee53 * ee2/ee29)) * ee4 + 5 *
        (ee50/ee26) + 6 * (left * ee51 * ee2/ee27)) * ee4 + 3 * (ee48/
          ee24) + 4 * (left * ee49 * ee2/ee25)) * ee4 + 2 * (left *
            ee47 * ee2/ee23) + ee46/ee22)) * ee31 - ee68 * ee70) * ee2) *
            ee41/ee56 + 1) * ee1);
      out(j, 4) += -(ee2 * (left * (ee39 * ((1 + left * (ee58/ee3 -
        1) * ee2) * ee1 - left * (ee17 - ee3) * ee2) * ee61 + (ee62 +
        3 * ((ee4 + 2 * (R_pow(left, 2) * R_pow(ee2, 2)))/ee24) +
        left * (((((100/ee19 + 121 * ee34 + 144 * ee32) * ee4 + 64/
          ee29 + 81 * ee33) * ee4 + 36/ee27 + 49 * ee36) * ee4 + 16/
            ee25 + 25 * ee35) * ee4 + 4/ee23) * ee2) * ee31 + left * (R_pow(ee72, 2) *
              ee41/ee56 + 2 * (ee65 * ee70)) * ee2) * ee41/
                ee56 - y));
      
}

}

return out;

}

// //' @rdname ltgammad0
// [[Rcpp::export]]
arma::mat ltgammad34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat leftmat)
{
  
  int nobs = nhere.size();
  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, left, p1, p2;
  
  double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
  double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee37, ee38, ee39;
  double ee40, ee41, ee42, ee43, ee44, ee45, ee46, ee47, ee48, ee49;
  double ee50, ee52, ee55, ee56, ee57, ee59;
  double ee60, ee61, ee62, ee63, ee64, ee65, ee66, ee67, ee68, ee69;
  double ee70, ee71, ee72, ee73, ee74, ee75, ee76, ee77, ee78, ee79;
  double ee80, ee81, ee82, ee83, ee84, ee85, ee86, ee87, ee88, ee89;
  double ee91, ee92, ee93, ee94, ee97;
  double ee100, ee103, ee106, ee107;
  double ee110, ee112, ee113, ee116, ee117, ee118, ee119;
  double ee120, ee121, ee122, ee123, ee124, ee125, ee126, ee127, ee128, ee129;
  double ee130, ee131, ee132, ee135, ee136, ee137, ee138, ee139;
  double ee140, ee141, ee142, ee143, ee144, ee145, ee146, ee147, ee148, ee149;
  double ee150, ee151, ee152, ee154, ee156, ee158, ee159;
  double ee160, ee162, ee163, ee165;
  double ee170, ee172, ee173, ee174, ee175, ee176, ee177, ee178;
  double ee184;
  double ee191, ee192, ee193, ee194, ee197;
  double ee200, ee203, ee206, ee209;
  double ee212, ee215, ee218;
  double ee221, ee224, ee225, ee227;
  double ee230, ee231, ee233, ee235, ee236, ee237, ee238, ee239;
  double ee240, ee241, ee242, ee243, ee244, ee245, ee246, ee247;
  double ee250, ee251, ee252, ee253, ee254, ee255, ee256, ee257, ee258, ee259;
  double ee260, ee261, ee262, ee263, ee264, ee265, ee266, ee267, ee268, ee269;
  double ee270, ee271, ee272, ee273, ee274, ee275, ee276, ee277, ee278, ee279;
  double ee280, ee281, ee282, ee283, ee284, ee285, ee286, ee287, ee288, ee289;
  double ee290, ee291, ee292, ee293, ee294, ee295, ee296, ee297, ee298, ee299;
  double ee300, ee301, ee302, ee303, ee304, ee305, ee306, ee307, ee308, ee309;
  double ee310, ee311, ee312, ee313, ee314, ee315, ee316, ee317, ee318, ee319;
  double ee320, ee321, ee322, ee323, ee324, ee325, ee326, ee327, ee328, ee329;
  double ee330, ee331, ee332, ee333, ee334, ee336, ee338;
  double ee340, ee342, ee344, ee346, ee348;
  double ee350, ee352, ee354, ee356, ee358, ee359;
  double ee360, ee361, ee362, ee363, ee364, ee365, ee366, ee367, ee369;
  double ee370, ee372, ee373, ee375, ee376, ee377, ee378, ee379;
  double ee380, ee381, ee382, ee383, ee384, ee385, ee386, ee387, ee388, ee389;
  double ee390, ee393, ee394, ee395, ee396, ee397, ee399;
  double ee400, ee401, ee402, ee404, ee405, ee407, ee408, ee409;
  double ee410, ee411, ee412, ee413, ee414, ee415, ee416, ee417, ee418, ee419;
  double ee420, ee421, ee424;
  double ee432, ee437;
  double ee441, ee444, ee446, ee449;
  double ee453, ee456, ee458, ee459;
  double ee461, ee462, ee465, ee466, ee467, ee468;
  double ee470, ee474, ee475, ee477;
  double ee487;
  double ee491, ee492, ee493, ee496, ee498, ee499;
  double ee500, ee502, ee504, ee505, ee506, ee507, ee508, ee509;
  double ee510, ee511, ee512, ee513, ee514, ee515, ee516, ee517, ee518, ee519;
  double ee520, ee521, ee522, ee523, ee524, ee525, ee529;
  double ee531, ee532, ee533, ee534, ee535, ee536, ee538, ee539;
  double ee540, ee541, ee543, ee544, ee545, ee546, ee547, ee548;
  double ee550, ee551, ee552, ee553;
  
  for (int j=0; j < nobs; j++) {
    
    p1 = p1vec[j];
    p2 = p2vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      left = leftmat(j, l);
      
      ee1 = exp(p1);
      ee2 = exp(p2);
      ee3 = left * ee2;
      ee4 = R_pow(ee3, 2);
      ee5 = 10 + ee1;
      ee6 = 11 + ee1;
      ee7 = 12 + ee1;
      ee8 = 13 + ee1;
      ee9 = 3 + ee1;
      ee10 = 4 + ee1;
      ee11 = 5 + ee1;
      ee12 = 6 + ee1;
      ee13 = 7 + ee1;
      ee14 = 8 + ee1;
      ee15 = 9 + ee1;
      ee16 = 2 + ee1;
      ee17 = 1 + ee1;
      ee18 = tgamma(ee5);
      ee19 = tgamma(ee6);
      ee20 = tgamma(ee7);
      ee21 = tgamma(ee8);
      ee22 = tgamma(ee16);
      ee23 = tgamma(ee9);
      ee24 = tgamma(ee10);
      ee25 = tgamma(ee11);
      ee26 = tgamma(ee12);
      ee27 = tgamma(ee13);
      ee28 = tgamma(ee14);
      ee29 = tgamma(ee15);
      ee30 = Rf_digamma(ee5);
      ee31 = Rf_digamma(ee6);
      ee32 = Rf_digamma(ee7);
      ee33 = Rf_digamma(ee8);
      ee34 = Rf_digamma(ee9);
      ee35 = Rf_digamma(ee10);
      ee36 = Rf_digamma(ee11);
      ee37 = Rf_digamma(ee12);
      ee38 = Rf_digamma(ee13);
      ee39 = Rf_digamma(ee14);
      ee40 = Rf_digamma(ee15);
      ee41 = Rf_digamma(ee16);
      ee42 = tgamma(ee17);
      ee43 = Rf_digamma(ee17);
      ee44 = R_pow(ee3, ee1);
      ee45 = ee4/ee21;
      ee46 = ee3/ee18;
      ee47 = ee3/ee20;
      ee48 = ee3/ee26;
      ee49 = ee3/ee28;
      ee50 = ee1 - 1;
      ee52 = log(left) + p2;
      ee55 = (((((ee45 + 1/ee19 + ee47) * ee4 + 1/ee29 + ee46) *  ee4 + 1/ee27 + ee49) * ee4 + 1/ee25 + ee48) * ee4 + 1/ee23 +  ee3/ee24) * ee4 + 1/ee42 + ee3/ee22;
      ee56 = R_pow(ee3, ee50);
      ee57 = 1/ee22;
      ee59 = exp(-ee3);
      ee60 = Rf_trigamma(ee5);
      ee61 = Rf_trigamma(ee6);
      ee62 = Rf_trigamma(ee7);
      ee63 = Rf_trigamma(ee8);
      ee64 = Rf_trigamma(ee16);
      ee65 = Rf_trigamma(ee9);
      ee66 = Rf_trigamma(ee10);
      ee67 = Rf_trigamma(ee11);
      ee68 = Rf_trigamma(ee12);
      ee69 = Rf_trigamma(ee13);
      ee70 = Rf_trigamma(ee14);
      ee71 = Rf_trigamma(ee15);
      ee72 = R_pow(ee30, 2);
      ee73 = R_pow(ee31, 2);
      ee74 = R_pow(ee32, 2);
      ee75 = R_pow(ee33, 2);
      ee76 = R_pow(ee34, 2);
      ee77 = R_pow(ee35, 2);
      ee78 = R_pow(ee36, 2);
      ee79 = R_pow(ee37, 2);
      ee80 = R_pow(ee38, 2);
      ee81 = R_pow(ee39, 2);
      ee82 = R_pow(ee40, 2);
      ee83 = ee3/ee19;
      ee84 = ee3/ee21;
      ee85 = ee3/ee27;
      ee86 = ee3/ee29;
      ee87 = Rf_trigamma(ee17);
      ee88 = ee34/ee23;
      ee89 = ee31/ee19;
      ee91 = ee33 * ee4/ee21;
      ee92 = ee36/ee25;
      ee93 = ee38/ee27;
      ee94 = ee40/ee29;
      ee97 = left * ee30 * ee2/ee18;
      ee100 = left * ee32 * ee2/ee20;
      ee103 = left * ee37 * ee2/ee26;
      ee106 = left * ee39 * ee2/ee28;
      ee107 = ee3/ee23;
      ee110 = (((((11/ee20 + 12 * ee84) * ee4 + 10 * ee83 + 9/ee18) *  ee4 + 7/ee28 + 8 * ee86) * ee4 + 5/ee26 + 6 * ee85) *  ee4 + 3/ee24 + 4 * (ee3/ee25)) * ee4 + ee57 + 2 * ee107;
      ee112 = ee1 * ee56 - ee44;
      ee113 = R_pow(ee41, 2);
      ee116 = (((((ee89 + ee91 + ee100) * ee4 + ee94 + ee97) *  ee4 + ee93 + ee106) * ee4 + ee92 + ee103) * ee4 + ee88 +  left * ee35 * ee2/ee24) * ee4 + ee43/ee42 + left * ee41 *  ee2/ee22;
      ee117 = R_pow(ee2, 2);
      ee118 = R_pow(left, 2);
      ee119 = ee118 * ee117;
      ee120 = R_pow(ee43, 2);
      ee121 = (ee72 + ee60) * ee1;
      ee122 = (ee73 + ee61) * ee1;
      ee123 = (ee74 + ee62) * ee1;
      ee124 = (ee75 + ee63) * ee1;
      ee125 = (ee113 + ee64) * ee1;
      ee126 = (ee76 + ee65) * ee1;
      ee127 = (ee77 + ee66) * ee1;
      ee128 = (ee78 + ee67) * ee1;
      ee129 = (ee79 + ee68) * ee1;
      ee130 = (ee80 + ee69) * ee1;
      ee131 = (ee81 + ee70) * ee1;
      ee132 = (ee82 + ee71) * ee1;
      ee135 = 1 - ee55 * ee59 * ee44;
      ee136 = ee1 * ee52;
      ee137 = (ee120 + ee87) * ee1;
      ee138 = ee4 + 2 * ee119;
      ee139 = ee41 * ee1;
      ee140 = ee50/ee3;
      ee141 = ee30 * ee1;
      ee142 = ee31 * ee1;
      ee143 = ee32 * ee1;
      ee144 = ee33 * ee1;
      ee145 = ee34 * ee1;
      ee146 = ee35 * ee1;
      ee147 = ee36 * ee1;
      ee148 = ee37 * ee1;
      ee149 = ee38 * ee1;
      ee150 = ee39 * ee1;
      ee151 = ee40 * ee1;
      ee152 = ee43 * ee1;
      ee154 = ee17 - ee3;
      ee156 = left * (ee140 - 1) * ee2;
      ee158 = (1 - 2 * ee139) * ee41 + ee125;
      ee159 = 1 + ee136;
      ee160 = ee41/ee22;
      ee162 = ee55 * ee112 + ee110 * ee44;
      ee163 = left * ee110;
      ee165 = ee55 * ee52 - ee116;
      ee170 = (1 + ee156) * ee1 - left * ee154 * ee2;
      ee172 = ee57 + 3 * (ee138/ee24) + left * (((((100/ee19 +  121 * ee47 + 144 * ee45) * ee4 + 64/ee29 + 81 * ee46) *  ee4 + 36/ee27 + 49 * ee49) * ee4 + 16/ee25 + 25 * ee48) *  ee4 + 4/ee23) * ee2;
      ee173 = ((1 - 2 * ee152) * ee43 + ee137)/ee42;
      ee174 = ((1 - 2 * ee142) * ee31 + ee122)/ee19;
      ee175 = ((1 - 2 * ee145) * ee34 + ee126)/ee23;
      ee176 = ((1 - 2 * ee147) * ee36 + ee128)/ee25;
      ee177 = ((1 - 2 * ee149) * ee38 + ee130)/ee27;
      ee178 = ((1 - 2 * ee151) * ee40 + ee132)/ee29;
      ee184 = (1 - 2 * ee144) * ee33 + ee124;
      ee191 = ee30/ee18;
      ee192 = ee32/ee20;
      ee193 = ee37/ee26;
      ee194 = ee39/ee28;
      ee197 = left * ((1 - 2 * ee141) * ee30 + ee121) * ee2/ee18;
      ee200 = left * ((1 - 2 * ee143) * ee32 + ee123) * ee2/ee20;
      ee203 = left * ee158 * ee2/ee22;
      ee206 = left * ((1 - 2 * ee146) * ee35 + ee127) * ee2/ee24;
      ee209 = left * ((1 - 2 * ee148) * ee37 + ee129) * ee2/ee26;
      ee212 = left * ((1 - 2 * ee150) * ee39 + ee131) * ee2/ee28;
      ee215 = left * ee31 * ee2/ee19;
      ee218 = left * ee33 * ee2/ee21;
      ee221 = left * ee38 * ee2/ee27;
      ee224 = left * ee40 * ee2/ee29;
      ee225 = (((((11 * ee192 + 12 * ee218) * ee4 + 10 * ee215 +  9 * ee191) * ee4 + 7 * ee194 + 8 * ee224) * ee4 + 5 *  ee193 + 6 * ee221) * ee4 + 3 * (ee35/ee24) + 4 * (left *  ee36 * ee2/ee25)) * ee4;
      ee227 = (ee1 - ee3) * ee52 + 1;
      ee230 = left * ee34 * ee2/ee23;
      ee231 = (((((ee174 + ee184 * ee4/ee21 + ee200) * ee4 + ee178 +  ee197) * ee4 + ee177 + ee212) * ee4 + ee176 + ee209) *  ee4 + ee175 + ee206) * ee4;
      ee233 = ee225 + 2 * ee230 + ee160;
      ee235 = ee231 + ee173 + ee203;
      ee236 = ee121 + ee30;
      ee237 = ee122 + ee31;
      ee238 = ee123 + ee32;
      ee239 = ee124 + ee33;
      ee240 = ee126 + ee34;
      ee241 = ee127 + ee35;
      ee242 = ee128 + ee36;
      ee243 = ee129 + ee37;
      ee244 = ee130 + ee38;
      ee245 = ee131 + ee39;
      ee246 = ee132 + ee40;
      ee247 = R_pow(ee3, (ee1 - 2));
      ee250 = ee55 * ee170 * ee56 + ee172 * ee44;
      ee251 = ee72 * ee1;
      ee252 = ee73 * ee1;
      ee253 = ee74 * ee1;
      ee254 = ee75 * ee1;
      ee255 = ee76 * ee1;
      ee256 = ee77 * ee1;
      ee257 = ee78 * ee1;
      ee258 = ee79 * ee1;
      ee259 = ee80 * ee1;
      ee260 = ee81 * ee1;
      ee261 = ee82 * ee1;
      ee262 = (3 + ee136) * ee1;
      ee263 = (ee55 * ee159 - 2 * (ee116 * ee1)) * ee52;
      ee264 = 8 * ee257;
      ee265 = Rf_psigamma(ee5, 2);
      ee266 = Rf_psigamma(ee6, 2);
      ee267 = Rf_psigamma(ee7, 2);
      ee268 = Rf_psigamma(ee8, 2);
      ee269 = Rf_psigamma(ee16, 2);
      ee270 = Rf_psigamma(ee9, 2);
      ee271 = Rf_psigamma(ee10, 2);
      ee272 = Rf_psigamma(ee11, 2);
      ee273 = Rf_psigamma(ee12, 2);
      ee274 = Rf_psigamma(ee13, 2);
      ee275 = Rf_psigamma(ee14, 2);
      ee276 = Rf_psigamma(ee15, 2);
      ee277 = (ee55 * ee227 + ee163 * ee2 * ee52) * ee44;
      ee278 = ee250 + 2 * (ee163 * ee112 * ee2);
      ee279 = ee116 * ee112;
      ee280 = 4 * ee242;
      ee281 = Rf_psigamma(ee17, 2);
      ee282 = ee263 - ee235;
      ee283 = ee17 * ee50;
      ee284 = (2 * ee60 + ee72) * ee1;
      ee285 = (2 * ee61 + ee73) * ee1;
      ee286 = (2 * ee62 + ee74) * ee1;
      ee287 = (2 * ee63 + ee75) * ee1;
      ee288 = (2 * ee64 + ee113) * ee1;
      ee289 = (2 * ee65 + ee76) * ee1;
      ee290 = (2 * ee66 + ee77) * ee1;
      ee291 = (2 * ee67 + ee78) * ee1;
      ee292 = (2 * ee68 + ee79) * ee1;
      ee293 = (2 * ee69 + ee80) * ee1;
      ee294 = (2 * ee70 + ee81) * ee1;
      ee295 = (2 * ee71 + ee82) * ee1;
      ee296 = (ee30 * ee60 + ee265) * ee1;
      ee297 = (ee31 * ee61 + ee266) * ee1;
      ee298 = (ee32 * ee62 + ee267) * ee1;
      ee299 = (ee33 * ee63 + ee268) * ee1;
      ee300 = (ee41 * ee64 + ee269) * ee1;
      ee301 = ee125 + ee41;
      ee302 = (ee34 * ee65 + ee270) * ee1;
      ee303 = (ee35 * ee66 + ee271) * ee1;
      ee304 = (ee36 * ee67 + ee272) * ee1;
      ee305 = (ee37 * ee68 + ee273) * ee1;
      ee306 = (ee38 * ee69 + ee274) * ee1;
      ee307 = (ee39 * ee70 + ee275) * ee1;
      ee308 = (ee40 * ee71 + ee276) * ee1;
      ee309 = 3 * ee60;
      ee310 = 3 * ee61;
      ee311 = 3 * ee62;
      ee312 = 3 * ee63;
      ee313 = 3 * ee64;
      ee314 = 3 * ee65;
      ee315 = 3 * ee66;
      ee316 = 3 * ee67;
      ee317 = 3 * ee68;
      ee318 = 3 * ee69;
      ee319 = 3 * ee70;
      ee320 = 3 * ee71;
      ee321 = ee119/ee4;
      ee322 = ee277 - left * (ee233 * ee44 + ee279) * ee2;
      ee323 = ((ee284 + 3 * ee30) * ee30 + ee296 + ee309) * ee1;
      ee324 = ((ee285 + 3 * ee31) * ee31 + ee297 + ee310) * ee1;
      ee325 = ((ee286 + 3 * ee32) * ee32 + ee298 + ee311) * ee1;
      ee326 = ((ee287 + 3 * ee33) * ee33 + ee299 + ee312) * ee1;
      ee327 = ((ee288 + 3 * ee41) * ee41 + ee300 + ee313) * ee1;
      ee328 = ((ee289 + 3 * ee34) * ee34 + ee302 + ee314) * ee1;
      ee329 = ((ee290 + 3 * ee35) * ee35 + ee303 + ee315) * ee1;
      ee330 = ((ee291 + 3 * ee36) * ee36 + ee304 + ee316) * ee1;
      ee331 = ((ee292 + 3 * ee37) * ee37 + ee305 + ee317) * ee1;
      ee332 = ((ee293 + 3 * ee38) * ee38 + ee306 + ee318) * ee1;
      ee333 = ((ee294 + 3 * ee39) * ee39 + ee307 + ee319) * ee1;
      ee334 = ((ee295 + 3 * ee40) * ee40 + ee308 + ee320) * ee1;
      ee336 = (1 + ee141) * ee30 + ee121;
      ee338 = (1 + ee142) * ee31 + ee122;
      ee340 = (1 + ee143) * ee32 + ee123;
      ee342 = (1 + ee144) * ee33 + ee124;
      ee344 = (1 + ee139) * ee41 + ee125;
      ee346 = (1 + ee145) * ee34 + ee126;
      ee348 = (1 + ee146) * ee35 + ee127;
      ee350 = (1 + ee147) * ee36 + ee128;
      ee352 = (1 + ee148) * ee37 + ee129;
      ee354 = (1 + ee149) * ee38 + ee130;
      ee356 = (1 + ee150) * ee39 + ee131;
      ee358 = (1 + ee151) * ee40 + ee132;
      ee359 = ee283 + 1;
      ee360 = (2 * ee87 + ee120) * ee1;
      ee361 = (ee43 * ee87 + ee281) * ee1;
      ee362 = ee137 + ee43;
      ee363 = ee56 + left * ee50 * ee2 * ee247;
      ee364 = 3 * ee87;
      ee365 = R_pow(ee162, 2);
      ee366 = ((ee360 + 3 * ee43) * ee43 + ee361 + ee364) * ee1;
      ee367 = ee363 * ee52;
      ee369 = (1 + ee152) * ee43 + ee137;
      ee370 = 2 * ee344;
      ee372 = 8 * ee4;
      ee373 = ee113 * ee1;
      ee375 = left * (3 - ee3) * ee2;
      ee376 = R_pow(ee3, (2 * ee1));
      ee377 = ee4 + 8 * ee119;
      ee378 = 1 - ee321;
      ee379 = 2 * ee369;
      ee380 = 2 * ee336;
      ee381 = 2 * ee338;
      ee382 = 2 * ee340;
      ee383 = 2 * ee342;
      ee384 = 2 * ee346;
      ee385 = 2 * ee348;
      ee386 = 2 * ee350;
      ee387 = 2 * ee352;
      ee388 = 2 * ee354;
      ee389 = 2 * ee356;
      ee390 = 2 * ee358;
      ee393 = 8 * ee373;
      ee394 = 8 * ee255;
      ee395 = ee372 + left * (3 * ee138 + 3 * ee4) * ee2;
      ee396 = ee120 * ee1;
      ee397 = left * ee233;
      ee399 = left * ee1 * ee2;
      ee400 = ee327 + (1 - (ee370 + 4 * ee301 - ee393) * ee1) *  ee41;
      ee401 = ee158/ee22;
      ee402 = ee378 * ee44;
      ee404 = ee262 * ee52 + 1;
      ee405 = 2 * ee240;
      ee407 = 3 * ee241 - 6 * ee256;
      ee408 = 4 * ee240;
      ee409 = 4 * ee255;
      ee410 = 8 * ee396;
      ee411 = 8 * ee251;
      ee412 = 8 * ee252;
      ee413 = 8 * ee253;
      ee414 = 8 * ee254;
      ee415 = 8 * ee256;
      ee416 = 8 * ee258;
      ee417 = 8 * ee259;
      ee418 = 8 * ee260;
      ee419 = 8 * ee261;
      ee420 = ((((((ee324 + (1 - (ee381 + 4 * ee237 - ee412) *  ee1) * ee31)/ee19 + (ee326 + (1 - (ee383 + 4 * ee239 -  ee414) * ee1) * ee33) * ee4/ee21 + left * (ee325 + (1 -  (ee382 + 4 * ee238 - ee413) * ee1) * ee32) * ee2/ee20) *  ee4 + (ee334 + (1 - (ee390 + 4 * ee246 - ee419) * ee1) *  ee40)/ee29 + left * (ee323 + (1 - (ee380 + 4 * ee236 -  ee411) * ee1) * ee30) * ee2/ee18) * ee4 + (ee332 + (1 -  (ee388 + 4 * ee244 - ee417) * ee1) * ee38)/ee27 + left *  (ee333 + (1 - (ee389 + 4 * ee245 - ee418) * ee1) * ee39) *  ee2/ee28) * ee4 + (ee330 + (1 - (ee386 + ee280 - ee264) *  ee1) * ee36)/ee25 + left * (ee331 + (1 - (ee387 + 4 *  ee243 - ee416) * ee1) * ee37) * ee2/ee26) * ee4 + (ee328 +  (1 - (ee384 + ee408 - ee394) * ee1) * ee34)/ee23 + left *  (ee329 + (1 - (ee385 + 4 * ee241 - ee415) * ee1) * ee35) *  ee2/ee24) * ee4;
      ee421 = ((((((11 * ee238 - 22 * ee253)/ee20 + left * (12 *  ee239 - 24 * ee254) * ee2/ee21) * ee4 + (9 * ee236 -  18 * ee251)/ee18 + left * (10 * ee237 - 20 * ee252) *  ee2/ee19) * ee4 + (7 * ee245 - 14 * ee260)/ee28 + left *  (8 * ee246 - 16 * ee261) * ee2/ee29) * ee4 + (5 * ee243 -  10 * ee258)/ee26 + left * (6 * ee244 - 12 * ee259) *  ee2/ee27) * ee4 + ee407/ee24 + left * (ee280 - ee264) *  ee2/ee25) * ee4;
      ee424 = (ee366 + (1 - (ee379 + 4 * ee362 - ee410) * ee1) *  ee43)/ee42;
      ee432 = (ee262 - left * ee159 * ee2) * ee52 + 1;
      ee437 = ee402 + left * ((ee367 + (2 - ee3 * ee52) * ee56) *  ee1 - (ee154 * ee52 + 2) * ee44) * ee2;
      ee441 = ee57 + 3 * (ee377/ee24) + 4 * (ee395/ee25) + left *  (((((2662/ee20 + 3456 * ee84) * ee4 + 1458/ee18 + 2000 *  ee83) * ee4 + 1024 * ee86 + 686/ee28) * ee4 + 250/ee26 +  432 * ee85) * ee4 + 8/ee23) * ee2;
      ee444 = 3 * (ee138 * ee35/ee24);
      ee446 = left * (((((100 * ee89 + 121 * ee100 + 144 * ee91) *  ee4 + 64 * ee94 + 81 * ee97) * ee4 + 36 * ee93 + 49 *  ee106) * ee4 + 16 * ee92 + 25 * ee103) * ee4 + 4 * ee88) *  ee2;
      ee449 = left * ee400 * ee2/ee22;
      ee453 = left * ((ee359 - ee399) * ee247 - (2 + ee156) *  ee56) * ee1 * ee2 - 2 * (((3 + left * (ee140 - 2) * ee2) *  ee1 + 1 - ee375) * ee56);
      ee456 = left * (ee405 - ee409) * ee2/ee23;
      ee458 = ee420 + ee424 + ee449;
      ee459 = R_pow(ee165, 2);
      ee461 = ee421 + ee401 + ee456;
      ee462 = left * ee162;
      ee465 = left * ee365 * ee59 * ee2;
      ee466 = ee235 * ee112;
      ee467 = (ee55 * ee432 + ee163 * ee159 * ee2 * ee52 - 2 *  ((ee116 * ee227 + ee397 * ee2 * ee52) * ee1)) * ee44;
      ee468 = ee278 * ee135;
      ee470 = (ee55 * ee404 - (3 * ee231 + 3 * (ee116 * ee159) +  3 * ee173 + 3 * ee203) * ee1) * ee52 - ee458;
      ee474 = ee459 * ee59 * ee1;
      ee475 = ee55 * ee437;
      ee477 = ee55 * ee453 + ee441 * ee44;
      ee487 = (ee172 * ee52 + 2 * (ee110 * ee227) - (ee116 * ee170/ee3 +  ee444 + ee160 + ee446)) * ee44;
      ee491 = 2 * (ee397 * ee112 * ee2);
      ee492 = 3 * (ee110 * ee170 * ee56);
      ee493 = 3 * (ee172 * ee112);
      ee496 = ee462 * ee165 * ee59 * ee2;
      ee498 = left * (2 * (ee365 * ee59/ee135) + 2 * (ee110 *  ee112)) * ee2;
      ee499 = ee322 * ee162;
      ee500 = ee322 * ee135;
      ee502 = ee467 - left * (ee466 + ee461 * ee44) * ee2;
      ee504 = ee278 * ee165 * ee44;
      ee505 = ee468 - ee465;
      ee506 = ee477 + left * (ee492 + ee493) * ee2;
      ee507 = ee323 + ee30;
      ee508 = ee324 + ee31;
      ee509 = ee325 + ee32;
      ee510 = ee326 + ee33;
      ee511 = ee328 + ee34;
      ee512 = ee329 + ee35;
      ee513 = ee330 + ee36;
      ee514 = ee331 + ee37;
      ee515 = ee332 + ee38;
      ee516 = ee333 + ee39;
      ee517 = ee334 + ee40;
      ee518 = R_pow(ee3, (3 * ee1));
      ee519 = 1 - ee3;
      ee520 = ee3/ee4;
      ee521 = ee502 * ee165;
      ee522 = ee470 * ee165;
      ee523 = ee499 + ee504;
      ee524 = ee282 * ee135;
      ee525 = ee263 + ((((((2 * (ee474 * R_pow(ee3, (ee1 - 12))/ee135) - ee184/ee21) * ee4 - (ee174 + ee200)) * ee4 - (ee178 + ee197)) * ee4 - (ee177 + ee212)) * ee4 - (ee176 + ee209)) * ee4 - (ee175 + ee206)) * ee4;
      ee529 = ee475 + left * (ee487 - ee491) * ee2;
      ee531 = ee250 + 2 * ee278 + ee498;
      ee532 = ee250 + ee498;
      ee533 = ee366 + ee43;
      ee534 = ee327 + ee41;
      ee535 = ee173 + ee203;
      ee536 = ((6 + ee136) * ee1 * ee52 + 7) * ee1;
      ee538 = 1 - ee375;
      ee539 = 16 * ee242;
      ee540 = 2 * ee505;
      ee541 = 2 * (ee162 * ee165 * ee59/ee135);
      ee543 = ee444 + ee160 + ee446;
      ee544 = 32 * ee257;
      ee545 = 4 * ee513;
      ee546 = 6 * ee348;
      ee547 = 8 * ee465;
      ee548 = Rf_digamma(ee1);
      ee550 = left * (ee367 + 2 * ee56) * ee2;
      ee551 = log(y);
      ee552 = Rf_psigamma(ee1, 2);
      ee553 = Rf_trigamma(ee1);
      
      out(j, 0) += -(((ee470 * ee44 + (ee525 + 2 * ee282 - ee535) *
        ee165 * ee59 * ee1 * ee376/ee135) * ee59/ee135 + ee551 + p2 -
        ((3 * ee553 + ee1 * ee552) * ee1 + ee548)) * ee1);
      out(j, 1) += -(((ee467 + ee165 * (2 * (ee322 * ee44) + 2 * (ee496 *
        ee376/ee135)) * ee59 * ee1/ee135 + left * ((ee282 *
        ee162 * ee59/ee135 - ee461) * ee44 - ee466) * ee2) * ee59/ee135 +
        1) * ee1);
      out(j, 2) += -((ee475 + left * ((ee504 + ee162 * (2 * ee277 +
        left * ((ee541 - (2 * ee225 + 2 * ee160 + 4 * ee230)) * ee44 -
        2 * ee279) * ee2)) * ee59/ee135 + ee487 - ee491) * ee2) *
        ee59 * ee1/ee135);
      out(j, 3) += -(ee2 * (left * (ee477 + left * (ee531 * ee162 *
        ee59/ee135 + ee492 + ee493) * ee2) * ee59/ee135 - y));
      out(j, 4) += -(((((ee282 * (ee525 - ee535) + 2 * (ee522 + R_pow(ee282, 2)) +
        2 * ee522) * ee376 + ee459 * ((2 * ((ee524 *
        ee44 - ee474 * ee376) * ee376) + 4 * (ee524 * ee518) + 8 *
        (ee474 * R_pow(ee3, (4 * ee1))))/ee135 + 4 * (ee282 * ee518)) *
        ee59 * ee1/ee135) * ee59 * ee1/ee135 + ((ee55 * (ee536 *
        ee52 + 1) - (4 * ee420 + 4 * (ee116 * ee404) + 4 * ee424 +
        4 * ee449 + 6 * (ee235 * ee159)) * ee1) * ee52 - (((((((((((ee285 +
        (6 - 64 * ee142) * ee31 + 4 * (ee381 + 4 * ee252) +
        8 * (ee381 + 2 * ee237)) * ee31 + 
        ee297 + 9 * ee61) *
        ee1 + 7 * ee31 - (2 * (ee324 + (1 + 3 * (ee237 * ee1)) * ee31) +
        4 * ee508)) * ee31 + ((3 * (ee1 * ee266 + ee61) + ee310) *
        ee31 + 3 * (ee237 * ee61) + 6 * ee266 + ee1 * Rf_psigamma(ee6, 3)) *
        ee1 + 7 * ee61 - ee237 * (6 * ee338 - ee412)) * ee1 +
        (1 - 2 * (ee508 * ee1)) * ee31)/ee19 + (((((ee287 + (6 -
        64 * ee144) * ee33 + 4 * (ee383 + 4 * ee254) + 8 * (ee383 +
        2 * ee239)) * ee33 + ee299 + 9 * ee63) * ee1 + 7 * ee33 - (2 *
        (ee326 + (1 + 3 * (ee239 * ee1)) * ee33) + 4 * ee510)) *
        
        ee33 + ((3 * (ee1 * ee268 + ee63) + ee312) * ee33 + 3 *
        (ee239 * ee63) + 6 * ee268 + ee1 * Rf_psigamma(ee8, 3)) * ee1 +
        7 * ee63 - ee239 * (6 * ee342 - ee414)) * ee1 + (1 - 2 *
        (ee510 * ee1)) * ee33) * ee4/ee21 + left * (((((ee286 + (6 -
        64 * ee143) * ee32 + 4 * (ee382 + 4 * ee253) + 8 * (ee382 +
        2 * ee238)) * ee32 + ee298 + 9 * ee62) * ee1 + 7 * ee32 - (2 *
        (ee325 + (1 + 3 * (ee238 * ee1)) * ee32) + 4 * ee509)) *
        ee32 + ((3 * (ee1 * ee267 + ee62) + ee311) * ee32 + 3 * (ee238 *
        ee62) + 6 * ee267 + 
        ee1 * Rf_psigamma(ee7, 3)) * ee1 +
        7 * ee62 - ee238 * (6 * ee340 - ee413)) * ee1 + (1 - 2 *
        (ee509 * ee1)) * ee32) * ee2/ee20) * ee4 + (((((ee295 + (6 -
        64 * ee151) * ee40 + 4 * (ee390 + 4 * ee261) + 8 * (ee390 +
        2 * ee246)) * ee40 + ee308 + 9 * ee71) * ee1 + 7 * ee40 - (2 *
        (ee334 + (1 + 3 * (ee246 * ee1)) * ee40) + 4 * ee517)) *
        ee40 + ((3 * (ee1 * ee276 + ee71) + ee320) * ee40 + 3 * (ee246 *
        ee71) + 6 * ee276 + ee1 * Rf_psigamma(ee15, 3)) * ee1 + 7 *
        ee71 - ee246 * (6 * ee358 - ee419)) * ee1 + (1 - 
        2 *
        (ee517 * ee1)) * ee40)/ee29 + left * (((((ee284 + (6 - 64 *
        ee141) * ee30 + 4 * (ee380 + 4 * ee251) + 8 * (ee380 + 2 *
        ee236)) * ee30 + ee296 + 9 * ee60) * ee1 + 7 * ee30 - (2 * (ee323 +
        (1 + 3 * (ee236 * ee1)) * ee30) + 4 * ee507)) * ee30 +
        ((3 * (ee1 * ee265 + ee60) + ee309) * ee30 + 3 * (ee236 *
        ee60) + 6 * ee265 + ee1 * Rf_psigamma(ee5, 3)) * ee1 + 7 * ee60 -
        ee236 * (6 * ee336 - ee411)) * ee1 + (1 - 2 * (ee507 * ee1)) *
        ee30) * ee2/ee18) * ee4 + (((((ee293 + (6 - 64 * ee149) *
        ee38 + 4 * 
        (ee388 + 4 * ee259) + 8 * (ee388 + 2 *
        ee244)) * ee38 + ee306 + 9 * ee69) * ee1 + 7 * ee38 - (2 * (ee332 +
        (1 + 3 * (ee244 * ee1)) * ee38) + 4 * ee515)) * ee38 +
        ((3 * (ee1 * ee274 + ee69) + ee318) * ee38 + 3 * (ee244 *
        ee69) + 6 * ee274 + ee1 * Rf_psigamma(ee13, 3)) * ee1 + 7 * ee69 -
        ee244 * (6 * ee354 - ee417)) * ee1 + (1 - 2 * (ee515 *
        ee1)) * ee38)/ee27 + left * (((((ee294 + (6 - 64 * ee150) *
        ee39 + 4 * (ee389 + 4 * ee260) + 8 * (ee389 + 2 * ee245)) *
        ee39 + ee307 + 9 * ee70) * ee1 + 7 * ee39 - 
        (2 * (ee333 +
        (1 + 3 * (ee245 * ee1)) * ee39) + 4 * ee516)) * ee39 + ((3 *
        (ee1 * ee275 + ee70) + ee319) * ee39 + 3 * (ee245 * ee70) +
        6 * ee275 + ee1 * Rf_psigamma(ee14, 3)) * ee1 + 7 * ee70 -
        ee245 * (6 * ee356 - ee418)) * ee1 + (1 - 2 * (ee516 * ee1)) *
        ee39) * ee2/ee28) * ee4 + (((((ee291 + (6 - 64 * ee147) *
        ee36 + 4 * (ee386 + 4 * ee257) + 8 * (ee386 + 2 * ee242)) *
        ee36 + ee304 + 9 * ee67) * ee1 + 7 * ee36 - (2 * (ee330 + (1 +
        3 * (ee242 * ee1)) * ee36) + ee545)) * ee36 + ((3 * (ee1 *
        ee272 + 
        ee67) + ee316) * ee36 + 3 * (ee242 * ee67) +
        6 * ee272 + ee1 * Rf_psigamma(ee11, 3)) * ee1 + 7 * ee67 - ee242 *
        (6 * ee350 - ee264)) * ee1 + (1 - 2 * (ee513 * ee1)) *
        ee36)/ee25 + left * (((((ee292 + (6 - 64 * ee148) * ee37 + 4 *
        (ee387 + 4 * ee258) + 8 * (ee387 + 2 * ee243)) * ee37 + ee305 +
        9 * ee68) * ee1 + 7 * ee37 - (2 * (ee331 + (1 + 3 * (ee243 *
        ee1)) * ee37) + 4 * ee514)) * ee37 + ((3 * (ee1 * ee273 +
        ee68) + ee317) * ee37 + 3 * (ee243 * ee68) + 6 * ee273 +
        ee1 * Rf_psigamma(ee12, 3)) * ee1 + 
        7 * ee68 - ee243 *
        (6 * ee352 - ee416)) * ee1 + (1 - 2 * (ee514 * ee1)) * ee37) *
        ee2/ee26) * ee4 + (((((ee289 + (6 - 64 * ee145) * ee34 +
        4 * (ee384 + ee409) + 8 * (ee384 + ee405)) * ee34 + ee302 +
        9 * ee65) * ee1 + 7 * ee34 - (2 * (ee328 + (1 + 3 * (ee240 *
        ee1)) * ee34) + 4 * ee511)) * ee34 + ((3 * (ee1 * ee270 + ee65) +
        ee314) * ee34 + 3 * (ee240 * ee65) + 6 * ee270 + ee1 *
        Rf_psigamma(ee9, 3)) * ee1 + 7 * ee65 - ee240 * (6 * ee346 -
        ee394)) * ee1 + (1 - 2 * (ee511 * ee1)) * ee34)/ee23 + left *
        
        (((((ee290 + (6 - 64 * ee146) * ee35 + 4 * (ee385 + 4 *
        ee256) + 8 * (ee385 + 2 * ee241)) * ee35 + ee303 + 9 * ee66) *
        ee1 + 7 * ee35 - (2 * (ee329 + (1 + 3 * (ee241 * ee1)) *
        ee35) + 4 * ee512)) * ee35 + ((3 * (ee1 * ee271 + ee66) +
        ee315) * ee35 + 3 * (ee241 * ee66) + 6 * ee271 + ee1 * Rf_psigamma(ee10, 3)) *
        ee1 + 7 * ee66 - ee241 * (ee546 - ee415)) *
        ee1 + (1 - 2 * (ee512 * ee1)) * ee35) * ee2/ee24) * ee4 + (((((ee360 +
        (6 - 64 * ee152) * ee43 + 4 * (ee379 + 4 * ee396) +
        8 * (ee379 + 2 * ee362)) * 
        ee43 + ee361 + 9 * ee87) *
        ee1 + 7 * ee43 - (2 * (ee366 + (1 + 3 * (ee362 * ee1)) *
        ee43) + 4 * ee533)) * ee43 + ((3 * (ee1 * ee281 + ee87) + ee364) *
        ee43 + 3 * (ee362 * ee87) + 6 * ee281 + ee1 * Rf_psigamma(ee17, 3)) *
        ee1 + 7 * ee87 - ee362 * (6 * ee369 - ee410)) *
        ee1 + (1 - 2 * (ee533 * ee1)) * ee43)/ee42 + left * (((((ee288 +
        (6 - 64 * ee139) * ee41 + 4 * (ee370 + 4 * ee373) + 8 *
        (ee370 + 2 * ee301)) * ee41 + ee300 + 9 * ee64) * ee1 + 7 *
        ee41 - (2 * (ee327 + (1 + 3 * (ee301 * ee1)) * ee41) + 
        4 *
        ee534)) * ee41 + ((3 * (ee1 * ee269 + ee64) + ee313) *
        ee41 + 3 * (ee301 * ee64) + 6 * ee269 + ee1 * Rf_psigamma(ee16, 3)) *
        ee1 + 7 * ee64 - ee301 * (6 * ee344 - ee393)) * ee1 +
        (1 - 2 * (ee534 * ee1)) * ee41) * ee2/ee22)) * ee44) * ee59/
          ee135 + ee551 + p2 - (((6 * ee552 + ee1 * Rf_psigamma(ee1, 3)) *
            ee1 + 7 * ee553) * ee1 + ee548)) * ee1);
      out(j, 5) += -(((((ee521 + (ee277 + left * ((ee541 - ee233) *
        ee44 - ee279) * ee2) * ee282 + 2 * (ee521 + ee322 * ee282)) *
        ee44 + (ee165 * ((2 * (ee500 - ee496 * ee44) + 4 * ee500) *
        ee376 + 8 * (ee496 * ee518)) * ee1/ee135 + 4 * (left * ee282 *
        ee162 * ee2 * ee376)) * ee165 * ee59/ee135) * ee59 *
        ee1/ee135 + (ee55 * ((ee536 - left * ee404 * ee2) * ee52 + 1) +
        ee163 * ee404 * ee2 * ee52 - (3 * (ee235 * ee227) + 3 *
        (ee116 * ee432) + left * (3 * ee421 + 3 * (ee233 * ee159) +
        3 * ee401 + 3 * ee456) * ee2 * ee52) * ee1) * ee44 + left *
        ((ee470 * ee162 * ee59/ee135 - (((((((11 * ee509 - (22 * ee340 +
        44 * ee238 - 88 * ee253) * ee32 * ee1)/ee20 + left * (12 *
        ee510 - (24 * ee342 + 48 * ee239 - 96 * ee254) * ee33 *
        ee1) * ee2/ee21) * ee4 + (9 * ee507 - (18 * ee336 + 36 * ee236 -
        72 * ee251) * ee30 * ee1)/ee18 + left * (10 * ee508 - (20 *
        ee338 + 40 * ee237 -   80 * ee252) * ee31 * ee1) * ee2/
          ee19) * ee4 + (7 * ee516 - (14 * ee356 + 28 * ee245 - 56 *   ee260) *
            ee39 * ee1)/ee28 + left * (8 * ee517 - (16 * ee358 +
            32 * ee246 - 64 * ee261) * ee40 *   ee1) * ee2/ee29) * ee4 +
            (5 * ee514 - (10 * ee352 + 20 * ee243 - 40 * ee258) * ee37 *
            ee1)/ee26 + left * (6 * ee515 - (12 * ee354 + 24 * ee244 -
            48 * ee259) * ee38 * ee1) * ee2/ee27) * ee4 + (3 * ee512 -
            (12 * ee241 + ee546 - 24 * ee256) *   ee35 * ee1)/ee24 +
            left * (ee545 - (ee539 + 8 * ee350 - ee544) * ee36 * ee1) *
            ee2/ee25) * ee4 + ee400/ee22 + left * (2 * ee511 - (4 * ee346 +
            8 * ee240 - 16 * ee255) * ee34 * ee1) * ee2/ee23)) * ee44 -
            ee458 * ee112) * ee2) * ee59/ee135 + 1) * ee1);
      out(j, 6) += -((ee55 * (((ee402 + left * (ee367 + 3 * ee56) *
        ee1 * ee2) * ee52 + ((ee17 - ee321) * ee52 + 2) * ee44 + ee550) *
        ee1 + (1 - left * ((ee519 * ee159 + 2 * ee262) * ee52 +
        2 + ee520) * ee2) * ee44) + ((2 * (R_pow(ee322, 2) + ee529 *
        ee165 * ee44) + left * ((ee162 * (4 * (ee500/ee44) + 8 *
        ee496) + 2 * (ee505 * ee165)) * ee376/ee135 + 4 * (ee499 *
        ee44)) * ee165 * ee59 * ee2/ee135) * ee59/ee135 - 2 * (ee116 *
        ee437 + left * (ee543 * ee52 + 2 * (ee233 * ee227)) * ee2 *
        ee44)) * ee1 + left * ((ee282 * 
        ee532 * ee44 + 2 *
        (ee502 * ee162)) * ee59/ee135 + (ee159 * ee172 * ee52 + 2 *
        (ee110 * ee432) - (ee235 * ee170/ee3 + ee401 + ee138 * ee407/
          ee24 + left * ((((((100 * ee237 - 200 * ee252)/ee19 + (144 *
            ee239 - 288 * ee254) * ee4/ee21 + left * (121 * ee238 - 242 *
            ee253) * ee2/ee20) * ee4 + (64 * ee246 - 128 * ee261)/ee29 +
            left * (81 * ee236 - 162 * ee251) * ee2/ee18) * ee4 + (36 *
            ee244 - 72 * ee259)/ee27 + left * (49 * ee245 - 98 * ee260) *
            ee2/ee28) * ee4 + (ee539 - ee544)/ee25 + left * (25 *
            
            ee243 - 50 * ee258) * ee2/ee26) * ee4 + (ee408 - ee394)/
              ee23) * ee2)) * ee44 - 2 * (left * ee461 * ee112 * ee2)) *
                ee2) * ee59 * ee1/ee135);
      out(j, 7) += -(left * ((ee322 * ee531 + (ee475 + 2 * ee529 +
        left * (((ee165 * (ee540 + ee547) * ee44 + 2 * (ee523 * ee135))/
          ee135 + 2 * ee523) * ee59/ee135 + ee487 - ee491) * ee2) *
            ee162 + ee506 * ee165 * ee44) * ee59/ee135 + ee55 * ((((3 +
            left * (ee140 - 2 * ee520) * ee2) * ee1 + 1 - ee118 * (3 -
            2 * ee321) * ee117/ee4)/ee3 + ee321 - (ee538 * ee52 + 1 +
            2 * ee378)) * ee44 + ((2 + left * (ee140 - ee520) * ee2) * ee56 +
            2 * ((ee359 * ee52 + ee1) * ee247) - 3 * ee550) * ee1 -
            3 * ((ee44 + ee399 * ee56 * ee52) * ee519)) + (ee441 * ee52 +
            3 * (ee227 * ee172) - (3 * (ee377 * ee35/ee24) + 4 * (ee395 *
            ee36/ee25) + ee160 + left * (((((2662 * ee192 + 3456 *
            ee218) * ee4 + 1458 * ee191 + 2000 * ee215) * ee4 + 1024 *
            ee224 + 686 * ee194) * ee4 + 250 * ee193 + 432 * ee221) * ee4 +
            8 * ee88) * ee2)) * ee44 + 3 * (ee110 * ee437) - (ee116 *
            ee453 + left * (3 * (ee233 * ee170 * ee56) + 3 * (ee543 *
            ee112)) * ee2)) * ee59 * ee1 * ee2/ee135);
      out(j, 8) += -(ee2 * (left * (ee55 * (((2 * ((2 + ee399) * ee50) +
        2 * (left * ee359 * ee2)) * R_pow(ee3, (ee1 - 3)) - left *
        ((1 + 2 * ee538 + left * (ee359/ee3 + ee3 - 3) * ee2) *
        ee56 + 6 * (ee363 * ee519) + left * (ee283 + 2 + 2 * ee359 -
        ee375) * ee2 * ee247) * ee2) * ee1 - (1 - left * (7 - left *
        (6 - ee3) * ee2) * ee2) * ee44) + ((((10 * ((3 * (1638 +
        81 * ee3) + 400 * ee3)/ee19) + 11 * (left * (3 * (100 * ee3 +
        2220) + 484 * ee3) * ee2/ee20) + 12 * ((3 * (121 * ee3 +
        2926) + 576 * ee3) * ee4/ee21)) * 
        ee4 + 8 * ((256 * ee3 +
        3 * (49 * ee3 + 798))/ee29) + 9 * (left * (3 * (1168 + 64 *
        ee3) + 324 * ee3) * ee2/ee18)) * ee4 + 6 * ((144 * ee3 +
        3 * (25 * ee3 + 310))/ee27) + 7 * (left * (196 * ee3 + 3 * (36 *
        ee3 + 516)) * ee2/ee28)) * ee4 + ee57 + 16 * ee107 + 3 *
        ((ee4 + 26 * ee119)/ee24) + 4 * ((ee372 + left * (3 * ee377 +
        6 * ee4 + 9 * ee138) * ee2)/ee25) + 5 * ((100 * ee4 + 2 *
        ((16 * ee3 + 40) * ee4 + 4 * ee395))/ee26)) * ee44 + left *
        ((ee278 * ee532 + ee162 * (2 * ee506 + ee462 * ((ee540 + 
        4 *
        ee468 + ee547)/ee135 + 4 * ee278) * ee59 * ee2/ee135) +
        2 * (R_pow(ee278, 2) + ee162 * ee506)) * ee59/ee135 + 4 *
        (ee110 * ee453) + 4 * (ee441 * ee112) + 6 * (ee170 * ee172 *
        ee56)) * ee2) * ee59/ee135 - y));
      
      
    }
    
  }
  
  return out;
  
}

// //' Left-truncated fixed-alpha Gamma distribution negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each GEV parameter
// //' @param X1 a design matrix for the shape parameter
// //' @param ymat a matrix
// //' @param leftmat a matrix
// //' @param alphamat a matrix
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return ltgammabd0 a scalar, the negative log-likelihood
// //' @return ltgammab12 a matrix, first then second derivatives w.r.t. parameters
// //' @return ltgammab34 a matrix, third then fourth derivatives w.r.t. parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double ltgammabd0(Rcpp::List pars, arma::mat X1, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat leftmat, arma::mat alphamat)
{
  
  int nobs = nhere.size();
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
  }
  
  double y, left, p1, alpha, beta, lnum;
  
  double nllh = 0.0;
  
  for (int j=0; j < nobs; j++) {
    
    p1 = p1vec[j];
    beta = exp(p1);
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      left = leftmat(j, l);
      alpha = alphamat(j, l);
      
      lnum = (alpha - 1) * log(y) - beta * y + alpha * p1 - lgamma(alpha);
      nllh += ldenom(left, alpha, beta) - lnum;
      
    }
    
  }
  
  return(nllh);
  
}

// //' @rdname ltgammabd0
// [[Rcpp::export]]
arma::mat ltgammabd12(Rcpp::List pars, arma::mat X1, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat leftmat, arma::mat alphamat)
{
  
  int nobs = nhere.size();
  arma::mat out = arma::mat(nobs, 2, arma::fill::zeros);
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
  }
  
  double y, left, alpha, p1;
  
  double ee1, ee2, ee3;
  double ee16, ee17, ee18, ee19;
  double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
  double ee30, ee31, ee32, ee33, ee34, ee38, ee39;
  double ee40, ee41, ee42, ee47, ee48;
  double ee50, ee52;
  
  for (int j=0; j < nobs; j++) {
    
    p1 = p1vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      left = leftmat(j, l);
      alpha = alphamat(j, l);
      
      ee1 = exp(p1);
      ee2 = left * ee1;
      ee3 = R_pow(ee2, 2);
      ee16 = tgamma(10 + alpha);
      ee17 = tgamma(11 + alpha);
      ee18 = tgamma(12 + alpha);
      ee19 = tgamma(13 + alpha);
      ee20 = tgamma(2 + alpha);
      ee21 = tgamma(3 + alpha);
      ee22 = tgamma(4 + alpha);
      ee23 = tgamma(5 + alpha);
      ee24 = tgamma(6 + alpha);
      ee25 = tgamma(7 + alpha);
      ee26 = tgamma(8 + alpha);
      ee27 = tgamma(9 + alpha);
      ee28 = R_pow(ee2, alpha);
      ee29 = ee3/ee19;
      ee30 = 1 + alpha;
      ee31 = ee2/ee16;
      ee32 = ee2/ee18;
      ee33 = ee2/ee24;
      ee34 = ee2/ee26;
      ee38 = (((((ee29 + 1/ee17 + ee32) * ee3 + 1/ee27 + ee31) *  ee3 + 1/ee25 + ee34) * ee3 + 1/ee23 + ee33) * ee3 + 1/ee21 +  ee2/ee22) * ee3 + 1/tgamma(ee30) + ee2/ee20;
      ee39 = exp(-ee2);
      ee40 = alpha - 1;
      ee41 = R_pow(ee2, ee40);
      ee42 = 1/ee20;
      ee47 = (((((11/ee18 + 12 * (ee2/ee19)) * ee3 + 10 * (ee2/ee17) +  9/ee16) * ee3 + 7/ee26 + 8 * (ee2/ee27)) * ee3 + 5/ee24 +  6 * (ee2/ee25)) * ee3 + 3/ee22 + 4 * (ee2/ee23)) * ee3 +  ee42 + 2 * (ee2/ee21);
      ee48 = 1 - ee38 * ee39 * ee28;
      ee50 = alpha * ee41 - ee28;
      ee52 = ee38 * ee50 + ee47 * ee28;
      
      out(j, 0) += -(alpha + ee1 * (left * ee52 * ee39/ee48 - y));
      out(j, 1) += -(ee1 * (left * (ee38 * (alpha * (1 + left * (ee40/
        ee2 - 1) * ee1) - left * (ee30 - ee2) * ee1) * ee41 + (ee42 +
          3 * ((ee3 + 2 * (R_pow(left, 2) * R_pow(ee1, 2)))/ee22) +
          left * (((((100/ee17 + 121 * ee32 + 144 * ee29) * ee3 +
          64/ee27 + 81 * ee31) * ee3 + 36/ee25 + 49 * ee34) * ee3 + 16/
            ee23 + 25 * ee33) * ee3 + 4/ee21) * ee1) * ee28 + left * (R_pow(ee52, 2) *
              ee39/ee48 + 2 * (ee47 * ee50)) * ee1) * ee39/
                ee48 - y));
      
    }
    
  }
  
  return out;
  
}

// //' @rdname ltgammabd0
// [[Rcpp::export]]
arma::mat ltgammabd34(Rcpp::List pars, arma::mat X1, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat leftmat, arma::mat alphamat)
{
  
  int nobs = nhere.size();
  arma::mat out = arma::mat(nobs, 2, arma::fill::zeros);
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
  }
  
  double y, left, alpha, p1;
  
  double ee1, ee2, ee3, ee4;
  double ee17, ee18, ee19;
  double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
  double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee37;
  double ee40, ee42, ee43, ee44, ee45, ee47, ee48;
  double ee51, ee52, ee54, ee56, ee57;
  double ee60, ee67, ee69;
  double ee72, ee73, ee76, ee78, ee79;
  double ee80, ee81, ee82, ee84, ee86, ee87;
  double ee90, ee91, ee92;
  double ee102, ee106, ee108;
  double ee112, ee113, ee114, ee117;
  double ee123, ee125;
  
  
  for (int j=0; j < nobs; j++) {
    
    p1 = p1vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      left = leftmat(j, l);
      alpha = alphamat(j, l);
      
      ee1 = exp(p1);
      ee2 = left * ee1;
      ee3 = R_pow(ee2, 2);
      ee4 = alpha - 1;
      ee17 = tgamma(10 + alpha);
      ee18 = tgamma(11 + alpha);
      ee19 = tgamma(12 + alpha);
      ee20 = tgamma(13 + alpha);
      ee21 = tgamma(2 + alpha);
      ee22 = tgamma(3 + alpha);
      ee23 = tgamma(4 + alpha);
      ee24 = tgamma(5 + alpha);
      ee25 = tgamma(6 + alpha);
      ee26 = tgamma(7 + alpha);
      ee27 = tgamma(8 + alpha);
      ee28 = tgamma(9 + alpha);
      ee29 = R_pow(ee2, alpha);
      ee30 = 1 + alpha;
      ee31 = R_pow(ee2, ee4);
      ee32 = ee3/ee20;
      ee33 = ee2/ee17;
      ee34 = ee2/ee19;
      ee35 = ee2/ee25;
      ee36 = ee2/ee27;
      ee37 = 1/ee21;
      ee40 = (((((ee32 + 1/ee18 + ee34) * ee3 + 1/ee28 + ee33) *  ee3 + 1/ee26 + ee36) * ee3 + 1/ee24 + ee35) * ee3 + 1/ee22 +  ee2/ee23) * ee3 + 1/tgamma(ee30) + ee2/ee21;
      ee42 = ee2/ee18;
      ee43 = ee2/ee20;
      ee44 = ee2/ee26;
      ee45 = ee2/ee28;
      ee47 = R_pow(left, 2) * R_pow(ee1, 2);
      ee48 = ee2/ee22;
      ee51 = (((((11/ee19 + 12 * ee43) * ee3 + 10 * ee42 + 9/ee17) *  ee3 + 7/ee27 + 8 * ee45) * ee3 + 5/ee25 + 6 * ee44) *  ee3 + 3/ee23 + 4 * (ee2/ee24)) * ee3 + ee37 + 2 * ee48;
      ee52 = ee4/ee2;
      ee54 = alpha * ee31 - ee29;
      ee56 = exp(-ee2);
      ee57 = ee3 + 2 * ee47;
      ee60 = left * (ee52 - 1) * ee1;
      ee67 = ee37 + 3 * (ee57/ee23) + left * (((((100/ee18 + 121 *  ee34 + 144 * ee32) * ee3 + 64/ee28 + 81 * ee33) * ee3 +  36/ee26 + 49 * ee36) * ee3 + 16/ee24 + 25 * ee35) * ee3 +  4/ee22) * ee1;
      ee69 = alpha * (1 + ee60) - left * (ee30 - ee2) * ee1;
      ee72 = 1 - ee40 * ee56 * ee29;
      ee73 = alpha * left;
      ee76 = ee40 * ee69 * ee31 + ee67 * ee29;
      ee78 = ee40 * ee54 + ee51 * ee29;
      ee79 = ee30 * ee4;
      ee80 = ee79 + 1;
      ee81 = ee76 + 2 * (left * ee51 * ee54 * ee1);
      ee82 = R_pow(ee2, (alpha - 2));
      ee84 = 8 * ee3;
      ee86 = left * (3 - ee2) * ee1;
      ee87 = ee3 + 8 * ee47;
      ee90 = ee84 + left * (3 * ee57 + 3 * ee3) * ee1;
      ee91 = ee73 * ee1;
      ee92 = R_pow(ee78, 2);
      ee102 = ee37 + 3 * (ee87/ee23) + 4 * (ee90/ee24) + left *  (((((2662/ee19 + 3456 * ee43) * ee3 + 1458/ee17 + 2000 *  ee42) * ee3 + 1024 * ee45 + 686/ee27) * ee3 + 250/ee25 +  432 * ee44) * ee3 + 8/ee22) * ee1;
      ee106 = ee73 * ((ee80 - ee91) * ee82 - (2 + ee60) * ee31) *  ee1 - 2 * ((1 + alpha * (3 + left * (ee52 - 2) * ee1) -  ee86) * ee31);
      ee108 = ee40 * ee106 + ee102 * ee29;
      ee112 = 3 * (ee51 * ee69 * ee31);
      ee113 = 3 * (ee67 * ee54);
      ee114 = ee81 * ee72;
      ee117 = ee108 + left * (ee112 + ee113) * ee1;
      ee123 = left * ee92 * ee56 * ee1;
      ee125 = left * (2 * (ee92 * ee56/ee72) + 2 * (ee51 * ee54)) *  ee1;
      
      out(j, 0) += -(ee1 * (left * (ee108 + left * ((ee76 + 2 * ee81 +
        ee125) * ee78 * ee56/ee72 + ee112 + ee113) * ee1) * ee56/
          ee72 - y));
      out(j, 1) += -(ee1 * (left * (ee40 * (alpha * ((2 * ((2 + ee91) *
        ee4) + 2 * (left * ee80 * ee1)) * R_pow(ee2, (alpha -
        3)) - left * ((1 + 2 * (1 - ee86) + left * (ee80/ee2 + ee2 -
        3) * ee1) * ee31 + 6 * ((ee31 + left * ee4 * ee1 * ee82) *
        (1 - ee2)) + left * (ee79 + 2 + 2 * ee80 - ee86) * ee1 * ee82) *
        ee1) - (1 - left * (7 - left * (6 - ee2) * ee1) * ee1) *
        ee29) + ((((10 * ((3 * (1638 + 81 * ee2) + 400 * ee2)/ee18) +
        11 * (left * (3 * (100 * ee2 + 2220) + 484 * ee2) * ee1/
          ee19) + 12 * ((3 * (121 * ee2 + 2926) + 
            576 * ee2) * ee3/
              ee20)) * ee3 + 8 * ((256 * ee2 + 3 * (49 * ee2 + 798))/ee28) +
                9 * (left * (3 * (1168 + 64 * ee2) + 324 * ee2) * ee1/
                  ee17)) * ee3 + 6 * ((144 * ee2 + 3 * (25 * ee2 + 310))/ee26) +
                    7 * (left * (196 * ee2 + 3 * (36 * ee2 + 516)) * ee1/ee27)) *
                    ee3 + ee37 + 16 * ee48 + 3 * ((ee3 + 26 * ee47)/ee23) +
                    4 * ((ee84 + left * (3 * ee87 + 6 * ee3 + 9 * ee57) * ee1)/
                      ee24) + 5 * ((100 * ee3 + 2 * ((16 * ee2 + 40) * ee3 + 4 * ee90))/
                        ee25)) * ee29 + left * ((ee81 * (ee76 + ee125) + ee78 *
                          (2 * 
                          ee117 + left * ee78 * ((2 * (ee114 - ee123) + 4 *
                          ee114 + 8 * ee123)/ee72 + 4 * ee81) * ee56 * ee1/ee72) + 2 *
                          (R_pow(ee81, 2) + ee78 * ee117)) * ee56/ee72 + 4 * (ee51 *
                          ee106) + 4 * (ee102 * ee54) + 6 * (ee67 * ee69 * ee31)) * ee1) *
                          ee56/ee72 - y));
      
      
    }
    
  }
  
  return out;
  
}
