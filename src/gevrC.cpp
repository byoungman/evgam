// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "shared.h"

// //' Right-censored generalized extreme value (GEV) distribution negative 
// //' log-likelihood with constrained shape parameter
// //'
// //' @param pars a list of vectors of coefficients for each GEV parameter
// //' @param X1 a design matrix for the GEV location parameter
// //' @param X2 a design matrix for the GEV log scale parameter
// //' @param X3 a design matrix for the GEV transformed shape parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return gevd0 a scalar, the negative log-likelihood
// //' @return gevd12 a matrix, first then second derivatives w.r.t. GEV parameters
// //' @return gevd34 a matrix, third then fourth derivatives w.r.t. GEV parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double gevrd0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = nhere.size();

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    txivec = txivec.elem(dupid);
}

double y, mu, lpsi, txi, xi;
double ee1, ee2;
double nllh = 0.0;

for (int j=0; j < nobs; j++) {

mu = muvec[j];
lpsi = lpsivec[j];
txi = txivec[j];
xi = 1.5 / (1.0 + exp(-txi)) - 1.0;

xi = xi_from_zero(xi, xieps);  
txi = xi2txi(xi);

for (int l=0; l < nhere[j]; l++) {

y = ymat(j, l);

ee1 = xi * (y - mu) / exp(lpsi);

if (ee1 <= -1.0) {
    nllh = 1e20;
    break;
} else {

ee2 = 1.0 / xi;

nllh -= log(1.0 - exp(-R_pow(1.0 + ee1, -ee2)));

}

} 

}

return(nllh);

}

// //' @rdname gevrd0
// [[Rcpp::export]]
arma::mat gevrd12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = nhere.size();

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    txivec = txivec.elem(dupid);
}

double y, mu, lpsi, txi, xi;

arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);

double ee2, ee3, ee5, ee6, ee7, ee9;
double ee10, ee11, ee13, ee14, ee15, ee16, ee17, ee18;
double ee21, ee22, ee23, ee26, ee27, ee28, ee29;
double ee34, ee35, ee36, ee37;
double ee40, ee41, ee43, ee44, ee45, ee46, ee47;

for (int j=0; j < nobs; j++) {

mu = muvec[j];
lpsi = lpsivec[j];
txi = txivec[j];
xi = 1.5 / (1.0 + exp(-txi)) - 1.0;

xi = xi_from_zero(xi, xieps);  
txi = xi2txi(xi);

for (int l=0; l < nhere[j]; l++) {

y = ymat(j, l);

  ee2 = exp(-txi);
  ee3 = 1 + ee2;
  ee5 = 1.5/ee3 - 1;
  ee6 = exp(lpsi);
  ee7 = y - mu;
  ee9 = ee5 * ee7/ee6;
  ee10 = 1/ee5;
  ee11 = ee9 + 1;
  ee13 = exp(-R_pow(ee11, -ee10));
  ee14 = 1 + ee10;
  ee15 = R_pow(ee11, ee14);
  ee16 = 1 - ee13;
  ee17 = log1p(ee9);
  ee18 = R_pow(ee11, ee10);
  ee21 = 1.5 * (ee7/(ee15 * ee6));
  ee22 = R_pow(ee3, 2);
  ee23 = 2 * ee14;
  ee26 = 1.5 * (ee17/(ee18 * ee5)) - ee21;
  ee27 = ee16 * ee22;
  ee28 = ee26 * ee13;
  ee29 = ee10 + 2;
  ee34 = ee14 * ee7;
  ee35 = (1 + ee13/ee16 - ee14 * ee5/R_pow(ee11, (ee29 - ee23)))/R_pow(ee11, ee23);
  ee36 = ee28/ee16;
  ee37 = 1.5/ee18;
  ee40 = ((ee36 + (ee37 - 1.5) * ee17/ee5 - ee21)/ee5 + 1.5 *  (ee34/(ee11 * ee6)))/ee15 * ee13 * ee2;
  ee41 = (ee35 * ee7/ee6 + 1/ee15) * ee13;
  ee43 = ee15 * ee16 * ee6;
  ee44 = ee27 * ee5;
  ee45 = ee27 * ee6;
  ee46 = ee16 * ee6;
  ee47 = ee3 * ee5;
  
  out(j, 0) += -(ee13/ee43);
  out(j, 1) += -(ee13 * ee7/ee43);
  out(j, 2) += -(ee28 * ee2/ee44);
  out(j, 3) += ee35 * ee13/(ee16 * R_pow(ee6, 2));
  out(j, 4) += ee41/ee46;
  out(j, 5) += ee40/ee45;
  out(j, 6) += ee41 * ee7/ee46;
  out(j, 7) += ee40 * ee7/ee45;
  out(j, 8) += -(((((1.5 - ee37) * ee17/ee5 + ee21 - ee36) * ee26 *
    ee2/ee22 + (2.25 * (ee2 * ee7/(ee11 * ee22 * ee6)) - ((4.5/
      ee47 - 3) * ee2/ee3 + 1.5) * ee17)/ee18)/ee5 + (((2.25/
        ee47 - 3) * ee2/ee3 + 1.5)/ee15 - 1.5 * ((1.5 * (ee17/(ee15 *
          R_pow(ee5, 2))) - 1.5 * (ee34/(R_pow(ee11, ee29) * ee6))) *
          ee2/ee22)) * ee7/ee6) * ee13 * ee2/ee44);
}

}

return out;

}

// //' @rdname gevrd0
// [[Rcpp::export]]
arma::mat gevrd34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = nhere.size();

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    txivec = txivec.elem(dupid);
}

double y, mu, lpsi, txi, xi;

arma::mat out = arma::mat(nobs, 25, arma::fill::zeros);

double ee2, ee3, ee5, ee6, ee7, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
double ee21, ee22, ee24, ee25, ee26, ee27, ee28, ee29;
double ee30, ee31, ee33, ee36, ee37;
double ee40, ee42, ee43, ee44, ee45, ee47;
double ee50, ee51, ee52, ee53, ee54, ee55, ee56, ee57;
double ee60, ee63;
double ee70, ee72, ee73, ee75, ee77, ee78, ee79;
double ee80, ee81, ee82, ee85, ee86, ee87, ee88;
double ee90, ee91, ee92, ee93, ee94, ee95, ee96, ee97, ee98, ee99;
double ee100, ee105, ee107, ee108, ee109;
double ee111, ee113, ee114, ee115, ee116, ee117, ee118;
double ee120, ee121, ee122, ee123, ee124, ee125, ee126, ee127, ee128;
double ee131, ee132, ee133, ee134, ee137, ee138, ee139;
double ee142, ee144, ee146, ee147, ee148, ee149;
double ee152, ee156, ee157, ee158;
double ee161, ee162, ee164, ee165, ee167, ee169;
double ee170, ee172, ee175;
double ee182, ee183, ee184, ee189;
double ee191, ee192, ee193, ee195, ee196, ee198, ee199;
double ee200, ee204, ee205, ee206, ee207, ee208, ee209;
double ee210, ee211, ee212, ee214, ee215, ee216, ee217;
double ee220, ee221, ee222, ee223, ee226, ee227, ee228, ee229;
double ee232, ee235, ee237, ee238, ee239;
double ee240, ee241, ee243, ee244, ee245, ee246, ee247, ee248, ee249;
double ee250, ee251, ee252, ee253, ee254, ee255, ee258, ee259;
double ee260, ee262, ee264, ee267, ee268;
double ee270, ee272, ee273, ee274, ee278;
double ee281, ee282, ee285, ee286;
double ee291, ee292, ee293, ee295, ee297, ee298, ee299;
double ee300, ee301, ee302, ee304, ee306, ee308, ee309;
double ee310, ee314;
double ee324, ee325, ee329;
double ee331, ee332, ee333, ee334, ee337, ee339;
double ee345, ee346, ee347, ee348, ee349;
double ee350, ee351, ee352, ee353;
double ee360, ee364, ee366, ee367, ee368;
double ee370, ee372, ee374, ee375, ee376, ee377, ee378, ee379;
double ee380, ee381, ee382, ee383, ee384, ee386;
double ee393, ee394, ee395, ee397, ee398, ee399;
double ee400, ee401, ee402, ee404, ee405, ee407, ee409;
double ee410, ee411, ee412, ee413, ee414, ee415, ee416, ee417, ee418;

for (int j=0; j < nobs; j++) {

mu = muvec[j];
lpsi = lpsivec[j];
txi = txivec[j];
xi = 1.5 / (1.0 + exp(-txi)) - 1.0;

xi = xi_from_zero(xi, xieps);  
txi = xi2txi(xi);

  for (int l=0; l < nhere[j]; l++) {

  y = ymat(j, l);

    ee2 = exp(-txi);
    ee3 = 1 + ee2;
    ee5 = 1.5/ee3 - 1;
    ee6 = exp(lpsi);
    ee7 = y - mu;
    ee9 = ee5 * ee7/ee6;
    ee10 = 1/ee5;
    ee11 = ee9 + 1;
    ee12 = 1 + ee10;
    ee13 = log1p(ee9);
    ee14 = R_pow(ee11, ee12);
    ee15 = ee10 + 2;
    ee16 = R_pow(ee3, 2);
    ee17 = R_pow(ee11, ee15);
    ee18 = R_pow(ee11, ee10);
    ee21 = 1.5 * (ee7/(ee14 * ee6));
    ee22 = R_pow(ee5, 2);
    ee24 = exp(-R_pow(ee11, -ee10));
    ee25 = ee3 * ee5;
    ee26 = ee12 * ee7;
    ee27 = ee18 * ee5;
    ee28 = 1.5 * (ee13/ee27);
    ee29 = ee28 - ee21;
    ee30 = ee17 * ee6;
    ee31 = 2 * ee12;
    ee33 = 1.5 * (ee13/(ee14 * ee22));
    ee36 = ee33 - 1.5 * (ee26/ee30);
    ee37 = 1 - ee24;
    ee40 = (4.5/ee25 - 3) * ee2/ee3 + 1.5;
    ee42 = ee11 * ee16 * ee6;
    ee43 = ee12 * ee5;
    ee44 = R_pow(ee11, ee31);
    ee45 = 1.5/ee18;
    ee47 = ee2 * ee7/ee42;
    ee50 = 2.25 * ee47 - ee40 * ee13;
    ee51 = 3 * (ee2/ee3);
    ee52 = ee11 * ee6;
    ee53 = 1.5 - ee51;
    ee54 = ee10 + 3;
    ee55 = R_pow(ee11, ee54);
    ee56 = ee36 * ee2;
    ee57 = 2.25/ee25;
    ee60 = 1/ee44 - ee43/ee17;
    ee63 = (ee57 - 3) * ee2/ee3 + 1.5;
    ee70 = 1.5 * (ee13/(ee17 * ee22)) - 1.5 * (ee15 * ee7/(ee55 *  ee6));
    ee72 = (ee63/ee14 - 1.5 * (ee56/ee16)) * ee7/ee6;
    ee73 = ee50/ee18;
    ee75 = (ee45 - 1.5) * ee13/ee5;
    ee77 = 1.5 * (ee26/ee52);
    ee78 = (ee75 - ee21)/ee5;
    ee79 = 1/ee14;
    ee80 = ee78 + ee77;
    ee81 = ee60 * ee7;
    ee82 = ee81/ee6;
    ee85 = (1.5 - ee45) * ee13/ee5 + ee21;
    ee86 = ee2/ee16;
    ee87 = ee16 * ee22;
    ee88 = 1 + 2 * ee2;
    ee90 = (ee85 * ee29 * ee2/ee16 + ee73)/ee5 + ee72;
    ee91 = ee53 * ee5;
    ee92 = ee43 * ee7;
    ee93 = ee82 + ee79;
    ee94 = 2.25 * (ee2/ee87);
    ee95 = 3 * ee88;
    ee96 = 3 * ee3;
    ee97 = ee12 * ee53;
    ee98 = 12 * ee2;
    ee99 = ee12 * ee70;
    ee100 = ee97 + ee94;
    ee105 = ee5 * ee15 * ee7/ee52;
    ee107 = ee50/ee14 + 1.5 * (ee56 * ee13/ee16);
    ee108 = 1.5/ee5;
    ee109 = 3 * ee53;
    ee111 = (ee100/ee17 - 1.5 * (ee99 * ee2/ee16)) * ee7/ee6;
    ee113 = ee108 - 1.5 * ee12;
    ee114 = 2/ee5;
    ee115 = ee16 * ee5;
    ee116 = ee70 * ee5;
    ee117 = ee29 * ee24;
    ee118 = ee92/ee30;
    ee120 = 2.25 * ee86 - ee91;
    ee121 = ee95 + ee96;
    ee122 = ee12 * (3 - ee105);
    ee123 = 3 * ee91;
    ee124 = 3 * ee120;
    ee125 = ee92/ee52;
    ee126 = (ee124 - (27 * ee86 + ee123))/ee5;
    ee127 = ee114 + 3;
    ee128 = (ee116 + 1.5/ee17) * ee12;
    ee131 = (((ee126 - ee109)/ee5 + ee95 + ee96 - ee98)/ee3 +  3) * ee2/ee3 - 1.5;
    ee132 = ee111 + ee107/ee22;
    ee133 = R_pow(ee11, ee127);
    ee134 = ee113/ee17;
    ee137 = ee14 * ee5;
    ee138 = ee107/ee5;
    ee139 = ee37 * ee16;
    ee142 = (2.25 * (ee7/(ee11 * ee3 * ee6)) - 3) * ee2/ee3 +  1.5;
    ee144 = ee131 * ee13 + (1.5 * ee142 + 3 * ee40) * ee2 *  ee7/ee42;
    ee146 = 1.5 * ee40 + ee109;
    ee147 = ee72 + (ee73 + 1.5 * (ee29 * ee2 * ee13/ee115))/ee5;
    ee148 = ee36/ee14;
    ee149 = ee13/ee137;
    ee152 = ((ee146/ee5 + ee98 - ee121)/ee3 - 3) * ee2/ee3 +  1.5;
    ee156 = 2 * ee125;
    ee157 = 2 * (ee36 * ee29 * ee2/ee16);
    ee158 = ee144/ee18;
    ee161 = ee80 * ee37;
    ee162 = ee80/ee44;
    ee164 = ((2 * (ee63 * ee36) - 1.5 * ee132) * ee2/ee16 -  ee152/ee14) * ee7/ee6;
    ee165 = ee118 - ee79;
    ee167 = ee117/ee5;
    ee169 = 1.5 * (ee147 * ee13) + 2 * (ee29 * ee50);
    ee170 = ee24 * ee7;
    ee172 = ee90/ee14 + ee157;
    ee175 = (ee122 - 1) * ee5 * ee7/ee6;
    ee182 = (1 - ee122) * ee5 * ee7/ee6 + 1;
    ee183 = ee139 * ee6;
    ee184 = ee37 * ee5;
    ee189 = 1.5 * ee149;
    ee191 = 1.5 * (ee13/(ee55 * ee22)) - 1.5 * (ee54 * ee7/(R_pow(ee11, (ee10 + 4)) * ee6));
    ee192 = 2 - ee105;
    ee193 = 2 * (ee117/ee37);
    ee195 = (((4.5 - ee45) * ee13/ee5 + ee21) * ee29 * ee2/ee16 +  3 * ee73)/ee5;
    ee196 = ee175 - 1;
    ee198 = (ee148 + ee134 - (ee162 + ee128)) * ee7/ee6;
    ee199 = R_pow(ee11, (ee15 - ee31));
    ee200 = ee169/ee5;
    ee204 = R_pow(ee29, 2) * ee24 * ee2/ee115;
    ee205 = ee79 - ee118;
    ee206 = 3 * ee72;
    ee207 = ee44 * ee6;
    ee208 = ee93 * ee37;
    ee209 = ee37 * ee6;
    ee210 = ee43/ee199;
    ee211 = ee70 * ee2;
    ee212 = ee192/ee17;
    ee214 = (ee172 - ee138)/ee5 - ee111;
    ee215 = ee90 * ee37;
    ee216 = (ee195 + ee206) * ee29;
    ee217 = ee196/ee11;
    ee220 = ee80 * ee29 * ee2/ee16;
    ee221 = R_pow(ee11, (ee12 + ee31));
    ee222 = ee27 * ee15;
    ee223 = ee133 * ee6;
    ee226 = ee29/ee17;
    ee227 = 1.5 - ((ee121 - ee98)/ee3 + 3) * ee2/ee3;
    ee228 = 4 * ee161;
    ee229 = 8 * ee167;
    ee232 = ((ee85 - ee193) * ee29 * ee2/ee16 + ee73)/ee5 +  ee72;
    ee235 = ((ee53 * ee15 + ee94)/ee55 - 1.5 * (ee191 * ee15 *  ee2/ee16)) * ee7/ee6 + (ee50/ee17 + 1.5 * (ee211 * ee13/ee16))/ee22;
    ee237 = ((ee200 - ee216) * ee2/ee16 - ee158)/ee5 + ee164;
    ee238 = ee182/ee11;
    ee239 = (ee125 - 2)/ee14;
    ee240 = (ee226 - ee116) * ee12;
    ee241 = R_pow(ee11, (ee31 - ee15));
    ee243 = (ee212 + ee7/ee223) * ee12 * ee5;
    ee244 = ee29 * ee205;
    ee245 = ee60/ee14;
    ee246 = (3 - ee156)/ee14;
    ee247 = 2 - ee222;
    ee248 = 2 * (ee215 + ee204);
    ee249 = ee24/ee37;
    ee250 = (ee239 - ee82)/ee14;
    ee251 = ee217 - (ee82 + ee246) * ee7/ee6;
    ee252 = ee165 * ee29;
    ee253 = R_pow(ee11, (ee12 - ee31));
    ee254 = R_pow(ee11, (2 + ee31 + ee114));
    ee255 = R_pow(ee37, 2);
    ee258 = ee43 * ee247/ee133 - ee245;
    ee259 = (ee156 - 3)/ee14;
    ee260 = R_pow(ee6, 2);
    ee262 = ee250 + ee243;
    ee264 = (ee252 + ee189)/ee5 + ee198;
    ee267 = ee80/ee241;
    ee268 = ee198 + (ee189 - ee244)/ee5;
    ee270 = ee240 + ee148 + ee134;
    ee272 = (ee75 + ee193 - ee21)/ee5 + ee77;
    ee273 = ee11 * ee36;
    ee274 = ee221 * ee37;
    ee278 = ee93 * ee24/ee37;
    ee281 = (ee259 - ee82) * ee7/ee6 - ee238;
    ee282 = ee191 * ee5;
    ee285 = ee29 * ee60 * ee24/ee184;
    ee286 = ee117/ee184;
    ee291 = ee113 * ee15/ee55;
    ee292 = (ee108 - 1.5 * ee15)/ee55;
    ee293 = (2 * (ee161 - ee167) + ee228 + ee229) * ee24;
    ee295 = (2 * (ee208 - ee170/ee207) + 4 * ee208)/ee44 + 8 *  (ee170/(ee254 * ee6));
    ee297 = 1 + 2 * ee249 - ee210;
    ee298 = 2 * ee80;
    ee299 = 3 - ee5 * ee54 * ee7/ee52;
    ee300 = 4 * ee2;
    ee301 = 4.5 * ee12;
    ee302 = 4.5/ee5;
    ee304 = ee144/ee14;
    ee306 = ee262/ee14;
    ee308 = ee235 * ee5;
    ee309 = (((ee292 - ee282 * ee15) * ee12 + ee291) * ee7/ee6 +  2 * ee99) * ee5;
    ee310 = ((ee36/R_pow(ee11, (ee10 - ee31)) + ee113/ee253 - 1.5 * (ee26/ee6))/ee11 - (ee78 + (ee298 + 2 * ee286) * ee24/ee37))/ee44;
    ee314 = ee80 * (2 - ee210);
    ee324 = ee44 * ee37;
    ee325 = ee93 * ee36;
    ee329 = (ee15 * ee299 + ee10) * ee5 * ee7/ee6;
    ee331 = ((2 * (ee100 * ee70) - 1.5 * (ee235 * ee12)) * ee2/ee16 -  (ee12 * ee227 + ee146 * ee2/ee87)/ee17) * ee7/ee6;
    ee332 = ee295/ee37;
    ee333 = ((2 * ee138 - ee172)/ee5 + 2 * ee111 - 2 * (ee80 *  ee36 * ee2/ee16))/ee14;
    ee334 = ee139 * ee260;
    ee337 = ee297 * ee7/ee207 + ee79;
    ee339 = ee70 * ee113;
    ee345 = ee293 * ee7/(ee221 * ee255 * ee6);
    ee346 = (2 * (ee93/ee14) + 2 * (ee170/(ee274 * ee6))) *  ee24;
    ee347 = (4 * ee220 - (ee248 - ee29 * (ee228 + ee229) * ee2/ee16)/ee37) *  ee24;
    ee348 = 1.5 * (ee132 * ee13);
    ee349 = 2 * (R_pow(ee80, 2) * ee2/(ee44 * ee16) - ee214/ee14);
    ee350 = 2 * (ee36 * ee50);
    ee351 = 8 * ee204;
    ee352 = R_pow(ee6, 3);
    ee353 = ee13/ee5;
    ee360 = (((ee90 - (ee272 * ee29 * ee2/ee16 - ee90) * ee24/ee37)/ee14 +  ee29 * (2 * ee36 - ee80 * ee24/(ee14 * ee37)) * ee2/ee16 -  ee138)/ee5 - ee111) * ee24 * ee2;
    ee364 = ((ee237/ee14 + (ee29 * (3 * ee132 - ee214 * ee24/ee37) +  3 * (ee90 * ee36) - (ee214 * ee29 + ee232 * ee80/ee14) *  ee24/ee37) * ee2/ee16 - ((((((ee90 - ((2 * ((ee220 -  ee90) * ee37) + ee351 - ee248)/ee37 + 2 * (ee220 - ee90)) *  ee24/ee37)/ee14 + ee157 - ee138)/ee5 - ee111) * ee29 +  ee90 * ee80/ee14) * ee2/ee16 + (ee90 * ee80 * ee2/ee16 -  ee237)/ee14) * ee24/ee37 + (((ee126 + 1.5 * ee53 + 4.5 *  ee40)/ee14 + ee348 + ee350) * ee2/ee16 - ee304)/ee5))/ee5 -  ee331) * ee24 * ee2;
    ee366 = ee251/R_pow(ee11, (ee31 - ee127));
    ee367 = ee251/ee44;
    ee368 = ee232 * ee93;
    ee370 = (ee310 - ee128) * ee7/ee6;
    ee372 = ((ee273 + ee302 - (ee267 + ee301))/ee17 - 3 * ee128) *  ee7/ee6;
    ee374 = ((((ee57 + 3) * ee2/ee3 - 1.5)/ee5 + 2 * ee97)/ee17 +  ee333 + (ee347/(ee324 * ee5) + ee349) * ee24/ee37 + 2 *  (ee339 * ee2/ee16) - (ee308 + 1.5 * (ee211/ee16)) * ee12) *  ee7/ee6;
    ee375 = ee80 * ee165;
    ee376 = ee80 * ee93;
    ee377 = ee314/ee44;
    ee378 = ee272 * ee93;
    ee379 = (ee332 + 4 * (ee93/ee44)) * ee24;
    ee380 = ee281/ee44;
    ee381 = ee258 * ee24;
    ee382 = ee258/ee14;
    ee383 = ee165 * ee93;
    ee384 = ee337 * ee60;
    ee386 = (ee169 * ee2/ee115 - ee158)/ee5 + ee164;
    ee393 = ee278 + ee82;
    ee394 = ee278 + ee79;
    ee395 = R_pow(ee93, 2);
    ee397 = ((2 - ee329)/ee55 + ee212) * ee12 * ee5;
    ee398 = ((6 * (2 * ee88 + ee300) + 8 * ee121 - 96 * ee2)/ee3 +  12) * ee2;
    ee399 = ee139 * ee5;
    ee400 = ee37 * ee60;
    ee401 = ee37 * ee260;
    ee402 = ee37 * ee352;
    ee404 = ee227 * ee5;
    ee405 = ee134 - ee99 * ee5;
    ee407 = ee346 * ee7/ee209;
    ee409 = 1.5 * ee353;
    ee410 = 1/ee241;
    ee411 = 2 * ((ee134 - ee128) * ee7/ee6 + ee33);
    ee412 = 2 * ee285;
    ee413 = 2 * ee93;
    ee414 = 3 * (1 + ee300);
    ee415 = 3 * ee149;
    ee416 = 4 * (ee93 * ee29 * ee24/ee184);
    ee417 = 6 * ee3;
    ee418 = 9 * ee88;
    
    out(j, 0) += (((ee43 * (3 - ee222) - ee410)/ee17 - ee60 * ee24/
      ee37)/ee14 - (2 * ee245 + 2 * (ee24/ee274)) * ee24/ee37) *
        ee24/ee402;
    out(j, 1) += ((ee239 - (1 + ee249) * ee60 * ee7/ee6)/ee14 +
      ee243 - ee346/ee37) * ee24/ee401;
    out(j, 2) += (ee310 + ee240 - ee285) * ee24 * ee2/ee334;
    out(j, 3) += (((ee259 - ee393) * ee7/ee6 - ee238)/ee14 - ee407) *
      ee24/ee209;
    out(j, 4) += (ee370 + ((ee118 - ee394) * ee29 + ee189)/ee5) *
      ee24 * ee2/ee183;
    out(j, 5) += ee360/ee183;
    out(j, 6) += ((ee217 - (ee393 + ee246) * ee7/ee6)/ee14 - ee407) *
      ee24 * ee7/ee209;
    out(j, 7) += (ee370 + (ee189 - (ee394 - ee118) * ee29)/ee5) *
      ee24 * ee2 * ee7/ee183;
    out(j, 8) += ee360 * ee7/ee183;
    out(j, 9) += -((((ee200 - ((ee232 + 2 * ee90) * ee24/ee37 +
      ee195 + ee206) * ee29) * ee2/ee16 - ee158)/ee5 + ee164) * ee24 *
      ee2/ee399);
    out(j, 10) += (((((2 * (ee400 - ee24/ee44) + 4 * ee400)/ee44 +
      8 * (ee24/ee254))/ee37 + 4 * (ee60/ee44)) * ee24/ee37 + ee297 *
      ee60/ee44 + 2 * (R_pow(ee60, 2) - ee382) - ee382) * ee24/
        ee37 - ((((ee410 - (ee15 * (2 - ee27 * ee54) + ee114 + 3) *
          ee5)/ee17 + ee247/R_pow(ee11, (2 + ee114)) + 2 * ee60) *
          ee12 * ee5 - ee60/ee18)/ee11 + ee381/ee37)/ee14) * ee24/(ee37 *
          R_pow(ee6, 4));
    out(j, 11) += (((ee332 + 4 * (ee81/ee207)) * ee24/ee37 + ee384 +
      2 * (ee93 * ee60 - ee306) - ee306) * ee24/ee37 - ((ee250 +
      ee381 * ee7/ee209 + ee243)/ee14 + ee165 * ee60 + ((ee413 +
      2 * (ee192/ee14))/ee17 - (ee299/ee55 + ee7/(R_pow(ee11, (ee114 +
      4)) * ee6)) * ee5 * ee15) * ee12 * ee5)) * ee24/ee402;
    out(j, 12) += ((((ee314 + ee293/ee255 - (ee36/ee253 + ee113/
      ee199))/ee44 + 2 * (ee377 - ee270) + 6 * ee285 - ee240)/ee14 -
        ee258 * ee29/ee5) * ee24/ee37 - ((((ee80 * ee5 - 4.5) * ee12 +
        ee273 + ee302 - ee267)/ee17 + (ee226 - 3 * ee116) * ee12)/
          ee14 + ((ee80/ee133 - ((ee29/ee55 - ee282) * ee15 + ee292)) *
            ee12 - ee291) * ee5 + ee36 * ee60)) * ee24 * ee2/(ee139 *
            ee352);
    out(j, 13) += (((ee379/ee37 + ee384 - 2 * ee306) * ee7/ee6 +
      2 * (ee395 - ee380)) * ee24/ee37 + ee182/ee133 - (((ee14 *
      ee93 + 4 - 2 * ee105) * ee7/ee223 - ((ee329 - 2)/ee55 + (ee105 -
      2)/ee17)) * ee12 * ee5 + ee380 + 2 * ee383)) * ee24/ee401;
    out(j, 14) += ((((ee377 + ee412 - ee270) * ee24/ee37 + 2 * ee128 -
      (ee273 + 3/ee5 - (ee267 + 3 * ee12))/ee17) * ee7/ee6 -
      ((ee252 + ee415)/ee5 + ee375 + (ee12 * ee29 * ee192 + 2 *
      ee113)/ee11))/ee14 + ee309 + ((2 * (ee376 - ee264) + ee416)/
        ee14 + ee345 - ee262 * ee29/ee5) * ee24/ee37 - ((ee80 * ee12 *
          ee5/ee133 + ee405/ee14) * ee7/ee6 + ee325)) * ee24 * ee2/
            ee334;
    out(j, 15) += (((ee347/ee324 - (ee90 * ee60 + (ee270 - (ee162 +
      ee412)) * ee29 * ee2/ee16))/ee5 + ee349) * ee24/ee37 + (ee97 +
      (ee51 - 1.5)/ee5)/ee17 + ee333 - ((ee308 - ee90/ee17) *
      ee12 + (((ee270 - ee162) * ee24/ee37 + 2 * ee405) * ee29/
        ee5 - 2 * ee339) * ee2/ee16)) * ee24 * ee2/ee334;
    out(j, 16) += ((((((1 + 6 * ee249 - ee210) * ee7/ee6 + 1/ee253)/
      ee44 + ee413) * ee93 + ee295 * ee24 * ee7/(ee255 * ee6) -
        ((ee5 * (2 * ee122 - 2) * ee7/ee6 - 2)/ee11 + ((4 * ee125 -
        6)/ee14 - 2 * ee82) * ee7/ee6 + 2 * ee281)/ee44) * ee24/ee37 -
        ((ee366 - (ee182 + 2 * ee182))/ee133 + ee397 + 3 * ee383)) *
        ee7/ee6 - ee196/ee17) * ee24/ee209;
    out(j, 17) += ((((ee378 - ee264) * ee24/ee37 - (ee372 + ((ee28 +
      4.5)/ee5 - ee301)/ee11 + (ee415 - ee244)/ee5 + 2 * ee375))/
        ee14 + ee309 + ((ee93 * (ee298 + 4 * ee286) - ((ee29 * (ee156 -
          2) + 3 * ee353)/ee137 + 2 * ee198))/ee14 + ee345) *
          ee24/ee37 + ee128 - ee325) * ee7/ee6 - ((ee281 * ee24/ee37 -
          ee238) * ee29 + ee409)/ee137) * ee24 * ee2/ee183;
    out(j, 18) += (ee374 - ((ee264 * ee29 * ee2/ee16 + ee368) *
      ee24/ee37 + (ee264 * ee24/ee37 + ee411) * ee29 * ee2/ee16 -
      (ee90 * ee165 + ee138))/ee5) * ee24 * ee2/ee183;
    out(j, 19) += ee364/ee183;
    out(j, 20) += (((ee379 * ee7/ee209 + ee337 * ee93 + 2 * (ee395 -
      ee367) - 2 * ee367) * ee24/ee37 + 3 * (ee93 * ee205) -
      ((ee366 + ee175 + 2 * ee196 - 1)/ee133 + ee397)) * ee7/ee6 +
      ee182/ee17) * ee24 * ee7/ee209;
    out(j, 21) += ((((ee378 - ee268) * ee24/ee37 + 2 * (ee80 * ee205) -
      (ee372 + (4.5 * ee149 - ee244)/ee5 + (ee302 - ee301)/
        ee11))/ee14 + ee309 + ee128 + ((2 * (ee376 - ee268) + ee416)/
          ee14 + ee345) * ee24/ee37 - ee325) * ee7/ee6 - ((ee251 * ee24/
            ee37 + ee217) * ee29 + ee409)/ee137) * ee24 * ee2 * ee7/
              ee183;
    out(j, 22) += (ee374 + (ee138 - ((ee368 + ee268 * ee29 * ee2/
      ee16) * ee24/ee37 + ee90 * ee205 + (ee268 * ee24/ee37 + ee411) *
        ee29 * ee2/ee16))/ee5) * ee24 * ee2 * ee7/ee183;
    out(j, 23) += ee364 * ee7/ee183;
    out(j, 24) += -(((((1.5 * (ee386 * ee13) + 3 * (ee147 * ee50) -
      3 * (ee144 * ee29))/ee5 - ((ee232 * ee90 + (((ee200 - (ee195 +
      ((ee248 + 4 * ee215 - ee351)/ee37 + 4 * ee90) * ee24/
        ee37 + ee206) * ee29) * ee2/ee16 - ee158)/ee5 + ee164) * ee29 +
          2 * (R_pow(ee90, 2) + ee237 * ee29)) * ee24/ee37 + (ee237 *
          ee24/ee37 + ((2 * ee200 - ee216) * ee2/ee16 - 2 * ee158)/
            ee5 + 2 * ee386 + 2 * ee164) * ee29 + 3 * (ee90 * ee147))) *
              ee2/ee16 - ((4.5 * ee131 - (1.5 * ((((4.5 * ee53 + 6.75 *
              ee47) * ee7/ee52 + ee98 - ee121)/ee3 - 
              3) * ee2/ee3 + 1.5) +
              3 * (ee142 * ee40))) * ee2 * ee7/ee42 - (((((ee53 * (6 *
              ee120 - 18 * ee86) + (12 * (ee124 - ee123) + 9 * (2 * ee120 +
              9 * ee86) - 324 * ee86) * ee2/ee115 + 3 * (4.5 * (ee53 *
              ee2/ee16) - ee404) - 6 * ee404)/ee5 - 3 * ee227)/ee5 + ee414 +
              ee417 + ee418 - ee398)/ee3 + 3) * ee2/ee3 - 1.5) * ee13)/
                ee18)/ee5 + ((3 * (ee132 * ee63) - (1.5 * (((ee348 + ee350) *
                  ee2/ee16 - ee304)/ee22 + ee331) + 3 * (ee152 * ee36))) *
                  ee2/ee16 - ((((1.5 * ee131 - (3 * (ee40 * ee53) + 4.5 * ee227))/
                    ee5 + 
                      ee414 + ee417 + ee418 - ee398)/ee3 + 3) * ee2/
                        ee3 - 1.5)/ee14) * ee7/ee6) * ee24 * ee2/ee399);
    
  }

}

return out;

}

