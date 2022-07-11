// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double xieps = 0.0;

// //' Extended generalized Pareto distribution of type 1 (eGPD1) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each eGPD parameter
// //' @param X1 a design matrix for the eGPD log scale parameter
// //' @param X2 a design matrix for the eGPD shape parameter
// //' @param X3 a design matrix for the eGPD log kappa
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return egpd1d0 a scalar, the negative log-liklihood
// //' @return egpd1d12 a matrix, first then second derivatives w.r.t. eGPD1 parameters
// //' @return egpd1d34 a matrix, third then fourth derivatives w.r.t. eGPD1 parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double egpd1d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, arma::vec yvec, const arma::uvec& dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec lkappavec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();

if (dcate == 1) {
  lpsivec = lpsivec.elem(dupid);
  xivec = xivec.elem(dupid);
  lkappavec = lkappavec.elem(dupid);
}

double y, lpsi, xi, lkappa;
double ee1, ee4;
double nllh=0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
lkappa = lkappavec[j];

ee1 = 1.0 / xi;
ee4 = xi * y / exp(lpsi);

if (ee4 <= -1.0) {
    nllh = 1e20;
    break;
}

nllh += (1.0 - exp(lkappa)) * log(1.0 - 1.0 / R_pow(1.0 + ee4, ee1)) + 
          (1.0 + ee1) * log1p(ee4) + lpsi - lkappa;
    
}

return(nllh);

}

// //' @rdname egpd1d0
// [[Rcpp::export]]
arma::mat egpd1d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec, const arma::uvec dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec lkappavec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 9);

if (dcate == 1) {
  lpsivec = lpsivec.elem(dupid);
  xivec = xivec.elem(dupid);
  lkappavec = lkappavec.elem(dupid);
}

double y, lpsi, xi, lkappa;
double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee14, ee16, ee17, ee18, ee19;
double ee20, ee21, ee23, ee25, ee26, ee27, ee28, ee29;
double ee31;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
lkappa = lkappavec[j];

ee1 = exp(lpsi);
ee2 = xi * y;
ee3 = ee2/ee1;
ee4 = 1 + ee3;
ee5 = 1/xi;
ee6 = 1 + ee5;
ee7 = R_pow(ee4, ee5);
ee8 = 1 - 1/ee7;
ee9 = R_pow(ee4, ee6);
ee10 = log1p(ee3);
ee11 = exp(lkappa);
ee12 = ee4 * ee1;
ee14 = 1 - ee11;
ee16 = ee10/(xi * ee7) - y/(ee9 * ee1);
ee17 = ee10/xi;
ee18 = y * ee6;
ee19 = R_pow(ee4, ee5 + 2);
ee20 = xi * ee6;
ee21 = ee8 * ee9;
ee23 = 1/ee9;
ee25 = ee11 * log(ee8);
ee26 = ee10/(xi * ee9);
ee27 = ee2/ee12;
ee28 = ee18/ee12;
ee29 = ee18/(ee19 * ee1);
ee31 = y/ee12 - 2 * ee17;

out(j, 0) = 1 - y * (ee14/(ee8 * ee7) + ee20)/ee12;
out(j, 1) = ee28 - (ee14 * ee16/ee8 + ee17)/xi;
out(j, 2) = -(1 + ee25);
out(j, 3) = y * (ee14 * (ee23 - y * (1/(ee8 * R_pow(ee4, 2 * ee6)) +
   ee20/ee19)/ee1)/ee8 - ee20 * (ee27 - 1)/ee4)/ee1;
out(j, 4) = -(y * (((ee16/ee21 + ee26)/xi - ee29) * ee14/ee8 +
   ((1 - ee27) * ee6 - ee5)/ee4)/ee1);
out(j, 5) = y * ee11/(ee21 * ee1);
out(j, 6) =  - (((((ee16/ee8 + ee17) * ee16 + ee31/ee7)/xi +
   y * ((ee23 - ee26)/xi + ee29)/ee1) * ee14/ee8 + ee31/xi)/xi +
   y * (1/R_pow(xi, 2) + ee28)/ee12);
out(j, 7) = ee11 * ee16/(xi * ee8);
out(j, 8) = -ee25;

}

return out;

}

// //' @rdname egpd1d0
// [[Rcpp::export]]
arma::mat egpd1d34(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec, const arma::uvec dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec lkappavec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 25);

if (dcate == 1) {
  lpsivec = lpsivec.elem(dupid);
  xivec = xivec.elem(dupid);
  lkappavec = lkappavec.elem(dupid);
}

double y, lpsi, xi, lkappa;
double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee19;
double ee20, ee21, ee22, ee24, ee25, ee26, ee27, ee28, ee29;
double ee30, ee31, ee34, ee35, ee37, ee39;
double ee40, ee41, ee42, ee43, ee45, ee46, ee47, ee48, ee49;
double ee50, ee51, ee52, ee53, ee54, ee55, ee57, ee59;
double ee60, ee62, ee63, ee64, ee65, ee66, ee67, ee68, ee69;
double ee70, ee71, ee72, ee73, ee74, ee77, ee78, ee79;
double ee80, ee81, ee82, ee84, ee85, ee87, ee88, ee89;
double ee90, ee91, ee93, ee94, ee95, ee96, ee97, ee98, ee99;
double ee101, ee102, ee104, ee106, ee107, ee108;
double ee110, ee111, ee113, ee114, ee115, ee116, ee117, ee118, ee119;
double ee122, ee123, ee124, ee125, ee127, ee129;
double ee131, ee132, ee134, ee138, ee139;
double ee142, ee144, ee148, ee149;
double ee150, ee152, ee154, ee155, ee157, ee158;
double ee160, ee162, ee163, ee166, ee168;
double ee171, ee172, ee173;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
lkappa = lkappavec[j];

ee1 = exp(lpsi);
ee2 = xi * y;
ee3 = ee2/ee1;
ee4 = 1/xi;
ee5 = 1 + ee3;
ee6 = 1 + ee4;
ee7 = log1p(ee3);
ee8 = R_pow(ee5, ee6);
ee9 = ee4 + 2;
ee10 = R_pow(ee5, ee4);
ee11 = R_pow(ee5, ee9);
ee12 = ee5 * ee1;
ee13 = R_pow(xi, 2);
ee14 = y/ee12;
ee15 = ee11 * ee1;
ee16 = 1 - 1/ee10;
ee19 = ee7/(xi * ee10) - y/(ee8 * ee1);
ee20 = ee7/xi;
ee21 = y * ee6;
ee22 = ee21/ee15;
ee24 = ee14 - 2 * ee20;
ee25 = 1/ee8;
ee26 = ee4 + 3;
ee27 = xi * ee8;
ee28 = R_pow(ee5, ee26);
ee29 = ee7/ee27;
ee30 = ee13 * ee11;
ee31 = ee28 * ee1;
ee34 = ee24/ee10;
ee35 = ee7/(ee13 * ee8);
ee37 = y * ((ee25 - ee29)/xi + ee22)/ee1;
ee39 = y * ee9/ee31;
ee40 = exp(lkappa);
ee41 = ee35 - ee22;
ee42 = ee7/ee30;
ee43 = ee42 - ee39;
ee45 = (ee34 + ee7 * ee19/xi)/xi + ee37;
ee46 = 1/ee11;
ee47 = 2 * ee14;
ee48 = 2 * ee6;
ee49 = 2/xi;
ee50 = ee24/ee8;
ee51 = ee2 * ee6;
ee52 = ee6 * ee43;
ee53 = 2/ee8;
ee54 = R_pow(ee5, ee48);
ee55 = ee25 - ee51/ee15;
ee57 = 1/ee30 - ee52;
ee59 = y * ee57/ee1;
ee60 = ee49 + ee14;
ee62 = (ee47 - 6 * ee20)/xi + y * ee60/ee12;
ee63 = ee7 * ee41;
ee64 = ee2/ee12;
ee65 = ee16 * ee1;
ee66 = ee50 + ee63;
ee67 = 1 - ee40;
ee68 = 2/ee11;
ee69 = xi * ee16;
ee70 = 2 * ee22;
ee71 = xi * ee6;
ee72 = ee16 * ee8;
ee73 = 2 * ee64;
ee74 = ee13 * ee28;
ee77 = (ee45 * ee7 + 2 * (ee19 * ee24))/xi;
ee78 = ee62/ee10;
ee79 = ee6 * (ee46 + xi * ee43);
ee80 = ee19/ee16;
ee81 = ee46 + ee68;
ee82 = ee2 * ee9;
ee84 = y * (((((ee53 - ee29)/xi + ee22) * ee7 - (ee50 +  ee53))/xi - ee70)/xi - ee59)/ee1;
ee85 = ee66/xi;
ee87 = R_pow(ee5, ee4 + 4) * ee1;
ee88 = ee19/ee72;
ee89 = R_pow(ee19, 2);
ee90 = ee81 - ee82/ee31;
ee91 = ee47 + 6/xi;
ee93 = (ee77 - ee78)/xi + ee84;
ee94 = ee45 * ee16;
ee95 = ee45/ee8;
ee96 = ((2 * ee80 + ee20) * ee19 + ee34)/xi;
ee97 = ee19 * ee41;
ee98 = ee89/xi;
ee99 = ee47 + ee49;
ee101 = ee7/ee74 - y * ee26/ee87;
ee102 = ee94 - ee98;
ee104 = ee66/ee13 + ee59;
ee106 = (ee24/ee11 + ee7 * ee43)/ee13 + y * (1/ee74 - ee9 *  ee101)/ee1;
ee107 = ee16 * ee54;
ee108 = ee16 * ee41;
ee110 = R_pow(ee5, ee6 + ee48);
ee111 = ee55 * ee19;
ee113 = ee91/xi + y * ee99/ee12;
ee114 = 2 - ee73;
ee115 = 2 * ee88;
ee116 = 2 * ee39;
ee117 = 4 * ee5;
ee118 = ee35 - y * (ee79 + ee46)/ee1;
ee119 = ee71/ee11;
ee122 = ee51 * ee90/ee1 - ee25;
ee123 = y * ee40;
ee124 = -(ee40 * log(ee16));
ee125 = -(y * (ee25 - y * (1/ee107 + ee119)/ee1) * ee40/ee65);
ee127 = (ee77 + (ee96 + 2 * ee45 + ee37) * ee19/ee16 - ee78)/xi +  ee84;
ee129 = ((ee95 + ((ee115 + 2 * ee29)/xi - ee70) * ee19)/ee16 +  ee85)/xi + ee59;
ee131 = (((ee80 + ee20) * ee19 + ee34)/xi + ee37) * ee40/ee69;
ee132 = ee104 * ee7;
ee134 = ee95 + ee97;
ee138 = (ee111/ee16 + ee29)/xi - y * (ee79 + (2 * (ee19/(ee69 *  ee54)) + 2 * (ee41/ee8))/ee16 + ee46)/ee1;
ee139 = ee96 + ee37;
ee142 = ee62/ee8;
ee144 = ((24 * ee20 - 6 * ee14)/xi - y * ee91/ee12)/xi -  y * ee113/ee12;
ee148 = ee16 * ee55;
ee149 = ee118/ee8;
ee150 = ee122/ee8;
ee152 = 1 + ee2 * (ee73 - 3)/ee12;
ee154 = 2 * (1 + 2 * ee3) + ee117;
ee155 = 2/ee28;
ee157 = 6 * ee3;
ee158 = 8 * ee3;
ee160 = ee40 * ee19/ee69;
ee162 = ee2 * (3 - ee73)/ee12;
ee163 = ee64 - 1;
ee166 = y * ((ee55 * (ee25 + ee53) - 2 * (y/(ee16 * ee110 *  ee1)))/ee16 + ee71 * ee90)/ee1 - ee25;
ee168 = y * (((2 * (ee7/(xi * ee11)) - ee68)/xi - ee116)/ee13 -  ee106 * ee6)/ee1;
ee171 = y * ((ee88 + ee29)/xi - ee22) * ee40/ee65;
ee172 = ee21/ee12;
ee173 = ee123/(ee72 * ee1);

out(j, 0) = y * (ee67 * ee166/ee16 - ee71 * ee152/ee5)/ee1;
out(j, 1) = y * (ee138 * ee67/ee16 - (ee6 * (ee162 - 1) - ee163/xi)/ee5)/ee1;
out(j, 2) = ee125;
out(j, 3) = y * (2 * ((ee6 * ee114 - ee49)/R_pow(ee5, 2)) -
   ee129 * ee67/ee16)/ee1;
out(j, 4) = ee171;
out(j, 5) = ee173;
out(j, 6) = -((ee127 * ee67/ee16 - ee62/xi)/xi - y * (ee60/ee13 +
   y * (1/ee13 + 2 * ee172)/ee12)/ee12);
out(j, 7) = ee131;
out(j, 8) = ee160;
out(j, 9) = ee124;
out(j, 10) = y * (ee67 * (ee25 + y * ((2 * ee150 + y * ((2 * ((ee148 +
   y/(ee54 * ee1))/ee54) + 4 * (ee148/ee54) - 8 * (y/(R_pow(ee5, 2 +
   ee48 + ee49) * ee1)))/ee16 + 4 * (ee55/ee54))/ee65 -
   (ee55 * (ee25 - y * (2/ee107 + ee119)/ee1) + 2 * (R_pow(ee55, 2) -
   ee150)))/ee16 - ee71 * (ee81 + 4/ee11 - ee82 * (ee155 +
   4/ee28 - ee2 * ee26/ee87)/ee1))/ee1)/ee16 - ee71 * (ee2 * (7 +
   ee2 * ((ee158 - ee154)/ee5 - 6)/ee12)/ee12 -
   1)/ee5)/ee1;
out(j, 11) = y * (((ee19 * ee122/ee16 - ee29)/xi + y * ((((ee115 +
   ee29)/xi - ee22) * ee55 + (4 * (ee111/ee27) + y * (2 * ((ee19/ee27 -
   ee108)/ee54) - (4 * (ee108/ee54) + 8 * (ee19/(xi * ee110))))/ee65)/ee16 +
   ee149 + 2 * (ee55 * ee41 + ee149))/ee16 +
   ee79 + ee46 + ee68 + xi * (2 * ee52 - y * (ee6 * (ee155 +
   xi * ee9 * ee101) + ee9/ee28)/ee1))/ee1) * ee67/ee16 -
   (ee6 * (1 + ee2 * (ee2 * ((ee154 - ee158)/ee5 + 6)/ee12 -
   7)/ee12) - ee152/xi)/ee5)/ee1;
out(j, 12) = -(ee123 * ee166/ee65);
out(j, 13) = y * ((((ee139 * ee55 + 2 * (ee19 * ee118))/ee16 +
   ee85)/xi + y * ((ee46 - 2 * (ee7/ee11))/ee13 + ee116 - (((((4 * (ee108/ee8) +
   8 * (ee19/(xi * ee54))) * ee19 + 2 * (ee102/ee54))/ee16 +
   4 * (ee97/ee8))/ee69 + 2 * (ee104/ee8 + R_pow(ee41, 2)))/ee16 +
   ee6 * (ee42 + xi * ee106 - ee39)))/ee1) * ee67/ee16 -
   (ee21 * (4 - ee2 * ((ee117 - ee157)/ee5 + 6)/ee12)/ee12 -
   (2 * ee162 - (2 + 2 * ee163))/ee13)/ee5)/ee1;
out(j, 14) = -(y * ee138 * ee40/ee65);
out(j, 15) = ee125;
out(j, 16) = -(y * ((((ee93/ee8 + ((((2 * (ee102/ee8) + 2 * (ee134 * ee16) +
   8 * (ee89/ee27))/ee16 + 2 * ee134)/ee16 + ee85)/xi +
   (2 * ee85 + 2 * (ee97/ee16))/xi + 3 * y * ee57/ee1) * ee19 +
   3 * (ee45 * ee41))/ee16 + (ee132 + 2 * (ee41 * ee24) -
   ee142)/xi)/xi + ee168) * ee67/ee16 + (((6 * (1 - ee64) -
   6)/xi + 2 * (y * ee114/ee12))/ee13 + y * (ee114/ee13 + y * ((2 * ee5 -
   ee157)/ee5 + 4) * ee6/ee12)/ee12)/ee5)/ee1);
out(j, 17) = y * ee129 * ee40/ee65;
out(j, 18) = ee171;
out(j, 19) = ee173;
out(j, 20) =  - (((((ee93 * ee7 + 3 * (ee45 * ee24) - 3 * (ee62 * ee19))/xi +
   ((((2 * ee102 + 4 * ee94 + 8 * ee98)/ee16 +
   4 * ee45) * ee19/ee69 + 2 * ee93) * ee19 + ee139 * ee45 +
   2 * (ee93 * ee19 + R_pow(ee45, 2)))/ee16 - ee144/ee10)/xi +
   y * (((ee142 + (ee50 + 2 * ee66 + 6/ee8 + ee63)/xi - (ee132 +
   (2 * ee24 + 6) * ee41))/xi + 3 * ee59)/xi - ee168)/ee1) * ee67/ee16 -
   ee144/xi)/xi + y * (ee113/ee13 + y * (ee99/ee13 +
   y * (2/ee13 + 6 * ee172)/ee12)/ee12)/ee12);
out(j, 21) = ee127 * ee40/ee69;
out(j, 22) = ee131;
out(j, 23) = ee160;
out(j, 24) = ee124;

}

return out;

}

// //' Extended generalized Pareto distribution of type 2 (eGPD2) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each eGPD parameter
// //' @param X1 a design matrix for the eGPD log scale parameter
// //' @param X2 a design matrix for the eGPD shape parameter
// //' @param X3 a design matrix for the eGPD log kappa1 parameter
// //' @param X4 a design matrix for the eGPD log kappa2 parameter
// //' @param X5 a design matrix for the eGPD logit p parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return egpd3d0 a scalar, the negative log-liklihood
// //' @return egpd3d12 a matrix, first then second derivatives w.r.t. eGPD1 parameters
// //' @return egpd3d34 a matrix, third then fourth derivatives w.r.t. eGPD1 parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double egpd2d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, const arma::mat& X4, const arma::mat& X5, arma::vec yvec, const arma::uvec& dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec lkappa1vec = X3 * Rcpp::as<arma::vec>(pars[2]);
arma::vec lkappa2vec = X4 * Rcpp::as<arma::vec>(pars[3]);
arma::vec logitpvec = X5 * Rcpp::as<arma::vec>(pars[4]);
int nobs = yvec.size();

if (dcate == 1) {
  lpsivec = lpsivec.elem(dupid);
  xivec = xivec.elem(dupid);
  lkappa1vec = lkappa1vec.elem(dupid);
  lkappa2vec = lkappa2vec.elem(dupid);
  logitpvec = logitpvec.elem(dupid);
}

double y, lpsi, xi, lkappa1, lkappa2, logitp;
// double ee1, ee2, ee4, ee6, ee7, ee8, ee9;
double ee1, ee4, ee7, ee8, ee9, ee10;
double nllh=0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
lkappa1 = lkappa1vec[j];
lkappa2 = lkappa2vec[j];
logitp = logitpvec[j];

ee1 = 1/xi;
ee4 = xi * y/exp(lpsi);

if (ee4 <= -1.0) {
    nllh = 1e20;
    break;
}

ee7 = 1 - 1/R_pow(1 + ee4, ee1);
ee8 = 1 + exp(-logitp);
ee9 = exp(lkappa1);
ee10 = exp(lkappa2);

nllh += (1 + ee1) * log1p(ee4) + lpsi - log((1 - 1/ee8) * R_pow(ee7, ee10 - 1) * 
ee10 + R_pow(ee7, ee9 - 1) * ee9/ee8);

// ee1 = 1 + xi * y/exp(lpsi);
// ee2 = 1/xi;
// ee4 = R_pow(ee1, 1 + ee2);
// ee6 = 1 - 1/R_pow(ee1, ee2);
// ee7 = 1 + exp(-logitp);
// ee8 = exp(lkappa1);
// ee9 = exp(lkappa2);
// 
// nllh -= log((1 - 1/ee7) * R_pow(ee6, ee9 - 1) * ee9/ee4 + R_pow(ee6, ee8 - 1) * ee8/(ee7 * ee4)) - lpsi;
    
}

return(nllh);

}

// //' @rdname egpd2d0
// [[Rcpp::export]]
arma::mat egpd2d12(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, const arma::mat& X4, const arma::mat& X5, arma::vec yvec, const arma::uvec& dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec lkappa1vec = X3 * Rcpp::as<arma::vec>(pars[2]);
arma::vec lkappa2vec = X4 * Rcpp::as<arma::vec>(pars[3]);
arma::vec logitpvec = X5 * Rcpp::as<arma::vec>(pars[4]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 20);

if (dcate == 1) {
  lpsivec = lpsivec.elem(dupid);
  xivec = xivec.elem(dupid);
  lkappa1vec = lkappa1vec.elem(dupid);
  lkappa2vec = lkappa2vec.elem(dupid);
  logitpvec = logitpvec.elem(dupid);
}

double y, lpsi, xi, lkappa1, lkappa2, logitp;

// double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
// double ee10, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
// double ee20, ee21, ee22, ee23, ee24, ee25, ee27, ee28;
// double ee31, ee35, ee36, ee37, ee38, ee39;
// double ee40, ee41, ee47, ee48, ee49;
// double ee50, ee51, ee52, ee53, ee54, ee55, ee56, ee57, ee58, ee59;
// double ee60, ee61, ee63, ee64, ee65, ee66, ee68, ee69;
// double ee73, ee76, ee78, ee79;
// double ee80, ee81, ee83, ee85, ee88, ee89;
// double ee90, ee91, ee92, ee93, ee94, ee95, ee96, ee98;
// double ee101, ee103, ee104, ee105, ee106, ee107, ee108;
// double ee114, ee115, ee116, ee117, ee118, ee119;
// double ee121, ee122, ee126, ee127;
// double ee130, ee134, ee136, ee138, ee139;
// double ee140, ee143;

double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee20, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
double ee31, ee32, ee35, ee38, ee39;
double ee40, ee42, ee44, ee45, ee46, ee48, ee49;
double ee51, ee52, ee53, ee54, ee55, ee56, ee57, ee58, ee59;
double ee60, ee61, ee62, ee63, ee65, ee67, ee68, ee69;
double ee71, ee72, ee73, ee74, ee76, ee77, ee79;
double ee80, ee82, ee83, ee84, ee85;


for (int j=0; j < nobs; j++) {
    
y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
lkappa1 = lkappa1vec[j];
lkappa2 = lkappa2vec[j];
logitp = logitpvec[j];

// ee1 = exp(lpsi);
// ee2 = xi * y;
// ee3 = ee2/ee1;
// ee4 = 1/xi;
// ee5 = 1 + ee3;
// ee6 = R_pow(ee5, ee4);
// ee7 = 1 + ee4;
// ee8 = 1 - 1/ee6;
// ee9 = exp(lkappa1);
// ee10 = exp(lkappa2);
// ee12 = exp(-logitp);
// ee13 = 1 + ee12;
// ee14 = R_pow(ee5, ee7);
// ee15 = ee9 - 1;
// ee16 = ee10 - 1;
// ee17 = R_pow(ee8, ee15);
// ee18 = R_pow(ee8, ee16);
// ee19 = log1p(ee3);
// ee20 = ee13 * ee14;
// ee21 = 1 - 1/ee13;
// ee22 = ee4 - 1;
// ee23 = R_pow(ee5, ee22);
// ee24 = 2 * ee7;
// ee25 = ee17 * ee9;
// ee27 = y * ee23/ee1;
// ee28 = R_pow(ee20, 2);
// ee31 = ee27 - ee6 * ee19/xi;
// ee35 = ee21 * ee18 * ee10/ee14 + ee25/ee20;
// ee36 = ee9 - 2;
// ee37 = ee10 - 2;
// ee38 = R_pow(ee5, ee24);
// ee39 = R_pow(xi, 2);
// ee40 = R_pow(ee8, ee36);
// ee41 = R_pow(ee8, ee37);
// ee47 = y * ee7 * ee6/ee1 - ee14 * ee19/ee39;
// ee48 = ee40 * ee15;
// ee49 = 2/xi;
// ee50 = log(ee8);
// ee51 = ee41 * ee16;
// ee52 = 3/xi;
// ee53 = 1 + ee52;
// ee54 = R_pow(ee5, ee53);
// ee55 = R_pow(ee13, 2);
// ee56 = R_pow(ee5, ee49);
// ee57 = ee18 * ee10;
// ee58 = xi * ee54;
// ee59 = ee55 * ee14;
// ee60 = R_pow(ee5, ee4 - ee24);
// ee61 = xi * ee17;
// ee63 = xi * ee18 * ee7;
// ee64 = ee35 * ee1;
// ee65 = ee17 * ee13;
// ee66 = ee17 * ee14;
// ee68 = ee17 + ee25 * ee50;
// ee69 = ee48 * ee31;
// ee73 = ee18 + ee57 * ee50;
// ee76 = ee51 * ee31/ee58 - ee18 * ee47/ee38;
// ee78 = ee51/ee38 - ee63 * ee60;
// ee79 = ee61 * ee7;
// ee80 = xi * ee13;
// ee81 = ee35 * ee13;
// ee83 = (ee69/(ee80 * ee54) - ee65 * ee47/ee28) * ee9 + ee76 *  ee21 * ee10;
// ee85 = (ee48/(ee13 * ee38) - ee79 * ee13 * ee6/ee28) * ee9 +  ee78 * ee21 * ee10;
// ee88 = ee66 * ee9/ee28 - ee57/ee59;
// ee89 = R_pow(ee5, ee4 - 2);
// ee90 = ee5 * ee1;
// ee91 = ee81 * ee14;
// ee92 = ee35 * ee14;
// ee93 = ee7 * ee31;
// ee94 = ee23 * ee19;
// ee95 = R_pow(ee31, 2);
// ee96 = 1 + ee49;
// ee98 = xi * ee7;
// ee101 = y * ee89 * ee22/ee1;
// ee103 = y/ee90 - 2 * (ee19/xi);
// ee104 = ee17 * ee55;
// ee105 = R_pow(ee8, ee9 - 3);
// ee106 = R_pow(ee8, ee10 - 3);
// ee107 = R_pow(ee5, ee96);
// ee108 = 4 * ee7;
// ee114 = ee40 * ee13 * ee47;
// ee115 = ee48 * ee50;
// ee116 = ee105 * ee36;
// ee117 = ee51 * ee50;
// ee118 = ee41 * ee47;
// ee119 = ee106 * ee37;
// ee121 = ee93 + ee6;
// ee122 = R_pow(ee5, ee7 + ee24);
// ee126 = ee6 + ee27;
// ee127 = ee38 * ee1;
// ee130 = R_pow(ee47, 2);
// ee134 = (ee101 - ee94/ee39)/ee56 - 2 * (ee31/(xi * ee107));
// ee136 = (y * (ee101 - (ee23 + ee94/xi)/xi)/ee1 - (ee6 * ee103 + ee19 * ee31/xi)/xi)/ee56 - 2 * (ee95/(xi * R_pow(ee5, ee52)));
// ee138 = 2 * (y/(R_pow(ee5, ee4 + 2) * ee1)) - (ee23 + ee2 * ee89 * ee22/ee1)/ee56;
// ee139 = ee98 * ee13;
// ee140 = xi * R_pow(ee5, 4/xi);
// ee143 = y * (ee93 - ee6/xi)/ee1 - (ee14 * ee103 + ee19 *  ee47)/xi;
// 
// out(j, 0) = 1 + y * ee85/ee64;
// out(j, 1) = -(ee83/ee35);
// out(j, 2) = -(ee68 * ee9/ee91);
// out(j, 3) = -(ee73 * ee21 * ee10/ee92);
// out(j, 4) = -(ee88 * ee12/ee35);
// out(j, 5) = y * ((((ee40 * ee138 - y * ee105 * ee36/ee127)/(ee13 * ee6) +
//    ee2 * ee40 * ee7 * ee13/(ee28 * ee1)) * ee15/ee5 -
//    ee139 * (y * (2 * (ee79 * ee55 * ee54/ee28) - ee48/ee5)/ee1 -
//    ee126 * ee17)/ee28) * ee9 + (((ee41 * ee138 - y * ee106 * ee37/ee127)/ee14 +
//    ee2 * ee41 * ee7/(R_pow(ee5, 1 + ee24) * ee1)) * ee16 -
//    ee98 * (2 * (ee2 * ee18 * ee7 * R_pow(ee5, ee53 -
//    ee108)/ee1) - (ee126 * ee18 + y * ee41 * ee16/ee90)/ee38)) * ee21 * ee10 +
//    y * R_pow(ee85, 2)/ee64)/ee64;
// out(j, 6) = y * ((((ee134 * ee40 + ee116 * ee31/ee58)/ee20 -
//    ee114/(ee28 * ee14)) * ee15 - ((ee69/ee6 - 2 * (ee61 * ee55 * ee107 * ee47/ee28)) * ee7 +
//    ee121 * ee17) * ee13/ee28) * ee9 +
//    (((ee134 * ee41 + ee119 * ee31/ee58)/ee14 - ee118/ee122) * ee16 +
//    2 * (ee63 * R_pow(ee5, ee96 - ee108) * ee47) - (ee121 * ee18 +
//    ee41 * ee7 * ee16 * ee31/ee6)/ee38) * ee21 * ee10 -
//    ee83 * ee85/ee35)/ee64;
// out(j, 7) = y * (((ee115/ee14 + ee40/ee14) * ee9 + ee48/ee14)/ee20 -
//    (ee85/ee91 + ee139 * ee6/ee28) * ee68) * ee9/ee64;
// out(j, 8) = y * (((ee117/ee14 + ee41/ee14) * ee10 + ee51/ee14)/ee14 -
//    (ee85/ee92 + ee98 * ee60) * ee73) * ee21 * ee10/ee64;
// out(j, 9) = y * ((ee48 + xi * (ee17 * ee6 - 2 * (ee104 * R_pow(ee5, ee4 +
//    ee24)/ee28)) * ee7) * ee9/ee28 - (ee85 * ee88/ee35 +
//    ee78 * ee10/ee55)) * ee12/ee64;
// out(j, 10) =  - ((((((ee136 * ee41 + ee119 * ee95/ee140)/ee14 -
//    ee118 * ee31/R_pow(ee5, ee24 + ee49)) * ee16 - (ee18 * ee143 +
//    ee51 * ee47 * ee31/ee56)/ee38)/xi + 2 * (ee18 * R_pow(ee5, ee7 -
//    ee108) * ee130)) * ee21 * ee10 + (((ee136 * ee40 +
//    ee116 * ee95/ee140)/ee20 - ee114 * ee31/(ee28 * ee56)) * ee15/xi -
//    ((ee17 * ee143 + ee48 * ee47 * ee31/ee56)/xi - 2 * (ee104 * ee14 * ee130/ee28)) * ee13/ee28) * ee9 -
//    R_pow(ee83, 2)/ee35)/ee35);
// out(j, 11) = -((((ee115/ee56 + ee40/ee56) * ee9 + ee48/ee56) * ee31/(ee80 * ee14) -
//    (ee83/ee91 + ee13 * ee47/ee28) * ee68) * ee9/ee35);
// out(j, 12) = -((((ee117/ee56 + ee41/ee56) * ee10 + ee51/ee56) * ee31/(xi * ee14) -
//    (ee83/ee92 + ee47/ee38) * ee73) * ee21 * ee10/ee35);
// out(j, 13) =  - ((((ee17 - 2 * (ee104 * ee38/ee28)) * ee47 +
//    ee40 * R_pow(ee5, 1 - ee4) * ee15 * ee31/xi) * ee9/ee28 - (ee83 * ee88/ee35 +
//    ee76 * ee10/ee55)) * ee12/ee35);
// out(j, 14) =  - ((((ee68 + 2 * ee17) * ee9 * ee50 + ee17)/ee14 -
//    R_pow(ee68, 2) * ee9/(ee81 * ee38)) * ee9/ee81);
// out(j, 15) = ee68 * ee73 * ee21 * ee9 * ee10/(R_pow(ee35, 2) * ee13 * ee38);
// out(j, 16) = -(ee68 * (ee14/ee28 - ee88/ee91) * ee12 * ee9/ee35);
// out(j, 17) =  - ((((ee73 + 2 * ee18) * ee10 * ee50 + ee18)/ee14 -
//    R_pow(ee73, 2) * ee21 * ee10/(ee35 * ee38)) * ee21 * ee10/ee35);
// out(j, 18) = (ee88 * ee21/ee92 + 1/ee59) * ee73 * ee12 * ee10/ee35;
// out(j, 19) = ((ee66 - 2 * (ee65 * ee122 * ee12/ee28)) * ee9/ee28 +
//    R_pow(ee88, 2) * ee12/ee35 - ee18 * (1 - 2 * (ee12/ee13)) * ee10/ee59) * ee12/ee35;   
// 
ee1 = exp(lpsi);
ee2 = xi * y;
ee3 = ee2/ee1;
ee4 = 1/xi;
ee5 = 1 + ee3;
ee6 = R_pow(ee5, ee4);
ee7 = 1 - 1/ee6;
ee8 = exp(lkappa1);
ee9 = exp(lkappa2);
ee11 = exp(-logitp);
ee12 = 1 + ee11;
ee13 = ee8 - 1;
ee14 = ee9 - 1;
ee15 = R_pow(ee7, ee13);
ee16 = R_pow(ee7, ee14);
ee17 = 1 - 1/ee12;
ee18 = 1 + ee4;
ee19 = ee15 * ee8;
ee20 = R_pow(ee5, ee18);
ee23 = ee17 * ee16 * ee9 + ee19/ee12;
ee24 = ee8 - 2;
ee25 = ee9 - 2;
ee26 = R_pow(ee7, ee24);
ee27 = R_pow(ee7, ee25);
ee28 = log1p(ee3);
ee29 = log(ee7);
ee31 = ee26 * ee8 * ee13;
ee32 = ee16 * ee9;
ee35 = ee17 * ee27 * ee9 * ee14;
ee38 = ee28/(xi * ee6) - y/(ee20 * ee1);
ee39 = ee5 * ee1;
ee40 = ee19 - ee32;
ee42 = ee15 + ee19 * ee29;
ee44 = ee16 + ee32 * ee29;
ee45 = ee23 * ee12;
ee46 = ee35 + ee31/ee12;
ee48 = ee35/ee20 + ee31/(ee12 * ee20);
ee49 = R_pow(ee12, 2);
ee51 = R_pow(ee5, ee4 + 2) * ee1;
ee52 = y * ee18;
ee53 = ee23 * ee49;
ee54 = xi * ee23;
ee55 = ee26 * ee13;
ee56 = ee27 * ee14;
ee57 = 1/ee20;
ee58 = ee28/xi;
ee59 = xi * ee20;
ee60 = ee52/ee51;
ee61 = R_pow(ee7, ee8 - 3);
ee62 = R_pow(ee7, ee9 - 3);
ee63 = R_pow(ee38, 2);
ee65 = R_pow(xi, 2);
ee67 = y/ee39 - 2 * ee58;
ee68 = ee23 * ee1;
ee69 = ee40/ee45;
ee71 = (ee67/ee6 + ee28 * ee38/xi)/xi + y * ((ee57 - ee28/ee59)/xi +  ee60)/ee1;
ee72 = ee55 * ee29;
ee73 = ee61 * ee24;
ee74 = ee56 * ee29;
ee76 = ee27 * ee9 * ee14;
ee77 = ee62 * ee25;
ee79 = R_pow(ee5, 2 * ee18) * ee1;
ee80 = ee57 - ee2 * ee18/ee51;
ee82 = ee28/(ee65 * ee20) - ee60;
ee83 = xi * ee18;
ee84 = ee2/ee39;
ee85 = ee52/ee39;

out(j, 0) = 1 + y * (ee48/ee23 - ee83/ee5)/ee1;
out(j, 1) = (ee46 * ee38/ee23 - ee58)/xi + ee85;
out(j, 2) = -(ee42 * ee8/ee45);
out(j, 3) = -(ee44 * ee17 * ee9/ee23);
out(j, 4) = -(ee40 * ee11/ee53);
out(j, 5) =  - (y * (((ee26 * ee80 + y * ee61 * ee24/ee79) * ee8 * ee13/ee12 +
   (ee27 * ee80 + y * ee62 * ee25/ee79) * ee17 * ee9 * ee14 -
   y * R_pow(ee48, 2)/ee68)/ee23 + ee83 * (ee84 -
   1)/ee5)/ee1);
out(j, 6) = -(y * (((ee73 * ee38/ee59 - ee26 * ee82) * ee8 * ee13/ee12 +
   (ee77 * ee38/ee59 - ee27 * ee82) * ee17 * ee9 * ee14 -
   ee46 * ee48 * ee38/ee54)/ee23 + ((1 - ee84) * ee18 -
   ee4)/ee5)/ee1);
out(j, 7) = -(y * (ee48 * ee42/ee23 - ((ee72/ee20 + ee26/ee20) * ee8 +
   ee55/ee20)) * ee8/(ee45 * ee1));
out(j, 8) = -(y * (ee48 * ee44/ee23 - ((ee74/ee20 + ee27/ee20) * ee9 +
   ee56/ee20)) * ee17 * ee9/ee68);
out(j, 9) = -(y * (ee48 * ee40/ee23 + ee76/ee20 - ee31/ee20) * ee11/(ee53 * ee1));
out(j, 10) = (((ee71 * ee26 - ee73 * ee63/xi) * ee8 * ee13/ee12 +
   (ee71 * ee27 - ee77 * ee63/xi) * ee17 * ee9 * ee14 + R_pow(ee46, 2) * ee63/ee54)/ee23 -
   ee67/xi)/xi - y * (1/ee65 +
   ee85)/ee39;
out(j, 11) = ((ee26 + ee72) * ee8 + ee55 - ee46 * ee42/ee23) * ee8 * ee38/(ee54 * ee12);
out(j, 12) = ((ee27 + ee74) * ee9 + ee56 - ee46 * ee44/ee23) * ee17 * ee9 * ee38/ee54;
out(j, 13) = (ee31 - (ee46 * ee40/ee23 + ee76)) * ee11 * ee38/(ee54 * ee49);
out(j, 14) =  - ((((ee42 + 2 * ee15) * ee29 - R_pow(ee42, 2)/ee45) * ee8 +
   ee15) * ee8/ee45);
out(j, 15) = ee42 * ee44 * ee17 * ee8 * ee9/(R_pow(ee23, 2) * ee12);
out(j, 16) = -(ee42 * (1 - ee69) * ee11 * ee8/ee53);
out(j, 17) =  - ((((ee44 + 2 * ee16) * ee29 - R_pow(ee44, 2) * ee17/ee23) * ee9 +
   ee16) * ee17 * ee9/ee23);
out(j, 18) = (ee40 * ee17/ee23 + 1) * ee44 * ee11 * ee9/ee53;
out(j, 19) = ((ee69 - 2) * ee11/ee12 + 1) * ee40 * ee11/ee53;
    

}

return out;

}

// //' @rdname egpd2d0
// [[Rcpp::export]]
arma::mat egpd2d34(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, const arma::mat& X4, const arma::mat& X5, arma::vec yvec, const arma::uvec& dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec lkappa1vec = X3 * Rcpp::as<arma::vec>(pars[2]);
arma::vec lkappa2vec = X4 * Rcpp::as<arma::vec>(pars[3]);
arma::vec logitpvec = X5 * Rcpp::as<arma::vec>(pars[4]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 105);

if (dcate == 1) {
  lpsivec = lpsivec.elem(dupid);
  xivec = xivec.elem(dupid);
  lkappa1vec = lkappa1vec.elem(dupid);
  lkappa2vec = lkappa2vec.elem(dupid);
  logitpvec = logitpvec.elem(dupid);
}

double y, lpsi, xi, lkappa1, lkappa2, logitp;

double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28;
double ee31, ee32, ee33, ee35, ee36, ee37, ee38, ee39;
double ee40, ee41, ee42, ee43, ee44, ee45, ee46, ee47, ee48, ee49;
double ee50, ee51, ee54, ee56, ee57, ee58, ee59;
double ee60, ee62, ee63, ee64, ee65, ee66, ee67, ee68, ee69;
double ee70, ee71, ee72, ee73, ee74, ee75, ee76, ee77, ee78, ee79;
double ee80, ee81, ee82, ee83, ee84, ee85, ee86, ee87, ee88;
double ee90, ee91, ee92, ee93, ee94, ee95, ee96, ee97, ee98, ee99;
double ee100, ee101, ee102, ee103, ee104, ee105, ee106, ee107, ee109;
double ee111, ee115, ee116, ee118;
double ee120, ee121, ee123, ee127;
double ee131, ee132, ee133, ee134, ee135, ee136, ee137, ee138, ee139;
double ee140, ee142, ee143, ee144, ee145, ee146, ee149;
double ee151, ee152, ee153, ee154, ee155, ee157, ee158, ee159;
double ee161, ee162, ee163, ee164, ee165, ee166, ee169;
double ee170, ee171, ee172, ee173, ee176, ee177, ee178, ee179;
double ee180, ee182, ee183, ee184, ee187;
double ee190, ee191, ee192, ee193, ee194, ee195, ee197, ee199;
double ee200, ee201, ee202, ee203, ee205, ee207;
double ee210, ee211, ee212, ee214, ee216, ee217, ee218, ee219;
double ee220, ee221, ee222, ee223, ee224, ee225, ee226, ee227, ee228, ee229;
double ee230, ee231, ee232, ee233, ee234, ee238, ee239;
double ee241, ee243, ee244, ee245, ee246, ee247, ee248, ee249;
double ee250, ee251, ee252, ee254, ee255, ee256, ee257, ee258, ee259;
double ee260, ee263, ee264, ee267, ee268, ee269;
double ee270, ee271, ee272, ee274, ee275, ee276, ee277, ee279;
double ee280, ee281, ee282, ee286, ee289;
double ee293, ee294, ee295, ee297, ee298;
double ee300, ee301, ee303, ee304, ee305, ee306, ee307, ee308, ee309;
double ee310, ee312, ee314, ee317;
double ee320, ee321, ee322, ee323, ee324, ee326, ee327, ee328;
double ee332, ee337, ee339;
double ee340, ee342, ee344, ee345, ee346, ee347, ee348, ee349;
double ee351, ee353, ee355, ee357, ee359;
double ee361, ee362, ee364, ee365, ee367, ee368;
double ee371, ee373, ee374, ee376, ee377;
double ee380, ee383, ee384, ee385, ee386, ee387, ee388;
double ee390, ee391, ee393, ee395, ee399;
double ee400, ee404, ee405, ee406, ee407, ee409;
double ee411, ee412, ee413;
double ee420, ee421, ee423;
double ee430, ee433, ee435, ee436, ee437, ee438, ee439;
double ee440, ee441, ee442, ee445, ee446, ee448;
double ee450, ee451, ee452, ee453, ee455, ee456, ee457, ee458, ee459;
double ee461, ee463, ee465, ee467, ee468, ee469;
double ee470, ee471, ee472, ee473, ee474, ee475, ee477, ee479;
double ee480, ee482, ee483, ee485, ee486, ee487;
double ee490, ee492, ee493, ee494, ee497, ee499;
double ee500, ee501, ee502, ee503, ee505, ee506, ee508;
double ee510, ee512, ee514, ee515, ee516, ee518, ee519;
double ee520, ee521, ee522, ee523, ee524, ee525, ee526, ee527, ee529;
double ee530, ee531, ee533, ee535, ee537, ee538, ee539;
double ee540, ee541, ee545, ee547, ee548;
double ee550, ee553, ee554, ee555, ee556, ee557, ee559;
double ee560, ee561, ee562, ee563, ee564, ee565, ee566, ee567, ee569;
double ee570, ee571, ee572, ee574, ee575, ee576, ee577, ee578, ee579;
double ee581, ee583, ee584, ee585, ee586, ee587, ee589;
double ee590, ee591, ee592, ee594, ee595, ee596, ee597, ee598, ee599;
double ee601, ee602, ee603, ee604, ee606, ee607, ee608, ee609;
double ee611, ee613, ee614, ee615, ee616, ee617, ee618, ee619;
double ee620, ee621, ee622, ee623, ee625, ee626, ee627, ee628, ee629;
double ee630, ee632, ee633, ee634, ee635, ee637, ee638;
double ee640, ee642, ee643, ee644, ee645, ee647, ee648;
double ee650, ee651, ee654, ee655, ee658;
double ee660, ee661, ee662, ee663, ee664, ee665, ee668;
double ee670, ee671, ee672, ee673, ee674, ee675, ee677, ee678;
double ee681, ee682, ee685, ee686, ee687, ee688, ee689;
double ee690, ee691, ee692, ee693, ee694, ee695, ee696, ee697, ee698, ee699;
double ee700, ee701, ee702, ee703, ee704, ee705, ee706, ee707, ee708, ee709;
double ee710, ee711, ee712, ee713, ee717, ee718;
double ee722, ee723, ee724, ee725, ee726, ee727, ee728, ee729;
double ee730, ee731, ee732, ee733, ee734, ee735, ee736, ee737, ee738, ee739;
double ee740, ee741, ee742, ee744, ee746, ee748, ee749;
double ee750, ee752, ee753, ee754, ee755, ee756, ee757, ee758, ee759;
double ee760, ee761, ee762, ee763, ee764, ee765, ee766, ee767, ee768, ee769;
double ee770, ee771, ee772, ee773, ee775, ee776, ee777;
double ee780, ee781;

for (int j=0; j < nobs; j++) {
    
y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
lkappa1 = lkappa1vec[j];
lkappa2 = lkappa2vec[j];
logitp = logitpvec[j];

ee1 = exp(lpsi);
ee2 = xi * y;
ee3 = ee2/ee1;
ee4 = 1/xi;
ee5 = 1 + ee3;
ee6 = R_pow(ee5, ee4);
ee7 = 1 - 1/ee6;
ee8 = exp(lkappa2);
ee9 = exp(lkappa1);
ee10 = 1 + ee4;
ee12 = exp(-logitp);
ee13 = R_pow(ee5, ee10);
ee14 = 1 + ee12;
ee15 = ee8 - 1;
ee16 = ee9 - 1;
ee17 = log1p(ee3);
ee18 = R_pow(ee7, ee15);
ee19 = ee8 - 2;
ee20 = ee9 - 2;
ee21 = R_pow(ee7, ee16);
ee22 = 1 - 1/ee14;
ee23 = ee4 + 2;
ee24 = R_pow(ee7, ee19);
ee25 = R_pow(ee7, ee20);
ee26 = ee21 * ee9;
ee27 = R_pow(ee5, ee23);
ee28 = log(ee7);
ee31 = ee17/(xi * ee6) - y/(ee13 * ee1);
ee32 = ee27 * ee1;
ee33 = ee26/ee14;
ee35 = ee22 * ee18 * ee8;
ee36 = ee35 + ee33;
ee37 = y * ee10;
ee38 = ee37/ee32;
ee39 = R_pow(xi, 2);
ee40 = ee5 * ee1;
ee41 = 1/ee13;
ee42 = ee25 * ee9;
ee43 = ee9 - 3;
ee44 = ee8 - 3;
ee45 = ee18 * ee8;
ee46 = ee42 * ee16;
ee47 = y/ee40;
ee48 = R_pow(ee7, ee43);
ee49 = R_pow(ee7, ee44);
ee50 = xi * ee13;
ee51 = ee17/xi;
ee54 = ee22 * ee24 * ee8 * ee15;
ee56 = ee47 - 2 * ee51;
ee57 = ee17/(ee39 * ee13);
ee58 = ee17/ee50;
ee59 = ee57 - ee38;
ee60 = ee2 * ee10;
ee62 = (ee56/ee6 + ee17 * ee31/xi)/xi + y * ((ee41 - ee58)/xi +  ee38)/ee1;
ee63 = R_pow(ee31, 2);
ee64 = ee24 * ee15;
ee65 = ee26 - ee45;
ee66 = ee41 - ee60/ee32;
ee67 = ee45 * ee28;
ee68 = ee25 * ee16;
ee69 = ee18 + ee67;
ee70 = ee26 * ee28;
ee71 = ee48 * ee20;
ee72 = ee49 * ee19;
ee73 = ee21 + ee70;
ee74 = ee46/ee14;
ee75 = ee4 + 3;
ee76 = ee54/ee13;
ee77 = ee46/(ee14 * ee13);
ee78 = ee54 + ee74;
ee79 = ee39 * ee27;
ee80 = R_pow(ee5, 2 * ee10);
ee81 = ee76 + ee77;
ee82 = R_pow(ee5, ee75);
ee83 = ee80 * ee1;
ee84 = ee82 * ee1;
ee85 = ee24 * ee8;
ee86 = R_pow(ee14, 2);
ee87 = ee64 * ee28;
ee88 = ee68 * ee28;
ee90 = y * ee23/ee84;
ee91 = 2 * ee18;
ee92 = ee17/ee79;
ee93 = 2 * ee21;
ee94 = ee92 - ee90;
ee95 = ee24/ee13;
ee96 = ee36 * ee14;
ee97 = ee85 * ee15;
ee98 = ee25/ee13;
ee99 = 1/ee27;
ee100 = ee9 - 4;
ee101 = ee8 - 4;
ee102 = R_pow(ee7, ee100);
ee103 = R_pow(ee7, ee101);
ee104 = ee24 + ee87;
ee105 = ee25 + ee88;
ee106 = ee56/ee13;
ee107 = ee10 * ee94;
ee109 = ee62 * ee25 - ee71 * ee63/xi;
ee111 = ee62 * ee24 - ee72 * ee63/xi;
ee115 = y * (1/ee79 - ee107)/ee1;
ee116 = 2/ee13;
ee118 = ee25 * ee66 + y * ee48 * ee20/ee83;
ee120 = ee24 * ee66 + y * ee49 * ee19/ee83;
ee121 = ee12/ee14;
ee123 = ee109 * ee9 * ee16;
ee127 = ee71 * ee31/ee50 - ee25 * ee59;
ee131 = ee72 * ee31/ee50 - ee24 * ee59;
ee132 = R_pow(ee36, 2);
ee133 = ee64/ee13;
ee134 = 1 - 2 * ee121;
ee135 = ee65/ee96;
ee136 = ee68/ee13;
ee137 = 2 * ee47;
ee138 = ee46/ee13;
ee139 = ee69 + ee91;
ee140 = ee97/ee13;
ee142 = ee118 * ee9 * ee16;
ee143 = ee102 * ee43;
ee144 = ee103 * ee44;
ee145 = xi * ee36;
ee146 = ee123/ee14;
ee149 = ee111 * ee22 * ee8 * ee15;
ee151 = ee127 * ee9 * ee16;
ee152 = ee73 + ee93;
ee153 = ee146 + ee149;
ee154 = ee104 * ee8;
ee155 = 2/xi;
ee157 = (ee87/ee13 + ee95) * ee8 + ee133;
ee158 = ee154 + ee64;
ee159 = ee105 * ee9;
ee161 = (ee88/ee13 + ee98) * ee9 + ee136;
ee162 = ee159 + ee68;
ee163 = R_pow(ee69, 2);
ee164 = ee17 * ee59;
ee165 = ee142/ee14;
ee166 = ee163 * ee22;
ee169 = ee120 * ee22 * ee8 * ee15;
ee170 = ee71 * ee28;
ee171 = ee72 * ee28;
ee172 = ee165 + ee169;
ee173 = ee151/ee14;
ee176 = ee131 * ee22 * ee8 * ee15;
ee177 = ee106 + ee164;
ee178 = 2/ee27;
ee179 = ee155 + ee47;
ee180 = ee36 * ee86;
ee182 = (ee137 - 6 * ee51)/xi + y * ee179/ee40;
ee183 = ee173 + ee176;
ee184 = ee177/ee39;
ee187 = ee152 * ee9 * ee28 + ee21;
ee190 = ee139 * ee8 * ee28 + ee18;
ee191 = R_pow(ee73, 2);
ee192 = ee184 + ee115;
ee193 = ee48 * ee66;
ee194 = ee49 * ee66;
ee195 = ee65 * ee22;
ee197 = ee65 * ee12/ee86;
ee199 = ee10 * (ee99 + xi * ee94) + ee99;
ee200 = ee99 + ee178;
ee201 = ee2 * ee23;
ee202 = ee46 - ee97;
ee203 = ee140 - ee138;
ee205 = ((ee62 * ee17 + 2 * (ee31 * ee56))/xi - ee182/ee6)/xi +  y * (((((ee116 - ee58)/xi + ee38) * ee17 - (ee106 + ee116))/xi -  2 * ee38)/xi - ee115)/ee1;
ee207 = ee57 - y * ee199/ee1;
ee210 = ee60 * (ee200 - ee201/ee84)/ee1 - ee41;
ee211 = ee48 + ee170;
ee212 = ee49 + ee171;
ee214 = ee143 * ee63/xi;
ee216 = ee144 * ee63/xi;
ee217 = 2 * ee135;
ee218 = 2 * ee25;
ee219 = 2 * ee24;
ee220 = ee42 * ee28;
ee221 = ee48 * ee59;
ee222 = ee85 * ee28;
ee223 = ee49 * ee59;
ee224 = ee195/ee36;
ee225 = ee166 * ee8;
ee226 = ee25 + ee220;
ee227 = ee24 + ee222;
ee228 = ee180 * ee1;
ee229 = ee36 * ee1;
ee230 = ee211 * ee16;
ee231 = ee212 * ee15;
ee232 = ee36 * ee134;
ee233 = ee48/ee13;
ee234 = ee49/ee13;
ee238 = (ee217 - 2) * ee12/ee14 + 1;
ee239 = 2 * (ee197 - ee232);
ee241 = ee62 * ee48 - ee214;
ee243 = ee62 * ee49 - ee216;
ee244 = ee78 * ee65;
ee245 = 1 + 2 * ee12;
ee246 = 8 * ee197;
ee247 = ee39 * ee82;
ee248 = (ee135 - 2) * ee12;
ee249 = ee132 * ee14;
ee250 = ee78 * ee69;
ee251 = ee191 * ee9;
ee252 = ee139 * ee28;
ee254 = 2 * (ee166/ee36);
ee255 = 2 * ee245;
ee256 = ee145 * ee86;
ee257 = ee190 * ee36;
ee258 = ee78 * ee81;
ee259 = ee251/ee14;
ee260 = ee193 + y * ee102 * ee43/ee83;
ee263 = ee143 * ee31/ee50 - ee221;
ee264 = ee194 + y * ee103 * ee44/ee83;
ee267 = ee144 * ee31/ee50 - ee223;
ee268 = 2 * ee14;
ee269 = ee81 * ee65;
ee270 = ee230 + ee71;
ee271 = ee231 + ee72;
ee272 = ee226 + ee218;
ee274 = ee170/ee13 + ee233;
ee275 = ee48 * ee31;
ee276 = ee48 + 2 * ee48;
ee277 = ee227 + ee219;
ee279 = ee171/ee13 + ee234;
ee280 = ee49 * ee31;
ee281 = ee49 + 2 * ee49;
ee282 = ee255 + ee268;
ee286 = (ee272 * ee16 + 2 * ee42) * ee28 + ee25 + ee25 +  ee25;
ee289 = (ee252 - ee254) * ee8 + ee18;
ee293 = (ee277 * ee15 + 2 * ee85) * ee28 + ee24 + ee24 +  ee24;
ee294 = ee78 * ee73;
ee295 = ee81 * ee69;
ee297 = (ee26 - 2 * ee45) * ee22 - ee33;
ee298 = ee152 * ee28;
ee300 = ee274 * ee16;
ee301 = ee279 * ee15;
ee303 = R_pow(ee5, ee4 + 4) * ee1;
ee304 = 2 * (ee257 + ee225);
ee305 = 2 * ee224;
ee306 = 2 * (ee191/ee96);
ee307 = 8 * ee225;
ee308 = 8 * ee12;
ee309 = ee12/ee86;
ee310 = ee2/ee40;
ee312 = ee205 * ee25 - (ee62 * ee276 - ee214) * ee20 * ee31/xi;
ee314 = ee205 * ee24 - (ee62 * ee281 - ee216) * ee19 * ee31/xi;
ee317 = (ee241/ee13 + 2 * (ee275 * ee59)) * ee20/xi - ee192 *  ee25;
ee320 = (ee243/ee13 + 2 * (ee280 * ee59)) * ee19/xi - ee192 *  ee24;
ee321 = ee109 * ee16;
ee322 = ee111 * ee15;
ee323 = ee248/ee14;
ee324 = ee96 * ee1;
ee326 = R_pow(ee78, 2) * ee63;
ee327 = ee224 + 1;
ee328 = ee118 * ee16;
ee332 = ee120 * ee15;
ee337 = ee25 * ee207 + ee20 * (y * (ee221/ee13 - ee263/ee13)/ee1 -  ee193 * ee31/xi);
ee339 = ee25 * ee210 - y * (ee260/ee13 + 2 * (ee193/ee13)) *  ee20/ee1;
ee340 = ee71/ee13;
ee342 = ee24 * ee207 + ee19 * (y * (ee223/ee13 - ee267/ee13)/ee1 -  ee194 * ee31/xi);
ee344 = ee24 * ee210 - y * (ee264/ee13 + 2 * (ee194/ee13)) *  ee19/ee1;
ee345 = ee72/ee13;
ee346 = 2 * ee297;
ee347 = 2 * ee46;
ee348 = 8 * ee195;
ee349 = y * ee81;
ee351 = ee312 * ee9 * ee16;
ee353 = ee317 * ee9 * ee16;
ee355 = (ee62 * ee105 - ee270 * ee63/xi) * ee9 + ee321;
ee357 = (ee62 * ee104 - ee271 * ee63/xi) * ee8 + ee322;
ee359 = ee111 * ee8 * ee15;
ee361 = (ee298 - ee306) * ee9 + ee21;
ee362 = ee161 * ee73;
ee364 = (ee105 * ee66 + y * (ee300/ee13 + ee71/ee80)/ee1) *  ee9 + ee328;
ee365 = ee162 * ee73;
ee367 = (ee104 * ee66 + y * (ee301/ee13 + ee72/ee80)/ee1) *  ee8 + ee332;
ee368 = ee158 * ee36;
ee371 = ee337 * ee9 * ee16;
ee373 = ee339 * ee9 * ee16;
ee374 = ee127 * ee16;
ee376 = ee120 * ee8 * ee15;
ee377 = ee131 * ee15;
ee380 = ee70 + ee93 + ee93 + ee93;
ee383 = ee67 + ee91 + ee91 + ee91;
ee384 = ee239 - ee246;
ee385 = 2 * ee138;
ee386 = 2 * ee97;
ee387 = 2 + ee305;
ee388 = 8 * (ee65/ee14);
ee390 = ee17/ee247 - y * ee75/ee303;
ee391 = xi * ee132;
ee393 = ee351/ee14 + ee314 * ee22 * ee8 * ee15;
ee395 = ee353/ee14 + ee320 * ee22 * ee8 * ee15;
ee399 = ((ee226/ee13 + ee98 + ee98) * ee16 + 2 * (ee42/ee13)) *  ee28 + ee98 + ee98 + ee98;
ee400 = ((ee230/ee13 + ee340) * ee31/xi - ee105 * ee59) *  ee9;
ee404 = ((ee227/ee13 + ee95 + ee95) * ee15 + 2 * (ee85/ee13)) *  ee28 + ee95 + ee95 + ee95;
ee405 = ((ee231/ee13 + ee345) * ee31/xi - ee104 * ee59) *  ee8;
ee406 = ee187 * ee36;
ee407 = ee157 * ee36;
ee409 = ee157 * ee69 * ee22;
ee411 = ee158 * ee69 * ee22;
ee412 = ee36 * ee202;
ee413 = ee81 * ee73;
ee420 = ee380 * ee9 * ee28 + ee21 + ee21 + ee21 + ee93 +  ee93;
ee421 = ee371/ee14;
ee423 = ee373/ee14 + ee344 * ee22 * ee8 * ee15;
ee430 = ee383 * ee8 * ee28 + ee18 + ee18 + ee18 + ee91 +  ee91;
ee433 = ee342 * ee22 * ee8 * ee15;
ee435 = ee131 * ee8 * ee15;
ee436 = ee54 + (ee347 - ee97)/ee14;
ee437 = ee387/ee14;
ee438 = 1 - ee135;
ee439 = 2 * ee140;
ee440 = 2 * ee310;
ee441 = ee137 + 6/xi;
ee442 = ee153 * ee81;
ee445 = ee362 * ee9/ee14;
ee446 = ee162 * ee81;
ee448 = ee365 * ee9/ee14;
ee450 = ee183 * ee78 * ee31;
ee451 = ee158 * ee81;
ee452 = ee36 * ee203;
ee453 = ee258 * ee31;
ee455 = ee250/ee36;
ee456 = ee78 * ee134;
ee457 = ee81 * ee202;
ee458 = ee81 * ee134;
ee459 = ee421 + ee433;
ee461 = (ee140 - ee385)/ee14 - ee76;
ee463 = (ee56/ee27 + ee17 * ee94)/ee39 + y * (1/ee247 -  ee23 * ee390)/ee1;
ee465 = ee22 * ee134;
ee467 = ee22 * (ee439 - ee138) + ee77;
ee468 = ee48 * ee9;
ee469 = ee49 * ee8;
ee470 = 2 * (ee406 + ee259);
ee471 = 2 * (ee244/ee36);
ee472 = 2 * (ee35 + (2 * ee26 - ee45)/ee14);
ee473 = 2 * ee134;
ee474 = 2 * ee90;
ee475 = 8 * ee259;
ee477 = ee399 * ee9 + ee136;
ee479 = ee404 * ee8 + ee133;
ee480 = ((ee65 * (ee238 + ee473)/ee36 + ee308 - ee282)/ee14 -  2) * ee12;
ee482 = ee286 * ee9 + ee68;
ee483 = ee400 + ee374;
ee485 = ee293 * ee8 + ee64;
ee486 = ee405 + ee377;
ee487 = ee323 + 1;
ee490 = ee420 * ee9 * ee28 + ee21;
ee492 = ee161 * ee78 + ee446;
ee493 = ee162 * ee36;
ee494 = ee450 - ee442;
ee497 = ee430 * ee8 * ee28 + ee18;
ee499 = ee157 * ee78 + ee451;
ee500 = ee409 * ee8;
ee501 = ee411 * ee8;
ee502 = ee249 * ee1;
ee503 = ee132 * ee86;
ee505 = ee78 * ee203 - ee457;
ee506 = ee326/xi;
ee508 = ee295/ee36;
ee510 = ee202 * ee12/ee86;
ee512 = (ee46 - ee386) * ee22 - ee74;
ee514 = ee203 * ee12/ee86;
ee515 = R_pow(ee7, ee9 - 5);
ee516 = R_pow(ee7, ee8 - 5);
ee518 = (ee282 - ee308)/ee14 + 2;
ee519 = 1 - ee437;
ee520 = ee304 - ee307;
ee521 = 2 * ee455;
ee522 = 2 * (ee269/ee36);
ee523 = ee346 - ee348;
ee524 = ee137 + ee155;
ee525 = ee145 * ee14;
ee526 = ee153 * ee36;
ee527 = ee153 * ee65;
ee529 = ee480/ee14 + 1;
ee530 = ee123 - ee359;
ee531 = ee505 * ee36;
ee533 = ee323 + (ee238 * ee22 + ee309) * ee65/ee36 + 1;
ee535 = ee187 * ee78 + ee448;
ee537 = ee187 * ee81 + ee445;
ee538 = ee172 * ee65;
ee539 = ee161 * ee36;
ee540 = ee183 * ee65;
ee541 = ee192 * ee17;
ee545 = ee458 + ee514;
ee547 = ee142 - ee376;
ee548 = ee510 - ee456;
ee550 = ee518 * ee12/ee14;
ee553 = ee182/ee13;
ee554 = ee143 * ee28;
ee555 = ee515 * ee100;
ee556 = ee144 * ee28;
ee557 = ee516 * ee101;
ee559 = ee441/xi + y * ee524/ee40;
ee560 = ee41 + ee116;
ee561 = 2 * (ee294/ee36);
ee562 = 2 * (ee326/ee145);
ee563 = 2 * ee508;
ee564 = 2/ee82;
ee565 = 4 * ee36;
ee566 = ee246 - ee239;
ee567 = ee391 * ee14;
ee569 = y * (((2 * (ee17/(xi * ee27)) - ee178)/xi - ee474)/ee39 -  ee463 * ee10)/ee1;
ee570 = y * R_pow(ee81, 2);
ee571 = ee526 - ee506;
ee572 = ee153 * ee69;
ee574 = (ee248 + ee65 * ((3 - ee217) * ee12/ee14 - 1)/ee36)/ee14 +  1;
ee575 = ee499 * ee36;
ee576 = ee172 * ee69;
ee577 = ee493 + ee294;
ee578 = ee183 * ee36;
ee579 = ee183 * ee69;
ee581 = ee190 * ee78 + ee501;
ee583 = ee190 * ee81 + ee500;
ee584 = ee461 * ee36;
ee585 = ee368 + ee250;
ee586 = ee36 * ee436;
ee587 = ee412 + ee244;
ee589 = ee453/xi;
ee590 = ee78 * ee520;
ee591 = ee78 * ee384;
ee592 = ee78 * ee523;
ee594 = ee151 - ee435;
ee595 = (ee48 + ee468 * ee28) * ee20;
ee596 = (ee102 + ee554) * ee20;
ee597 = ee158 - ee521;
ee598 = (ee49 + ee469 * ee28) * ee19;
ee599 = (ee103 + ee556) * ee19;
ee601 = ((24 * ee51 - 6 * ee47)/xi - y * ee441/ee40)/xi -  y * ee559/ee40;
ee602 = ee465 + ee309;
ee603 = ee22 * ee384;
ee604 = ee102 * ee66;
ee606 = ee555 * ee63/xi;
ee607 = ee97 + ee471;
ee608 = ee140 + ee522;
ee609 = ee103 * ee66;
ee611 = ee557 * ee63/xi;
ee613 = (ee254 - ee252) * ee8 - ee18;
ee614 = 1 - ee550;
ee615 = 1 - 2/ee14;
ee616 = 2 - ee440;
ee617 = ee470 - ee475;
ee618 = 2 * (ee453/ee145);
ee619 = 2 * (ee413/ee36);
ee620 = ee563 - ee157;
ee621 = ee472 - ee388;
ee622 = 2 * ee98;
ee623 = 2 * ee95;
ee625 = 3 * ee121 - 1;
ee626 = 4 * ee407;
ee627 = 4 * ee5;
ee628 = 4 * ee121;
ee629 = ee348 - ee346;
ee630 = ee307 - ee304;
ee632 = ((ee205 * ee17 + 3 * (ee62 * ee56) - 3 * (ee182 *  ee31))/xi - ee601/ee6)/xi + y * (((ee553 + (ee106 + 2 *  ee177 + 6/ee13 + ee164)/xi - (ee541 + (2 * ee56 + 6) *  ee59))/xi + 3 * ee115)/xi - ee569)/ee1;
ee633 = ee205 * ee48;
ee634 = ee205 * ee49;
ee635 = ee483 * ee73;
ee637 = ee486 * ee69 * ee22;
ee638 = ee355 * ee73;
ee640 = ee357 * ee69 * ee22;
ee642 = ee153 * ee73;
ee643 = ee492 * ee36;
ee644 = ee364 * ee73;
ee645 = ((ee595 + 2 * ee468) * ee28 + ee48) * ee16;
ee647 = ee367 * ee69 * ee22;
ee648 = ((ee598 + 2 * ee469) * ee28 + ee49) * ee15;
ee650 = (ee541 + 2 * (ee59 * ee56) - ee553)/ee39 + ee569;
ee651 = (ee62 * ee211 - (ee596 + ee143) * ee63/xi) * ee16;
ee654 = (ee62 * (ee102 + 2 * ee102) - ee606) * ee43 * ee31/xi;
ee655 = (ee62 * ee212 - (ee599 + ee144) * ee63/xi) * ee15;
ee658 = (ee62 * (ee103 + 2 * ee103) - ee611) * ee44 * ee31/xi;
ee660 = ee153 + 2 * ee153 + ee562;
ee661 = ee153 + ee562;
ee662 = ee241 * ee20;
ee663 = ee243 * ee19;
ee664 = ee420 * ee28;
ee665 = ee172 * ee36;
ee668 = ee172 * ee78 * ee31/xi;
ee670 = ee172 * ee73;
ee671 = ee362/ee96;
ee672 = ee161 * ee69;
ee673 = ee512 * ee36;
ee674 = ee365/ee96;
ee675 = ee162 * ee69;
ee677 = ee183 * ee73;
ee678 = ee430 * ee28;
ee681 = ee407 * ee69 * ee22 * ee8;
ee682 = ee409/ee36;
ee685 = ee368 * ee69 * ee22 * ee8;
ee686 = ee411/ee36;
ee687 = ee192 * ee48;
ee688 = ee192 * ee49;
ee689 = ee36 * ee467;
ee690 = ee36 * ee615;
ee691 = ee78 * ee617;
ee692 = ee78 * ee621;
ee693 = ee81 * ee629;
ee694 = ee81 * ee566;
ee695 = ee81 * ee630;
ee696 = ee135 - 1;
ee697 = ee162 - ee561;
ee698 = ee226 * ee560;
ee699 = ee300 + ee340;
ee700 = ee227 * ee560;
ee701 = ee301 + ee345;
ee702 = ee184 + y * ((ee99 - 2 * (ee17/ee27))/ee39 + ee474 -  ee10 * (ee92 + xi * ee463 - ee90))/ee1;
ee703 = ee465 + 2 * ee309;
ee704 = ee46 - ee607;
ee705 = ee48 * ee207;
ee706 = ee233 + 2 * ee233;
ee707 = ee48/ee80;
ee708 = ee102 * ee59;
ee709 = ee608 - ee138;
ee710 = ee49 * ee207;
ee711 = ee234 + 2 * ee234;
ee712 = ee49/ee80;
ee713 = ee103 * ee59;
ee717 = (ee306 - ee298) * ee9 - ee21;
ee718 = ee437 - 1;
ee722 = 1 + ee2 * (ee440 - 3)/ee40;
ee723 = ee41 - ee60 * (ee200 + 4/ee27 - ee201 * (ee564 +  4/ee82 - ee2 * ee75/ee303)/ee1)/ee1;
ee724 = 2 * ee571;
ee725 = 2 * (ee527/ee36);
ee726 = 2 * ee575;
ee727 = 2 * (ee538/ee36);
ee728 = 2 * ee577;
ee729 = 2 * (ee540/ee36);
ee730 = 2 * ee585;
ee731 = 2 * ee587;
ee732 = ee471 + ee386;
ee733 = ee619 - ee161;
ee734 = 2 * ee545;
ee735 = ee305 + 4;
ee736 = 2 * ee548;
ee737 = 2 * ee151;
ee738 = 2 * ee190;
ee739 = 2 * ee157;
ee740 = 2 * ee158;
ee741 = 2 * ee435;
ee742 = 2 * ee602;
ee744 = 2 * (ee16/ee13) + 2 * (ee9/ee13);
ee746 = 2 * (ee15/ee13) + 2 * (ee8/ee13);
ee748 = 2 * (1 + 2 * ee3) + ee627;
ee749 = 2 * ee625;
ee750 = 2 * (ee570/ee229);
ee752 = 4 * ee539;
ee753 = 4 * ee578;
ee754 = 4 * ee368;
ee755 = 4 * ee452;
ee756 = ee565 - ee388;
ee757 = ee628 - 1;
ee758 = 4 * ee12;
ee759 = 6 * ee250;
ee760 = 6 * ee295;
ee761 = 6 * ee3;
ee762 = 8 * ee589;
ee763 = 8 * ee244;
ee764 = 8 * ee294;
ee765 = 8 * ee250;
ee766 = 8 * ee506;
ee767 = 8 * ee269;
ee768 = 8 * ee413;
ee769 = 8 * ee295;
ee770 = ee388 - ee472;
ee771 = ee475 - ee470;
ee772 = 8 * ee3;
ee773 = xi * ee10;
ee775 = ee2 * (3 - ee440)/ee40;
ee776 = ee310 - 1;
ee777 = ee570/ee1;
ee780 = y * (ee199 + ee178 + xi * (2 * ee107 - y * (ee10 *  (ee564 + xi * ee23 * ee390) + ee23/ee82)/ee1))/ee1 -  ee57;
ee781 = ee37/ee40;

out(j, 0) = -(y * ((ee423 + y * (ee172 + 2 * ee172 - ee750) * ee81/ee229)/ee36 +
   ee773 * ee722/ee5)/ee1);
out(j, 1) = -(y * (((ee668 - ee349 * (ee618 - 2 * ee183)/ee1)/ee36 +
   ee421 + ee433)/ee36 + (ee10 * (ee775 - 1) - ee776/xi)/ee5)/ee1);
out(j, 2) = -(y * (ee364 - (ee670 + ee349 * (2 * ee161 - ee619)/ee1)/ee36) * ee9/ee324);
out(j, 3) = -(y * (ee367 - (ee576 + ee349 * (ee739 - ee563)/ee1)/ee36) * ee22 * ee8/ee229);
out(j, 4) = -(y * (ee142 + (ee349 * (ee522 + 2 * ee203)/ee1 -
   ee538)/ee36 - ee376) * ee12/ee228);
out(j, 5) =  - (y * ((ee395 + (ee78 * (2 * ee173 + 2 * ee176 -
   ee618) * ee31 - ee442)/ee145)/ee36 - 2 * ((ee10 * ee616 -
   ee155)/R_pow(ee5, 2)))/ee1);
out(j, 6) = -(y * (((ee78 * ee733 - ee446) * ee31/xi - ee677)/ee36 +
   ee400 + ee374) * ee9/ee324);
out(j, 7) = -(y * (((ee78 * ee620 - ee451) * ee31/xi - ee579)/ee36 +
   ee405 + ee377) * ee22 * ee8/ee229);
out(j, 8) = -(y * (((ee78 * ee709 - ee457) * ee31/xi - ee540)/ee36 +
   ee151 - ee435) * ee12/ee228);
out(j, 9) = -(y * ((ee671 - ee399) * ee9 - ((ee81 * ee717 -
   ee445)/ee36 + ee136)) * ee9/ee324);
out(j, 10) = y * (ee73 * ee620 - ee672) * ee22 * ee9 * ee8/ee502;
out(j, 11) = -(y * (ee696 * ee161 - ((ee608 - ee385)/ee14 -
   ee76) * ee73/ee36) * ee12 * ee9/ee228);
out(j, 12) = -(y * ((ee682 - ee404) * ee8 - ((ee81 * ee613 -
   ee500)/ee36 + ee133)) * ee22 * ee8/ee229);
out(j, 13) = y * ((ee22 * (ee522 + ee439 - ee138) + ee77) * ee69/ee36 -
   ee327 * ee157) * ee12 * ee8/ee228;
out(j, 14) = y * (ee487 * ee203 + (ee81 * ee238 + ee514) * ee65/ee36) * ee12/ee228;
out(j, 15) = ((ee393 + ee660 * ee78 * ee31/ee145)/ee36 + ee182/xi)/xi +
   y * (ee179/ee39 + y * (1/ee39 + 2 * ee781)/ee40)/ee40;
out(j, 16) = (ee355 + (ee78 * (2 * ee162 - ee561) * ee63/xi -
   ee642)/ee36) * ee9/ee525;
out(j, 17) = (ee357 + (ee78 * (ee740 - ee521) * ee63/xi - ee572)/ee36) * ee22 * ee8/ee145;
out(j, 18) = (ee123 + (ee78 * (2 * ee202 - ee471) * ee63/xi -
   ee527)/ee36 - ee359) * ee12/ee256;
out(j, 19) = ((ee286 - ee674) * ee9 + ee68 - (ee361 * ee78 +
   ee448)/ee36) * ee9 * ee31/ee525;
out(j, 20) = -((ee675 + ee597 * ee73) * ee22 * ee9 * ee8 * ee31/ee567);
out(j, 21) = (ee162 * ee438 - (ee54 + (ee347 - ee607)/ee14) * ee73/ee36) * ee12 * ee9 * ee31/ee256;
out(j, 22) = ((ee293 - ee686) * ee8 + ee64 - (ee289 * ee78 +
   ee501)/ee36) * ee22 * ee8 * ee31/ee145;
out(j, 23) = -((ee327 * ee158 + ((ee46 - ee732) * ee22 - ee74) * ee69/ee36) * ee12 * ee8 * ee31/ee256);
out(j, 24) = -((ee487 * ee202 + (ee510 - ee78 * ee238) * ee65/ee36) * ee12 * ee31/ee256);
out(j, 25) = -(((ee664 - (ee361 + 2 * ee187) * ee73/ee96) * ee9 +
   ee21) * ee9/ee96);
out(j, 26) = ee361 * ee69 * ee22 * ee9 * ee8/ee249;
out(j, 27) = -((ee187 * ee438 - ee191 * (2 - ee217) * ee9/ee96) * ee12 * ee9/ee180);
out(j, 28) = ee289 * ee73 * ee22 * ee9 * ee8/ee249;
out(j, 29) = ee73 * ee69 * ee519 * ee12 * ee9 * ee8/ee503;
out(j, 30) = ee574 * ee73 * ee12 * ee9/ee180;
out(j, 31) = -(((ee678 - (ee289 + ee738) * ee69 * ee22/ee36) * ee8 +
   ee18) * ee22 * ee8/ee36);
out(j, 32) = (ee327 * ee190 - ee166 * ee387 * ee8/ee36) * ee12 * ee8/ee180;
out(j, 33) = -(ee533 * ee69 * ee12 * ee8/ee180);
out(j, 34) = -(ee529 * ee65 * ee12/ee180);
out(j, 35) =  - (y * (((ee25 * ee723 + y * (3 * (ee260 * ee66) -
   ((ee48 * ee210 - y * ((ee604 + y * ee515 * ee100/ee83)/ee13 +
   2 * (ee604/ee13)) * ee43/ee1)/ee13 + ee706 * ee210)) * ee20/ee1) * ee9 * ee16/ee14 +
   (ee24 * ee723 + y * (3 * (ee264 * ee66) -
   ((ee49 * ee210 - y * ((ee609 + y * ee516 * ee101/ee83)/ee13 +
   2 * (ee609/ee13)) * ee44/ee1)/ee13 + ee711 * ee210)) * ee19/ee1) * ee22 * ee8 * ee15 +
   y * (ee81 * (2 * ee423 +
   ee349 * ((2 * (ee665 + ee777) + 4 * ee665 - 8 * ee777)/ee36 +
   4 * ee172)/ee229) - (ee172 * (ee172 - ee750) + 2 * (R_pow(ee172, 2) -
   ee423 * ee81)))/ee229)/ee36 + ee773 * (ee2 * (7 +
   ee2 * ((ee772 - ee748)/ee5 - 6)/ee40)/ee40 - 1)/ee5)/ee1);
out(j, 36) = -(y * (((ee423 * ee78 * ee31/xi - y * (ee172 * (ee183 -
   ee618) + 2 * (ee172 * ee183 - ee459 * ee81) - (ee459 +
   (4 * ee668 + ee349 * (2 * (ee578 + ee589) +   ee753 - ee762)/ee229)/ee36) * ee81)/ee1)/ee36 +
   (ee25 * ee780 + ee20 * (y * (2 * (ee263 * ee66) -
     (ee260 * ee59 + (ee705 + ee43 * (y * (ee708/ee13 -
   (ee555 * ee31/ee50 - ee708)/ee13)/ee1 -
   ee604 * ee31/xi))/ee13 + 2 * (ee705/ee13)))/ee1 -   ee275 * ee210/xi)) * ee9 * ee16/ee14 +
   (ee24 * ee780 + ee19 * (y * (2 * (ee267 * ee66) -
   (ee264 * ee59 + (ee710 + ee44 * (y * (ee713/ee13 -
   (ee557 * ee31/ee50 - ee713)/ee13)/ee1 - ee609 * ee31/xi))/ee13 +
   2 * (ee710/ee13)))/ee1 - ee280 * ee210/xi)) * ee22 * ee8 * ee15)/ee36 +
   (ee10 * (1 + ee2 * (ee2 * ((ee748 -
   ee772)/ee5 + 6)/ee40 - 7)/ee40) - ee722/xi)/ee5)/ee1);
out(j, 37) = -(y * ((ee105 * ee210 - y * ((ee211 * ee66 + y * ((ee554/ee13 +
   ee102/ee13) * ee20/ee13 + ee143/ee80)/ee1) * ee16/ee13 +
   ee260 * ee20/ee13 + 2 * (ee699 * ee66))/ee1) * ee9 +
   ee339 * ee16 - (ee423 * ee73 + y * (ee172 * ee733 + ee81 * ((4 * ee670 -
   ee349 *   (ee768 - (2 * (ee539 + ee413) +
   ee752))/ee229)/ee36 -   ee364) - 2 * (ee364 * ee81 + ee172 * ee161))/ee1)/ee36) * ee9/ee324);
out(j, 38) = -(y * ((ee104 * ee210 - y * ((ee212 * ee66 + y * ((ee556/ee13 +
   ee103/ee13) * ee19/ee13 + ee144/ee80)/ee1) * ee15/ee13 +
   ee264 * ee19/ee13 + 2 * (ee701 * ee66))/ee1) * ee8 +
   ee344 * ee15 - (ee423 * ee69 + y * (ee172 * ee620 + ee81 * ((4 * ee576 -
   ee349 * (ee769 - (2 * (ee407 + ee295) +
   ee626))/ee229)/ee36 - ee367) - 2 * (ee367 * ee81 + ee172 * ee157))/ee1)/ee36) * ee22 * ee8/ee229);
out(j, 39) = -(y * (ee373 - ((ee423 * ee65 + y * (ee172 * ee709 +
   (ee376 + (4 * ee538 - ee349 * (2 * (ee452 - ee269) + ee755 +
     ee767)/ee229)/ee36 - ee142) * ee81 + 2 * (ee172 * ee203 -
   ee547 * ee81))/ee1)/ee36 + ee344 * ee8 * ee15)) * ee12/ee228);
out(j, 40) =  - (y * ((((ee661 * ee172 + 2 * (ee459 * ee78 * ee31))/xi -
   y * (((ee78 * (ee762 - ee753) * ee31 + 2 * (ee571 * ee81))/ee36 -
   4 * ee450) * ee81/ee145 + 2 * (R_pow(ee183, 2) -
   ee395 * ee81))/ee1)/ee36 + (ee702 * ee25 + ee20 * (y * (ee687/ee13 -
   ((((ee62 * ee102 - ee606)/ee13 + 2 * (ee102 * ee31 * ee59)) * ee43/xi -
   ee687)/ee13 + 2 * (ee263 * ee59)))/ee1 -
   (ee241 * ee66 + 2 * (ee275 * ee207))/xi)) * ee9 * ee16/ee14 +
   (ee702 * ee24 + ee19 * (y * (ee688/ee13 - ((((ee62 * ee103 -
   ee611)/ee13 + 2 * (ee103 * ee31 * ee59)) * ee44/xi -
   ee688)/ee13 + 2 * (ee267 * ee59)))/ee1 - (ee243 * ee66 +
   2 * (ee280 * ee207))/xi)) * ee22 * ee8 * ee15)/ee36 + (ee37 * (4 -
   ee2 * ((ee627 - ee761)/ee5 + 6)/ee40)/ee40 - (2 * ee775 -
   (2 + 2 * ee776))/ee39)/ee5)/ee1);
out(j, 41) = -(y * (((ee364 * ee78 + ee172 * ee697) * ee31/xi -
   (ee459 * ee73 + y * (((ee78 * (ee752 - ee768) + 2 * (ee577 * ee81)) * ee31/ee145 +
   4 * ee677) * ee81/ee36 - 2 * (ee483 * ee81 +
   ee161 * ee183))/ee1))/ee36 + (ee105 * ee207 + y * (ee699 * ee59 -
   (((ee596/ee13 +   ee143/ee13) * ee31/xi -
   ee211 * ee59) * ee16/ee13 +   ee263 * ee20/ee13))/ee1 - ee270 * ee66 * ee31/xi) *   ee9 +
   ee337 * ee16) * ee9/ee324);
out(j, 42) = -(y * (((ee367 * ee78 + ee172 * ee597) * ee31/xi -
   (ee459 * ee69 + y * (((ee78 * (ee626 - ee769) + 2 * (ee585 * ee81)) * ee31/ee145 +
   4 * ee579) * ee81/ee36 - 2 * (ee486 * ee81 +
   ee183 * ee157))/ee1))/ee36 + (ee104 * ee207 + y * (ee701 * ee59 -
   (((ee599/ee13 +   ee144/ee13) * ee31/xi -
   ee212 * ee59) * ee15/ee13 +   ee267 * ee19/ee13))/ee1 - ee271 * ee66 * ee31/xi) *   ee8 +
   ee342 * ee15) * ee22 * ee8/ee229);
out(j, 43) = -(y * (((ee547 * ee78 + ee172 * ee704) * ee31/xi -
   (ee459 * ee65 + y * (ee81 * ((2 * (ee587 * ee81) - ee78 * (ee755 +
   ee767)) * ee31/ee145 + 4 * ee540)/ee36 + 2 * (ee183 * ee203 -
   ee594 * ee81))/ee1))/ee36 + ee371 - ee342 * ee8 * ee15) * ee12/ee228);
out(j, 44) =  - (y * ((ee286 * ee66 + y * (((ee595/ee13 + 2 * (ee468/ee13)) * ee28 +
   ee233) * ee16/ee13 + ee274 * ee744 +
   (ee707 + 2 * ee707) * ee20)/ee1 - ee644/ee96) * ee9 + ee328 -
   (ee361 * ee172 + ee644 * ee9/ee14 + y * (2 * (ee477 * ee81 +
   R_pow(ee161, 2) * ee9/ee14) - ((ee81 * ee617 + 4 * (ee539 * ee73 * ee9/ee14))/ee36 +
   4 * ee445) * ee81/ee36)/ee1)/ee36) * ee9/ee324);
out(j, 45) = y * (ee364 * ee69 + (ee367 - 2 * (ee576/ee36)) * ee73 +
   y * (ee161 * (ee739 - 4 * ee508) - ee413 * (ee626 -
   ee760)/ee132)/ee1) * ee22 * ee9 * ee8/ee502;
out(j, 46) = -(y * (ee364 * ee438 - ((ee169 + (2 * ee142 - (ee376 +
   ee727))/ee14) * ee73 + y * (ee161 * (2 * (ee76 + (ee385 -
   ee140)/ee14) - 4 * (ee269/ee96)) - (ee81 * ee621 - 4 * (ee452/ee14)) * ee81 * ee73/ee132)/ee1)/ee36) * ee12 * ee9/ee228);
out(j, 47) =  - (y * ((ee293 * ee66 + y * (((ee598/ee13 + 2 * (ee469/ee13)) * ee28 +
   ee234) * ee15/ee13 + ee279 * ee746 +
   (ee712 + 2 * ee712) * ee19)/ee1 - ee647/ee36) * ee8 + ee332 -
   (ee289 * ee172 + ee647 * ee8 + y * (2 * (ee479 * ee81 +
   R_pow(ee157, 2) * ee22 * ee8) - ((ee81 * ee520 + 4 * ee681)/ee36 +
   4 * ee500) * ee81/ee36)/ee1)/ee36) * ee22 * ee8/ee229);
out(j, 48) = y * ((((ee142 - (ee727 + 2 * ee376)) * ee22 - ee165) * ee69 -
   y * ((ee81 * ee523 - 4 * (ee452 * ee22)) * ee81 * ee69/ee132 +
   ee157 * (2 * ee467 + 4 * (ee269 * ee22/ee36)))/ee1)/ee36 +
   ee367 * ee327) * ee12 * ee8/ee228;
out(j, 49) = y * (((ee547 * ee12/ee86 - ee172 * ee238) * ee65 +
   y * (ee203 * (ee734 + 4 * (ee269 * ee12/ee180)) - (ee81 * ee384 -
   4 * (ee452 * ee12/ee86)) * ee81 * ee65/ee132)/ee1)/ee36 +
   ee487 * ee547) * ee12/ee228;
out(j, 50) = -(y * (((((ee633 - ee654)/ee13 + ee192 * ee276 * ee31 +
   3 * (ee241 * ee59)) * ee20/xi - ee650 * ee25) * ee9 * ee16/ee14 +
   (((ee634 - ee658)/ee13 + ee192 * ee281 * ee31 +
   3 * (ee243 * ee59)) * ee19/xi - ee650 * ee24) * ee22 * ee8 * ee15 +
   ((ee395 + ((2 * (ee494 * ee36) - ee81 * (ee724 +
   ee766))/ee36 +   2 * ee494)/ee145 + 2 * ee395) * ee78 * ee31 +
   ee660 * ee183 - ee393 * ee81)/ee145)/ee36 + (((6 * (1 - ee310) -
   6)/xi + 2 * (y * ee616/ee40))/ee39 + y * (ee616/ee39 +
   y * ((2 * ee5 - ee761)/ee5 +   4) * ee10/ee40)/ee40)/ee5)/ee1);
out(j, 51) = -(y * ((((((ee81 * (ee764 - ee728) - 2 * ee643) * ee31/ee391 +
   2 * ee400 + 2 * ee374) * ee78 + ee183 * (2 * ee159 +
   2 * ee68 - ee561)) * ee31 - (ee355 * ee81 + ee153 * ee161 +
   2 * (ee494 * ee73/ee36)))/xi - ee395 * ee73)/ee36 +
   ((ee651/ee13 + ee662/ee13 + 2 * (ee270 * ee31 * ee59))/xi -
   ee192 * ee105) * ee9 + ee317 * ee16) * ee9/ee324);
out(j, 52) = -(y * ((((((ee81 * (ee765 - ee730) - ee726) * ee31/ee391 +
   2 * ee405 + 2 * ee377) * ee78 + ee183 * (2 * ee154 +
   2 * ee64 - ee521)) * ee31 - (ee357 * ee81 + ee153 * ee157 +
   2 * (ee494 * ee69/ee36)))/xi - ee395 * ee69)/ee36 + ((ee655/ee13 +
   ee663/ee13 + 2 * (ee271 * ee31 * ee59))/xi - ee192 * ee104) * ee8 +
   ee320 * ee15) * ee22 * ee8/ee229);
out(j, 53) = -(y * ((((((ee81 * (ee763 - ee731) + 2 * ee531) * ee31/ee391 +
   ee737 - ee741) * ee78 + ee183 * (ee347 - ee732)) * ee31 +
   ee153 * ee203 - (ee530 * ee81 + 2 * (ee494 * ee65/ee36)))/xi -
   ee395 * ee65)/ee36 + ee353 - ee320 * ee8 * ee15) * ee12/ee228);
out(j, 54) = -(y * (((ee645/ee13 + ee211 * ee744 + ee706 * ee20) * ee31/xi -
   (ee635/ee96 + ee286 * ee59)) * ee9 + (((2 * (ee492 * ee73/ee36) -
   2 * (ee161 * ee162)) * ee9/ee14 - (ee477 * ee78 +
   ee482 * ee81 + (ee258 * ee771 -   2 * (ee643 * ee73 * ee9/ee14))/ee132)) * ee31/xi -
   (ee635 * ee9/ee14 + ee361 * ee183))/ee36 +
   ee374) * ee9/ee324);
out(j, 55) = -(y * (((2 * (ee492 * ee69) - ee73 * (6 * (ee258 * ee69) -
   ee726)/ee36)/ee36 - (ee161 * ee158 + ee162 * ee157)) * ee31/xi -
   (ee483 * ee69 + (ee486 - 2 * (ee579/ee36)) * ee73)) * ee22 * ee9 * ee8/ee502);
out(j, 56) = -(y * (((ee162 * ee461 + (2 * (ee492 * ee65/ee14) -
   (ee258 * ee770 + 2 * (ee531/ee14)) * ee73/ee36)/ee36 -
   ee161 * ee436) * ee31/xi - (ee176 + (ee737 - (ee435 + ee729))/ee14) * ee73)/ee36 +
   ee483 * ee438) * ee12 * ee9/ee228);
out(j, 57) = -(y * (((ee648/ee13 + ee212 * ee746 + ee711 * ee19) * ee31/xi -
   (ee637/ee36 + ee293 * ee59)) * ee8 + ((ee22 * (2 * (ee499 * ee69/ee36) -
   2 * (ee157 * ee158)) * ee8 - (ee479 * ee78 +
   ee485 * ee81 + (ee258 * ee630 - 2 * (ee575 * ee69 * ee22 * ee8))/ee132)) * ee31/xi -
   (ee637 * ee8 + ee289 * ee183))/ee36 +
   ee377) * ee22 * ee8/ee229);
out(j, 58) = -(y * (((ee157 * (ee22 * (ee386 - ee46) + ee74) +
   ee158 * ee467 + (2 * (ee499 * ee65 * ee22) - (ee258 * ee629 +
   2 * (ee531 * ee22)) * ee69/ee36)/ee36) * ee31/xi - ((ee151 -
   (ee729 + ee741)) * ee22 - ee173) * ee69)/ee36 - ee486 * ee327) * ee12 * ee8/ee228);
out(j, 59) = -(y * (((ee545 * ee202 + ee548 * ee203 - ((ee258 * ee566 +
   2 * (ee531 * ee12/ee86))/ee36 + 2 * (ee505 * ee12/ee86)) * ee65/ee36) * ee31/xi -
   (ee594 * ee12/ee86 - ee183 * ee238) * ee65)/ee36 -
   ee487 * ee594) * ee12/ee228);
out(j, 60) = -(y * (((2 * (ee477 * ee73) - ee161 * ee717)/ee96 -
   ((((ee272 * ee9 * ee28 + ee25)/ee13 + ee698 + ee98 + ee98 +
   ee98) * ee16 + (ee698 + ee98 + ee98 + ee98 + ee622 + ee622 +
   ee622) * ee9) * ee28 + ee98 + ee98 + ee98 + ee98 + ee98 +
   ee98 + ee98)) * ee9 - (((ee73 * ((2 * ee537 - (ee81 * ee771 -
   2 * (ee537 *   ee36))/ee36)/ee36 - ee477) - 2 * (ee187 * ee161)) *   ee9/ee14 -
   ee490 * ee81)/ee36 + ee136)) * ee9/ee324);
out(j, 61) = y * (ee69 * ((2 * ee671 - ee399) * ee9 - ee136) +
   (2 * (ee537 * ee69) - ee191 * (ee760 - 2 * ee407) * ee9/ee96)/ee36 -
   ee187 * ee157) * ee22 * ee9 * ee8/ee502;
out(j, 62) = -(y * (ee477 * ee696 - (((ee161 * (ee217 - 4) -
   (ee81 * ee770 + 2 * ee584) * ee73/ee132) * ee73 * ee9 + 2 * (ee537 * ee65/ee36))/ee14 +
   ee187 * ee461)/ee36) * ee12 * ee9/ee228);
out(j, 63) = y * (ee161 * ee613 + ee73 * ((2 * ee682 - ee404) * ee8 -
   ((ee695 - 2 * ee681)/ee132 + ee133))) * ee22 * ee9 * ee8/ee502;
out(j, 64) = y * ((ee157 * ee718 - (ee693/ee14 + 2 * (ee584 * ee22)) * ee69/ee132) * ee73 +
   ee672 * ee718) * ee12 * ee9 * ee8/(ee503 * ee1);
out(j, 65) = y * (((ee65 * ((ee217 - 3) * ee12/ee14 + 1)/ee36 -
   ee248)/ee14 - 1) * ee161 + ((ee203 * ee757 - ((ee694 + 2 * (ee584 * ee12/ee14))/ee36 +
   2 * (ee461 * ee12/ee14)) * ee65/ee36)/ee14 +
   ee458) * ee73/ee36) * ee12 * ee9/ee228;
out(j, 66) = -(y * ((ee22 * (2 * (ee479 * ee69) - ee157 * ee613)/ee36 -
   ((((ee277 * ee8 * ee28 + ee24)/ee13 + ee700 + ee95 +
   ee95 + ee95) * ee15 + (ee700 + ee95 + ee95 + ee95 + ee623 +
   ee623 + ee623) * ee8) * ee28 + ee95 + ee95 + ee95 + ee95 +
   ee95 + ee95 + ee95)) * ee8 - (((ee69 * ((2 * ee583 - (ee695 -
   2 * (ee583 * ee36))/ee36)/ee36 - ee479) - 2 * (ee190 * ee157)) * ee22 * ee8 -
   ee497 * ee81)/ee36 + ee133)) * ee22 * ee8/ee229);
out(j, 67) = y * ((((ee157 * ee735 - (ee693 + 2 * ee689) * ee69/ee132) * ee69 * ee8 +
   2 * (ee583 * ee65/ee36)) * ee22 +
   ee190 * ee467)/ee36 - ee479 * ee327) * ee12 * ee8/ee228;
out(j, 68) = -(y * ((((ee81 * ee22 * ee566 + 2 * (ee689 * ee12/ee86))/ee36 +
   2 * (ee467 * ee12/ee86)) * ee65/ee36 + ee458 +
   ee703 * ee203) * ee69/ee36 - ee533 * ee157) * ee12 * ee8/ee228);
out(j, 69) = -(y * (((((ee694 + 2 * (ee545 * ee36))/ee36 + ee734) * ee65/ee36 +
   3 * (ee203 * ee134)) * ee12/ee86 + ee81 * ee614) * ee65/ee36 +
   ee529 * ee203) * ee12/ee228);
out(j, 70) = (((ee632 * ee25 - ((4 * ee633 - ee654) * ee31 +
   3 * (ee241 * ee62)) * ee20/xi) * ee9 * ee16/ee14 + (ee632 * ee24 -
   ((4 * ee634 - ee658) * ee31 + 3 * (ee243 * ee62)) * ee19/xi) * ee22 * ee8 * ee15 +
   (ee153 * ee661 + (ee78 * ((ee724 +
   4 * ee526 + ee766)/ee36 + 4 * ee153) * ee31/ee145 + 2 * ee393) * ee78 * ee31 +
   2 * (ee393 * ee78 * ee31 + R_pow(ee153, 2)))/ee145)/ee36 +
   ee601/xi)/xi - y * (ee559/ee39 + y * (ee524/ee39 +
   y * (2/ee39 + 6 * ee781)/ee40)/ee40)/ee40;
out(j, 71) = ((((ee355 + (ee78 * (ee728 + 4 * ee493 - ee764) * ee63/ee145 -
   4 * ee642)/ee36) * ee78 + ee153 * ee697 + 2 * (ee355 * ee78 +
   ee153 * ee162)) * ee31/xi - ee393 * ee73)/ee36 +
   (ee205 * ee105 - (ee651 + ee662 + 2 * (ee270 * ee62)) * ee31/xi) * ee9 +
   ee312 * ee16) * ee9/ee525;
out(j, 72) = ((((ee357 + (ee78 * (ee730 + ee754 - ee765) * ee63/ee145 -
   4 * ee572)/ee36) * ee78 + ee153 * ee597 + 2 * (ee357 * ee78 +
   ee153 * ee158)) * ee31/xi - ee393 * ee69)/ee36 +
   (ee205 * ee104 - (ee655 + ee663 + 2 * (ee271 * ee62)) * ee31/xi) * ee8 +
   ee314 * ee15) * ee22 * ee8/ee145;
out(j, 73) = ((((ee123 + (ee78 * (ee731 + 4 * ee412 - ee763) * ee63/ee145 -
   4 * ee527)/ee36 - ee359) * ee78 + ee153 * ee704 +
   2 * (ee530 * ee78 + ee153 * ee202)) * ee31/xi - ee393 * ee65)/ee36 +
   ee351 - ee314 * ee8 * ee15) * ee12/ee256;
out(j, 74) = ((ee286 * ee62 - ((ee645 + ee211 * (4 * ee9 - 2) +
   ee276 * ee20) * ee63/xi + ee638/ee96)) * ee9 + ee321 + ((2 * (ee482 * ee78 +
   R_pow(ee162, 2) * ee9/ee14) - ((ee691 +
   4 * (ee493 * ee73 * ee9/ee14))/ee36 + 4 * ee448) * ee78/ee36) * ee63/xi -
   (ee638 * ee9/ee14 + ee153 * ee361))/ee36) * ee9/ee525;
out(j, 75) = ((ee162 * (ee740 - 4 * ee455) - ee294 * (ee754 -
   ee759)/ee132) * ee63/xi - (ee355 * ee69 + (ee357 - 2 * (ee572/ee36)) * ee73)) * ee22 * ee9 * ee8/ee567;
out(j, 76) = (((ee162 * (2 * ee436 - 4 * (ee244/ee96)) - (ee692 +
   4 * (ee412/ee14)) * ee78 * ee73/ee132) * ee63/xi - (ee149 +
   (2 * ee123 - (ee359 + ee725))/ee14) * ee73)/ee36 + ee355 * ee438) * ee12 * ee9/ee256;
out(j, 77) = ((ee293 * ee62 - ((ee648 + ee212 * (4 * ee8 - 2) +
   ee281 * ee19) * ee63/xi + ee640/ee36)) * ee8 + ee322 + ((2 * (ee485 * ee78 +
   R_pow(ee158, 2) * ee22 * ee8) - ((ee590 +
   4 * ee685)/ee36 + 4 * ee501) * ee78/ee36) * ee63/xi - (ee640 * ee8 +
   ee153 * ee289))/ee36) * ee22 * ee8/ee145;
out(j, 78) = (((ee158 * (2 * ee512 - 4 * (ee244 * ee22/ee36)) -
   (ee592 + 4 * (ee412 * ee22)) * ee78 * ee69/ee132) * ee63/xi -
   ((ee123 - (ee725 + 2 * ee359)) * ee22 - ee146) * ee69)/ee36 -
   ee357 * ee327) * ee12 * ee8/ee256;
out(j, 79) = (((ee202 * (ee736 - 4 * (ee244 * ee12/ee180)) -
   (ee591 + 4 * (ee412 * ee12/ee86)) * ee78 * ee65/ee132) * ee63/xi -
   (ee530 * ee12/ee86 - ee153 * ee238) * ee65)/ee36 - ee530 * ee487) * ee12/ee256;
out(j, 80) = (((((ee220 + ee218 + ee218 + ee218) * ee9 * ee28 +
   ee25 + ee25 + ee25 + ee218 + ee218) * ee16 + ((ee25 + ee218) * ee9 * ee28 +
   ee25 + ee25 + ee218 + ee218 + ee218 + ee218 +
   ee218) * ee9) * ee28 + ee25 + ee25 + ee25 + ee25 + ee25 +
   ee25 + ee25 - (ee361 * ee162 + 2 * (ee482 * ee73))/ee96) * ee9 +
   ee68 - (((ee482 - ((ee691 + 2 * (ee535 * ee36))/ee36 +
   2 * ee535)/ee36) * ee73 + 2 * (ee187 * ee162)) * ee9/ee14 +
   ee490 * ee78)/ee36) * ee9 * ee31/ee525;
out(j, 81) = -((((ee286 - 2 * ee674) * ee9 + ee68) * ee69 +
   ee187 * ee158 - (ee191 * (2 * ee368 - ee759) * ee9/ee96 + 2 * (ee535 * ee69))/ee36) * ee22 * ee9 * ee8 * ee31/ee567);
out(j, 82) = (ee482 * ee438 - (((ee162 * (4 - ee217) - (ee692 +
   2 * ee586) * ee73/ee132) * ee73 * ee9 - 2 * (ee535 * ee65/ee36))/ee14 +
   ee187 * ee436)/ee36) * ee12 * ee9 * ee31/ee256;
out(j, 83) = -((((ee293 - 2 * ee686) * ee8 + ee64 - (ee590 +
   2 * ee685)/ee132) * ee73 + ee289 * ee162) * ee22 * ee9 * ee8 * ee31/ee567);
out(j, 84) = -(((ee158 * ee519 - (ee592/ee14 + 2 * (ee586 * ee22)) * ee69/ee132) * ee73 +
   ee675 * ee519) * ee12 * ee9 * ee8 * ee31/(ee391 * ee86));
out(j, 85) = -((ee574 * ee162 + ((ee202 * ee757 - ((ee591 +
   2 * (ee586 * ee12/ee14))/ee36 + 2 * (ee436 * ee12/ee14)) * ee65/ee36)/ee14 -
   ee456) * ee73/ee36) * ee12 * ee9 * ee31/ee256);
out(j, 86) = (((((ee222 + ee219 + ee219 + ee219) * ee8 * ee28 +
   ee24 + ee24 + ee24 + ee219 + ee219) * ee15 + ((ee24 + ee219) * ee8 * ee28 +
   ee24 + ee24 + ee219 + ee219 + ee219 + ee219 +
   ee219) * ee8) * ee28 + ee24 + ee24 + ee24 + ee24 + ee24 +
   ee24 + ee24 - (ee289 * ee158 + 2 * (ee485 * ee69)) * ee22/ee36) * ee8 +
   ee64 - (((ee485 - ((ee590 + 2 * (ee581 * ee36))/ee36 +
   2 * ee581)/ee36) * ee69 + 2 * (ee190 * ee158)) * ee22 * ee8 +
   ee497 * ee78)/ee36) * ee22 * ee8 * ee31/ee145;
out(j, 87) = -((ee485 * ee327 + (ee512 * ee190 - (((ee592 +
   2 * ee673) * ee69/ee132 + ee158 * ee735) * ee69 * ee8 + 2 * (ee581 * ee65/ee36)) * ee22)/ee36) * ee12 * ee8 * ee31/ee256);
out(j, 88) = ((((ee78 * ee22 * ee384 + 2 * (ee673 * ee12/ee86))/ee36 +
   2 * (ee512 * ee12/ee86)) * ee65/ee36 + ee703 * ee202 -
   ee456) * ee69/ee36 + ee533 * ee158) * ee12 * ee8 * ee31/ee256;
out(j, 89) = (((((ee591 + 2 * (ee548 * ee36))/ee36 + ee736) * ee65/ee36 +
   3 * (ee202 * ee134)) * ee12/ee86 - ee78 * ee614) * ee65/ee36 +
   ee529 * ee202) * ee12 * ee31/ee256;
out(j, 90) =  - ((((((ee380 + ee93 + ee93) * ee9 * ee28 + ee21 +
   ee21 + ee21 + ee21 + ee21 + ee93 + ee93 + ee93 + ee93 +
   ee93 + ee93 + ee93 + ee93 + ee93 + ee93) * ee9 * ee28 + ee21 +
   ee21 + ee21 + ee21 + ee21 + ee21 + ee21 + ee21 + ee21 +
   ee93 + ee93 + ee93) * ee28 - (ee361 * ee187 + ee73 * (2 * ee490 -
   ee73 * ((ee470 + 4 * ee406 - ee475)/ee36 + 4 * ee187) * ee9/ee96) +
   2 * (ee490 * ee73 + R_pow(ee187, 2)))/ee96) * ee9 +
   ee21) * ee9/ee96);
out(j, 91) = ((ee664 - ee73 * (6 * ee187 - 6 * (ee251/ee96))/ee96) * ee9 +
   ee21) * ee69 * ee22 * ee9 * ee8/ee249;
out(j, 92) = -((ee490 * ee438 - (ee187 * (6 - 6 * ee135) - ee191 * (ee472 +
   ee565 - ee388) * ee9/ee249) * ee73 * ee9/ee96) * ee12 * ee9/ee180);
out(j, 93) = (ee289 * ee187 - ee191 * ee520 * ee9/ee249) * ee22 * ee9 * ee8/ee249;
out(j, 94) = (ee187 * ee519 - (ee22 * ee756 + 2 * (ee297/ee14)) * ee191 * ee9/ee249) * ee69 * ee12 * ee9 * ee8/ee503;
out(j, 95) = (ee574 * ee187 + ee191 * (ee749 - ((ee756 * ee12/ee14 +
   ee239)/ee36 + ee628) * ee65/ee96) * ee9/ee96) * ee12 * ee9/ee180;
out(j, 96) = ((ee678 - ee69 * ((ee304 + 2 * ee257 - ee307)/ee36 +
   ee738) * ee22/ee36) * ee8 + ee18) * ee73 * ee22 * ee9 * ee8/ee249;
out(j, 97) = (ee190 * ee519 - ee163 * (ee523/ee14 + 2 * ee690) * ee22 * ee8/ee132) * ee73 * ee12 * ee9 * ee8/ee503;
out(j, 98) = -((((ee603 + 2 * (ee690 * ee12/ee14)) * ee65/ee132 +
   (2 * (ee65 * ee615/ee36) + 4) * ee12/ee14 - 1)/ee14 +
   ee465) * ee73 * ee69 * ee12 * ee9 * ee8/ee503);
out(j, 99) = -((((((ee65 * ((2 * (ee36 * ee625) + ee239 - ee246)/ee36 +
   ee749)/ee36 + ee255 + ee268 - 10 * ee12)/ee14 +
   ee473 + 3) * ee12/ee14 - 1) * ee65/ee36 + ee480)/ee14 + 1) * ee73 * ee12 * ee9/ee180);
out(j, 100) =  - ((((((ee383 + ee91 + ee91) * ee8 * ee28 + ee18 +
   ee18 + ee18 + ee18 + ee18 + ee91 + ee91 + ee91 + ee91 +
   ee91 + ee91 + ee91 + ee91 + ee91 + ee91) * ee8 * ee28 + ee18 +
   ee18 + ee18 + ee18 + ee18 + ee18 + ee18 + ee18 + ee18 +
   ee91 + ee91 + ee91) * ee28 - (ee289 * ee190 + ee69 * (2 * ee497 -
   ee69 * ((ee304 + 4 * ee257 - ee307)/ee36 + 4 * ee190) * ee22 * ee8/ee36) +
   2 * (ee497 * ee69 + R_pow(ee190, 2))) * ee22/ee36) * ee8 +
   ee18) * ee22 * ee8/ee36);
out(j, 101) = (ee497 * ee327 - (ee190 * (6 + 6 * ee224) + ee166 * (ee346 -
   (ee565 + ee348)) * ee8/ee132) * ee69 * ee22 * ee8/ee36) * ee12 * ee8/ee180;
out(j, 102) = -(((((ee603 - 4 * (ee36 * ee12/ee86))/ee36 - 4 * ee309) * ee65 * ee22/ee36 -
   ee742) * ee163 * ee8/ee36 + ee533 * ee190) * ee12 * ee8/ee180);
out(j, 103) = -((((((ee603 - 2 * (ee36 * ee602))/ee36 - ee742) * ee65/ee36 -
   3 * ee134) * ee12/ee86 - ee614 * ee22) * ee65/ee36 -
   ee529) * ee69 * ee12 * ee8/ee180);
out(j, 104) =  - ((((((ee65 * ((ee239 - (4 * ee232 + ee246))/ee36 -
   4 * ee134)/ee96 + 2 * ee518) * ee12/ee14 - (ee238 * ee134 +
   2 + 2 * (R_pow(ee134, 2) + 1 - ee550))) * ee65/ee36 +
   2 * (1 + ee758) + 4 * ee14 + 6 * ee245 - ((4 * (ee255 + ee758) +
   8 * ee282 - 64 * ee12)/ee14 + 8) * ee12)/ee14 + 2) * ee12/ee14 -
   1) * ee65 * ee12/ee180);

}

return out;

}

// //' Extended generalized Pareto distribution of type 3 (eGPD3) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each eGPD parameter
// //' @param X1 a design matrix for the eGPD log scale parameter
// //' @param X2 a design matrix for the eGPD shape parameter
// //' @param X3 a design matrix for the eGPD log kappa
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return egpd3d0 a scalar, the negative log-liklihood
// //' @return egpd3d12 a matrix, first then second derivatives w.r.t. eGPD1 parameters
// //' @return egpd3d34 a matrix, third then fourth derivatives w.r.t. eGPD1 parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double egpd3d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, arma::vec yvec, const arma::uvec& dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec ldeltavec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();

if (dcate == 1) {
  lpsivec = lpsivec.elem(dupid);
  xivec = xivec.elem(dupid);
  ldeltavec = ldeltavec.elem(dupid);
}

double y, lpsi, xi, ldelta;
double ee1, ee4;
double nllh=0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
ldelta = ldeltavec[j];

// ee1 = exp(ldelta);
// ee2 = exp(lpsi);
// ee3 = 1 + xi * y/ee2;
// ee4 = 1/xi;
// ee6 = 1 + ee1;
// ee7 = 1/R_pow(ee3, ee4);
// 
// nllh -= log(-((R_pow(ee7, ee1 - 1)/R_pow(ee3, 1 + 2/xi) - (1 - R_pow(ee7, ee1)/ee6) * 
// ee6/(R_pow(ee3, 1 + ee4) * ee1))/ee2));

ee1 = exp(ldelta);
ee4 = xi * y/exp(lpsi);

if (ee4 <= -1.0) {
    nllh = 1e20;
    break;
}

nllh += (1 + 1/xi) * log1p(ee4) + ldelta + lpsi - (log(1 - 1/R_pow(1 + ee4, ee1/xi)) + log1p(ee1));
    
}

return(nllh);

}

// //' @rdname egpd3d0
// [[Rcpp::export]]
arma::mat egpd3d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec, const arma::uvec dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec ldeltavec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 9);

if (dcate == 1) {
  lpsivec = lpsivec.elem(dupid);
  xivec = xivec.elem(dupid);
  ldeltavec = ldeltavec.elem(dupid);
}

double y, lpsi, xi, ldelta;

// double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
// double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
// double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
// double ee30, ee32, ee33, ee34, ee35, ee36, ee37, ee38, ee39;
// double ee40, ee49;
// double ee53, ee54, ee55, ee56, ee58, ee59;
// double ee60, ee61, ee62, ee63, ee65, ee66, ee67, ee68, ee69;
// double ee70, ee71, ee74, ee75, ee78;
// double ee80, ee81, ee82, ee83, ee84, ee85;
// double ee92, ee93, ee94, ee97, ee98;
// double ee100, ee101, ee103, ee104, ee105, ee109;
// double ee113, ee116, ee117, ee118;
// double ee122, ee124, ee126, ee127, ee128, ee129;
// double ee130, ee131;

double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
double ee21, ee22, ee23, ee25, ee26, ee27, ee28;
double ee30, ee31, ee32, ee33, ee34, ee36, ee37, ee38;
double ee40;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
ldelta = ldeltavec[j];

// ee1 = exp(lpsi);
// ee2 = xi * y;
// ee3 = ee2/ee1;
// ee4 = 1 + ee3;
// ee5 = 1/xi;
// ee6 = exp(ldelta);
// ee7 = R_pow(ee4, ee5);
// ee8 = 1/ee7;
// ee9 = 2/xi;
// ee10 = 1 + ee5;
// ee11 = 1 + ee9;
// ee12 = ee6 - 1;
// ee13 = 1 + ee6;
// ee14 = R_pow(ee4, ee10);
// ee15 = log1p(ee3);
// ee16 = R_pow(ee8, ee12);
// ee17 = R_pow(ee8, ee6);
// ee18 = ee14 * ee6;
// ee19 = R_pow(ee4, ee11);
// ee20 = ee17/ee13;
// ee21 = 1 - ee20;
// ee22 = ee5 - 1;
// ee23 = R_pow(ee4, ee22);
// ee24 = ee21 * ee13;
// ee25 = R_pow(ee4, ee9);
// ee26 = R_pow(ee18, 2);
// ee27 = ee6 - 2;
// ee28 = R_pow(ee8, ee27);
// ee29 = R_pow(xi, 2);
// ee30 = y * ee23;
// ee32 = ee7 * ee15/xi;
// ee33 = ee30/ee1;
// ee34 = ee33 - ee32;
// ee35 = ee24/ee18;
// ee36 = ee16/ee19;
// ee37 = 2 * ee11;
// ee38 = ee36 - ee35;
// ee39 = ee28 * ee12;
// ee40 = 3/xi;
// ee49 = y * ee10 * ee7/ee1 - ee14 * ee15/ee29;
// ee53 = y * ee11 * ee25/ee1 - 2 * (ee19 * ee15/ee29);
// ee54 = 2 * ee10;
// ee55 = R_pow(ee4, 1 + ee40);
// ee56 = R_pow(ee4, ee37);
// ee58 = ee17 * ee15/xi;
// ee59 = R_pow(ee6, 2);
// ee60 = xi * ee19;
// ee61 = ee58 + 1;
// ee62 = 4/xi;
// ee63 = R_pow(ee4, ee54);
// ee65 = R_pow(ee38/ee1, 2) * R_pow(ee1, 2);
// ee66 = ee24 * ee14;
// ee67 = R_pow(ee4, ee5 - 2);
// ee68 = R_pow(ee4, ee9 - ee37);
// ee69 = ee16 * ee15;
// ee70 = 1 + ee62;
// ee71 = (ee66/ee26 - ee69/ee60) * ee6;
// ee74 = (ee16/ee55 + ee39/R_pow(ee4, ee70)) * ee34/xi + ee16 * ee53/ee56;
// ee75 = ee61/ee14;
// ee78 = ee24 * ee6 * ee49/ee26;
// ee80 = ee11 * ee68 * ee16;
// ee81 = ee4 * ee1;
// ee82 = ee16 * ee6;
// ee83 = ee15/xi;
// ee84 = ee71 - ee75;
// ee85 = ee74 - ee78;
// ee92 = ee35 + y * (ee16/ee63 + ee39/R_pow(ee4, 2 + ee40) + xi * (ee80 - ee21 * ee10 * ee13 * ee7 * ee6/ee26))/ee1 - ee36;
// ee93 = ee23 * ee15;
// ee94 = R_pow(ee34, 2);
// ee97 = y * ee67 * ee22/ee1;
// ee98 = y/ee81;
// ee100 = ee82 * ee15;
// ee101 = R_pow(ee8, ee6 - 3);
// ee103 = 2 * (y * R_pow(ee4, ee9 - 1)/ee1);
// ee104 = 4 * ee11;
// ee105 = ee98 - 2 * ee83;
// ee109 = ee11 * (ee103 - 2 * (ee25 * ee15/xi));
// ee113 = ee63 * ee1;
// ee116 = ee39 * ee15;
// ee117 = ee28 * ee53;
// ee118 = ee101 * ee27;
// ee122 = (ee97 - ee93/ee29)/ee25 - 2 * (ee34/ee60);
// ee124 = (y * (ee97 - (ee23 + ee93/xi)/xi)/ee1 - (ee7 * ee105 + ee15 * ee34/xi)/xi)/ee25 - 2 * (ee94/(xi * R_pow(ee4, ee40)));
// ee126 = 2 * (y/(R_pow(ee4, ee5 + 2) * ee1)) - (ee23 + ee2 * ee67 * ee22/ee1)/ee25;
// ee127 = xi * ee14;
// ee128 = xi * ee55;
// ee129 = xi * ee25;
// ee130 = xi * R_pow(ee4, ee62);
// ee131 = ee2 * ee11;
// 
// out(j, 0) = -(ee92/ee38);
// out(j, 1) = ee85/ee38;
// out(j, 2) = -(ee84/ee38);
// out(j, 3) =  - (((y * (((ee16 * ee126 + y * ee28 * ee12/ee113)/ee7 +
//    ee2 * ee10 * ee16 * ee59/(ee26 * ee1) - (2 * (ee16/R_pow(ee4, ee54 -
//    1)) + 2 * (ee39/ee55)))/ee4 + ((ee28 * ee126 +
//    y * ee101 * ee27/ee113)/ee19 + ee131 * R_pow(ee4, ee5 -
//    (1 + ee37)) * ee28/ee1) * ee12 + xi * (((ee7 + 2 * ee7 - y * (2 * (xi * ee10 * ee55 * ee59/ee26) -
//    ee23)/ee1) * ee21 * ee13 +
//    y * ee16 * ee6/ee81) * ee10 * ee6/ee26 + ((ee30 * ee28 * ee12/ee1 -
//    (ee25 + ee103) * ee16)/ee56 + 2 * (ee131 * R_pow(ee4, 1 +
//    6/xi - ee104) * ee16/ee1) - 2 * (ee68 * ee16)) * ee11)) -
//    ee38 * ee1)/ee1 + 2 * ee38)/ee38 - R_pow(ee92, 2)/ee65);
// out(j, 4) =  - (ee85 * ee92/ee65 + (ee74 + y * (((ee122 * ee28 -
//    ee118 * ee34/ee128)/ee19 - ee117/R_pow(ee4, ee10 + ee37)) * ee12 +
//    ((ee109 + ee25) * ee16 - ee11 * ee28 * ee12 * ee34)/ee56 +
//    (ee122 * ee16 - ee39 * ee34/ee128)/ee14 - (((ee10 * ee16 * ee34/ee7 +
//    ee16 * ee49/ee14) * ee6 + (ee10 * (ee33 -
//    (ee32 + 2 * (ee60 * ee59 * ee49/ee26))) + ee7) * ee21 * ee13) * ee6/ee26 +
//    2 * (xi * ee11 * R_pow(ee4, ee70 - ee104) * ee16 * ee53)))/ee1 -
//    ee78)/ee38);
// out(j, 5) =  - ((ee75 + y * (((ee28/ee14 - ee116/ee127)/ee19 -
//    ((ee82 + xi * ((ee7 - 2 * (R_pow(ee4, ee5 + ee54) * ee59/ee26)) * ee21 * ee13 +
//    ee61 * ee7 * ee6) * ee10)/ee26 + ee80 * ee15)) * ee6 +
//    (ee16/ee14 - ee100/ee127)/ee14)/ee1 - ee71)/ee38 -
//    ee84 * ee92/ee65);
// out(j, 6) = ((((ee124 * ee28 - ee118 * ee94/ee130)/ee19 - ee117 * ee34/R_pow(ee4, ee37 +
//    ee9)) * ee12 + (ee124 * ee16 -
//    ee39 * ee94/ee130)/ee14 + (ee16 * (y * (ee109 - 2 * (ee25/xi))/ee1 -
//    (ee19 * (2 * ee98 - 4 * ee83) + 2 * (ee15 * ee53))/xi) -
//    ee39 * ee53 * ee34/ee25)/ee56 - ee16 * ee59 * ee49 * ee34/(ee26 * ee25))/xi -
//    (((ee24 * (y * (ee10 * ee34 - ee7/xi)/ee1 -
//    (ee14 * ee105 + ee15 * ee49)/xi) + ee82 * ee49 * ee34/ee25)/xi -
//    2 * (ee66 * ee59 * R_pow(ee49, 2)/ee26)) * ee6/ee26 +
//    2 * (R_pow(ee4, ee11 - ee104) * ee16 * R_pow(ee53, 2))))/ee38 +
//    R_pow(ee85, 2)/ee65;
// out(j, 7) = ((((ee28/ee25 - ee116/ee129) * ee34/ee19 - ee69 * ee53/ee56)/xi -
//    (ee61 * ee6 + ee21 * (1 - 2 * (ee63 * ee59/ee26)) * ee13) * ee49/ee26) * ee6 +
//    ((ee16/ee25 - ee100/ee129)/ee14 -
//    R_pow(ee4, 1 - ee5) * ee16 * ee59/ee26) * ee34/xi)/ee38 -
//    ee84 * ee85/ee65;
// out(j, 8) =  - (((((ee14 - 2 * (R_pow(ee4, ee10 + ee54) * ee59/ee26)) * ee21 * ee13 +
//    2 * (ee61 * ee14 * ee6))/ee26 + (ee100/xi -
//    ee16) * ee15/ee60) * ee6 - (((ee17 - 2 * ee17) * ee6/ee13 +
//    ee17)/ee13 + ((ee58 + ee20) * ee6 - ee17) * (1/ee13 -
//    ee83) + 1)/ee14)/ee38 - R_pow(ee84, 2)/ee65);

ee1 = exp(lpsi);
ee2 = xi * y;
ee3 = ee2/ee1;
ee4 = 1 + ee3;
ee5 = exp(ldelta);
ee6 = 1/xi;
ee7 = ee5/xi;
ee8 = R_pow(ee4, ee7);
ee9 = 1 + ee6;
ee10 = 1/ee8;
ee11 = log1p(ee3);
ee12 = 1 - ee10;
ee13 = ee5 - 1;
ee14 = R_pow(ee4, ee9);
ee15 = ee13/xi;
ee16 = ee4 * ee1;
ee17 = R_pow(ee4, ee6);
ee18 = R_pow(ee4, ee15);
ee21 = ee11/(xi * ee17) - y/(ee14 * ee1);
ee22 = xi * ee12;
ee23 = y * ee9;
ee25 = R_pow(ee4, ee6 + 2) * ee1;
ee26 = 1 + ee5;
ee27 = 1 + ee7;
ee28 = R_pow(xi, 2);
ee30 = (2 * ee5 - 1)/xi;
ee31 = 1/ee18;
ee32 = 1/ee14;
ee33 = 1/R_pow(ee4, ee27);
ee34 = 2 * ee7;
ee36 = xi * ee9;
ee37 = ee2/ee16;
ee38 = ee23/ee25;
ee40 = y/ee16 - 2 * (ee11/xi);

out(j, 0) = 1 + y * (ee5/(ee12 * ee8) - ee36)/ee16;
out(j, 1) = (ee5 * ee21/(ee22 * R_pow(ee4, ee15 - 1)) + ee23/ee1)/ee4 -
   ee11/ee28;
out(j, 2) = 1 - (1/ee26 + ee11/(ee22 * ee8)) * ee5;
out(j, 3) =  - (y * (((ee32 - ee2 * ee9/ee25)/ee18 - y * (ee13/R_pow(ee4, 2 +
   ee7) + ee5/(ee12 * R_pow(ee4, 2 * ee27)))/ee1) * ee5/ee12 +
   ee36 * (ee37 - 1)/ee4)/ee1);
out(j, 4) =  - (y * (((1 - ee37) * ee9 - ee6)/ee4 - ((ee13/R_pow(ee4, ee15 +
   1) + ee5/(ee12 * R_pow(ee4, ee30 + 1))) * ee21/xi +
   (ee11/(ee28 * ee14) - ee38)/ee18) * ee5/ee12)/ee1);
out(j, 5) =  - (y * ((1/(ee12 * R_pow(ee4, 1 + ee34)) + ee33) * ee5 * ee11/xi -
   ee33) * ee5/(ee12 * ee1));
out(j, 6) = ((((ee40/ee17 + ee11 * ee21/xi)/xi + y * ((ee32 -
   ee11/(xi * ee14))/xi + ee38)/ee1)/ee18 + (ee13/R_pow(ee4, (ee5 -
   2)/xi) + ee5/(ee12 * R_pow(ee4, 2 * ee15))) * R_pow(ee21, 2)/xi) * ee5/ee12 -
   ee40/xi)/xi - y * (1/ee28 + ee23/ee16)/ee16;
out(j, 7) = (ee31 - (1/(ee12 * R_pow(ee4, ee30)) + ee31) * ee5 * ee11/xi) * ee5 * ee21/ee22;
out(j, 8) = (((1/(ee12 * R_pow(ee4, ee34)) + ee10) * ee5 * ee11/xi -
   ee10) * ee11/ee22 - (1 - ee5/ee26)/ee26) * ee5;
    
}

return out;

}

// //' @rdname egpd3d0
// [[Rcpp::export]]
arma::mat egpd3d34(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec, const arma::uvec dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec ldeltavec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 25);

if (dcate == 1) {
  lpsivec = lpsivec.elem(dupid);
  xivec = xivec.elem(dupid);
  ldeltavec = ldeltavec.elem(dupid);
}

double y, lpsi, xi, ldelta;

// double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
// double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
// double ee20, ee21, ee22, ee24, ee25, ee26, ee27, ee28, ee29;
// double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee37, ee38, ee39;
// double ee40, ee41;
// double ee50, ee54, ee55, ee58, ee59;
// double ee60, ee61, ee62, ee63, ee64, ee66, ee67, ee68, ee69;
// double ee71, ee72, ee73, ee74, ee75, ee76, ee77, ee78, ee79;
// double ee80, ee81, ee84, ee88, ee89;
// double ee90, ee91, ee92, ee93, ee97, ee98, ee99;
// double ee100, ee101, ee102, ee103, ee104, ee105, ee108;
// double ee113, ee114, ee115, ee116, ee119;
// double ee120, ee121, ee122, ee123, ee124, ee125, ee126, ee127, ee128, ee129;
// double ee130, ee131, ee133, ee134, ee136, ee139;
// double ee141, ee142, ee144, ee146, ee149;
// double ee150;
// double ee161, ee162, ee165, ee166, ee167, ee168, ee169;
// double ee170, ee171, ee173, ee174, ee175, ee176, ee177, ee178, ee179;
// double ee180, ee181, ee182, ee183, ee184, ee185, ee186, ee187, ee189;
// double ee190, ee191, ee193, ee194, ee195, ee197, ee198, ee199;
// double ee200, ee202, ee203, ee205, ee207, ee208, ee209;
// double ee210, ee211, ee212, ee216, ee217, ee218, ee219;
// double ee220, ee221, ee223, ee224, ee225, ee226, ee227, ee229;
// double ee231, ee232, ee233, ee234, ee235, ee238, ee239;
// double ee240, ee241, ee243, ee244, ee245, ee246, ee248, ee249;
// double ee250, ee251, ee254, ee256, ee257, ee258, ee259;
// double ee260, ee261, ee265, ee266, ee267, ee268;
// double ee271, ee272, ee273, ee274, ee275, ee276, ee277, ee278, ee279;
// double ee280, ee289;
// double ee290, ee291, ee293, ee294, ee295, ee296, ee297, ee298;
// double ee300, ee301, ee302, ee303, ee304, ee305, ee308, ee309;
// double ee310, ee313, ee315, ee316, ee317, ee318;
// double ee324, ee325;
// double ee330, ee331, ee334, ee335, ee336, ee337;
// double ee340, ee341, ee342, ee343, ee344, ee345, ee346, ee347, ee348;
// double ee350, ee351, ee353, ee354, ee355, ee357, ee358, ee359;
// double ee360, ee364;
// double ee372, ee373, ee374, ee375, ee377;
// double ee381, ee389;
// double ee391, ee393, ee394, ee395, ee396, ee397;
// double ee402, ee403, ee404, ee407, ee408;
// double ee411, ee412, ee413, ee417, ee418, ee419;
// double ee420, ee422, ee423, ee424, ee429;
// double ee431, ee432, ee433, ee434, ee437, ee438, ee439;
// double ee440, ee441, ee442, ee443, ee444;

double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17;
double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
double ee30, ee31, ee32, ee33, ee34, ee35, ee37, ee38, ee39;
double ee40, ee41, ee42, ee43, ee44, ee46, ee47, ee48, ee49;
double ee50, ee51, ee52, ee53, ee54, ee55, ee57, ee58, ee59;
double ee60, ee61, ee62, ee63, ee64, ee65, ee66, ee67, ee68, ee69;
double ee70, ee71, ee72, ee73, ee74, ee75, ee76, ee77, ee78, ee79;
double ee80, ee84, ee85, ee86, ee87, ee89;
double ee90, ee91, ee92, ee93, ee95, ee98, ee99;
double ee100, ee101, ee103, ee104, ee105, ee106, ee108, ee109;
double ee110, ee111, ee112, ee113, ee114, ee115, ee116, ee117, ee118, ee119;
double ee120, ee122, ee123, ee124, ee125, ee126, ee127, ee128, ee129;
double ee130, ee132, ee134, ee135, ee136, ee138, ee139;
double ee142, ee144, ee146, ee148, ee149;
double ee152, ee153, ee154, ee155, ee156, ee157, ee158, ee159;
double ee160, ee161, ee163, ee164, ee167, ee168, ee169;
double ee170, ee171, ee172, ee173, ee174, ee176, ee178;
double ee180, ee181, ee183, ee184, ee185, ee186, ee187, ee188, ee189;
double ee191, ee192, ee193, ee194, ee195, ee196, ee197, ee198, ee199;
double ee200, ee201, ee202, ee203, ee204, ee208, ee209;
double ee210, ee211, ee212, ee213, ee214, ee215, ee217, ee219;
double ee220, ee221, ee222, ee224, ee226, ee228;
double ee230, ee231, ee233, ee237, ee239;
double ee240, ee241, ee242, ee243, ee244, ee245, ee246, ee247, ee248, ee249;
double ee250, ee251, ee252, ee253, ee256, ee257, ee258, ee259;
double ee260, ee263, ee264, ee265, ee266, ee267, ee268;
double ee271, ee273, ee274, ee275, ee276;
double ee283, ee285, ee286, ee287, ee288, ee289;
double ee290, ee291, ee293, ee294, ee296, ee297, ee298, ee299;
double ee301, ee302, ee303, ee304, ee305, ee306, ee307, ee308, ee309;
double ee311, ee312, ee314, ee315;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
ldelta = ldeltavec[j];

// ee1 = exp(lpsi);
// ee2 = xi * y;
// ee3 = ee2/ee1;
// ee4 = 1 + ee3;
// ee5 = 1/xi;
// ee6 = exp(ldelta);
// ee7 = R_pow(ee4, ee5);
// ee8 = 2/xi;
// ee9 = 1/ee7;
// ee10 = log1p(ee3);
// ee11 = 1 + ee5;
// ee12 = ee5 - 1;
// ee13 = 1 + ee8;
// ee14 = ee6 - 1;
// ee15 = R_pow(ee4, ee11);
// ee16 = R_pow(ee4, ee12);
// ee17 = 1 + ee6;
// ee18 = R_pow(ee9, ee14);
// ee19 = R_pow(ee4, ee8);
// ee20 = R_pow(ee9, ee6);
// ee21 = y * ee16;
// ee22 = ee21/ee1;
// ee24 = ee7 * ee10/xi;
// ee25 = ee22 - ee24;
// ee26 = R_pow(ee4, ee13);
// ee27 = ee6 - 2;
// ee28 = ee15 * ee6;
// ee29 = R_pow(xi, 2);
// ee30 = R_pow(ee9, ee27);
// ee31 = ee20/ee17;
// ee32 = 1 - ee31;
// ee33 = R_pow(ee28, 2);
// ee34 = ee5 - 2;
// ee35 = ee32 * ee17;
// ee36 = R_pow(ee4, ee34);
// ee37 = 3/xi;
// ee38 = ee4 * ee1;
// ee39 = 2 * ee13;
// ee40 = ee10/xi;
// ee41 = y/ee38;
// ee50 = y * ee11 * ee7/ee1 - ee15 * ee10/ee29;
// ee54 = y * ee13 * ee19/ee1 - 2 * (ee26 * ee10/ee29);
// ee55 = ee30 * ee14;
// ee58 = y * ee36 * ee12/ee1;
// ee59 = ee16 * ee10;
// ee60 = 2 * ee11;
// ee61 = R_pow(ee6, 2);
// ee62 = ee35/ee28;
// ee63 = ee18/ee26;
// ee64 = 4/xi;
// ee66 = ee8 - 1;
// ee67 = ee41 - 2 * ee40;
// ee68 = R_pow(ee4, ee66);
// ee69 = R_pow(ee25, 2);
// ee71 = ee20 * ee10/xi;
// ee72 = R_pow(ee4, ee39);
// ee73 = ee6 - 3;
// ee74 = R_pow(ee4, 1 + ee37);
// ee75 = ee63 - ee62;
// ee76 = xi * ee26;
// ee77 = R_pow(ee4, ee60);
// ee78 = R_pow(ee9, ee73);
// ee79 = 2 * (y * ee68/ee1);
// ee80 = ee71 + 1;
// ee81 = ee59/xi;
// ee84 = ee18 * ee6;
// ee88 = y * (ee58 - (ee16 + ee81)/xi)/ee1 - (ee7 * ee67 +  ee10 * ee25/xi)/xi;
// ee89 = ee59/ee29;
// ee90 = R_pow(ee4, ee64);
// ee91 = ee58 - ee89;
// ee92 = ee16 + ee2 * ee36 * ee12/ee1;
// ee93 = R_pow(ee4, ee37);
// ee97 = ee79 - 2 * (ee19 * ee10/xi);
// ee98 = 4 * ee13;
// ee99 = xi * ee93;
// ee100 = ee13 * ee97;
// ee101 = R_pow(ee4, ee8 - ee39);
// ee102 = ee84 * ee10;
// ee103 = xi * ee19;
// ee104 = ee18 * ee10;
// ee105 = ee78 * ee27;
// ee108 = ee88/ee19 - 2 * (ee69/ee99);
// ee113 = 2 * (y/(R_pow(ee4, ee5 + 2) * ee1)) - ee92/ee19;
// ee114 = ee35 * ee15;
// ee115 = 1 + ee64;
// ee116 = xi * ee90;
// ee119 = ee91/ee19 - 2 * (ee25/ee76);
// ee120 = ee11 * ee25;
// ee121 = 2 * ee41;
// ee122 = R_pow(ee75/ee1, 2);
// ee123 = ee55 * ee10;
// ee124 = ee77 * ee1;
// ee125 = ee18/ee15;
// ee126 = xi * ee15;
// ee127 = xi * ee74;
// ee128 = R_pow(ee4, ee115);
// ee129 = ee122 * R_pow(ee1, 2);
// ee130 = (ee114/ee33 - ee104/ee76) * ee6;
// ee131 = ee80/ee15;
// ee133 = ee13 * ee101 * ee18;
// ee134 = ee121 - 4 * ee40;
// ee136 = (ee18/ee74 + ee55/ee128) * ee25/xi;
// ee139 = ee35 * ee6 * ee50/ee33;
// ee141 = ee18 * ee54/ee72;
// ee142 = ee18/ee19;
// ee144 = ee136 + ee141;
// ee146 = R_pow(ee4, 2 + ee37);
// ee149 = y * (ee120 - ee7/xi)/ee1 - (ee15 * ee67 + ee10 *  ee50)/xi;
// ee150 = ee130 - ee131;
// ee161 = ee18/ee77 + ee55/ee146 + xi * (ee133 - ee32 * ee11 *  ee17 * ee7 * ee6/ee33);
// ee162 = 2 * ee20;
// ee165 = y * (ee100 - 2 * (ee19/xi))/ee1 - (ee26 * ee134 +  2 * (ee10 * ee54))/xi;
// ee166 = ee144 - ee139;
// ee167 = ee100 + ee19;
// ee168 = ee19 + ee79;
// ee169 = ee102/xi;
// ee170 = y * ee161;
// ee171 = ee108 * ee30;
// ee173 = ee62 + ee170/ee1 - ee63;
// ee174 = ee30 * ee113;
// ee175 = ee119 * ee30;
// ee176 = ee171 - ee105 * ee69/ee116;
// ee177 = ee2 * ee13;
// ee178 = ee174 + y * ee78 * ee27/ee124;
// ee179 = ee30/ee19;
// ee180 = R_pow(ee50, 2);
// ee181 = 1/ee17;
// ee182 = ee168 * ee18;
// ee183 = ee80 * ee6;
// ee184 = ee175 - ee105 * ee25/ee127;
// ee185 = ee108 * ee18;
// ee186 = R_pow(ee4, ee13 - ee98);
// ee187 = ee169 - ee18;
// ee189 = ee30/ee15;
// ee190 = ee179 - ee123/ee103;
// ee191 = xi * ee11;
// ee193 = ((ee20 - ee162) * ee6/ee17 + ee20)/ee17 + ((ee71 +  ee31) * ee6 - ee20) * (ee181 - ee40);
// ee194 = ee167 * ee18;
// ee195 = ee80 * ee7;
// ee197 = ee35 * ee149 + ee84 * ee50 * ee25/ee19;
// ee198 = R_pow(ee4, ee11 + ee60);
// ee199 = R_pow(ee4, ee5 - 3);
// ee200 = R_pow(ee4, ee60 - 1);
// ee202 = ee18 * ee113 + y * ee30 * ee14/ee124;
// ee203 = ee55 * ee25;
// ee205 = ee55 * ee69/ee116;
// ee207 = ee30 * ee54;
// ee208 = ee189 - ee123/ee126;
// ee209 = R_pow(ee54, 2);
// ee210 = 2 * ee7;
// ee211 = 6/xi;
// ee212 = xi * ee13;
// ee216 = ee21 * ee30 * ee14/ee1 - ee182;
// ee217 = ee193 + 1;
// ee218 = ee194 - ee13 * ee30 * ee14 * ee25;
// ee219 = ee119 * ee18;
// ee220 = ee185 - ee205;
// ee221 = ee120 + ee7;
// ee223 = ee198 * ee61/ee33;
// ee224 = R_pow(ee4, 1 + ee211 - ee98);
// ee225 = ee7 + ee210;
// ee226 = ee7 + ee22;
// ee227 = ee101 * ee18;
// ee229 = ee18 * ee165 - ee55 * ee54 * ee25/ee19;
// ee231 = ee102/ee103;
// ee232 = ee125 - ee102/ee126;
// ee233 = ee18/ee200;
// ee234 = ee203/ee127;
// ee235 = 2 * ee125;
// ee238 = y * ee18 * ee6/ee38;
// ee239 = ee33 * ee1;
// ee240 = ee187 * ee10;
// ee241 = ee75 * ee1;
// ee243 = ee80 * ee15 * ee6;
// ee244 = ee195 * ee6;
// ee245 = ee219 - ee234;
// ee246 = ee11 * ee18;
// ee248 = R_pow(ee4, 1 - ee5) * ee18;
// ee249 = ee186 * ee18;
// ee250 = R_pow(ee4, ee115 - ee98);
// ee251 = R_pow(ee4, ee5 - (1 + ee39));
// ee254 = ee7 * ee88;
// ee256 = ee18 * ee61 * ee50;
// ee257 = ee142 - ee231;
// ee258 = 2 * ee223;
// ee259 = 2 * ee142;
// ee260 = 2 * ee16;
// ee261 = ee2 * ee11;
// ee265 = y * ee199 * ee34/ee1 - ee36 * ee10/ee29;
// ee266 = ((ee176/ee26 - ee207 * ee25/R_pow(ee4, ee39 + ee8)) * ee14 + ee220/ee15 + ee229/ee72 - ee256 * ee25/(ee33 * ee19))/xi;
// ee267 = (((ee15 - ee258) * ee32 * ee17 + 2 * ee243)/ee33 +  ee240/ee76) * ee6;
// ee268 = ee217/ee15;
// ee271 = (ee197/xi - 2 * (ee114 * ee61 * ee180/ee33)) * ee6/ee33 +  2 * (ee249 * ee209);
// ee272 = ee190 * ee25;
// ee273 = ee183 + ee35;
// ee274 = ee248 * ee61;
// ee275 = R_pow(ee4, ee5 + ee60);
// ee276 = ee77 * ee61;
// ee277 = ee55/ee74;
// ee278 = 2 * ee227;
// ee279 = ((ee272/ee26 - ee104 * ee54/ee72)/xi - (ee183 +  ee32 * (1 - 2 * (ee276/ee33)) * ee17) * ee50/ee33) *  ee6;
// ee280 = ee218/ee72;
// ee289 = (ee257/ee15 - ee274/ee33) * ee25/xi;
// ee290 = ee245/ee15;
// ee291 = ee184/ee26;
// ee293 = (ee7 - 2 * (ee275 * ee61/ee33)) * ee32 * ee17;
// ee294 = ee232/ee15;
// ee295 = ee208/ee26;
// ee296 = ee35 * ee7;
// ee297 = R_pow(ee4, ee11 + ee39);
// ee298 = ee16 * ee67;
// ee300 = ee254 + ee69/xi;
// ee301 = R_pow(ee4, ee8 - 2);
// ee302 = R_pow(ee4, ee98);
// ee303 = ee30 * ee6;
// ee304 = R_pow(ee9, ee6 - 4);
// ee305 = ee20 * ee6;
// ee308 = (y * ((ee202/ee7 + ee261 * ee18 * ee61/ee239 - (2 *  ee233 + 2 * ee277))/ee4 + (ee178/ee26 + ee177 * ee251 *  ee30/ee1) * ee14 + xi * (((ee225 - y * (2 * (ee191 *  ee74 * ee61/ee33) - ee16)/ee1) * ee32 * ee17 + ee238) *  ee11 * ee6/ee33 + (ee216/ee72 + 2 * (ee177 * ee224 *  ee18/ee1) - ee278) * ee13)) - ee241)/ee1;
// ee309 = 2 * ee18;
// ee310 = 6 * ee13;
// ee313 = ee76 * ee61 * ee50/ee33;
// ee315 = y * (ee12 * ee265 - ee36/ee29)/ee1;
// ee316 = ee266 - ee271;
// ee317 = ee267 - ee268;
// ee318 = ee279 + ee289;
// ee324 = ee150 * ee173 * ee75;
// ee325 = ee122 * ee1;
// ee330 = ee221 * ee32 * ee17 + ee246 * ee6 * ee25/ee7;
// ee331 = ee167 * ee30;
// ee334 = ee226 * ee32 * ee17 + ee238;
// ee335 = ee244 + ee296;
// ee336 = ee119 * ee78;
// ee337 = ee108 * ee78;
// ee340 = ee18 - ee169;
// ee341 = ee18 * ee149;
// ee342 = ee18 + ee309;
// ee343 = ee18/ee7;
// ee344 = ee30 - ee123/xi;
// ee345 = ee30 * ee165;
// ee346 = ee78 * ee113;
// ee347 = ee105 * ee10;
// ee348 = ee304 * ee73;
// ee350 = ee305 * ee10/xi;
// ee351 = ee20 + ee162;
// ee353 = (ee121 - 6 * ee40)/xi + y * (ee8 + ee41)/ee38;
// ee354 = ee308 + 2 * ee75;
// ee355 = 1 + ee60;
// ee357 = 2 * ee71;
// ee358 = 4 * ee11;
// ee359 = y * ((ee291 - ee207/ee297) * ee14 + ee280 + ee290 -  (((ee246 * ee25/ee7 + ee18 * ee50/ee15) * ee6 + (ee11 *  (ee22 - (ee24 + 2 * ee313)) + ee7) * ee32 * ee17) *  ee6/ee33 + 2 * (ee212 * ee250 * ee18 * ee54)));
// ee360 = y * ((ee295 - ((ee84 + xi * (ee293 + ee244) * ee11)/ee33 +  ee133 * ee10)) * ee6 + ee294);
// ee364 = y * ee301 * ee66/ee1 - 2 * (ee68 * ee10/ee29);
// ee372 = ee166 * ee173 * ee75;
// ee373 = ee176 * ee54;
// ee374 = ee197 * ee15;
// ee375 = ee221 * ee18;
// ee377 = ee331/ee15 + xi * ee184 * ee13 * ee19;
// ee381 = (ee16 + ee2 * (ee36 + 2 * ee36 + ee2 * ee199 * ee34/ee1) * ee12/ee1)/ee19 + y * (8 * (y/(R_pow(ee4, ee5 + 3) * ee1)) - ((2 * ((ee92 * ee7 + y * R_pow(ee4, 2 * ee12)/ee1) * ee16) + 2 * (ee92 * ee68))/ee90 + 2 * (ee92/ee26)))/ee1;
// ee389 = ee168 * ee30/ee15;
// ee391 = ee232 * ee7 + ee18/ee4;
// ee393 = ee144 + ee359/ee1 - ee139;
// ee394 = ee178 * ee19;
// ee395 = (ee303 * ee10/xi - ee30) * ee14;
// ee396 = ee190 * ee6;
// ee397 = (ee78/ee19 - ee347/ee103) * ee14;
// ee402 = ee337 - ee348 * ee69/ee116;
// ee403 = ee11 * ee88;
// ee404 = ee13 * (y * (2 * ee364 - 2 * (ee68/xi))/ee1 - (ee19 *  ee134 + 2 * (ee97 * ee10/xi))/xi);
// ee407 = ee15 * ee149/xi + ee180;
// ee408 = ee15 + 2 * ee15;
// ee411 = ee26 * ee165/xi + ee209;
// ee412 = ee224 * ee18;
// ee413 = ee251 * ee30;
// ee417 = ee341 * ee25/ee103;
// ee418 = ee104/ee126;
// ee419 = ee104/ee103;
// ee420 = ee125 + ee235;
// ee422 = ee345 * ee25/ee103;
// ee423 = ee30/ee74;
// ee424 = ee20 - ee350;
// ee429 = (ee315 - (ee298 + ee10 * ee91)/ee29)/ee19 - ((2 *  (ee300 * ee16) + 2 * (ee7 * ee25 * ee91))/ee90 + (2 *  (ee91/ee93) - 8 * (ee25/ee127)) * ee25)/xi;
// ee431 = (y * (ee315 - ((ee298 + ee10 * (ee58 - (ee81 + ee260)/xi) -  ee260)/xi + 2 * ee58)/xi)/ee1 - ((2 * (ee25 * ee67) +  ee10 * ee88)/xi - ee353 * ee7)/xi)/ee19 - ((2 * ee300 +  2 * ee254)/ee90 + 2 * (ee88/ee93) - 8 * (ee69/ee116)) *  ee25/xi;
// ee432 = ee181 + ee40;
// ee433 = 2 * (ee150 * ee166 * ee75/ee129);
// ee434 = 2 * (ee324/ee325);
// ee437 = 2 * (ee92 * ee25/ee99) + y * ((2 * ((ee16 * ee25/xi + ee7 * ee91) * ee16) + 2 * (ee68 * ee91))/ee90 - 8 * (ee25/(xi * R_pow(ee4, 2 + ee8))))/ee1 - (y * ((ee36 + xi * ee265) * ee12 - ee36)/ee1 - ee89)/ee19;
// ee438 = 2 * ee195;
// ee439 = 2 * ee31;
// ee440 = 2 * ee68;
// ee441 = ee162 + 4 * ee20;
// ee442 = 2 + ee64;
// ee443 = 4 * ee74;
// ee444 = 5/xi;
// 
// out(j, 0) =  - (((3 * ee241 + 4 * (ee241 + ee170) + y * (((ee381 * ee18 +
//    y * (ee178/ee15 + 2 * (ee174/ee15)) * ee14/ee1)/ee7 -
//    (ee202 * (ee9 + 2/ee7) + (ee423 + 2 * ee423) * ee14 +
//    ee233 + ee233 + ee233 + ee2 * ee342 * ee11 * ee61/ee239))/ee4 +
//    ((ee381 * ee30 + y * ((ee346 + y * ee304 * ee73/ee124)/ee15 +
//    2 * (ee346/ee15)) * ee27/ee1)/ee26 + ee2 * ((ee394 -
//    ee389)/ee72 + ee178 * ee101 + 2 * (ee177 * R_pow(ee4, ee444 -
//    ee98) * ee30/ee1) - (ee413 + 2 * ee413)) * ee13/ee1 - ee178 * (1/ee26 +
//    2/ee26)) * ee14 + 2 * ee161 + xi * ((((ee19 +
//    y * (ee440 + 2 * (ee2 * ee301 * ee66/ee1) + 4 * ee68)/ee1) * ee18 +
//    y * (ee394 - 2 * ee389) * ee14/ee1)/ee72 + ee2 * ((2 * (ee128 * ee216) -
//    2 * ((ee168 * ee26 + ee177 * ee90/ee1) * ee19 * ee18))/ee302 +
//    2 * (ee250 * ee216) + 8 * (ee177 * R_pow(ee4, 10/xi +
//    2 - ee310) * ee18/ee1) - (2 * ee412 + 4 * ee412)) * ee13/ee1 -
//    (ee227 + (1/ee72 + 2/ee72) * ee216 + ee278)) * ee13 +
//    ((ee225 + y * (xi * (ee11 * (2 * ee74 + ee443) * ee61/ee33 -
//    ee58) - (ee260 + 4 * ee16))/ee1 - 4 * ee7) * ee32 * ee17 +
//    y * (ee225 * ee202 + ee191 * (2 * (y * ee19 * ee18 * ee6/ee1) +
//    4 * (ee334 * ee26) - ee35 * (8 * (ee261 * R_pow(ee4, 2 +
//    ee444) * ee61/ee239) - 2 * ((ee226 * ee15 +
//    ee261 * ee19/ee1) * ee7))) * ee6/ee33 - (ee226 * ee420 + ee342/ee4)) * ee6/ee1) * ee11 * ee6/ee33)))/ee1 -
//    8 * ee75)/ee75 -
//    (ee75 * (2 - 2 * (R_pow(ee173, 2)/ee129)) + ee308 + 2 * ee354) * ee173/ee129);
// out(j, 1) =  - (((ee166 * ee1 + y * (((((ee226 * ee18/ee19 +
//    ee343 + ee343) * ee25 - xi * (2 * (ee334 * ee15 * ee50) + y * (ee35 * (2 * ((ee221 * ee15 +
//    ee191 * ee7 * ee50) * ee7) -
//    8 * (ee191 * R_pow(ee4, ee442) * ee61 * ee50/ee33)) + 2 * (ee330 * ee26))/ee1) * ee6/ee33) * ee11 +
//    (ee235 - ee202) * ee50 +
//    y * (2 * (ee375/ee15) + ee191 * (2 * (ee245 * ee7) -
//    2 * (ee7 * ee18 * ee61 * ee50/ee33)))/ee1) * ee6 + (ee11 * (y * (ee16 +
//    ee260)/ee1 - (ee225 * ee10/xi + 4 * ee313)) + ee7 +
//    ee7 + ee7 + y * (ee16 + ee191 * ee91)/ee1) * ee32 * ee17) * ee6/ee33 +
//    ((ee168 * ee13 * ee30 * ee25/ee19 + y * ee377/ee1) * ee14 -
//    (ee167 + y * (ee440 + 2 * (ee212 * ee364))/ee1) * ee18)/ee72 +
//    ((ee30 * ee437 + ee27 * (y * ((ee336 - ee348 * ee25/ee127)/ee15 +
//    ee336/ee15)/ee1 - ee346 * ee25/ee103))/ee26 +
//    (2 * (ee30/ee297) - ee178/ee72) * ee54 + y * (ee377/ee72 -
//    2 * (ee212 * R_pow(ee4, ee37 - ee98) * ee30 * ee54))/ee1 -
//    2 * ee291) * ee14 + (ee18 * ee437 + ee14 * (y * (ee184/ee15 +
//    ee175/ee15)/ee1 - ee174 * ee25/ee103))/ee15 + xi * ((4 * (ee250 * ee18) -
//    2 * (ee186 * ee216)) * ee54 + y * ((2 * ((ee167 * ee26 +
//    ee212 * ee19 * ee54) * ee19 * ee18) + 2 * (ee218 * ee128))/ee302 -
//    8 * (ee212 * R_pow(ee4, 2 + 8/xi -
//    ee310) * ee18 * ee54))/ee1) * ee13 - (2 * ee280 + 2 * ee290)))/ee1 -
//    2 * ee166)/ee75 + (ee166 * ee354 - ee173 * (2 * (ee372/ee129) +
//    2 * ee393))/ee129);
// out(j, 2) =  - (((y * ((((ee182 * ee10 + ee2 * ee208 * ee19/ee1)/ee72 +
//    2 * (ee227 * ee10) + ee2 * (ee208 * ee101 - 2 * (ee13 * ee224 * ee18 * ee10))/ee1) * ee13 +
//    (ee344 * ee113 +
//    y * ((ee78/ee15 - ee347/ee126) * ee14/ee15 + ee105/ee77)/ee1)/ee26 +
//    ((ee309 + ee2 * (ee391 - 2 * (ee200 * ee18 * ee61/ee33)) * ee11/ee1 -
//    ee202 * ee15) * ee6 + xi * (ee273 * ee226 +
//    (ee438 + y * ee391/ee1 - (2 * (ee334 * ee77) + ee2 * (ee35 * (ee443 -
//    8 * (R_pow(ee4, ee355 + ee37) * ee61/ee33)) +
//    2 * (ee335 * ee26)) * ee11/ee1) * ee6/ee33) * ee6 + 2 * ee293) * ee11)/ee33 -
//    2 * ee295) * ee6 + (ee340 * ee113 + y * (ee208 * ee6/ee15 +
//    ee55/ee77)/ee1)/ee15 - 2 * ee294) - ee150 * ee1)/ee1 +
//    2 * ee150)/ee75 - (ee150 * ee354 + ee173 * (2 * (ee131 +
//    ee360/ee1 - ee130) - 2 * (ee324/ee129)))/ee129);
// out(j, 3) =  - ((ee316 * ee173 + ee166 * ((2 * (ee372/ee325) +
//    2 * ee359)/ee1 + 2 * ee136 + 2 * ee141 - 2 * ee139))/ee129 +
//    (ee266 + y * (((ee429 * ee30 - (ee402/ee15 + 2 * (ee336 * ee25/ee19)) * ee27/xi)/ee26 +
//    (2 * (R_pow(ee4, ee5 - ee98) * ee30 * ee54) -
//    ee184/ee72) * ee54 - (ee184 * ee54 + ee345/ee126)/ee72) * ee14 +
//    ((ee404 + 2 * (ee97/xi)) * ee18 - (ee176 * ee13 * ee19 +
//    2 * (ee331 * ee25/ee103)) * ee14)/ee72 +
//    (ee429 * ee18 - (ee176/ee15 + 2 * (ee175 * ee25/ee19)) * ee14/xi)/ee15 +
//    (8 * (ee212 * R_pow(ee4, 2 + ee211 - ee310) * ee18 * ee54) -
//    2 * (ee218 * ee186)) * ee54 - (((ee220 * ee11 * ee7 +
//    (ee341/ee15 + 2 * (ee375 * ee25/ee19))/xi + (2 * ee219 -
//    (2 * (ee256/ee33) + 2 * ee234)) * ee50 - (4 * (ee330 * ee15 * ee50) +
//    xi * ee32 * ee11 * ee17 * (2 * (ee407 * ee7) -
//    8 * (ee146 * ee61 * ee180/ee33))) * ee6/ee33) * ee6 + (ee403 +
//    2 * (ee25/xi)) * ee32 * ee17) * ee6/ee33 + (2 * (ee218 * ee26 * ee54) +
//    2 * (xi * ee411 * ee13 * ee19 * ee18))/ee302))/ee1 -
//    ee271)/ee75);
// out(j, 4) =  - ((ee318 * ee173 + ee166 * (ee131 + (ee360 - ee434)/ee1 -
//    ee130) - ee150 * ee393)/ee129 + (ee318 + y * (((ee344 * ee119 -
//    (ee397/ee15 + ee105/ee74) * ee25/xi)/ee26 +
//    (2 * (ee13 * ee250 * ee18 * ee10) - ee208/ee72) * ee54 - (((((ee257 * ee7 +
//    ee343) * ee25 - xi * (ee35 * (4 * ee26 - 8 * (R_pow(ee4, ee355 +
//    ee8) * ee61/ee33)) + 2 * (ee335 * ee15)) * ee6 * ee50/ee33) * ee11 +
//    ee245 * ee15 + (ee235 - (ee418 +
//    2 * (R_pow(ee4, ee60 - ee11) * ee18 * ee6/ee33)) * ee6) * ee50 -
//    2 * (ee330 * ee77 * ee6/ee33)) * ee6 + ee273 * ee221)/ee33 +
//    (ee194 * ee10/xi + ee190 * ee13 * ee19 * ee25)/ee72)) * ee6 +
//    (ee340 * ee119 - (ee396/ee15 + ee277) * ee25/xi)/ee15)/ee1)/ee75);
// out(j, 5) =  - ((ee268 + y * (((((ee432 * (ee235 - ee125) -
//    (ee418 + ee18/(ee17 * ee15))) * ee6 + ee125 + ee235 - ee420)/ee17 +
//    (ee187/ee15 - ee235) * ee10/xi) * ee6 + ee125)/ee15 +
//    (((ee395/ee15 - 2 * (ee303/ee15)) * ee10/xi + ee189)/ee26 +
//    ee187 * ee13 * ee101 * ee10 - ((ee18 + 2 * (ee232 * ee15) -
//    2 * (ee77 * ee18 * ee61/ee33)) * ee6 + xi * ((ee217 * ee7 +
//    ee438 - (ee35 * (4 * ee275 - 8 * (R_pow(ee4, ee5 + ee358) * ee61/ee33)) +
//    4 * (ee335 * ee77)) * ee6/ee33) * ee6 + ee296) * ee11)/ee33) * ee6)/ee1 -
//    ee267)/ee75 - (ee317 * ee173 +
//    ee150 * ((2 * ee360 - ee434)/ee1 + 2 * ee131 - 2 * ee130))/ee129);
// out(j, 6) = ((((ee431 * ee30 - (ee402/ee19 + 2 * (ee337/ee19)) * ee27 * ee25/xi)/ee26 +
//    (2 * (R_pow(ee4, 1 - ee98) * ee30 * ee54 * ee25) -
//    ee176/ee72) * ee54 - (ee373 + ee422)/ee72) * ee14 +
//    (ee431 * ee18 - (ee176/ee19 + 2 * (ee171/ee19)) * ee14 * ee25/xi)/ee15 +
//    (ee18 * (y * (ee404 - (4 * ee97 - 4 * ee19)/ee29)/ee1 -
//    (2 * (ee134 * ee54) + 2 * (ee10 * ee165/xi) -
//    ((4 * ee41 - 12 * ee40)/xi + y * (ee121 + ee64)/ee38) * ee26)/xi) -
//    (ee373 + 2 * ee422) * ee14)/ee72 - ((ee417 + (2 * ee185 -
//    (2 * (ee274 * ee50/ee33) + 2 * (ee203/ee116)) * ee25) * ee50) * ee61/ee33 +
//    2 * (ee229 * ee186 * ee54)))/xi -
//    (((((ee185 - (ee205 + 2 * (ee374 * ee6/ee33))) * ee50 + 2 * ee417) * ee6 +
//    ee35 * (y * (ee403 - (2 * ee22 - (2 * ee24 +
//    ee210))/ee29)/ee1 - (2 * (ee50 * ee67) + ee10 * ee149/xi -
//    ee353 * ee15)/xi))/xi - (ee35 * (2 * ee407 - 8 * (ee276 * ee180/ee33)) +
//    2 * (ee374/xi)) * ee61 * ee50/ee33) * ee6/ee33 +
//    ((2 * (ee411 * ee18) + 2 * (ee229 * ee26/xi))/ee302 - 8 * (R_pow(ee4, ee442 -
//    ee310) * ee18 * ee209)) * ee54))/ee75 +
//    (ee266 + 2 * ee316 + 2 * (R_pow(ee166, 2) * ee75/ee129) - ee271) * ee166/ee129;
// out(j, 7) = ((((ee344 * ee108 - (ee397/ee19 + ee105/ee90) * ee69/xi)/ee26 +
//    (2 * (ee249 * ee10 * ee54) - ee272/ee72) * ee54 -
//    (ee190 * ee54 * ee25 + ee104 * ee165/xi)/ee72)/xi - ((ee273 * ee149 +
//    ((ee259 - ee231) * ee50 * ee25 - 2 * (ee197 * ee77 * ee6/ee33)) * ee6)/xi -
//    (ee35 * (4 * ee15 - 8 * ee223) +
//    2 * (ee273 * ee15)) * ee61 * ee180/ee33)/ee33) * ee6 +
//    ((ee340 * ee108 - (ee396/ee19 + ee55/ee90) * ee69/xi)/ee15 -
//    (ee220 * ee15 + (ee259 - (ee419 + 2 * (R_pow(ee4, ee60 - ee8) * ee18 * ee6/ee33)) * ee6) * ee50 * ee25) * ee61/ee33)/xi)/ee75 +
//    (ee166 * (2 * ee318 - ee433) - ee316 * ee150)/ee129;
// out(j, 8) = ((((((ee432 * (ee259 - ee142) - (ee419 + ee18/(ee17 * ee19))) * ee6 +
//    ee142 + ee259 - (ee142 + ee259))/ee17 +
//    (ee187/ee19 - ee259) * ee10/xi) * ee6 + ee142)/ee15 - (ee248 +
//    2 * (ee257 * ee15) - 2 * (R_pow(ee4, ee355 - ee5) * ee18 * ee61/ee33)) * ee61/ee33) * ee25/xi +
//    ((((ee395/ee19 - 2 * (ee303/ee19)) * ee10/xi +
//    ee179) * ee25/ee26 + ee240 * ee54/ee72)/xi -
//    ((ee193 + ee357 + 3 - (ee35 * (4 * ee77 - 8 * (R_pow(ee4, ee358) * ee61/ee33)) +
//    4 * (ee273 * ee77)) * ee6/ee33) * ee6 +
//    ee35) * ee50/ee33) * ee6)/ee75 - (ee317 * ee166 +
//    ee150 * (2 * ee279 + 2 * ee289 - ee433))/ee129;
// out(j, 9) =  - (((((ee217 * ee408 + (ee408 - ee258) * ee80 -
//    (ee35 * (4 * ee198 - 8 * (R_pow(ee4, ee11 + ee358) * ee61/ee33)) +
//    4 * ((ee243 + ee114) * ee77)) * ee6/ee33) * ee6 + ee114)/ee33 +
//    ((ee342 - ee169) * ee6 * ee10/xi - ee18) * ee10/ee76) * ee6 -
//    ((((ee441 - ((ee351/ee17 + ee357) * ee6 + ee20 +
//    ee162)) * ee10/xi + (ee441 - ((2 * (ee424 * ee17) + 2 * ((1 +
//    2 * ee6) * ee20) - 8 * ee305)/ee17 + 2 * ee424 + 3 * ((ee71 +
//    ee439) * ee6)))/ee17) * ee6 + ee20 - ee20)/ee17 + 1 -
//    (((ee351 - ee350) * ee10/xi - ((ee357 + ee439) * ee6 - ee351)/ee17) * ee6 -
//    ee20) * ee10/xi)/ee15)/ee75 - (ee267 + 2 * ee317 -
//    (ee268 + 2 * (R_pow(ee150, 2) * ee75/ee129))) * ee150/ee129);

ee1 = exp(lpsi);
ee2 = xi * y;
ee3 = ee2/ee1;
ee4 = 1 + ee3;
ee5 = 1/xi;
ee6 = exp(ldelta);
ee7 = 1 + ee5;
ee8 = log1p(ee3);
ee9 = R_pow(ee4, ee7);
ee10 = ee6 - 1;
ee11 = ee6/xi;
ee12 = ee5 + 2;
ee13 = ee10/xi;
ee14 = R_pow(ee4, ee12);
ee15 = R_pow(ee4, ee5);
ee16 = R_pow(ee4, ee11);
ee17 = R_pow(ee4, ee13);
ee20 = ee8/(xi * ee15) - y/(ee9 * ee1);
ee21 = ee14 * ee1;
ee22 = ee6 - 2;
ee23 = ee4 * ee1;
ee24 = 1/ee16;
ee25 = R_pow(xi, 2);
ee26 = ee22/xi;
ee27 = y * ee7;
ee28 = ee27/ee21;
ee29 = 1 - ee24;
ee30 = y/ee23;
ee31 = R_pow(ee4, ee26);
ee32 = ee6 * ee8;
ee33 = 1/ee9;
ee34 = ee8/xi;
ee35 = 1 + ee11;
ee37 = ee30 - 2 * ee34;
ee38 = R_pow(ee4, ee35);
ee39 = ee8/(ee25 * ee9);
ee40 = ee8/(xi * ee9);
ee41 = 1/ee17;
ee42 = ee39 - ee28;
ee43 = ee5 + 3;
ee44 = ee2 * ee7;
ee46 = (ee37/ee15 + ee8 * ee20/xi)/xi + y * ((ee33 - ee40)/xi +  ee28)/ee1;
ee47 = R_pow(ee4, ee43);
ee48 = R_pow(ee4, ee13 + 1);
ee49 = R_pow(ee20, 2);
ee50 = ee25 * ee14;
ee51 = xi * ee17;
ee52 = ee33 - ee44/ee21;
ee53 = ee32/ee51;
ee54 = ee47 * ee1;
ee55 = xi * ee31;
ee57 = y * ee12/ee54;
ee58 = ee6 - 3;
ee59 = ee8/ee50;
ee60 = ee59 - ee57;
ee61 = ee58/xi;
ee62 = R_pow(ee4, ee61);
ee63 = ee41 - ee53;
ee64 = 1/ee14;
ee65 = 1/ee31;
ee66 = 1/ee38;
ee67 = xi * ee29;
ee68 = xi * ee16;
ee69 = ee10 * ee8;
ee70 = 2 * ee11;
ee71 = ee32/ee68;
ee72 = xi * ee48;
ee73 = R_pow(ee4, 2 + ee11);
ee74 = 2 * ee6;
ee75 = ee46/ee17;
ee76 = ee37/ee9;
ee77 = 2 * ee30;
ee78 = ee7 * ee60;
ee79 = (ee74 - 1)/xi;
ee80 = 2 * ee35;
ee84 = y * (1/ee50 - ee78)/ee1;
ee85 = ee75 + ee10 * ee49/ee55;
ee86 = 2/ee9;
ee87 = 2/xi;
ee89 = ee32/(xi * ee38) - ee66;
ee90 = ee71 - ee24;
ee91 = ee52/ee17;
ee92 = ee42/ee17;
ee93 = xi * ee62;
ee95 = ee91 - y * ee10/(ee73 * ee1);
ee98 = ee10 * ee20/ee72 + ee92;
ee99 = ee69/ee55;
ee100 = ee65 - ee99;
ee101 = R_pow(ee4, ee70);
ee103 = ee22 * ee49/ee93;
ee104 = 2/ee17;
ee105 = ee87 + ee30;
ee106 = R_pow(ee4, ee26 + 1);
ee108 = (ee77 - 6 * ee34)/xi + y * ee105/ee23;
ee109 = 2 * ee13;
ee110 = 2/ee14;
ee111 = ee8 * ee42;
ee112 = ee29 * ee1;
ee113 = R_pow(ee4, ee109);
ee114 = R_pow(ee4, ee80);
ee115 = ee100 * ee6;
ee116 = ee10/ee31;
ee117 = ee76 + ee111;
ee118 = ee117/ee25;
ee119 = ee41 + ee104;
ee120 = ee118 + ee84;
ee122 = ee7 * (ee64 + xi * ee60) + ee64;
ee123 = R_pow(ee4, ee79);
ee124 = ee79 + 1;
ee125 = 1/ee48;
ee126 = ee64 + ee110;
ee127 = 2/ee16;
ee128 = ee53 - ee41;
ee129 = ee2 * ee12;
ee130 = ee2/ee23;
ee132 = ((ee46 * ee8 + 2 * (ee20 * ee37))/xi - ee108/ee15)/xi +  y * (((((ee86 - ee40)/xi + ee28) * ee8 - (ee76 + ee86))/xi -  2 * ee28)/xi - ee84)/ee1;
ee134 = ee46/ee31 + ee103;
ee135 = R_pow(ee4, ee124);
ee136 = R_pow(ee4, ee13 + 2);
ee138 = ee39 - y * ee122/ee1;
ee139 = xi * ee106;
ee142 = ee44 * (ee126 - ee129/ee54)/ee1 - ee33;
ee144 = ee115 + ee116;
ee146 = ee52/ee31 - y * ee22/(ee136 * ee1);
ee148 = ee69/ee72 - ee125;
ee149 = ee10/ee48;
ee152 = ee22 * ee20/ee139 + ee42/ee31;
ee153 = 1 + ee6;
ee154 = ee24 + ee127;
ee155 = ee148 * ee6;
ee156 = ee29 * ee90;
ee157 = R_pow(ee4, 1 + ee70);
ee158 = ee10/ee73;
ee159 = ee65 + 2/ee31;
ee160 = 2 * ee130;
ee161 = ee25 * ee47;
ee163 = ee156 - ee32/(xi * ee101);
ee164 = ee52 * ee20;
ee167 = (ee24 + 2/(ee29 * ee101)) * ee6 * ee8/xi;
ee168 = 2/ee38;
ee169 = ee6 * ee49;
ee170 = ee132/ee17;
ee171 = (ee46 * ee159 + ee103) * ee10;
ee172 = ee46 * ee63;
ee173 = ee144 * ee49;
ee174 = (ee164/ee55 - y * (ee152/ee9 + ee42/ee48)/ee1) *  ee10;
ee176 = ee146/ee9 + 2 * (ee52/ee48);
ee178 = ee29 * ee63;
ee180 = R_pow(ee4, ee5 + 4) * ee1;
ee181 = ee63 * ee52;
ee183 = ee20 * ee42;
ee184 = ee138/ee17;
ee185 = ee142/ee17;
ee186 = ee41 + ee6 * (ee53 - ee119) * ee8/xi;
ee187 = 1/ee62;
ee188 = ee77 + 6/xi;
ee189 = ee6 * ee20;
ee191 = y * (ee158 - ee155/ee9)/ee1;
ee192 = ee170 + ee171 * ee20/xi;
ee193 = (ee134/ee9 + 2 * (ee183/ee31)) * ee10;
ee194 = ee85 * ee29;
ee195 = ee85/ee38;
ee196 = ee120/ee17;
ee197 = ee172 + ee173/xi;
ee198 = (ee115/ee9 + ee149) * ee20;
ee199 = ee174 + ee184;
ee200 = ee98 * ee20;
ee201 = ee29 * ee89;
ee202 = ee181 - ee191;
ee203 = ee63/ee38;
ee204 = ee167 - ee24;
ee208 = (ee154 - ee71) * ee6 * ee8/xi - ee24;
ee209 = ee22 * ee8;
ee210 = ee128 * ee42;
ee211 = ee128/ee9;
ee212 = ee89/ee17;
ee213 = ee185 + y * ee176 * ee10/ee1;
ee214 = ee77 + ee87;
ee215 = 4 * ee6;
ee217 = ee8/ee161 - y * ee43/ee180;
ee219 = ee193/xi + ee196;
ee220 = ee194 - ee169/(xi * ee113);
ee221 = ee195 + ee200/ee17;
ee222 = ee75 + (ee116 + 2 * (ee6/(ee29 * ee113))) * ee49/xi;
ee224 = ee98 * ee29;
ee226 = (ee37/ee14 + ee8 * ee60)/ee25 + y * (1/ee161 - ee12 *  ee217)/ee1;
ee228 = ee178 + ee32/(xi * ee123);
ee230 = ee63/ee16;
ee231 = (ee187 - ee209/ee93) * ee10;
ee233 = ee188/xi + y * ee214/ee23;
ee237 = (ee168 - ee211) * ee6 * ee8/xi - ee66;
ee239 = ee58 * ee49/(xi * R_pow(ee4, (ee6 - 4)/xi));
ee240 = ee210 - ee198/xi;
ee241 = ee212 - ee203;
ee242 = ee89/ee16;
ee243 = ee90/ee17;
ee244 = ee90/ee38;
ee245 = 2 - ee160;
ee246 = 2 * ee57;
ee247 = 2/(ee29 * ee123);
ee248 = 3 * ee6;
ee249 = 4 * ee4;
ee250 = ee32/ee55;
ee251 = y * ee6;
ee252 = ee132/ee31;
ee253 = ee197/ee17;
ee256 = (ee46 * (ee187 + 2/ee62) + ee239) * ee22 * ee20/xi;
ee257 = ee85 * ee8;
ee258 = ee134 * ee10;
ee259 = ee199/ee38;
ee260 = ee120 * ee8;
ee263 = ee46 * ee100 + (ee231 + ee22/ee62) * ee49/xi;
ee264 = ee202/ee38;
ee265 = ee95 * ee29;
ee266 = ee95 * ee20;
ee267 = ee95 * ee8;
ee268 = ee208/ee16;
ee271 = ee108/ee9;
ee273 = ((24 * ee34 - 6 * ee30)/xi - y * ee188/ee23)/xi -  y * ee233/ee23;
ee274 = ee98 * ee8;
ee275 = ee241 * ee29;
ee276 = ee213/ee38;
ee283 = (ee119 - ee53) * ee6 * ee8/xi - ee41;
ee285 = ee167 + 2 * ee90 - ee24;
ee286 = (ee248 - 1)/xi;
ee287 = (ee215 - 1)/xi;
ee288 = ee10 * (ee250 - ee65);
ee289 = ee242 + ee244;
ee290 = ee243 - ee230;
ee291 = 1 + ee80;
ee293 = 1 + ee2 * (ee160 - 3)/ee23;
ee294 = ee41 - (ee41 + ee247) * ee6 * ee8/xi;
ee296 = 2 * (1 + 2 * ee3) + ee249;
ee297 = 2 * ee10;
ee298 = 2 * (ee6/ee31);
ee299 = 2/ee47;
ee301 = 3 * ee11;
ee302 = 4 * ee11;
ee303 = 4/ee17;
ee304 = 4/ee16;
ee305 = 6 * ee3;
ee306 = 8 * ee3;
ee307 = xi * ee7;
ee308 = xi * ee135;
ee309 = xi * ee157;
ee311 = ee2 * (3 - ee160)/ee23;
ee312 = ee130 - 1;
ee314 = y * (((2 * (ee8/(xi * ee14)) - ee110)/xi - ee246)/ee25 -  ee226 * ee7)/ee1;
ee315 = ee27/ee23;

out(j, 0) =  - (y * ((ee185 + y * ((ee95 * (ee66 + ee168) -
   2 * (ee251/(ee29 * R_pow(ee4, ee291 + ee11) * ee1))) * ee6/ee29 +
   ee176 * ee10)/ee1) * ee6/ee29 + ee307 * ee293/ee4)/ee1);
out(j, 1) =  - (y * (((ee266/ee51 - y * (2 * (ee98/ee38) + 2 * (ee189/(ee67 * R_pow(ee4, ee13 +
   ee80))))/ee1) * ee6/ee29 +
   ee174 + ee184) * ee6/ee29 + (ee7 * (ee311 - 1) - ee312/xi)/ee4)/ee1);
out(j, 2) =  - (y * (ee181 + ee6 * (y * (2 * (ee89/ee38) + 2 * (ee32/(ee67 * R_pow(ee4, ee80 +
   ee11))))/ee1 - ee267/ee68)/ee29 -
   ee191) * ee6/ee112);
out(j, 3) = y * ((((ee195 + (2 * (ee98/ee17) + 2 * (ee189/(ee67 * R_pow(ee4, (ee248 -
   2)/xi + 1)))) * ee20) * ee6/ee29 +
   ee193)/xi + ee196) * ee6/ee29 + 2 * (ee7 * ee245 - ee87)/R_pow(ee4, 2))/ee1;
out(j, 4) =  - (y * (((ee274/ee16 + (ee212 + 2 * (ee32/(ee67 * R_pow(ee4, ee286 +
   1))) - ee203) * ee20) * ee6/ee29 - ee198)/xi +
   ee210) * ee6/ee112);
out(j, 5) =  - (y * ((ee168 - ((ee244 + 2 * ee242 + 2 * (ee32/(ee67 * R_pow(ee4, 1 +
   ee301))))/ee29 + ee211)) * ee6 * ee8/xi -
   ee66) * ee6/ee112);
out(j, 6) = ((ee170 + ((ee85 * ee119 + 2 * (ee169/(ee67 * R_pow(ee4, 3 * ee13)))) * ee6/ee29 +
   ee171) * ee20/xi) * ee6/ee29 +
   ee108/xi)/xi + y * (ee105/ee25 + y * (1/ee25 + 2 * ee315)/ee23)/ee23;
out(j, 7) = ((ee173 + ((2 * (ee63/ee17) - 2 * (ee32/(ee67 * R_pow(ee4, (ee297 +
   ee6)/xi)))) * ee49 - ee257/ee16) * ee6/ee29)/xi +
   ee172) * ee6/ee67;
out(j, 8) = (((ee243 + 2 * (ee32/(ee67 * R_pow(ee4, ee286))) -
   2 * ee230)/ee29 + ee53 - ee119) * ee6 * ee8/xi + ee41) * ee6 * ee20/ee67;
out(j, 9) = (((ee154 - ((ee154 * ee90 + 2 * (ee32/(ee67 * R_pow(ee4, ee301))))/ee29 +
   ee71)) * ee6 * ee8/xi - ee24) * ee8/ee67 -
   (1 - (3 - 2 * (ee6/ee153)) * ee6/ee153)/ee153) * ee6;
out(j, 10) =  - (y * (((ee33 - ee44 * (ee126 + 4/ee14 - ee129 * (ee299 +
   4/ee47 - ee2 * ee43/ee180)/ee1)/ee1)/ee17 + y * (((ee142/ee31 +
   y * ((ee52/ee62 - y * ee58/(R_pow(ee4, ee26 +
   2) * ee1))/ee9 + 2 * (ee52/ee106)) * ee22/ee1)/ee9 + (ee125 +
   2/ee48) * ee142 - 3 * (ee146 * ee52)) * ee10 + (2 * ee276 +
   y * ((2 * ((ee265 + ee251/(ee114 * ee1))/ee114) + 4 * (ee265/ee114) -
   8 * (ee251/(R_pow(ee4, 2 + ee80 + ee70) * ee1)))/ee29 +
   4 * (ee95/ee114)) * ee6/ee112 - ((ee91 - y * (ee158 +
   2 * (ee6/(ee29 * ee114)))/ee1) * ee95 + 2 * (R_pow(ee95, 2) -
   ee276))) * ee6/ee29)/ee1) * ee6/ee29 + ee307 * (ee2 * (7 +
   ee2 * ((ee306 - ee296)/ee4 - 6)/ee23)/ee23 - 1)/ee4)/ee1);
out(j, 11) =  - (y * (((ee213 * ee20/ee51 + y * (ee259 + ((ee149 +
   2 * (ee6/(ee29 * ee135))) * ee20/xi + ee92) * ee95 +
   (4 * (ee266/ee308) + y * (2 * ((ee189/ee308 - ee224)/ee114) -
   (4 * (ee224/ee114) + 8 * (ee189/(xi * R_pow(ee4, ee124 + ee80)))))/ee112) * ee6/ee29 +
   2 * (ee259 + ee95 * ee98))/ee1) * ee6/ee29 +
   (ee20 * ee142/ee55 + y * (((ee164/ee93 - y * ((ee58 * ee20/(xi * R_pow(ee4, ee61 +
   1)) + ee42/ee62)/ee9 +
   ee42/ee106)/ee1) * ee22 + ee138/ee31)/ee9 + ee146 * ee42 + 2 * (ee152 * ee52) +
   2 * (ee138/ee48))/ee1) * ee10 + (y * (ee122 +
   ee110 + xi * (2 * ee78 - y * (ee7 * (ee299 + xi * ee12 * ee217) +
   ee12/ee47)/ee1))/ee1 - ee39)/ee17) * ee6/ee29 +
   (ee7 * (1 + ee2 * (ee2 * ((ee296 - ee306)/ee4 + 6)/ee23 - 7)/ee23) -
   ee293/xi)/ee4)/ee1);
out(j, 12) =  - (y * (ee63 * ee142 - ((ee213 * ee8/ee68 + y * (ee95 * ((ee66 +
   2/(ee29 * ee157)) * ee6 * ee8/xi - ee66) +
   (4 * (ee267/ee309) - y * (2 * ((ee201 - ee32/ee309)/ee114) +
   4 * (ee201/ee114) + 8 * (ee32/(xi * R_pow(ee4, ee291 + ee70))))/ee112) * ee6/ee29 +
   2 * (ee95 * ee89 - ee264) - ee264)/ee1) * ee6/ee29 +
   y * (2 * ((ee155 - ee149) * ee52) - ((ee100 * ee52 -
   y * (ee22/ee136 - (ee209/ee139 - 1/ee106) * ee10/ee9)/ee1) * ee6/ee9 +
   ee146 * ee10/ee9))/ee1)) * ee6/ee112);
out(j, 13) =  - (y * ((((ee222 * ee95 + 2 * (ee199 * ee20/ee17))/xi -
   y * ((((4 * (ee224/ee135) + 8 * (ee189/(xi * R_pow(ee4, ee109 +
   ee80)))) * ee20 + 2 * (ee220/ee114))/ee29 + 4 * (ee200/ee135)) * ee6/ee67 +
   2 * (ee219/ee38 + R_pow(ee98, 2)))/ee1) * ee6/ee29 +
   ((ee134 * ee52 + 2 * (ee20 * ee138/ee31))/xi -
   y * ((((ee46/ee62 + ee239)/ee9 + 2 * (ee183/ee62)) * ee22/xi +
   ee120/ee31)/ee9 + ee120/ee48 + 2 * (ee152 * ee42))/ee1) * ee10 +
   (ee118 + y * ((ee64 - 2 * (ee8/ee14))/ee25 +
   ee246 - ee7 * (ee59 + xi * ee226 - ee57))/ee1)/ee17) * ee6/ee29 +
   (ee27 * (4 - ee2 * ((ee249 - ee305)/ee4 + 6)/ee23)/ee23 -
   (2 * ee311 - (2 + 2 * ee312))/ee25)/ee4)/ee1);
out(j, 14) =  - (y * ((((ee202/ee17 + ee95 * ee294) * ee20 -
   ee199 * ee8/ee16) * ee6/ee29 + ee144 * ee52 * ee20)/xi + ee63 * ee138 -
   y * (((((2 * (ee228/ee114) - (4 * (ee201/ee135) +
   8 * (ee32/(xi * R_pow(ee4, ee79 + ee80))))) * ee20/ee29 -
   4 * (ee274/ee157)) * ee6/ee67 - 2 * (ee98 * ee89 + ee240/ee38))/ee29 -
   ((ee99 - ee65) * ee42 - (ee231/ee9 + ee22/ee106) * ee20/xi)/ee9) * ee6 +
   (ee149 - ee155) * ee42 + ee152 * ee10/ee9)/ee1) * ee6/ee112);
out(j, 15) =  - (y * (((ee95 * ee204 - 2 * (ee202/ee16)) * ee8/xi -
   y * (((2 * (ee163/ee114) + 4 * (ee201/ee157) + 8 * (ee32/(xi * R_pow(ee4, ee80 +
   ee70))))/ee29 + 4 * (ee89/ee157)) * ee6 * ee8/ee67 +
   2 * (R_pow(ee89, 2) - ee237/ee38))/ee1) * ee6/ee29 +
   ee186 * ee52 - y * (ee158 - (((2 * (ee6/ee48) -
   ee288/ee9) * ee8/xi - ee125)/ee9 + 2 * (ee148/ee9)) * ee6)/ee1) * ee6/ee112);
out(j, 16) =  - (y * ((((6 * (1 - ee130) - 6)/xi + 2 * (y * ee245/ee23))/ee25 +
   y * (ee245/ee25 + y * ((2 * ee4 - ee305)/ee4 +
   4) * ee7/ee23)/ee23)/ee4 - (((ee192/ee38 + (ee219 * ee119 +
   ((2 * (ee220/ee135) + 2 * (ee221 * ee29/ee17) + 8 * (ee169/(xi * R_pow(ee4, (ee215 -
   3)/xi + 1))))/ee29 + 2 * (ee221/ee17)) * ee6/ee67) * ee20 +
   (ee222 + 2 * ee85) * ee98) * ee6/ee29 +
   ((ee252 + ee256)/ee9 + ee120 * ee159 * ee20 + 3 * (ee134 * ee42)) * ee10)/xi +
   ((ee260 + 2 * (ee42 * ee37) -
   ee271)/ee25 + ee314)/ee17) * ee6/ee29)/ee1);
out(j, 17) =  - (y * ((((ee219 * ee8/ee16 + ee85 * ee89 + (((ee247 +
   ee104) * ee6 * ee8/xi - ee104) * ee98 + 2 * (ee240/ee17)) * ee20 +
   ((2 * (ee275/ee17) + 8 * (ee32/(xi * R_pow(ee4, (ee215 -
   2)/xi + 1))) - 2 * (ee228/ee135)) * ee49/ee29 +
   2 * (ee221 * ee8/ee16)) * ee6/ee67 - ee197/ee38)/ee29 - ee263/ee9) * ee6 -
   (ee258/ee9 + 2 * (ee144 * ee20 * ee42)))/xi +
   ee120 * ee128) * ee6/ee112);
out(j, 18) =  - (y * ((((ee237/ee17 + 2 * (ee63 * ee89) - (((2 * (ee163/ee135) +
   2 * (ee275/ee16) + 8 * (ee32/(xi * R_pow(ee4, ee287 +
   1))))/ee29 + 2 * (ee241/ee16)) * ee6 * ee8/ee67 +
   ee186/ee38)) * ee20 - (ee204 * ee98 + 2 * (ee240/ee16)) * ee8) * ee6/ee29 -
   ((((ee288 - ee298) * ee8/xi + ee65)/ee9 +
   2 * (ee100/ee9)) * ee6 + ee149) * ee20)/xi + ee283 * ee42) * ee6/ee112);
out(j, 19) =  - (y * (((ee285 * ee89 + ((2 * (ee163/ee157) +
   2 * (ee289 * ee29/ee16) + 8 * (ee32/(xi * R_pow(ee4, 1 + ee302))))/ee29 +
   2 * (ee289/ee16)) * ee6 * ee8/ee67 - (ee208/ee38 +
   ee237 * ee154))/ee29 + ee66 + ee168 - (ee283/ee9 + (ee33 +
   ee86) * ee128)) * ee6 * ee8/xi - ee66) * ee6/ee112);
out(j, 20) = (((((ee132 * ee8 + 3 * (ee46 * ee37) - 3 * (ee108 * ee20))/xi -
   ee273/ee15)/xi + y * (((ee271 + (ee76 + 2 * ee117 +
   6/ee9 + ee111)/xi - (ee260 + (2 * ee37 + 6) * ee42))/xi +
   3 * ee84)/xi - ee314)/ee1)/ee17 + (((ee256 + 4 * ee252) * ee20 +
   3 * (ee134 * ee46)) * ee10 + (ee222 * ee85 + (((2 * (ee220/ee113) +
   4 * (ee194/ee113) + 8 * (ee169/(xi * R_pow(ee4, 4 * ee13))))/ee29 +
   4 * (ee85/ee113)) * ee6 * ee20/ee67 +
   2 * (ee192/ee17)) * ee20 + 2 * (ee192 * ee20/ee17 + R_pow(ee85, 2))) * ee6/ee29)/xi) * ee6/ee29 +
   ee273/xi)/xi - y * (ee233/ee25 +
   y * (ee214/ee25 + y * (2/ee25 + 6 * ee315)/ee23)/ee23)/ee23;
out(j, 21) = ((((ee253 + ee85 * ee294 + ((2 * (ee228/ee113) +
   4 * (ee178/ee113) - 8 * (ee32/(xi * R_pow(ee4, (3 * ee10 +
   ee6)/xi)))) * ee49/ee29 - 4 * (ee257/ee123)) * ee6/ee67 + 2 * (ee253 +
   ee85 * ee63)) * ee20 - ee192 * ee8/ee16) * ee6/ee29 +
   (ee263 * ee6 + ee258 + 2 * (ee46 * ee144)) * ee20)/xi +
   ee132 * ee63) * ee6/ee67;
out(j, 22) = ((((ee85 * ee204 - 2 * (ee197/ee16)) * ee8 + (((2 * (ee163/ee113) +
   8 * (ee32/(xi * R_pow(ee4, (ee297 + ee74)/xi))) -
   4 * (ee178/ee123))/ee29 - 4 * (ee63/ee123)) * ee6 * ee8/ee67 +
   2 * (R_pow(ee63, 2) + ee186/ee17)) * ee49) * ee6/ee29 +
   (((ee10 * (ee250 - ee159) - ee298) * ee8/xi + ee65 +
   ee65 + ee65) * ee6 + ee116) * ee49)/xi + ee46 * ee186) * ee6/ee67;
out(j, 23) = (((ee285 * ee63 + ee208/ee17 - (((2 * (ee163/ee123) +
   2 * (ee290 * ee29/ee16) + 8 * (ee32/(xi * R_pow(ee4, ee287))))/ee29 +
   2 * (ee290/ee16)) * ee6 * ee8/ee67 + ee186 * ee154))/ee29 +
   (ee104 + ee303 - ee53) * ee6 * ee8/xi - (ee119 +
   ee303)) * ee6 * ee8/xi + ee41) * ee6 * ee20/ee67;
out(j, 24) = ((((ee204 * ee90 + ((2 * (ee163/ee101) + 4 * (ee156/ee101) +
   8 * (ee32/(xi * R_pow(ee4, ee302))))/ee29 + 4 * (ee90/ee101)) * ee6 * ee8/ee67 +
   2 * (R_pow(ee90, 2) - ee268) -
   2 * ee268)/ee29 + ee24 + ee127 + ee304 - (ee127 + ee304 -
   ee71) * ee6 * ee8/xi) * ee6 * ee8/xi - ee24) * ee8/ee67 -
   (1 - (7 - ((2 * (1 + ee74) + 4 * ee153 - 8 * ee6)/ee153 +
   6) * ee6/ee153) * ee6/ee153)/ee153) * ee6;
    
}

return out;

}

// //' Extended generalized Pareto distribution of type 4 (eGPD4) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each eGPD parameter
// //' @param X1 a design matrix for the eGPD log scale parameter
// //' @param X2 a design matrix for the eGPD shape parameter
// //' @param X3 a design matrix for the eGPD log delta
// //' @param X4 a design matrix for the eGPD log kappa
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return egpd4d0 a scalar, the negative log-liklihood
// //' @return egpd4d12 a matrix, first then second derivatives w.r.t. eGPD4 parameters
// //' @return egpd4d34 a matrix, third then fourth derivatives w.r.t. eGPD4 parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double egpd4d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, const arma::mat& X4, arma::vec yvec, const arma::uvec& dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec ldeltavec = X3 * Rcpp::as<arma::vec>(pars[2]);
arma::vec lkappavec = X4 * Rcpp::as<arma::vec>(pars[3]);
int nobs = yvec.size();

if (dcate == 1) {
  lpsivec = lpsivec.elem(dupid);
  xivec = xivec.elem(dupid);
  ldeltavec = ldeltavec.elem(dupid);
  lkappavec = lkappavec.elem(dupid);
}

double y, lpsi, xi, ldelta, lkappa;
double ee1, ee4, ee5, ee6, ee7, ee8;
double nllh=0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
ldelta = ldeltavec[j];
lkappa = lkappavec[j];
/*
ee1 = exp(ldelta);
ee2 = exp(lpsi);
ee3 = 1 + xi * y/ee2;
ee4 = 1/xi;
ee5 = R_pow(ee3, ee4);
ee7 = 1/ee5;
ee8 = (1 - R_pow(ee7, ee1)/(1 + ee1)) * (1 + ee1);
ee9 = exp(lkappa);

nllh -= log(-((R_pow(ee7, ee1 - 1)/R_pow(ee3, 1 + 2/xi) - ee8/(R_pow(ee3, 1 + ee4) * 
ee1)) * R_pow(1 - ee8/(ee5 * ee1), ee9/2 - 1) * ee9/(2 * ee2)));
 */   

// ee1 = exp(ldelta);
// ee4 = xi * y/exp(lpsi);
// ee5 = 1 + ee4;
// ee6 = R_pow(ee5, ee1/xi);
// ee7 = 1 + ee1;
// ee8 = 1/xi;
// 
// nllh += (1 - 0.5 * exp(lkappa)) * log(1 - (1 - 1/(ee7 * ee6)) * ee7/(R_pow(ee5, ee8) * 
// ee1)) + (1 + ee8) * log1p(ee4) + 0.693147180559945 + 
// ldelta + lpsi - (lkappa + log(1 - 1/ee6) + log1p(ee1));

ee1 = exp(ldelta);
ee4 = xi * y/exp(lpsi);

if (ee4 <= -1.0) {
    nllh = 1e20;
    break;
}

ee5 = 1 + ee4;
ee6 = R_pow(ee5, ee1/xi);
ee7 = 1 + ee1;
ee8 = 1/xi;

nllh += (1 - 0.5 * exp(lkappa)) * log(1 - (1 - 1/(ee7 * ee6)) * ee7/(R_pow(ee5, ee8) * 
ee1)) + (1 + ee8) * log1p(ee4) + 0.693147180559945 + 
ldelta + lpsi - (lkappa + log(1 - 1/ee6) + log1p(ee1));

}

return(nllh);

}

// //' @rdname egpd4d0
// [[Rcpp::export]]
arma::mat egpd4d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, arma::vec yvec, const arma::uvec dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec ldeltavec = X3 * Rcpp::as<arma::vec>(pars[2]);
arma::vec lkappavec = X4 * Rcpp::as<arma::vec>(pars[3]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 14);

if (dcate == 1) {
  lpsivec = lpsivec.elem(dupid);
  xivec = xivec.elem(dupid);
  ldeltavec = ldeltavec.elem(dupid);
  lkappavec = lkappavec.elem(dupid);
}

double y, lpsi, xi, ldelta, lkappa;

// double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
// double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
// double ee20, ee21, ee22, ee23, ee25, ee26, ee27, ee28, ee29;
// double ee30, ee31, ee32, ee34, ee35, ee36, ee38, ee39;
// double ee40, ee41, ee42, ee43, ee44, ee45, ee46, ee47, ee48;
// double ee50, ee51, ee52, ee53, ee54, ee55, ee56;
// double ee61, ee63, ee67;
// double ee71, ee72, ee73, ee74, ee75, ee79;
// double ee80, ee81, ee82, ee84, ee87, ee88;
// double ee90, ee91, ee97;
// double ee100, ee103, ee104, ee106;
// double ee110, ee111, ee112, ee113, ee114, ee115, ee116, ee117, ee118, ee119;
// double ee122, ee123, ee124, ee127;
// double ee130, ee131, ee136, ee138, ee139;
// double ee140, ee142, ee144, ee145, ee146, ee147, ee148;
// double ee151, ee152, ee156, ee157, ee159;
// double ee160, ee161, ee162, ee165, ee169;
// double ee173, ee177, ee179;
// double ee183, ee185, ee186, ee187, ee189;
// double ee190, ee191, ee192, ee193, ee196, ee197, ee199;
// double ee201, ee203, ee204, ee205, ee207;
// double ee210, ee212, ee214, ee215, ee216, ee217, ee218, ee219;

// double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
// double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
// double ee20, ee22, ee23, ee25, ee28, ee29;
// double ee30, ee31, ee32, ee33, ee34, ee36, ee37, ee38;
// double ee40, ee41, ee43, ee44, ee45, ee46, ee47, ee49;
// double ee50, ee51, ee54, ee56, ee58, ee59;
// double ee60, ee61, ee62, ee64, ee66, ee67, ee68, ee69;
// double ee70, ee71, ee72, ee73, ee74, ee75, ee76, ee77, ee78, ee79;
// double ee81, ee82, ee85, ee86, ee87, ee88, ee89;
// double ee90, ee91, ee92;

double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee20, ee21, ee22, ee23, ee24, ee25, ee27, ee28, ee29;
double ee30, ee31, ee33, ee34, ee35, ee36, ee38, ee39;
double ee41, ee42, ee43, ee44, ee45, ee46, ee48, ee49;
double ee50, ee53, ee55, ee57, ee59;
double ee61, ee62, ee64, ee66, ee67, ee69;
double ee71, ee72, ee73, ee74, ee75, ee77, ee78, ee79;
double ee80, ee81, ee82, ee83, ee84, ee85, ee86, ee87, ee88;
double ee90, ee91, ee92, ee93, ee94, ee95;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
ldelta = ldeltavec[j];
lkappa = lkappavec[j];

// ee1 = exp(lpsi);
// ee2 = exp(ldelta);
// ee3 = xi * y;
// ee4 = ee3/ee1;
// ee5 = 1 + ee4;
// ee6 = 1/xi;
// ee7 = R_pow(ee5, ee6);
// ee8 = 1 + ee2;
// ee9 = 1/ee7;
// ee10 = R_pow(ee9, ee2);
// ee11 = ee10/ee8;
// ee12 = 1 - ee11;
// ee13 = ee12 * ee8;
// ee14 = 2/xi;
// ee15 = 1 + ee6;
// ee16 = ee2 - 1;
// ee17 = ee7 * ee2;
// ee18 = exp(lkappa);
// ee19 = 1 + ee14;
// ee20 = R_pow(ee9, ee16);
// ee21 = R_pow(ee5, ee15);
// ee22 = log1p(ee4);
// ee23 = ee18/2;
// ee25 = 1 - ee13/ee17;
// ee26 = R_pow(ee5, ee19);
// ee27 = ee21 * ee2;
// ee28 = ee23 - 1;
// ee29 = ee6 - 1;
// ee30 = R_pow(ee5, ee29);
// ee31 = R_pow(ee25, ee28);
// ee32 = ee20/ee26;
// ee34 = ee32 - ee13/ee27;
// ee35 = R_pow(ee17, 2);
// ee36 = y * ee30;
// ee38 = ee7 * ee22/xi;
// ee39 = ee36/ee1;
// ee40 = ee39 - ee38;
// ee41 = R_pow(ee27, 2);
// ee42 = ee23 - 2;
// ee43 = 3/xi;
// ee44 = ee2 - 2;
// ee45 = R_pow(ee25, ee42);
// ee46 = R_pow(ee5, ee14);
// ee47 = R_pow(ee9, ee44);
// ee48 = R_pow(xi, 2);
// ee50 = ee10 * ee22/xi;
// ee51 = 2 * ee19;
// ee52 = ee47 * ee16;
// ee53 = ee50 + 1;
// ee54 = ee34 * ee31;
// ee55 = ee13 * ee2;
// ee56 = 2 * ee1;
// ee61 = R_pow(ee5, ee43);
// ee63 = 2 * ee15;
// ee67 = y * ee15 * ee7/ee1 - ee21 * ee22/ee48;
// ee71 = y * ee19 * ee46/ee1 - 2 * (ee26 * ee22/ee48);
// ee72 = R_pow(ee2, 2);
// ee73 = R_pow(ee5, ee63);
// ee74 = R_pow(ee1, 2);
// ee75 = xi * ee26;
// ee79 = R_pow(ee5, 1 + ee43);
// ee80 = R_pow(ee5, ee51);
// ee81 = ee32 - ee13 * ee30 * ee2/ee35;
// ee82 = 4/xi;
// ee84 = ee53/ee7 - ee13 * ee7 * ee2/ee35;
// ee87 = ee20/ee61 - ee55/ee35;
// ee88 = log(ee25);
// ee90 = ee19 * R_pow(ee5, ee14 - ee51) * ee20;
// ee91 = R_pow(ee5, ee6 - 2);
// ee97 = ee13 * ee21;
// ee100 = ee20 * ee22;
// ee103 = ee20/ee73 + ee52/R_pow(ee5, 2 + ee43) + xi * (ee90 - ee12 * ee15 * ee8 * ee7 * ee2/ee41);
// ee104 = 1 + ee82;
// ee106 = (ee97/ee41 - ee100/ee75) * ee2 - ee53/ee21;
// ee110 = (ee20/ee79 + ee52/R_pow(ee5, ee104)) * ee40/xi + ee20 * ee71/ee80 - ee55 * ee67/ee41;
// ee111 = R_pow(ee56, 2);
// ee112 = R_pow(ee54 * ee18/ee56, 2);
// ee113 = R_pow(ee18, 2);
// ee114 = ee103 * ee31;
// ee115 = ee20 * ee2;
// ee116 = 2 * ee74;
// ee117 = ee22/xi;
// ee118 = ee106 * ee31;
// ee119 = ee110 * ee31;
// ee122 = ee81 * ee34 * ee45 * ee28;
// ee123 = ee30 * ee22;
// ee124 = R_pow(ee40, 2);
// ee127 = y * ee91 * ee29/ee1;
// ee130 = ee84 * ee34 * ee45 * ee28;
// ee131 = ee122 + ee114;
// ee136 = ee34 * ee87 * ee45 * ee28 * ee40/xi;
// ee138 = ee31 * ee18 * ee88;
// ee139 = ee5 * ee1;
// ee140 = ee118 - ee130;
// ee142 = ee119 + ee136;
// ee144 = ee31 + ee138/2;
// ee145 = R_pow(ee25, ee23 - 3);
// ee146 = 4 * (ee112 * ee74);
// ee147 = y * ee131;
// ee148 = y/ee139;
// ee151 = ee123/xi;
// ee152 = ee115 * ee22;
// ee156 = ee147/ee116 - 2 * (ee54 * ee1/ee111);
// ee157 = ee148 - 2 * ee117;
// ee159 = ee7 * ee157;
// ee160 = R_pow(ee5, ee14 - 1);
// ee161 = xi * ee21;
// ee162 = xi * ee46;
// ee165 = y * (ee127 - (ee30 + ee151)/xi)/ee1;
// ee169 = ee31/2;
// ee173 = ee45 * ee28 * ee88/2 + ee45/2;
// ee177 = ee73 * ee1;
// ee179 = R_pow(ee9, ee2 - 3);
// ee183 = (ee127 - ee123/ee48)/ee46 - 2 * (ee40/ee75);
// ee185 = (ee165 - (ee159 + ee22 * ee40/xi)/xi)/ee46 - 2 *  (ee124/(xi * ee61));
// ee186 = 2 * (ee112 * ee1);
// ee187 = 2 * (y * ee160/ee1);
// ee189 = 2 * (y/(R_pow(ee5, ee6 + 2) * ee1)) - (ee30 + ee3 * ee91 * ee29/ee1)/ee46;
// ee190 = 4 * ee19;
// ee191 = xi * ee79;
// ee192 = xi * R_pow(ee5, ee82);
// ee193 = y * ee20;
// ee196 = ((ee10 - 2 * ee10) * ee2/ee8 + ee10)/ee8 + ((ee50 +  ee11) * ee2 - ee10) * (1/ee8 - ee117) + 1;
// ee197 = ee144 * ee34;
// ee199 = ee53 * ee7 * ee2;
// ee201 = ee183 * ee20 - ee52 * ee40/ee191;
// ee203 = ee185 * ee20 - ee52 * ee124/ee192;
// ee204 = ee19 * (ee187 - 2 * (ee46 * ee22/xi));
// ee205 = R_pow(ee5, 2);
// ee207 = ee20 * ee189 + y * ee47 * ee16/ee177;
// ee210 = ee20 * ee72;
// ee212 = ee20/ee21 - ee152/ee161;
// ee214 = ee20/ee46 - ee152/ee162;
// ee215 = ee52 * ee22;
// ee216 = ee47 * ee71;
// ee217 = ee179 * ee44;
// ee218 = ee3 * ee19;
// ee219 = ee193 * ee2;
// 
// out(j, 0) = -(2 * (ee1 * ee156/ee54));
// out(j, 1) = ee142/ee54;
// out(j, 2) = -(ee140/ee54);
// out(j, 3) = -(ee144/ee31);
// out(j, 4) =  - (2 * (ee1 * (y * ((((((ee207/ee91 + ee193 * ee72/(ee35 * ee1))/ee205 -
//    (ee13 * (y * (2 * (R_pow(ee5, ee43 -
//    2) * ee72/ee35) - xi * ee91 * ee29)/ee1 - ee30) - ee219/(ee205 * ee1)) * ee2/ee35) * ee45 +
//    y * R_pow(ee81, 2) * ee145 * ee42/ee1) * ee34 +
//    2 * (y * ee81 * ee103 * ee45/ee1)) * ee28 +
//    ((ee207/ee7 + ee3 * ee15 * ee20 * ee72/(ee41 * ee1))/ee5 +
//    ((ee47 * ee189 + y * ee179 * ee44/ee177)/ee26 + ee218 * R_pow(ee5, ee6 -
//    (1 + ee51)) * ee47/ee1) * ee16 + xi * (((ee36 * ee47 * ee16/ee1 -
//    (ee46 + ee187) * ee20)/ee80 + 2 * (ee218 * R_pow(ee5, 1 +
//    6/xi - ee190) * ee20/ee1)) * ee19 - (ee13 * (y * (2 * (xi * ee15 * ee79 * ee72/ee41) -
//    ee30)/ee1 -
//    ee7) - ee219/ee139) * ee15 * ee2/ee41)) * ee31)/ee116 - 2 * (ee131/ee111)) -
//    (ee34 * (2 * ee31 - 16 * (ee31 * ee74/ee111)) * ee1 +
//    2 * ee147)/ee111)/ee54) - ee113 * R_pow(ee156, 2)/ee112);
// out(j, 5) =  - (ee142 * ee113 * ee156/ee186 + 2 * ((2 * (ee142 * ee1/ee111) +
//    y * ((((ee201/ee7 - (ee13 * (ee127 - (ee151 +
//    2 * (ee160 * ee72 * ee40/ee35))/xi) + 2 * (ee115 * ee40/ee161)) * ee2/ee35) * ee45 -
//    ee81 * ee87 * ee145 * ee42 * ee40/xi) * ee34 -
//    (ee110 * ee81 * ee45 + ee103 * ee87 * ee45 * ee40/xi)) * ee28 +
//    (((ee183 * ee47 - ee217 * ee40/ee191)/ee26 -
//    ee216/R_pow(ee5, ee15 + ee51)) * ee16 + ((ee204 + ee46) * ee20 -
//    ee19 * ee47 * ee16 * ee40)/ee80 + ee201/ee21 - (((ee15 * ee20 * ee40/ee7 +
//    ee20 * ee67/ee21) * ee2 + (ee15 * (ee39 -
//    (ee38 + 2 * (ee75 * ee72 * ee67/ee41))) + ee7) * ee12 * ee8) * ee2/ee41 +
//    2 * (xi * ee19 * R_pow(ee5, ee104 - ee190) * ee20 * ee71))) * ee31)/ee116) * ee1/ee54));
// out(j, 6) =  - (2 * (ee1 * (y * (((((ee212/ee30 - ee210/ee35)/ee5 -
//    ((ee30 - 2 * (R_pow(ee5, ee43 - 1) * ee72/ee35)) * ee12 * ee8 +
//    ee53 * ee30 * ee2) * ee2/ee35) * ee45 - ee84 * ee81 * ee145 * ee42) * ee34 +
//    ee106 * ee81 * ee45 - ee84 * ee103 * ee45) * ee28 +
//    (((ee47/ee21 - ee215/ee161)/ee26 - ((ee115 +
//    xi * ((ee7 - 2 * (R_pow(ee5, ee6 + ee63) * ee72/ee41)) * ee12 * ee8 +
//    ee199) * ee15)/ee41 + ee90 * ee22)) * ee2 +
//    ee212/ee21) * ee31)/ee116 - 2 * (ee140 * ee1/ee111))/ee54) -
//    ee140 * ee113 * ee156/ee186);
// out(j, 7) = -(2 * (ee1 * (y * ((ee173 * ee81 * ee34 + ee114 * ee88/2) * ee18 +
//    ee122 + ee114)/ee116 - 2 * (ee197 * ee1/ee111))/ee54) -
//    ee197 * ee113 * ee156/ee186);
// out(j, 8) = ((((ee203/ee7 - (ee13 * (ee165 - (ee159 + (2 * (ee7 * ee72 * ee40/ee35) +
//    ee117) * ee40)/xi) + 2 * (ee115 * ee124/ee162)) * ee2/ee35) * ee45 -
//    R_pow(ee87, 2) * ee145 * ee42 * ee124/xi) * ee34 -
//    2 * (ee110 * ee87 * ee45 * ee40)) * ee28/xi +
//    ((((ee185 * ee47 - ee217 * ee124/ee192)/ee26 - ee216 * ee40/R_pow(ee5, ee51 +
//    ee14)) * ee16 + ee203/ee21 + (ee20 * (y * (ee204 -
//    2 * (ee46/xi))/ee1 - (ee26 * (2 * ee148 -
//    4 * ee117) + 2 * (ee22 * ee71))/xi) - ee52 * ee71 * ee40/ee46)/ee80 -
//    ee210 * ee67 * ee40/(ee41 * ee46))/xi - (((ee13 * (y * (ee15 * ee40 -
//    ee7/xi)/ee1 - (ee21 * ee157 + ee22 * ee67)/xi) +
//    ee115 * ee67 * ee40/ee46)/xi - 2 * (ee97 * ee72 * R_pow(ee67, 2)/ee41)) * ee2/ee41 +
//    2 * (R_pow(ee5, ee19 -
//    ee190) * ee20 * R_pow(ee71, 2)))) * ee31)/ee54 + R_pow(ee142, 2) * ee113/ee146;
// out(j, 9) = (((((ee214/ee7 - ((ee20/ee7 + ee50 + 1) * ee2 +
//    ee12 * (1 - 2 * (ee46 * ee72/ee35)) * ee8) * ee2/ee35) * ee45 -
//    ee84 * ee87 * ee145 * ee42) * ee34 + ee106 * ee87 * ee45) * ee40/xi -
//    ee110 * ee84 * ee45) * ee28 + ((((ee47/ee46 -
//    ee215/ee162) * ee40/ee26 - ee100 * ee71/ee80)/xi - (ee53 * ee2 +
//    ee12 * (1 - 2 * (ee73 * ee72/ee41)) * ee8) * ee67/ee41) * ee2 +
//    (ee214/ee21 - R_pow(ee5, 1 - ee6) * ee20 * ee72/ee41) * ee40/xi) * ee31)/ee54 -
//    ee140 * ee142 * ee113/ee146;
// out(j, 10) = ((ee119 * ee88/2 + ee173 * ee34 * ee87 * ee40/xi) * ee18 +
//    ee119 + ee136)/ee54 - ee142 * ee144 * ee34 * ee113/ee146;
// out(j, 11) =  - ((((((ee21 - 2 * (R_pow(ee5, ee15 + ee63) * ee72/ee41)) * ee12 * ee8 +
//    2 * (ee53 * ee21 * ee2))/ee41 + (ee152/xi -
//    ee20) * ee22/ee75) * ee2 - ee196/ee21) * ee31 - (((ee196/ee7 -
//    ((ee7 - 2 * (ee61 * ee72/ee35)) * ee12 * ee8 +
//    2 * ee199) * ee2/ee35) * ee45 - R_pow(ee84, 2) * ee145 * ee42) * ee34 +
//    2 * (ee106 * ee84 * ee45)) * ee28)/ee54 - R_pow(ee140, 2) * ee113/ee146);
// out(j, 12) = -(((ee118 * ee88/2 - ee84 * ee173 * ee34) * ee18 +
//    ee118 - ee130)/ee54 - ee140 * ee144 * ee34 * ee113/ee146);
// out(j, 13) =  - (((ee138/4 + ee169 + ee169 + ee169) * ee18 * ee88 +
//    ee31)/ee31 - R_pow(ee144, 2) * R_pow(ee34, 2) * ee113/ee146);

// ee1 = exp(lpsi);
// ee2 = exp(ldelta);
// ee3 = xi * y;
// ee4 = ee3/ee1;
// ee5 = 1 + ee4;
// ee6 = 1/xi;
// ee7 = 1 + ee2;
// ee8 = ee2/xi;
// ee9 = R_pow(ee5, ee8);
// ee10 = R_pow(ee5, ee6);
// ee11 = log1p(ee4);
// ee12 = ee7 * ee9;
// ee13 = 1 + ee6;
// ee14 = 1/ee12;
// ee15 = 1 - ee14;
// ee16 = R_pow(ee5, ee13);
// ee17 = ee2 - 1;
// ee18 = ee15 * ee7;
// ee19 = 1/ee9;
// ee20 = ee17/xi;
// ee22 = R_pow(ee5, ee20);
// ee23 = 1 - ee18/(ee10 * ee2);
// ee25 = 1/ee2 - 1;
// ee28 = ee11/(xi * ee10) - y/(ee16 * ee1);
// ee29 = ee25 * ee2;
// ee30 = 1 - ee19;
// ee31 = exp(lkappa);
// ee32 = xi * ee9;
// ee33 = ee29/xi;
// ee34 = R_pow(ee5, ee33);
// ee36 = R_pow(ee5, ee6 + 2) * ee1;
// ee37 = 1 + ee8;
// ee38 = ee11/ee32;
// ee40 = ee5 * ee1;
// ee41 = R_pow(ee5, ee37);
// ee43 = 1 - ee7/ee2;
// ee44 = 1 - 0.5 * ee31;
// ee45 = ee14 + ee38;
// ee46 = y * ee13;
// ee47 = ee7/xi;
// ee49 = R_pow(ee5, ee47 + 1);
// ee50 = 1/ee16;
// ee51 = ee46/ee36;
// ee54 = ee43 * ee15/ee10 + ee45/ee10;
// ee56 = ee18/(ee16 * ee2) - 1/ee49;
// ee58 = ee18/ee2 - ee19;
// ee59 = R_pow(ee28, 2);
// ee60 = 1/ee22;
// ee61 = 1/ee41;
// ee62 = ee11/xi;
// ee64 = R_pow(xi, 2);
// ee66 = y/ee40 - 2 * ee62;
// ee67 = ((ee66/ee10 + ee11 * ee28/xi)/xi + y * ((ee50 - ee11/(xi *  ee16))/xi + ee51)/ee1)/ee22;
// ee68 = R_pow(ee5, ee20 + 1);
// ee69 = R_pow(ee5, (ee2 - 2)/xi);
// ee70 = R_pow(ee5, 2 + ee8);
// ee71 = (1/ee34 - ee29/ee34) * ee11;
// ee72 = (ee50 - ee3 * ee13/ee36)/ee22;
// ee73 = (ee11/(ee64 * ee16) - ee51)/ee22;
// ee74 = 2 * ee2;
// ee75 = ee2 * ee11;
// ee76 = xi * ee22;
// ee77 = xi * ee41;
// ee78 = ee67 + ee17 * ee59/(xi * ee69);
// ee79 = ee43 * ee2;
// ee81 = ee72 - y * ee17/(ee70 * ee1);
// ee82 = (ee74 - 1)/xi;
// ee85 = ee17 * ee28/(xi * ee68) + ee73;
// ee86 = 0.5 * (ee31 * log(ee23));
// ee87 = 2 * ee8;
// ee88 = xi * ee23;
// ee89 = xi * ee30;
// ee90 = xi * ee13;
// ee91 = ee3/ee40;
// ee92 = ee46/ee40;
// 
// out(j, 0) = 1 - y * (ee56 * ee44/ee23 + (ee90 - ee2/(ee30 * ee9))/ee5)/ee1;
// out(j, 1) = ee92 - ((ee58 * ee44/ee23 - ee2/(ee30 * ee22)) * ee28 +
//    ee62)/xi;
// out(j, 2) = 1 - (ee54 * ee44/ee23 + (1/ee7 + ee11/(ee89 * ee9)) * ee2);
// out(j, 3) = -(ee86 + 1);
// out(j, 4) =  - (y * ((ee81/ee10 + y * (R_pow(ee56, 2)/ee23 -
//    2/R_pow(ee5, ee47 + 2))/ee1 - (ee81/ee34 - y * ee25 * ee2/ee36) * ee15 * ee7/ee2) * ee44/ee23 +
//    (ee72 - y * (ee17/ee70 +
//    ee2/(ee30 * R_pow(ee5, 2 * ee37)))/ee1) * ee2/ee30 + ee90 * (ee91 -
//    1)/ee5)/ee1);
// out(j, 5) =  - (y * ((ee56 * ee58 * ee28/ee88 - (ee85/R_pow(ee5, ee6 -
//    1) + 2 * (ee28/ee32) - (ee85/R_pow(ee5, ee33 - 1) +
//    ee29 * ee28/xi) * ee15 * ee7/ee2)/ee5) * ee44/ee23 + ((1 -
//    ee91) * ee13 - ee6)/ee5 - ((ee17/ee68 + ee2/(ee30 * R_pow(ee5, ee82 +
//    1))) * ee28/xi + ee73) * ee2/ee30)/ee1);
// out(j, 6) =  - (y * ((ee54 * ee56/ee23 + ((1/(ee7 * ee41) +
//    ee11/ee77) * ee2 - ee61)/ee10 + ee45/ee16 - ((((ee61 + ee75/ee77 -
//    ee61)/ee34 - ee71/ee77) * ee7/ee2 - ee43/ee16) * ee15 +
//    ee79/(ee7 * ee49))) * ee44/ee23 + ((1/(ee30 * R_pow(ee5, 1 +
//    ee87)) + ee61) * ee2 * ee11/xi - ee61) * ee2/ee30)/ee1);
// out(j, 7) = 0.5 * (y * ee56 * ee31/(ee23 * ee1));
// out(j, 8) =  - ((((ee78/ee34 + ee10 * ee25 * ee2 * ee59/xi) * ee15 * ee7/ee2 +
//    (R_pow(ee58, 2)/ee23 - 2/ee22) * ee59/xi -
//    ee78/ee10) * ee44/ee23 + ee66/xi - (ee67 + (ee17/ee69 + ee2/(ee30 * R_pow(ee5, 2 * ee20))) * ee59/xi) * ee2/ee30)/xi +
//    y * (1/ee64 + ee92)/ee40);
// out(j, 9) =  - (((((ee71/ee76 + (ee60 - (ee60 + ee75/ee76))/ee34 -
//    1) * ee7/ee2 + 1) * ee15 + ee54 * ee58/ee23 + (ee19 -
//    ee79/ee9)/ee7 + ee38 - (ee60 - (1/(ee7 * ee22) + ee11/ee76) * ee2)/ee10) * ee44/ee23 -
//    (ee60 - (1/(ee30 * R_pow(ee5, ee82)) +
//    ee60) * ee2 * ee11/xi) * ee2/ee30) * ee28/xi);
// out(j, 10) = 0.5 * (ee58 * ee31 * ee28/ee88);
// out(j, 11) = (((1/(ee30 * R_pow(ee5, ee87)) + ee19) * ee2 * ee11/xi -
//    ee19) * ee11/ee89 - (1 - ee2/ee7)/ee7) * ee2 - ((((ee19 +
//    ee75/ee32 - ee19)/ee34 - ee71/ee32) * ee7 * ee11/xi -
//    (1 + ee74 - 2 * ee7)/ee10) * ee15/ee2 + R_pow(ee54, 2)/ee23 +
//    2 * (ee43 * ee45 * ee2/(ee7 * ee10)) - ((ee45 * ee2 - ee19) * ee11/xi -
//    (ee19 - (2/ee12 + ee38) * ee2)/ee7)/ee10) * ee44/ee23;
// out(j, 12) = 0.5 * (ee54 * ee31/ee23);
// out(j, 13) = -ee86;

ee1 = exp(ldelta);
ee2 = exp(lpsi);
ee3 = xi * y;
ee4 = ee3/ee2;
ee5 = 1 + ee4;
ee6 = ee1/xi;
ee7 = 1 + ee1;
ee8 = R_pow(ee5, ee6);
ee9 = 1/xi;
ee10 = log1p(ee4);
ee11 = ee7 * ee8;
ee12 = R_pow(ee5, ee9);
ee13 = 1/ee11;
ee14 = 1 - ee13;
ee15 = 1 + ee6;
ee16 = R_pow(ee5, ee15);
ee17 = ee14 * ee7;
ee18 = 1/ee1;
ee19 = ee18 - 1;
ee20 = ee19 * ee1;
ee21 = xi * ee8;
ee22 = 1 + ee9;
ee23 = 1/ee8;
ee24 = ee10/ee21;
ee25 = ee20/xi;
ee27 = R_pow(ee5, ee25);
ee28 = 1 - ee17/(ee12 * ee1);
ee29 = R_pow(ee5, ee22);
ee30 = 1 - ee23;
ee31 = exp(lkappa);
ee33 = 1/ee16;
ee34 = ee1 * ee10;
ee35 = ee24 - y/(ee16 * ee2);
ee36 = ee1 - 1;
ee38 = 1 - ee7/ee1;
ee39 = ee5 * ee2;
ee41 = 1 - 0.5 * ee31;
ee42 = ee13 + ee24;
ee43 = ee36/xi;
ee44 = ee7/xi;
ee45 = R_pow(ee5, 2 + ee6);
ee46 = xi * ee16;
ee48 = R_pow(ee5, ee44 + 1);
ee49 = R_pow(ee5, ee43);
ee50 = ee45 * ee2;
ee53 = ee38 * ee14/ee12 + ee42/ee12;
ee55 = ee17/(ee27 * ee1) - 1/ee12;
ee57 = ee17/(ee29 * ee1) - 1/ee48;
ee59 = ee34/ee46;
ee61 = ee10/(xi * ee12) - y/(ee29 * ee2);
ee62 = ee10/xi;
ee64 = R_pow(ee5, ee9 + 2) * ee2;
ee66 = R_pow(xi, 2);
ee67 = y * ee22;
ee69 = y * ee15/ee50;
ee71 = y/ee39 - 2 * ee62;
ee72 = (ee20/ee27 - 1/ee27) * ee10;
ee73 = ee33 - ee59;
ee74 = 2 * ee1;
ee75 = ee34/ee21;
ee77 = (ee71/ee8 + ee34 * ee35/xi)/xi + y * (ee73/xi + ee69)/ee2;
ee78 = ee38 * ee1;
ee79 = ee7 * ee12;
ee80 = (ee74 - 1)/xi;
ee81 = R_pow(ee35, 2);
ee82 = 0.5 * (ee31 * log(ee28));
ee83 = 1/ee49;
ee84 = 1/ee29;
ee85 = ee33 - ee3 * ee15/ee50;
ee86 = ee33 + ee59;
ee87 = ee23 - (ee23 + ee75);
ee88 = 2 * ee6;
ee90 = ee34/(ee66 * ee16) - ee69;
ee91 = xi * ee30;
ee92 = xi * ee22;
ee93 = ee3/ee39;
ee94 = ee67/ee39;
ee95 = ee67/ee64;

out(j, 0) = 1 - y * (ee57 * ee41/ee28 + (ee92 - ee1/(ee30 * ee8))/ee5)/ee2;
out(j, 1) = ee94 - (ee55 * ee41 * ee35/ee28 + ee62 - ee1 * ee61/(ee30 * ee49))/xi;
out(j, 2) = 1 - (ee53 * ee41/ee28 + (1/ee7 + ee10/(ee91 * ee8)) * ee1);
out(j, 3) = -(ee82 + 1);
out(j, 4) =  - (y * (((ee84 - ee3 * ee22/ee64)/ee49 - y * (ee36/ee45 +
   ee1/(ee30 * R_pow(ee5, 2 * ee15)))/ee2) * ee1/ee30 +
   (ee85/ee12 + y * (R_pow(ee57, 2)/ee28 - 2/R_pow(ee5, ee44 +
   2))/ee2 - (ee85/ee27 - y * ee19 * ee1/ee64) * ee14 * ee7/ee1) * ee41/ee28 +
   ee92 * (ee93 - 1)/ee5)/ee2);
out(j, 5) =  - (y * (((ee55 * ee57/ee28 - 2/ee29) * ee35/xi +
   (ee20 * ee35/(xi * R_pow(ee5, ee25 + 1)) + ee90/ee27) * ee14 * ee7/ee1 -
   ee90/ee12) * ee41/ee28 + ((1 - ee93) * ee22 -
   ee9)/ee5 - ((ee36/R_pow(ee5, ee43 + 1) + ee1/(ee30 * R_pow(ee5, ee80 +
   1))) * ee61/xi + (ee10/(ee66 * ee29) - ee95)/ee49) * ee1/ee30)/ee2);
out(j, 6) =  - (y * (((((ee33 - ee86)/ee27 - ee72/ee46) * ee7/ee1 +
   ee38/ee29) * ee14 + ee53 * ee57/ee28 + ee42/ee29 - (ee78/(ee7 * ee48) +
   (ee33 - (1/(ee7 * ee16) + ee10/ee46) * ee1)/ee12)) * ee41/ee28 +
   ((1/(ee30 * R_pow(ee5, 1 + ee88)) +
   ee33) * ee1 * ee10/xi - ee33) * ee1/ee30)/ee2);
out(j, 7) = 0.5 * (y * ee57 * ee31/(ee28 * ee2));
out(j, 8) =  - ((((ee77/ee27 + ee20 * ee81/(xi * R_pow(ee5, (ee18 -
   2) * ee1/xi))) * ee14 * ee7/ee1 + (R_pow(ee55, 2)/ee28 -
   2/ee27) * ee81/xi - ee77/ee12) * ee41/ee28 + ee71/xi -
   (((ee71/ee12 + ee10 * ee61/xi)/xi + y * ((ee84 - ee10/(xi * ee29))/xi +
   ee95)/ee2)/ee49 + (ee36/R_pow(ee5, (ee1 - 2)/xi) +
   ee1/(ee30 * R_pow(ee5, 2 * ee43))) * R_pow(ee61, 2)/xi) * ee1/ee30)/xi +
   y * (1/ee66 + ee94)/ee39);
out(j, 9) =  - ((((((ee87 * ee10/xi + y * (ee86 - ee33)/ee2)/ee27 -
   ee72 * ee35/xi) * ee7/ee1 + ee38 * ee35/ee27) * ee14 +
   (ee53 * ee55/ee28 + ee42/ee27 - ee78/ee79) * ee35 - ((ee23 -
   ee75) * ee10/xi - (ee1 * ee35/ee7 + y * ee73/ee2))/ee12) * ee41/ee28 -
   (ee83 - (1/(ee30 * R_pow(ee5, ee80)) + ee83) * ee1 * ee10/xi) * ee1 * ee61/ee30)/xi);
out(j, 10) = 0.5 * (ee55 * ee31 * ee35/(xi * ee28));
out(j, 11) = (((1/(ee30 * R_pow(ee5, ee88)) + ee23) * ee1 * ee10/xi -
   ee23) * ee10/ee91 - (1 - ee1/ee7)/ee7) * ee1 - (R_pow(ee53, 2)/ee28 +
   ((ee23 - ee42 * ee1) * ee10/xi + (ee23 -
   (2/ee11 + ee24) * ee1)/ee7)/ee12 + 2 * (ee38 * ee42 * ee1/ee79) -
   ((ee87/ee27 - ee72/ee21) * ee7 * ee10/xi + (1 + ee74 -
   2 * ee7)/ee12) * ee14/ee1) * ee41/ee28;
out(j, 12) = 0.5 * (ee53 * ee31/ee28);
out(j, 13) = -ee82;
    

}

return out;

}

// //' @rdname egpd4d0
// [[Rcpp::export]]
arma::mat egpd4d34(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, arma::vec yvec, const arma::uvec dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec ldeltavec = X3 * Rcpp::as<arma::vec>(pars[2]);
arma::vec lkappavec = X4 * Rcpp::as<arma::vec>(pars[3]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 55);

if (dcate == 1) {
  lpsivec = lpsivec.elem(dupid);
  xivec = xivec.elem(dupid);
  ldeltavec = ldeltavec.elem(dupid);
  lkappavec = lkappavec.elem(dupid);
}

double y, lpsi, xi, ldelta, lkappa;

double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28;
double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee37, ee38, ee39;
double ee40, ee41, ee42, ee43, ee44, ee45, ee46, ee47, ee48;
double ee50, ee51, ee53, ee54, ee55, ee56, ee57, ee58, ee59;
double ee60, ee61, ee63, ee64, ee65, ee66, ee67, ee68, ee69;
double ee72, ee74, ee75, ee76, ee77, ee78;
double ee80, ee81, ee82, ee83, ee84, ee85, ee86, ee87, ee88, ee89;
double ee90, ee91, ee92, ee93, ee94, ee95, ee96, ee97, ee98, ee99;
double ee100, ee101, ee102, ee104, ee106, ee107, ee108, ee109;
double ee110, ee111, ee112, ee113, ee116, ee118;
double ee120, ee121, ee124, ee125, ee126, ee129;
double ee131, ee132, ee134, ee135, ee136, ee138, ee139;
double ee141, ee142, ee143, ee144, ee145, ee147, ee148, ee149;
double ee150, ee151, ee152, ee153, ee154, ee155, ee158, ee159;
double ee160, ee161, ee164, ee165, ee167, ee169;
double ee173, ee175, ee176, ee178;
double ee180, ee181, ee182, ee183, ee184, ee186;
double ee190, ee192, ee195, ee196, ee197, ee198, ee199;
double ee200, ee201, ee204, ee205, ee206, ee209;
double ee212, ee213, ee214, ee215, ee216, ee217, ee219;
double ee222, ee224, ee225, ee226, ee227, ee228;
double ee230, ee231, ee233, ee234, ee235, ee239;
double ee240, ee241, ee242, ee243, ee244, ee245, ee246, ee247, ee248, ee249;
double ee251, ee254, ee255, ee258;
double ee262, ee265, ee266, ee268, ee269;
double ee270, ee271, ee272, ee276, ee279;
double ee281, ee282, ee283, ee284, ee285, ee286, ee287, ee288, ee289;
double ee290, ee291, ee292, ee293, ee294, ee296, ee297, ee298, ee299;
double ee302, ee303, ee305, ee307, ee308, ee309;
double ee310, ee316, ee317;
double ee320, ee321, ee324, ee325, ee327, ee328, ee329;
double ee330, ee331, ee332, ee335, ee336, ee338, ee339;
double ee340, ee341, ee343, ee344, ee345, ee346, ee347, ee348, ee349;
double ee350, ee351, ee352, ee353, ee354, ee355, ee358, ee359;
double ee360, ee363, ee364, ee366, ee367;
double ee370, ee371, ee372, ee373, ee374, ee375, ee376, ee377;
double ee380, ee383, ee385, ee386, ee387, ee388, ee389;
double ee390, ee392, ee399;
double ee402, ee404, ee407;
double ee413, ee414, ee416, ee417, ee418, ee419;
double ee420, ee421, ee422, ee426, ee427;
double ee430, ee431, ee432, ee433, ee434, ee438, ee439;
double ee440, ee441, ee442, ee443, ee444, ee446;
double ee450, ee451, ee454, ee456, ee457, ee459;
double ee461, ee462, ee465, ee466, ee467, ee468, ee469;
double ee470, ee472, ee476;
double ee483, ee487, ee489;
double ee490, ee494, ee499;
double ee501, ee502, ee504, ee506, ee508, ee509;
double ee512, ee513, ee516, ee517, ee518, ee519;
double ee520, ee521, ee522, ee523, ee524, ee525, ee526;
double ee530, ee532, ee535, ee537, ee539;
double ee542, ee546;
double ee551, ee553, ee555, ee556, ee558;
double ee561, ee563, ee568, ee569;
double ee570, ee572, ee573, ee574, ee575, ee576, ee577, ee579;
double ee580, ee581, ee584, ee587, ee588, ee589;
double ee590, ee592, ee594, ee596;
double ee600, ee601, ee604, ee605, ee607, ee608;
double ee610, ee611, ee614, ee615, ee618, ee619;
double ee621, ee623, ee624, ee625, ee626;
double ee631, ee632, ee633, ee634, ee635, ee636, ee637, ee638, ee639;
double ee640, ee643, ee644, ee648, ee649;
double ee650, ee652, ee654, ee655, ee656, ee658;
double ee661, ee662, ee663, ee664, ee665, ee666, ee667, ee668, ee669;
double ee670, ee671, ee672, ee674, ee675, ee676, ee677;
double ee680, ee683, ee686, ee688, ee689;
double ee690, ee693, ee694, ee695, ee696, ee697, ee698, ee699;
double ee700, ee701, ee702, ee703, ee704, ee705, ee707, ee708, ee709;
double ee710, ee712, ee713, ee714, ee715, ee716, ee717, ee718, ee719;
double ee721, ee722, ee723, ee724, ee726, ee727, ee728;
double ee731, ee732, ee734, ee735, ee737, ee739;
double ee741, ee742, ee743, ee744, ee745, ee746, ee747, ee748, ee749;
double ee750, ee752, ee753, ee754, ee755, ee756, ee759;
double ee760, ee761, ee762, ee763, ee764, ee766, ee768, ee769;
double ee771, ee773, ee774, ee775, ee776, ee777, ee778, ee779;
double ee780, ee781, ee784, ee787;
double ee791, ee793, ee794, ee795, ee796, ee797;
double ee801, ee802, ee803, ee806, ee807, ee808, ee809;
double ee810, ee811, ee812, ee813, ee814, ee815, ee816, ee818;
double ee820, ee822, ee823, ee825, ee827, ee828;
double ee831, ee832, ee833, ee834, ee836, ee837, ee838, ee839;
double ee842, ee843, ee845, ee846, ee847, ee848;
double ee852, ee853, ee854, ee856, ee857, ee858, ee859;
double ee863, ee864, ee866, ee868, ee869;
double ee872, ee873, ee875, ee879;
double ee883, ee888, ee889;
double ee890, ee891, ee892, ee893, ee895, ee896, ee897, ee898, ee899;
double ee900, ee901, ee902, ee903, ee904, ee905, ee906, ee907, ee908, ee909;
double ee910, ee911, ee912, ee913, ee914, ee916, ee917;
double ee921, ee923, ee926;
double ee931, ee934, ee938, ee939;
double ee940, ee943, ee947, ee948, ee949;
double ee953, ee954, ee957, ee958;
double ee961, ee962, ee963, ee967;
double ee978, ee979;
double ee986, ee987, ee988, ee989;
double ee992, ee995;
double ee1011, ee1012, ee1013, ee1014, ee1015, ee1016, ee1017, ee1018, ee1019;
double ee1020, ee1021, ee1022, ee1024, ee1025, ee1026, ee1029;
double ee1030, ee1037, ee1038, ee1039;
double ee1040, ee1041, ee1042;
double ee1054, ee1058;
double ee1066, ee1068;
double ee1073, ee1074, ee1075, ee1076, ee1077, ee1078;
double ee1080, ee1081, ee1082, ee1083, ee1084, ee1085, ee1086, ee1087, ee1088, ee1089;
double ee1090, ee1092, ee1093, ee1094, ee1095, ee1096, ee1097, ee1098, ee1099;
double ee1100, ee1101, ee1102, ee1103, ee1104, ee1105, ee1106, ee1107, ee1108, ee1109;
double ee1110, ee1112, ee1113, ee1114, ee1115, ee1116, ee1117;
double ee1120, ee1121, ee1122, ee1123, ee1124, ee1125, ee1126, ee1128;
double ee1130, ee1131, ee1132, ee1133, ee1134, ee1135, ee1136, ee1137, ee1139;
double ee1140, ee1142, ee1144, ee1146, ee1149;
double ee1150, ee1152, ee1154, ee1156;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
ldelta = ldeltavec[j];
lkappa = lkappavec[j];

ee1 = exp(lpsi);
ee2 = exp(ldelta);
ee3 = xi * y;
ee4 = ee3/ee1;
ee5 = 1 + ee4;
ee6 = ee2/xi;
ee7 = log1p(ee4);
ee8 = 1 + ee6;
ee9 = R_pow(ee5, ee6);
ee10 = R_pow(ee5, ee8);
ee11 = 1/ee2;
ee12 = 1 + ee2;
ee13 = 1/xi;
ee14 = ee11 - 1;
ee15 = ee2 * ee7;
ee16 = xi * ee9;
ee17 = ee14 * ee2;
ee18 = ee17/xi;
ee19 = 2 + ee6;
ee20 = R_pow(ee5, ee18);
ee21 = 1/ee10;
ee22 = ee7/ee16;
ee23 = R_pow(ee5, ee19);
ee24 = R_pow(ee5, ee13);
ee25 = 1/ee9;
ee26 = 1 + ee13;
ee27 = ee12 * ee9;
ee28 = 1/ee27;
ee30 = y/(ee10 * ee1);
ee31 = ee22 - ee30;
ee32 = xi * ee10;
ee33 = R_pow(ee5, ee26);
ee34 = 1 - ee28;
ee35 = ee23 * ee1;
ee36 = ee15/ee32;
ee37 = ee11 - 2;
ee38 = ee5 * ee1;
ee39 = R_pow(xi, 2);
ee40 = ee37 * ee2;
ee41 = y/ee38;
ee42 = ee15/ee16;
ee43 = y * ee8;
ee44 = ee43/ee35;
ee45 = ee40/xi;
ee46 = ee7/xi;
ee47 = ee2 - 1;
ee48 = R_pow(ee5, ee45);
ee50 = ee41 - 2 * ee46;
ee51 = ee13 + 2;
ee53 = 1 - ee12/ee2;
ee54 = ee34 * ee12;
ee55 = R_pow(ee5, ee51);
ee56 = 1/ee20;
ee57 = ee47/xi;
ee58 = 2/ee10;
ee59 = ee21 - ee36;
ee60 = ee28 + ee22;
ee61 = ee55 * ee1;
ee63 = ee17/ee20 - ee56;
ee64 = ee39 * ee10;
ee65 = ee15/ee64;
ee66 = R_pow(ee5, ee57);
ee67 = ee21 + ee36;
ee68 = ee24 * ee2;
ee69 = ee65 - ee44;
ee72 = 1 - ee54/ee68;
ee74 = ee7/(xi * ee24) - y/(ee33 * ee1);
ee75 = 1/ee23;
ee76 = ee25 + ee42;
ee77 = 2/ee9;
ee78 = ee3 * ee8;
ee80 = (ee50/ee9 + ee15 * ee31/xi)/xi + y * (ee59/xi + ee44)/ee1;
ee81 = ee2 - 2;
ee82 = ee63 * ee7;
ee83 = ee25 - ee76;
ee84 = 3 + ee6;
ee85 = ee25 - ee42;
ee86 = R_pow(ee5, ee84);
ee87 = ee78/ee35;
ee88 = ee21 - ee87;
ee89 = ee39 * ee23;
ee90 = ee81/xi;
ee91 = y * ee26;
ee92 = ee91/ee61;
ee93 = 1 - ee25;
ee94 = ee86 * ee1;
ee95 = 2 * ee2;
ee96 = R_pow(ee5, ee90);
ee97 = ee12/xi;
ee98 = 1/ee33;
ee99 = R_pow(ee31, 2);
ee100 = 1/ee48;
ee101 = ee21 + ee58;
ee102 = 1 + ee95;
ee104 = y * ee19/ee94;
ee106 = R_pow(ee5, ee97 + 1);
ee107 = R_pow(ee5, ee18 + 1);
ee108 = ee25 + ee77;
ee109 = ee15/ee89;
ee110 = ee109 - ee104;
ee111 = ee21 - ee67;
ee112 = 1/ee24;
ee113 = xi * ee33;
ee116 = ee53 * ee34/ee24 + ee60/ee24;
ee118 = ee40/ee48 - ee100;
ee120 = ee54/(ee20 * ee2) - ee112;
ee121 = ee11 - 3;
ee124 = ee121 * ee2;
ee125 = ee75 - ee8 * ee7/ee23;
ee126 = ee7/ee32;
ee129 = ee83 * ee7/xi + y * (ee67 - ee21)/ee1;
ee131 = ee85 * ee7/xi;
ee132 = ee102 - 2 * ee12;
ee134 = y * ee59/ee1;
ee135 = ee12 * ee24;
ee136 = ee83/ee20;
ee138 = ee54/(ee33 * ee2) - 1/ee106;
ee139 = ee7/(ee39 * ee33);
ee141 = y * ee125/ee1;
ee142 = ee124/xi;
ee143 = xi * ee48;
ee144 = R_pow(ee5, ee142);
ee145 = ee7/ee113;
ee147 = ee136 - ee82/ee16;
ee148 = ee69/ee20;
ee149 = 1/ee66;
ee150 = 2/ee27;
ee151 = ee139 - ee92;
ee152 = ee12 * ee10;
ee153 = ee88/ee20;
ee154 = ee13 + 3;
ee155 = 2 * ee41;
ee158 = y * ee14 * ee2/ee61;
ee159 = ee50/ee10;
ee160 = xi * ee107;
ee161 = ee8 * ee110;
ee164 = ee17 * ee31/ee160 + ee148;
ee165 = ee3 * ee26;
ee167 = ee129/ee20 - ee82 * ee31/xi;
ee169 = (ee50/ee24 + ee7 * ee74/xi)/xi + y * ((ee98 - ee145)/xi +  ee92)/ee1;
ee173 = y * (ee2/ee89 - ee161)/ee1;
ee175 = ee53 * ee2;
ee176 = R_pow(ee5, ee154);
ee178 = ee111/ee20 - ee82/ee32;
ee180 = ee2 * ee31/ee12;
ee181 = R_pow(ee5, ee57 + 1);
ee182 = R_pow(ee74, 2);
ee183 = ee58 - ee36;
ee184 = 2/ee23;
ee186 = ee80/ee20 + ee17 * ee99/ee143;
ee190 = ee147 * ee12 * ee7/xi + ee132/ee24;
ee192 = ee153 - ee158;
ee195 = (ee25 - ee60 * ee2) * ee7/xi + (ee25 - (ee150 +  ee22) * ee2)/ee12;
ee196 = 2/xi;
ee197 = exp(lkappa);
ee198 = ee39 * ee55;
ee199 = xi * ee66;
ee200 = ee98 - ee165/ee61;
ee201 = ee15/ee199;
ee204 = ee167 * ee12/ee2 + ee53 * ee31/ee20;
ee205 = ee131 - (ee180 + ee134);
ee206 = ee196 + ee41;
ee209 = ee178 * ee12/ee2 + ee53/ee33;
ee212 = (ee155 - 6 * ee46)/xi + y * ee206/ee38;
ee213 = ee21 - (1/ee152 + ee126) * ee2;
ee214 = ee176 * ee1;
ee215 = xi * ee96;
ee216 = ee75 + ee184;
ee217 = ee108 - ee42;
ee219 = ((ee56 - ee118/ee9) * ee14 * ee2 - ee56) * ee7/xi;
ee222 = ee85 * ee14 * ee2 - ee25;
ee224 = y * ee51/ee214;
ee225 = ee219 + ee222/ee48;
ee226 = ee101 - ee36;
ee227 = ee2 - 3;
ee228 = ee7/ee198;
ee230 = ee69 * ee7 + ee159;
ee231 = ee228 - ee224;
ee233 = ee190 * ee34/ee2;
ee234 = ee195/ee24;
ee235 = ee230 * ee2;
ee239 = ee131 - ee134;
ee240 = ee227/xi;
ee241 = 2 * (ee53 * ee60 * ee2/ee135);
ee242 = 2 * ee44;
ee243 = R_pow(ee5, ee240);
ee244 = ee149 - ee201;
ee245 = 1/ee55;
ee246 = 2/ee33;
ee247 = 4/ee10;
ee248 = ee235/ee39;
ee249 = 1/ee96;
ee251 = ee8 * (ee75 + xi * ee110) + ee75;
ee254 = ee186 * ee34 * ee12/ee2;
ee255 = ee80/ee24;
ee258 = ee192 * ee34 * ee12/ee2;
ee262 = ee239 * ee14 * ee2 + ee30 - ee22;
ee265 = ee164 * ee34 * ee12/ee2;
ee266 = ee248 + ee173;
ee268 = R_pow(ee5, ee97 + 2);
ee269 = ee88/ee24;
ee270 = ee69/ee24;
ee271 = xi * ee144;
ee272 = xi * ee23;
ee276 = ((ee80 * ee7 + 2 * (ee31 * ee50)) * ee2/xi - ee212/ee9)/xi +  y * (((((ee183/xi + ee44) * ee7 - ee159) * ee2 - ee58)/xi -  ee242)/xi - ee173)/ee1;
ee279 = ee118 * ee14 * ee2 * ee7;
ee281 = ee40 * ee99/ee271;
ee282 = ee67 - ee58;
ee283 = ee58 + ee36;
ee284 = ee204 * ee34;
ee285 = ee205/ee24;
ee286 = ee175/ee135;
ee287 = ee60/ee20;
ee288 = ee65 - y * ee251/ee1;
ee289 = ee3 * ee19;
ee290 = ee209 * ee34;
ee291 = ee175/(ee12 * ee106);
ee292 = ee60/ee33;
ee293 = ee213/ee24;
ee294 = ee47 * ee7;
ee296 = 1 - 0.5 * ee197;
ee297 = 2 * ee6;
ee298 = xi * ee181;
ee299 = ee262/ee48;
ee302 = ee234 + ee241 - ee233;
ee303 = ee269 - (ee258 + 2 * (y/(ee268 * ee1)));
ee305 = ee226 * ee7/xi;
ee307 = ee270 + 2 * (ee31/ee113) - ee265;
ee308 = ee216 - ee289/ee94;
ee309 = ee254 - (ee255 + 2 * (ee99/(xi * ee20)));
ee310 = ee291 + ee293;
ee316 = (ee77 + ee42 - ee108) * ee2 * ee7/xi + ee25 + ee25 -  ee77;
ee317 = xi * ee93;
ee320 = ee78 * ee308/ee1 - ee21;
ee321 = ee169/ee66;
ee324 = ee217 * ee2 * ee7/xi;
ee325 = ee50/ee33;
ee327 = ee284 + (ee287 - ee286) * ee31 - ee285;
ee328 = ee290 + ee292;
ee329 = ee26 * ee231;
ee330 = (ee95 - 1)/xi;
ee331 = 2 * ee8;
ee332 = ee100 + 2/ee48;
ee335 = 4/ee9;
ee336 = ee15/ee272;
ee338 = y * (1/ee198 - ee329)/ee1;
ee339 = ee328 - ee310;
ee340 = ee321 + ee47 * ee182/ee215;
ee341 = ee111 * ee14;
ee343 = ee217 * ee7/xi;
ee344 = ee36 - ee21;
ee345 = ee42 - ee25;
ee346 = xi * ee2;
ee347 = ee200/ee66;
ee348 = ee151/ee66;
ee349 = ee25 - ee324;
ee350 = ee77 - ee42;
ee351 = 4 * ee2;
ee352 = xi * ee243;
ee353 = R_pow(ee116, 2);
ee354 = R_pow(ee5, ee45 + 1);
ee355 = ee347 - y * ee47/ee35;
ee358 = ee47 * ee74/ee298 + ee348;
ee359 = ee294/ee215;
ee360 = ee72 * ee1;
ee363 = ee59 * ee14 * ee2 - ee21;
ee364 = ee249 - ee359;
ee366 = ee80/ee48 + ee281;
ee367 = R_pow(ee5, ee297);
ee370 = ee350 * ee7/xi - y * ee183/ee1;
ee371 = ee67 - ee101;
ee372 = 2/ee20;
ee373 = ee15/xi;
ee374 = ee299 - ee279 * ee31/xi;
ee375 = ee116 * ee120;
ee376 = (ee305 - ee141) * ee2;
ee377 = (ee371 * ee7/xi + ee141) * ee2;
ee380 = ee9 * ee85 + ee373 - 1;
ee383 = ee40 * ee31/(xi * ee354) + ee69/ee48;
ee385 = ee81 * ee182/ee352;
ee386 = 2 * ee36;
ee387 = 2/ee66;
ee388 = ee77 - ee76;
ee389 = xi * ee72;
ee390 = ee39 * ee86;
ee392 = ee363/ee48 - ee279/ee32;
ee399 = (ee388 * ee7/xi + y * ee282/ee1) * ee2 * ee7/xi +  ee83 * ee50 + y * (ee377 + ee21 - ee21)/ee1;
ee402 = R_pow(ee5, ee90 + 1);
ee404 = ee88/ee48 - y * ee37 * ee2/(R_pow(ee5, ee18 + 2) * ee1);
ee407 = ee282 * ee7/xi + ee141;
ee413 = (ee283 - ee101) * ee2 * ee7/xi + ee21 + ee21 - ee58;
ee414 = ee101 - ee283;
ee416 = ee21 + ee2 * (ee126 + ee141) - ee21;
ee417 = 2 * ee57;
ee418 = 2 * ee42;
ee419 = 2/ee55;
ee420 = 4 * ee12;
ee421 = ee7 * ee151;
ee422 = xi * ee8;
ee426 = (ee80 * ee332 + ee281) * ee14 * ee2 * ee31/xi;
ee427 = ee116 * ee138;
ee430 = ee316 * ee7/xi + y * (ee414 * ee2 * ee7/xi + ee21 +  ee21 - ee58)/ee1;
ee431 = R_pow(ee5, ee417);
ee432 = R_pow(ee5, ee331);
ee433 = ee364 * ee2;
ee434 = ee88 * ee31;
ee438 = ee47/ee96;
ee439 = ee288/ee20;
ee440 = ee325 + ee421;
ee441 = 1/ee144;
ee442 = ee21 - ee226 * ee2 * ee7/xi;
ee443 = 2 * ee22;
ee444 = 4 * ee102;
ee446 = ee276/ee20 + ee426;
ee450 = (ee370 * ee2 * ee7/xi + ee85 * ee50)/xi + y * ((ee21 -  ee376)/xi + ee44)/ee1;
ee451 = ee380/ee24;
ee454 = (ee434/ee143 - y * (ee383/ee10 + ee69/ee107)/ee1) *  ee14 * ee2 + ee439;
ee456 = (ee136 - (ee225/ee9 + ee63 * (2 * ee83 + ee22))) *  ee7/xi;
ee457 = ee440/ee39;
ee459 = ee125 * ee2;
ee461 = ee349 * ee7/xi;
ee462 = ee69 * ee31;
ee465 = 1 + ee420 + ee351 - ee444;
ee466 = ee149 + ee387;
ee467 = 2 * (ee102/ee9);
ee468 = 2 * ee126;
ee469 = ee2/ee9;
ee470 = ee422/ee23;
ee472 = y * ee442/ee1;
ee476 = (ee366/ee10 + 2 * (ee462/ee48)) * ee14 * ee2/xi +  ee266/ee20;
ee483 = ((ee129/ee48 - ee118 * ee7 * ee31/xi) * ee14 * ee2 +  ee299) * ee31 + ee399/ee20 - ee80 * ee63 * ee7;
ee487 = ee430/ee20 - (ee225 * ee31 + 2 * (ee129 * ee63)) *  ee7/xi;
ee489 = ee80 * ee2/ee12;
ee490 = R_pow(ee120, 2);
ee494 = ee456 + ee316/ee20 - ee451;
ee499 = ee63 * ee88 * ee7/xi + ee416/ee20 + y * (ee392/ee10 +  ee341 * ee2/ee107)/ee1;
ee501 = ee413/ee20 - (ee225/ee10 + 2 * (ee63 * ee111)) *  ee7/xi;
ee502 = ee457 + ee338;
ee504 = ee26 * (ee245 + xi * ee231) + ee245;
ee506 = R_pow(ee5, ee330);
ee508 = R_pow(ee5, 4 + ee6) * ee1;
ee509 = ee88/ee12;
ee512 = ee124/ee144 - ee441;
ee513 = ee330 + 1;
ee516 = ee183 * ee7/xi - ee141;
ee517 = ee320/ee20;
ee518 = 1/ee181;
ee519 = ee245 + ee419;
ee520 = ee108 - ee418;
ee521 = ee467 + 2 * (ee12 * ee85);
ee522 = 2 * ee85;
ee523 = 8 * ee469;
ee524 = ee201 - ee149;
ee525 = ee3 * ee51;
ee526 = ee3/ee38;
ee530 = y * (ee404/ee10 + 2 * (ee88/ee107)) * ee14 * ee2/ee1;
ee532 = y * (ee459 + ee470)/ee1;
ee535 = ee483 * ee12/ee346 + ee186 * ee53;
ee537 = (ee487 * ee12 - ee132 * ee31/ee20)/ee2 + 2 * (ee167 *  ee53);
ee539 = ((ee169 * ee7 + 2 * (ee74 * ee50))/xi - ee212/ee24)/xi +  y * (((((ee246 - ee145)/xi + ee92) * ee7 - (ee325 + ee246))/xi -  2 * ee92)/xi - ee338)/ee1;
ee542 = (ee494 * ee12/ee2 + 3 * (ee147 * ee53)) * ee7/xi +  ee465/ee68;
ee546 = (ee407/ee20 - ee341 * ee31/ee48) * ee2 + ee63 *  ee69 * ee7 - ee374/ee10;
ee551 = (((ee468 + 2/ee152) * ee2 - ee101)/ee12 - ee305) *  ee2 + ee21;
ee553 = ee450 - ee489;
ee555 = ee169/ee96 + ee385;
ee556 = ee120 * ee138;
ee558 = (ee132/ee33 - ee501 * ee12)/ee2 - 2 * (ee178 * ee53);
ee561 = (ee509 + ee126) * ee2 + ee532 - ee21;
ee563 = ee192 * ee53 - ee499 * ee12/ee2;
ee568 = ee516/xi;
ee569 = R_pow(ee5, ee513);
ee570 = R_pow(ee5, ee57 + 2);
ee572 = (ee25 - (((ee521 - ee523)/ee12 + ee522)/ee12 + ee343) *  ee2)/ee12 + (ee25 - ((ee108 - (ee443 + ee150) * ee2)/ee12 +  ee343) * ee2) * ee7/xi;
ee573 = ee461 - ((ee520 * ee7/xi - (2 * ee180 + y * (ee101 -  ee386)/ee1)) * ee2/ee12 + ee472);
ee574 = ee69/ee12;
ee575 = ee517 + ee530;
ee576 = ee56 + ee372;
ee577 = ee98 + ee246;
ee579 = ee155 + 6/xi;
ee580 = ee139 - y * ee504/ee1;
ee581 = xi * ee402;
ee584 = ee165 * (ee519 - ee525/ee214)/ee1 - ee98;
ee587 = ee546 * ee12/ee346 - ee164 * ee53;
ee588 = ee375 * ee31;
ee589 = ee427/ee72;
ee590 = ee353/ee72;
ee592 = (ee574 - ee568) * ee2 + ee44;
ee594 = ee433 + ee438;
ee596 = ee200/ee96 - y * ee81/(ee570 * ee1);
ee600 = ee294/ee298 - ee518;
ee601 = ee47/ee181;
ee604 = ee81 * ee74/ee581 + ee151/ee96;
ee605 = ee112 + 2/ee24;
ee607 = 1/ee86 - ee19 * ee7/ee86;
ee608 = 2 * ee336;
ee610 = ee15/ee390 - y * ee84/ee508;
ee611 = ee2/ee272;
ee614 = ee446 * ee34 * ee12/ee2;
ee615 = ee276/ee24;
ee618 = ee454 * ee34 * ee12/ee2;
ee619 = ee302 * ee72;
ee621 = ee120 * ee307;
ee623 = ee164/ee10 + ee69/ee33;
ee624 = ee600 * ee2;
ee625 = ee93 * ee345;
ee626 = R_pow(ee5, 1 + ee297);
ee631 = ee332 * ee14 * ee2 * ee99/xi + 2 * (ee80 * ee576);
ee632 = ee47/ee23;
ee633 = ee288/ee24;
ee634 = ee249 + 2/ee96;
ee635 = ee101 + ee247;
ee636 = 2 * ee153;
ee637 = 2 * ee526;
ee638 = ee39 * ee176;
ee639 = ee535 * ee34;
ee640 = ee537 * ee34;
ee643 = ee476 * ee34 * ee12/ee2;
ee644 = ee542 * ee34;
ee648 = ee204 * ee2/ee12 + ee205/ee20;
ee649 = (ee190 * ee31 + (2 * (ee204 * ee60) - 2 * (ee205 *  ee53/ee24)) * ee2)/ee12;
ee650 = (ee190/ee10 + (2 * (ee209 * ee60) - 2 * (ee53 *  ee213/ee24)) * ee2)/ee12;
ee652 = ee551/ee24 + ee558 * ee34;
ee654 = ee553/ee24 + ee80 * ee53 * ee2/ee135;
ee655 = ee186 * ee60;
ee656 = ee186/ee10;
ee658 = (((ee100 - ee512/ee9) * ee37 * ee2 - ee100) * ee7/xi +  (ee85 * ee37 * ee2 - ee25)/ee144) * ee14 * ee2;
ee661 = ee209 * ee2/ee152 + ee213/ee33;
ee662 = ee561/ee24;
ee663 = ee563 * ee34;
ee664 = (ee195 * ee53 * ee605 * ee2 - 3 * (ee190 * ee60))/ee12;
ee665 = ee619 - ee353;
ee666 = ee266/ee24;
ee667 = ee80/ee33;
ee668 = R_pow(ee138, 2);
ee669 = ee192 * ee60;
ee670 = ee192 * ee101;
ee671 = ((ee75 - ee336) * ee7/xi - y * ee607/ee1) * ee8;
ee672 = ee572/ee24;
ee674 = ee195 * ee31/ee20;
ee675 = ee195/ee33;
ee676 = ee234 + 2 * ee590;
ee677 = ee573/ee24;
ee680 = (ee110 * ee7 + ee50/ee23) * ee2/ee39 + y * (ee2/ee390 -  ee19 * ee610)/ee1;
ee683 = ee575 * ee34 * ee12/ee2;
ee686 = ee53 * ee88 * ee2/ee135;
ee688 = ee625 - ee15/(xi * ee367);
ee689 = ee200 * ee74;
ee690 = ee577 * ee88;
ee693 = (ee25 + 2/(ee93 * ee367)) * ee2 * ee7/xi;
ee694 = ee320/ee24;
ee695 = 2 * ee589;
ee696 = 2 * ee611;
ee697 = 2 * ee104;
ee698 = ee155 + ee196;
ee699 = ee58 + ee247;
ee700 = ee77 + ee335;
ee701 = ee2 * ee182;
ee702 = ee539/ee66;
ee703 = ee587 * ee34;
ee704 = ee327 * ee138;
ee705 = ee327 * ee72;
ee707 = (ee204/ee10 + ee209 * ee31) * ee2/ee12;
ee708 = ee309 * ee138;
ee709 = (ee169 * ee634 + ee385) * ee47;
ee710 = ee618 + 2 * (y * ee623/ee1);
ee712 = ee662 + ee663 + ee669;
ee713 = ee592/ee24;
ee714 = ee169 * ee244;
ee715 = ee375/ee72;
ee716 = ee116 * ee303;
ee717 = ee116 * ee307;
ee718 = ee556 * ee31;
ee719 = ee120 * ee303;
ee721 = ee490 * ee99/xi;
ee722 = ee490/ee72;
ee723 = ee594 * ee182;
ee724 = (ee689/ee215 - y * (ee604/ee33 + ee151/ee181)/ee1) *  ee47;
ee726 = ee596/ee33 + 2 * (ee200/ee181);
ee727 = ee205/ee33;
ee728 = ee164 * ee60;
ee731 = ee175 * ee69/ee135;
ee732 = ee93 * ee244;
ee734 = R_pow(ee5, ee13 + 4) * ee1;
ee735 = ee244 * ee200;
ee737 = ee213 * ee31/ee20;
ee739 = ee579/xi + y * ee698/ee38;
ee741 = ee74 * ee151;
ee742 = ee580/ee66;
ee743 = ee584/ee66;
ee744 = ee149 + ee2 * (ee201 - ee466) * ee7/xi;
ee745 = 1/ee243;
ee746 = ee21 + ee247;
ee747 = ee108 + ee335;
ee748 = 2 * ee164;
ee749 = 2 * ee148;
ee750 = ee2 * ee74;
ee752 = y * (ee632 - ee624/ee33)/ee1;
ee753 = ee614 - (ee615 + ee631 * ee31/xi);
ee754 = (ee707 + ee727 + ee737)/xi;
ee755 = ee639 + ee655;
ee756 = ee702 + ee709 * ee74/xi;
ee759 = (ee555/ee33 + 2 * (ee741/ee96)) * ee47;
ee760 = ee309 * ee72;
ee761 = ee339 * ee72;
ee762 = ee654 + 2 * (ee648 * ee31/xi);
ee763 = ee340 * ee93;
ee764 = ee340/ee10;
ee766 = ee664 + ee672 - ee644;
ee768 = ee302 * ee120 * ee31;
ee769 = ee302 * ee138;
ee771 = ee266 * ee7 + 2 * (ee69 * ee50);
ee773 = ee502/ee66;
ee774 = ee714 + ee723/xi;
ee775 = ee716/ee72;
ee776 = ee556/ee72;
ee777 = ee719/ee72;
ee778 = ee621 * ee31;
ee779 = ee668/ee72;
ee780 = (ee433/ee33 + ee601) * ee74;
ee781 = ee724 + ee742;
ee784 = ee676 + 2 * ee302 + ee241 - ee233;
ee787 = ee512 * ee37 * ee2 * ee7;
ee791 = ee212/ee10;
ee793 = ((24 * ee46 - 6 * ee41)/xi - y * ee579/ee38)/xi -  y * ee739/ee38;
ee794 = ee358 * ee74;
ee795 = ee93 * ee344;
ee796 = ee735 - ee752;
ee797 = ee244/ee10;
ee801 = ee693 - ee25;
ee802 = ee324 - ee25;
ee803 = ee19/ee86;
ee806 = (ee700 - ee42) * ee2 * ee7/xi;
ee807 = ee81 * ee7;
ee808 = ee524 * ee151;
ee809 = ee524/ee33;
ee810 = ee344/ee66;
ee811 = ee743 + y * ee726 * ee47/ee1;
ee812 = 2 * ee102;
ee813 = 2/ee268;
ee814 = 2/ee86;
ee815 = 8 * ee2;
ee816 = 8/ee9;
ee818 = ee7/ee638 - y * ee154/ee734;
ee820 = y * ee661/ee1;
ee822 = y * (((ee608 - ee184)/xi - ee697) * ee2/ee39 - ee680 *  ee8)/ee1;
ee823 = y * ee138;
ee825 = y * ((ee75 - ee608)/xi + ee104 - ee671)/ee1;
ee827 = ee754 + ee703 + ee731;
ee828 = ee755 - ee762;
ee831 = ee640 + ee649 + ee674 - ee677;
ee832 = ee327 * ee116;
ee833 = ee704 + ee339 * ee120 * ee31;
ee834 = ee705 - ee588;
ee836 = ee759/xi + ee773;
ee837 = ee309 * ee116;
ee838 = ee708 - ee778;
ee839 = ee760 - ee721;
ee842 = (ee656 + ee667 + (ee748 + ee749) * ee31)/xi + ee666 -  ee643;
ee843 = ee339 * ee116;
ee845 = ee650 + ee675 - ee652;
ee846 = ee763 - ee701/(xi * ee431);
ee847 = ee764 + ee794/ee66;
ee848 = ee254 + (2 * ee722 - ee372) * ee99/xi;
ee852 = (ee219 + ee222 * (ee100 - 2 * ee118) - ee658/ee9) *  ee7/xi + ((ee349 * ee14 + ee443) * ee2 + ee25 - ee77)/ee48 -  ee380/ee20;
ee853 = ee712 + 2 * ee820;
ee854 = (ee376/xi - ee242) * ee7;
ee856 = (ee239 * ee37 * ee2 + ee30 - ee22)/ee144 - ee787 *  ee31/xi;
ee857 = ee321 + (ee438 + 2 * (ee2/(ee93 * ee431))) * ee182/xi;
ee858 = ee718/xi;
ee859 = ee621/ee72;
ee863 = ee358 * ee93;
ee864 = ee307 * ee72;
ee866 = (ee50/ee55 + ee7 * ee231)/ee39 + y * (1/ee638 -  ee51 * ee818)/ee1;
ee868 = ee732 + ee15/(xi * ee506);
ee869 = ee93 * ee1;
ee872 = ee244/ee9;
ee873 = (ee745 - ee807/ee352) * ee47;
ee875 = ee124 * ee99/(xi * R_pow(ee5, (ee11 - 4) * ee2/xi));
ee879 = (ee636 - ee158) * ee31/xi + ee633 - ee710;
ee883 = (ee58 - ee809) * ee2 * ee7/xi - ee21;
ee888 = ee227 * ee182/(xi * R_pow(ee5, (ee2 - 4)/xi));
ee889 = ee808 - ee780/xi;
ee890 = ee810 - ee797;
ee891 = ee344/ee9;
ee892 = ee345/ee66;
ee893 = ee345/ee10;
ee895 = ee694 + y * (ee670 + ee690)/ee1 - ee683;
ee896 = ee101 + ee58;
ee897 = ee635 - (ee699 - ee36) * ee2 * ee7/xi;
ee898 = ee216 - (ee8 * (ee216 - ee336) + ee696) * ee7;
ee899 = ee747 - ee806;
ee900 = 2 - ee637;
ee901 = 2 * ee665;
ee902 = 2 * ee715;
ee903 = 2 * ee224;
ee904 = 2/(ee93 * ee506);
ee905 = 3 * ee2;
ee906 = 4 * ee761;
ee907 = 4 * ee5;
ee908 = 8 * ee588;
ee909 = 8 * ee427;
ee910 = 8 * ee353;
ee911 = 8/ee10;
ee912 = ee15/ee215;
ee913 = y * ee2;
ee914 = -(0.5 * (ee197 * log(ee72)));
ee916 = ((ee276 * ee7 + 3 * (ee80 * ee50) - 3 * (ee212 *  ee31)) * ee2/xi - ee793/ee9)/xi + y * ((((ee791 + 2 *  (ee230/xi) - ee771) * ee2 + (ee235 + 6/ee10)/xi - 6 *  ee69)/xi + 3 * ee173)/xi - ee822)/ee1;
ee917 = ee753 * ee120;
ee921 = ((ee704 + (ee328 + ee695 - ee310) * ee120 * ee31)/xi -  ee717)/ee72 + ee713 + ee728 - ee827;
ee923 = ee614 + ((ee848 + 2 * ee309 - ee255) * ee120/ee72 -  ee631) * ee31/xi - ee615;
ee926 = ((ee708 + (ee265 + (2 * ee776 - ee246) * ee31/xi -  ee270) * ee120 * ee31)/ee72 - (ee656 + (ee859 + ee748 +  ee749) * ee31 + ee667))/xi + ee643 - ee666;
ee931 = ee639 + ee837/ee72 + ee655 + (ee120 * (2 * ee327 +  2 * (ee588/ee72))/ee72 - 2 * ee648) * ee31/xi - ee654;
ee934 = ((ee856 * ee14 + ee262 * ee37/ee144) * ee31 - ee80 *  ee118 * ee14 * ee7) * ee2/xi + (ee450 * ee14 * ee2 -  ee80)/ee48;
ee938 = ee640 + (ee768 + ee116 * ((ee902 + 2 * ee287 - 2 *  ee286) * ee31 + 2 * ee284 - 2 * ee285))/ee72 + ee649 +  ee674 - ee677;
ee939 = ee539/ee96;
ee940 = ee276/ee48;
ee943 = ee832 + ee768;
ee947 = (ee658 * ee31 + 2 * (ee262 * ee118)) * ee7/xi;
ee948 = ee843 + ee769;
ee949 = ee339 * ee303;
ee953 = (ee769 + ee116 * (2 * ee290 + ee695 + 2 * ee292 -  (2 * ee291 + 2 * ee293)))/ee72 + ee650 + ee675 - ee652;
ee954 = ee665 * ee120;
ee957 = (ee771 - ee791) * ee2/ee39 + ee822;
ee958 = ee774/ee66;
ee961 = (ee169 * (ee745 + 2/ee243) + ee888) * ee81 * ee74/xi;
ee962 = ee340 * ee7;
ee963 = ee555 * ee47;
ee967 = (ee80 * (ee441 + 2/ee144) + ee875) * ee37 * ee2 *  ee31/xi;
ee978 = (ee777 + ee636 - ee158) * ee31/xi + ee633 + y *  (ee138 * (2 * (ee718/ee389) - 2 * ee307)/ee72 - 2 * ee623)/ee1 -  ee618;
ee979 = ee781/ee10;
ee986 = ee664 + ee784 * ee116/ee72 + ee672 - ee644;
ee987 = ((ee461 - ee472) * ee14 * ee2 + y * (ee58 - (ee21 +  ee386))/ee1 - (ee77 - (ee25 + ee418)) * ee7/xi)/ee48;
ee988 = ((ee183 * (3 * ee46 - ee41) - (ee854 + ee141))/xi -  ee825) * ee2;
ee989 = ee502 * ee7;
ee992 = ee169 * ee364 + (ee873 + ee81/ee243) * ee182/xi;
ee995 = ee775 + ee686 + y * (ee138 * (2 * ee339 + ee695)/ee72 -  2 * ee661)/ee1 - ee712;
ee1011 = (ee9 + 2 * ee9) * ee85;
ee1012 = ee796/ee10;
ee1013 = ee355 * ee93;
ee1014 = ee355 * ee74;
ee1015 = ee355 * ee7;
ee1016 = ee363 * ee37;
ee1017 = ee303 * ee307;
ee1018 = ee303 * ee72;
ee1019 = R_pow(ee303, 2);
ee1020 = ee407 * ee63;
ee1021 = ee407 * ee2;
ee1022 = ((ee75 - ee611) * ee7/xi + ee104 - ee671) * ee2;
ee1024 = ee676 + ee241 - ee233;
ee1025 = ee147 * ee132;
ee1026 = ee802/ee9;
ee1029 = ee212/ee33;
ee1030 = ee413 * ee14;
ee1037 = ((ee700 - (ee108 + ee42)) * ee2 * ee7/xi + ee25 +  ee25 + ee25 + ee77 + ee77 + ee77 - (ee77 + ee816)) *  ee2 * ee7/xi + ee25 + ee25 + ee77 - ee335;
ee1038 = ((ee247 - ee36) * ee7/xi - ee141) * ee2;
ee1039 = ee358 * ee7;
ee1040 = ee890 * ee93;
ee1041 = ee248 + y * ((ee75 - 2 * (ee7/ee23)) * ee2/ee39 +  ee697 - ee8 * (ee109 + xi * ee680 - ee104))/ee1;
ee1042 = ee811/ee10;
ee1054 = (ee466 - ee201) * ee2 * ee7/xi - ee149;
ee1058 = ee897 * ee7/xi;
ee1066 = ee693 + 2 * ee345 - ee25;
ee1068 = ee899 * ee7/xi;
ee1073 = (2 * ee777 + ee636 - ee158) * ee31/xi;
ee1074 = (ee905 - 1)/xi;
ee1075 = (ee351 - 1)/xi;
ee1076 = ee47 * (ee912 - ee249);
ee1077 = ee891 + ee893;
ee1078 = ee892 - ee872;
ee1080 = ee694 + y * (ee138 * (ee269 + 2 * ee303 + y * (2 *  ee779 - ee813)/ee1 - ee258)/ee72 + ee670 + ee690)/ee1 -  ee683;
ee1081 = 0.5 * ((ee284 + (ee715 + ee287 - ee286) * ee31 -  ee285) * ee197/ee389);
ee1082 = 0.5 * ((ee254 + (ee722 - ee372) * ee99/xi - ee255) *  ee197/ee389);
ee1083 = 0.5 * ((ee590 + ee234 + ee241 - ee233) * ee197/ee72);
ee1084 = 0.5 * (ee116 * ee197/ee72);
ee1085 = 0.5 * (ee120 * ee197 * ee31/ee389);
ee1086 = 0.5 * (y * ((ee776 - ee246) * ee31/xi + ee265 -  ee270) * ee197/ee360);
ee1087 = 0.5 * (y * (ee290 + ee589 + ee292 - ee310) * ee197/ee360);
ee1088 = 0.5 * (ee823 * ee197/ee360);
ee1089 = 0.5 * (y * (ee269 + y * (ee779 - ee813)/ee1 - ee258) *  ee197/ee360);
ee1090 = 1 + ee331;
ee1092 = 1 + ee351;
ee1093 = 1 + ee3 * (ee637 - 3)/ee38;
ee1094 = ee149 - (ee149 + ee904) * ee2 * ee7/xi;
ee1095 = ee21 - ee78 * (ee216 + 4/ee23 - ee289 * (ee814 +  4/ee86 - ee3 * ee84/ee508)/ee1)/ee1;
ee1096 = ee896 - ee386;
ee1097 = ee101 + ee36;
ee1098 = ee216 - (ee8 * ee216 * ee7 + ee3 * (ee8 * ee607 +  ee803)/ee1);
ee1099 = ee25 + ee335;
ee1100 = 2 * (ee833 * ee72);
ee1101 = 2 * ee703;
ee1102 = 2 * ee834;
ee1103 = 2 * ee839;
ee1104 = ee901 + ee910;
ee1105 = 2 * ee713;
ee1106 = 2 * ee775;
ee1107 = 2 * (ee717/ee72);
ee1108 = 2 * ee728;
ee1109 = 2 * ee731;
ee1110 = 2 * ee370;
ee1112 = 2 * (1 + 2 * ee4) + ee907;
ee1113 = ee812 + ee420;
ee1114 = 2 * ee47;
ee1115 = 2 * (ee2/ee96);
ee1116 = 2 * ee141;
ee1117 = 2/ee176;
ee1120 = 3 * ee6;
ee1121 = 4 * ee705;
ee1122 = ee906 + ee909;
ee1123 = 4 * ee864;
ee1124 = 4 * ee6;
ee1125 = 4/ee66;
ee1126 = ee247 - (ee746 - ee36) * ee2 * ee7/xi;
ee1128 = ee247 + ee911 - ee101 * ee2 * ee7/xi;
ee1130 = ee335 + ee816 - ee108 * ee2 * ee7/xi;
ee1131 = 6 * ee4;
ee1132 = 8 * ee858;
ee1133 = 8 * ee721;
ee1134 = 8 * ee4;
ee1135 = xi * ee26;
ee1136 = xi * ee569;
ee1137 = xi * ee626;
ee1139 = ee3 * (3 - ee637)/ee38;
ee1140 = ee526 - 1;
ee1142 = y * ((ee8 * (ee184 + ee336 - ee216) - ee696) *  ee7 + ee75 + ee184 - ee184)/ee1;
ee1144 = y * (((2 * (ee7/(xi * ee55)) - ee419)/xi - ee903)/ee39 -  ee866 * ee26)/ee1;
ee1146 = y * ee668/ee1;
ee1149 = y * (ee251 + ee184 + xi * (2 * ee161 - y * (ee8 *  (ee814 + xi * ee19 * ee610) + ee803)/ee1))/ee1 - ee65;
ee1150 = ee91/ee38;
ee1152 = ee43 * ee216/ee1;
ee1154 = ee43 * (ee184 - ee75)/ee1;
ee1156 = y * ee898/ee1;

out(j, 0) =  - (y * ((ee743 + y * ((ee355 * ee101 - 2 * (ee913/(ee93 * R_pow(ee5, ee1090 +
   ee6) * ee1))) * ee2/ee93 + ee726 * ee47)/ee1) * ee2/ee93 +
   ee1080 * ee296/ee72 + ee1135 * ee1093/ee5)/ee1);
out(j, 1) =  - (y * (ee978 * ee296/ee72 + ((ee1014/ee199 - y * (2 * (ee358/ee10) +
   2 * (ee750/(ee317 * R_pow(ee5, ee57 +
   ee331))))/ee1) * ee2/ee93 + ee724 + ee742) * ee2/ee93 + (ee26 * (ee1139 -
   1) - ee1140/xi)/ee5)/ee1);
out(j, 2) =  - (y * (ee995 * ee296/ee72 + (ee735 + ee2 * (y * (2 * (ee344/ee10) +
   2 * (ee15/(ee317 * R_pow(ee5, ee331 +
   ee6))))/ee1 - ee1015/ee16)/ee93 - ee752) * ee2/ee93)/ee1);
out(j, 3) = ee1089;
out(j, 4) =  - (y * (ee926 * ee296/ee72 - ((((ee764 + (2 * (ee358/ee66) +
   2 * (ee750/(ee317 * R_pow(ee5, (ee905 - 2)/xi +
   1)))) * ee74) * ee2/ee93 + ee759)/xi + ee773) * ee2/ee93 +
   2 * ((ee26 * ee900 - ee196)/R_pow(ee5, 2))))/ee1);
out(j, 5) =  - (y * (ee921 * ee296/ee72 + (((ee1039/ee9 + (ee810 +
   2 * (ee15/(ee317 * R_pow(ee5, ee1074 + 1))) - ee797) * ee74) * ee2/ee93 -
   ee780)/xi + ee808) * ee2/ee93)/ee1);
out(j, 6) = ee1086;
out(j, 7) =  - (y * (ee953 * ee296/ee72 + ((ee58 - ((ee893 +
   2 * ee891 + 2 * (ee15/(ee317 * R_pow(ee5, 1 + ee1120))))/ee93 +
   ee809)) * ee2 * ee7/xi - ee21) * ee2/ee93)/ee1);
out(j, 8) = ee1087;
out(j, 9) = ee1088;
out(j, 10) =  - ((ee923 * ee296/ee72 - ((ee702 + ((ee340 * ee466 +
   2 * (ee701/(ee317 * R_pow(ee5, 3 * ee57)))) * ee2/ee93 +
   ee709) * ee74/xi) * ee2/ee93 + ee212/xi))/xi - y * (ee206/ee39 +
   y * (1/ee39 + 2 * ee1150)/ee38)/ee38);
out(j, 11) =  - ((ee931 * ee296/ee72 - ((ee723 + ((2 * (ee244/ee66) -
   2 * (ee15/(ee317 * R_pow(ee5, (ee1114 + ee2)/xi)))) * ee182 -
   ee962/ee9) * ee2/ee93)/xi + ee714) * ee2/ee93)/xi);
out(j, 12) = ee1082;
out(j, 13) =  - ((ee938 * ee296/ee72 - (((ee892 + 2 * (ee15/(ee317 * R_pow(ee5, ee1074))) -
   2 * ee872)/ee93 + ee201 - ee466) * ee2 * ee7/xi +
   ee149) * ee2 * ee74/ee93)/xi);
out(j, 14) = ee1081;
out(j, 15) = ee1085;
out(j, 16) = (((ee108 - ((ee108 * ee345 + 2 * (ee15/(ee317 * R_pow(ee5, ee1120))))/ee93 +
   ee42)) * ee2 * ee7/xi - ee25) * ee7/ee317 -
   (1 - (3 - 2 * (ee2/ee12)) * ee2/ee12)/ee12) * ee2 -
   ee986 * ee296/ee72;
out(j, 17) = ee1083;
out(j, 18) = ee1084;
out(j, 19) = ee914;
out(j, 20) =  - (y * (((ee98 - ee165 * (ee519 + 4/ee55 - ee525 * (ee1117 +
   4/ee176 - ee3 * ee154/ee734)/ee1)/ee1)/ee66 +
   y * (((ee584/ee96 + y * ((ee200/ee243 - y * ee227/(R_pow(ee5, ee90 +
   2) * ee1))/ee33 + 2 * (ee200/ee402)) * ee81/ee1)/ee33 +
   (ee518 + 2/ee181) * ee584 - 3 * (ee596 * ee200)) * ee47 +
   (2 * ee1042 + y * ((2 * ((ee1013 + ee913/(ee432 * ee1))/ee432) +
   4 * (ee1013/ee432) - 8 * (ee913/(R_pow(ee5, 2 + ee331 +
   ee297) * ee1)))/ee93 + 4 * (ee355/ee432)) * ee2/ee869 -
   ((ee347 - y * (ee632 + 2 * (ee2/(ee93 * ee432)))/ee1) * ee355 +
   2 * (R_pow(ee355, 2) - ee1042))) * ee2/ee93)/ee1) * ee2/ee93 +
   (ee1095/ee24 + y * ((ee138 * (2 * ee694 + y * (ee138 * ((2 * (ee1018 -
   ee1146) + 4 * ee1018 + 8 * ee1146)/ee72 +
   6 * ee303)/ee72 + 2 * ee670 + 2 * ee690)/ee1 - 2 * ee683) +
   ee1019 + 2 * (ee138 * ee895 + ee1019))/ee72 + 4 * (ee575/ee10) +
   4 * (ee320/ee33) - 6 * (ee192 * ee88))/ee1 - (ee1095/ee20 +
   y * ((ee320/ee48 + y * ((ee88/ee144 - y * ee121 * ee2/(R_pow(ee5, ee45 +
   2) * ee1))/ee10 + 2 * (ee88/ee354)) * ee37 * ee2/ee1)/ee10 +
   (1/ee107 + 2/ee107) * ee320 - 3 * (ee404 * ee88)) * ee14 * ee2/ee1) * ee34 * ee12/ee2) * ee296/ee72 +
   ee1135 * (ee3 * (7 + ee3 * ((ee1134 - ee1112)/ee5 - 6)/ee38)/ee38 -
   1)/ee5)/ee1);
out(j, 21) =  - (y * (((ee120 * ee895/ee72 + 2 * ee517 + ee530) * ee31/xi +
   ee1149/ee24 + y * ((ee138 * (ee1073 + (4 * (ee719 * ee31/xi) +
   ee823 * (ee1132 - (2 * (ee858 + ee864) +
   ee1123))/ee360)/ee72 + ee633 - ee710) + 2 * (ee138 * ee879 -
   ee1017) - ee1017)/ee72 + ee454 * ee101 + ee577 * ee288 + 3 * (ee192 * ee69) +
   3 * (ee164 * ee88))/ee1 - ((ee31 * ee320/ee143 +
   y * (((ee434/ee271 - y * ((ee124 * ee31/(xi * R_pow(ee5, ee142 +
   1)) + ee69/ee144)/ee10 + ee69/ee354)/ee1) * ee37 * ee2 +
   ee288/ee48)/ee10 + ee404 * ee69 + 2 * (ee383 * ee88) +
   2 * (ee288/ee107))/ee1) * ee14 * ee2 + ee1149/ee20) * ee34 * ee12/ee2) * ee296/ee72 +
   ((ee811 * ee74/ee199 + y * (ee979 +
   ((ee601 + 2 * (ee2/(ee93 * ee569))) * ee74/xi + ee348) * ee355 +
   (4 * (ee1014/ee1136) + y * (2 * ((ee750/ee1136 -
   ee863)/ee432) - (4 * (ee863/ee432) + 8 * (ee750/(xi * R_pow(ee5, ee513 +
   ee331)))))/ee869) * ee2/ee93 + 2 * (ee979 + ee355 * ee358))/ee1) * ee2/ee93 +
   (ee74 * ee584/ee215 + y * (((ee689/ee352 -
   y * ((ee227 * ee74/(xi * R_pow(ee5, ee240 + 1)) +
   ee151/ee243)/ee33 + ee151/ee402)/ee1) * ee81 + ee580/ee96)/ee33 +
   ee596 * ee151 + 2 * (ee604 * ee200) + 2 * (ee580/ee181))/ee1) * ee47 +
   (y * (ee504 + ee419 + xi * (2 * ee329 -
   y * (ee26 * (ee1117 + xi * ee51 * ee818) + ee51/ee176)/ee1))/ee1 -
   ee139)/ee66) * ee2/ee93 + (ee26 * (1 + ee3 * (ee3 * ((ee1112 -
   ee1134)/ee5 + 6)/ee38 - 7)/ee38) - ee1093/xi)/ee5)/ee1);
out(j, 22) =  - (y * ((ee116 * ee895/ee72 + ee175 * ee320/ee135 +
   y * ((ee949 + (ee686 + (6 * ee716 + ee823 * (2 * (ee761 -
   ee427) + ee906 + ee909)/ee360)/ee72 - ee853) * ee138 + 2 * (ee949 +
   (ee686 - ee853) * ee138))/ee72 + (ee563 * ee101 +
   3 * (ee209 * ee88)) * ee2/ee12 + 3 * (ee192 * ee213) - ee561 * ee577)/ee1 -
   (((ee320/ee12 - ee126) * ee2 + ee21 - y * (ee1098 * ee2 +
   ee422 * ee308)/ee1)/ee24 + (ee575 * ee53 - (ee82 * ee320/xi +
   (ee21 - (ee21 + ee2 * (ee126 + y * ee1098/ee1)))/ee20 -
   y * ((ee404 * ee111 - 2 * (ee416/ee107)) * ee14 * ee2 +
   2 * (ee392 * ee88) - ((ee118 * ee88 * ee14 * ee7/xi +
   y * (((ee59 * ee37 * ee2 - ee21)/ee144 - ee787/ee32) * ee14/ee10 +
   ee1016/ee354)/ee1) * ee2 + (ee17 * (ee36 + ee532 -
   ee21) + ee21 - ee87)/ee48)/ee10)/ee1) * ee12/ee2) * ee34 +
   ee575 * ee60)) * ee296/ee72 + (ee244 * ee584 - ((ee811 * ee7/ee16 +
   y * (ee355 * ((ee21 + 2/(ee93 * ee626)) * ee2 * ee7/xi -
   ee21) + (4 * (ee1015/ee1137) - y * (2 * ((ee795 - ee15/ee1137)/ee432) +
   4 * (ee795/ee432) + 8 * (ee15/(xi * R_pow(ee5, ee1090 +
   ee297))))/ee869) * ee2/ee93 + 2 * (ee355 * ee344 -
   ee1012) - ee1012)/ee1) * ee2/ee93 + y * (2 * ((ee624 -
   ee601) * ee200) - ((ee364 * ee200 - y * (ee81/ee570 - (ee807/ee581 -
   1/ee402) * ee47/ee33)/ee1) * ee2/ee33 + ee596 * ee47/ee33))/ee1)) * ee2/ee93)/ee1);
out(j, 23) = 0.5 * (y * ee1080 * ee197/ee360);
out(j, 24) =  - (y * ((((ee309 * ee303 + ee120 * (ee1073 + ee633 -
   ee710) * ee31)/ee72 + ee186 * ee88 + (ee120 * ee879/ee72 +
   2 * ee454 + 2 * ee439) * ee31 + ee80 * ee192)/xi + ee1041/ee24 +
   y * ((((ee120 * (ee1132 - ee1123) * ee31 + 2 * (ee839 * ee138))/ee72 -
   4 * ee778) * ee138/ee389 + 2 * (R_pow(ee307, 2) -
   ee842 * ee138))/ee72 - 2 * (ee476/ee10 + ee266/ee33 +
   2 * (ee164 * ee69)))/ee1 - (((ee366 * ee88 + 2 * (ee288 * ee31/ee48))/xi -
   y * ((((ee80/ee144 + ee875)/ee10 + 2 * (ee462/ee144)) * ee37 * ee2/xi +
   ee266/ee48)/ee10 + ee266/ee107 +
   2 * (ee383 * ee69))/ee1) * ee14 * ee2 + ee1041/ee20) * ee34 * ee12/ee2) * ee296/ee72 +
   (((ee857 * ee355 + 2 * (ee781 * ee74/ee66))/xi -
   y * ((((4 * (ee863/ee569) + 8 * (ee750/(xi * R_pow(ee5, ee417 +
   ee331)))) * ee74 + 2 * (ee846/ee432))/ee93 +
   4 * (ee794/ee569)) * ee2/ee317 + 2 * (ee836/ee10 +
   R_pow(ee358, 2)))/ee1) * ee2/ee93 + ((ee555 * ee200 + 2 * (ee74 * ee580/ee96))/xi -
   y * ((((ee169/ee243 + ee888)/ee33 +
   2 * (ee741/ee243)) * ee81/xi + ee502/ee96)/ee33 + ee502/ee181 +
   2 * (ee604 * ee151))/ee1) * ee47 + (ee457 + y * ((ee245 -
   2 * (ee7/ee55))/ee39 + ee903 - ee26 * (ee228 + xi * ee866 -
   ee224))/ee1)/ee66) * ee2/ee93 + (ee91 * (4 - ee3 * ((ee907 -
   ee1131)/ee5 + 6)/ee38)/ee38 - (2 * ee1139 - (2 + 2 * ee1140))/ee39)/ee5)/ee1);
out(j, 25) =  - (y * ((((ee327 * ee303 + (ee686 + ee1106 - ee853) * ee120 * ee31)/xi +
   ee116 * ee879)/ee72 + ((ee204 * ee88 +
   ee563 * ee31) * ee2/ee12 + ee192 * ee205 - ee561 * ee31/ee20)/xi +
   ((ee568 - ee288/ee12) * ee2 - y * (ee251 - ee1022)/ee1)/ee24 +
   ee175 * ee288/ee135 + y * ((((ee120 * ee1122 * ee31 +
   2 * (ee834 * ee138))/ee389 - 4 * ee717) * ee138/ee72 +
   2 * ((ee713 + ee728 - ee827) * ee138 - ee339 * ee307))/ee72 -
   2 * ((ee209 * ee69 - ee587/ee10) * ee2/ee12 + ee164 * ee213 -
   ee592/ee33))/ee1 - ((((ee374 * ee88 - (ee63 * ee288 * ee7 +
   ee416 * ee14 * ee2 * ee31/ee48))/xi + (((ee58 - ee67) * ee7/xi -
   ee141) * ee2/xi + y * (ee1022 + ee75 - ee75)/ee1)/ee20 -
   y * (ee392 * ee69 + (ee383 * ee111 - ee1021/ee160) * ee14 * ee2 -
   (((ee118 * ee69 * ee7 - ee856/ee10) * ee14 -
   ee1016 * ee31/ee144) * ee2/xi + (ee2 * (ee7/ee64 - (ee516 * ee2/xi -
   ee44) * ee14) - ee44)/ee48)/ee10)/ee1) * ee12/ee2 +
   ee454 * ee53) * ee34 + ee454 * ee60)) * ee296/ee72 + ((((ee796/ee66 +
   ee355 * ee1094) * ee74 - ee781 * ee7/ee9) * ee2/ee93 +
   ee594 * ee200 * ee74)/xi + ee244 * ee580 - y * (((((2 * (ee868/ee432) -
   (4 * (ee795/ee569) + 8 * (ee15/(xi * R_pow(ee5, ee330 +
   ee331))))) * ee74/ee93 - 4 * (ee1039/ee626)) * ee2/ee317 -
   2 * (ee358 * ee344 + ee889/ee10))/ee93 - ((ee359 -
   ee249) * ee151 - (ee873/ee33 + ee81/ee402) * ee74/xi)/ee33) * ee2 +
   (ee601 - ee624) * ee151 + ee604 * ee47/ee33)/ee1) * ee2/ee93)/ee1);
out(j, 26) = 0.5 * (y * ee978 * ee197/ee360);
out(j, 27) =  - (y * ((((((ee225 * ee88 - 2 * (ee63 * ee416)) * ee7/xi +
   ((ee414 * ee7/xi + ee1142) * ee2 + ee21 + ee21 -
   ee58)/ee20 + y * ((((ee442 * ee14 + ee468) * ee2 + ee21 -
   ee58)/ee48 - (ee658/ee10 + 2 * (ee363 * ee118)) * ee7/xi)/ee10 +
   ee1030 * ee2/ee107 + 2 * (ee392 * ee111))/ee1) * ee12 +
   ee192 * ee132)/ee2 + 2 * (ee499 * ee53)) * ee34 + (ee302 * ee303 +
   ee116 * (ee1106 + 2 * ee686 - (2 * ee662 + 2 * ee663 +
   2 * ee669 + 4 * ee820)))/ee72 + y * ((((ee116 * ee1122 +
   2 * (ee665 * ee138))/ee72 + 4 * ee843) * ee138/ee72 + 2 * (ee845 * ee138 +
   R_pow(ee339, 2)))/ee72 - 2 * (ee551/ee33 + (2 * (ee209 * ee213) -
   ee558/ee10) * ee2/ee12))/ee1 - ((ee190 * ee88 +
   (2 * (ee561 * ee53/ee24) + 2 * (ee563 * ee60)) * ee2)/ee12 +
   (((ee101 - ((2 * ee509 + ee468) * ee2 + y * (2 * ee459 +
   ee422 * ee216)/ee1))/ee12 + ee305) * ee2 + y * (ee898 * ee2 +
   ee470)/ee1 - ee21)/ee24 + ee192 * ee195)) * ee296/ee72 +
   (((ee355 * ee801 - 2 * (ee796/ee9)) * ee7/xi - y * (((2 * (ee688/ee432) +
   4 * (ee795/ee626) + 8 * (ee15/(xi * R_pow(ee5, ee331 +
   ee297))))/ee93 + 4 * (ee344/ee626)) * ee2 * ee7/ee317 +
   2 * (R_pow(ee344, 2) - ee883/ee10))/ee1) * ee2/ee93 +
   ee744 * ee200 - y * (ee632 - (((2 * (ee2/ee181) - ee1076/ee33) * ee7/xi -
   ee518)/ee33 + 2 * (ee600/ee33)) * ee2)/ee1) * ee2/ee93)/ee1);
out(j, 28) = 0.5 * (y * ee995 * ee197/ee360);
out(j, 29) = ee1089;
out(j, 30) =  - (y * ((((ee753 * ee138 + ((((ee138 * (ee1103 +
   ee1133) + 2 * (ee838 * ee72))/ee72 + 2 * ee838)/ee72 - ((2 * ee859 +
   4 * ee164 + 4 * ee148) * ee31 + 2 * ee656 + 2 * ee667))/xi +
   2 * ee643 - 2 * ee666) * ee120 * ee31 - 3 * (ee309 * ee307))/ee72 -
   (ee446/ee10 + (ee842 * ee120/ee72 + ee266 * ee576 +
   3 * ee476) * ee31 + ee276/ee33 + 3 * (ee186 * ee69) +
   3 * (ee80 * ee164)))/xi + (((ee940 + ee967)/ee10 + ee266 * ee332 * ee31 +
   3 * (ee366 * ee69)) * ee14 * ee2/xi + ee957/ee20) * ee34 * ee12/ee2 -
   ee957/ee24) * ee296/ee72 + (((6 * (1 -
   ee526) - 6)/xi + 2 * (y * ee900/ee38))/ee39 + y * (ee900/ee39 +
   y * ((2 * ee5 - ee1131)/ee5 + 4) * ee26/ee38)/ee38)/ee5 -
   (((ee756/ee10 + (ee836 * ee466 + ((2 * (ee846/ee569) +
   2 * (ee847 * ee93/ee66) + 8 * (ee701/(xi * R_pow(ee5, (ee351 -
   3)/xi + 1))))/ee93 + 2 * (ee847/ee66)) * ee2/ee317) * ee74 +
   (ee857 + 2 * ee340) * ee358) * ee2/ee93 + ((ee939 +
   ee961)/ee33 + ee502 * ee634 * ee74 + 3 * (ee555 * ee151)) * ee47)/xi +
   ((ee989 + 2 * (ee151 * ee50) - ee1029)/ee39 +
   ee1144)/ee66) * ee2/ee93)/ee1);
out(j, 31) =  - (y * ((((ee828 * ee138 + (((ee138 * (ee1102 +
   ee908) + ee1100)/R_pow(ee72, 2) - (2 * ee707 + 2 * ee727 +
   2 * ee737))/xi + ee1105 + ee1108 - (ee1101 + ee1107 + ee1109)) * ee120 * ee31 +
   ee309 * ee339 + 2 * (ee838 * ee116/ee72) -
   2 * (ee327 * ee307))/xi - ee842 * ee116)/ee72 + ee476 * ee60 +
   ((ee266/ee12 - ((ee854 + ee183 * ee50)/xi + ee825)/xi) * ee2 -
   ee173)/ee24 + ((2 * (ee587 * ee31) - (ee535/ee10 +
   ee209 * ee80 + 2 * (ee204 * ee69))) * ee2/ee12 + 2 * (ee592 * ee31/ee20) -
   (ee553/ee33 + ee186 * ee213 + 2 * (ee205 * ee164)))/xi -
   (((((((ee377/xi + ee1154) * ee7 + ee282 * ee50)/xi -
   ee825)/ee20 + ee14 * (2 * (ee1021 * ee31/ee143) - ee366 * ee111)) * ee2 +
   ee266 * ee63 * ee7 - (ee934/ee10 + 2 * (ee374 * ee69))) * ee12/ee346 -
   ee476 * ee53) * ee34 + ee266 * ee53 * ee2/ee135)) * ee296/ee72 +
   ((((ee836 * ee7/ee9 + ee340 * ee344 +
   (((ee904 + ee387) * ee2 * ee7/xi - ee387) * ee358 +
   2 * (ee889/ee66)) * ee74 + ((2 * (ee1040/ee66) + 8 * (ee15/(xi * R_pow(ee5, (ee351 -
   2)/xi + 1))) - 2 * (ee868/ee569)) * ee182/ee93 +
   2 * (ee847 * ee7/ee9)) * ee2/ee317 - ee774/ee10)/ee93 -
   ee992/ee33) * ee2 - (ee963/ee33 + 2 * (ee594 * ee74 * ee151)))/xi +
   ee502 * ee524) * ee2/ee93)/ee1);
out(j, 32) = 0.5 * (y * ee926 * ee197/ee360);
out(j, 33) =  - (y * ((((ee831 * ee138 + ee845 * ee120 * ee31 +
   ((ee116 * (ee1100 + 8 * (ee375 * ee138 * ee31)) + 2 * (ee954 * ee138 * ee31))/ee72 +
   2 * (ee833 * ee116))/ee72 + 2 * (ee327 * ee339))/xi +
   ee116 * (ee1105 + ee1108 - (2 * ee754 +
   ee1101 + ee1107 + ee1109)) - ee302 * ee307)/ee72 + (((((ee1096 * ee7/xi -
   ee1116)/xi - 2 * ee574) * ee2 - ee1152)/ee12 -
   (ee1126 * ee7/xi - ee1156)/xi) * ee2 + ee44)/ee24 + ((ee558 * ee31 -
   (ee537/ee10 + 2 * (ee204 * ee213) + 2 * (ee209 * ee205))) * ee2/ee12 -
   (ee551 * ee31/ee20 + ee573/ee33))/xi +
   (ee190 * ee69 + (2 * (ee592 * ee53/ee24) - 2 * (ee587 * ee60)) * ee2)/ee12 +
   ee195 * ee164 - (((((((ee746 - ee283) * ee2 * ee7/xi +
   ee58 + ee58 - ee746) * ee7/xi + ee1142)/ee20 -
   (ee1020 * ee7/xi + ee1030 * ee31/ee48)) * ee2 - ((ee987 -
   ee947)/ee10 + (ee1020 * ee2/xi - ee225 * ee69) * ee7 + 2 * (ee374 * ee111))) * ee12/ee2 +
   2 * (ee546 * ee53))/xi + ee164 * ee132/ee2) * ee34) * ee296/ee72 +
   ((((ee883/ee66 + 2 * (ee244 * ee344) -
   (((2 * (ee688/ee569) + 2 * (ee1040/ee9) + 8 * (ee15/(xi * R_pow(ee5, ee1075 +
   1))))/ee93 + 2 * (ee890/ee9)) * ee2 * ee7/ee317 +
   ee744/ee10)) * ee74 - (ee801 * ee358 +
   2 * (ee889/ee9)) * ee7) * ee2/ee93 - ((((ee1076 - ee1115) * ee7/xi +
   ee249)/ee33 + 2 * (ee364/ee33)) * ee2 + ee601) * ee74)/xi +
   ee1054 * ee151) * ee2/ee93)/ee1);
out(j, 34) = 0.5 * (y * ee921 * ee197/ee360);
out(j, 35) = ee1086;
out(j, 36) =  - (y * ((((ee542/ee10 + 3 * (ee209 * ee195) -
   (ee551 * ee53 * ee605 + 3 * (ee558 * ee60))) * ee2 + 3 * (ee190 * ee213))/ee12 +
   ((((ee138 * ee1104 + 2 * (ee948 * ee72))/ee72 +
   2 * ee948)/ee72 + ee650 + ee675 + 2 * ee845 - ee652) * ee116 +
   ee339 * ee784 + ee766 * ee138)/ee72 + ee572/ee33 -
   ((((((2 * ee59 + 2 * ee183 - (8 * (ee2/ee10) - (2 * (ee102/ee10) +
   2 * (ee12 * ee183)))/ee12)/ee12 + ee1128 * ee7/xi) * ee2 -
   ee635)/ee12 - ee1058) * ee2 + ee21)/ee24 + ((ee465/ee33 +
   3 * (ee178 * ee132) - ((((ee699 - ee1097) * ee2 * ee7/xi +
   ee21 + ee21 + ee21 + ee58 + ee58 + ee58 - (ee58 + ee911)) * ee2 * ee7/xi +
   ee21 + ee21 + ee58 - ee247)/ee20 - (ee852/ee10 +
   3 * (ee225 * ee111) + 3 * (ee63 * ee413)) * ee7/xi) * ee12)/ee2 -
   3 * (ee501 * ee53)) * ee34)) * ee296/ee72 +
   (((ee1066 * ee344 + ((2 * (ee688/ee626) + 2 * (ee1077 * ee93/ee9) +
   8 * (ee15/(xi * R_pow(ee5, 1 + ee1124))))/ee93 + 2 * (ee1077/ee9)) * ee2 * ee7/ee317 -
   (ee802/ee10 + ee883 * ee108))/ee93 +
   ee21 + ee58 - (ee1054/ee33 + ee577 * ee524)) * ee2 * ee7/xi -
   ee21) * ee2/ee93)/ee1);
out(j, 37) = 0.5 * (y * ee953 * ee197/ee360);
out(j, 38) = ee1087;
out(j, 39) = ee1088;
out(j, 40) =  - ((((ee916/ee20 + ((ee967 + 4 * ee940) * ee31 +
   3 * (ee366 * ee80)) * ee14 * ee2/xi) * ee34 * ee12/ee2 +
   ((ee917/ee72 - (ee276 * (ee372 + 4/ee20) + 2 * ee446 + 2 * ee426)) * ee31 +
   ((ee614 + (ee120 * ((ee1103 + 4 * ee760 + ee1133)/ee72 +
   4 * ee309)/ee72 - ee631) * ee31/xi - ee615) * ee120 * ee31 +
   ee309 * (ee848 - ee255) + 2 * (ee917 * ee31 +
   R_pow(ee309, 2)))/ee72 - 6 * (ee186 * ee80))/xi - ee916/ee24) * ee296/ee72 -
   (((((ee539 * ee7 + 3 * (ee169 * ee50) - 3 * (ee212 * ee74))/xi -
   ee793/ee24)/xi + y * (((ee1029 + (ee325 +
   2 * ee440 + 6/ee33 + ee421)/xi - (ee989 + (2 * ee50 + 6) * ee151))/xi +
   3 * ee338)/xi - ee1144)/ee1)/ee66 + (((ee961 +
   4 * ee939) * ee74 + 3 * (ee555 * ee169)) * ee47 + (ee857 * ee340 +
   (((2 * (ee846/ee431) + 4 * (ee763/ee431) + 8 * (ee701/(xi * R_pow(ee5, 4 * ee57))))/ee93 +
   4 * (ee340/ee431)) * ee2 * ee74/ee317 +
   2 * (ee756/ee66)) * ee74 + 2 * (ee756 * ee74/ee66 +
   R_pow(ee340, 2))) * ee2/ee93)/xi) * ee2/ee93 +
   ee793/xi))/xi + y * (ee739/ee39 + y * (ee698/ee39 + y * (2/ee39 +
   6 * ee1150)/ee38)/ee38)/ee38);
out(j, 41) =  - ((((((ee755 + (ee120 * (ee1102 + ee1121 + ee908) * ee31/ee389 +
   4 * ee837)/ee72 - ee762) * ee120 * ee31 +
   (ee284 + (ee287 + ee902 - ee286) * ee31 - ee285) * ee309 +
   2 * (ee828 * ee120 * ee31 + ee327 * ee309))/xi + ee753 * ee116)/ee72 +
   ((ee934 * ee31 + ((((((ee108 - ee76) * ee7/xi +
   y * ee371/ee1) * ee2 * ee7/xi + ee388 * ee50)/xi + y * ((((ee67 -
   ee247) * ee7/xi + ee141) * ee2 + ee58 - ee21)/xi + ee1154)/ee1) * ee7 +
   (ee1110 - 2 * ee31) * ee50) * ee2/xi + ee212 * (ee76 -
   ee25) + y * (((((ee247 - ee58)/xi - ee44) * ee7 +
   ee159 - ee141) * ee2 + ee58 - ee58)/xi + ee988)/ee1)/ee20 +
   (ee366 * ee129 + 2 * (ee399 * ee31/ee143)) * ee14 * ee2 +
   2 * (ee374 * ee80) - ee276 * ee63 * ee7) * ee12/ee346 + ee446 * ee53) * ee34 +
   ee446 * ee60 - (((((((ee343 - y * ee226/ee1) * ee2 * ee7/xi +
   ee350 * ee50)/xi + y * ((ee58 - ee1038)/xi +
   ee242)/ee1) * ee7 + 2 * (ee370 * ee50)) * ee2/xi - ee212 * ee85)/xi +
   y * ((ee988 + (ee1038 - ee58)/xi - ee242)/xi -
   ee173)/ee1 - ee276 * ee2/ee12)/ee24 + ee276 * ee53 * ee2/ee135 +
   (ee553 * ee576 * ee31 + (3 * (ee535 * ee31) + 3 * (ee204 * ee80)) * ee2/ee12 +
   3 * (ee186 * ee205))/xi)) * ee296/ee72 -
   ((((ee958 + ee340 * ee1094 + ((2 * (ee868/ee431) +
   4 * (ee732/ee431) - 8 * (ee15/(xi * R_pow(ee5, (3 * ee47 + ee2)/xi)))) * ee182/ee93 -
   4 * (ee962/ee506)) * ee2/ee317 + 2 * (ee958 +
   ee340 * ee244)) * ee74 - ee756 * ee7/ee9) * ee2/ee93 +
   (ee992 * ee2 + ee963 + 2 * (ee169 * ee594)) * ee74)/xi +
   ee539 * ee244) * ee2/ee93)/xi);
out(j, 42) = 0.5 * (ee923 * ee197/ee389);
out(j, 43) =  - ((((((((ee987 + ee430 * ee14 * ee2/ee48 - ee947) * ee31 +
   2 * (ee374 * ee129) - (ee225 * ee80 + 2 * (ee399 * ee63/xi)) * ee7)/xi +
   (2 * ((((ee76 - ee1099) * ee2 * ee7/xi +
   ee335 - ee77) * ee7/xi + y * ((ee746 - ee67) * ee2 * ee7/xi +
   ee58 - ee247)/ee1) * ee2 * ee7/xi + ((ee76 - ee108) * ee2 * ee7/xi +
   ee25 - ee25) * ee50 + y * ((((ee699 - ee67) * ee2 * ee7/xi +
   ee21 + ee58 - ee635) * ee7/xi + y * ((ee8 * (ee75 +
   ee336 - ee216) - ee696) * ee7 + ee75 + ee184 - ee75)/ee1) * ee2 +
   ee21 - ee21)/ee1) - 2 * ee399)/ee20) * ee12 -
   ee186 * ee132)/ee2 + 2 * (ee483 * ee53/xi)) * ee34 + ((((ee116 * (ee1121 +
   ee908) + 2 * (ee954 * ee31))/ee72 + 4 * ee832) * ee120 * ee31/ee72 +
   2 * (ee831 * ee120 * ee31 + R_pow(ee327, 2)))/xi +
   ee309 * ee1024 + 2 * (ee828 * ee116))/ee72 +
   (ee190 * ee80 + (2 * (ee535 * ee60) - 2 * (ee553 * ee53/ee24)) * ee2)/ee12 +
   ee186 * ee195 - (((((ee335 - (ee1099 -
   ee42) * ee2 * ee7/xi) * ee7/xi - y * ee1126/ee1) * ee2 * ee7/xi +
   ee349 * ee50)/xi + y * ((ee21 - (ee1058 - ee1156) * ee2)/xi +
   ee44)/ee1 - ((((ee108 + ee77 - ee418) * ee7/xi - y * ee1096/ee1) * ee2 * ee7/xi +
   ee520 * ee50)/xi + y * ((ee101 -
   ((ee635 - ee386) * ee7/xi - ee1116) * ee2)/xi + ee1152)/ee1 -
   2 * ee489) * ee2/ee12)/ee24 + 2 * (((ee537 * ee31 + 2 * (ee204 * ee205)) * ee2/ee12 +
   ee573 * ee31/ee20)/xi))) * ee296/ee72 -
   ((((ee340 * ee801 - 2 * (ee774/ee9)) * ee7 + (((2 * (ee688/ee431) +
   8 * (ee15/(xi * R_pow(ee5, (ee1114 + ee95)/xi))) -
   4 * (ee732/ee506))/ee93 - 4 * (ee244/ee506)) * ee2 * ee7/ee317 +
   2 * (R_pow(ee244, 2) + ee744/ee66)) * ee182) * ee2/ee93 +
   (((ee47 * (ee912 - ee634) - ee1115) * ee7/xi +
   ee249 + ee249 + ee249) * ee2 + ee438) * ee182)/xi + ee169 * ee744) * ee2/ee93)/xi);
out(j, 44) = 0.5 * (ee931 * ee197/ee389);
out(j, 45) = ee1082;
out(j, 46) =  - (((((ee640 + ((ee120 * ee1104 * ee31 + 2 * (ee943 * ee72))/ee72 +
   2 * ee943)/ee72 + ee649 + ee674 + 2 * ee831 -
   ee677) * ee116 + ee327 * ee784 + ee766 * ee120 * ee31)/ee72 +
   ((ee542 * ee31 + 3 * (ee537 * ee60) + 3 * (ee204 * ee195) -
   ee573 * ee53 * ee605) * ee2 + 3 * (ee190 * ee205))/ee12 +
   ((((ee1037 * ee7/xi + y * (((ee1097 - ee699) * ee2 * ee7/xi +
   ee21 + ee21 + ee58 + ee58 + ee247 - (ee896 + ee247)) * ee2 * ee7/xi +
   ee21 + ee21 + ee58 - ee247)/ee1)/ee20 -
   (ee852 * ee31 + 3 * (ee225 * ee129) + 3 * (ee430 * ee63)) * ee7/xi) * ee12 -
   (ee465 * ee31/ee20 + 3 * (ee167 * ee132)))/ee2 +
   3 * (ee487 * ee53)) * ee34 + ee572 * ee31/ee20 - ((ee25 -
   ee899 * ee2 * ee7/xi) * ee7/xi - (((ee747 - ee1130 * ee2 * ee7/xi) * ee7/xi -
   ((((ee812 - ee815) * ee31 + 2 * (ee370 * ee12))/ee12 +
   2 * ee239 + ee1110) * ee2/ee12 + y * (ee635 -
   ee1128 * ee2 * ee7/xi)/ee1)) * ee2/ee12 + y * (ee21 - ee897 * ee2 * ee7/xi)/ee1))/ee24) * ee296/ee72 -
   (((ee1066 * ee244 +
   ee802/ee66 - (((2 * (ee688/ee506) + 2 * (ee1078 * ee93/ee9) +
   8 * (ee15/(xi * R_pow(ee5, ee1075))))/ee93 + 2 * (ee1078/ee9)) * ee2 * ee7/ee317 +
   ee744 * ee108))/ee93 + (ee387 +
   ee1125 - ee201) * ee2 * ee7/xi - (ee466 + ee1125)) * ee2 * ee7/xi +
   ee149) * ee2 * ee74/ee93)/xi);
out(j, 47) = 0.5 * (ee938 * ee197/ee389);
out(j, 48) = ee1081;
out(j, 49) = ee1085;
out(j, 50) = ((((ee801 * ee345 + ((2 * (ee688/ee367) + 4 * (ee625/ee367) +
   8 * (ee15/(xi * R_pow(ee5, ee1124))))/ee93 +
   4 * (ee345/ee367)) * ee2 * ee7/ee317 + 2 * (R_pow(ee345, 2) -
   ee1026) - 2 * ee1026)/ee93 + ee25 + ee77 + ee335 - ee806) * ee2 * ee7/xi -
   ee25) * ee7/ee317 - (1 - (7 - ((ee1113 - ee815)/ee12 +
   6) * ee2/ee12) * ee2/ee12)/ee12) * ee2 - (((ee116 * ((ee901 +
   4 * ee619 + ee910)/ee72 + 4 * ee302)/ee72 + 2 * ee766) * ee116 +
   ee1024 * ee302 + 2 * (ee766 * ee116 + R_pow(ee302, 2)))/ee72 +
   ((ee1025 * ee7/xi - (1 + 8 * ee1113 +
   ee815 - (ee444 + 40 * ee12 + 6 * ee1092))/ee24)/ee2 - ((((ee456 +
   ee316 * (ee56 - 3 * ee63) - (ee852/ee9 + ee451 + 3 * (ee225 * ee83))) * ee7/xi +
   ee1037/ee20 - ((ee1011 + 2 * ee373 -
   3) * ee2 * ee7/xi + ee9 * ee349 + 2 - ee1011)/ee24) * ee12 -
   5 * ee1025)/ee2 + 4 * (ee494 * ee53)) * ee7/xi) * ee34 +
   ((ee25 - (((ee85 * (6 * ee102 - ee815) + 2 * (ee1092/ee9) +
   4 * (ee12 * ee349) - (4 * ((ee812 + ee351)/ee9) + 8 * ee521 -
   64 * ee469) * ee2/ee12)/ee12 + 2 * ee349)/ee12 + ee1068) * ee2)/ee12 +
   (ee25 - ((ee747 - (((ee467 + 2 * (ee12 * ee350) -
   ee523)/ee12 + ee522 + 2 * ee350)/ee12 + ee1130 * ee7/xi) * ee2)/ee12 +
   ee1068) * ee2) * ee7/xi)/ee24 + ((4 * (ee572 * ee53/ee24) -
   4 * (ee542 * ee60)) * ee2 - 6 * (ee190 * ee195))/ee12) * ee296/ee72;
out(j, 51) = 0.5 * (ee986 * ee197/ee72);
out(j, 52) = ee1083;
out(j, 53) = ee1084;
out(j, 54) = ee914;

}

return out;

}
