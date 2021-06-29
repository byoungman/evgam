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
double ee1, ee2, ee4, ee6, ee7, ee8, ee9;
double nllh=0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
lkappa1 = lkappa1vec[j];
lkappa2 = lkappa2vec[j];
logitp = logitpvec[j];

ee1 = 1 + xi * y/exp(lpsi);
ee2 = 1/xi;
ee4 = R_pow(ee1, 1 + ee2);
ee6 = 1 - 1/R_pow(ee1, ee2);
ee7 = 1 + exp(-logitp);
ee8 = exp(lkappa1);
ee9 = exp(lkappa2);

nllh -= log((1 - 1/ee7) * R_pow(ee6, ee9 - 1) * ee9/ee4 + R_pow(ee6, ee8 - 1) * ee8/(ee7 * ee4)) - lpsi;
    
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
double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee20, ee21, ee22, ee23, ee24, ee25, ee27, ee28;
double ee31, ee35, ee36, ee37, ee38, ee39;
double ee40, ee41, ee47, ee48, ee49;
double ee50, ee51, ee52, ee53, ee54, ee55, ee56, ee57, ee58, ee59;
double ee60, ee61, ee63, ee64, ee65, ee66, ee68, ee69;
double ee73, ee76, ee78, ee79;
double ee80, ee81, ee83, ee85, ee88, ee89;
double ee90, ee91, ee92, ee93, ee94, ee95, ee96, ee98;
double ee101, ee103, ee104, ee105, ee106, ee107, ee108;
double ee114, ee115, ee116, ee117, ee118, ee119;
double ee121, ee122, ee126, ee127;
double ee130, ee134, ee136, ee138, ee139;
double ee140, ee143;

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
ee7 = 1 + ee4;
ee8 = 1 - 1/ee6;
ee9 = exp(lkappa1);
ee10 = exp(lkappa2);
ee12 = exp(-logitp);
ee13 = 1 + ee12;
ee14 = R_pow(ee5, ee7);
ee15 = ee9 - 1;
ee16 = ee10 - 1;
ee17 = R_pow(ee8, ee15);
ee18 = R_pow(ee8, ee16);
ee19 = log1p(ee3);
ee20 = ee13 * ee14;
ee21 = 1 - 1/ee13;
ee22 = ee4 - 1;
ee23 = R_pow(ee5, ee22);
ee24 = 2 * ee7;
ee25 = ee17 * ee9;
ee27 = y * ee23/ee1;
ee28 = R_pow(ee20, 2);
ee31 = ee27 - ee6 * ee19/xi;
ee35 = ee21 * ee18 * ee10/ee14 + ee25/ee20;
ee36 = ee9 - 2;
ee37 = ee10 - 2;
ee38 = R_pow(ee5, ee24);
ee39 = R_pow(xi, 2);
ee40 = R_pow(ee8, ee36);
ee41 = R_pow(ee8, ee37);
ee47 = y * ee7 * ee6/ee1 - ee14 * ee19/ee39;
ee48 = ee40 * ee15;
ee49 = 2/xi;
ee50 = log(ee8);
ee51 = ee41 * ee16;
ee52 = 3/xi;
ee53 = 1 + ee52;
ee54 = R_pow(ee5, ee53);
ee55 = R_pow(ee13, 2);
ee56 = R_pow(ee5, ee49);
ee57 = ee18 * ee10;
ee58 = xi * ee54;
ee59 = ee55 * ee14;
ee60 = R_pow(ee5, ee4 - ee24);
ee61 = xi * ee17;
ee63 = xi * ee18 * ee7;
ee64 = ee35 * ee1;
ee65 = ee17 * ee13;
ee66 = ee17 * ee14;
ee68 = ee17 + ee25 * ee50;
ee69 = ee48 * ee31;
ee73 = ee18 + ee57 * ee50;
ee76 = ee51 * ee31/ee58 - ee18 * ee47/ee38;
ee78 = ee51/ee38 - ee63 * ee60;
ee79 = ee61 * ee7;
ee80 = xi * ee13;
ee81 = ee35 * ee13;
ee83 = (ee69/(ee80 * ee54) - ee65 * ee47/ee28) * ee9 + ee76 *  ee21 * ee10;
ee85 = (ee48/(ee13 * ee38) - ee79 * ee13 * ee6/ee28) * ee9 +  ee78 * ee21 * ee10;
ee88 = ee66 * ee9/ee28 - ee57/ee59;
ee89 = R_pow(ee5, ee4 - 2);
ee90 = ee5 * ee1;
ee91 = ee81 * ee14;
ee92 = ee35 * ee14;
ee93 = ee7 * ee31;
ee94 = ee23 * ee19;
ee95 = R_pow(ee31, 2);
ee96 = 1 + ee49;
ee98 = xi * ee7;
ee101 = y * ee89 * ee22/ee1;
ee103 = y/ee90 - 2 * (ee19/xi);
ee104 = ee17 * ee55;
ee105 = R_pow(ee8, ee9 - 3);
ee106 = R_pow(ee8, ee10 - 3);
ee107 = R_pow(ee5, ee96);
ee108 = 4 * ee7;
ee114 = ee40 * ee13 * ee47;
ee115 = ee48 * ee50;
ee116 = ee105 * ee36;
ee117 = ee51 * ee50;
ee118 = ee41 * ee47;
ee119 = ee106 * ee37;
ee121 = ee93 + ee6;
ee122 = R_pow(ee5, ee7 + ee24);
ee126 = ee6 + ee27;
ee127 = ee38 * ee1;
ee130 = R_pow(ee47, 2);
ee134 = (ee101 - ee94/ee39)/ee56 - 2 * (ee31/(xi * ee107));
ee136 = (y * (ee101 - (ee23 + ee94/xi)/xi)/ee1 - (ee6 * ee103 + ee19 * ee31/xi)/xi)/ee56 - 2 * (ee95/(xi * R_pow(ee5, ee52)));
ee138 = 2 * (y/(R_pow(ee5, ee4 + 2) * ee1)) - (ee23 + ee2 * ee89 * ee22/ee1)/ee56;
ee139 = ee98 * ee13;
ee140 = xi * R_pow(ee5, 4/xi);
ee143 = y * (ee93 - ee6/xi)/ee1 - (ee14 * ee103 + ee19 *  ee47)/xi;

out(j, 0) = 1 + y * ee85/ee64;
out(j, 1) = -(ee83/ee35);
out(j, 2) = -(ee68 * ee9/ee91);
out(j, 3) = -(ee73 * ee21 * ee10/ee92);
out(j, 4) = -(ee88 * ee12/ee35);
out(j, 5) = y * ((((ee40 * ee138 - y * ee105 * ee36/ee127)/(ee13 * ee6) +
   ee2 * ee40 * ee7 * ee13/(ee28 * ee1)) * ee15/ee5 -
   ee139 * (y * (2 * (ee79 * ee55 * ee54/ee28) - ee48/ee5)/ee1 -
   ee126 * ee17)/ee28) * ee9 + (((ee41 * ee138 - y * ee106 * ee37/ee127)/ee14 +
   ee2 * ee41 * ee7/(R_pow(ee5, 1 + ee24) * ee1)) * ee16 -
   ee98 * (2 * (ee2 * ee18 * ee7 * R_pow(ee5, ee53 -
   ee108)/ee1) - (ee126 * ee18 + y * ee41 * ee16/ee90)/ee38)) * ee21 * ee10 +
   y * R_pow(ee85, 2)/ee64)/ee64;
out(j, 6) = y * ((((ee134 * ee40 + ee116 * ee31/ee58)/ee20 -
   ee114/(ee28 * ee14)) * ee15 - ((ee69/ee6 - 2 * (ee61 * ee55 * ee107 * ee47/ee28)) * ee7 +
   ee121 * ee17) * ee13/ee28) * ee9 +
   (((ee134 * ee41 + ee119 * ee31/ee58)/ee14 - ee118/ee122) * ee16 +
   2 * (ee63 * R_pow(ee5, ee96 - ee108) * ee47) - (ee121 * ee18 +
   ee41 * ee7 * ee16 * ee31/ee6)/ee38) * ee21 * ee10 -
   ee83 * ee85/ee35)/ee64;
out(j, 7) = y * (((ee115/ee14 + ee40/ee14) * ee9 + ee48/ee14)/ee20 -
   (ee85/ee91 + ee139 * ee6/ee28) * ee68) * ee9/ee64;
out(j, 8) = y * (((ee117/ee14 + ee41/ee14) * ee10 + ee51/ee14)/ee14 -
   (ee85/ee92 + ee98 * ee60) * ee73) * ee21 * ee10/ee64;
out(j, 9) = y * ((ee48 + xi * (ee17 * ee6 - 2 * (ee104 * R_pow(ee5, ee4 +
   ee24)/ee28)) * ee7) * ee9/ee28 - (ee85 * ee88/ee35 +
   ee78 * ee10/ee55)) * ee12/ee64;
out(j, 10) =  - ((((((ee136 * ee41 + ee119 * ee95/ee140)/ee14 -
   ee118 * ee31/R_pow(ee5, ee24 + ee49)) * ee16 - (ee18 * ee143 +
   ee51 * ee47 * ee31/ee56)/ee38)/xi + 2 * (ee18 * R_pow(ee5, ee7 -
   ee108) * ee130)) * ee21 * ee10 + (((ee136 * ee40 +
   ee116 * ee95/ee140)/ee20 - ee114 * ee31/(ee28 * ee56)) * ee15/xi -
   ((ee17 * ee143 + ee48 * ee47 * ee31/ee56)/xi - 2 * (ee104 * ee14 * ee130/ee28)) * ee13/ee28) * ee9 -
   R_pow(ee83, 2)/ee35)/ee35);
out(j, 11) = -((((ee115/ee56 + ee40/ee56) * ee9 + ee48/ee56) * ee31/(ee80 * ee14) -
   (ee83/ee91 + ee13 * ee47/ee28) * ee68) * ee9/ee35);
out(j, 12) = -((((ee117/ee56 + ee41/ee56) * ee10 + ee51/ee56) * ee31/(xi * ee14) -
   (ee83/ee92 + ee47/ee38) * ee73) * ee21 * ee10/ee35);
out(j, 13) =  - ((((ee17 - 2 * (ee104 * ee38/ee28)) * ee47 +
   ee40 * R_pow(ee5, 1 - ee4) * ee15 * ee31/xi) * ee9/ee28 - (ee83 * ee88/ee35 +
   ee76 * ee10/ee55)) * ee12/ee35);
out(j, 14) =  - ((((ee68 + 2 * ee17) * ee9 * ee50 + ee17)/ee14 -
   R_pow(ee68, 2) * ee9/(ee81 * ee38)) * ee9/ee81);
out(j, 15) = ee68 * ee73 * ee21 * ee9 * ee10/(R_pow(ee35, 2) * ee13 * ee38);
out(j, 16) = -(ee68 * (ee14/ee28 - ee88/ee91) * ee12 * ee9/ee35);
out(j, 17) =  - ((((ee73 + 2 * ee18) * ee10 * ee50 + ee18)/ee14 -
   R_pow(ee73, 2) * ee21 * ee10/(ee35 * ee38)) * ee21 * ee10/ee35);
out(j, 18) = (ee88 * ee21/ee92 + 1/ee59) * ee73 * ee12 * ee10/ee35;
out(j, 19) = ((ee66 - 2 * (ee65 * ee122 * ee12/ee28)) * ee9/ee28 +
   R_pow(ee88, 2) * ee12/ee35 - ee18 * (1 - 2 * (ee12/ee13)) * ee10/ee59) * ee12/ee35;   
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
double ee1, ee2, ee3, ee4, ee6, ee7;
double nllh=0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
ldelta = ldeltavec[j];

ee1 = exp(ldelta);
ee2 = exp(lpsi);
ee3 = 1 + xi * y/ee2;
ee4 = 1/xi;
ee6 = 1 + ee1;
ee7 = 1/R_pow(ee3, ee4);

nllh -= log(-((R_pow(ee7, ee1 - 1)/R_pow(ee3, 1 + 2/xi) - (1 - R_pow(ee7, ee1)/ee6) * 
ee6/(R_pow(ee3, 1 + ee4) * ee1))/ee2));
    
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
double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
double ee30, ee32, ee33, ee34, ee35, ee36, ee37, ee38, ee39;
double ee40, ee49;
double ee53, ee54, ee55, ee56, ee58, ee59;
double ee60, ee61, ee62, ee63, ee65, ee66, ee67, ee68, ee69;
double ee70, ee71, ee74, ee75, ee78;
double ee80, ee81, ee82, ee83, ee84, ee85;
double ee92, ee93, ee94, ee97, ee98;
double ee100, ee101, ee103, ee104, ee105, ee109;
double ee113, ee116, ee117, ee118;
double ee122, ee124, ee126, ee127, ee128, ee129;
double ee130, ee131;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
ldelta = ldeltavec[j];


ee1 = exp(lpsi);
ee2 = xi * y;
ee3 = ee2/ee1;
ee4 = 1 + ee3;
ee5 = 1/xi;
ee6 = exp(ldelta);
ee7 = R_pow(ee4, ee5);
ee8 = 1/ee7;
ee9 = 2/xi;
ee10 = 1 + ee5;
ee11 = 1 + ee9;
ee12 = ee6 - 1;
ee13 = 1 + ee6;
ee14 = R_pow(ee4, ee10);
ee15 = log1p(ee3);
ee16 = R_pow(ee8, ee12);
ee17 = R_pow(ee8, ee6);
ee18 = ee14 * ee6;
ee19 = R_pow(ee4, ee11);
ee20 = ee17/ee13;
ee21 = 1 - ee20;
ee22 = ee5 - 1;
ee23 = R_pow(ee4, ee22);
ee24 = ee21 * ee13;
ee25 = R_pow(ee4, ee9);
ee26 = R_pow(ee18, 2);
ee27 = ee6 - 2;
ee28 = R_pow(ee8, ee27);
ee29 = R_pow(xi, 2);
ee30 = y * ee23;
ee32 = ee7 * ee15/xi;
ee33 = ee30/ee1;
ee34 = ee33 - ee32;
ee35 = ee24/ee18;
ee36 = ee16/ee19;
ee37 = 2 * ee11;
ee38 = ee36 - ee35;
ee39 = ee28 * ee12;
ee40 = 3/xi;
ee49 = y * ee10 * ee7/ee1 - ee14 * ee15/ee29;
ee53 = y * ee11 * ee25/ee1 - 2 * (ee19 * ee15/ee29);
ee54 = 2 * ee10;
ee55 = R_pow(ee4, 1 + ee40);
ee56 = R_pow(ee4, ee37);
ee58 = ee17 * ee15/xi;
ee59 = R_pow(ee6, 2);
ee60 = xi * ee19;
ee61 = ee58 + 1;
ee62 = 4/xi;
ee63 = R_pow(ee4, ee54);
ee65 = R_pow(ee38/ee1, 2) * R_pow(ee1, 2);
ee66 = ee24 * ee14;
ee67 = R_pow(ee4, ee5 - 2);
ee68 = R_pow(ee4, ee9 - ee37);
ee69 = ee16 * ee15;
ee70 = 1 + ee62;
ee71 = (ee66/ee26 - ee69/ee60) * ee6;
ee74 = (ee16/ee55 + ee39/R_pow(ee4, ee70)) * ee34/xi + ee16 * ee53/ee56;
ee75 = ee61/ee14;
ee78 = ee24 * ee6 * ee49/ee26;
ee80 = ee11 * ee68 * ee16;
ee81 = ee4 * ee1;
ee82 = ee16 * ee6;
ee83 = ee15/xi;
ee84 = ee71 - ee75;
ee85 = ee74 - ee78;
ee92 = ee35 + y * (ee16/ee63 + ee39/R_pow(ee4, 2 + ee40) + xi * (ee80 - ee21 * ee10 * ee13 * ee7 * ee6/ee26))/ee1 - ee36;
ee93 = ee23 * ee15;
ee94 = R_pow(ee34, 2);
ee97 = y * ee67 * ee22/ee1;
ee98 = y/ee81;
ee100 = ee82 * ee15;
ee101 = R_pow(ee8, ee6 - 3);
ee103 = 2 * (y * R_pow(ee4, ee9 - 1)/ee1);
ee104 = 4 * ee11;
ee105 = ee98 - 2 * ee83;
ee109 = ee11 * (ee103 - 2 * (ee25 * ee15/xi));
ee113 = ee63 * ee1;
ee116 = ee39 * ee15;
ee117 = ee28 * ee53;
ee118 = ee101 * ee27;
ee122 = (ee97 - ee93/ee29)/ee25 - 2 * (ee34/ee60);
ee124 = (y * (ee97 - (ee23 + ee93/xi)/xi)/ee1 - (ee7 * ee105 + ee15 * ee34/xi)/xi)/ee25 - 2 * (ee94/(xi * R_pow(ee4, ee40)));
ee126 = 2 * (y/(R_pow(ee4, ee5 + 2) * ee1)) - (ee23 + ee2 * ee67 * ee22/ee1)/ee25;
ee127 = xi * ee14;
ee128 = xi * ee55;
ee129 = xi * ee25;
ee130 = xi * R_pow(ee4, ee62);
ee131 = ee2 * ee11;

out(j, 0) = -(ee92/ee38);
out(j, 1) = ee85/ee38;
out(j, 2) = -(ee84/ee38);
out(j, 3) =  - (((y * (((ee16 * ee126 + y * ee28 * ee12/ee113)/ee7 +
   ee2 * ee10 * ee16 * ee59/(ee26 * ee1) - (2 * (ee16/R_pow(ee4, ee54 -
   1)) + 2 * (ee39/ee55)))/ee4 + ((ee28 * ee126 +
   y * ee101 * ee27/ee113)/ee19 + ee131 * R_pow(ee4, ee5 -
   (1 + ee37)) * ee28/ee1) * ee12 + xi * (((ee7 + 2 * ee7 - y * (2 * (xi * ee10 * ee55 * ee59/ee26) -
   ee23)/ee1) * ee21 * ee13 +
   y * ee16 * ee6/ee81) * ee10 * ee6/ee26 + ((ee30 * ee28 * ee12/ee1 -
   (ee25 + ee103) * ee16)/ee56 + 2 * (ee131 * R_pow(ee4, 1 +
   6/xi - ee104) * ee16/ee1) - 2 * (ee68 * ee16)) * ee11)) -
   ee38 * ee1)/ee1 + 2 * ee38)/ee38 - R_pow(ee92, 2)/ee65);
out(j, 4) =  - (ee85 * ee92/ee65 + (ee74 + y * (((ee122 * ee28 -
   ee118 * ee34/ee128)/ee19 - ee117/R_pow(ee4, ee10 + ee37)) * ee12 +
   ((ee109 + ee25) * ee16 - ee11 * ee28 * ee12 * ee34)/ee56 +
   (ee122 * ee16 - ee39 * ee34/ee128)/ee14 - (((ee10 * ee16 * ee34/ee7 +
   ee16 * ee49/ee14) * ee6 + (ee10 * (ee33 -
   (ee32 + 2 * (ee60 * ee59 * ee49/ee26))) + ee7) * ee21 * ee13) * ee6/ee26 +
   2 * (xi * ee11 * R_pow(ee4, ee70 - ee104) * ee16 * ee53)))/ee1 -
   ee78)/ee38);
out(j, 5) =  - ((ee75 + y * (((ee28/ee14 - ee116/ee127)/ee19 -
   ((ee82 + xi * ((ee7 - 2 * (R_pow(ee4, ee5 + ee54) * ee59/ee26)) * ee21 * ee13 +
   ee61 * ee7 * ee6) * ee10)/ee26 + ee80 * ee15)) * ee6 +
   (ee16/ee14 - ee100/ee127)/ee14)/ee1 - ee71)/ee38 -
   ee84 * ee92/ee65);
out(j, 6) = ((((ee124 * ee28 - ee118 * ee94/ee130)/ee19 - ee117 * ee34/R_pow(ee4, ee37 +
   ee9)) * ee12 + (ee124 * ee16 -
   ee39 * ee94/ee130)/ee14 + (ee16 * (y * (ee109 - 2 * (ee25/xi))/ee1 -
   (ee19 * (2 * ee98 - 4 * ee83) + 2 * (ee15 * ee53))/xi) -
   ee39 * ee53 * ee34/ee25)/ee56 - ee16 * ee59 * ee49 * ee34/(ee26 * ee25))/xi -
   (((ee24 * (y * (ee10 * ee34 - ee7/xi)/ee1 -
   (ee14 * ee105 + ee15 * ee49)/xi) + ee82 * ee49 * ee34/ee25)/xi -
   2 * (ee66 * ee59 * R_pow(ee49, 2)/ee26)) * ee6/ee26 +
   2 * (R_pow(ee4, ee11 - ee104) * ee16 * R_pow(ee53, 2))))/ee38 +
   R_pow(ee85, 2)/ee65;
out(j, 7) = ((((ee28/ee25 - ee116/ee129) * ee34/ee19 - ee69 * ee53/ee56)/xi -
   (ee61 * ee6 + ee21 * (1 - 2 * (ee63 * ee59/ee26)) * ee13) * ee49/ee26) * ee6 +
   ((ee16/ee25 - ee100/ee129)/ee14 -
   R_pow(ee4, 1 - ee5) * ee16 * ee59/ee26) * ee34/xi)/ee38 -
   ee84 * ee85/ee65;
out(j, 8) =  - (((((ee14 - 2 * (R_pow(ee4, ee10 + ee54) * ee59/ee26)) * ee21 * ee13 +
   2 * (ee61 * ee14 * ee6))/ee26 + (ee100/xi -
   ee16) * ee15/ee60) * ee6 - (((ee17 - 2 * ee17) * ee6/ee13 +
   ee17)/ee13 + ((ee58 + ee20) * ee6 - ee17) * (1/ee13 -
   ee83) + 1)/ee14)/ee38 - R_pow(ee84, 2)/ee65);
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
arma::mat out = arma::mat(nobs, 25, arma::fill::zeros);

if (dcate == 1) {
  lpsivec = lpsivec.elem(dupid);
  xivec = xivec.elem(dupid);
  ldeltavec = ldeltavec.elem(dupid);
}

double y, lpsi, xi, ldelta;
double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee20, ee21, ee22, ee24, ee25, ee26, ee27, ee28, ee29;
double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee37, ee38, ee39;
double ee40, ee41;
double ee50, ee54, ee55, ee58, ee59;
double ee60, ee61, ee62, ee63, ee64, ee66, ee67, ee68, ee69;
double ee71, ee72, ee73, ee74, ee75, ee76, ee77, ee78, ee79;
double ee80, ee81, ee84, ee88, ee89;
double ee90, ee91, ee92, ee93, ee97, ee98, ee99;
double ee100, ee101, ee102, ee103, ee104, ee105, ee108;
double ee113, ee114, ee115, ee116, ee119;
double ee120, ee121, ee122, ee123, ee124, ee125, ee126, ee127, ee128, ee129;
double ee130, ee131, ee133, ee134, ee136, ee139;
double ee141, ee142, ee144, ee146, ee149;
double ee150;
double ee161, ee162, ee165, ee166, ee167, ee168, ee169;
double ee170, ee171, ee173, ee174, ee175, ee176, ee177, ee178, ee179;
double ee180, ee181, ee182, ee183, ee184, ee185, ee186, ee187, ee189;
double ee190, ee191, ee193, ee194, ee195, ee197, ee198, ee199;
double ee200, ee202, ee203, ee205, ee207, ee208, ee209;
double ee210, ee211, ee212, ee216, ee217, ee218, ee219;
double ee220, ee221, ee223, ee224, ee225, ee226, ee227, ee229;
double ee231, ee232, ee233, ee234, ee235, ee238, ee239;
double ee240, ee241, ee243, ee244, ee245, ee246, ee248, ee249;
double ee250, ee251, ee254, ee256, ee257, ee258, ee259;
double ee260, ee261, ee265, ee266, ee267, ee268;
double ee271, ee272, ee273, ee274, ee275, ee276, ee277, ee278, ee279;
double ee280, ee289;
double ee290, ee291, ee293, ee294, ee295, ee296, ee297, ee298;
double ee300, ee301, ee302, ee303, ee304, ee305, ee308, ee309;
double ee310, ee313, ee315, ee316, ee317, ee318;
double ee324, ee325;
double ee330, ee331, ee334, ee335, ee336, ee337;
double ee340, ee341, ee342, ee343, ee344, ee345, ee346, ee347, ee348;
double ee350, ee351, ee353, ee354, ee355, ee357, ee358, ee359;
double ee360, ee364;
double ee372, ee373, ee374, ee375, ee377;
double ee381, ee389;
double ee391, ee393, ee394, ee395, ee396, ee397;
double ee402, ee403, ee404, ee407, ee408;
double ee411, ee412, ee413, ee417, ee418, ee419;
double ee420, ee422, ee423, ee424, ee429;
double ee431, ee432, ee433, ee434, ee437, ee438, ee439;
double ee440, ee441, ee442, ee443, ee444;

for (int j=0; j < nobs; j++) {

y = yvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
ldelta = ldeltavec[j];

ee1 = exp(lpsi);
ee2 = xi * y;
ee3 = ee2/ee1;
ee4 = 1 + ee3;
ee5 = 1/xi;
ee6 = exp(ldelta);
ee7 = R_pow(ee4, ee5);
ee8 = 2/xi;
ee9 = 1/ee7;
ee10 = log1p(ee3);
ee11 = 1 + ee5;
ee12 = ee5 - 1;
ee13 = 1 + ee8;
ee14 = ee6 - 1;
ee15 = R_pow(ee4, ee11);
ee16 = R_pow(ee4, ee12);
ee17 = 1 + ee6;
ee18 = R_pow(ee9, ee14);
ee19 = R_pow(ee4, ee8);
ee20 = R_pow(ee9, ee6);
ee21 = y * ee16;
ee22 = ee21/ee1;
ee24 = ee7 * ee10/xi;
ee25 = ee22 - ee24;
ee26 = R_pow(ee4, ee13);
ee27 = ee6 - 2;
ee28 = ee15 * ee6;
ee29 = R_pow(xi, 2);
ee30 = R_pow(ee9, ee27);
ee31 = ee20/ee17;
ee32 = 1 - ee31;
ee33 = R_pow(ee28, 2);
ee34 = ee5 - 2;
ee35 = ee32 * ee17;
ee36 = R_pow(ee4, ee34);
ee37 = 3/xi;
ee38 = ee4 * ee1;
ee39 = 2 * ee13;
ee40 = ee10/xi;
ee41 = y/ee38;
ee50 = y * ee11 * ee7/ee1 - ee15 * ee10/ee29;
ee54 = y * ee13 * ee19/ee1 - 2 * (ee26 * ee10/ee29);
ee55 = ee30 * ee14;
ee58 = y * ee36 * ee12/ee1;
ee59 = ee16 * ee10;
ee60 = 2 * ee11;
ee61 = R_pow(ee6, 2);
ee62 = ee35/ee28;
ee63 = ee18/ee26;
ee64 = 4/xi;
ee66 = ee8 - 1;
ee67 = ee41 - 2 * ee40;
ee68 = R_pow(ee4, ee66);
ee69 = R_pow(ee25, 2);
ee71 = ee20 * ee10/xi;
ee72 = R_pow(ee4, ee39);
ee73 = ee6 - 3;
ee74 = R_pow(ee4, 1 + ee37);
ee75 = ee63 - ee62;
ee76 = xi * ee26;
ee77 = R_pow(ee4, ee60);
ee78 = R_pow(ee9, ee73);
ee79 = 2 * (y * ee68/ee1);
ee80 = ee71 + 1;
ee81 = ee59/xi;
ee84 = ee18 * ee6;
ee88 = y * (ee58 - (ee16 + ee81)/xi)/ee1 - (ee7 * ee67 +  ee10 * ee25/xi)/xi;
ee89 = ee59/ee29;
ee90 = R_pow(ee4, ee64);
ee91 = ee58 - ee89;
ee92 = ee16 + ee2 * ee36 * ee12/ee1;
ee93 = R_pow(ee4, ee37);
ee97 = ee79 - 2 * (ee19 * ee10/xi);
ee98 = 4 * ee13;
ee99 = xi * ee93;
ee100 = ee13 * ee97;
ee101 = R_pow(ee4, ee8 - ee39);
ee102 = ee84 * ee10;
ee103 = xi * ee19;
ee104 = ee18 * ee10;
ee105 = ee78 * ee27;
ee108 = ee88/ee19 - 2 * (ee69/ee99);
ee113 = 2 * (y/(R_pow(ee4, ee5 + 2) * ee1)) - ee92/ee19;
ee114 = ee35 * ee15;
ee115 = 1 + ee64;
ee116 = xi * ee90;
ee119 = ee91/ee19 - 2 * (ee25/ee76);
ee120 = ee11 * ee25;
ee121 = 2 * ee41;
ee122 = R_pow(ee75/ee1, 2);
ee123 = ee55 * ee10;
ee124 = ee77 * ee1;
ee125 = ee18/ee15;
ee126 = xi * ee15;
ee127 = xi * ee74;
ee128 = R_pow(ee4, ee115);
ee129 = ee122 * R_pow(ee1, 2);
ee130 = (ee114/ee33 - ee104/ee76) * ee6;
ee131 = ee80/ee15;
ee133 = ee13 * ee101 * ee18;
ee134 = ee121 - 4 * ee40;
ee136 = (ee18/ee74 + ee55/ee128) * ee25/xi;
ee139 = ee35 * ee6 * ee50/ee33;
ee141 = ee18 * ee54/ee72;
ee142 = ee18/ee19;
ee144 = ee136 + ee141;
ee146 = R_pow(ee4, 2 + ee37);
ee149 = y * (ee120 - ee7/xi)/ee1 - (ee15 * ee67 + ee10 *  ee50)/xi;
ee150 = ee130 - ee131;
ee161 = ee18/ee77 + ee55/ee146 + xi * (ee133 - ee32 * ee11 *  ee17 * ee7 * ee6/ee33);
ee162 = 2 * ee20;
ee165 = y * (ee100 - 2 * (ee19/xi))/ee1 - (ee26 * ee134 +  2 * (ee10 * ee54))/xi;
ee166 = ee144 - ee139;
ee167 = ee100 + ee19;
ee168 = ee19 + ee79;
ee169 = ee102/xi;
ee170 = y * ee161;
ee171 = ee108 * ee30;
ee173 = ee62 + ee170/ee1 - ee63;
ee174 = ee30 * ee113;
ee175 = ee119 * ee30;
ee176 = ee171 - ee105 * ee69/ee116;
ee177 = ee2 * ee13;
ee178 = ee174 + y * ee78 * ee27/ee124;
ee179 = ee30/ee19;
ee180 = R_pow(ee50, 2);
ee181 = 1/ee17;
ee182 = ee168 * ee18;
ee183 = ee80 * ee6;
ee184 = ee175 - ee105 * ee25/ee127;
ee185 = ee108 * ee18;
ee186 = R_pow(ee4, ee13 - ee98);
ee187 = ee169 - ee18;
ee189 = ee30/ee15;
ee190 = ee179 - ee123/ee103;
ee191 = xi * ee11;
ee193 = ((ee20 - ee162) * ee6/ee17 + ee20)/ee17 + ((ee71 +  ee31) * ee6 - ee20) * (ee181 - ee40);
ee194 = ee167 * ee18;
ee195 = ee80 * ee7;
ee197 = ee35 * ee149 + ee84 * ee50 * ee25/ee19;
ee198 = R_pow(ee4, ee11 + ee60);
ee199 = R_pow(ee4, ee5 - 3);
ee200 = R_pow(ee4, ee60 - 1);
ee202 = ee18 * ee113 + y * ee30 * ee14/ee124;
ee203 = ee55 * ee25;
ee205 = ee55 * ee69/ee116;
ee207 = ee30 * ee54;
ee208 = ee189 - ee123/ee126;
ee209 = R_pow(ee54, 2);
ee210 = 2 * ee7;
ee211 = 6/xi;
ee212 = xi * ee13;
ee216 = ee21 * ee30 * ee14/ee1 - ee182;
ee217 = ee193 + 1;
ee218 = ee194 - ee13 * ee30 * ee14 * ee25;
ee219 = ee119 * ee18;
ee220 = ee185 - ee205;
ee221 = ee120 + ee7;
ee223 = ee198 * ee61/ee33;
ee224 = R_pow(ee4, 1 + ee211 - ee98);
ee225 = ee7 + ee210;
ee226 = ee7 + ee22;
ee227 = ee101 * ee18;
ee229 = ee18 * ee165 - ee55 * ee54 * ee25/ee19;
ee231 = ee102/ee103;
ee232 = ee125 - ee102/ee126;
ee233 = ee18/ee200;
ee234 = ee203/ee127;
ee235 = 2 * ee125;
ee238 = y * ee18 * ee6/ee38;
ee239 = ee33 * ee1;
ee240 = ee187 * ee10;
ee241 = ee75 * ee1;
ee243 = ee80 * ee15 * ee6;
ee244 = ee195 * ee6;
ee245 = ee219 - ee234;
ee246 = ee11 * ee18;
ee248 = R_pow(ee4, 1 - ee5) * ee18;
ee249 = ee186 * ee18;
ee250 = R_pow(ee4, ee115 - ee98);
ee251 = R_pow(ee4, ee5 - (1 + ee39));
ee254 = ee7 * ee88;
ee256 = ee18 * ee61 * ee50;
ee257 = ee142 - ee231;
ee258 = 2 * ee223;
ee259 = 2 * ee142;
ee260 = 2 * ee16;
ee261 = ee2 * ee11;
ee265 = y * ee199 * ee34/ee1 - ee36 * ee10/ee29;
ee266 = ((ee176/ee26 - ee207 * ee25/R_pow(ee4, ee39 + ee8)) * ee14 + ee220/ee15 + ee229/ee72 - ee256 * ee25/(ee33 * ee19))/xi;
ee267 = (((ee15 - ee258) * ee32 * ee17 + 2 * ee243)/ee33 +  ee240/ee76) * ee6;
ee268 = ee217/ee15;
ee271 = (ee197/xi - 2 * (ee114 * ee61 * ee180/ee33)) * ee6/ee33 +  2 * (ee249 * ee209);
ee272 = ee190 * ee25;
ee273 = ee183 + ee35;
ee274 = ee248 * ee61;
ee275 = R_pow(ee4, ee5 + ee60);
ee276 = ee77 * ee61;
ee277 = ee55/ee74;
ee278 = 2 * ee227;
ee279 = ((ee272/ee26 - ee104 * ee54/ee72)/xi - (ee183 +  ee32 * (1 - 2 * (ee276/ee33)) * ee17) * ee50/ee33) *  ee6;
ee280 = ee218/ee72;
ee289 = (ee257/ee15 - ee274/ee33) * ee25/xi;
ee290 = ee245/ee15;
ee291 = ee184/ee26;
ee293 = (ee7 - 2 * (ee275 * ee61/ee33)) * ee32 * ee17;
ee294 = ee232/ee15;
ee295 = ee208/ee26;
ee296 = ee35 * ee7;
ee297 = R_pow(ee4, ee11 + ee39);
ee298 = ee16 * ee67;
ee300 = ee254 + ee69/xi;
ee301 = R_pow(ee4, ee8 - 2);
ee302 = R_pow(ee4, ee98);
ee303 = ee30 * ee6;
ee304 = R_pow(ee9, ee6 - 4);
ee305 = ee20 * ee6;
ee308 = (y * ((ee202/ee7 + ee261 * ee18 * ee61/ee239 - (2 *  ee233 + 2 * ee277))/ee4 + (ee178/ee26 + ee177 * ee251 *  ee30/ee1) * ee14 + xi * (((ee225 - y * (2 * (ee191 *  ee74 * ee61/ee33) - ee16)/ee1) * ee32 * ee17 + ee238) *  ee11 * ee6/ee33 + (ee216/ee72 + 2 * (ee177 * ee224 *  ee18/ee1) - ee278) * ee13)) - ee241)/ee1;
ee309 = 2 * ee18;
ee310 = 6 * ee13;
ee313 = ee76 * ee61 * ee50/ee33;
ee315 = y * (ee12 * ee265 - ee36/ee29)/ee1;
ee316 = ee266 - ee271;
ee317 = ee267 - ee268;
ee318 = ee279 + ee289;
ee324 = ee150 * ee173 * ee75;
ee325 = ee122 * ee1;
ee330 = ee221 * ee32 * ee17 + ee246 * ee6 * ee25/ee7;
ee331 = ee167 * ee30;
ee334 = ee226 * ee32 * ee17 + ee238;
ee335 = ee244 + ee296;
ee336 = ee119 * ee78;
ee337 = ee108 * ee78;
ee340 = ee18 - ee169;
ee341 = ee18 * ee149;
ee342 = ee18 + ee309;
ee343 = ee18/ee7;
ee344 = ee30 - ee123/xi;
ee345 = ee30 * ee165;
ee346 = ee78 * ee113;
ee347 = ee105 * ee10;
ee348 = ee304 * ee73;
ee350 = ee305 * ee10/xi;
ee351 = ee20 + ee162;
ee353 = (ee121 - 6 * ee40)/xi + y * (ee8 + ee41)/ee38;
ee354 = ee308 + 2 * ee75;
ee355 = 1 + ee60;
ee357 = 2 * ee71;
ee358 = 4 * ee11;
ee359 = y * ((ee291 - ee207/ee297) * ee14 + ee280 + ee290 -  (((ee246 * ee25/ee7 + ee18 * ee50/ee15) * ee6 + (ee11 *  (ee22 - (ee24 + 2 * ee313)) + ee7) * ee32 * ee17) *  ee6/ee33 + 2 * (ee212 * ee250 * ee18 * ee54)));
ee360 = y * ((ee295 - ((ee84 + xi * (ee293 + ee244) * ee11)/ee33 +  ee133 * ee10)) * ee6 + ee294);
ee364 = y * ee301 * ee66/ee1 - 2 * (ee68 * ee10/ee29);
ee372 = ee166 * ee173 * ee75;
ee373 = ee176 * ee54;
ee374 = ee197 * ee15;
ee375 = ee221 * ee18;
ee377 = ee331/ee15 + xi * ee184 * ee13 * ee19;
ee381 = (ee16 + ee2 * (ee36 + 2 * ee36 + ee2 * ee199 * ee34/ee1) * ee12/ee1)/ee19 + y * (8 * (y/(R_pow(ee4, ee5 + 3) * ee1)) - ((2 * ((ee92 * ee7 + y * R_pow(ee4, 2 * ee12)/ee1) * ee16) + 2 * (ee92 * ee68))/ee90 + 2 * (ee92/ee26)))/ee1;
ee389 = ee168 * ee30/ee15;
ee391 = ee232 * ee7 + ee18/ee4;
ee393 = ee144 + ee359/ee1 - ee139;
ee394 = ee178 * ee19;
ee395 = (ee303 * ee10/xi - ee30) * ee14;
ee396 = ee190 * ee6;
ee397 = (ee78/ee19 - ee347/ee103) * ee14;
ee402 = ee337 - ee348 * ee69/ee116;
ee403 = ee11 * ee88;
ee404 = ee13 * (y * (2 * ee364 - 2 * (ee68/xi))/ee1 - (ee19 *  ee134 + 2 * (ee97 * ee10/xi))/xi);
ee407 = ee15 * ee149/xi + ee180;
ee408 = ee15 + 2 * ee15;
ee411 = ee26 * ee165/xi + ee209;
ee412 = ee224 * ee18;
ee413 = ee251 * ee30;
ee417 = ee341 * ee25/ee103;
ee418 = ee104/ee126;
ee419 = ee104/ee103;
ee420 = ee125 + ee235;
ee422 = ee345 * ee25/ee103;
ee423 = ee30/ee74;
ee424 = ee20 - ee350;
ee429 = (ee315 - (ee298 + ee10 * ee91)/ee29)/ee19 - ((2 *  (ee300 * ee16) + 2 * (ee7 * ee25 * ee91))/ee90 + (2 *  (ee91/ee93) - 8 * (ee25/ee127)) * ee25)/xi;
ee431 = (y * (ee315 - ((ee298 + ee10 * (ee58 - (ee81 + ee260)/xi) -  ee260)/xi + 2 * ee58)/xi)/ee1 - ((2 * (ee25 * ee67) +  ee10 * ee88)/xi - ee353 * ee7)/xi)/ee19 - ((2 * ee300 +  2 * ee254)/ee90 + 2 * (ee88/ee93) - 8 * (ee69/ee116)) *  ee25/xi;
ee432 = ee181 + ee40;
ee433 = 2 * (ee150 * ee166 * ee75/ee129);
ee434 = 2 * (ee324/ee325);
ee437 = 2 * (ee92 * ee25/ee99) + y * ((2 * ((ee16 * ee25/xi + ee7 * ee91) * ee16) + 2 * (ee68 * ee91))/ee90 - 8 * (ee25/(xi * R_pow(ee4, 2 + ee8))))/ee1 - (y * ((ee36 + xi * ee265) * ee12 - ee36)/ee1 - ee89)/ee19;
ee438 = 2 * ee195;
ee439 = 2 * ee31;
ee440 = 2 * ee68;
ee441 = ee162 + 4 * ee20;
ee442 = 2 + ee64;
ee443 = 4 * ee74;
ee444 = 5/xi;

out(j, 0) =  - (((3 * ee241 + 4 * (ee241 + ee170) + y * (((ee381 * ee18 +
   y * (ee178/ee15 + 2 * (ee174/ee15)) * ee14/ee1)/ee7 -
   (ee202 * (ee9 + 2/ee7) + (ee423 + 2 * ee423) * ee14 +
   ee233 + ee233 + ee233 + ee2 * ee342 * ee11 * ee61/ee239))/ee4 +
   ((ee381 * ee30 + y * ((ee346 + y * ee304 * ee73/ee124)/ee15 +
   2 * (ee346/ee15)) * ee27/ee1)/ee26 + ee2 * ((ee394 -
   ee389)/ee72 + ee178 * ee101 + 2 * (ee177 * R_pow(ee4, ee444 -
   ee98) * ee30/ee1) - (ee413 + 2 * ee413)) * ee13/ee1 - ee178 * (1/ee26 +
   2/ee26)) * ee14 + 2 * ee161 + xi * ((((ee19 +
   y * (ee440 + 2 * (ee2 * ee301 * ee66/ee1) + 4 * ee68)/ee1) * ee18 +
   y * (ee394 - 2 * ee389) * ee14/ee1)/ee72 + ee2 * ((2 * (ee128 * ee216) -
   2 * ((ee168 * ee26 + ee177 * ee90/ee1) * ee19 * ee18))/ee302 +
   2 * (ee250 * ee216) + 8 * (ee177 * R_pow(ee4, 10/xi +
   2 - ee310) * ee18/ee1) - (2 * ee412 + 4 * ee412)) * ee13/ee1 -
   (ee227 + (1/ee72 + 2/ee72) * ee216 + ee278)) * ee13 +
   ((ee225 + y * (xi * (ee11 * (2 * ee74 + ee443) * ee61/ee33 -
   ee58) - (ee260 + 4 * ee16))/ee1 - 4 * ee7) * ee32 * ee17 +
   y * (ee225 * ee202 + ee191 * (2 * (y * ee19 * ee18 * ee6/ee1) +
   4 * (ee334 * ee26) - ee35 * (8 * (ee261 * R_pow(ee4, 2 +
   ee444) * ee61/ee239) - 2 * ((ee226 * ee15 +
   ee261 * ee19/ee1) * ee7))) * ee6/ee33 - (ee226 * ee420 + ee342/ee4)) * ee6/ee1) * ee11 * ee6/ee33)))/ee1 -
   8 * ee75)/ee75 -
   (ee75 * (2 - 2 * (R_pow(ee173, 2)/ee129)) + ee308 + 2 * ee354) * ee173/ee129);
out(j, 1) =  - (((ee166 * ee1 + y * (((((ee226 * ee18/ee19 +
   ee343 + ee343) * ee25 - xi * (2 * (ee334 * ee15 * ee50) + y * (ee35 * (2 * ((ee221 * ee15 +
   ee191 * ee7 * ee50) * ee7) -
   8 * (ee191 * R_pow(ee4, ee442) * ee61 * ee50/ee33)) + 2 * (ee330 * ee26))/ee1) * ee6/ee33) * ee11 +
   (ee235 - ee202) * ee50 +
   y * (2 * (ee375/ee15) + ee191 * (2 * (ee245 * ee7) -
   2 * (ee7 * ee18 * ee61 * ee50/ee33)))/ee1) * ee6 + (ee11 * (y * (ee16 +
   ee260)/ee1 - (ee225 * ee10/xi + 4 * ee313)) + ee7 +
   ee7 + ee7 + y * (ee16 + ee191 * ee91)/ee1) * ee32 * ee17) * ee6/ee33 +
   ((ee168 * ee13 * ee30 * ee25/ee19 + y * ee377/ee1) * ee14 -
   (ee167 + y * (ee440 + 2 * (ee212 * ee364))/ee1) * ee18)/ee72 +
   ((ee30 * ee437 + ee27 * (y * ((ee336 - ee348 * ee25/ee127)/ee15 +
   ee336/ee15)/ee1 - ee346 * ee25/ee103))/ee26 +
   (2 * (ee30/ee297) - ee178/ee72) * ee54 + y * (ee377/ee72 -
   2 * (ee212 * R_pow(ee4, ee37 - ee98) * ee30 * ee54))/ee1 -
   2 * ee291) * ee14 + (ee18 * ee437 + ee14 * (y * (ee184/ee15 +
   ee175/ee15)/ee1 - ee174 * ee25/ee103))/ee15 + xi * ((4 * (ee250 * ee18) -
   2 * (ee186 * ee216)) * ee54 + y * ((2 * ((ee167 * ee26 +
   ee212 * ee19 * ee54) * ee19 * ee18) + 2 * (ee218 * ee128))/ee302 -
   8 * (ee212 * R_pow(ee4, 2 + 8/xi -
   ee310) * ee18 * ee54))/ee1) * ee13 - (2 * ee280 + 2 * ee290)))/ee1 -
   2 * ee166)/ee75 + (ee166 * ee354 - ee173 * (2 * (ee372/ee129) +
   2 * ee393))/ee129);
out(j, 2) =  - (((y * ((((ee182 * ee10 + ee2 * ee208 * ee19/ee1)/ee72 +
   2 * (ee227 * ee10) + ee2 * (ee208 * ee101 - 2 * (ee13 * ee224 * ee18 * ee10))/ee1) * ee13 +
   (ee344 * ee113 +
   y * ((ee78/ee15 - ee347/ee126) * ee14/ee15 + ee105/ee77)/ee1)/ee26 +
   ((ee309 + ee2 * (ee391 - 2 * (ee200 * ee18 * ee61/ee33)) * ee11/ee1 -
   ee202 * ee15) * ee6 + xi * (ee273 * ee226 +
   (ee438 + y * ee391/ee1 - (2 * (ee334 * ee77) + ee2 * (ee35 * (ee443 -
   8 * (R_pow(ee4, ee355 + ee37) * ee61/ee33)) +
   2 * (ee335 * ee26)) * ee11/ee1) * ee6/ee33) * ee6 + 2 * ee293) * ee11)/ee33 -
   2 * ee295) * ee6 + (ee340 * ee113 + y * (ee208 * ee6/ee15 +
   ee55/ee77)/ee1)/ee15 - 2 * ee294) - ee150 * ee1)/ee1 +
   2 * ee150)/ee75 - (ee150 * ee354 + ee173 * (2 * (ee131 +
   ee360/ee1 - ee130) - 2 * (ee324/ee129)))/ee129);
out(j, 3) =  - ((ee316 * ee173 + ee166 * ((2 * (ee372/ee325) +
   2 * ee359)/ee1 + 2 * ee136 + 2 * ee141 - 2 * ee139))/ee129 +
   (ee266 + y * (((ee429 * ee30 - (ee402/ee15 + 2 * (ee336 * ee25/ee19)) * ee27/xi)/ee26 +
   (2 * (R_pow(ee4, ee5 - ee98) * ee30 * ee54) -
   ee184/ee72) * ee54 - (ee184 * ee54 + ee345/ee126)/ee72) * ee14 +
   ((ee404 + 2 * (ee97/xi)) * ee18 - (ee176 * ee13 * ee19 +
   2 * (ee331 * ee25/ee103)) * ee14)/ee72 +
   (ee429 * ee18 - (ee176/ee15 + 2 * (ee175 * ee25/ee19)) * ee14/xi)/ee15 +
   (8 * (ee212 * R_pow(ee4, 2 + ee211 - ee310) * ee18 * ee54) -
   2 * (ee218 * ee186)) * ee54 - (((ee220 * ee11 * ee7 +
   (ee341/ee15 + 2 * (ee375 * ee25/ee19))/xi + (2 * ee219 -
   (2 * (ee256/ee33) + 2 * ee234)) * ee50 - (4 * (ee330 * ee15 * ee50) +
   xi * ee32 * ee11 * ee17 * (2 * (ee407 * ee7) -
   8 * (ee146 * ee61 * ee180/ee33))) * ee6/ee33) * ee6 + (ee403 +
   2 * (ee25/xi)) * ee32 * ee17) * ee6/ee33 + (2 * (ee218 * ee26 * ee54) +
   2 * (xi * ee411 * ee13 * ee19 * ee18))/ee302))/ee1 -
   ee271)/ee75);
out(j, 4) =  - ((ee318 * ee173 + ee166 * (ee131 + (ee360 - ee434)/ee1 -
   ee130) - ee150 * ee393)/ee129 + (ee318 + y * (((ee344 * ee119 -
   (ee397/ee15 + ee105/ee74) * ee25/xi)/ee26 +
   (2 * (ee13 * ee250 * ee18 * ee10) - ee208/ee72) * ee54 - (((((ee257 * ee7 +
   ee343) * ee25 - xi * (ee35 * (4 * ee26 - 8 * (R_pow(ee4, ee355 +
   ee8) * ee61/ee33)) + 2 * (ee335 * ee15)) * ee6 * ee50/ee33) * ee11 +
   ee245 * ee15 + (ee235 - (ee418 +
   2 * (R_pow(ee4, ee60 - ee11) * ee18 * ee6/ee33)) * ee6) * ee50 -
   2 * (ee330 * ee77 * ee6/ee33)) * ee6 + ee273 * ee221)/ee33 +
   (ee194 * ee10/xi + ee190 * ee13 * ee19 * ee25)/ee72)) * ee6 +
   (ee340 * ee119 - (ee396/ee15 + ee277) * ee25/xi)/ee15)/ee1)/ee75);
out(j, 5) =  - ((ee268 + y * (((((ee432 * (ee235 - ee125) -
   (ee418 + ee18/(ee17 * ee15))) * ee6 + ee125 + ee235 - ee420)/ee17 +
   (ee187/ee15 - ee235) * ee10/xi) * ee6 + ee125)/ee15 +
   (((ee395/ee15 - 2 * (ee303/ee15)) * ee10/xi + ee189)/ee26 +
   ee187 * ee13 * ee101 * ee10 - ((ee18 + 2 * (ee232 * ee15) -
   2 * (ee77 * ee18 * ee61/ee33)) * ee6 + xi * ((ee217 * ee7 +
   ee438 - (ee35 * (4 * ee275 - 8 * (R_pow(ee4, ee5 + ee358) * ee61/ee33)) +
   4 * (ee335 * ee77)) * ee6/ee33) * ee6 + ee296) * ee11)/ee33) * ee6)/ee1 -
   ee267)/ee75 - (ee317 * ee173 +
   ee150 * ((2 * ee360 - ee434)/ee1 + 2 * ee131 - 2 * ee130))/ee129);
out(j, 6) = ((((ee431 * ee30 - (ee402/ee19 + 2 * (ee337/ee19)) * ee27 * ee25/xi)/ee26 +
   (2 * (R_pow(ee4, 1 - ee98) * ee30 * ee54 * ee25) -
   ee176/ee72) * ee54 - (ee373 + ee422)/ee72) * ee14 +
   (ee431 * ee18 - (ee176/ee19 + 2 * (ee171/ee19)) * ee14 * ee25/xi)/ee15 +
   (ee18 * (y * (ee404 - (4 * ee97 - 4 * ee19)/ee29)/ee1 -
   (2 * (ee134 * ee54) + 2 * (ee10 * ee165/xi) -
   ((4 * ee41 - 12 * ee40)/xi + y * (ee121 + ee64)/ee38) * ee26)/xi) -
   (ee373 + 2 * ee422) * ee14)/ee72 - ((ee417 + (2 * ee185 -
   (2 * (ee274 * ee50/ee33) + 2 * (ee203/ee116)) * ee25) * ee50) * ee61/ee33 +
   2 * (ee229 * ee186 * ee54)))/xi -
   (((((ee185 - (ee205 + 2 * (ee374 * ee6/ee33))) * ee50 + 2 * ee417) * ee6 +
   ee35 * (y * (ee403 - (2 * ee22 - (2 * ee24 +
   ee210))/ee29)/ee1 - (2 * (ee50 * ee67) + ee10 * ee149/xi -
   ee353 * ee15)/xi))/xi - (ee35 * (2 * ee407 - 8 * (ee276 * ee180/ee33)) +
   2 * (ee374/xi)) * ee61 * ee50/ee33) * ee6/ee33 +
   ((2 * (ee411 * ee18) + 2 * (ee229 * ee26/xi))/ee302 - 8 * (R_pow(ee4, ee442 -
   ee310) * ee18 * ee209)) * ee54))/ee75 +
   (ee266 + 2 * ee316 + 2 * (R_pow(ee166, 2) * ee75/ee129) - ee271) * ee166/ee129;
out(j, 7) = ((((ee344 * ee108 - (ee397/ee19 + ee105/ee90) * ee69/xi)/ee26 +
   (2 * (ee249 * ee10 * ee54) - ee272/ee72) * ee54 -
   (ee190 * ee54 * ee25 + ee104 * ee165/xi)/ee72)/xi - ((ee273 * ee149 +
   ((ee259 - ee231) * ee50 * ee25 - 2 * (ee197 * ee77 * ee6/ee33)) * ee6)/xi -
   (ee35 * (4 * ee15 - 8 * ee223) +
   2 * (ee273 * ee15)) * ee61 * ee180/ee33)/ee33) * ee6 +
   ((ee340 * ee108 - (ee396/ee19 + ee55/ee90) * ee69/xi)/ee15 -
   (ee220 * ee15 + (ee259 - (ee419 + 2 * (R_pow(ee4, ee60 - ee8) * ee18 * ee6/ee33)) * ee6) * ee50 * ee25) * ee61/ee33)/xi)/ee75 +
   (ee166 * (2 * ee318 - ee433) - ee316 * ee150)/ee129;
out(j, 8) = ((((((ee432 * (ee259 - ee142) - (ee419 + ee18/(ee17 * ee19))) * ee6 +
   ee142 + ee259 - (ee142 + ee259))/ee17 +
   (ee187/ee19 - ee259) * ee10/xi) * ee6 + ee142)/ee15 - (ee248 +
   2 * (ee257 * ee15) - 2 * (R_pow(ee4, ee355 - ee5) * ee18 * ee61/ee33)) * ee61/ee33) * ee25/xi +
   ((((ee395/ee19 - 2 * (ee303/ee19)) * ee10/xi +
   ee179) * ee25/ee26 + ee240 * ee54/ee72)/xi -
   ((ee193 + ee357 + 3 - (ee35 * (4 * ee77 - 8 * (R_pow(ee4, ee358) * ee61/ee33)) +
   4 * (ee273 * ee77)) * ee6/ee33) * ee6 +
   ee35) * ee50/ee33) * ee6)/ee75 - (ee317 * ee166 +
   ee150 * (2 * ee279 + 2 * ee289 - ee433))/ee129;
out(j, 9) =  - (((((ee217 * ee408 + (ee408 - ee258) * ee80 -
   (ee35 * (4 * ee198 - 8 * (R_pow(ee4, ee11 + ee358) * ee61/ee33)) +
   4 * ((ee243 + ee114) * ee77)) * ee6/ee33) * ee6 + ee114)/ee33 +
   ((ee342 - ee169) * ee6 * ee10/xi - ee18) * ee10/ee76) * ee6 -
   ((((ee441 - ((ee351/ee17 + ee357) * ee6 + ee20 +
   ee162)) * ee10/xi + (ee441 - ((2 * (ee424 * ee17) + 2 * ((1 +
   2 * ee6) * ee20) - 8 * ee305)/ee17 + 2 * ee424 + 3 * ((ee71 +
   ee439) * ee6)))/ee17) * ee6 + ee20 - ee20)/ee17 + 1 -
   (((ee351 - ee350) * ee10/xi - ((ee357 + ee439) * ee6 - ee351)/ee17) * ee6 -
   ee20) * ee10/xi)/ee15)/ee75 - (ee267 + 2 * ee317 -
   (ee268 + 2 * (R_pow(ee150, 2) * ee75/ee129))) * ee150/ee129);

}

return out;

}
