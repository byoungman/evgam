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
arma::vec lkappa = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();

if (dcate == 1) {
  lpsivec = lpsivec.elem(dupid);
  xivec = xivec.elem(dupid);
  lkappavec = lkappavec.elem(dupid);
}

double y, lpsi, xi, lkappa;
double ee1, ee2;
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
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 9);

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    xivec = xivec.elem(dupid);
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
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 9);

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    xivec = xivec.elem(dupid);
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
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];

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
