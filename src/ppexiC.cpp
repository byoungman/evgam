// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// added with evgam_0.1.2 (05/04/2020)

// [[Rcpp::export]]
double ppexi1d0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, arma::vec yvec, arma::vec wvec)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
arma::vec thetavec = X4 * Rcpp::as<arma::vec>(pars[3]);
int nobs = yvec.size();

double y, mu, lpsi, xi, theta;
double ee1;
double nllh=0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
theta = 1 / (1 + exp(-thetavec[j]));

ee1 = xi * (y - mu) / exp(lpsi);

// if (ee1 <= -1.0) {
//     nllh = 1e20;
//     break;
// } else {
// 
// nllh += wvec[j] * R_pow(1.0 + ee1, -1/xi);
//     
// }

if (ee1 > -1.0) {
  nllh += wvec[j] * theta * R_pow(1.0 + ee1, -1/xi);
}

}

return(nllh);

}

// [[Rcpp::export]]
arma::mat ppexi1d12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, arma::vec yvec, arma::vec wvec)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
arma::vec thetavec = X4 * Rcpp::as<arma::vec>(pars[3]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 14);
out.zeros();

double y, w, mu, lpsi, xi, theta;
double ee1, ee2, ee4, ee5, ee6, ee8, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
double ee20, ee23, ee25, ee27, ee29;
double ee30, ee33, ee36, ee37;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
theta = thetavec[j];
w = wvec[j];

if (xi * (y - mu) / exp(lpsi) > -1.0) {
    
ee1 = exp(lpsi);
ee2 = y - mu;
ee4 = xi * ee2/ee1;
ee5 = 1/xi;
ee6 = 1 + ee4;
ee8 = exp(-theta);
ee9 = 1 + ee5;
ee10 = 1 + ee8;
ee11 = R_pow(ee6, ee9);
ee12 = log1p(ee4);
ee13 = R_pow(ee6, ee5 + 2);
ee14 = R_pow(ee6, ee5);
ee15 = R_pow(ee10, 2);
ee16 = ee13 * ee1;
ee17 = ee10 * ee1;
ee18 = w * ee8;
ee20 = ee9 * ee2/ee16;
ee23 = 1/ee11;
ee25 = ee12/(xi * ee14) - ee2/(ee11 * ee1);
ee27 = ee10 * ee11 * ee1;
ee29 = ee15 * ee11 * ee1;
ee30 = ee15 * ee14;
ee33 = w * (ee12/(R_pow(xi, 2) * ee11) - ee20);
ee36 = xi * ee9 * ee2/ee16;
ee37 = xi * ee10;

out(j, 0) = w/ee27;
out(j, 1) = w * ee2/ee27;
out(j, 2) = w * ee25/ee37;
out(j, 3) = ee18/ee30;
out(j, 4) = w * xi * ee9/(ee10 * ee13 * R_pow(ee1, 2));
out(j, 5) = w * (ee36 - ee23)/ee17;
out(j, 6) = ee33/ee17;
out(j, 7) = ee18/ee29;
out(j, 8) = -(w * (ee23 - ee36) * ee2/ee17);
out(j, 9) = ee33 * ee2/ee17;
out(j, 10) = ee18 * ee2/ee29;
out(j, 11) = w * (((ee2/(ee6 * ee1) - 2 * (ee12/xi))/ee14 +
   ee12 * ee25/xi)/xi + (ee20 + (ee23 - ee12/(xi * ee11))/xi) 
  * ee2/ee1)/ee37;
out(j, 12) = ee18 * ee25/(xi * ee15);
out(j, 13) = -(w * (1 - 2 * (ee8/ee10)) * ee8/ee30);

    
}

}

return out;

}

// [[Rcpp::export]]
arma::mat ppexi1d34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, arma::vec yvec, arma::vec wvec)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
arma::vec thetavec = X4 * Rcpp::as<arma::vec>(pars[3]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 20);
out.zeros();

double y, w, mu, lpsi, xi, theta;
double ee1, ee2, ee4, ee5, ee6, ee7, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee20, ee21, ee23, ee24, ee25, ee26, ee27, ee28;
double ee31, ee32, ee33, ee35, ee36, ee37;
double ee41, ee42, ee43, ee46, ee48, ee49;
double ee50, ee52, ee53, ee56;
double ee60, ee62, ee63, ee65, ee68, ee69;
double ee70, ee72, ee75, ee77, ee78;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
theta = thetavec[j];
w = wvec[j];

if (xi * (y - mu) / exp(lpsi) > -1.0) {
    
ee1 = exp(lpsi);
ee2 = y - mu;
ee4 = xi * ee2/ee1;
ee5 = 1/xi;
ee6 = 1 + ee4;
ee7 = 1 + ee5;
ee9 = exp(-theta);
ee10 = ee5 + 2;
ee11 = log1p(ee4);
ee12 = R_pow(ee6, ee10);
ee13 = R_pow(ee6, ee7);
ee14 = 1 + ee9;
ee15 = R_pow(xi, 2);
ee16 = R_pow(ee14, 2);
ee17 = ee12 * ee1;
ee18 = R_pow(ee6, ee5 + 3);
ee19 = ee6 * ee1;
ee20 = ee18 * ee1;
ee21 = ee15 * ee12;
ee23 = ee7 * ee2/ee17;
ee24 = R_pow(ee6, ee5);
ee25 = ee2/ee19;
ee26 = 1/ee12;
ee27 = ee11/xi;
ee28 = ee14 * ee1;
ee31 = ee25 - 2 * ee27;
ee32 = 1/ee13;
ee33 = ee11/(ee15 * ee13);
ee35 = ee11/ee21 - ee10 * ee2/ee20;
ee36 = xi * ee7;
ee37 = ee16 * ee1;
ee41 = ee11/(xi * ee24) - ee2/(ee13 * ee1);
ee42 = ee33 - ee23;
ee43 = ee31/ee13;
ee46 = (1/ee21 - ee7 * ee35) * ee2/ee1;
ee48 = 2/ee12;
ee49 = R_pow(ee1, 2);
ee50 = ee11/(xi * ee13);
ee52 = w * (1 - 2 * (ee9/ee14)) * ee9;
ee53 = w * ee9;
ee56 = xi * ee10 * ee2/ee20;
ee60 = (ee31/ee24 + ee11 * ee41/xi)/xi + (ee23 + (ee32 -  ee50)/xi) * ee2/ee1;
ee62 = (ee7 * (ee26 + xi * ee35) + ee26) * ee2/ee1;
ee63 = ee14 * ee49;
ee65 = ee16 * ee13 * ee1;
ee68 = 2/ee13;
ee69 = w * ((ee43 + ee11 * ee42)/ee15 + ee46);
ee70 = ee53 * ee42;
ee72 = w * xi * ee7;
ee75 = ee36 * (ee26 + ee48 - ee56) * ee2/ee1;
ee77 = ee36 * ee2/ee17;
ee78 = xi * ee16;

out(j, 0) = w * ee15 * (1 + ee7) * ee7/(ee14 * ee18 * R_pow(ee1, 3));
out(j, 1) = -(ee72 * (ee48 - ee56)/ee63);
out(j, 2) = w * (ee26 + ee36 * ee35)/ee63;
out(j, 3) = ee72 * ee9/(ee16 * ee12 * ee49);
out(j, 4) = w * (ee32 - ee75)/ee28;
out(j, 5) = w * (ee62 - ee33)/ee28;
out(j, 6) = ee53 * (ee77 - ee32)/ee37;
out(j, 7) = ee69/ee28;
out(j, 8) = ee70/ee37;
out(j, 9) = -(ee52/ee65);
out(j, 10) = -(w * (ee75 - ee32) * ee2/ee28);
out(j, 11) = -(w * (ee33 - ee62) * ee2/ee28);
out(j, 12) = -(w * (ee32 - ee77) * ee9 * ee2/ee37);
out(j, 13) = ee69 * ee2/ee28;
out(j, 14) = ee70 * ee2/ee37;
out(j, 15) = -(ee52 * ee2/ee65);
out(j, 16) = w * (((ee60 * ee11 + 2 * (ee31 * ee41))/xi - ((ee25 +
   2/xi) * ee2/ee19 + (2 * ee25 - 6 * ee27)/xi)/ee24)/xi +
   ((((ee23 + (ee68 - ee50)/xi) * ee11 - (ee43 + ee68))/xi -
   2 * ee23)/xi - ee46) * ee2/ee1)/(xi * ee14);
out(j, 17) = w * ee60 * ee9/ee78;
out(j, 18) = -(ee52 * ee41/ee78);
out(j, 19) = w * (1 - ((2 * (1 + 2 * ee9) + 2 * ee14 - 8 * ee9)/ee14 +
   2) * ee9/ee14) * ee9/(ee16 * ee24);

}

}

return out;

}

// [[Rcpp::export]]
double ppexi2d0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, arma::vec yvec, arma::vec wvec)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
arma::vec thetavec = X4 * Rcpp::as<arma::vec>(pars[3]);
double y, mu, lpsi, xi, theta;

int nobs = yvec.size();

double ee2, ee3, ee6;
double nllh=0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
theta = thetavec[j];

ee2 = 1/xi;
ee3 = exp(-theta);
ee6 = xi * (y - mu)/exp(lpsi);

if (ee6 <= -1.0) {
    nllh = 1e20;
    break;
} else {

nllh += wvec[j] * ((1 + ee2) * log1p(ee6) + 1/((1 + ee3) * R_pow(1 + ee6, ee2)) + log1p(ee3) + lpsi);
    
}
 
}

return(nllh);

}

// [[Rcpp::export]]
arma::mat ppexi2d12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, arma::vec yvec, arma::vec wvec)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
arma::vec thetavec = X4 * Rcpp::as<arma::vec>(pars[3]);
double y, mu, lpsi, xi, theta;

int nobs = yvec.size();

arma::mat out = arma::mat(nobs, 14);

double ee1, ee2, ee4, ee5, ee6, ee7, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee17, ee18, ee19;
double ee21, ee22, ee23, ee24, ee26;
double ee31, ee32, ee34, ee35, ee37, ee38, ee39;
double ee41, ee42, ee43, ee44, ee46, ee49;
double ee51, ee52, ee53, ee55, ee56, ee59;
double ee60, ee61, ee65, ee66, ee67, ee69;
double ee70, ee72, ee75, ee77, ee78, ee79;
double ee80, ee82, ee83, ee84, ee85;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
theta = thetavec[j];

ee1 = exp(lpsi);
ee2 = y - mu;
ee4 = xi * ee2/ee1;
ee5 = 1/xi;
ee6 = 1 + ee4;
ee7 = 1 + ee5;
ee9 = exp(-theta);
ee10 = R_pow(ee6, ee7);
ee11 = 1 + ee9;
ee12 = ee5 + 2;
ee13 = R_pow(ee6, ee12);
ee14 = R_pow(ee6, ee5);
ee15 = log1p(ee4);
ee17 = 1/(ee11 * ee14);
ee18 = 1/ee10;
ee19 = xi * ee7;
ee21 = exp(-ee17);
ee22 = ee11 * ee10;
ee23 = ee13 * ee1;
ee24 = R_pow(ee6, 2 * ee7);
ee26 = ee7 * ee2/ee23;
ee31 = 1/(ee11 * ee24) - ee19/ee13;
ee32 = 2/xi;
ee34 = ee15/(xi * ee14) - ee2/(ee10 * ee1);
ee35 = R_pow(ee1, 2);
ee37 = R_pow( - (ee21/(ee22 * ee1)), 2);
ee38 = R_pow(ee11, 2);
ee39 = R_pow(ee21, 2);
ee41 = ee37 * ee38 * ee35;
ee42 = ee15/(xi * ee10);
ee43 = R_pow(xi, 2);
ee44 = ee11 * R_pow(ee6, 1 + ee32);
ee46 = (ee34/ee22 - ee42)/xi + ee26;
ee49 = ee31 * ee2/ee1 + ee18;
ee51 = 1/ee44 - ee18;
ee52 = R_pow(ee6, ee5 + 3);
ee53 = 1/ee13;
ee55 = ee15/(ee43 * ee10) - ee26;
ee56 = ee52 * ee1;
ee59 = ee19 * ee2/ee23 - ee18;
ee60 = ee43 * ee13;
ee61 = R_pow(ee6, ee32 + 3);
ee65 = ee15/ee60 - ee12 * ee2/ee56;
ee66 = ee15/xi;
ee67 = xi * ee12;
ee69 = ee49/ee10;
ee70 = (ee55/ee10 - ee46/ee10)/ee11;
ee72 = ee11 * ee1;
ee75 = ee59/ee10;
ee77 = ee2/(ee6 * ee1) - 2 * ee66;
ee78 = 1 - ee17;
ee79 = 1 - 2 * (ee9/ee11);
ee80 = 1/(ee11 * ee61);
ee82 = 1/ee24 - ee51/ee10;
ee83 = 2/ee13;
ee84 = ee67/ee52;
ee85 = xi * ee65;

out(j, 0) = ee10 * ee31/ee1;
out(j, 1) = ee49 * ee10;
out(j, 2) = ee46 * ee10;
out(j, 3) = ee10 * ee51 * ee9/ee11;
out(j, 4) = (((ee19/ee61 - ee31/ee10)/ee11 - ee19 * (ee84 -
   ee80)) * ee10 + R_pow(ee31, 2) * ee39/ee41)/ee35;
out(j, 5) = (((ee75 - ee69)/ee11 + xi * ((ee80 - ee84) * ee2/ee1 +
   ee83) * ee7) * ee10 + ee49 * ee31 * ee39/ee41)/ee1;
out(j, 6) = (ee46 * ee31 * ee39/ee41 + (ee70 - (ee7 * (ee85 -
   ee34/(ee11 * ee13)) + ee53)) * ee10)/ee1;
out(j, 7) = ((ee82/ee11 - ee19 * (ee53 - 1/(ee11 * R_pow(ee6, 2 +
   ee32)))) * ee10 + ee51 * ee31 * ee39/ee41) * ee9/ee72;
out(j, 8) = (((2 * ee75 - ee69)/ee11 + ee19 * (ee53 + ee83 -
   ee67 * ee2/ee56)) * ee2/ee1 - ee18) * ee10 + R_pow(ee49, 2) * ee39/ee41;
out(j, 9) = ((ee70 - (ee7 * (ee53 + ee85) + ee53)) * ee2/ee1 +
   (ee34 * ee59/ee11 + ee42)/xi) * ee10 + ee46 * ee49 * ee39/ee41;
out(j, 10) = (ee49 * ee51 * ee39/ee41 + (ee82 * ee2/ee72 - ee78 * ee59) * ee10) * ee9/ee11;
out(j, 11) = (((((ee77/ee14 + ee34 * (ee66 - ee34/ee11))/xi +
   (ee26 + (ee18 - ee42)/xi) * ee2/ee1)/ee10 + 2 * (ee34 * ee55))/ee11 -
   (ee77/ee10 + ee15 * ee55)/xi)/xi - (1/ee60 - ee7 * ee65) * ee2/ee1) * ee10 +
   R_pow(ee46, 2) * ee39/ee41;
out(j, 12) = (((ee78/ee10 + ee18) * ee34/(xi * ee11) + (ee17 -
   1) * ee55) * ee10 + ee46 * ee51 * ee39/ee41) * ee9/ee11;
out(j, 13) = (R_pow(ee51, 2) * ee39 * ee9/(ee37 * R_pow(ee11, 3) * ee35) - (((ee79/ee14 + ee9/(ee38 * R_pow(ee6, ee32)))/ee10 -
   2 * (ee9/ee44))/ee11 - ee79/ee10) * ee10) * ee9/ee11;

}

return out;

}

// [[Rcpp::export]]
arma::mat ppexi2d34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, arma::vec yvec, arma::vec wvec)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
arma::vec thetavec = X4 * Rcpp::as<arma::vec>(pars[3]);
double y, mu, lpsi, xi, theta;

int nobs = yvec.size();

arma::mat out = arma::mat(nobs, 20);

double ee1, ee2, ee4, ee5, ee6, ee7, ee8;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee21, ee22, ee24, ee28, ee29;
double ee30, ee32, ee33, ee34, ee37, ee38, ee39;
double ee40, ee41, ee42, ee43, ee44, ee45, ee47, ee49;
double ee50, ee51, ee52, ee53, ee54, ee55, ee56, ee58;
double ee61, ee62, ee63, ee65, ee66, ee67, ee68, ee69;
double ee70, ee71, ee72, ee73, ee74, ee75, ee76, ee78, ee79;
double ee81, ee82, ee83, ee84, ee85, ee86, ee87, ee89;
double ee90, ee92, ee94, ee95, ee96, ee97, ee98, ee99;
double ee101, ee103, ee104, ee106, ee107, ee108, ee109;
double ee110, ee111, ee112, ee113, ee114, ee115, ee116, ee117, ee118, ee119;
double ee121, ee122, ee124, ee125, ee126, ee129;
double ee130, ee131, ee132, ee135, ee136, ee137, ee138, ee139;
double ee140, ee144, ee145, ee146, ee147, ee149;
double ee151, ee152, ee153, ee157, ee159;
double ee162, ee163, ee164, ee166, ee167, ee169;
double ee172, ee174, ee175, ee176, ee177, ee178, ee179;
double ee180, ee181, ee182, ee184, ee185, ee186, ee189;
double ee193, ee195, ee196, ee197, ee198, ee199;
double ee200, ee201, ee202, ee203, ee205, ee207, ee208;
double ee210, ee211, ee213, ee215, ee216, ee218;
double ee222, ee224, ee225, ee226, ee227, ee229;
double ee234, ee236, ee237, ee238, ee239;
double ee241, ee242, ee243, ee244, ee246, ee247;
double ee250, ee251, ee252, ee253, ee254, ee257, ee258;
double ee260, ee261;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
theta = thetavec[j];

ee1 = exp(lpsi);
ee2 = y - mu;
ee4 = xi * ee2/ee1;
ee5 = 1/xi;
ee6 = 1 + ee4;
ee7 = 1 + ee5;
ee8 = R_pow(ee6, ee7);
ee10 = exp(-theta);
ee11 = 1 + ee10;
ee12 = ee5 + 2;
ee13 = R_pow(ee6, ee12);
ee14 = log1p(ee4);
ee15 = R_pow(ee6, ee5);
ee16 = 1/ee8;
ee17 = xi * ee7;
ee18 = ee13 * ee1;
ee19 = R_pow(xi, 2);
ee21 = 1/(ee11 * ee15);
ee22 = ee5 + 3;
ee24 = ee7 * ee2/ee18;
ee28 = ee14/(xi * ee15) - ee2/(ee8 * ee1);
ee29 = R_pow(ee6, ee22);
ee30 = 2/xi;
ee32 = exp(-ee21);
ee33 = ee11 * ee8;
ee34 = R_pow(ee6, 2 * ee7);
ee37 = 1/(ee11 * ee34) - ee17/ee13;
ee38 = ee29 * ee1;
ee39 = ee14/(xi * ee8);
ee40 = ee19 * ee13;
ee41 = ee14/(ee19 * ee8);
ee42 = 1/ee13;
ee43 = ee41 - ee24;
ee44 = ee14/xi;
ee45 = R_pow(ee6, 1 + ee30);
ee47 = ee12 * ee2/ee38;
ee49 = R_pow( - (ee32/(ee33 * ee1)), 2);
ee50 = R_pow(ee32, 2);
ee51 = R_pow(ee1, 2);
ee52 = ee11 * ee45;
ee53 = R_pow(ee11, 2);
ee54 = ee14/ee40;
ee55 = ee54 - ee47;
ee56 = ee6 * ee1;
ee58 = (ee28/ee33 - ee39)/xi + ee24;
ee61 = ee37 * ee2/ee1 + ee16;
ee62 = ee2/ee56;
ee63 = 1/ee52;
ee65 = ee17 * ee2/ee18;
ee66 = ee62 - 2 * ee44;
ee67 = ee63 - ee16;
ee68 = ee65 - ee16;
ee69 = xi * ee12;
ee70 = ee49 * ee53;
ee71 = 1 - 2 * (ee10/ee11);
ee72 = 2/ee13;
ee73 = ee70 * ee51;
ee74 = R_pow(ee6, ee30 + 3);
ee75 = 1 - ee21;
ee76 = xi * ee55;
ee78 = ee69 * ee2/ee38;
ee79 = ee7 * ee55;
ee81 = (ee24 + (ee16 - ee39)/xi) * ee2/ee1;
ee82 = ee66/ee8;
ee83 = R_pow(ee6, ee30);
ee84 = ee61/ee8;
ee85 = ee82 + ee14 * ee43;
ee86 = ee66/ee15;
ee87 = ee71/ee15;
ee89 = (1/ee40 - ee79) * ee2/ee1;
ee90 = ee58/ee8;
ee92 = ee7 * (ee42 + ee76) + ee42;
ee94 = ee11 * ee1;
ee95 = ee28/ee11;
ee96 = ee43/ee8;
ee97 = 1/(ee11 * ee74);
ee98 = ee42 + ee72;
ee99 = 1/ee34;
ee101 = ee49 * ee11 * ee8;
ee103 = (ee86 + ee28 * (ee44 - ee95))/xi + ee81;
ee104 = ee68/ee8;
ee106 = 2/ee8;
ee107 = ee17 * (ee98 - ee78);
ee108 = ee69/ee29;
ee109 = ee101 * ee51;
ee110 = (ee96 - ee90)/ee11;
ee111 = ee87 + ee10/(ee53 * ee83);
ee112 = R_pow(ee6, ee5 + 4);
ee113 = ee70 * ee8;
ee114 = xi * ee11;
ee115 = ee19 * ee29;
ee116 = ee103/ee8;
ee117 = ee85/xi;
ee118 = ee71/ee8;
ee119 = ee112 * ee1;
ee121 = ee99 - ee67/ee8;
ee122 = ee113 * ee51;
ee124 = ((ee116 + 2 * (ee28 * ee43))/ee11 - ee117)/xi -  ee89;
ee125 = (ee110 - ee92) * ee2;
ee126 = ee111/ee8;
ee129 = ee21 - 1;
ee130 = 2 * (ee37 * ee50/ee109);
ee131 = 2 * ee104;
ee132 = ee125/ee1;
ee135 = (ee75/ee8 + ee16) * ee28/ee114 + ee129 * ee43;
ee136 = ((ee131 - ee84)/ee11 + ee107) * ee2;
ee137 = (ee28 * ee68/ee11 + ee39)/xi;
ee138 = ee11 * R_pow(ee6, 2 + ee30);
ee139 = ee37/ee8;
ee140 = ee121 * ee2;
ee144 = ee14/ee115 - ee22 * ee2/ee119;
ee145 = ee17 * (ee108 - ee97);
ee146 = ee17/ee74;
ee147 = ee132 + ee137;
ee149 = ee85/ee19 + ee89;
ee151 = ee136/ee1 - ee16;
ee152 = ee75 * ee68;
ee153 = ee7 * (ee76 - ee28/(ee11 * ee13));
ee157 = ee140/ee94;
ee159 = (ee146 - ee139)/ee11 - ee145;
ee162 = ee42 + ee17 * ee55;
ee163 = ee72 - ee78;
ee164 = 2/ee29;
ee166 = xi * ((ee97 - ee108) * ee2/ee1 + ee72) * ee7;
ee167 = ee17 * (ee42 - 1/ee138);
ee169 = (ee126 - 2 * (ee10/ee52))/ee11 - ee118;
ee172 = ee92 * ee2/ee1 - ee41;
ee174 = (ee104 - ee84)/ee11 + ee166;
ee175 = ee153 + ee42;
ee176 = R_pow(ee6, ee30 + 4);
ee177 = ee157 - ee152;
ee178 = ee16 - ee107 * ee2/ee1;
ee179 = ee130 - ee16;
ee180 = 2/ee45;
ee181 = 2/ee15;
ee182 = xi * ee22;
ee184 = (ee66/ee13 + ee14 * ee55)/ee19 + (1/ee115 - ee12 *  ee144) * ee2/ee1;
ee185 = ee61/ee13;
ee186 = ee110 - ee175;
ee189 = ee67/ee13 - 1/ee74;
ee193 = ee12/ee29;
ee195 = (2 * (R_pow(ee67, 2) * ee50/ee109) + ee180) * ee10/ee11;
ee196 = ee163/ee8;
ee197 = 1/(ee11 * ee83);
ee198 = ee16 - ee65;
ee199 = 1/ee29;
ee200 = 2 - ee21;
ee201 = ee130 - ee106;
ee202 = 2 * ee47;
ee203 = ee113 * ee1;
ee205 = ee124/ee8 + 2 * (ee58 * ee43);
ee207 = ee169/ee8 - (ee71/ee34 + 2 * (ee67 * ee10/ee33));
ee208 = ee149/ee8;
ee210 = ee135/ee8 + ee90;
ee211 = ee172/ee8;
ee213 = (((ee181 - ee197)/ee11 - 2) * ee10 - ee87)/ee11 +  1;
ee215 = ee58 * ee61 * ee50;
ee216 = ee58 * ee68;
ee218 = (ee86 + ee14 * ee28/xi)/xi + ee81;
ee222 = (ee67 * ee179 + ee99)/ee11 - ee167;
ee224 = ee61 * ee67 * ee50;
ee225 = ee61 * ee43;
ee226 = ee61 * ee68;
ee227 = ee185 + ee196;
ee229 = (ee195 - ee126)/ee11 + ee118;
ee234 = (ee62 + ee30) * ee2/ee56 + (2 * ee62 - 6 * ee44)/xi;
ee236 = (ee63 - ee106) * ee68;
ee237 = ee37/ee13;
ee238 = ee162 * ee28;
ee239 = ee162/ee8;
ee241 = ee121/ee11 - ee167;
ee242 = ee200 * ee28;
ee243 = (ee106 - ee63) * ee43;
ee244 = 1 - ((2 * (1 + 2 * ee10) + 2 * ee11 - 8 * ee10)/ee11 +  2) * ee10/ee11;
ee246 = 1/(ee11 * ee176) - ee182/ee112;
ee247 = ee16 + ee106;
ee250 = 2 * ee177;
ee251 = 2 * (ee14/(xi * ee13));
ee252 = 4/ee13;
ee253 = 4/ee29;
ee254 = xi * ee184;
ee257 = xi * ee58 * ee7/ee13;
ee258 = ee69/ee176;
ee260 = ee182 * ee2/ee119;
ee261 = xi * (2 * ee79 - (ee7 * (ee164 + ee69 * ee144) +  ee193) * ee2/ee1);

out(j, 0) =  - ((((ee159/ee8 + ee17 * (2 * ee237 - ee258))/ee11 +
   xi * ((ee237 - ee258)/ee11 - xi * ee246 * ee12) * ee7) * ee8 -
   ((ee37 * ee179 + ee146)/ee11 + 2 * ee159 - ee145) * ee37 * ee50/ee73)/R_pow(ee1, 3));
out(j, 1) = -((((ee174/ee8 + ee37 * ee68 + xi * ee227 * ee7)/ee11 +
   xi * (ee227/ee11 - xi * (ee246 * ee2/ee1 + ee199 + ee164) * ee12) * ee7) * ee8 -
   (ee61 * ee159 + ee37 * (2 * (ee61 * ee37 * ee50/ee122) +
   2 * ee174)) * ee50/ee73)/ee51);
out(j, 2) = -((((ee186/ee8 + ee37 * ee43 + 2 * ee257 - 2 * ee239)/ee11 -
   xi * (((ee28/(ee11 * ee29) - xi * ee144) * ee12 -
   ee164) * ee7 - ee193)) * ee8 - (ee58 * ee159 + ee37 * (2 * (ee58 * ee37 * ee50/ee122) +
   2 * ee186)) * ee50/ee73)/ee51);
out(j, 3) = -((((ee241/ee8 + ee139 + xi * ee189 * ee7)/ee11 +
   xi * (ee189/ee11 - xi * (ee97 - ee199) * ee12) * ee7) * ee8 -
   (ee159 * ee67 + ee37 * (2 * (ee67 * ee37 * ee50/ee122) +
   2 * ee241)) * ee50/ee73) * ee10/(ee11 * ee51));
out(j, 4) = -((((ee151/ee8 + 2 * ee226 - ee178/ee8)/ee11 + xi * (((ee185 +
   2 * ee196)/ee11 - xi * (ee199 + ee253 - ee260) * ee12) * ee2/ee1 +
   ee252) * ee7) * ee8 - (((ee61 * ee201 +
   ee131)/ee11 + 2 * ee166) * ee61 + ee151 * ee37) * ee50/ee73)/ee1);
out(j, 5) = -((((ee147/ee8 + ee216 + ee225 + ee7 * ee163 * ee28 +
   (ee257 - ee239) * ee2/ee1 - ee211)/ee11 - (ee72 + ee261)) * ee8 -
   (((ee58 * ee179 + ee96)/ee11 - ee175) * ee61 + ee147 * ee37 +
   ee58 * ee174) * ee50/ee73)/ee1);
out(j, 6) = -((((ee84 + ee177/ee8 + ee236)/ee11 + xi * (ee189 * ee2/ee94 -
   ee75 * ee163) * ee7) * ee8 - (ee222 * ee61 +
   ee174 * ee67 + ee177 * ee37) * ee50/ee73) * ee10/ee94);
out(j, 7) = ((ee124 * ee37 + ((ee58 * ee201 + 2 * ee96)/ee11 -
   (2 * ee153 + ee72)) * ee58) * ee50/ee73 + ((ee208 + ee103 * ee7/ee13 +
   ee238/xi - ee205)/ee11 + ee202 - ((ee251 - ee238/ee11)/xi +
   ee254 * ee7)) * ee8)/ee1;
out(j, 8) = ((ee135 * ee37 + ee222 * ee58 + ee186 * ee67) * ee50/ee73 +
   (((ee75/ee13 + ee42) * ee7 * ee28 + ee243 - ee210)/ee11 -
   ee75 * ee162) * ee8) * ee10/ee94;
out(j, 9) = ((ee207/ee11 + xi * (ee71/ee13 + (2 * (ee10/ee138) -
   ee111/ee13)/ee11) * ee7) * ee8 + (((ee67 * ee201 + ee99 +
   ee99)/ee11 - 2 * ee167) * ee67 * ee10/ee11 - ee169 * ee37) * ee50/ee73) * ee10/ee94;
out(j, 10) =  - ((((((ee107 - (ee84 + 2 * (ee198/ee8))/ee11) * ee2/ee1 -
   ee16)/ee8 + 3 * ee226 - ee178 * ee247)/ee11 + ee17 * (ee98 +
   ee252 - ee69 * (ee164 + ee253 - ee260) * ee2/ee1)) * ee2/ee1 -
   ee16) * ee8 - ((ee136 + 2 * (R_pow(ee61, 2) * ee50/ee203))/ee1 +
   2 * ee151 - ee16) * ee61 * ee50/ee73);
out(j, 11) = -(((((ee132 + (ee39 - ee198 * ee28/ee11)/xi)/ee8 +
   ee225 + 2 * ee216 - 2 * ee211)/ee11 - (ee92 + ee72 + ee261)) * ee2/ee1 +
   (ee39 - ee178 * ee28/ee11)/xi) * ee8 - (ee151 * ee58 +
   ee61 * (2 * (ee215/ee122) + 2 * ee147)) * ee50/ee73);
out(j, 12) = -(((((ee75 * ee198 + ee157)/ee8 + ee84 + 2 * ee236) * ee2/ee94 +
   ee75 * ee178) * ee8 - (ee151 * ee67 + ee61 * (2 * (ee224/ee122) +
   ee250)) * ee50/ee73) * ee10/ee11);
out(j, 13) = (ee124 * ee61 + ((ee125 + 2 * (ee215/ee203))/ee1 +
   ee137) * ee58) * ee50/ee73 + ((((ee208 - ee205)/ee11 + ee202 -
   (ee7 * (ee54 + ee254 - ee47) + (2 * (ee14/ee13) - ee42)/ee19)) * ee2/ee1 +
   ((ee103 * ee68 + 2 * (ee172 * ee28))/ee11 +
   ee117)/xi) * ee8 + ee147 * ee58 * ee50/ee73);
out(j, 14) = ((ee147 * ee67 + ee135 * ee61 + ((ee140 + 2 * (ee224/(ee101 * ee1)))/ee94 -
   ee152) * ee58) * ee50/ee73 + (((ee243 -
   ee210) * ee2/ee1 + ee242 * ee68/xi)/ee11 - ee172 * ee75) * ee8) * ee10/ee11;
out(j, 15) = ((ee207 * ee2/ee94 + ee213 * ee68) * ee8 + (ee229 * ee61 +
   ee250 * ee67 * ee10/ee11) * ee50/ee73) * ee10/ee11;
out(j, 16) = ((((((ee218 * ee14 + 2 * (ee66 * ee28))/xi - (((ee66 * (1/ee15 +
   ee181) + (3 * ee44 - ee95) * ee28)/xi + 3 * ee81) * ee28/ee11 +
   ee234/ee15))/xi + ((((ee24 + (ee106 -
   ee39)/xi) * ee14 - (ee82 + ee106))/xi - 2 * ee24)/xi - ee89) * ee2/ee1)/ee8 +
   3 * (ee149 * ee28) + 3 * (ee103 * ee43))/ee11 -
   (ee149 * ee14 + 2 * (ee66 * ee43) - ee234/ee8)/xi)/xi -
   (((ee251 - ee72)/xi - ee202)/ee19 - ee184 * ee7) * ee2/ee1) * ee8 +
   ee124 * ee58 * ee50/ee73 + ee58 * (2 * ee124 + 2 * (R_pow(ee58, 2) * ee50/ee122)) * ee50/ee73;
out(j, 17) = ((ee124 * ee67 + ee58 * (2 * (ee58 * ee67 * ee50/ee122) +
   2 * ee135)) * ee50/ee73 + (((ee218 * ee75 - ee200 * R_pow(ee28, 2)/ee114)/ee8 +
   ee116 + 2 * (ee242 * ee43))/ee114 +
   ee149 * ee129) * ee8) * ee10/ee11;
out(j, 18) = ((ee229 * ee58 + 2 * (ee135 * ee67 * ee10/ee11)) * ee50/ee73 -
   ((ee213/ee8 + ee118 - 2 * (ee75 * ee10/ee33)) * ee28/ee114 +
   ((((ee197 - ee181)/ee11 + 2) * ee10 + ee87)/ee11 -
     1) * ee43) * ee8) * ee10/ee11;
out(j, 19) =  - ((((ee126 - ee195)/ee11 + 2 * ee169 - ee118) * ee67 * ee50 * ee10/(ee49 * R_pow(ee11, 3) * ee51) -
   ((((ee111/ee15 +
   2 * (ee71/ee83)) * ee10/ee53 + ee244/ee15)/ee8 -
   (ee111 * ee247 + ee71 * (1/ee45 + ee180)) * ee10/ee11)/ee11 -
   ee244/ee8) * ee8) * ee10/ee11);
    
}

return out;

}
