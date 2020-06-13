// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// start with the integral bit

// [[Rcpp::export]]
double pp1d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, arma::vec yvec, arma::vec wvec)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();

double y, mu, lpsi, xi;
double ee1;
double nllh=0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];

ee1 = xi * (y - mu) / exp(lpsi);

if (ee1 > -1.0) {
  nllh += wvec[j] * R_pow(1.0 + ee1, -1/xi);
}

}

return(nllh);

}

// [[Rcpp::export]]
arma::mat pp1d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec, arma::vec wvec)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 9);
out.zeros();

double y, w, mu, lpsi, xi;
double ee1, ee2, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee14, ee15, ee16, ee19;
double ee22, ee25;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
w = wvec[j];

ee1 = exp(lpsi);
ee2 = y - mu;
ee4 = xi * ee2/ee1;

if (ee4 > -1.0) {

ee5 = 1/xi;
ee6 = 1 + ee4;
ee7 = 1 + ee5;
ee8 = R_pow(ee6, ee7);
ee9 = log1p(ee4);
ee10 = R_pow(ee6, ee5 + 2);
ee11 = ee10 * ee1;
ee12 = ee8 * ee1;
ee14 = ee7 * ee2/ee11;
ee15 = R_pow(ee6, ee5);
ee16 = 1/ee8;
ee19 = ee9/(xi * ee15) - ee2/ee12;
ee22 = w * (ee9/(R_pow(xi, 2) * ee8) - ee14);
ee25 = xi * ee7 * ee2/ee11;

out(j, 0) = w/ee12;
out(j, 1) = w * ee2/ee12;
out(j, 2) = w * ee19/xi;
out(j, 3) = w * xi * ee7/(ee10 * R_pow(ee1, 2));
out(j, 4) = w * (ee25 - ee16)/ee1;
out(j, 5) = ee22/ee1;
out(j, 6) = -(w * (ee16 - ee25) * ee2/ee1);
out(j, 7) = ee22 * ee2/ee1;
out(j, 8) = w * (((ee2/(ee6 * ee1) - 2 * (ee9/xi))/ee15 + ee9 * ee19/xi)/xi +
   (ee14 + (ee16 - ee9/(xi * ee8))/xi) * ee2/ee1)/xi;
            
}

}

return out;

}

// [[Rcpp::export]]
arma::mat pp1d34(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec, arma::vec wvec)
{

arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 25);
out.zeros();

double y, w, mu, lpsi, xi;
double ee1, ee2, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee21, ee22, ee23, ee24, ee26, ee27, ee28, ee29;
double ee30, ee31, ee32, ee33, ee34, ee35, ee36;
double ee40, ee41, ee42, ee43, ee44, ee46, ee47, ee48, ee49;
double ee50, ee52, ee55, ee58;
double ee60, ee61, ee62, ee63, ee64, ee66, ee68, ee69;
double ee70, ee74, ee75, ee78;
double ee80, ee81, ee86, ee87, ee88, ee89;
double ee90, ee92, ee93, ee95, ee99;
double ee101, ee102, ee106, ee107, ee108, ee109;
double ee110, ee111, ee114, ee117;
double ee120;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
w = wvec[j];

ee1 = exp(lpsi);
ee2 = y - mu;
ee4 = xi * ee2/ee1;

if (ee4 > -1.0) {

ee5 = 1/xi;
ee6 = 1 + ee4;
ee7 = ee5 + 2;
ee8 = log1p(ee4);
ee9 = 1 + ee5;
ee10 = R_pow(ee6, ee7);
ee11 = R_pow(xi, 2);
ee12 = ee5 + 3;
ee13 = R_pow(ee6, ee12);
ee14 = R_pow(ee6, ee9);
ee15 = ee6 * ee1;
ee16 = ee2/ee15;
ee17 = ee13 * ee1;
ee18 = ee8/xi;
ee19 = ee11 * ee10;
ee21 = ee7 * ee2/ee17;
ee22 = ee16 - 2 * ee18;
ee23 = ee8/ee19;
ee24 = ee23 - ee21;
ee26 = ee9 * ee2/(ee10 * ee1);
ee27 = ee8/(ee11 * ee14);
ee28 = ee11 * ee13;
ee29 = R_pow(ee6, ee5 + 4);
ee30 = 1/ee10;
ee31 = ee29 * ee1;
ee32 = R_pow(ee6, ee5);
ee33 = ee27 - ee26;
ee34 = ee22/ee14;
ee35 = ee9 * ee24;
ee36 = 2/ee10;
ee40 = ee8/ee28 - ee12 * ee2/ee31;
ee41 = 2 * ee16;
ee42 = ee8 * ee33;
ee43 = xi * ee7;
ee44 = ee34 + ee42;
ee46 = (1/ee19 - ee35) * ee2/ee1;
ee47 = ee44/ee11;
ee48 = 1/ee14;
ee49 = 2/ee13;
ee50 = 2/xi;
ee52 = (ee22/ee10 + ee8 * ee24)/ee11 + (1/ee28 - ee7 * ee40) *  ee2/ee1;
ee55 = (ee16 + ee50) * ee2/ee15 + (ee41 - 6 * ee18)/xi;
ee58 = 2 * ee21;
ee60 = ee8/(xi * ee32) - ee2/(ee14 * ee1);
ee61 = ee47 + ee46;
ee62 = R_pow(ee1, 2);
ee63 = ee8/(xi * ee14);
ee64 = xi * ee9;
ee66 = ee9 * (ee30 + xi * ee24) + ee30;
ee68 = ee9 * (ee49 + ee43 * ee40) + ee7/ee13;
ee69 = 1 + ee9;
ee70 = ee30 + ee36;
ee74 = xi * ee12 * ee2/ee31;
ee75 = ee61 * ee8;
ee78 = (((2 * (ee8/(xi * ee10)) - ee36)/xi - ee58)/ee11 -  ee52 * ee9) * ee2/ee1;
ee80 = (ee22/ee32 + ee8 * ee60/xi)/xi + (ee26 + (ee48 -  ee63)/xi) * ee2/ee1;
ee81 = ee55/ee14;
ee86 = 4/ee10;
ee87 = 4/ee13;
ee88 = R_pow(ee1, 3);
ee89 = w * xi;
ee90 = xi * ee52;
ee92 = ee43 * ee2/ee17;
ee93 = xi * (2 * ee35 - ee68 * ee2/ee1);
ee95 = ((ee80 * ee8 + 2 * (ee22 * ee60))/xi - ee55/ee32)/xi +  ((((ee26 + (2/ee14 - ee63)/xi) * ee8 - (ee34 + 2/ee14))/xi -  2 * ee26)/xi - ee46) * ee2/ee1;
ee99 = (ee66 + ee36 + ee93) * ee2/ee1;
ee101 = ee66 * ee2/ee1;
ee102 = ee9 * (ee23 + ee90 - ee21);
ee106 = 1/ee13;
ee107 = ee41 + 6/xi;
ee108 = 2 * (ee8/ee10);
ee109 = w * ((ee75 + 2 * (ee22 * ee33) - ee81)/ee11 + ee78);
ee110 = w * ee61;
ee111 = ee89 * ee9;
ee114 = w * ee11 * ee69 * ee9;
ee117 = ee64 * (ee70 - ee92) * ee2/ee1;
ee120 = ee64 * (ee70 + ee86 - ee43 * (ee49 + ee87 - ee74) *  ee2/ee1) * ee2/ee1;

out(j, 0) = ee114/(ee13 * ee88);
out(j, 1) = -(ee111 * (ee36 - ee92)/ee62);
out(j, 2) = w * (ee30 + ee64 * ee24)/ee62;
out(j, 3) = w * (ee48 - ee117)/ee1;
out(j, 4) = w * (ee101 - ee27)/ee1;
out(j, 5) = ee110/ee1;
out(j, 6) = -(w * (ee117 - ee48) * ee2/ee1);
out(j, 7) = -(w * (ee27 - ee101) * ee2/ee1);
out(j, 8) = ee110 * ee2/ee1;
out(j, 9) = w * ee95/xi;
out(j, 10) = w * R_pow(xi, 3) * (1 + ee69) * ee69 * ee9/(ee29 * R_pow(ee1, 4));
out(j, 11) = -(ee114 * (ee106 + ee49 - ee74)/ee88);
out(j, 12) = ee89 * ee68/ee88;
out(j, 13) = -(ee111 * (xi * (ee106 + ee87 - ee74) * ee7 * ee2/ee1 -
   ee86)/ee62);
out(j, 14) = -(w * (ee36 + ee93)/ee62);
out(j, 15) = w * (2 * ee23 + ee90 * ee9 - ee58)/ee62;
out(j, 16) = w * (ee120 - ee48)/ee1;
out(j, 17) = w * (ee27 - ee99)/ee1;
out(j, 18) = w * ((ee102 + (ee108 - ee30)/ee11 - ee58) * ee2/ee1 -
   ee47)/ee1;
out(j, 19) = ee109/ee1;
out(j, 20) = -(w * (ee48 - ee120) * ee2/ee1);
out(j, 21) = -(w * (ee99 - ee27) * ee2/ee1);
out(j, 22) = -(w * (ee47 + ((ee30 - ee108)/ee11 + ee58 - ee102) * ee2/ee1) * ee2/ee1);
out(j, 23) = ee109 * ee2/ee1;
out(j, 24) = w * (((ee95 * ee8 + 3 * (ee80 * ee22) - 3 * (ee55 * ee60))/xi -
   (((24 * ee18 - 6 * ee16)/xi - ee107 * ee2/ee15)/xi -
     ((ee41 + ee50) * ee2/ee15 + ee107/xi) * ee2/ee15)/ee32)/xi +
   ((((ee34 + 2 * ee44 + 6/ee14 + ee42)/xi + ee81 -
     (ee75 + (2 * ee22 + 6) * ee33))/xi + 3 * ee46)/xi -   ee78) * ee2/ee1)/xi;

}

}

return out;

}

// the density bit

// [[Rcpp::export]]
double pp2d0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, arma::vec yvec, arma::vec wvec)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();

double y, mu, lpsi, xi;
double ee1;
double nllh=0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];

ee1 = xi * (y - mu) / exp(lpsi);

if (ee1 <= -1.0) {
    nllh = 1e20;
    break;
} else {

nllh += wvec[j] * (lpsi + (1 / xi + 1) * log1p(ee1));
    
}
 
}

return(nllh);

}

// [[Rcpp::export]]
arma::mat pp2d12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 9);

double y, mu, lpsi, xi;
double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee14, ee16, ee17;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];

ee1 = exp(lpsi);
ee2 = y - mu;
ee3 = xi * ee2;
ee4 = ee3/ee1;
ee5 = 1 + ee4;
ee6 = ee5 * ee1;
ee7 = 1/xi;
ee8 = 1 + ee7;
ee9 = ee3/ee6;
ee10 = R_pow(xi, 2);
ee11 = 1 - ee9;
ee12 = xi * ee8;
ee14 = ee11 * ee8 - ee7;
ee16 = ee8 * ee2/ee6;
ee17 = log1p(ee4);

out(j, 0) = -(ee12/ee6);
out(j, 1) = 1 - ee12 * ee2/ee6;
out(j, 2) = ee16 - ee17/ee10;
out(j, 3) =  - (ee10 * ee8/(R_pow(ee5, 2) * R_pow(ee1, 2)));
out(j, 4) = xi * ee11 * ee8/ee6;
out(j, 5) = -(ee14/ee6);
out(j, 6) = -(ee12 * (ee9 - 1) * ee2/ee6);
out(j, 7) = -(ee14 * ee2/ee6);
out(j, 8) = -((ee16 + 1/ee10) * ee2/ee6 + (ee2/ee6 - 2 * (ee17/xi))/ee10);

}

return out;

}

// [[Rcpp::export]]
arma::mat pp2d34(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec)
{

arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 25);

double y, mu, lpsi, xi;
double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee15, ee17, ee18, ee19;
double ee20, ee22, ee25, ee26, ee27, ee29;
double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee37;
double ee40, ee42, ee45, ee47, ee49;
double ee50, ee53, ee55, ee56, ee57, ee59;
double ee62, ee63, ee64, ee66;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];

ee1 = exp(lpsi);
ee2 = y - mu;
ee3 = xi * ee2;
ee4 = ee3/ee1;
ee5 = 1 + ee4;
ee6 = ee5 * ee1;
ee7 = 1 + 1/xi;
ee8 = ee3/ee6;
ee9 = 2 * ee8;
ee10 = R_pow(xi, 2);
ee11 = ee2/ee6;
ee12 = 4 * ee5;
ee13 = 2 - ee9;
ee15 = 2 * ee11;
ee17 = 2 * (1 + 2 * ee4) + ee12;
ee18 = 8 * ee4;
ee19 = 2/xi;
ee20 = 1 - ee8;
ee22 = xi * ee7;
ee25 = xi * (3 - ee9) * ee2/ee6;
ee26 = ee17 - ee18;
ee27 = 6 * ee4;
ee29 = R_pow(ee5, 2) * R_pow(ee1, 2);
ee30 = ee20/xi;
ee31 = ee26/ee5;
ee32 = (ee12 - ee27)/ee5;
ee33 = (ee18 - ee17)/ee5;
ee34 = ee15 + ee19;
ee35 = ee15 + 6/xi;
ee36 = R_pow(ee1, 3);
ee37 = ee25 - 1;
ee40 = (((2 * ee5 - ee27)/ee5 + 4) * ee7 * ee2/ee6 + ee13/ee10) *  ee2/ee6 + ((4 * ee20 - 6)/xi + 2 * ee30 + 2 * ee13 *  ee2/ee6)/ee10;
ee42 = (ee7 * ee13 - ee19) * ee2/ee6;
ee45 = ee7 * (4 - xi * (ee32 + 6) * ee2/ee6) * ee2/ee6;
ee47 = ee7 * ee2/ee6;
ee49 = R_pow(ee5, 3) * ee36;
ee50 = R_pow(ee5, 4);
ee53 = ee34 * ee2/ee6 + ee35/xi;
ee55 = ee11 + ee19;
ee56 = 1 + xi * (ee9 - 3) * ee2/ee6;
ee57 = 2 * ee25;
ee59 = log1p(ee4)/xi;
ee62 = xi * (ee31 + 6) * ee2/ee6;
ee63 = ee8 - 1;
ee64 = ee10 * ee7;
ee66 = R_pow(xi, 3) * ee7;

out(j, 0) = -(2 * (ee66/ee49));
out(j, 1) = ee64 * ee13/ee29;
out(j, 2) = -((ee22 * ee13 - 1)/ee29);
out(j, 3) = ee22 * ee37/ee6;
out(j, 4) = ((1 - ee25) * ee7 - ee30)/ee6;
out(j, 5) = ee42/ee6;
out(j, 6) = -(ee22 * ee56 * ee2/ee6);
out(j, 7) = -((ee7 * ee37 - ee63/xi) * ee2/ee6);
out(j, 8) = ee42 * ee2/ee6;
out(j, 9) = (ee55 * ee2/ee6 + (ee15 - 6 * ee59)/xi)/ee10 + (ee55/ee10 +
   (1/ee10 + 2 * ee47) * ee2/ee6) * ee2/ee6;
out(j, 10) =  - (6 * (R_pow(xi, 4) * ee7/(ee50 * R_pow(ee1, 4))));
out(j, 11) = ee66 * ee26/(ee50 * ee36);
out(j, 12) = -(xi * (ee22 * ee26/ee5 - 2)/ee49);
out(j, 13) = ee64 * (xi * (4 - ee33) * ee2/ee6 - 4)/ee29;
out(j, 14) = (xi * (ee7 * (4 - xi * (ee31 + 4) * ee2/ee6) +
   ee15) - 2)/ee29;
out(j, 15) = -((ee7 * (2 - xi * (ee32 + 4) * ee2/ee6) + 4 * ee2 -
   ee13/xi)/ee6/ee6);
out(j, 16) = ee22 * (1 + xi * (xi * (6 - ee33) * ee2/ee6 - 7) * ee2/ee6)/ee6;
out(j, 17) = (ee7 * (xi * (7 - ee62) * ee2/ee6 - 1) - ee37/xi)/ee6;
out(j, 18) = -((ee45 + (2 - (2 * ee20 + ee57))/ee10)/ee6);
out(j, 19) = -(ee40/ee6);
out(j, 20) = -(ee22 * (xi * (7 + xi * (ee33 - 6) * ee2/ee6) * ee2/ee6 -
   1) * ee2/ee6);
out(j, 21) = -((ee7 * (1 + xi * (ee62 - 7) * ee2/ee6) - ee56/xi) * ee2/ee6);
out(j, 22) = -((ee45 - (ee57 - (2 + 2 * ee63))/ee10) * ee2/ee6);
out(j, 23) = -(ee40 * ee2/ee6);
out(j, 24) = (((24 * ee59 - 6 * ee11)/xi - ee35 * ee2/ee6)/xi -
   ee53 * ee2/ee6)/ee10 - (ee53/ee10 + (ee34/ee10 + (2/ee10 +
   6 * ee47) * ee2/ee6) * ee2/ee6) * ee2/ee6;

}

return out;

}

// [[Rcpp::export]]
double ppcd0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, arma::vec yvec, arma::vec wvec)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();

double y, mu, lpsi, xi;
double ee1;
double nllh=0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];

ee1 = xi * (y - mu) / exp(lpsi);

if (ee1 <= -1.0) {
    nllh = 1e20;
    break;
} else {

nllh += (wvec[j] / xi) * log1p(ee1);
    
}

}

return(nllh);

}

// [[Rcpp::export]]
arma::mat ppcd12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 9);

double y, mu, lpsi, xi;
double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee10; 
double ee12; 

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];

ee1 = exp(lpsi);
ee2 = y - mu;
ee3 = xi * ee2;
ee4 = ee3/ee1;
ee5 = 1 + ee4;
ee6 = ee5 * ee1;
ee7 = ee2/ee6;
ee8 = ee3/ee6;
ee10 = log1p(ee4)/xi;
ee12 = xi * ee5 * ee1;

out(j, 0) = -(1/ee6);
out(j, 1) = -ee7;
out(j, 2) = (ee7 - ee10)/xi;

out(j, 3) = -(xi/(ee5 * ee5 * ee1 * ee1));
out(j, 4) = (1 - ee8)/ee6;
out(j, 5) = ee8/ee12;
out(j, 6) = -((ee8 - 1) * ee2/ee6);
out(j, 7) = ee8 * ee2/ee12;
out(j, 8) = -(((ee7 - 2 * ee10)/xi + (ee7 + 1/xi) * ee2/ee6)/xi);

}

return out;

}

// [[Rcpp::export]]
arma::mat ppcd34(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec)
{

arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 25);

double y, mu, lpsi, xi;
double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9; 
double ee10, ee12, ee14, ee16, ee17, ee19; 
double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29; 
double ee31, ee34, ee37, ee39; 
double ee40, ee43, ee46, ee47, ee49; 
double ee50, ee51, ee53, ee56, ee57, ee58; 

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];

ee1 = exp(lpsi);
ee2 = y - mu;
ee3 = xi * ee2;
ee4 = ee3/ee1;
ee5 = 1 + ee4;
ee6 = ee5 * ee1;
ee7 = ee3/ee6;
ee8 = ee2/ee6;
ee9 = 2 * ee7;
ee10 = 4 * ee5;
ee12 = xi * ee5 * ee1;
ee14 = 2 * ee8;
ee16 = 2 * (1 + 2 * ee4) + ee10;
ee17 = 8 * ee4;
ee19 = ee5 * ee5 * ee1 * ee1;
ee20 = 2 - ee9;
ee21 = 2/xi;
ee22 = 1 - ee7;
ee23 = ee16 - ee17;
ee24 = 6 * ee4;
ee25 = ee23/ee5;
ee26 = (ee10 - ee24)/ee5;
ee27 = (ee17 - ee16)/ee5;
ee28 = ee14 + ee21;
ee29 = ee14 + 6/xi;
ee31 = R_pow(ee1, 3);
ee34 = xi * (3 - ee9) * ee2/ee6;
ee37 = (((2 * ee5 - ee24)/ee5 + 4) * ee2/ee6 + ee20/xi) * 
  ee2/ee6 + ((4 * ee22 - 6)/xi + 2 * (ee22/xi) + 2 * ee20 * 
  ee2/ee6)/xi;
ee39 = R_pow(ee5, 3) * ee31;
ee40 = R_pow(ee5, 4);
ee43 = ee28 * ee2/ee6 + ee29/xi;
ee46 = (4 - xi * (ee26 + 6) * ee2/ee6) * ee2/ee6;
ee47 = ee8 + ee21;
ee49 = 2 * ee34;
ee50 = 2 * (xi * ee2 * ee2/ee19);
ee51 = 6 * ee8;
ee53 = log1p(ee4)/xi;
ee56 = xi * (ee25 + 4) * ee2/ee6;
ee57 = xi * ee20;
ee58 = xi * xi;

// third derivatives
// 1=location, 2=log(scale), 3=shape
// order: 111, 112, 113, 122, 123, 133, 222, 223, 233, 333

out(j, 0) = -(2 * (ee58/ee39));
out(j, 1) = ee57/ee19;
out(j, 2) = -((1 - ee9)/ee19);
out(j, 3) = (ee34 - 1)/ee6;
out(j, 4) = xi * (ee9 - 2) * ee2/ee6/ee12;
out(j, 5) = -(ee50/ee12);
out(j, 6) = -((1 + xi * (ee9 - 3) * ee2/ee6) * ee2/ee6);
out(j, 7) = -(ee57 * ee2/ee6 * ee2/ee12);
out(j, 8) = -(ee50 * ee2/ee12);
out(j, 9) = ((ee47 * ee2/ee6 + (ee14 - 6 * ee53)/xi)/xi + 
  (ee47/xi + (1/xi + ee14) * ee2/ee6) * ee2/ee6)/xi;

// fourth derivatives
// 1=location, 2=log(scale), 3=shape
// order: 1111, 1112, 1113, 1122, 1123, 1133, 1222, 1223, 1233, 1333
//  2222, 2223, 2233, 2333, 3333

out(j, 10) =  -(6 * (R_pow(xi, 3)/(ee40 * R_pow(ee1, 4))));
out(j, 11) = ee58 * ee23/(ee40 * ee31);
out(j, 12) = -(xi * (ee25 - 2)/ee39);
out(j, 13) = xi * (xi * (4 - ee27) * ee2/ee6 - 4)/ee19;
out(j, 14) = (2 - xi * (ee25 + 2) * ee2/ee6)/ee19;
out(j, 15) = -((4 - 2 * (ee26 + 2)) * ee2/ee6/ee6);
out(j, 16) = (1 + xi * (xi * (6 - ee27) * ee2/ee6 - 7) * ee2/ee6)/ee6;
out(j, 17) = xi * (4 - ee56) * ee2/ee6/ee12;
out(j, 18) = -(((2 - (2 * ee22 + ee49))/xi + ee46)/ee12);
out(j, 19) = -(ee37/ee12);
out(j, 20) = -((xi * (7 + xi * (ee27 - 6) * ee2/ee6) * ee2/ee6 - 1) * ee2/ee6);
out(j, 21) = -(xi * (ee56 - 4) * ee2/ee6 * ee2/ee12);
out(j, 22) = -((ee46 - (ee49 - (2 + 2 * (ee7 - 1)))/xi) * ee2/ee12);
out(j, 23) = -(ee37 * ee2/ee12);
out(j, 24) = ((((24 * ee53 - ee51)/xi - 
  ee29 * ee2/ee6)/xi - ee43 * ee2/ee6)/xi - (ee43/xi + 
  (ee28/xi + (ee21 + ee51) * ee2/ee6) * ee2/ee6) * 
  ee2/ee6)/xi;

}

return out;

}
