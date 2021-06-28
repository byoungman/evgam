// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double xieps = 0.0001;

// [[Rcpp::export]]
arma::mat mean_arr(arma::mat x, double n)
{
    
int nr = x.n_rows;
int nc = x.n_cols;
int nr2 = nr / n;

arma::mat x2(nc, nr2);
arma::mat xi;

arma::rowvec pre(n);
pre.fill(1.0 / n);

for (int i=0; i < nc; i++) {
    
xi = x.col(i);
xi.reshape(n, nr2);
x2.row(i) = pre * xi;

}

return x2.t();

}

// [[Rcpp::export]]
arma::vec ragged_mean_vec(arma::vec x, arma::uvec n)
{

int m = n.size();    
arma::vec out(m, arma::fill::zeros);
int i = 0;
int nj;

for (int j = 0; j < m; j++) {
  nj = n[j];
  for (int k = 0; k < nj; k++) {
    out[j] += x[i];
    i += 1;
  }
  out[j] /= nj;
}

return out;

}

// [[Rcpp::export]]
arma::mat ragged_mean_mat(arma::mat x, arma::uvec n)
{

int m = n.size();
int p = x.n_cols;
arma::mat out(m, p, arma::fill::zeros);
int i, nj;

for (int l = 0; l < p; l++) {
  i = 0;
  for (int j = 0; j < m; j++) {
    nj = n[j];
    for (int k = 0; k < nj; k++) {
      out(j, l) += x(i, l);
      i += 1;
    }
    out(j, l) /= nj;
  }
}

return out;

}


// // [[Rcpp::export]]
// arma::vec gevd0(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec)
// {
//     
// int nobs = yvec.size();
// arma::vec nllh(nobs);
// 
// double y, mu, lpsi, xi;
// double ee1, ee2;
// 
// for (int j=0; j < nobs; j++) {
// 
// y = yvec[j];
// mu = muvec[j];
// lpsi = lpsivec[j];
// xi = xivec[j];
// 
// if (fabs(xi) >= xieps) {
// 
// ee1 = xi * (y - mu) / exp(lpsi);
// 
// if (ee1 <= -1.0) {
//     nllh[j] = 1e20;
//     break;
// } else {
// 
// ee2 = 1.0 / xi;
// 
// nllh[j] = lpsi + (ee2 + 1.0) * log1p(ee1) + R_pow(1.0 + ee1, -ee2);
// }
// 
// } else {
// 
// ee1 = (y - mu) / exp(lpsi);
// nllh[j] = lpsi + ee1 + exp(-ee1);
//     
// }
// 
// }
// 
// return(nllh);
// 
// }

// [[Rcpp::export]]
double ldgev(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec)
{
    
int nobs = yvec.size();

double y, mu, lpsi, xi;
double ee1, ee2;
double nllh = 0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];

if (fabs(xi) >= xieps) {

ee1 = xi * (y - mu) / exp(lpsi);

if (ee1 <= -1.0) {
    nllh = 1e20;
    break;
} else {

ee2 = 1.0 / xi;

nllh += lpsi + (ee2 + 1.0) * log1p(ee1) + R_pow(1.0 + ee1, -ee2);
}

} else {

ee1 = (y - mu) / exp(lpsi);
nllh += lpsi + ee1 + exp(-ee1);
    
}

}

return nllh;

}

// // [[Rcpp::export]]
// arma::mat ldgev1(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec)
// {
//     
// int nobs = yvec.size();
// arma::mat out(nobs, 4);
// 
// double y, mu, lpsi, xi;
// double ee1, ee2, ee4, ee5, ee6, ee8, ee9, ee10, ee11, ee12;
// 
// for (int j=0; j < nobs; j++) {
// 
// y = yvec[j];
// mu = muvec[j];
// lpsi = lpsivec[j];
// xi = xivec[j];
// 
// if (fabs(xi) >= xieps) {
// 
// ee1 = exp(lpsi);
// ee2 = y - mu;
// ee4 = xi * ee2/ee1;
// ee5 = 1 + ee4;
// ee6 = 1/xi;
// ee8 = (1 - 1/R_pow(ee5, ee6))/xi + 1;
// ee9 = ee5 * ee1;
// ee10 = log1p(ee4);
// ee11 = xi * ee8;
// ee12 = xi * xi;
// 
// out(j, 0) = -(ee11/ee9);
// out(j, 1) = 1 - ee11 * ee2/ee9;
// out(j, 2) = (ee8 * ee2/ee1 + ee10/(ee12 * R_pow(ee5, ee6 - 1)))/ee5 - ee10/ee12;
// out(j, 3) = 0;
// 
// } else {
//     
// ee1 = exp(lpsi);
// ee2 = y - mu;
// ee4 = 1 - exp(-(ee2/ee1));
// 
// out(j, 0) = -(ee4/ee1);
// out(j, 1) = 1 - ee4 * ee2/ee1;
// out(j, 2) = 0;
// out(j, 3) = 0;
//     
// }
// 
// }
// 
// return(out);
// 
// }

// [[Rcpp::export]]
arma::mat ldgev12(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec)
{
    
int nobs = yvec.size();
arma::mat out(nobs, 14);

double y, mu, lpsi, xi;
double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee16, ee17, ee18, ee19;
double ee20, ee22, ee23;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];

if (fabs(xi) >= xieps) {

ee1 = exp(lpsi);
ee2 = y - mu;
ee4 = xi * ee2/ee1;
ee5 = 1 + ee4;
ee6 = 1/xi;
ee7 = 1 + ee6;
ee8 = R_pow(ee5, ee6);
ee9 = ee5 * ee1;
ee10 = 1/ee8;
ee11 = log1p(ee4);
ee12 = R_pow(ee5, ee7);
ee13 = ee7 * ee2;
ee16 = (ee10 - xi) * ee2/ee9 + 1;
ee17 = ee11/(xi * ee8);
ee18 = xi * ee7;
ee19 = (ee16 * ee7 - (1 + ee17)/xi)/ee5;
ee20 = ee13/ee9;
ee22 = ee2/(ee12 * ee1);
ee23 = xi - ee10;

out(j, 0) = -((ee18 - ee10)/ee9);
out(j, 1) = (ee10 - ee18) * ee2/ee9 + 1;
out(j, 2) = ((ee10 - 1) * ee11/xi - ee22)/xi + ee20;
out(j, 3) = 0;
out(j, 4) =  - (ee18 * ee23/(R_pow(ee5, 2) * R_pow(ee1, 2)));
out(j, 5) = (xi * ee16 * ee7 - ee10)/ee5/ee1;
out(j, 6) = -(ee19/ee1);
out(j, 7) = 0;
out(j, 8) = -((ee10 + xi * (ee23 * ee2/ee9 - 1) * ee7)/ee5 * ee2/ee1);
out(j, 9) = -(ee19 * ee2/ee1);
out(j, 10) = 0;
out(j, 11) = ((((ee2/ee9 - 2 * (ee11/xi))/R_pow(ee5, ee6 - 1) -
   ee2/ee1)/ee5 + (2 + ee17 - ee22) * ee11/xi)/xi + (ee13/(R_pow(ee5, ee6 +
   2) * ee1) + (1/ee12 - ee11/(xi * ee12))/xi) * ee2/ee1)/xi -
   (ee20 + 1/R_pow(xi, 2)) * ee2/ee9;
out(j, 12) = 0;
out(j, 13) = 0;

} else {
    
ee1 = exp(lpsi);
ee2 = y - mu;
ee3 = ee2/ee1;
ee5 = exp(-ee3);
ee7 = (ee3 - 1) * ee5 + 1;
ee8 = ee5 - 1;

out(j, 0) = ee8/ee1;
out(j, 1) = ee8 * ee2/ee1 + 1;
out(j, 2) = 0;
out(j, 3) = 0;
out(j, 4) = ee5/R_pow(ee1, 2);
out(j, 5) = ee7/ee1;
out(j, 6) = 0;
out(j, 7) = 0;
out(j, 8) = ee7 * ee2/ee1;
out(j, 9) = 0;
out(j, 10) = 0;
out(j, 11) = 0;
out(j, 12) = 0;
out(j, 13) = 0;
    
}

}

return(out);

}

// [[Rcpp::export]]
double ldgevagg(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec, arma::vec thetavec)
{
    
int nobs = yvec.size();

double y, mu, lpsi, psi, xi, theta;
double ee1, ee2;
double nllh = 0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
theta = thetavec[j];
psi = exp(lpsi);

if (fabs(xi) >= xieps) {

mu = mu - psi * (1 - exp(theta * xi)) / xi;
lpsi = lpsi + xi * theta;

ee1 = xi * (y - mu) / exp(lpsi);

if (ee1 <= -1.0) {
    nllh = 1e20;
    break;
} else {

ee2 = 1.0 / xi;

nllh += lpsi + (ee2 + 1.0) * log1p(ee1) + R_pow(1.0 + ee1, -ee2);
}

} else {

mu = mu + psi * theta;

ee1 = (y - mu) / exp(lpsi);
nllh += lpsi + ee1 + exp(-ee1);
    
}

}

return nllh;

}

// // [[Rcpp::export]]
// arma::mat ldgevagg1(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec, arma::vec thetavec)
// {
//     
// int nobs = yvec.size();
// arma::mat out(nobs, 4);
// 
// double y, mu, lpsi, xi, theta;
// double ee1, ee2, ee3, ee4, ee7, ee8, ee9;
// double ee10, ee11, ee12, ee13, ee14, ee16, ee17, ee20, ee21, ee22;
// 
// for (int j=0; j < nobs; j++) {
// 
// y = yvec[j];
// mu = muvec[j];
// lpsi = lpsivec[j];
// xi = xivec[j];
// theta = thetavec[j];
// 
// if (fabs(xi) >= xieps) {
// 
// ee1 = log(theta);
// ee2 = R_pow(theta, xi);
// ee3 = exp(lpsi);
// ee4 = 1 - ee2;
// ee7 = ee4 * ee3/xi + y;
// ee8 = exp(lpsi + xi * ee1);
// ee9 = ee7 - mu;
// ee10 = xi * ee9;
// ee11 = ee10/ee8;
// ee12 = 1 + ee11;
// ee13 = 1/xi;
// ee14 = R_pow(ee12, ee13);
// ee16 = (1 - 1/ee14)/xi + 1;
// ee17 = ee12 * ee8;
// ee20 = ee7 - (mu + xi * (ee9 * ee1 + (ee4/xi + ee2 * ee1) * ee3/xi));
// ee21 = log1p(ee11);
// ee22 = xi * ee16;
// 
// out(j, 0) = -(ee22/ee17);
// out(j, 1) = 1 - ee22 * (y - mu)/ee17;
// out(j, 2) = (ee20 * (1 + ee13)/ee8 - (ee20 * R_pow(ee12, ee13 - 1)/ee8 - ee14 * ee21/xi)/(xi * R_pow(ee12, 2/xi - 1)))/ee12 + ee1 - ee21/ (xi * xi);
// out(j, 3) = xi * (1/theta - ee16 * (R_pow(theta, xi - 1) * ee3 + ee10/theta)/ee17);
// 
// } else {
// 
// ee1 = exp(lpsi);
// ee2 = y - mu;
// ee4 = 1 - exp(-(ee2/ee1));
// 
// out(j, 0) = -(ee4/ee1);
// out(j, 1) = 1 - ee4 * ee2/ee1;
// out(j, 2) = 0;
// out(j, 3) = 0;
// 
//     
// }
// 
// }
// 
// return(out);
// 
// }

// [[Rcpp::export]]
arma::mat ldgevagg12(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec, arma::vec thetavec)
{
    
int nobs = yvec.size();
arma::mat out(nobs, 14);

double y, mu, lpsi, xi, theta;
// double ee1, ee2, ee3, ee4, ee5, ee6, ee8, ee9;
// double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
// double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
// double ee30, ee31, ee33, ee34, ee35, ee36, ee37, ee38, ee39;
// double ee40, ee41, ee42, ee43, ee44, ee45, ee46, ee48, ee49;
// double ee50, ee51, ee52, ee53, ee54, ee55, ee56, ee57, ee58;
// double ee60, ee61;
// double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
// double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
// double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
// double ee30, ee31, ee32, ee33, ee34, ee35, ee37, ee38, ee39;
// double ee41, ee42, ee43, ee44, ee45, ee46, ee47, ee48;
// double ee51, ee52, ee53, ee55, ee56, ee57, ee59;
// double ee60, ee63, ee64, ee67, ee68, ee69;
// double ee70, ee71, ee72;
double ee1, ee2, ee3, ee4, ee5, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee18, ee19;
double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee28, ee29;
double ee30, ee31, ee32, ee33, ee34, ee35, ee37, ee39;
double ee41, ee42, ee43, ee44, ee49;
double ee51, ee52, ee53;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
theta = thetavec[j];

if (fabs(xi) >= xieps) {

// ee1 = R_pow(theta, xi);
// ee2 = exp(lpsi);
// ee3 = (1 - ee1) * ee2;
// ee6 = ee3/xi + y - mu;
// ee7 = xi * ee6;
// ee8 = ee1 * ee2;
// ee9 = log(theta);
// ee10 = ee7/ee8;
// ee11 = 1 + ee10;
// ee12 = 1/xi;
// ee13 = 1 + ee12;
// ee15 = y - (mu + ee8 * ee9);
// ee16 = xi - 1;
// ee17 = ee7 * ee9;
// ee18 = ee15/ee1;
// ee19 = ee17/ee1;
// ee20 = 2 * xi;
// ee21 = R_pow(ee11, ee13);
// ee22 = ee18 - ee19;
// ee23 = R_pow(ee11, ee12);
// ee24 = R_pow(theta, ee20);
// ee25 = R_pow(theta, ee16);
// ee26 = 1 - (ee1 + ee7/ee2);
// ee27 = log1p(ee10);
// ee29 = 1/ee1;
// ee30 = R_pow(theta, 1 + xi);
// ee31 = R_pow(xi, 2);
// ee32 = ee22 * ee13;
// ee33 = R_pow(ee11, ee12 + 2);
// ee34 = 2 * ee16;
// ee35 = ee30 * ee21;
// ee37 = ee25 * xi * ee9;
// ee38 = R_pow(theta, 2);
// ee39 = ee1 * xi;
// ee41 = ee11 * ee2;
// ee42 = 1/ee23;
// ee44 = ee9/ee1;
// ee46 = ee27/(ee31 * ee21) - ee32/(ee33 * ee2);
// ee47 = R_pow(theta, ee20 - 1);
// ee48 = (((ee1 * ee15 + ee39 * ee6 * ee9)/ee24 + ee18 - 2 *  ee19)/ee2 + ee9) * ee9;
// ee51 = ee6 * (ee25 + ee37)/ee2 - ee47 * ee9;
// ee52 = (ee15/ee2 + ee1 * ee9)/ee1;
// ee53 = ee22/(ee21 * ee2);
// ee55 = ee22/(ee1 * ee11 * ee2) + ee44;
// ee56 = R_pow(ee26, 2);
// ee57 = ee3 - ee7;
// ee59 = ee13/ee35 - 1/ee35;
// ee60 = R_pow(ee11, ee12 - 1);
// ee63 = ee42 - 1;
// ee64 = 1/(ee1 * ee23);
// ee67 = 2 * ee25 + ee37 - ee25;
// ee68 = R_pow(ee2, 2);
// ee69 = ee27/(xi * ee23);
// ee70 = ee24 * ee11;
// ee71 = xi * ((R_pow(theta, xi - 2) * ee6 * ee16/ee2 - R_pow(theta, ee34))/ee24 - (2 * (ee7/(R_pow(theta, 2 + xi) * ee2)) + R_pow(theta, ee34 - ee20)));
// ee72 = xi * ee13;
// 
// out(j, 0) = -((ee72/ee1 - ee64)/ee41);
// out(j, 1) = (ee13/ee1 - 1/(ee39 * ee23)) * ee26/ee11 + 1;
// out(j, 2) = (ee63 * ee27/xi - ee53)/xi + ee32/ee41 + ee9;
// out(j, 3) = ee63/theta;
// out(j, 4) =  - (ee72 * (xi/ee24 - 1/(ee24 * ee23))/(R_pow(ee11, 2) * ee68));
// out(j, 5) = -(((ee26/(ee24 * ee21) - xi * (ee26/ee70 + ee29)) * ee13 +
//    ee64)/ee11/ee2);
// out(j, 6) = -(((((ee29 + xi * ee9/ee1 - ee29)/ee23 - ee29)/xi +
//    ee13 * (ee29 - xi * ee55))/ee11 - ee46/ee1)/ee2);
// out(j, 7) = xi * ee59/ee2;
// out(j, 8) = (ee57/(ee1 * ee21 * ee2) + ee56 * ee13/(ee24 * ee33))/xi -
//    (ee57/ee8 + ee56/ee70) * ee13/ee11;
// out(j, 9) = ((ee52 + ee26 * (1/ee39 + ee44))/ee21 - ee26 * ee46/ee1)/xi -
//    ((ee55 * ee26 + ee52) * ee13 + ee26/(ee1 * ee31))/ee11;
// out(j, 10) = -(ee59 * ee26);
// out(j, 11) = ((((ee22/ee41 - 2 * (ee27/xi))/ee60 - ee22/ee2)/ee11 +
//    (2 + ee69 - ee53) * ee27/xi)/xi + (ee48 + ee22/(xi * ee2))/ee21 -
//    ee22 * ee46/ee2)/xi - ((ee48 + R_pow(ee22, 2)/(ee11 * ee68)) * ee13 +
//    ee22/(ee31 * ee2))/ee11;
// out(j, 12) = ((1 + ee69 - ee42)/xi + 1 - (ee22/(ee23 * ee2) +
//    ee67/ee25 + xi * (ee51/ee47 - ee7 * (2/ee1 - ee29) * ee9/ee2)) * ee13/ee11)/theta +
//    (ee67/ee1 + xi * (ee51/ee24 + (ee15/ee30 -
//    2 * (ee17/ee30))/ee2))/(xi * ee21);
// out(j, 13) = -(xi * ((((xi - (1 + 1/ee60))/ee38 + ee71)/ee11 +
//    xi/ee38) * ee13 + 1/ee38) - (ee16/ee38 + ee71)/ee21);

ee1 = theta * xi;
ee2 = exp(lpsi);
ee3 = exp(ee1);
ee4 = (1 - ee3) * ee2;
ee7 = ee4/xi + y - mu;
ee8 = xi * ee7;
ee9 = exp(lpsi + ee1);
ee10 = ee8/ee9;
ee11 = 1 + ee10;
ee12 = 1/xi;
ee13 = 1 + ee12;
ee14 = ee2 * ee3;
ee15 = ee14 + ee8;
ee16 = theta * ee15;
ee18 = y - (mu + ee16);
ee19 = ee11 * ee9;
ee20 = R_pow(ee11, ee13);
ee21 = R_pow(ee11, ee12);
ee22 = ee4 - ee8;
ee23 = ee20 * ee9;
ee24 = 1/ee21;
ee25 = log1p(ee10);
ee26 = ee13 * ee18;
ee28 = R_pow(ee11, ee12 + 2) * ee9;
ee29 = ee15/ee23;
ee30 = R_pow(xi, 2);
ee31 = ee26/ee28;
ee32 = ee15/ee19;
ee33 = ee18/ee19;
ee34 = ee18/ee23;
ee35 = ((ee29 + xi * (1 - ee32)) * ee13 - ee24)/ee11;
ee37 = ee22/ee19 + 1;
ee39 = ee7 * (1 - ee1) + y;
ee41 = (1 + ee1) * ee2 * ee3;
ee42 = ee33 + theta;
ee43 = 1/ee20;
ee44 = 2 * (theta * ee2 * ee3);
ee49 = ee25/(xi * ee21);
ee51 = ee25/(ee30 * ee20) - ee31;
ee52 = theta * (2 * y - (2 * mu + ee16));
ee53 = xi * ee13;

out(j, 0) = -((ee53 - ee24)/ee19);
out(j, 1) = ((1 - ee24)/xi + 1) * ee22/ee19 + 1;
out(j, 2) = ((ee24 - 1) * ee25/xi - ee34)/xi + ee26/ee19 + theta;
out(j, 3) = ee29 + xi * (1 - ee13 * ee15/ee19);
out(j, 4) =  - (ee53 * (xi - ee24)/(R_pow(ee11, 2) * R_pow(ee9, 2)));
out(j, 5) = -(((ee22/ee23 - xi * ee37) * ee13 + ee24)/ee11/ee9);
out(j, 6) = -(((ee34 + 1 - xi * ee42) * ee13 + theta/ee21 -
   (1 + ee49)/xi)/ee11/ee9);
out(j, 7) = xi * ee35/ee9;
out(j, 8) = ((ee22 * ee13/ee28 + ee43)/xi - ee37 * ee13/ee11) * ee22/ee9;
out(j, 9) = (((ee22 * (ee12 + theta) + y - mu)/ee20 - ee22 * ee51)/xi -
   ((ee22 * ee42 + y - mu) * ee13 + ee22/ee30)/ee11)/ee9;
out(j, 10) = -(ee35 * ee22/ee9);
out(j, 11) = ((((ee33 - 2 * (ee25/xi))/R_pow(ee11, ee12 - 1) -
   ee18/ee9)/ee11 + (2 + ee49 - ee34) * ee25/xi)/xi + ((ee18/xi +
   ee52)/ee20 - ee51 * ee18)/ee9)/xi - ((R_pow(ee18, 2)/ee19 +
   ee52) * ee13 + ee18/ee30)/ee19;
out(j, 12) = 1 - (((ee41 + xi * (ee39 - (ee15 * ee18/ee19 +
   ee44 + mu))) * ee13 - ((ee41 + xi * (ee39 - (ee44 + mu)))/ee21 +
   ee14 + ee8)/xi)/ee11 + (ee31 + (ee43 - ee25/(xi * ee20))/xi) * ee15)/ee9;
out(j, 13) = -(xi * ((ee13 * (xi * (ee32 - 1) - ee29) + ee24)/ee11) * ee15/ee9);

// ee1 = exp(lpsi);
// ee2 = R_pow(theta, xi);
// ee3 = log(theta);
// ee4 = 1 - ee2;
// ee5 = ee4 * ee1;
// ee8 = ee5/xi + y - mu;
// ee9 = xi * ee8;
// ee10 = xi * ee3;
// ee11 = exp(lpsi + ee10);
// ee12 = ee9/ee11;
// ee13 = 1 + ee12;
// ee14 = 1/xi;
// ee15 = 1 + ee14;
// ee16 = xi - 1;
// ee17 = R_pow(theta, ee16);
// ee18 = ee2 * ee1;
// ee19 = ee3 * (ee18 + ee9);
// ee20 = ee13 * ee11;
// ee22 = y - (ee19 + mu);
// ee23 = R_pow(ee13, ee15);
// ee24 = ee17 * ee1;
// ee25 = ee9/theta;
// ee26 = ee5 - ee9;
// ee27 = R_pow(ee13, ee14);
// ee28 = ee24 + ee25;
// ee29 = log1p(ee12);
// ee30 = ee23 * ee11;
// ee31 = R_pow(theta, xi - 2);
// ee33 = R_pow(ee13, ee14 + 2) * ee11;
// ee34 = 1/ee27;
// ee35 = ee15 * ee22;
// ee36 = R_pow(xi, 2);
// ee37 = ee35/ee33;
// ee38 = ee22/ee30;
// ee39 = 1/theta;
// ee40 = R_pow(theta, 2);
// ee41 = xi * ee15;
// ee42 = (ee8 * (1 - ee10) + y - (mu + ee18 * ee3))/theta;
// ee43 = ee26 * ee28;
// ee45 = ee26/ee20 + 1;
// ee48 = ee8 * (1 + xi)/ee40 + ee31 * ee1;
// ee49 = (ee4/theta + ee17 - ee17) * ee1;
// ee50 = ((1 + ee2 - 2 * ee2) * ee1 - ee9) * ee3;
// ee51 = (2 * ee17 + ee17 * xi * ee3 - ee17) * ee1;
// ee52 = (2 * y - (2 * mu + ee19)) * ee3;
// ee53 = ee28/ee30;
// ee54 = R_pow(ee28, 2);
// ee55 = ee22/ee20;
// ee56 = 1/ee23;
// ee57 = ee1 * (ee31 * ee16 - ee31 * xi);
// ee58 = ee29/(xi * ee27);
// ee60 = ee29/(ee36 * ee23) - ee37;
// ee61 = ee24 * ee3;
// 
// out(j, 0) = -((ee41 - ee34)/ee20);
// out(j, 1) = ((1 - ee34)/xi + 1) * ee26/ee20 + 1;
// out(j, 2) = ((ee34 - 1) * ee29/xi - ee38)/xi + ee35/ee20 + ee3;
// out(j, 3) = ee53 + xi * (ee39 - ee15 * ee28/ee20);
// out(j, 4) =  - (ee41 * (xi - ee34)/(R_pow(ee13, 2) * R_pow(ee11, 2)));
// out(j, 5) = -(((ee26/ee30 - xi * ee45) * ee15 + ee34)/ee13/ee11);
// out(j, 6) = -(((ee38 + 1 - xi * (ee55 + ee3)) * ee15 + ee3/ee27 -
//    (1 + ee58)/xi)/ee13/ee11);
// out(j, 7) = xi * (((ee53 + xi * (ee39 - ee28/ee20)) * ee15 -
//    1/(theta * ee27))/ee13)/ee11;
// out(j, 8) = ((ee26 * ee15/ee33 + ee56)/xi - ee45 * ee15/ee13) * ee26/ee11;
// out(j, 9) = (((ee26/xi + ee50 + y - mu)/ee23 - ee26 * ee60)/xi -
//    ((ee26 * ee22/ee20 + ee50 + y - mu) * ee15 + ee26/ee36)/ee13)/ee11;
// out(j, 10) = -(((ee43/ee30 + xi * (ee49 - (ee43/ee20 + ee25))) * ee15 -
//    (ee49 - ee25)/ee27)/ee13/ee11);
// out(j, 11) = ((((ee55 - 2 * (ee29/xi))/R_pow(ee13, ee14 - 1) -
//    ee22/ee11)/ee13 + (2 + ee58 - ee38) * ee29/xi)/xi + ((ee52 +
//    ee22/xi)/ee23 - ee60 * ee22)/ee11)/xi - ((ee52 + R_pow(ee22, 2)/ee20) * ee15 +
//    ee22/ee36)/ee20;
// out(j, 12) = ee39 - (((ee51 + xi * (ee42 - (ee28 * ee22/ee20 +
//    ee61))) * ee15 - ((ee51 + xi * (ee42 - ee61))/ee27 + ee24 +
//    ee25)/xi)/ee13 + (ee37 + (ee56 - ee29/(xi * ee23))/xi) * ee28)/ee11;
// out(j, 13) = -(xi * (ee15 * (ee57 + xi * (ee54/ee20 - ee48))/ee20 +
//    1/ee40) - ((ee57 - xi * ee48)/ee23 + ee41 * ee54/ee33)/ee11);

} else {
/*
ee1 = log(theta);
ee2 = R_pow(theta, xi);
ee3 = exp(lpsi);
ee4 = 1 - ee2;
ee5 = xi * ee1;
ee6 = exp(lpsi + ee5);
ee8 = ee4 * ee3/xi;
ee10 = ee8 + y - mu;
ee11 = xi - 1;
ee12 = R_pow(theta, ee11);
ee14 = exp(-(ee10/ee6));
ee15 = ee4/xi;
ee16 = ee15 + ee2 * ee1;
ee18 = ee10 * ee1 + ee16 * ee3/xi;
ee20 = xi * ee10/theta;
ee21 = R_pow(theta, xi - 2);
ee22 = mu - y;
ee24 = ee12 * ee3 + ee20;
ee26 = ee18/ee6 - ee1;
ee27 = ee14 - 1;
ee28 = R_pow(theta, 2);
ee29 = xi/theta;
ee31 = ee26 * ee14 + ee1;
ee33 = ee10 * (1 - ee5)/theta;
ee37 = (ee4/theta + ee12 - ee12) * ee3;
ee41 = ((1 + ee2 + ee2 - ee2) * ee1 + 2 * ee15) * ee3/R_pow(xi, 2);
ee42 = ((ee12 + ee12 * xi * ee1 - ee12)/xi - (ee16/theta +  ee12 * ee1)) * ee3;
ee43 = 2 * ee10;
ee44 = 2 * ee8;
ee45 = ee3 * (ee21 * ee11 - ee21 * xi);
ee46 = xi * (ee10 * (1 + xi)/ee28 + ee21 * ee3);

out(j, 0) = ee27/ee6;
out(j, 1) = (1 - ee14) * ee22/ee6 + 1;
out(j, 2) = ee18 * ee27/ee6 + ee1;
out(j, 3) = ee27 * ee24/ee6 + ee29;
out(j, 4) = ee14/R_pow(ee6, 2);
out(j, 5) = -(((ee22/ee6 + 1) * ee14 - 1)/ee6);
out(j, 6) = ee31/ee6;
out(j, 7) = ((ee24/ee6 - ee29) * ee14 + ee29)/ee6;
out(j, 8) = ((R_pow(ee22, 2)/ee6 + ee44 + y - (ee43 + mu)) * ee14 +
   ee43 + mu - (ee44 + y))/ee6;
out(j, 9) = -(ee31 * ee22/ee6);
out(j, 10) = -((ee37 + (ee22 * ee24/ee6 + ee20 - ee37) * ee14 -
   ee20)/ee6);
out(j, 11) = ((ee26 * ee18 - ee41) * ee14 + ee18 * ee1 + ee41)/ee6;
out(j, 12) = ((ee18 * ee24/ee6 + ee33 + ee42) * ee14 - (ee33 +
   ee42))/ee6 + 1/theta;
out(j, 13) = ((R_pow(ee24, 2)/ee6 + ee45 - ee46) * ee14 + ee46 -
   ee45)/ee6 - xi/ee28;
   */
ee1 = exp(lpsi);
ee2 = theta * ee1;
ee4 = y - (mu + ee2);
ee5 = ee4/ee1;
ee7 = exp(-ee5);
ee8 = ee5 + theta;
ee9 = ee8 * ee7;
ee10 = (y - (2 * ee2 + 2 * ee4 + mu))/ee1;
ee11 = ee7 - 1;

out(j, 0) = ee11/ee1;
out(j, 1) = ee9 + 1 - ee8;
out(j, 2) = 0;
out(j, 3) = ee11;
out(j, 4) = ee7/R_pow(ee1, 2);
out(j, 5) = ((ee8 - 1) * ee7 + 1)/ee1;
out(j, 6) = 0;
out(j, 7) = ee7/ee1;
out(j, 8) = (R_pow(ee8, 2) + ee10) * ee7 - ee10;
out(j, 9) = 0;
out(j, 10) = ee9;
out(j, 11) = 0;
out(j, 12) = 0;
out(j, 13) = ee7;
    
}

}

return(out);

}

// // [[Rcpp::export]]
// arma::mat ldgevagg12v2(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec, arma::vec thetavec)
// {
//     
// int nobs = yvec.size();
// arma::mat out(nobs, 14);
// 
// double y, mu, lpsi, xi, theta;
// double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
// double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee18;
// double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
// double ee31, ee33, ee37;
// double ee41, ee42, ee43, ee44, ee45, ee46, ee48;
// double ee50, ee53;
// double ee64, ee66, ee68;
// double ee70, ee72, ee73, ee74, ee76, ee78;
// double ee81, ee82, ee83, ee85;
// double ee91, ee105, ee106, ee107, ee109;
// double ee111, ee113, ee114, ee115, ee116, ee117, ee119;
// double ee121, ee124, ee135, ee138, ee141, ee149, ee152, ee172, ee176, ee198, ee202, ee213, ee217;
// double ee220, ee222, ee224, ee231, ee233, ee235, ee243, ee252, ee258, ee260, ee280, ee307, ee313;
// double ee320, ee321, ee335, ee337, ee349, ee358;
// 
// 
// for (int j=0; j < nobs; j++) {
// 
// y = yvec[j];
// mu = muvec[j];
// lpsi = lpsivec[j];
// xi = xivec[j];
// theta = thetavec[j];
// 
// if (fabs(xi) >= xieps) {
// 
// ee2 = 1/xi + 1;
// ee3 = log(theta);
// ee6 = exp(lpsi + xi * ee3);
// ee7 = R_pow(theta, xi);
// ee8 = 1 - ee7;
// ee9 = ee6 * ee8;
// ee10 = ee9/xi;
// ee12 = y - (mu - ee10);
// ee13 = xi * ee12;
// ee14 = ee13/ee6;
// ee15 = log1p(ee14);
// ee18 = 1 + ee14;
// ee20 = -1/xi;
// ee21 = R_pow(ee18, ee20);
// ee23 = ee20 - 1;
// ee24 = R_pow(ee18, ee23);
// ee25 = xi/ee6;
// ee26 = ee20 * ee25;
// ee28 = ee25/ee18;
// ee33 = ee18 * ee18;
// ee37 = R_pow(ee18, ee23 - 1);
// ee43 = xi * ee10;
// ee44 = ee43/ee6;
// ee45 = ee13 * ee6;
// ee46 = ee6 * ee6;
// ee48 = ee44 - ee45/ee46;
// ee50 = ee37 * (ee23 * ee48);
// ee53 = xi * ee6/ee46;
// ee64 = ee6 * ee3;
// ee66 = ee7 * ee3;
// ee68 = ee64 * ee8 - ee6 * ee66;
// ee70 = xi * xi;
// ee72 = ee68/xi - ee9/ee70;
// ee73 = xi * ee72;
// ee74 = ee12 + ee73;
// ee76 = ee13 * ee64;
// ee78 = ee74/ee6 - ee76/ee46;
// ee81 = log(ee18);
// ee82 = 1/ee70;
// ee83 = ee81 * ee82;
// ee85 = ee37 * (ee23 * ee78) + ee24 * ee83;
// ee91 = 1/ee6 - xi * ee64/ee46;
// ee105 = 1/theta;
// ee106 = xi * ee105;
// ee107 = ee6 * ee106;
// ee109 = xi - 1;
// ee111 = R_pow(theta, ee109) * xi;
// ee113 = ee107 * ee8 - ee6 * ee111;
// ee114 = ee113/xi;
// ee115 = xi * ee114;
// ee116 = ee115/ee6;
// ee117 = ee13 * ee107;
// ee119 = ee116 - ee117/ee46;
// ee121 = ee37 * (ee23 * ee119);
// ee124 = xi * ee107/ee46;
// ee135 = ee48/ee18;
// ee138 = ee20 * ee48;
// ee141 = ee43 * ee6;
// ee149 = ee46 * ee46;
// ee152 = ee44 - ee141/ee46 - ((ee141 + ee45)/ee46 - 
// ee45 * (2 * (ee6 * ee6))/ee149);
// ee172 = 2 * (ee64 * ee6);
// ee176 = (ee10 + ee73)/ee6 - ee43 * ee64/ee46 - 
// ((ee74 * ee6 + ee76)/ee46 - ee45 * ee172/ee149);
// ee198 = 2 * (ee107 * ee6);
// ee202 = ee116 - ee43 * ee107/ee46 - ((ee115 * 
// ee6 + ee117)/ee46 - ee45 * ee198/ee149);
// ee213 = ee78/ee18;
// ee217 = ee20 * ee78;
// ee220 = ee24 * ee217 + ee21 * ee83;
// ee222 = ee64 * ee3;
// ee224 = ee64 * ee66;
// ee231 = ee68/ee70;
// ee233 = 2 * xi;
// ee235 = ee70 * ee70;
// ee243 = ee74 * ee64;
// ee252 = (ee72 + (ee72 + xi * ((ee222 * ee8 - 
// ee224 - (ee224 + ee6 * (ee66 * ee3)))/xi - 
// ee231 - (ee231 - ee9 * ee233/ee235))))/ee6 - 
// ee243/ee46 - ((ee243 + ee13 * ee222)/ee46 - 
// ee76 * ee172/ee149);
// ee258 = ee82 * ee213;
// ee260 = ee233/ee235;
// ee280 = ee107 * ee3 + ee6 * ee105;
// ee307 = (ee114 + xi * ((ee280 * ee8 - ee64 * 
// ee111 - (ee107 * ee66 + ee6 * (ee111 * 
// ee3 + ee7 * ee105)))/xi - ee113/ee70))/ee6 - 
// ee74 * ee107/ee46 - ((ee115 * ee64 + ee13 * 
// ee280)/ee46 - ee76 * ee198/ee149);
// ee313 = ee119/ee18;
// ee320 = ee20 * ee119;
// ee321 = ee24 * ee320;
// ee335 = ee107 * ee106 - ee6 * (xi * (1/R_pow(theta, 2)));
// ee337 = ee107 * ee111;
// ee349 = ee115 * ee107;
// ee358 = xi * ((ee335 * ee8 - ee337 - (ee337 + 
// ee6 * (R_pow(theta, ee109 - 1) * ee109 * xi)))/xi)/ee6 - 
// ee349/ee46 - ((ee349 + ee13 * ee335)/ee46 - 
// ee117 * ee198/ee149);
// 
// out(j, 0) = -(ee24 * ee26 + ee2 * ee28);
// out(j, 1) = 1 + ee2 * ee135 + ee24 * ee138;
// out(j, 2) = ee2 * ee213 - ee82 * ee15 + ee220;
// out(j, 3) = ee2 * ee313 + ee321;
// 
// out(j, 4) =  -(ee2 * (ee25 * ee25/ee33) - ee37 * (ee23 * ee25) * ee26);
// out(j, 5) = -(ee50 * 
// ee26 - ee24 * (ee20 * ee53) - ee2 * (ee53/ee18 + 
// ee25 * ee48/ee33));
// out(j, 6) = -(ee85 * 
// ee26 + ee24 * (ee82 * ee25 + ee20 * ee91) + 
// (ee2 * (ee91/ee18 - ee25 * ee78/ee33) - 
// ee82 * ee28));
// out(j, 7) = -(ee121 * 
// ee26 - ee24 * (ee20 * ee124) - ee2 * (ee124/ee18 + 
// ee25 * ee119/ee33));
// 
// out(j, 8) = ee2 * (ee152/ee18 - 
// ee48 * ee48/ee33) + (ee50 * ee138 + ee24 * 
// (ee20 * ee152));
// out(j, 9) = ee2 * 
// (ee176/ee18 - ee48 * ee78/ee33) - ee82 * 
// ee135 + (ee85 * ee138 + ee24 * (ee82 * 
// ee48 + ee20 * ee176));
// out(j, 10) = ee2 * 
// (ee202/ee18 - ee48 * ee119/ee33) + (ee121 * 
// ee138 + ee24 * (ee20 * ee202));
// 
// out(j, 11) = ee2 * (ee252/ee18 - ee78 * 
// ee78/ee33) - ee258 - (ee258 - ee260 * 
// ee15) + (ee85 * ee217 + ee24 * (ee82 * 
// ee78 + ee20 * ee252) + (ee220 * ee83 + 
// ee21 * (ee213 * ee82 - ee81 * ee260)));
// out(j, 12) = ee2 * 
// (ee307/ee18 - ee78 * ee119/ee33) - ee82 * 
// ee313 + (ee121 * ee217 + ee24 * (ee20 * 
// ee307) + (ee321 * ee83 + ee21 * (ee313 * 
// ee82)));
// out(j, 13) = ee2 * (ee358/ee18 - 
// ee119 * ee119/ee33) + (ee121 * ee320 + 
// ee24 * (ee20 * ee358));
// 
// } else {
// 
// ee1 = log(theta);
// ee2 = R_pow(theta, xi);
// ee3 = exp(lpsi);
// ee4 = 1 - ee2;
// ee5 = xi * ee1;
// ee6 = exp(lpsi + ee5);
// ee8 = ee4 * ee3/xi;
// ee10 = ee8 + y - mu;
// ee11 = xi - 1;
// ee12 = R_pow(theta, ee11);
// ee14 = exp(-(ee10/ee6));
// ee15 = ee4/xi;
// ee16 = ee15 + ee2 * ee1;
// ee18 = ee10 * ee1 + ee16 * ee3/xi;
// ee20 = xi * ee10/theta;
// ee21 = R_pow(theta, xi - 2);
// ee22 = mu - y;
// ee24 = ee12 * ee3 + ee20;
// ee26 = ee18/ee6 - ee1;
// ee27 = ee14 - 1;
// ee28 = R_pow(theta, 2);
// ee29 = xi/theta;
// ee31 = ee26 * ee14 + ee1;
// ee33 = ee10 * (1 - ee5)/theta;
// ee37 = (ee4/theta + ee12 - ee12) * ee3;
// ee41 = ((1 + ee2 + ee2 - ee2) * ee1 + 2 * ee15) * ee3/R_pow(xi, 2);
// ee42 = ((ee12 + ee12 * xi * ee1 - ee12)/xi - (ee16/theta +  ee12 * ee1)) * ee3;
// ee43 = 2 * ee10;
// ee44 = 2 * ee8;
// ee45 = ee3 * (ee21 * ee11 - ee21 * xi);
// ee46 = xi * (ee10 * (1 + xi)/ee28 + ee21 * ee3);
// 
// out(j, 0) = ee27/ee6;
// out(j, 1) = (1 - ee14) * ee22/ee6 + 1;
// out(j, 2) = ee18 * ee27/ee6 + ee1;
// out(j, 3) = ee27 * ee24/ee6 + ee29;
// out(j, 4) = ee14/R_pow(ee6, 2);
// out(j, 5) = -(((ee22/ee6 + 1) * ee14 - 1)/ee6);
// out(j, 6) = ee31/ee6;
// out(j, 7) = ((ee24/ee6 - ee29) * ee14 + ee29)/ee6;
// out(j, 8) = ((R_pow(ee22, 2)/ee6 + ee44 + y - (ee43 + mu)) * ee14 +
//    ee43 + mu - (ee44 + y))/ee6;
// out(j, 9) = -(ee31 * ee22/ee6);
// out(j, 10) = -((ee37 + (ee22 * ee24/ee6 + ee20 - ee37) * ee14 -
//    ee20)/ee6);
// out(j, 11) = ((ee26 * ee18 - ee41) * ee14 + ee18 * ee1 + ee41)/ee6;
// out(j, 12) = ((ee18 * ee24/ee6 + ee33 + ee42) * ee14 - (ee33 +
//    ee42))/ee6 + 1/theta;
// out(j, 13) = ((R_pow(ee24, 2)/ee6 + ee45 - ee46) * ee14 + ee46 -
//    ee45)/ee6 - xi/ee28;
//     
// }
// 
// }
// 
// return(out);
// 
// }

// // [[Rcpp::export]]
// arma::vec lpsi(arma::vec yvec, arma::vec muvec, arma::vec psivec, arma::vec xivec, arma::vec thetavec)
// {
//     
// int nobs = yvec.size();
// arma::vec nllh(nobs);
// 
// double y, mu, psi, xi, theta;
// double ee1, ee2;
// 
// for (int j=0; j < nobs; j++) {
// 
// y = yvec[j];
// mu = muvec[j];
// psi = psivec[j];
// xi = xivec[j];
// theta = thetavec[j];
// 
// if (fabs(xi) >= xieps) {
// 
// mu = mu - psi * (1 - R_pow(theta, xi)) / xi;
// psi = psi * R_pow(theta, xi);
// 
// ee1 = xi * (y - mu) / psi;
// 
// if (ee1 <= -1.0) {
//     nllh[j] = 1e20;
//     break;
// } else {
// 
// ee2 = 1.0 / xi;
// 
// nllh[j] = log(psi) + (ee2 + 1.0) * log1p(ee1) + R_pow(1.0 + ee1, -ee2);
// }
// 
// } else {
// 
// mu = mu + psi * log(theta);
// 
// ee1 = (y - mu) / psi;
// nllh[j] = log(psi) + ee1 + exp(-ee1);
//     
// }
// 
// }
// 
// return(nllh);
// 
// }
// 
// // [[Rcpp::export]]
// arma::mat lpsi12(arma::vec yvec, arma::vec muvec, arma::vec psivec, arma::vec xivec, arma::vec thetavec)
// {
//     
// int nobs = yvec.size();
// arma::mat out(nobs, 14);
// 
// double y, mu, psi, xi, theta;
// double ee1, ee2, ee4, ee5, ee6, ee7, ee8, ee9;
// double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
// double ee22, ee23, ee24, ee25, ee26, ee29;
// double ee30, ee31, ee32, ee33, ee34, ee35, ee37, ee39;
// double ee40, ee41, ee42, ee43, ee44, ee46, ee47, ee48, ee49;
// double ee50, ee51, ee52, ee55;
// double ee63, ee67, ee68, ee69;
// double ee70, ee74;
// 
// for (int j=0; j < nobs; j++) {
// 
// y = yvec[j];
// mu = muvec[j];
// psi = psivec[j];
// xi = xivec[j];
// theta = thetavec[j];
// 
// if (fabs(xi) >= xieps) {
// 
// ee1 = R_pow(theta, xi);
// ee2 = psi * ee1;
// ee4 = mu - psi * (1 - ee1);
// ee6 = y - ee4/xi;
// ee7 = log(theta);
// ee9 = xi * ee6/ee2;
// ee10 = R_pow(ee2, 2);
// ee11 = 1 + ee9;
// ee12 = 1/xi;
// ee13 = xi - 1;
// ee14 = R_pow(theta, ee13);
// ee15 = 1 + ee12;
// ee16 = ee2 * ee7;
// ee17 = y - ee16;
// ee18 = ee2 * xi;
// ee19 = ee17/ee2;
// ee22 = ee18 * ee7 * ee6/ee10;
// ee23 = ee19 - ee22;
// ee24 = R_pow(ee11, ee15);
// ee25 = psi * ee14;
// ee26 = 1/theta;
// ee29 = ee25 * xi * ee6/ee10;
// ee30 = ee26 + ee29;
// ee31 = log1p(ee9);
// ee32 = R_pow(ee11, ee12);
// ee33 = 2 * xi;
// ee34 = 3 * xi;
// ee35 = R_pow(psi, 2);
// ee37 = R_pow(theta, 2 * ee13);
// ee39 = ee14 * xi * ee7;
// ee40 = ee23 * ee15;
// ee41 = R_pow(ee11, ee12 + 2);
// ee42 = R_pow(xi, 2);
// ee43 = ee40/ee41;
// ee44 = ee13/R_pow(theta, 2);
// ee46 = R_pow(theta, ee33 - 1);
// ee47 = R_pow(theta, ee33);
// ee48 = ee14 + ee39;
// ee49 = (((2 - ee12) * ee4 + ee16)/xi + psi * ee7 * (ee1 *  ee7 - ee1/xi))/ee2;
// ee50 = ee23/ee24;
// ee51 = R_pow(ee30, 2);
// ee52 = (2 * ee14 + ee39 - ee14)/ee1;
// ee55 = ee48 * ee6 + ee14 * ee17 - psi * (2 * (psi * R_pow(theta, ee34 - 1) * xi * ee6/ee10) + ee46) * ee7;
// ee63 = ee31/(ee42 * ee24) - ee43;
// ee67 = psi * (2 * (ee1 * ee17) + xi * ee7 * (ee1 - 2 * (ee35 * R_pow(theta, ee34)/ee10)) * ee6) * ee7/ee10;
// ee68 = ee2 * ee11;
// ee69 = ee16/ee10;
// ee70 = psi * xi;
// ee74 = R_pow(theta, xi - 2) * ee13 * ee6 - psi * (2 * (psi * R_pow(theta, ee34 - 2) * xi * ee6/ee10) + ee37 + ee37);
// 
// out(j, 0) = -((ee15/ee1 - 1/(ee1 * xi * ee32))/(psi * ee11));
// out(j, 1) = 0;
// out(j, 2) = ((1/ee32 - 1) * ee31/xi - ee50)/xi + ee40/ee11 +
//    ee7;
// out(j, 3) = ee30/ee24 + xi * (ee26 - ee15 * ee30/ee11);
// out(j, 4) =  - (ee15 * (1/ee47 - 1/(ee47 * xi * ee32))/(ee35 * R_pow(ee11, 2)));
// out(j, 5) = 0;
// out(j, 6) = -(((1/ee18 + ee69)/ee24 - ee63/ee2)/xi - ((ee23/ee68 +
//    ee69) * ee15 + 1/(ee2 * ee42))/ee11);
// out(j, 7) = ((ee30/(ee2 * ee24) + xi * (ee25/ee10 - ee30/ee68)) * ee15 -
//    ee25/(ee32 * ee10))/ee11;
// out(j, 8) = 0;
// out(j, 9) = 0;
// out(j, 10) = 0;
// out(j, 11) = ((((ee23/ee11 - 2 * (ee31/xi))/R_pow(ee11, ee12 -
//    1) + ee22 - ee19)/ee11 + (2 + ee31/(xi * ee32) - ee50) * ee31/xi)/xi +
//    (ee49 + ee23/xi + ee67)/ee24 - ee23 * ee63)/xi +
//    (1 - ee35 * ee47/ee10) * R_pow(ee7, 2) - ((ee49 + R_pow(ee23, 2)/ee11 +
//    ee67) * ee15 + ee23/ee42)/ee11;
// out(j, 12) = (((ee52 + ee70 * ee55/ee10)/ee32 + ee26 + ee29)/xi -
//    (ee52 + xi * (psi * ee55/ee10 - ee23 * ee30/ee11)) * ee15)/ee11 +
//    ee48/ee1 - ((ee43 + (1/ee24 - ee31/(xi * ee24))/xi) * ee30 +
//    ee35 * ee46 * xi * ee7/ee10);
// out(j, 13) = (ee44 + ee70 * ee74/ee10)/ee24 + xi * ee15 * ee51/ee41 +
//    xi * (ee44 - ((ee44 + xi * (ee51/ee11 + psi * ee74/ee10)) * ee15/ee11 +
//    ee35 * ee37 * xi/ee10));
// 
// } else {
// 
// ee1 = log(theta);
// ee2 = (y - (mu + psi * ee1))/psi;
// ee4 = exp(-ee2);
// ee5 = ee2 + ee1;
// ee6 = R_pow(psi, 2);
// ee7 = ee5 * ee4;
// ee8 = 2 * ee2;
// ee9 = 2 * ee1;
// ee10 = ee4 - 1;
// ee11 = psi * theta;
// 
// out(j, 0) = ee10/psi;
// out(j, 1) = (ee7 + 1 - ee5)/psi;
// out(j, 2) = 0;
// out(j, 3) = ee10/theta;
// out(j, 4) = ee4/ee6;
// out(j, 5) = ((ee5 - 1) * ee4 + 1)/ee6;
// out(j, 6) = 0;
// out(j, 7) = ee4/ee11;
// out(j, 8) = ((R_pow(ee5, 2) - (ee8 + ee9)) * ee4 + ee8 + ee9 -
//    1)/ee6;
// out(j, 9) = 0;
// out(j, 10) = ee7/ee11;
// out(j, 11) = 0;
// out(j, 12) = 0;
// out(j, 13) = 1/R_pow(theta, 2);
// 
// }
// 
// }
// 
// return(out);
// 
// }
