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

// shape = 1.5 / (1 + exp(-xi)) - .5

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

// shape = 1.5 / (1 + exp(-xi)) - .5

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
out(j, 3) =  0.0;
out(j, 4) =  - (ee18 * ee23/(R_pow(ee5, 2) * R_pow(ee1, 2)));
out(j, 5) = (xi * ee16 * ee7 - ee10)/ee5/ee1;
out(j, 6) = -(ee19/ee1);
out(j, 7) = 0.0;
out(j, 8) = -((ee10 + xi * (ee23 * ee2/ee9 - 1) * ee7)/ee5 * ee2/ee1);
out(j, 9) = -(ee19 * ee2/ee1);
out(j, 10) = 0.0;
out(j, 11) = ((((ee2/ee9 - 2 * (ee11/xi))/R_pow(ee5, ee6 - 1) -
   ee2/ee1)/ee5 + (2 + ee17 - ee22) * ee11/xi)/xi + (ee13/(R_pow(ee5, ee6 +
   2) * ee1) + (1/ee12 - ee11/(xi * ee12))/xi) * ee2/ee1)/xi -
   (ee20 + 1/R_pow(xi, 2)) * ee2/ee9;
out(j, 12) = 0.0;
out(j, 13) = 0.0;

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

// shape = 1.5 / (1 + exp(-xi)) - .5
// dep = 1 / (1 + exp(-theta))

// [[Rcpp::export]]
double ldgevagg_logit(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec, arma::vec thetavec)
{
    
int nobs = yvec.size();

double y, mu, lpsi, psi, xi, lgttheta, theta, theta2;
double ee1, ee2;
double nllh = 0.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
psi = exp(lpsi);
xi = xivec[j];
lgttheta = thetavec[j];
theta = 1 / (1 + exp(-lgttheta));

if (fabs(xi) >= xieps) {

theta2 = R_pow(theta, xi);
mu = mu - exp(lpsi) * (1 - theta2) / xi;
lpsi = lpsi + xi * log(theta);

ee1 = xi * (y - mu) / exp(lpsi);

if (ee1 <= -1.0) {
    nllh = 1e20;
    break;
} else {

ee2 = 1.0 / xi;
nllh += lpsi + (ee2 + 1.0) * log1p(ee1) + R_pow(1.0 + ee1, -ee2);
    
}

} else {

mu = mu + psi * log(theta2);

ee1 = (y - mu) / exp(lpsi);
nllh += lpsi + ee1 + exp(-ee1);
    
}

}

return nllh;

}

// shape = 1.5 / (1 + exp(-xi)) - .5
// dep = 1 / (1 + exp(-theta))

// [[Rcpp::export]]
arma::mat ldgevagg12_logit(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec, arma::vec thetavec)
{
    
int nobs = yvec.size();
arma::mat out(nobs, 14);

double y, mu, lpsi, xi, lgttheta, theta;

// double ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
// double ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
// double ee20, ee23, ee24, ee25, ee26, ee28, ee29;
// double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee37, ee38, ee39;
// double ee40, ee41, ee42, ee44, ee45, ee47, ee48, ee49;
// double ee51, ee52, ee53, ee54, ee55, ee57, ee59;
// double ee62, ee63, ee64, ee65, ee66, ee68, ee69;

double ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee20, ee21, ee23, ee24, ee25, ee27, ee28, ee29;
double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee37, ee38, ee39;
double ee40, ee41, ee42, ee44, ee45, ee46, ee48, ee49;
double ee51, ee53, ee54, ee57, ee58;
double ee60, ee61, ee62, ee63, ee65, ee66;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];
lgttheta = thetavec[j];
theta = 1.0 / (1.0 + exp(-lgttheta));

if (fabs(xi) >= xieps) {

// ee2 = exp(-lgttheta);
// ee3 = 1 + ee2;
// ee4 = 1/ee3;
// ee5 = exp(lpsi);
// ee6 = R_pow(ee4, xi);
// ee7 = log1p(ee2);
// ee8 = (1 - ee6) * ee5;
// ee11 = ee8/xi + y - mu;
// ee12 = xi * ee11;
// ee13 = xi * ee7;
// ee14 = exp(lpsi + ee13);
// ee15 = ee12/ee14;
// ee16 = 1 + ee15;
// ee17 = 1/xi;
// ee18 = 1 + ee17;
// ee19 = xi - 1;
// ee20 = R_pow(ee4, ee19);
// ee23 = (ee6 * ee5 - ee12) * ee7 + y - mu;
// ee24 = ee16 * ee14;
// ee25 = R_pow(ee16, ee18);
// ee26 = ee8 - ee12;
// ee28 = ee20 * ee5/ee3;
// ee29 = R_pow(ee16, ee17);
// ee30 = ee28 - ee12;
// ee31 = log1p(ee15);
// ee32 = ee25 * ee14;
// ee33 = 1/ee29;
// ee34 = R_pow(ee16, ee17 + 2);
// ee35 = ee23 * ee18;
// ee36 = R_pow(ee3, 2);
// ee37 = ee34 * ee14;
// ee38 = 2 * ee6;
// ee39 = xi * ee20;
// ee40 = R_pow(xi, 2);
// ee41 = ee35/ee37;
// ee42 = ee23/ee32;
// ee44 = 1/(ee36 * R_pow(ee4, 2));
// ee45 = ((ee20 - (ee20 + ee39 * ee7))/ee3 + ee6) * ee5;
// ee47 = ((R_pow(ee4, xi - 2) * ee19/ee3 + ee39) * ee2/ee3 - (1 - 2 * (ee2/ee3)) * ee20) * ee5/ee3;
// ee48 = ee23/ee24;
// ee49 = ee26 * ee30;
// ee51 = ee26/ee24 + 1;
// ee52 = ee11 * ((ee44 - (2 + xi)) * ee2/ee3 + 1);
// ee53 = ee11 * (1 - ee13);
// ee54 = ee30/ee32;
// ee55 = R_pow(ee30, 2);
// ee57 = (ee20/ee3 + ee6) * ee5 * ee7;
// ee59 = (ee6 - 1) * ee5 + ee12;
// ee62 = ((ee38 - (ee6 + 1)) * ee5 + ee12) * ee7 + mu;
// ee63 = ee3 * ee14;
// ee64 = (2 * mu + ee7 * (ee12 - (ee6 + ee38) * ee5) - 2 *  y) * ee7;
// ee65 = 1/ee25;
// ee66 = ee31/(xi * ee29);
// ee68 = ee31/(ee40 * ee25) - ee41;
// ee69 = xi * ee18;
// 
// out(j, 0) = -((ee69 - ee33)/ee24);
// out(j, 1) = ee26 * ((1 - ee33)/xi + 1)/ee24 + 1;
// out(j, 2) = ee35/ee24 + ((ee33 - 1) * ee31/xi - ee42)/xi + ee7;
// out(j, 3) = -(ee2 * (xi * (ee30 * ee18/ee24 + 1) - ee54)/ee3);
// out(j, 4) =  - (ee69 * (xi - ee33)/(R_pow(ee16, 2) * R_pow(ee14, 2)));
// out(j, 5) = -(((ee26/ee32 - xi * ee51) * ee18 + ee33)/ee16/ee14);
// out(j, 6) = -(((ee42 + 1 - xi * (ee48 + ee7)) * ee18 + ee7/ee29 -
//    (1 + ee66)/xi)/ee16/ee14);
// out(j, 7) = -(xi * ((ee18 * (xi * (ee30/ee24 + 1) - ee54) -
//    ee33)/ee16) * ee2/ee63);
// out(j, 8) = ((ee26 * ee18/ee37 + ee65)/xi - ee51 * ee18/ee16) * ee26/ee14;
// out(j, 9) = (((ee26/xi + y - ee62)/ee25 - ee26 * ee68)/xi +
//    ((ee62 - (ee23 * ee26/ee24 + y)) * ee18 - ee26/ee40)/ee16)/ee14;
// out(j, 10) = -(((ee49/ee32 + xi * (ee59 - ee49/ee24)) * ee18 -
//    ee59/ee29)/ee16 * ee2/ee63);
// out(j, 11) = ((((ee48 - 2 * (ee31/xi))/R_pow(ee16, ee17 - 1) -
//    ee23/ee14)/ee16 + (2 + ee66 - ee42) * ee31/xi)/xi + ((ee23/xi -
//    ee64)/ee25 - ee23 * ee68)/ee14)/xi + ((ee64 - R_pow(ee23, 2)/ee24) * ee18 -
//    ee23/ee40)/ee24;
// out(j, 12) = ((((ee28 - ((xi * (ee53 + ee57 + y - mu) - ee45)/ee29 +
//    ee12))/xi + ee18 * (xi * (ee23 * ee30/ee24 + ee53 +
//    ee57 + y - mu) - ee45))/ee16 - (ee41 + (ee65 - ee31/(xi * ee25))/xi) * ee30)/ee14 -
//    1) * ee2/ee3;
// out(j, 13) = -(ee2 * (xi * ((ee47 + xi * ((ee55/ee24 + ee28) * ee2/ee3 +
//    ee52)) * ee18/ee24 - ((ee44 - 2) * ee2/ee3 + 1)) -
//    ((ee47 + xi * (ee52 + ee20 * ee2 * ee5/ee36))/ee25 + xi * ee55 * ee18 * ee2/(ee3 * ee34 * ee14))/ee14)/ee3);

ee2 = exp(-lgttheta);
ee3 = 1 + ee2;
ee4 = 1/ee3;
ee5 = exp(lpsi);
ee6 = R_pow(ee4, xi);
ee7 = log1p(ee2);
ee8 = (1 - ee6) * ee5;
ee11 = ee8/xi + y - mu;
ee12 = xi * ee11;
ee13 = xi * ee7;
ee14 = exp(lpsi - ee13);
ee15 = ee12/ee14;
ee16 = 1 + ee15;
ee17 = 1/xi;
ee18 = 1 + ee17;
ee19 = xi - 1;
ee20 = R_pow(ee4, ee19);
ee21 = (ee6 * ee5 + ee12) * ee7;
ee23 = ee21 + y - mu;
ee24 = ee16 * ee14;
ee25 = R_pow(ee16, ee18);
ee27 = ee20 * ee5/ee3;
ee28 = ee8 - ee12;
ee29 = R_pow(ee16, ee17);
ee30 = ee27 + ee12;
ee31 = log1p(ee15);
ee32 = ee25 * ee14;
ee33 = 1/ee29;
ee34 = R_pow(ee16, ee17 + 2);
ee35 = ee23 * ee18;
ee36 = R_pow(ee3, 2);
ee37 = ee34 * ee14;
ee38 = xi * ee20;
ee39 = R_pow(xi, 2);
ee40 = ee35/ee37;
ee41 = ee23/ee32;
ee42 = ee30/ee32;
ee44 = 1/(ee36 * R_pow(ee4, 2));
ee45 = ((ee42 + xi * (1 - ee30/ee24)) * ee18 - ee33)/ee16;
ee46 = ((ee20 - (ee20 + ee38 * ee7))/ee3 + ee6) * ee5;
ee48 = ((R_pow(ee4, xi - 2) * ee19/ee3 - ee38) * ee2/ee3 - (1 - 2 * (ee2/ee3)) * ee20) * ee5/ee3;
ee49 = ee23/ee24;
ee51 = ((ee6 + 1 - 2 * ee6) * ee5 - ee12) * ee7 + mu;
ee53 = ee28/ee24 + 1;
ee54 = ee11 * ((2 - (ee44 + xi)) * ee2/ee3 - 1);
ee57 = ee11 * (1 + ee13) + (ee20/ee3 + ee6) * ee5 * ee7 +  y;
ee58 = R_pow(ee30, 2);
ee60 = ee3 * ee14;
ee61 = (2 * mu - (ee21 + 2 * y)) * ee7;
ee62 = 1/ee25;
ee63 = ee31/(xi * ee29);
ee65 = ee31/(ee39 * ee25) - ee40;
ee66 = xi * ee18;

out(j, 0) = -((ee66 - ee33)/ee24);
out(j, 1) = ee28 * ((1 - ee33)/xi + 1)/ee24 + 1;
out(j, 2) = ee35/ee24 + ((ee33 - 1) * ee31/xi - ee41)/xi - ee7;
out(j, 3) = (ee42 + xi * (1 - ee30 * ee18/ee24)) * ee2/ee3;
out(j, 4) =  - (ee66 * (xi - ee33)/(R_pow(ee16, 2) * R_pow(ee14, 2)));
out(j, 5) = -(((ee28/ee32 - xi * ee53) * ee18 + ee33)/ee16/ee14);
out(j, 6) = -(((ee41 + 1 + xi * (ee7 - ee49)) * ee18 - ((1 +
   ee63)/xi + ee7/ee29))/ee16/ee14);
out(j, 7) = xi * ee45 * ee2/ee60;
out(j, 8) = ((ee28 * ee18/ee37 + ee62)/xi - ee53 * ee18/ee16) * ee28/ee14;
out(j, 9) = (((ee51 - (ee23 * ee28/ee24 + y)) * ee18 - ee28/ee39)/ee16 +
   ((ee28/xi + y - ee51)/ee25 - ee28 * ee65)/xi)/ee14;
out(j, 10) = -(ee45 * ee28 * ee2/ee60);
out(j, 11) = ((((ee49 - 2 * (ee31/xi))/R_pow(ee16, ee17 - 1) -
   ee23/ee14)/ee16 + (2 + ee63 - ee41) * ee31/xi)/xi + ((ee23/xi +
   ee61)/ee25 - ee23 * ee65)/ee14)/xi - ((R_pow(ee23, 2)/ee24 +
   ee61) * ee18 + ee23/ee39)/ee24;
out(j, 12) = (1 - (((ee46 + xi * (ee57 - (ee23 * ee30/ee24 +
   mu))) * ee18 - ((ee46 + xi * (ee57 - mu))/ee29 + ee27 + ee12)/xi)/ee16 +
   (ee40 + (ee62 - ee31/(xi * ee25))/xi) * ee30)/ee14) * ee2/ee3;
out(j, 13) = -(ee2 * (xi * ((ee48 + xi * ((ee58/ee24 - ee27) * ee2/ee3 +
   ee54)) * ee18/ee24 + (ee44 - 2) * ee2/ee3 + 1) -
   ((ee48 + xi * (ee54 - ee20 * ee2 * ee5/ee36))/ee25 + xi * ee58 * ee18 * ee2/(ee3 * ee34 * ee14))/ee14)/ee3);
   
} else {

ee2 = exp(-theta);
ee3 = exp(lpsi);
ee4 = log1p(ee2);
ee5 = ee3 * ee4;
ee7 = y - (ee5 + mu);
ee8 = ee7/ee3;
ee9 = 1 + ee2;
ee11 = exp(-ee8);
ee12 = ee8 + ee4;
ee13 = ee12 * ee11;
ee16 = (y - (2 * ee5 + 2 * ee7 + mu))/ee3;
ee17 = 1/(R_pow(ee9, 2) * R_pow(1/ee9, 2));

out(j, 0) = (ee11 - 1)/ee3;
out(j, 1) = ee13 + 1 - ee12;
out(j, 2) = 0;
out(j, 3) = (1 - ee11) * ee2/ee9;
out(j, 4) = ee11/R_pow(ee3, 2);
out(j, 5) = -(((1 - ee12) * ee11 - 1)/ee3);
out(j, 6) = 0;
out(j, 7) = -(ee11 * ee2/(ee9 * ee3));
out(j, 8) = (R_pow(ee12, 2) + ee16) * ee11 - ee16;
out(j, 9) = 0;
out(j, 10) = -(ee13 * ee2/ee9);
out(j, 11) = 0;
out(j, 12) = 0;
out(j, 13) = -(((ee17 - 2) * ee2/ee9 + 1 - ((ee17 - 1) * ee2/ee9 +
   1) * ee11) * ee2/ee9);

}

}

return(out);

}

// // shape = 1.5 / (1 + exp(-xi)) - .5
// // dep = exp(theta)
// 
// double ldgevagg_log(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec, arma::vec thetavec)
// {
//     
// int nobs = yvec.size();
// 
// double y, mu, lpsi, psi, xi, xi2, theta, theta2;
// double ee1, ee2;
// double nllh = 0.0;
// 
// for (int j=0; j < nobs; j++) {
// 
// y = yvec[j];
// mu = muvec[j];
// lpsi = lpsivec[j];
// xi = xivec[j];
// xi2 = 1.5 / (1 + exp(-xi)) - 0.5;
// theta = thetavec[j];
// theta2 = exp(theta);
// psi = exp(lpsi);
// 
// if (fabs(xi2) < xieps) {
//   if (xi2 < 0) {
//     xi2 = -xieps;
//   } else {
//     xi2 = xieps;
//   }
// }
// 
// ee1 = xi2 * (y - mu) / psi;
// 
// if (ee1 <= -1.0) {
//     nllh = 1e20;
//     break;
// } else {
// 
// ee2 = 1.0 / xi2;
// nllh += lpsi - theta + (ee2 + 1.0) * log1p(ee1) + theta2 * R_pow(1.0 + ee1, -ee2);
//     
// }
// 
// }
// 
// return nllh;
// 
// }
// 
// // shape = 1.5 / (1 + exp(-xi)) - .5
// // dep = exp(theta)
// 
// arma::mat ldgevagg12_log(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec, arma::vec thetavec)
// {
//     
// int nobs = yvec.size();
// arma::mat out(nobs, 14);
// 
// double y, mu, lpsi, xi, xi2, theta;
// 
// double ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
// double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
// double ee21, ee22, ee23, ee24, ee25, ee26;
// double ee30, ee35, ee37, ee38;
// double ee43, ee44, ee45, ee46, ee47, ee49;
// 
// for (int j=0; j < nobs; j++) {
// 
// y = yvec[j];
// mu = muvec[j];
// lpsi = lpsivec[j];
// xi = xivec[j];
// xi2 = 1.5 / (1 + exp(-xi)) - 0.5;
// theta = thetavec[j];
// 
// if (fabs(xi2) < xieps) {
//   if (xi2 < 0) {
//     xi2 = -xieps;
//   } else {
//     xi2 = xieps;
//   }
// }
// 
// ee2 = exp(-xi);
// ee3 = 1 + ee2;
// ee4 = 1.5/ee3;
// ee5 = ee4 - 0.5;
// ee6 = exp(lpsi);
// ee7 = y - mu;
// ee8 = ee5 * ee7;
// ee9 = ee8/ee6;
// ee10 = ee9 + 1;
// ee11 = 1/ee5;
// ee12 = 1 + ee11;
// ee13 = exp(theta);
// ee14 = R_pow(ee10, ee12);
// ee15 = R_pow(ee3, 2);
// ee16 = log1p(ee9);
// ee17 = R_pow(ee10, ee11);
// ee18 = ee10 * ee6;
// ee19 = ee14 * ee6;
// ee21 = R_pow(ee10, ee11 + 2) * ee6;
// ee22 = ee13/ee17;
// ee23 = ee12 * ee5;
// ee24 = ee12 * ee7;
// ee25 = ee8/ee18;
// ee26 = R_pow(ee5, 2);
// ee30 = ee3 * ee5;
// ee35 = 1.5 * (ee16/(ee14 * ee26)) - 1.5 * (ee24/ee21);
// ee37 = 1.5 * (ee16/(ee17 * ee5)) - 1.5 * (ee7/ee19);
// ee38 = ((ee12 * (1.5 - 1.5 * ee25) - 1.5/ee5)/ee10 - ee35 *  ee13) * ee2;
// ee43 = ee23 * ee7/ee21;
// ee44 = ee15 * ee5;
// ee45 = ee15 * ee6;
// ee46 = ee37 * ee2;
// ee47 = 1/ee14;
// ee49 = 2.25 * (ee2 * ee7/(ee10 * ee15 * ee6)) - ((4.5/ee30 -  3) * ee2/ee3 + 1.5) * ee16;
// 
// out(j, 0) = -((ee23 - ee22)/ee18);
// out(j, 1) = (ee22 - ee23) * ee7/ee18 + 1;
// out(j, 2) = ((ee37 * ee13 - 1.5 * (ee16/ee5))/ee5 + 1.5 * (ee24/ee18)) * ee2/ee15;
// out(j, 3) = ee22 - 1;
// out(j, 4) =  - (ee12 * (ee4 - (0.5 + ee22)) * ee5/(R_pow(ee10, 2) * R_pow(ee6, 2)));
// out(j, 5) = ((ee43 - ee47) * ee13 + (1 - ee25) * ee12 * ee5/ee10)/ee6;
// out(j, 6) = -(ee38/ee45);
// out(j, 7) = ee13/ee19;
// out(j, 8) = -(((ee25 - 1) * ee12 * ee5/ee10 + (ee47 - ee43) * ee13) * ee7/ee6);
// out(j, 9) = -(ee38 * ee7/ee45);
// out(j, 10) = ee13 * ee7/ee19;
// out(j, 11) = ((((((2.25/ee30 - 3) * ee2/ee3 + 1.5)/ee14 - 1.5 * (ee35 * ee2/ee15)) * ee7/ee6 +
//    (ee49/ee17 + 1.5 * (ee46 * ee16/ee44))/ee5) * ee13 -
//    ee49/ee5)/ee5 - (((2.25 * (ee7/(ee10 * ee3 * ee6)) -
//    3) * ee2/ee3 + 1.5) * ee12 + 2.25 * (ee2/(ee15 * ee26))) * ee7/ee18) * ee2/ee15;
// out(j, 12) = ee46 * ee13/ee44;
// out(j, 13) = ee22;
// 
// }
// 
// return(out);
// 
// }
// 
// // shape = xi
// // dep = exp(theta)
// 
// double ldgev_log_id(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec)
// {
//     
// int nobs = yvec.size();
// 
// double y, mu, lpsi, xi;
// double ee1, ee2;
// double nllh = 0.0;
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
//     nllh = 1e20;
//     break;
// } else {
// 
// ee2 = 1.0 / xi;
// 
// nllh += lpsi + (ee2 + 1.0) * log1p(ee1) + R_pow(1.0 + ee1, -ee2);
// }
// 
// } else {
// 
// ee1 = (y - mu) / exp(lpsi);
// nllh += lpsi + ee1 + exp(-ee1);
//     
// }
// 
// }
// 
// return nllh;
// 
// }
// 
// // shape = xi
// // dep = exp(theta)
// 
// arma::mat ldgev12_log_id(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec)
// {
//     
// int nobs = yvec.size();
// arma::mat out(nobs, 14);
// 
// double y, mu, lpsi, xi;
// double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
// double ee10, ee11, ee12, ee13, ee16, ee17, ee18, ee19;
// double ee20, ee22, ee23;
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
// ee7 = 1 + ee6;
// ee8 = R_pow(ee5, ee6);
// ee9 = ee5 * ee1;
// ee10 = 1/ee8;
// ee11 = log1p(ee4);
// ee12 = R_pow(ee5, ee7);
// ee13 = ee7 * ee2;
// ee16 = (ee10 - xi) * ee2/ee9 + 1;
// ee17 = ee11/(xi * ee8);
// ee18 = xi * ee7;
// ee19 = (ee16 * ee7 - (1 + ee17)/xi)/ee5;
// ee20 = ee13/ee9;
// ee22 = ee2/(ee12 * ee1);
// ee23 = xi - ee10;
// 
// out(j, 0) = -((ee18 - ee10)/ee9);
// out(j, 1) = (ee10 - ee18) * ee2/ee9 + 1;
// out(j, 2) = ((ee10 - 1) * ee11/xi - ee22)/xi + ee20;
// out(j, 3) = 0;
// out(j, 4) =  - (ee18 * ee23/(R_pow(ee5, 2) * R_pow(ee1, 2)));
// out(j, 5) = (xi * ee16 * ee7 - ee10)/ee5/ee1;
// out(j, 6) = -(ee19/ee1);
// out(j, 7) = 0;
// out(j, 8) = -((ee10 + xi * (ee23 * ee2/ee9 - 1) * ee7)/ee5 * ee2/ee1);
// out(j, 9) = -(ee19 * ee2/ee1);
// out(j, 10) = 0;
// out(j, 11) = ((((ee2/ee9 - 2 * (ee11/xi))/R_pow(ee5, ee6 - 1) -
//    ee2/ee1)/ee5 + (2 + ee17 - ee22) * ee11/xi)/xi + (ee13/(R_pow(ee5, ee6 +
//    2) * ee1) + (1/ee12 - ee11/(xi * ee12))/xi) * ee2/ee1)/xi -
//    (ee20 + 1/R_pow(xi, 2)) * ee2/ee9;
// out(j, 12) = 0;
// out(j, 13) = 0;
// 
// } else {
//     
// ee1 = exp(lpsi);
// ee2 = y - mu;
// ee3 = ee2/ee1;
// ee5 = exp(-ee3);
// ee7 = (ee3 - 1) * ee5 + 1;
// ee8 = ee5 - 1;
// 
// out(j, 0) = ee8/ee1;
// out(j, 1) = ee8 * ee2/ee1 + 1;
// out(j, 2) = 0;
// out(j, 3) = 0;
// out(j, 4) = ee5/R_pow(ee1, 2);
// out(j, 5) = ee7/ee1;
// out(j, 6) = 0;
// out(j, 7) = 0;
// out(j, 8) = ee7 * ee2/ee1;
// out(j, 9) = 0;
// out(j, 10) = 0;
// out(j, 11) = 0;
// out(j, 12) = 0;
// out(j, 13) = 0;
//     
// }
// 
// }
// 
// return(out);
// 
// }
// 
// double ldgevagg_log_id2(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec, arma::vec thetavec)
// {
//     
// int nobs = yvec.size();
// 
// double y, mu, lpsi, psi, xi, theta;
// double ee1, ee2;
// double nllh = 0.0;
// 
// for (int j=0; j < nobs; j++) {
// 
// y = yvec[j];
// mu = muvec[j];
// lpsi = lpsivec[j];
// xi = xivec[j];
// theta = thetavec[j];
// psi = exp(lpsi);
// 
// if (fabs(xi) >= xieps) {
// 
// mu = mu - psi * (1 - exp(theta * xi)) / xi;
// lpsi = lpsi + xi * theta;
// 
// ee1 = xi * (y - mu) / exp(lpsi);
// 
// if (ee1 <= -1.0) {
//     nllh = 1e20;
//     break;
// } else {
// 
// ee2 = 1.0 / xi;
// 
// nllh += lpsi + (ee2 + 1.0) * log1p(ee1) + R_pow(1.0 + ee1, -ee2);
// }
// 
// } else {
// 
// mu = mu + psi * theta;
// 
// ee1 = (y - mu) / exp(lpsi);
// nllh += lpsi + ee1 + exp(-ee1);
//     
// }
// 
// }
// 
// return nllh;
// 
// }
// 
// arma::mat ldgevagg12_log_id2(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec, arma::vec thetavec)
// {
//     
// int nobs = yvec.size();
// arma::mat out(nobs, 14);
// 
// double y, mu, lpsi, xi, theta;
// double ee1, ee2, ee3, ee4, ee5, ee7, ee8, ee9;
// double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee18, ee19;
// double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee28, ee29;
// double ee30, ee31, ee32, ee33, ee34, ee35, ee37, ee39;
// double ee41, ee42, ee43, ee44, ee49;
// double ee51, ee52, ee53;
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
// ee1 = theta * xi;
// ee2 = exp(lpsi);
// ee3 = exp(ee1);
// ee4 = (1 - ee3) * ee2;
// ee7 = ee4/xi + y - mu;
// ee8 = xi * ee7;
// ee9 = exp(lpsi + ee1);
// ee10 = ee8/ee9;
// ee11 = 1 + ee10;
// ee12 = 1/xi;
// ee13 = 1 + ee12;
// ee14 = ee2 * ee3;
// ee15 = ee14 + ee8;
// ee16 = theta * ee15;
// ee18 = y - (mu + ee16);
// ee19 = ee11 * ee9;
// ee20 = R_pow(ee11, ee13);
// ee21 = R_pow(ee11, ee12);
// ee22 = ee4 - ee8;
// ee23 = ee20 * ee9;
// ee24 = 1/ee21;
// ee25 = log1p(ee10);
// ee26 = ee13 * ee18;
// ee28 = R_pow(ee11, ee12 + 2) * ee9;
// ee29 = ee15/ee23;
// ee30 = R_pow(xi, 2);
// ee31 = ee26/ee28;
// ee32 = ee15/ee19;
// ee33 = ee18/ee19;
// ee34 = ee18/ee23;
// ee35 = ((ee29 + xi * (1 - ee32)) * ee13 - ee24)/ee11;
// ee37 = ee22/ee19 + 1;
// ee39 = ee7 * (1 - ee1) + y;
// ee41 = (1 + ee1) * ee2 * ee3;
// ee42 = ee33 + theta;
// ee43 = 1/ee20;
// ee44 = 2 * (theta * ee2 * ee3);
// ee49 = ee25/(xi * ee21);
// ee51 = ee25/(ee30 * ee20) - ee31;
// ee52 = theta * (2 * y - (2 * mu + ee16));
// ee53 = xi * ee13;
// 
// out(j, 0) = -((ee53 - ee24)/ee19);
// out(j, 1) = ((1 - ee24)/xi + 1) * ee22/ee19 + 1;
// out(j, 2) = ((ee24 - 1) * ee25/xi - ee34)/xi + ee26/ee19 + theta;
// out(j, 3) = ee29 + xi * (1 - ee13 * ee15/ee19);
// out(j, 4) =  - (ee53 * (xi - ee24)/(R_pow(ee11, 2) * R_pow(ee9, 2)));
// out(j, 5) = -(((ee22/ee23 - xi * ee37) * ee13 + ee24)/ee11/ee9);
// out(j, 6) = -(((ee34 + 1 - xi * ee42) * ee13 + theta/ee21 -
//    (1 + ee49)/xi)/ee11/ee9);
// out(j, 7) = xi * ee35/ee9;
// out(j, 8) = ((ee22 * ee13/ee28 + ee43)/xi - ee37 * ee13/ee11) * ee22/ee9;
// out(j, 9) = (((ee22 * (ee12 + theta) + y - mu)/ee20 - ee22 * ee51)/xi -
//    ((ee22 * ee42 + y - mu) * ee13 + ee22/ee30)/ee11)/ee9;
// out(j, 10) = -(ee35 * ee22/ee9);
// out(j, 11) = ((((ee33 - 2 * (ee25/xi))/R_pow(ee11, ee12 - 1) -
//    ee18/ee9)/ee11 + (2 + ee49 - ee34) * ee25/xi)/xi + ((ee18/xi +
//    ee52)/ee20 - ee51 * ee18)/ee9)/xi - ((R_pow(ee18, 2)/ee19 +
//    ee52) * ee13 + ee18/ee30)/ee19;
// out(j, 12) = 1 - (((ee41 + xi * (ee39 - (ee15 * ee18/ee19 +
//    ee44 + mu))) * ee13 - ((ee41 + xi * (ee39 - (ee44 + mu)))/ee21 +
//    ee14 + ee8)/xi)/ee11 + (ee31 + (ee43 - ee25/(xi * ee20))/xi) * ee15)/ee9;
// out(j, 13) = -(xi * ((ee13 * (xi * (ee32 - 1) - ee29) + ee24)/ee11) * ee15/ee9);
// 
// } else {
// 
// ee1 = exp(lpsi);
// ee2 = theta * ee1;
// ee4 = y - (mu + ee2);
// ee5 = ee4/ee1;
// ee7 = exp(-ee5);
// ee8 = ee5 + theta;
// ee9 = ee8 * ee7;
// ee10 = (y - (2 * ee2 + 2 * ee4 + mu))/ee1;
// ee11 = ee7 - 1;
// 
// out(j, 0) = ee11/ee1;
// out(j, 1) = ee9 + 1 - ee8;
// out(j, 2) = 0;
// out(j, 3) = ee11;
// out(j, 4) = ee7/R_pow(ee1, 2);
// out(j, 5) = ((ee8 - 1) * ee7 + 1)/ee1;
// out(j, 6) = 0;
// out(j, 7) = ee7/ee1;
// out(j, 8) = (R_pow(ee8, 2) + ee10) * ee7 - ee10;
// out(j, 9) = 0;
// out(j, 10) = ee9;
// out(j, 11) = 0;
// out(j, 12) = 0;
// out(j, 13) = ee7;
//     
// }
// 
// }
// 
// return(out);
// 
// }

// [[Rcpp::export]]
double ldgevagg(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec, arma::vec thetavec)
{
// double out = ldgevagg_log(yvec, muvec, lpsivec, xivec, thetavec);
double out = ldgevagg_logit(yvec, muvec, lpsivec, xivec, thetavec);
return out;
}

// [[Rcpp::export]]
arma::mat ldgevagg12(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec, arma::vec thetavec)
{
// arma::mat out = ldgevagg12_log(yvec, muvec, lpsivec, xivec, thetavec);
arma::mat out = ldgevagg12_logit(yvec, muvec, lpsivec, xivec, thetavec);
return out;
}
