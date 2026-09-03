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
arma::mat out(nobs, 14, arma::fill::zeros);

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

double y, mu, lpsi, xi, lgttheta, theta;
double ee1, ee2;
double nllh = 0.0;
double theta2 = 1.0;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
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

mu = mu + exp(lpsi) * log(theta2);

ee1 = (y - mu) / exp(lpsi);
nllh += lpsi + ee1 + exp(-ee1);
    
}

}

return nllh;

}

// [[Rcpp::export]]
arma::mat ldgevagg12_logit(arma::vec yvec, arma::vec muvec, arma::vec lpsivec, arma::vec xivec, arma::vec thetavec)
{
    
int nobs = yvec.size();
arma::mat out(nobs, 14, arma::fill::zeros);

double y, mu, lpsi, xi, lgttheta, theta;

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
