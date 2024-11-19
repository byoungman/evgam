// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

double xi_from_zero_v2(double xi, double eps) 
{
if (fabs(xi) <= eps) {
  if (xi >= 0.0) {
    xi = eps;
  } else {
    xi = -eps;
  }
}
return xi;
}

const double xieps2 = 0.0001;

double xi2txi_v2(double xi) 
{
return -log(1.5 / (1.0 + xi) - 1.0);
}

// //' R-largest generalized extreme value (GEV) distribution negative log-likelihood with constrained shape parameter
// //'
// //' @param pars a list of vectors of coefficients for each GEV parameter
// //' @param X1 a design matrix for the GEV location parameter
// //' @param X2 a design matrix for the GEV log scale parameter
// //' @param X3 a design matrix for the GEV transformed shape parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return rlarged0 a scalar, the negative log-likelihood
// //' @return rlarged12 a matrix, first then second derivatives w.r.t. GEV parameters
// //' @return rlarged34 a matrix, third then fourth derivatives w.r.t. GEV parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double rlarged0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = ymat.n_rows;
int r = ymat.n_cols;

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    txivec = txivec.elem(dupid);
}

double y, mu, lpsi, txi, xi;
double ee1;
double nllh = 0.0;

for (int j=0; j < nobs; j++) {

mu = muvec[j];
lpsi = lpsivec[j];
txi = txivec[j];
xi = 1.5 / (1.0 + exp(-txi)) - 1.0;

xi = xi_from_zero_v2(xi, xieps2);  
txi = xi2txi_v2(xi);

for (int l=0; l < r; l++) {

y = ymat(j, l);

ee1 = xi * (y - mu) / exp(lpsi);

if (ee1 <= -1.0) {
    nllh = 1e20;
    break;
} else {

nllh += (1 + 1/xi) * log1p(ee1) + lpsi;

}

if (l == r - 1) {
  
  nllh += R_pow(1 + ee1, -1.0 / xi);
  
}

} 

}

return(nllh);

}

// //' @rdname rlarged0
// [[Rcpp::export]]
arma::mat rlarged12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate)
{
    
  arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = ymat.n_rows;
  int r = ymat.n_cols;
  
  if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    txivec = txivec.elem(dupid);
  }
  
  double y, mu, lpsi, txi, xi;
  double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee17, ee18, ee19;
  double ee20;
  
  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    mu = muvec[j];
    lpsi = lpsivec[j];
    txi = txivec[j];
    xi = 1.5 / (1.0 + exp(-txi)) - 1.0;
    
    xi = xi_from_zero_v2(xi, xieps2);  
    txi = xi2txi_v2(xi);
    
    for (int l=0; l < r; l++) {
      
      y = ymat(j, l);
      
      ee2 = exp(-txi);
      ee3 = 1 + ee2;
      ee5 = 1.5/ee3 - 1;
      ee6 = exp(lpsi);
      ee7 = y - mu;
      ee8 = ee5 * ee7;
      ee9 = ee8/ee6;
      ee10 = ee9 + 1;
      ee11 = ee10 * ee6;
      ee12 = 1 + 1/ee5;
      ee13 = R_pow(ee3, 2);
      ee14 = ee8/ee11;
      ee15 = R_pow(ee5, 2);
      ee17 = ee10 * ee13 * ee6;
      ee18 = (ee12 * (1.5 - 1.5 * ee14) - 1.5/ee5) * ee2;
      ee19 = ee12 * ee5;
      ee20 = log1p(ee9);
      
      out(j, 0) += -(ee19/ee11);
      out(j, 1) += 1 - ee19 * ee7/ee11;
      out(j, 2) += (1.5 * (ee12 * ee7/ee11) - 1.5 * (ee20/ee15)) *
        ee2/ee13;
      out(j, 3) += -(ee12 * ee15/(R_pow(ee10, 2) * R_pow(ee6, 2)));
      out(j, 4) += (1 - ee14) * ee12 * ee5/ee11;
      out(j, 5) += -(ee18/ee17);
      out(j, 6) += -((ee14 - 1) * ee12 * ee5 * ee7/ee11);
      out(j, 7) += -(ee18 * ee7/ee17);
      out(j, 8) += -(((((2.25 * (ee7/(ee10 * ee3 * ee6)) - 3) * ee2/
        ee3 + 1.5) * ee12 + 2.25 * (ee2/(ee13 * ee15))) * ee7/ee11 +
          (2.25 * (ee2 * ee7/ee17) - ((4.5/(ee3 * ee5) - 3) * ee2/ee3 +
          1.5) * ee20)/ee15) * ee2/ee13);
      
      if (l == r - 1) {
        
        double ee2, ee3, ee5, ee6, ee7, ee9;
        double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
        double ee20, ee21, ee24, ee25, ee26, ee29;
        double ee30, ee31, ee32, ee33;
        
        ee2 = exp(-txi);
        ee3 = 1 + ee2;
        ee5 = 1.5/ee3 - 1;
        ee6 = exp(lpsi);
        ee7 = y - mu;
        ee9 = ee5 * ee7/ee6;
        ee10 = 1/ee5;
        ee11 = ee9 + 1;
        ee12 = 1 + ee10;
        ee13 = R_pow(ee11, ee12);
        ee14 = R_pow(ee3, 2);
        ee15 = log1p(ee9);
        ee16 = R_pow(ee11, (ee10 + 2));
        ee17 = ee16 * ee6;
        ee18 = ee13 * ee6;
        ee20 = R_pow(ee11, ee10);
        ee21 = ee12 * ee5;
        ee24 = ee14 * ee5;
        ee25 = (1.5 * (ee15/(ee13 * R_pow(ee5, 2))) - 1.5 * (ee12 * ee7/ee17)) * ee2;
        ee26 = ee7/ee18;
        ee29 = ee21 * ee7/ee17;
        ee30 = ee3 * ee5;
        ee31 = ee14 * ee6;
        ee32 = (1.5 * (ee15/(ee20 * ee5)) - 1.5 * ee26) * ee2;
        ee33 = 1/ee13;
        
        out(j, 0) += 1/ee18;
        out(j, 1) += ee26;
        out(j, 2) += ee32/ee24;
        out(j, 3) += ee21/(ee16 * R_pow(ee6, 2));
        out(j, 4) += (ee29 - ee33)/ee6;
        out(j, 5) += ee25/ee31;
        out(j, 6) += -((ee33 - ee29) * ee7/ee6);
        out(j, 7) += ee25 * ee7/ee31;
        out(j, 8) += ((((2.25/ee30 - 3) * ee2/ee3 + 1.5)/ee13 - 1.5 *
          (ee25/ee14)) * ee7/ee6 + ((2.25 * (ee2 * ee7/(ee11 * ee14 *
          ee6)) - ((4.5/ee30 - 3) * ee2/ee3 + 1.5) * ee15)/ee20 + 1.5 *
          (ee32 * ee15/ee24))/ee5) * ee2/ee24;
        
      }
      
    }
      
}

return out;

}

// //' @rdname rlarged0
// [[Rcpp::export]]
arma::mat rlarged34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate)
{
    
  arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = ymat.n_rows;
  int r = ymat.n_cols;
  
  if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    txivec = txivec.elem(dupid);
  }
  
  double y, mu, lpsi, txi, xi;
  
arma::mat out = arma::mat(nobs, 25, arma::fill::zeros);

double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee20, ee21, ee23, ee26, ee27, ee28, ee29;
double ee30, ee32, ee33, ee35, ee38, ee39;
double ee40, ee44, ee45, ee46, ee47, ee48;
double ee50, ee51, ee52, ee54, ee55, ee56, ee57, ee58, ee59;
double ee61, ee64, ee66, ee67, ee69;
double ee74, ee77, ee78, ee79;
double ee80, ee81, ee84, ee85, ee87, ee89;
double ee90, ee91, ee92, ee93, ee94, ee95, ee97, ee98, ee99;
double ee100, ee103, ee104, ee106, ee107, ee108;
double ee112, ee116, ee117, ee118, ee119;
double ee121, ee122, ee123, ee124, ee125, ee126, ee127, ee128, ee129;

double ee22, ee37, ee42, ee43, ee53, ee60, ee62, ee63, ee65;
double ee70, ee71, ee75, ee82, ee88;
double ee101, ee111, ee114, ee120, ee130, ee131, ee133, ee136, ee139;
double ee144, ee146, ee147, ee152, ee153, ee154, ee155, ee156, ee159;
double ee160, ee161, ee165, ee166, ee167, ee168, ee170;
double ee174, ee175, ee176, ee178, ee179, ee180;

for (int j=0; j < nobs; j++) {

mu = muvec[j];
lpsi = lpsivec[j];
txi = txivec[j];
xi = 1.5 / (1.0 + exp(-txi)) - 1.0;

xi = xi_from_zero_v2(xi, xieps2);  
txi = xi2txi_v2(xi);

for (int l=0; l < r; l++) {

y = ymat(j, l);

  ee2 = exp(-txi);
  ee3 = 1 + ee2;
  ee5 = 1.5/ee3 - 1;
  ee6 = exp(lpsi);
  ee7 = y - mu;
  ee8 = ee5 * ee7;
  ee9 = ee8/ee6;
  ee10 = ee9 + 1;
  ee11 = R_pow(ee3, 2);
  ee12 = ee10 * ee6;
  ee13 = 1.5 - 3 * (ee2/ee3);
  ee14 = ee2/ee11;
  ee15 = 1 + 1/ee5;
  ee16 = 1 + 2 * ee2;
  ee17 = ee13 * ee5;
  ee18 = 2.25 * ee14;
  ee19 = ee8/ee12;
  ee20 = 3 * ee16;
  ee21 = 3 * ee3;
  ee23 = ee10 * ee11 * ee6;
  ee26 = ee7/(ee10 * ee3 * ee6);
  ee27 = 12 * ee2;
  ee28 = ee2 * ee7;
  ee29 = ee18 - ee17;
  ee30 = ee20 + ee21;
  ee32 = ee28/(ee11 * ee6);
  ee33 = ee3 * ee5;
  ee35 = (2.25 * ee26 - 3) * ee2/ee3;
  ee38 = (4.5/ee33 - 3) * ee2/ee3 + 1.5;
  ee39 = 3 * ee29;
  ee40 = ee10 * ee13;
  ee44 = (4.5 * ee26 - 3) * ee2/ee3;
  ee45 = 1.5 - ((ee30 - ee27)/ee3 + 3) * ee2/ee3;
  ee46 = 4.5 * ee13;
  ee47 = R_pow(ee5, 2);
  ee48 = 2 * ee19;
  ee50 = 2.25 * ee32 - ee40;
  ee51 = 3 * ee17;
  ee52 = 3 * ee19;
  ee54 = (((ee46 + 6.75 * (ee28/ee23)) * ee7/ee12 + ee27 -  ee30)/ee3 - 3) * ee2/ee3;
  ee55 = R_pow(ee10, 2);
  ee56 = (ee39 - (27 * ee14 + ee51))/ee5;
  ee57 = ee44 + 1.5;
  ee58 = 3 * ee13;
  ee59 = R_pow(ee6, 2);
  ee61 = ee45 * ee5;
  ee64 = ee5 * (4.5 - ee52) * ee7/ee12;
  ee66 = ee35 + (ee18 - ee57 * ee5) * ee7/ee12 + 1.5;
  ee67 = ee35 + 1.5;
  ee69 = 1.5 - 1.5 * ee19;
  ee74 = 2 * (1 + 2 * ee9) + 4 * ee10;
  ee77 = 2 * (1.5 * ee10 + 1.5 * ee9) + 6 * ee10 - 12 * ee9;
  ee78 = 2 * ee29;
  ee79 = 4 * ee2;
  ee80 = 8 * ee9;
  ee81 = 9 * ee14;
  ee84 = (((ee56 - ee58)/ee5 + ee20 + ee21 - ee27)/ee3 + 3) *  ee2/ee3 - 1.5;
  ee85 = ee54 + 1.5;
  ee87 = (ee5 * (2 * ee50 - 18 * ee32) + 9 * (ee10 * ee2/ee11))/ee10 +  ee81;
  ee89 = ee55 * ee11 * ee59;
  ee90 = ee15 * ee5;
  ee91 = ee11 * ee5;
  ee92 = ee77/ee10;
  ee93 = (ee80 - ee74)/ee10;
  ee94 = R_pow(ee6, 3);
  ee95 = ((((((ee5 * (3 * ee50 - 27 * ee32) + 3 * (ee10 *  ee29))/ee10 + ee39) * ee7/ee12 + ee46) * ee2/ee11 - ee61) *  ee7/ee12 + ee54 + 1.5) * ee15 + ((((ee56 + 3 * (ee38 *  ee69) + 4.5 * ee66 - ee58)/ee5 + ee20 + ee21 - ee27)/ee3 +  3) * ee2/ee3 - 1.5)/ee5) * ee2;
  ee97 = ((((3 * ee69 - 4.5)/ee33 + 3) * ee2/ee3 - 1.5)/ee5 +  ee66 * ee15) * ee2;
  ee98 = ee55 * ee59;
  ee99 = R_pow(ee10, 3);
  ee100 = R_pow(ee10, 4);
  ee103 = (ee92 + 9) * ee5 * ee7/ee12;
  ee104 = ((ee78 + ee18 - (ee87 * ee7/ee12 + ee44 + 1.5) *  ee5) * ee7/ee12 + ee35 + 1.5) * ee15;
  ee106 = ((6 * (2 * ee16 + ee79) + 8 * ee30 - 96 * ee2)/ee3 +  12) * ee2;
  ee107 = ee15 * R_pow(ee5, 3);
  ee108 = ee11 * ee47;
  ee112 = ee5 * (ee48 - 3) * ee7/ee12 + 1;
  ee116 = ee5 * (3 - ee48) * ee7/ee12 - 1;
  ee117 = ee64 - 1.5;
  ee118 = ee19 - 1;
  ee119 = 1 - ee19;
  ee121 = 1.5 - ee64;
  ee122 = 2 - ee48;
  ee123 = 3 - ee52;
  ee124 = 3 * (ee67 * ee38);
  ee125 = 3 * (1 + ee79);
  ee126 = 3 * ee45;
  ee127 = 6 * ee3;
  ee128 = 9 * ee16;
  ee129 = log1p(ee9);
  
  out(j, 0) += -(2 * (ee107/(ee99 * ee94)));
  out(j, 1) += ee15 * ee47 * ee122/ee98;
  out(j, 2) += -((ee90 * ee123 - 1.5) * ee2/ee89);
  out(j, 3) += ee116 * ee15 * ee5/ee12;
  out(j, 4) += (ee15 * ee121 - 1.5 * (ee119/ee5)) * ee2/ee23;
  out(j, 5) += ee97/ee23;
  out(j, 6) += -(ee112 * ee15 * ee5 * ee7/ee12);
  out(j, 7) += -((ee117 * ee15 - 1.5 * (ee118/ee5)) * ee2 * ee7/
    ee23);
  out(j, 8) += ee97 * ee7/ee23;
  out(j, 9) += ((ee84 * ee129 + (1.5 * ee67 + 3 * ee38) * ee2 *
    ee7/ee23)/ee47 + (ee85 * ee15 + (1.5 * ee38 + 3 * ee67) * ee2/
      ee108) * ee7/ee12) * ee2/ee11;
  out(j, 10) += -(6 * (ee15 * R_pow(ee5, 4)/(ee100 * R_pow(ee6, 4))));
  out(j, 11) += ee107 * (ee74 - ee80)/(ee100 * ee94);
  out(j, 12) += -((ee90 * ee77/ee10 - 3) * ee5 * ee2/(ee99 * ee11 *
    ee94));
  out(j, 13) += (ee5 * (4 - ee93) * ee7/ee12 - 4) * ee15 * ee47/
    ee98;
  out(j, 14) += (ee90 * (6 - (ee92 + 6) * ee5 * ee7/ee12) - 1.5 *
    ee122) * ee2/ee89;
  out(j, 15) += -((ee15 * (ee78 - ee87 * ee5 * ee7/ee12) + 1.5 -
    ((3 * ee123 - 4.5)/ee33 + 3) * ee2/ee3) * ee2/ee89);
  out(j, 16) += ((ee5 * (6 - ee93) * ee7/ee12 - 7) * ee5 * ee7/
    ee12 + 1) * ee15 * ee5/ee12;
  out(j, 17) += ((ee5 * (10.5 - ee103) * ee7/ee12 - 1.5) * ee15 -
    1.5 * (ee116/ee5)) * ee2/ee23;
  out(j, 18) += -((ee104 + (3 * (ee121 * ee2/ee91) - ee38 * ee119)/
    ee5) * ee2/ee23);
  out(j, 19) += -(ee95/ee23);
  out(j, 20) += -((((ee93 - 6) * ee5 * ee7/ee12 + 7) * ee5 * ee7/
    ee12 - 1) * ee15 * ee5 * ee7/ee12);
  out(j, 21) += -((((ee103 - 10.5) * ee5 * ee7/ee12 + 1.5) * ee15 -
    1.5 * (ee112/ee5)) * ee2 * ee7/ee23);
  out(j, 22) += -((ee104 - (3 * (ee117 * ee2/ee91) - ee118 * ee38)/
    ee5) * ee2 * ee7/ee23);
  out(j, 23) += -(ee95 * ee7/ee23);
  out(j, 24) += ((((((((4.5 * ee50 - (40.5 * ee32 + 9 * ee40))/
    ee10 - 9 * ee13) * ee2 * ee7/ee23 - (ee57 * ee13 + 2 * (R_pow(ee13, 2) +
      1.5 * ee45) + ee126)) * ee7/ee12 + ee125 + ee127 +
      ee128 - ee106)/ee3 + 3) * ee2/ee3 - 1.5) * ee15 + (1.5 *
      ee84 - (ee124 + 4.5 * ee85)) * ee2/ee108) * ee7/ee12 + ((4.5 *
      ee84 - (1.5 * ee85 + ee124)) * ee2 * ee7/ee23 - (((((ee13 *
      (6 * ee29 - 18 * ee14) + (12 * (ee39 - ee51) + 9 * (ee78 +
      ee81) - 324 * ee14) * ee2/ee91 + 3 * (4.5 * (ee13 * ee2/
        ee11) - ee61) - 6 * ee61)/ee5 - ee126)/ee5 + 
          ee125 + ee127 +
          ee128 - ee106)/ee3 + 3) * ee2/ee3 - 1.5) * ee129)/ee47) *
          ee2/ee11;
  
  if (l == r - 1) {
    
    ee2 = exp(-txi);
    ee3 = 1 + ee2;
    ee5 = 1.5/ee3 - 1;
    ee6 = exp(lpsi);
    ee7 = y - mu;
    ee9 = ee5 * ee7/ee6;
    ee10 = 1/ee5;
    ee11 = ee9 + 1;
    ee12 = ee10 + 2;
    ee13 = R_pow(ee3, 2);
    ee14 = 1 + ee10;
    ee15 = log1p(ee9);
    ee16 = R_pow(ee5, 2);
    ee17 = R_pow(ee11, ee12);
    ee18 = ee10 + 3;
    ee19 = 3 * (ee2/ee3);
    ee20 = 1.5 - ee19;
    ee21 = R_pow(ee11, ee18);
    ee22 = R_pow(ee11, ee14);
    ee23 = ee3 * ee5;
    ee26 = (4.5/ee23 - 3) * ee2/ee3 + 1.5;
    ee28 = ee11 * ee13 * ee6;
    ee35 = 1.5 * (ee15/(ee17 * ee16)) - 1.5 * (ee12 * ee7/(ee21 *  ee6));
    ee37 = ee2 * ee7/ee28;
    ee38 = ee2/ee13;
    ee40 = 1 + 2 * ee2;
    ee42 = 2.25 * ee37 - ee26 * ee15;
    ee43 = ee20 * ee5;
    ee45 = 1.5 * (ee15/(ee22 * ee16));
    ee46 = 3 * ee40;
    ee47 = 3 * ee3;
    ee48 = ee13 * ee16;
    ee53 = ee45 - 1.5 * (ee14 * ee7/(ee17 * ee6));
    ee54 = 12 * ee2;
    ee55 = 2.25 * (ee2/ee48);
    ee56 = 1.5/ee5;
    ee57 = ee14 * ee20;
    ee58 = 3 * ee20;
    ee59 = ee11 * ee6;
    ee60 = ee53 * ee2;
    ee62 = 2.25 * ee38 - ee43;
    ee63 = ee46 + ee47;
    ee64 = R_pow(ee11, ee10);
    ee65 = ee14 * ee35;
    ee66 = R_pow(ee11, (ee10 + 4));
    ee67 = ee57 + ee55;
    ee69 = ee56 - 1.5 * ee14;
    ee70 = 3 * ee43;
    ee71 = 3 * ee62;
    ee74 = (ee42/ee22 + 1.5 * (ee60 * ee15/ee13))/ee16;
    ee75 = ee13 * ee6;
    ee78 = (ee71 - (27 * ee38 + ee70))/ee5;
    ee81 = 1.5 * (ee15/(ee21 * ee16)) - 1.5 * (ee18 * ee7/(ee66 *  ee6));
    ee82 = ee13 * ee5;
    ee85 = (((ee78 - ee58)/ee5 + ee46 + ee47 - ee54)/ee3 + 3) *  ee2/ee3 - 1.5;
    ee88 = (ee67/ee17 - 1.5 * (ee65 * ee2/ee13)) * ee7/ee6 +  ee74;
    ee89 = ee35 * ee2;
    ee92 = ee5 * ee12 * ee7/ee59;
    ee93 = 2.25/ee23;
    ee100 = (2.25 * (ee7/(ee11 * ee3 * ee6)) - 3) * ee2/ee3 +  1.5;
    ee101 = 1.5 - ((ee63 - ee54)/ee3 + 3) * ee2/ee3;
    ee103 = ee85 * ee15 + (1.5 * ee100 + 3 * ee26) * ee2 * ee7/ee28;
    ee106 = ((ee20 * ee12 + ee55)/ee21 - 1.5 * (ee81 * ee12 *  ee2/ee13)) * ee7/ee6 + (ee42/ee17 + 1.5 * (ee89 * ee15/ee13))/ee16;
    ee111 = (ee93 - 3) * ee2/ee3 + 1.5;
    ee114 = 1.5 * ee26 + ee58;
    ee117 = 1.5 * (ee15/(ee64 * ee5)) - 1.5 * (ee7/(ee22 * ee6));
    ee118 = R_pow(ee6, 2);
    ee119 = ee106 * ee14;
    ee120 = (ee35 * ee5 + 1.5/ee17) * ee14;
    ee122 = ((ee56 - 1.5 * ee12)/ee21 - ee81 * ee5 * ee12) *  ee14 + ee69 * ee12/ee21;
    ee123 = ee14 * (3 - ee92);
    ee127 = ee69/ee17;
    ee128 = 1 + ee14;
    ee129 = 3 - ee5 * ee18 * ee7/ee59;
    ee130 = 4 * ee2;
    ee131 = ee103/ee22;
    ee133 = (ee122 * ee7/ee6 + 2 * ee65) * ee5;
    ee136 = ((ee114/ee5 + ee54 - ee63)/ee3 - 3) * ee2/ee3 +  1.5;
    ee139 = (ee111/ee22 - 1.5 * (ee60/ee13)) * ee7/ee6 + (ee42/ee64 +  1.5 * (ee117 * ee2 * ee15/ee82))/ee5;
    ee144 = (ee12 * ee129 + ee10) * ee5 * ee7/ee6;
    ee146 = ((2 * (ee67 * ee35) - 1.5 * ee119) * ee2/ee13 -  (ee14 * ee101 + ee114 * ee2/ee48)/ee17) * ee7/ee6;
    ee147 = ee13 * ee118;
    ee152 = (2 - ee92)/ee17;
    ee153 = 1.5 * (ee88 * ee15);
    ee154 = 2 * (ee53 * ee42);
    ee155 = 2 * (ee35 * ee69 * ee2/ee13);
    ee156 = R_pow(ee6, 3);
    ee159 = (((((ee93 + 3) * ee2/ee3 - 1.5)/ee5 + 2 * ee57)/ee17 +  ee155 - (ee106 * ee5 + 1.5 * (ee89/ee13)) * ee14) * ee7/ee6 +  ee74) * ee2;
    ee160 = ((((ee78 + 1.5 * ee20 + 4.5 * ee26)/ee22 + ee153 +  ee154) * ee2/ee13 - ee131)/ee16 + ee146) * ee2;
    ee161 = ee88 * ee2;
    ee165 = ((ee123 - 1) * ee5 * ee7/ee6 - 1)/ee17;
    ee166 = ((ee127 - ee120) * ee7/ee6 + ee45) * ee2;
    ee167 = (((4.5/ee5 - 4.5 * ee14)/ee17 - (ee133 + ee120)) *  ee7/ee6 + ee45) * ee2;
    ee168 = ((1 - ee123) * ee5 * ee7/ee6 + 1)/ee17;
    ee170 = ((1.5 * (ee139 * ee15) + 2 * (ee117 * ee42)) * ee2/ee82 -  ee103/ee64)/ee5 + ((2 * (ee111 * ee53) - 1.5 * ee88) *  ee2/ee13 - ee136/ee22) * ee7/ee6;
    ee174 = ((2 - ee144)/ee21 + ee152) * ee14 * ee5 * ee7/ee6;
    ee175 = ((6 * (2 * ee40 + ee130) + 8 * ee63 - 96 * ee2)/ee3 +  12) * ee2;
    ee176 = ee101 * ee5;
    ee178 = 3 * (1 + ee130);
    ee179 = 6 * ee3;
    ee180 = 9 * ee40;
    
    out(j, 0) += ee128 * ee14 * ee16/(ee21 * ee156);
    out(j, 1) += -(ee152 * ee14 * ee5/ee118);
    out(j, 2) += -((ee127 - ee65 * ee5) * ee2/ee147);
    out(j, 3) += ee168/ee6;
    out(j, 4) += -(ee166/ee75);
    out(j, 5) += ee161/ee75;
    out(j, 6) += -(ee165 * ee7/ee6);
    out(j, 7) += -(ee166 * ee7/ee75);
    out(j, 8) += ee161 * ee7/ee75;
    out(j, 9) += ee170 * ee2/ee82;
    out(j, 10) += (1 + ee128) * ee128 * ee14 * R_pow(ee5, 3)/(ee66 *
      R_pow(ee6, 4));
    out(j, 11) += -(ee129/ee21 * ee128 * ee14 * ee16/ee156);
    out(j, 12) += -(ee122 * ee5 * ee2/(ee13 * ee156));
    out(j, 13) += -(((ee144 - 2)/ee21 + (ee92 - 2)/ee17) * ee14 *
      ee5/ee118);
    out(j, 14) += (2 * ee127 - ee133) * ee2/ee147;
    out(j, 15) += -(((ee57 + (ee19 - 1.5)/ee5)/ee17 + ee155 - ee119 *
      ee5) * ee2/ee147);
    out(j, 16) += (ee165 + ee174)/ee6;
    out(j, 17) += ee167/ee75;
    out(j, 18) += -(ee159/ee75);
    out(j, 19) += ee160/ee75;
    out(j, 20) += -((ee168 - ee174) * ee7/ee6);
    out(j, 21) += ee167 * ee7/ee75;
    out(j, 22) += -(ee159 * ee7/ee75);
    out(j, 23) += ee160 * ee7/ee75;
    out(j, 24) += (((1.5 * (ee170 * ee15) + 3 * (ee139 * ee42) -
      3 * (ee103 * ee117)) * ee2/ee82 - ((4.5 * ee85 - (1.5 * ((((4.5 *
      ee20 + 6.75 * ee37) * ee7/ee59 + ee54 - ee63)/ee3 - 3) *
      ee2/ee3 + 1.5) + 3 * (ee100 * ee26))) * ee2 * ee7/ee28 -
      (((((ee20 * (6 * ee62 - 18 * ee38) + (12 * (ee71 - ee70) + 9 *
      (2 * ee62 + 9 * ee38) - 324 * ee38) * ee2/ee82 + 3 * (4.5 *
      (ee20 * ee2/ee13) - ee176) - 6 * ee176)/ee5 - 3 * ee101)/
        ee5 + ee178 + ee179 + ee180 - ee175)/ee3 + 3) * ee2/ee3 - 1.5) *
          ee15)/ee64)/ee5 + ((3 * (ee88 * ee111) - (1.5 * (((ee153 +
          ee154) * ee2/ee13 - ee131)/ee16 + ee146) + 3 * (ee136 *
          ee53))) * ee2/ee13 - ((((1.5 * ee85 - (3 * (ee26 * ee20) + 4.5 *
          ee101))/ee5 + ee178 + ee179 + ee180 - ee175)/ee3 + 3) *
          ee2/ee3 - 1.5)/ee22) * ee7/ee6) * ee2/ee82;
    
  }
  
  }

}

return out;

}

