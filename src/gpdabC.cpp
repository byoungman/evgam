// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// //' Generalized Pareto distribution (GPD) negative log-likelihood
// //' with parameters constrained to interval [a, b]
// //'
// //' @param pars a list of vectors of coefficients for each GPD parameter
// //' @param X1 a design matrix for the GPD unlinked scale parameter
// //' @param X2 a design matrix for the GPD unlinked shape parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return gpdd0 a scalar, the negative log-likelihood
// //' @return gpdd12 a matrix, first then second derivatives w.r.t. GPD parameters
// //' @return gpdd34 a matrix, third then fourth derivatives w.r.t. GPD parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double gpdabd0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::vec yvec, arma::vec ab, arma::uvec dupid, int dcate)
{
    
arma::vec psi0vec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xi0vec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 5);
  
if (dcate == 1) {
  psi0vec = psi0vec.elem(dupid);
  xi0vec = xi0vec.elem(dupid);
}
  
double y, psi0, xi0, psi, xi;
  
double a1 = ab[0];
double b1 = ab[2];
double a2 = ab[1];
double b2 = ab[3];
  
double ee1;
double nllh = 0.0;

for (int j=0; j < nobs; j++) {

  y = yvec[j];
  psi0 = psi0vec[j];
  xi0 = xi0vec[j];
  psi = a1 + b1 / (1 + exp(-psi0));
  xi = a2 + b2 / (1 + exp(-xi0));

  ee1 = xi * y / psi;

if (ee1 <= -1.0) {
  nllh = 1e20;
  break;
} else {
  nllh += log(psi) + (1.0 / xi + 1.0) * log1p(ee1);
}

}

return(nllh);

}

// //' @rdname gpdabd0
// [[Rcpp::export]]
arma::mat gpdabd12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::vec yvec, arma::vec ab, arma::uvec dupid, int dcate)
{
    
arma::vec psi0vec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xi0vec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 5);

if (dcate == 1) {
    psi0vec = psi0vec.elem(dupid);
    xi0vec = xi0vec.elem(dupid);
}

double y, psi0, xi0;

double a1 = ab[0];
double b1 = ab[2];
double a2 = ab[1];
double b2 = ab[3];

double ee2, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee20, ee21, ee22, ee23, ee24;

for (int j=0; j < nobs; j++) {

  y = yvec[j];
  psi0 = psi0vec[j];
  xi0 = xi0vec[j];
  
  ee2 = exp(-xi0);
  ee4 = exp(-psi0);
  ee5 = 1 + ee2;
  ee6 = 1 + ee4;
  ee7 = a2 + b2/ee5;
  ee8 = a1 + b1/ee6;
  ee9 = y * ee7;
  ee10 = ee9/ee8;
  ee11 = 1 + ee10;
  ee12 = ee11 * ee8;
  ee13 = 1/ee7;
  ee14 = R_pow(ee5, 2);
  ee15 = 1 + ee13;
  ee16 = R_pow(ee6, 2);
  ee17 = R_pow(ee7, 2);
  ee18 = ee6 * ee8;
  ee19 = ee16 * ee8;
  ee20 = b2 * ee2;
  ee21 = b2 * y;
  ee22 = log1p(ee10);
  ee23 = y * ee15;
  ee24 = ee9/ee12;
  
  out(j, 0) = b1 * (1 - ee23 * ee7/ee12) * ee4/ee19;
  out(j, 1) = ee20 * (ee23/ee12 - ee22/ee17)/ee14;
  out(j, 2) = -(b1 * ((b1/ee18 - 2) * ee4/ee6 + 1 - y * ((b1 * (2 -
    ee24)/ee18 - 2) * ee4/ee6 + 1) * ee15 * ee7/ee12) * ee4/ee19);
  out(j, 3) =  - (b1 * b2 * y * ((1 - ee24) * ee15 - ee13) * ee4 * ee2/(ee16 * ee14 * ee11 * R_pow(ee8, 2)));
  out(j, 4) = -(b2 * ((ee21 * ee2/(ee14 * ee11 * ee8) - ((2 * (b2/(ee5 * ee7)) -
    2) * ee2/ee5 + 1) * ee22)/ee17 + y * (((ee21/(ee5 * ee11 * ee8) -
    2) * ee2/ee5 + 1) * ee15 + ee20/(ee14 * ee17))/ee12) * ee2/ee14);

}

return out;

}

// //' @rdname gpdabd0
// [[Rcpp::export]]
arma::mat gpdabd34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::vec yvec, arma::vec ab, arma::uvec dupid, int dcate)
{
    
arma::vec psi0vec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xi0vec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 9);
  
if (dcate == 1) {
  psi0vec = psi0vec.elem(dupid);
  xi0vec = xi0vec.elem(dupid);
}
  
double y, psi0, xi0;

double a1 = ab[0];
double b1 = ab[2];
double a2 = ab[1];
double b2 = ab[3];

double ee2, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee19;
double ee20, ee21, ee23, ee24, ee25, ee26, ee28, ee29;
double ee30, ee31, ee32, ee35, ee36, ee37, ee38, ee39;
double ee40, ee41, ee42, ee43, ee46, ee47, ee48, ee49;
double ee50, ee51, ee52, ee55, ee56, ee57, ee58, ee59;
double ee61, ee63, ee66, ee68, ee69;
double ee70, ee72, ee73, ee74, ee76, ee78, ee79;
double ee81, ee83, ee86, ee87, ee88, ee89;
double ee90, ee91, ee92, ee93, ee94, ee96, ee97, ee99;
double ee101, ee102, ee103, ee104, ee107, ee109;
double ee112, ee114, ee115, ee119;
double ee120, ee123, ee124, ee127, ee128, ee129;
double ee131, ee132, ee133, ee136, ee137;
double ee142, ee145, ee147, ee148;
double ee151, ee152, ee153;

for (int j=0; j < nobs; j++) {

  y = yvec[j];
  psi0 = psi0vec[j];
  xi0 = xi0vec[j];
  
  ee2 = exp(-xi0);
  ee4 = exp(-psi0);
  ee5 = 1 + ee4;
  ee6 = 1 + ee2;
  ee7 = a1 + b1/ee5;
  ee8 = a2 + b2/ee6;
  ee9 = y * ee8;
  ee10 = ee9/ee7;
  ee11 = 1 + ee10;
  ee12 = R_pow(ee6, 2);
  ee13 = 1 - 2 * (ee2/ee6);
  ee14 = R_pow(ee5, 2);
  ee15 = ee11 * ee7;
  ee16 = b2 * y;
  ee17 = 1 - 2 * (ee4/ee5);
  ee19 = b2 * ee2/ee12;
  ee20 = 1 + 2 * ee2;
  ee21 = ee5 * ee7;
  ee23 = b1 * ee4/ee14;
  ee24 = ee13 * ee8;
  ee25 = 1 + 2 * ee4;
  ee26 = ee17 * ee7;
  ee28 = ee6 * ee11 * ee7;
  ee29 = 2 * ee20;
  ee30 = 2 * ee25;
  ee31 = ee16/ee28;
  ee32 = ee29 + 2 * ee6;
  ee35 = (2 * (b1/ee21) - 2) * ee4/ee5 + 1;
  ee36 = R_pow(ee7, 2);
  ee37 = ee30 + 2 * ee5;
  ee38 = 8 * ee2;
  ee39 = ee9/ee15;
  ee40 = 8 * ee4;
  ee41 = ee19 - ee24;
  ee42 = 2 * ee13;
  ee43 = 2 * ee31;
  ee46 = (ee43 - 2) * ee2/ee6 + 1;
  ee47 = 1 + 1/ee8;
  ee48 = 2 * ee41;
  ee49 = ee23 - ee26;
  ee50 = 2 * ee49;
  ee51 = b1 * y;
  ee52 = 8 * ee23;
  ee55 = (ee32 - ee38)/ee6 + 2;
  ee56 = 2 * ee26;
  ee57 = 2 * ee24;
  ee58 = 2 * ee17;
  ee59 = 8 * ee19;
  ee61 = ee16 * ee2/(ee12 * ee7);
  ee63 = ee55 * ee2/ee6;
  ee66 = ee14 * ee12 * ee11 * ee36;
  ee68 = (ee37 - ee40)/ee5 + 2;
  ee69 = (ee48 - (ee57 + ee59))/ee8;
  ee70 = ee46 + ee42;
  ee72 = (ee31 - 2) * ee2/ee6;
  ee73 = ee56 + ee52;
  ee74 = 2 * ee39;
  ee76 = b1 * b2 * y;
  ee78 = ee68 * ee4/ee5;
  ee79 = ee13 * ee11;
  ee81 = ee14 * ee11 * ee36;
  ee83 = ee69 - ee42;
  ee86 = (2 * (b2/(ee6 * ee8)) - 2) * ee2/ee6 + 1;
  ee87 = R_pow(ee8, 2);
  ee88 = ee72 + 1;
  ee89 = 1 - ee63;
  ee90 = 1 - ee39;
  ee91 = 2 - ee74;
  ee92 = 4 * ee6;
  ee93 = 4 * ee4;
  ee94 = 4 * ee2;
  ee96 = ee89 * ee8;
  ee97 = ee14 * ee7;
  ee99 = (ee37 + b1 * ((ee50 - ee73)/ee7 - ee58)/ee7 - ee40)/ee5 +  2;
  ee101 = (ee38 + ee16 * ee70/ee15 - ee32)/ee6 - 2;
  ee102 = 1 - ee78;
  ee103 = 2 * (ee61 - ee79);
  ee104 = 8 * ee61;
  ee107 = y * (ee19 - ee46 * ee8)/ee15;
  ee109 = ee99 * ee4/ee5;
  ee112 = ((ee37 + b1 * ((ee50 + y * ((b1 * ee91/ee21 - 2) *  ee4/ee5 + 1 + 2 * ee35) * ee8/ee11 - ee73)/ee7 - ee58)/ee7 -  ee40)/ee5 + 2) * ee4/ee5 - 1;
  ee114 = ee35 * ee90 - ee51 * ee91 * ee8 * ee4/ee81;
  ee115 = ee35 * ee11;
  ee119 = ((4 * (ee30 + ee93) + 8 * ee37 - 64 * ee4)/ee5 +  8) * ee4;
  ee120 = ((4 * (ee29 + ee94) + 8 * ee32 - 64 * ee2)/ee6 +  8) * ee2;
  ee123 = ee101 * ee2/ee6 + 1;
  ee124 = ee102 * ee7;
  ee127 = ee12 * ee11 * ee7;
  ee128 = ee12 * ee8;
  ee129 = ee12 * ee87;
  ee131 = (ee32 + b2 * ee83/ee8 - ee38)/ee6 + 2;
  ee132 = (ee103 - ee104) * ee8;
  ee133 = (4 * ee20 + ee92 + b2 * (ee83/ee8 - y * ee70/ee15) -  16 * ee2)/ee6;
  ee136 = (b1 * (2 - ee39)/ee21 - 2) * ee4/ee5 + 1;
  ee137 = ee88 + ee107;
  ee142 = 2 * (1 + ee93) + 4 * ee5 + 6 * ee25;
  ee145 = 2 * (1 + ee94) + ee92 + 6 * ee20;
  ee147 = 3 * (ee86 * ee88) + 4;
  ee148 = 4 * ee19;
  ee151 = ee51 * ee8 * ee4/(ee14 * ee36);
  ee152 = log1p(ee10);
  ee153 = y/ee15;
  
  // third derivatives
  // 1=scale0, 2=shape0
  // order: 111, 112, 122, 222

  out(j, 0) = b1 * (((ee40 + b1 * (ee35 + ee58)/ee7 - ee37)/ee5 -
    2) * ee4/ee5 + 1 + y * ee112 * ee47 * ee8/ee15) * ee4/ee97;
  out(j, 1) = ee76 * (ee114 * ee47 - ee136/ee8) * ee4 * ee2/ee66;
  out(j, 2) = ee76 * (((2 - ee43) * ee2/ee6 - 1)/ee8 + ee137 * ee47) * ee4 * ee2/ee66;
  out(j, 3) = b2 * (((ee131 * ee2/ee6 - 1) * ee152 + ee16 * ((b2 * (4/ee8 +
    ee153)/ee6 - 6) * ee2/ee6 + 3) * ee2/ee127)/ee87 +
    y * (ee123 * ee47 + b2 * ((b2 * (2 * ee153 + 2/ee8)/ee6 -
    6) * ee2/ee6 + 3) * ee2/ee129)/ee15) * ee2/ee12;

  // fourth derivatives
  // 1=scale0, 2=shape0
  // order: 1111, 1112, 1122, 1222, 2222

  out(j, 4) = b1 * (((ee142 + b1 * ((2 * ee68 + b1 * ((ee50 -
    (4 * ee26 + ee52))/ee7 - 4 * ee17)/ee21) * ee4/ee5 - (ee35 * ee17 +
    2 + 2 * (R_pow(ee17, 2) + 1 - ee78)))/ee7 - ee119)/ee5 +
    2) * ee4/ee5 + y * (((ee119 + b1 * ((4 * ee124 + y * (ee35 * ((2 +
    b1 * (ee74 - 2)/ee21) * ee4/ee5 - 1) + (2 * ee99 +
    ee51 * ((2 * (ee115 + ee151) + 4 * ee115 - 8 * ee151)/ee11 +
    4 * ee35) * ee8/(ee5 * ee11 * ee36)) * ee4/ee5 + 2 * (ee109 -
    (R_pow(ee35, 2) + 1)) - 2) * ee8/ee11 - (ee17 * (6 * ee49 -
    ee52) + 2 * (3 * (b1 * ee17 * ee4/ee14) - ee124) + b1 * (4 * (ee50 +
    4 * ee23) + 8 * (ee50 - ee56) - 64 * ee23) * ee4/ee97))/ee7 +
    2 * ee102)/ee7 - ee142)/ee5 - 2) * ee4/ee5 +
    1) * ee47 * ee8/ee15 - 1) * ee4/ee97;
  out(j, 5) = ee76 * (((ee109 - 1) * ee90 + ee51 * (ee35 * (6 -
    6 * ee39) - ee51 * (2 * (1 + 2 * ee10) + 4 * ee11 - 8 * ee10) * ee8 * ee4/(ee14 * R_pow(ee11, 2) * ee36)) * ee8 * ee4/ee81) * ee47 -
    ee112/ee8) * ee4 * ee2/ee66;
  out(j, 6) = -(ee76 * ((ee35 * ee137 + ee51 * (ee48 - y * ((ee132 +
    4 * (b2 * ee11 * ee2/ee12))/ee11 + ee148) * ee8/ee15) * ee4/ee81) * ee47 +
    (2 * (b2 * ee114 * ee2/ee128) - ee86 * ee136)/ee8) * ee4 * ee2/ee66);
  out(j, 7) = -(ee76 * ((((ee32 + b2 * (ee69 + 3 + 3 * (ee86 * ee90) +
    3 * ee72 + 3 * ee107 - ee42)/ee8 - ee38)/ee6 + 2) * ee2/ee6 -
    1)/ee8 + (ee123 + y * (b2 * (3 * ee13 + y * ((ee132 +
    2 * (ee11 * ee41))/ee11 + ee48)/ee15) * ee2/ee12 - ee96)/ee15) * ee47) * ee4 * ee2/ee66);
  out(j, 8) = b2 * ((ee16 * ((ee133 + 2 * ee131 + 4) * ee2/ee6 -
    ee147) * ee2/ee127 - (((ee145 + b2 * ((ee13 * (6 * ee41 -
    ee59) + 2 * (3 * (b2 * ee13 * ee2/ee12) - ee96) + b2 * (4 * (ee48 +
    ee148) + 8 * (ee48 - ee57) - 64 * ee19) * ee2/ee128 -
    4 * ee96)/ee8 - 2 * ee89)/ee8 - ee120)/ee6 + 2) * ee2/ee6 -
    1) * ee152)/ee87 + y * ((((ee145 + ee16 * ((2 * ee55 + ee16 * ((ee103 -
    (4 * ee79 + ee104))/ee11 - 4 * ee13)/ee28) * ee2/ee6 -
    (ee46 * ee13 + 2 + 2 * (R_pow(ee13, 2) + 1 - ee63)))/ee15 -
    ee120)/ee6 + 2) * ee2/ee6 - 1) * ee47 + b2 * ((ee133 +
    4 - 2 * ee101) * ee2/ee6 - ee147) * ee2/ee129)/ee15) * ee2/ee12;
  
}

return out;

}
