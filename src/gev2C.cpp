// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

double xi_from_zero(double xi, double eps) 
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

const double xieps = 0.0001;

double xi2txi(double xi) 
{
return -log(1.5 / (1.0 + xi) - 1.0);
}

// //' Generalized extreme value (GEV) distribution negative log-likelihood with constrained shape parameter
// //'
// //' @param pars a list of vectors of coefficients for each GEV parameter
// //' @param X1 a design matrix for the GEV location parameter
// //' @param X2 a design matrix for the GEV log scale parameter
// //' @param X3 a design matrix for the GEV transformed shape parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return gevd0 a scalar, the negative log-likelihood
// //' @return gevd12 a matrix, first then second derivatives w.r.t. GEV parameters
// //' @return gevd34 a matrix, third then fourth derivatives w.r.t. GEV parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double gev2d0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate, arma::uvec nhere)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = nhere.size();

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    txivec = txivec.elem(dupid);
}

double y, mu, lpsi, txi, xi;
double ee1, ee2;
double nllh = 0.0;

for (int j=0; j < nobs; j++) {

mu = muvec[j];
lpsi = lpsivec[j];
txi = txivec[j];
xi = 1.5 / (1.0 + exp(-txi)) - 1.0;

xi = xi_from_zero(xi, xieps);  
txi = xi2txi(xi);

for (int l=0; l < nhere[j]; l++) {

y = ymat(j, l);

ee1 = xi * (y - mu) / exp(lpsi);

if (ee1 <= -1.0) {
    nllh = 1e20;
    break;
} else {

ee2 = 1.0 / xi;

nllh += lpsi + (ee2 + 1.0) * log1p(ee1) + R_pow(1.0 + ee1, -ee2);

}

} 

}

return(nllh);

}

// //' @rdname gev2d0
// [[Rcpp::export]]
arma::mat gev2d12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate, arma::uvec nhere)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = nhere.size();

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    txivec = txivec.elem(dupid);
}

double y, mu, lpsi, txi, xi;

arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);

double ee2, ee3, ee4, ee5, ee6, ee7, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
double ee20, ee21, ee22, ee23, ee25, ee27, ee28, ee29;
double ee30, ee31, ee33, ee34;

for (int j=0; j < nobs; j++) {

mu = muvec[j];
lpsi = lpsivec[j];
txi = txivec[j];
xi = 1.5 / (1.0 + exp(-txi)) - 1.0;

xi = xi_from_zero(xi, xieps);  
txi = xi2txi(xi);

for (int l=0; l < nhere[j]; l++) {

y = ymat(j, l);

ee2 = exp(-txi);
ee3 = 1 + ee2;
ee4 = 1.5/ee3;
ee5 = ee4 - 1;
ee6 = exp(lpsi);
ee7 = y - mu;
ee9 = ee5 * ee7/ee6;
ee10 = ee9 + 1;
ee11 = 1/ee5;
ee12 = R_pow(ee10, ee11);
ee13 = 1 + ee11;
ee14 = ee10 * ee6;
ee15 = R_pow(ee3, 2);
ee16 = log1p(ee9);
ee17 = 1/ee12;
ee18 = R_pow(ee10, ee13);
ee20 = ee3 * ee5;
ee21 = 1 + ee17;
ee22 = 1.5 * (ee16/(ee12 * ee5));
ee23 = 1.5/ee12;
ee25 = (((ee23 - 1.5 * ee5) * ee7/ee14 + 1.5) * ee13 - (1.5 +  ee22)/ee5)/ee10 * ee2;
ee27 = ((4.5/ee20 - 3) * ee2/ee3 + 1.5) * ee16;
ee28 = ee13 * ee5;
ee29 = ee13 * ee7;
ee30 = ee15 * ee6;
ee31 = R_pow(ee5, 2);
ee33 = 1.5 * (ee7/(ee18 * ee6));
ee34 = ee4 - ee21;

out(j, 0) = -((ee28 - ee17)/ee14);
out(j, 1) = (ee17 - ee28) * ee7/ee14 + 1;
out(j, 2) = (((ee23 - 1.5) * ee16/ee5 - ee33)/ee5 + 1.5 * (ee29/
  ee14)) * ee2/ee15;
out(j, 3) = -(ee13 * ee34 * ee5/(R_pow(ee10, 2) * R_pow(ee6, 2)));
out(j, 4) = (((ee21 - ee4) * ee7/ee14 + 1) * ee13 * ee5 - ee17)/
  ee10/ee6;
out(j, 5) = -(ee25/ee30);
out(j, 6) = -(((ee34 * ee7/ee14 - 1) * ee13 * ee5 + ee17)/ee10 *
   ee7/ee6);
out(j, 7) = -(ee25 * ee7/ee30);
out(j, 8) = (((((2.25/ee20 - 3) * ee2/ee3 + 1.5)/ee18 - 1.5 *
   ((1.5 * (ee16/(ee18 * ee31)) - 1.5 * (ee29/(R_pow(ee10, (ee11 +
   2)) * ee6))) * ee2/ee15)) * ee7/ee6 + (ee27 + (1.5 * ((ee22 -
   ee33) * ee16/ee5) - 2.25 * (ee7/ee14)) * ee2/ee15 +
   (2.25 * (ee2 * ee7/(ee10 * ee15 * ee6)) - ee27)/ee12)/ee5)/
  ee5 - (((2.25 * (ee7/(ee10 * ee3 * ee6)) - 3) * ee2/ee3 + 1.5) *
   ee13 + 2.25 * (ee2/(ee15 * ee31))) * ee7/ee14) * ee2/ee15;

}

}

return out;

}

// //' @rdname gev2d0
// [[Rcpp::export]]
arma::mat gev2d34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate, arma::uvec nhere)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = nhere.size();

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    txivec = txivec.elem(dupid);
}

double y, mu, lpsi, txi, xi;

arma::mat out = arma::mat(nobs, 25, arma::fill::zeros);

double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee20, ee21, ee22, ee23, ee24, ee25, ee28, ee29;
double ee30, ee31, ee32, ee33, ee35, ee36, ee37, ee38, ee39;
double ee45, ee46, ee47;
double ee50, ee53, ee54, ee55, ee56;
double ee61, ee64, ee65, ee66, ee67, ee68;
double ee70, ee71, ee72, ee73, ee74, ee76, ee79;
double ee80, ee81, ee82, ee83, ee84, ee87, ee88, ee89;
double ee92, ee97, ee98;
double ee100, ee101, ee102, ee103, ee104, ee105, ee106, ee108;
double ee110, ee111, ee112, ee113, ee116, ee118;
double ee120, ee123, ee125, ee127, ee129;
double ee130, ee131, ee132, ee133, ee134, ee135, ee136, ee137, ee138, ee139;
double ee140, ee142, ee143, ee146, ee149;
double ee150, ee151, ee153, ee156, ee159;
double ee161, ee164, ee167, ee168, ee169;
double ee170, ee171, ee172, ee174, ee177;
double ee180, ee184, ee186, ee188, ee189;
double ee190, ee195, ee196, ee197, ee198, ee199;
double ee200, ee201, ee202, ee203, ee204, ee205, ee206, ee207, ee208;
double ee210, ee217, ee218, ee219;
double ee221, ee223, ee226, ee227;
double ee237;
double ee241, ee242, ee243, ee244, ee246, ee248, ee249;
double ee251, ee252, ee253, ee255, ee256, ee259;
double ee260, ee261, ee263;

for (int j=0; j < nobs; j++) {

mu = muvec[j];
lpsi = lpsivec[j];
txi = txivec[j];
xi = 1.5 / (1.0 + exp(-txi)) - 1.0;

xi = xi_from_zero(xi, xieps);  
txi = xi2txi(xi);

for (int l=0; l < nhere[j]; l++) {

y = ymat(j, l);

ee2 = exp(-txi);
ee3 = 1 + ee2;
ee5 = 1.5/ee3 - 1;
ee6 = exp(lpsi);
ee7 = y - mu;
ee8 = ee5 * ee7;
ee9 = ee8/ee6;
ee10 = ee9 + 1;
ee11 = 1/ee5;
ee12 = R_pow(ee3, 2);
ee13 = 1 + ee11;
ee14 = ee11 + 2;
ee15 = log1p(ee9);
ee16 = 3 * (ee2/ee3);
ee17 = 1.5 - ee16;
ee18 = R_pow(ee5, 2);
ee19 = R_pow(ee10, ee14);
ee20 = ee11 + 3;
ee21 = R_pow(ee10, ee20);
ee22 = ee10 * ee6;
ee23 = ee3 * ee5;
ee24 = ee2/ee12;
ee25 = R_pow(ee10, ee13);
ee28 = (4.5/ee23 - 3) * ee2/ee3 + 1.5;
ee29 = 1 + 2 * ee2;
ee30 = ee17 * ee5;
ee31 = ee2 * ee7;
ee32 = 3 * ee29;
ee33 = 3 * ee3;
ee35 = ee10 * ee12 * ee6;
ee36 = 12 * ee2;
ee37 = 2.25 * ee24;
ee38 = ee21 * ee6;
ee39 = ee31/ee35;
ee45 = 1.5 * (ee15/(ee19 * ee18)) - 1.5 * (ee14 * ee7/ee38);
ee46 = ee37 - ee30;
ee47 = ee32 + ee33;
ee50 = 2.25 * ee39 - ee28 * ee15;
ee53 = ee12 * ee18;
ee54 = ee7/(ee10 * ee3 * ee6);
ee55 = R_pow(ee10, ee11);
ee56 = ee12 * ee6;
ee61 = ee8/ee22;
ee64 = 1.5 * (ee15/(ee25 * ee18)) - 1.5 * (ee13 * ee7/(ee19 *  ee6));
ee65 = 1.5/ee5;
ee66 = 2.25 * (ee2/ee53);
ee67 = 3 * ee17;
ee68 = 3 * ee46;
ee70 = (2.25 * ee54 - 3) * ee2/ee3;
ee71 = 3 * ee30;
ee72 = ee13 * ee17;
ee73 = (ee68 - (27 * ee24 + ee71))/ee5;
ee74 = 1.5 * ee13;
ee76 = R_pow(ee10, (ee11 + 4)) * ee6;
ee79 = 1.5 - ((ee47 - ee36)/ee3 + 3) * ee2/ee3;
ee80 = ee72 + ee66;
ee81 = ee64 * ee2;
ee82 = ee65 - ee74;
ee83 = ee13 * ee45;
ee84 = ee31/ee56;
ee87 = (((ee73 - ee67)/ee5 + ee32 + ee33 - ee36)/ee3 + 3) *  ee2/ee3 - 1.5;
ee88 = ee70 + 1.5;
ee89 = ee80/ee19;
ee92 = (ee50/ee25 + 1.5 * (ee81 * ee15/ee12))/ee18;
ee97 = 1.5 * (ee15/(ee21 * ee18)) - 1.5 * (ee20 * ee7/ee76);
ee98 = 1.5 * (ee15/(ee55 * ee5));
ee100 = (ee89 - 1.5 * (ee83 * ee2/ee12)) * ee7/ee6;
ee101 = ee12 * ee5;
ee102 = ee45 * ee2;
ee103 = 4 * ee2;
ee104 = 4.5 * ee17;
ee105 = ee87 * ee15;
ee106 = ee10 * ee17;
ee108 = (4.5 * ee54 - 3) * ee2/ee3;
ee110 = 1.5 * ee88 + 3 * ee28;
ee111 = 1.5 * ee28;
ee112 = ee105 + ee110 * ee2 * ee7/ee35;
ee113 = ee100 + ee92;
ee116 = ((ee17 * ee14 + ee66)/ee21 - 1.5 * (ee97 * ee14 *  ee2/ee12)) * ee7/ee6 + (ee50/ee19 + 1.5 * (ee102 * ee15/ee12))/ee18;
ee118 = (((ee104 + 6.75 * ee39) * ee7/ee22 + ee36 - ee47)/ee3 -  3) * ee2/ee3;
ee120 = ee79 * ee5;
ee123 = (2.25/ee23 - 3) * ee2/ee3 + 1.5;
ee125 = ee111 + ee67;
ee127 = ee98 - 1.5 * (ee7/(ee25 * ee6));
ee129 = 2.25 * ee84 - ee106;
ee130 = 3 * ee61;
ee131 = R_pow(ee10, 2);
ee132 = ee108 + 1.5;
ee133 = 2 * ee61;
ee134 = 2 * ee46;
ee135 = 9 * ee24;
ee136 = R_pow(ee6, 2);
ee137 = ee116 * ee13;
ee138 = ee118 + 1.5;
ee139 = (ee45 * ee5 + 1.5/ee19) * ee13;
ee140 = ((ee65 - 1.5 * ee14)/ee21 - ee97 * ee5 * ee14) *  ee13;
ee142 = ((6 * (2 * ee29 + ee103) + 8 * ee47 - 96 * ee2)/ee3 +  12) * ee2;
ee143 = ee5 * ee14;
ee146 = ee5 * ee20 * ee7/ee76;
ee149 = ee5 * (4.5 - ee130) * ee7/ee22;
ee150 = ee82 * ee14;
ee151 = ee82/ee19;
ee153 = ee70 + (ee37 - ee132 * ee5) * ee7/ee22 + 1.5;
ee156 = 1.5 - 1.5 * ee61;
ee159 = 1/ee25;
ee161 = 1/ee19 + 2/ee19;
ee164 = 2 * (1 + 2 * ee9) + 4 * ee10;
ee167 = 2 * (1.5 * ee10 + 1.5 * ee9) + 6 * ee10 - 12 * ee9;
ee168 = 3 * (1 + ee103);
ee169 = 6 * ee3;
ee170 = 8 * ee9;
ee171 = 9 * ee29;
ee172 = ee112/ee25;
ee174 = ((ee140 + ee150/ee21) * ee7/ee6 + 2 * ee83) * ee5;
ee177 = ((ee125/ee5 + ee36 - ee47)/ee3 - 3) * ee2/ee3 +  1.5;
ee180 = (ee123/ee25 - 1.5 * (ee81/ee12)) * ee7/ee6 + (ee50/ee55 +  1.5 * (ee127 * ee2 * ee15/ee101))/ee5;
ee184 = (ee5 * (2 * ee129 - 18 * ee84) + 9 * (ee10 * ee2/ee12))/ee10 +  ee135;
ee186 = ((2 * (ee80 * ee45) - 1.5 * ee137) * ee2/ee12 -  (ee13 * ee79 + ee125 * ee2/ee53)/ee19) * ee7/ee6;
ee188 = ee72 + (ee16 - 1.5)/ee5;
ee189 = ee13 * ee5;
ee190 = ee12 * ee136;
ee195 = ee167/ee10;
ee196 = (ee170 - ee164)/ee10;
ee197 = 1.5 * (ee113 * ee15);
ee198 = 2 * (ee64 * ee50);
ee199 = 2 * (ee45 * ee82 * ee2/ee12);
ee200 = 2/ee21;
ee201 = 3 * (ee88 * ee28);
ee202 = 3 * ee79;
ee203 = 4/ee19;
ee204 = 4/ee21;
ee205 = R_pow(ee6, 3);
ee206 = (((((((ee5 * (3 * ee129 - 27 * ee84) + 3 * (ee10 *  ee46))/ee10 + ee68) * ee7/ee22 + ee104) * ee2/ee12 -  ee120) * ee7/ee22 + ee118 + 1.5) * ee13 + ((((ee73 +  3 * (ee28 * ee156) + 4.5 * ee153 - ee67)/ee5 + ee32 +  ee33 - ee36)/ee3 + 3) * ee2/ee3 - 1.5)/ee5)/ee10 - ((((ee73 +  1.5 * ee17 + 4.5 * ee28)/ee25 + ee197 + ee198) * ee2/ee12 -  ee172)/ee18 + ee186)) * ee2;
ee207 = (((((ee17 * (6 * ee46 - 18 * ee24) + (12 * (ee68 -  ee71) + 9 * (ee134 + ee135) - 324 * ee24) * ee2/ee101 +  3 * (4.5 * (ee17 * ee2/ee12) - ee120) - 6 * ee120)/ee5 -  ee202)/ee5 + ee168 + ee169 + ee171 - ee142)/ee3 + 3) *  ee2/ee3 - 1.5) * ee15;
ee208 = ee112/ee55;
ee210 = (((((3 * ee156 - 4.5)/ee23 + 3) * ee2/ee3 - 1.5)/ee5 +  ee153 * ee13)/ee10 + ee100 + ee92) * ee2;
ee217 = (ee188/ee19 + ee89 + ee199 - (ee116 * ee5 + 1.5 *  (ee102/ee12)) * ee13) * ee7/ee6;
ee218 = R_pow(ee10, 3);
ee219 = R_pow(ee10, 4);
ee221 = (ee151 - ee139) * ee7/ee6;
ee223 = ((2 * (ee123 * ee64) - 1.5 * ee113) * ee2/ee12 -  ee177/ee25) * ee7/ee6;
ee226 = (ee195 + 9) * ee5 * ee7/ee22;
ee227 = ((ee134 + ee37 - (ee184 * ee7/ee22 + ee108 + 1.5) *  ee5) * ee7/ee22 + ee70 + 1.5) * ee13;
ee237 = ee5 * (ee133 - 3) * ee7/ee22 + 1;
ee241 = ee5 * (3 - ee133) * ee7/ee22 - 1;
ee242 = ee149 - 1.5;
ee243 = ee61 - 1;
ee244 = ee82/ee55;
ee246 = (ee161 - ee143 * ee7/ee38) * ee7/ee6;
ee248 = (ee161 + ee203 - ee143 * (ee200 + ee204 - ee146) *  ee7/ee6) * ee7/ee6;
ee249 = ee14/ee55;
ee251 = (2 * ee151 - (ee174 + ee139 + (ee74 - ee65)/ee19)) *  ee7/ee6;
ee252 = 1 - ee61;
ee253 = 1.5 - ee149;
ee255 = 1.5 * (ee180 * ee15) + 2 * (ee127 * ee50);
ee256 = 1.5 * ee87;
ee259 = 1/ee21;
ee260 = 2 * ee5;
ee261 = 3 - ee130;
ee263 = 4.5 * ee87 - (1.5 * ee138 + ee201);

out(j, 0) = -(ee13 * ee18 * (ee260 - ee249)/(ee218 * ee205));
out(j, 1) = (((ee249 - ee260) * ee7/ee22 + 2) * ee5 - 2/ee55)/
  ee131 * ee13 * ee5/ee136;
out(j, 2) = -(((ee189 * ee261 + ee244 - 1.5)/ee131 - ee83 *
   ee5) * ee2/ee190);
out(j, 3) = ((ee241/ee10 - ee246) * ee13 * ee5 + ee159)/ee6;
out(j, 4) = ((ee13 * ee253 - (1.5 * ee252 + ee98)/ee5)/ee10 -
   ee221) * ee2/ee56;
out(j, 5) = ee210/ee56;
out(j, 6) = -(((ee237/ee10 + ee246) * ee13 * ee5 - ee159) *
   ee7/ee6);
out(j, 7) = -(((ee242 * ee13 + (ee98 - 1.5 * ee243)/ee5)/ee10 +
   ee221) * ee2 * ee7/ee56);
out(j, 8) = ee210 * ee7/ee56;
out(j, 9) = (((ee105 + (ee255/ee5 + ee110 * ee7/ee22) * ee2/
  ee12 - ee208)/ee5 + ee223)/ee5 + (ee138 * ee13 + (ee111 + 3 *
   ee88) * ee2/ee53) * ee7/ee22) * ee2/ee12;
out(j, 10) = (ee14 * ee20/ee55 - 6 * ee5) * ee13 * R_pow(ee5, 3)/
  (ee219 * R_pow(ee6, 4));
out(j, 11) = (ee5 * (ee164 - ee170)/ee219 - (ee259 + ee200 -
   ee146) * ee14) * ee13 * ee18/ee205;
out(j, 12) = -((((ee189 * ee167 + ee150/R_pow(ee10, (ee11 -
   1)))/ee10 - 3)/ee218 + ee140) * ee5 * ee2/(ee12 * ee205));
out(j, 13) = (((ee5 * (4 - ee196) * ee7/ee22 - 4)/ee131 - (ee259 +
   ee204 - ee146) * ee14 * ee7/ee6) * ee5 + ee203) * ee13 *
   ee5/ee136;
out(j, 14) = ((ee189 * (6 - (ee195 + 6) * ee5 * ee7/ee22) +
   2 * ee244 - 1.5 * (2 - ee133))/ee131 - ee174) * ee2/ee190;
out(j, 15) = -(((ee188/ee55 + ee13 * (ee134 - ee184 * ee5 *
   ee7/ee22) + 1.5 - ((3 * ee261 - 4.5)/ee23 + 3) * ee2/ee3)/ee131 +
   ee199 - ee137 * ee5) * ee2/ee190);
out(j, 16) = ((((ee5 * (6 - ee196) * ee7/ee22 - 7) * ee5 * ee7/
  ee22 + 1)/ee10 + ee248) * ee13 * ee5 - ee159)/ee6;
out(j, 17) = (((ee5 * (10.5 - ee226) * ee7/ee22 - 1.5) * ee13 +
   (ee98 - 1.5 * ee241)/ee5)/ee10 + ee251) * ee2/ee56;
out(j, 18) = -((ee217 + (ee227 + (3 * (ee253 * ee2/ee101) -
   ee28 * ee252)/ee5)/ee10 + ee92) * ee2/ee56);
out(j, 19) = -(ee206/ee56);
out(j, 20) = -((((((ee196 - 6) * ee5 * ee7/ee22 + 7) * ee5 *
   ee7/ee22 - 1)/ee10 - ee248) * ee13 * ee5 + ee159) * ee7/ee6);
out(j, 21) = -(((((ee226 - 10.5) * ee5 * ee7/ee22 + 1.5) * ee13 -
   (1.5 * ee237 + ee98)/ee5)/ee10 - ee251) * ee2 * ee7/ee56);
out(j, 22) = -((ee217 + (ee227 - (3 * (ee242 * ee2/ee101) -
   ee243 * ee28)/ee5)/ee10 + ee92) * ee2 * ee7/ee56);
out(j, 23) = -(ee206 * ee7/ee56);
out(j, 24) = ((((((((4.5 * ee129 - (40.5 * ee84 + 9 * ee106))/
  ee10 - 9 * ee17) * ee2 * ee7/ee35 - (ee132 * ee17 + 2 * (R_pow(ee17, 2) +
   1.5 * ee79) + ee202)) * ee7/ee22 + ee168 + ee169 +
   ee171 - ee142)/ee3 + 3) * ee2/ee3 - 1.5) * ee13 + (ee256 -
   (ee201 + 4.5 * ee138)) * ee2/ee53) * ee7/ee22 + ((((1.5 *
   (((ee255 * ee2/ee101 - ee208)/ee5 + ee223) * ee15) + 3 *
   (ee180 * ee50) - 3 * (ee112 * ee127))/ee5 + ee263 * ee7/ee22) *
   ee2/ee12 - (ee207 + (ee263 * ee2 * ee7/ee35 - ee207)/ee55))/
  ee5 + ((3 * (ee113 * ee123) - 
    (1.5 * (((ee197 + ee198) *
   ee2/ee12 - ee172)/ee18 + ee186) + 3 * (ee177 * ee64))) *
   ee2/ee12 - ((((ee256 - (3 * (ee28 * ee17) + 4.5 * ee79))/
  ee5 + ee168 + ee169 + ee171 - ee142)/ee3 + 3) * ee2/ee3 - 1.5)/
  ee25) * ee7/ee6)/ee5) * ee2/ee12;
}

}

return out;

}

// //' Generalized extreme value (GEV) distribution negative log-likelihood with constrained shape parameter
// //' and sparse design matrices
// //'
// //' @param pars a list of vectors of coefficients for each GEV parameter
// //' @param X1 a sparse design matrix for the GEV location parameter
// //' @param X2 a sparse design matrix for the GEV log scale parameter
// //' @param X3 a sparse design matrix for the GEV transformed shape parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return gevspd0 a scalar, the negative log-likelihood
// //' @return gevspd12 a matrix, first then second derivatives w.r.t. GEV parameters
// //' @return gevspd34 a matrix, third then fourth derivatives w.r.t. GEV parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double gev2spd0(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::sp_mat X3, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = nhere.size();

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    txivec = txivec.elem(dupid);
}

double y, mu, lpsi, txi, xi;
double ee1, ee2;
double nllh = 0.0;

for (int j=0; j < nobs; j++) {

mu = muvec[j];
lpsi = lpsivec[j];
txi = txivec[j];
xi = 1.5 / (1.0 + exp(-txi)) - 1.0;

xi = xi_from_zero(xi, xieps);  
txi = xi2txi(xi);

for (int l=0; l < nhere[j]; l++) {

y = ymat(j, l);

ee1 = xi * (y - mu) / exp(lpsi);

if (ee1 <= -1.0) {
    nllh = 1e20;
    break;
} else {

ee2 = 1.0 / xi;

nllh += lpsi + (ee2 + 1.0) * log1p(ee1) + R_pow(1.0 + ee1, -ee2);

}

} 

}

return(nllh);

}

// //' @rdname gev2spd0
// [[Rcpp::export]]
arma::mat gev2spd12(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::sp_mat X3, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = nhere.size();

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    txivec = txivec.elem(dupid);
}

double y, mu, lpsi, txi, xi;

arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);

double ee2, ee3, ee4, ee5, ee6, ee7, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
double ee20, ee21, ee22, ee23, ee25, ee27, ee28, ee29;
double ee30, ee31, ee33, ee34;

for (int j=0; j < nobs; j++) {

mu = muvec[j];
lpsi = lpsivec[j];
txi = txivec[j];
xi = 1.5 / (1.0 + exp(-txi)) - 1.0;

xi = xi_from_zero(xi, xieps);  
txi = xi2txi(xi);

for (int l=0; l < nhere[j]; l++) {

y = ymat(j, l);

ee2 = exp(-txi);
ee3 = 1 + ee2;
ee4 = 1.5/ee3;
ee5 = ee4 - 1;
ee6 = exp(lpsi);
ee7 = y - mu;
ee9 = ee5 * ee7/ee6;
ee10 = ee9 + 1;
ee11 = 1/ee5;
ee12 = R_pow(ee10, ee11);
ee13 = 1 + ee11;
ee14 = ee10 * ee6;
ee15 = R_pow(ee3, 2);
ee16 = log1p(ee9);
ee17 = 1/ee12;
ee18 = R_pow(ee10, ee13);
ee20 = ee3 * ee5;
ee21 = 1 + ee17;
ee22 = 1.5 * (ee16/(ee12 * ee5));
ee23 = 1.5/ee12;
ee25 = (((ee23 - 1.5 * ee5) * ee7/ee14 + 1.5) * ee13 - (1.5 +  ee22)/ee5)/ee10 * ee2;
ee27 = ((4.5/ee20 - 3) * ee2/ee3 + 1.5) * ee16;
ee28 = ee13 * ee5;
ee29 = ee13 * ee7;
ee30 = ee15 * ee6;
ee31 = R_pow(ee5, 2);
ee33 = 1.5 * (ee7/(ee18 * ee6));
ee34 = ee4 - ee21;

out(j, 0) = -((ee28 - ee17)/ee14);
out(j, 1) = (ee17 - ee28) * ee7/ee14 + 1;
out(j, 2) = (((ee23 - 1.5) * ee16/ee5 - ee33)/ee5 + 1.5 * (ee29/
  ee14)) * ee2/ee15;
out(j, 3) = -(ee13 * ee34 * ee5/(R_pow(ee10, 2) * R_pow(ee6, 2)));
out(j, 4) = (((ee21 - ee4) * ee7/ee14 + 1) * ee13 * ee5 - ee17)/
  ee10/ee6;
out(j, 5) = -(ee25/ee30);
out(j, 6) = -(((ee34 * ee7/ee14 - 1) * ee13 * ee5 + ee17)/ee10 *
   ee7/ee6);
out(j, 7) = -(ee25 * ee7/ee30);
out(j, 8) = (((((2.25/ee20 - 3) * ee2/ee3 + 1.5)/ee18 - 1.5 *
   ((1.5 * (ee16/(ee18 * ee31)) - 1.5 * (ee29/(R_pow(ee10, (ee11 +
   2)) * ee6))) * ee2/ee15)) * ee7/ee6 + (ee27 + (1.5 * ((ee22 -
   ee33) * ee16/ee5) - 2.25 * (ee7/ee14)) * ee2/ee15 +
   (2.25 * (ee2 * ee7/(ee10 * ee15 * ee6)) - ee27)/ee12)/ee5)/
  ee5 - (((2.25 * (ee7/(ee10 * ee3 * ee6)) - 3) * ee2/ee3 + 1.5) *
   ee13 + 2.25 * (ee2/(ee15 * ee31))) * ee7/ee14) * ee2/ee15;
}

}

return out;

}

// //' @rdname gev2spd0
// [[Rcpp::export]]
arma::mat gev2spd34(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::sp_mat X3, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = nhere.size();

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    txivec = txivec.elem(dupid);
}

double y, mu, lpsi, txi, xi;

arma::mat out = arma::mat(nobs, 25, arma::fill::zeros);

double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee20, ee21, ee22, ee23, ee24, ee25, ee28, ee29;
double ee30, ee31, ee32, ee33, ee35, ee36, ee37, ee38, ee39;
double ee45, ee46, ee47;
double ee50, ee53, ee54, ee55, ee56;
double ee61, ee64, ee65, ee66, ee67, ee68;
double ee70, ee71, ee72, ee73, ee74, ee76, ee79;
double ee80, ee81, ee82, ee83, ee84, ee87, ee88, ee89;
double ee92, ee97, ee98;
double ee100, ee101, ee102, ee103, ee104, ee105, ee106, ee108;
double ee110, ee111, ee112, ee113, ee116, ee118;
double ee120, ee123, ee125, ee127, ee129;
double ee130, ee131, ee132, ee133, ee134, ee135, ee136, ee137, ee138, ee139;
double ee140, ee142, ee143, ee146, ee149;
double ee150, ee151, ee153, ee156, ee159;
double ee161, ee164, ee167, ee168, ee169;
double ee170, ee171, ee172, ee174, ee177;
double ee180, ee184, ee186, ee188, ee189;
double ee190, ee195, ee196, ee197, ee198, ee199;
double ee200, ee201, ee202, ee203, ee204, ee205, ee206, ee207, ee208;
double ee210, ee217, ee218, ee219;
double ee221, ee223, ee226, ee227;
double ee237;
double ee241, ee242, ee243, ee244, ee246, ee248, ee249;
double ee251, ee252, ee253, ee255, ee256, ee259;
double ee260, ee261, ee263;

for (int j=0; j < nobs; j++) {

mu = muvec[j];
lpsi = lpsivec[j];
txi = txivec[j];
xi = 1.5 / (1.0 + exp(-txi)) - 1.0;

xi = xi_from_zero(xi, xieps);  
txi = xi2txi(xi);

for (int l=0; l < nhere[j]; l++) {

y = ymat(j, l);

ee2 = exp(-txi);
ee3 = 1 + ee2;
ee5 = 1.5/ee3 - 1;
ee6 = exp(lpsi);
ee7 = y - mu;
ee8 = ee5 * ee7;
ee9 = ee8/ee6;
ee10 = ee9 + 1;
ee11 = 1/ee5;
ee12 = R_pow(ee3, 2);
ee13 = 1 + ee11;
ee14 = ee11 + 2;
ee15 = log1p(ee9);
ee16 = 3 * (ee2/ee3);
ee17 = 1.5 - ee16;
ee18 = R_pow(ee5, 2);
ee19 = R_pow(ee10, ee14);
ee20 = ee11 + 3;
ee21 = R_pow(ee10, ee20);
ee22 = ee10 * ee6;
ee23 = ee3 * ee5;
ee24 = ee2/ee12;
ee25 = R_pow(ee10, ee13);
ee28 = (4.5/ee23 - 3) * ee2/ee3 + 1.5;
ee29 = 1 + 2 * ee2;
ee30 = ee17 * ee5;
ee31 = ee2 * ee7;
ee32 = 3 * ee29;
ee33 = 3 * ee3;
ee35 = ee10 * ee12 * ee6;
ee36 = 12 * ee2;
ee37 = 2.25 * ee24;
ee38 = ee21 * ee6;
ee39 = ee31/ee35;
ee45 = 1.5 * (ee15/(ee19 * ee18)) - 1.5 * (ee14 * ee7/ee38);
ee46 = ee37 - ee30;
ee47 = ee32 + ee33;
ee50 = 2.25 * ee39 - ee28 * ee15;
ee53 = ee12 * ee18;
ee54 = ee7/(ee10 * ee3 * ee6);
ee55 = R_pow(ee10, ee11);
ee56 = ee12 * ee6;
ee61 = ee8/ee22;
ee64 = 1.5 * (ee15/(ee25 * ee18)) - 1.5 * (ee13 * ee7/(ee19 *  ee6));
ee65 = 1.5/ee5;
ee66 = 2.25 * (ee2/ee53);
ee67 = 3 * ee17;
ee68 = 3 * ee46;
ee70 = (2.25 * ee54 - 3) * ee2/ee3;
ee71 = 3 * ee30;
ee72 = ee13 * ee17;
ee73 = (ee68 - (27 * ee24 + ee71))/ee5;
ee74 = 1.5 * ee13;
ee76 = R_pow(ee10, (ee11 + 4)) * ee6;
ee79 = 1.5 - ((ee47 - ee36)/ee3 + 3) * ee2/ee3;
ee80 = ee72 + ee66;
ee81 = ee64 * ee2;
ee82 = ee65 - ee74;
ee83 = ee13 * ee45;
ee84 = ee31/ee56;
ee87 = (((ee73 - ee67)/ee5 + ee32 + ee33 - ee36)/ee3 + 3) *  ee2/ee3 - 1.5;
ee88 = ee70 + 1.5;
ee89 = ee80/ee19;
ee92 = (ee50/ee25 + 1.5 * (ee81 * ee15/ee12))/ee18;
ee97 = 1.5 * (ee15/(ee21 * ee18)) - 1.5 * (ee20 * ee7/ee76);
ee98 = 1.5 * (ee15/(ee55 * ee5));
ee100 = (ee89 - 1.5 * (ee83 * ee2/ee12)) * ee7/ee6;
ee101 = ee12 * ee5;
ee102 = ee45 * ee2;
ee103 = 4 * ee2;
ee104 = 4.5 * ee17;
ee105 = ee87 * ee15;
ee106 = ee10 * ee17;
ee108 = (4.5 * ee54 - 3) * ee2/ee3;
ee110 = 1.5 * ee88 + 3 * ee28;
ee111 = 1.5 * ee28;
ee112 = ee105 + ee110 * ee2 * ee7/ee35;
ee113 = ee100 + ee92;
ee116 = ((ee17 * ee14 + ee66)/ee21 - 1.5 * (ee97 * ee14 *  ee2/ee12)) * ee7/ee6 + (ee50/ee19 + 1.5 * (ee102 * ee15/ee12))/ee18;
ee118 = (((ee104 + 6.75 * ee39) * ee7/ee22 + ee36 - ee47)/ee3 -  3) * ee2/ee3;
ee120 = ee79 * ee5;
ee123 = (2.25/ee23 - 3) * ee2/ee3 + 1.5;
ee125 = ee111 + ee67;
ee127 = ee98 - 1.5 * (ee7/(ee25 * ee6));
ee129 = 2.25 * ee84 - ee106;
ee130 = 3 * ee61;
ee131 = R_pow(ee10, 2);
ee132 = ee108 + 1.5;
ee133 = 2 * ee61;
ee134 = 2 * ee46;
ee135 = 9 * ee24;
ee136 = R_pow(ee6, 2);
ee137 = ee116 * ee13;
ee138 = ee118 + 1.5;
ee139 = (ee45 * ee5 + 1.5/ee19) * ee13;
ee140 = ((ee65 - 1.5 * ee14)/ee21 - ee97 * ee5 * ee14) *  ee13;
ee142 = ((6 * (2 * ee29 + ee103) + 8 * ee47 - 96 * ee2)/ee3 +  12) * ee2;
ee143 = ee5 * ee14;
ee146 = ee5 * ee20 * ee7/ee76;
ee149 = ee5 * (4.5 - ee130) * ee7/ee22;
ee150 = ee82 * ee14;
ee151 = ee82/ee19;
ee153 = ee70 + (ee37 - ee132 * ee5) * ee7/ee22 + 1.5;
ee156 = 1.5 - 1.5 * ee61;
ee159 = 1/ee25;
ee161 = 1/ee19 + 2/ee19;
ee164 = 2 * (1 + 2 * ee9) + 4 * ee10;
ee167 = 2 * (1.5 * ee10 + 1.5 * ee9) + 6 * ee10 - 12 * ee9;
ee168 = 3 * (1 + ee103);
ee169 = 6 * ee3;
ee170 = 8 * ee9;
ee171 = 9 * ee29;
ee172 = ee112/ee25;
ee174 = ((ee140 + ee150/ee21) * ee7/ee6 + 2 * ee83) * ee5;
ee177 = ((ee125/ee5 + ee36 - ee47)/ee3 - 3) * ee2/ee3 +  1.5;
ee180 = (ee123/ee25 - 1.5 * (ee81/ee12)) * ee7/ee6 + (ee50/ee55 +  1.5 * (ee127 * ee2 * ee15/ee101))/ee5;
ee184 = (ee5 * (2 * ee129 - 18 * ee84) + 9 * (ee10 * ee2/ee12))/ee10 +  ee135;
ee186 = ((2 * (ee80 * ee45) - 1.5 * ee137) * ee2/ee12 -  (ee13 * ee79 + ee125 * ee2/ee53)/ee19) * ee7/ee6;
ee188 = ee72 + (ee16 - 1.5)/ee5;
ee189 = ee13 * ee5;
ee190 = ee12 * ee136;
ee195 = ee167/ee10;
ee196 = (ee170 - ee164)/ee10;
ee197 = 1.5 * (ee113 * ee15);
ee198 = 2 * (ee64 * ee50);
ee199 = 2 * (ee45 * ee82 * ee2/ee12);
ee200 = 2/ee21;
ee201 = 3 * (ee88 * ee28);
ee202 = 3 * ee79;
ee203 = 4/ee19;
ee204 = 4/ee21;
ee205 = R_pow(ee6, 3);
ee206 = (((((((ee5 * (3 * ee129 - 27 * ee84) + 3 * (ee10 *  ee46))/ee10 + ee68) * ee7/ee22 + ee104) * ee2/ee12 -  ee120) * ee7/ee22 + ee118 + 1.5) * ee13 + ((((ee73 +  3 * (ee28 * ee156) + 4.5 * ee153 - ee67)/ee5 + ee32 +  ee33 - ee36)/ee3 + 3) * ee2/ee3 - 1.5)/ee5)/ee10 - ((((ee73 +  1.5 * ee17 + 4.5 * ee28)/ee25 + ee197 + ee198) * ee2/ee12 -  ee172)/ee18 + ee186)) * ee2;
ee207 = (((((ee17 * (6 * ee46 - 18 * ee24) + (12 * (ee68 -  ee71) + 9 * (ee134 + ee135) - 324 * ee24) * ee2/ee101 +  3 * (4.5 * (ee17 * ee2/ee12) - ee120) - 6 * ee120)/ee5 -  ee202)/ee5 + ee168 + ee169 + ee171 - ee142)/ee3 + 3) *  ee2/ee3 - 1.5) * ee15;
ee208 = ee112/ee55;
ee210 = (((((3 * ee156 - 4.5)/ee23 + 3) * ee2/ee3 - 1.5)/ee5 +  ee153 * ee13)/ee10 + ee100 + ee92) * ee2;
ee217 = (ee188/ee19 + ee89 + ee199 - (ee116 * ee5 + 1.5 *  (ee102/ee12)) * ee13) * ee7/ee6;
ee218 = R_pow(ee10, 3);
ee219 = R_pow(ee10, 4);
ee221 = (ee151 - ee139) * ee7/ee6;
ee223 = ((2 * (ee123 * ee64) - 1.5 * ee113) * ee2/ee12 -  ee177/ee25) * ee7/ee6;
ee226 = (ee195 + 9) * ee5 * ee7/ee22;
ee227 = ((ee134 + ee37 - (ee184 * ee7/ee22 + ee108 + 1.5) *  ee5) * ee7/ee22 + ee70 + 1.5) * ee13;
ee237 = ee5 * (ee133 - 3) * ee7/ee22 + 1;
ee241 = ee5 * (3 - ee133) * ee7/ee22 - 1;
ee242 = ee149 - 1.5;
ee243 = ee61 - 1;
ee244 = ee82/ee55;
ee246 = (ee161 - ee143 * ee7/ee38) * ee7/ee6;
ee248 = (ee161 + ee203 - ee143 * (ee200 + ee204 - ee146) *  ee7/ee6) * ee7/ee6;
ee249 = ee14/ee55;
ee251 = (2 * ee151 - (ee174 + ee139 + (ee74 - ee65)/ee19)) *  ee7/ee6;
ee252 = 1 - ee61;
ee253 = 1.5 - ee149;
ee255 = 1.5 * (ee180 * ee15) + 2 * (ee127 * ee50);
ee256 = 1.5 * ee87;
ee259 = 1/ee21;
ee260 = 2 * ee5;
ee261 = 3 - ee130;
ee263 = 4.5 * ee87 - (1.5 * ee138 + ee201);

out(j, 0) = -(ee13 * ee18 * (ee260 - ee249)/(ee218 * ee205));
out(j, 1) = (((ee249 - ee260) * ee7/ee22 + 2) * ee5 - 2/ee55)/
  ee131 * ee13 * ee5/ee136;
out(j, 2) = -(((ee189 * ee261 + ee244 - 1.5)/ee131 - ee83 *
   ee5) * ee2/ee190);
out(j, 3) = ((ee241/ee10 - ee246) * ee13 * ee5 + ee159)/ee6;
out(j, 4) = ((ee13 * ee253 - (1.5 * ee252 + ee98)/ee5)/ee10 -
   ee221) * ee2/ee56;
out(j, 5) = ee210/ee56;
out(j, 6) = -(((ee237/ee10 + ee246) * ee13 * ee5 - ee159) *
   ee7/ee6);
out(j, 7) = -(((ee242 * ee13 + (ee98 - 1.5 * ee243)/ee5)/ee10 +
   ee221) * ee2 * ee7/ee56);
out(j, 8) = ee210 * ee7/ee56;
out(j, 9) = (((ee105 + (ee255/ee5 + ee110 * ee7/ee22) * ee2/
  ee12 - ee208)/ee5 + ee223)/ee5 + (ee138 * ee13 + (ee111 + 3 *
   ee88) * ee2/ee53) * ee7/ee22) * ee2/ee12;
out(j, 10) = (ee14 * ee20/ee55 - 6 * ee5) * ee13 * R_pow(ee5, 3)/
  (ee219 * R_pow(ee6, 4));
out(j, 11) = (ee5 * (ee164 - ee170)/ee219 - (ee259 + ee200 -
   ee146) * ee14) * ee13 * ee18/ee205;
out(j, 12) = -((((ee189 * ee167 + ee150/R_pow(ee10, (ee11 -
   1)))/ee10 - 3)/ee218 + ee140) * ee5 * ee2/(ee12 * ee205));
out(j, 13) = (((ee5 * (4 - ee196) * ee7/ee22 - 4)/ee131 - (ee259 +
   ee204 - ee146) * ee14 * ee7/ee6) * ee5 + ee203) * ee13 *
   ee5/ee136;
out(j, 14) = ((ee189 * (6 - (ee195 + 6) * ee5 * ee7/ee22) +
   2 * ee244 - 1.5 * (2 - ee133))/ee131 - ee174) * ee2/ee190;
out(j, 15) = -(((ee188/ee55 + ee13 * (ee134 - ee184 * ee5 *
   ee7/ee22) + 1.5 - ((3 * ee261 - 4.5)/ee23 + 3) * ee2/ee3)/ee131 +
   ee199 - ee137 * ee5) * ee2/ee190);
out(j, 16) = ((((ee5 * (6 - ee196) * ee7/ee22 - 7) * ee5 * ee7/
  ee22 + 1)/ee10 + ee248) * ee13 * ee5 - ee159)/ee6;
out(j, 17) = (((ee5 * (10.5 - ee226) * ee7/ee22 - 1.5) * ee13 +
   (ee98 - 1.5 * ee241)/ee5)/ee10 + ee251) * ee2/ee56;
out(j, 18) = -((ee217 + (ee227 + (3 * (ee253 * ee2/ee101) -
   ee28 * ee252)/ee5)/ee10 + ee92) * ee2/ee56);
out(j, 19) = -(ee206/ee56);
out(j, 20) = -((((((ee196 - 6) * ee5 * ee7/ee22 + 7) * ee5 *
   ee7/ee22 - 1)/ee10 - ee248) * ee13 * ee5 + ee159) * ee7/ee6);
out(j, 21) = -(((((ee226 - 10.5) * ee5 * ee7/ee22 + 1.5) * ee13 -
   (1.5 * ee237 + ee98)/ee5)/ee10 - ee251) * ee2 * ee7/ee56);
out(j, 22) = -((ee217 + (ee227 - (3 * (ee242 * ee2/ee101) -
   ee243 * ee28)/ee5)/ee10 + ee92) * ee2 * ee7/ee56);
out(j, 23) = -(ee206 * ee7/ee56);
out(j, 24) = ((((((((4.5 * ee129 - (40.5 * ee84 + 9 * ee106))/
  ee10 - 9 * ee17) * ee2 * ee7/ee35 - (ee132 * ee17 + 2 * (R_pow(ee17, 2) +
   1.5 * ee79) + ee202)) * ee7/ee22 + ee168 + ee169 +
   ee171 - ee142)/ee3 + 3) * ee2/ee3 - 1.5) * ee13 + (ee256 -
   (ee201 + 4.5 * ee138)) * ee2/ee53) * ee7/ee22 + ((((1.5 *
   (((ee255 * ee2/ee101 - ee208)/ee5 + ee223) * ee15) + 3 *
   (ee180 * ee50) - 3 * (ee112 * ee127))/ee5 + ee263 * ee7/ee22) *
   ee2/ee12 - (ee207 + (ee263 * ee2 * ee7/ee35 - ee207)/ee55))/
  ee5 + ((3 * (ee113 * ee123) - 
    (1.5 * (((ee197 + ee198) *
   ee2/ee12 - ee172)/ee18 + ee186) + 3 * (ee177 * ee64))) *
   ee2/ee12 - ((((ee256 - (3 * (ee28 * ee17) + 4.5 * ee79))/
  ee5 + ee168 + ee169 + ee171 - ee142)/ee3 + 3) * ee2/ee3 - 1.5)/
  ee25) * ee7/ee6)/ee5) * ee2/ee12;
}

}

return out;

}

