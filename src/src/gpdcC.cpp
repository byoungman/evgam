// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double xieps = 0.0;
const double xieps3 = 0.0;

// //' Generalized Pareto distribution (GPD) negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each GPD parameter
// //' @param X1 a design matrix for the GEV log scale parameter
// //' @param X2 a design matrix for the GEV shape parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return gpdd0 a scalar, the negative log-likelihood
// //' @return gpdd12 a matrix, first then second derivatives w.r.t. GPD parameters
// //' @return gpdd34 a matrix, third then fourth derivatives w.r.t. GPD parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double gpdcd0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xivec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = ymat.n_rows;

if (dcate == 1) {
    lpsivec = lpsivec.elem(dupid);
    xivec = xivec.elem(dupid);
}

double yl, yh, lpsi, xi;
double nllh = 0.0;

double ee1, ee2, ee3, ee4;

for (int j=0; j < nobs; j++) {

yl = ymat(j, 0);
yh = ymat(j, 1);
lpsi = lpsivec[j];
xi = xivec[j];

if (fabs(xi) >= xieps3) {

ee1 = 1/xi;
ee2 = exp(lpsi);
ee3 = xi * yl/ee2;
if (ee3 <= -1.0) {
    nllh = 1e20;
    break;
} else {

ee4 = xi * yh/ee2;

if (ee4 <= -1.0) {
    nllh = 1e20;
    break;
} else {

nllh -= log(R_pow(1 + ee3, -ee1) - R_pow(1 + ee4, -ee1));
  
}
}
} else {

nllh -= log(exp(-yl / exp(lpsi)) - exp(-yl / exp(lpsi)));
    
}

}

return(nllh);

}

// //' @rdname gpdd0
// [[Rcpp::export]]
arma::mat gpdcd12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xivec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = ymat.n_rows;
arma::mat out = arma::mat(nobs, 5);

if (dcate == 1) {
    lpsivec = lpsivec.elem(dupid);
    xivec = xivec.elem(dupid);
}

double yl, yh, lpsi, xi;

double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee16, ee17, ee18, ee19;
double ee21, ee22, ee24, ee26, ee28, ee29;
double ee30, ee31, ee32, ee34, ee36;

for (int j=0; j < nobs; j++) {

yl = ymat(j, 0);
yh = ymat(j, 1);
lpsi = lpsivec[j];
xi = xivec[j];

if (fabs(xi) >= xieps) {

ee1 = exp(lpsi);
ee2 = 1/xi;
ee3 = xi * yh;
ee4 = xi * yl;
ee5 = ee3/ee1;
ee6 = ee4/ee1;
ee7 = 1 + ee2;
ee8 = 1 + ee5;
ee9 = 1 + ee6;
ee10 = R_pow(ee8, ee2);
ee11 = R_pow(ee9, ee2);
ee12 = R_pow(ee8, ee7);
ee13 = R_pow(ee9, ee7);
ee16 = 1/ee11 - 1/ee10;
ee17 = log1p(ee5);
ee18 = log1p(ee6);
ee19 = ee2 + 2;
ee21 = yh/ee12 - yl/ee13;
ee22 = ee16 * ee1;
ee24 = R_pow(ee8, ee19) * ee1;
ee26 = R_pow(ee9, ee19) * ee1;
ee28 = (ee18/ee11 - ee17/ee10)/xi + ee21/ee1;
ee29 = xi * ee16;
ee30 = 1/ee12;
ee31 = 1/ee13;
ee32 = R_pow(xi, 2);
ee34 = yh * ee7/ee24;
ee36 = yl * ee7/ee26;

out(j, 0) = ee21/ee22;
out(j, 1) = -(ee28/ee29);
out(j, 2) = (R_pow(ee21, 2)/ee22 + yl * (ee31 - ee4 * ee7/ee26) -
   yh * (ee30 - ee3 * ee7/ee24))/ee22;
out(j, 3) = (yh * (ee17/(ee32 * ee12) - ee34) - (ee28 * ee21/ee29 +
   yl * (ee18/(ee32 * ee13) - ee36)))/ee22;
out(j, 4) =  - ((((ee18 * (ee18/(xi * ee11) - yl/(ee13 * ee1)) -
   ee17 * (ee17/(xi * ee10) - yh/(ee12 * ee1)))/xi + (yl/(ee9 * ee1) -
   2 * (ee18/xi))/ee11 - (R_pow(ee28, 2)/ee16 + (yh/(ee8 * ee1) -
   2 * (ee17/xi))/ee10))/xi + (yl * ((ee31 - ee18/(xi * ee13))/xi +
   ee36) - yh * ((ee30 - ee17/(xi * ee12))/xi +
   ee34))/ee1)/ee29);

} else {
    
ee1 = exp(lpsi);
ee2 = yh/ee1;
ee3 = yl/ee1;
ee6 = exp(-ee2);
ee7 = exp(-ee3);
ee8 = (ee7 - ee6) * ee1;
ee9 = yh * ee6;
ee10 = yl * ee7;
ee11 = ee10 - ee9;

out(j, 0) = -(ee11/ee8);
out(j, 1) = 0;
out(j, 2) =  - ((ee10 * (ee3 - 1) - (R_pow(ee11, 2)/ee8 + ee9 * (ee2 -
   1)))/ee8);
out(j, 3) = 0;
out(j, 4) = 0;
    
}

}

return out;

}

// //' @rdname gpdd0
// [[Rcpp::export]]
arma::mat gpdcd34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate)
{
    
arma::vec lpsivec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec xivec = X2 * Rcpp::as<arma::vec>(pars[1]);
int nobs = ymat.n_rows;
arma::mat out = arma::mat(nobs, 9);

if (dcate == 1) {
    lpsivec = lpsivec.elem(dupid);
    xivec = xivec.elem(dupid);
}

double yl, yh, lpsi, xi;

double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28;
double ee30, ee32, ee34, ee37, ee38;
double ee41, ee42, ee43, ee44, ee45, ee46, ee47, ee48, ee49;
double ee51, ee52, ee53, ee57, ee59;
double ee60, ee61, ee63, ee65, ee66, ee67;
double ee72, ee73, ee74, ee75, ee76, ee77, ee78, ee79;
double ee80, ee81, ee82, ee83, ee84, ee85, ee86, ee88, ee89;
double ee91, ee92, ee93, ee96, ee97, ee98, ee99;
double ee100, ee101, ee102, ee103, ee104, ee109;
double ee110, ee112, ee114, ee115, ee118;
double ee120, ee122, ee123, ee124, ee125, ee126, ee127, ee128, ee129;
double ee130, ee131, ee132, ee134, ee135, ee137;
double ee143, ee145, ee146, ee147, ee148, ee149;
double ee150, ee151, ee153, ee155, ee157, ee159;
double ee160, ee161, ee162, ee163, ee164, ee165, ee167, ee168, ee169;
double ee170, ee173, ee176, ee179;
double ee186, ee187, ee188, ee189;
double ee191, ee192, ee193, ee194, ee195, ee197, ee199;
double ee200, ee201, ee202, ee203, ee204, ee205, ee206, ee207, ee209;
double ee215, ee216, ee218;
double ee220, ee221, ee222, ee223, ee224, ee225, ee226, ee227, ee228, ee229;
double ee230, ee231, ee232, ee233, ee234, ee236, ee238, ee239;

for (int j=0; j < nobs; j++) {

yl = ymat(j, 0);
yh = ymat(j, 1);
lpsi = lpsivec[j];
xi = xivec[j];

if (fabs(xi) >= xieps) {

ee1 = exp(lpsi);
ee2 = 1/xi;
ee3 = xi * yh;
ee4 = xi * yl;
ee5 = ee3/ee1;
ee6 = ee4/ee1;
ee7 = 1 + ee2;
ee8 = 1 + ee5;
ee9 = 1 + ee6;
ee10 = ee2 + 2;
ee11 = log1p(ee5);
ee12 = log1p(ee6);
ee13 = R_pow(ee8, ee7);
ee14 = R_pow(ee9, ee7);
ee15 = R_pow(xi, 2);
ee16 = R_pow(ee8, ee10);
ee17 = R_pow(ee9, ee10);
ee18 = R_pow(ee8, ee2);
ee19 = R_pow(ee9, ee2);
ee20 = ee2 + 3;
ee21 = ee8 * ee1;
ee22 = ee9 * ee1;
ee23 = ee16 * ee1;
ee24 = ee17 * ee1;
ee25 = yh/ee21;
ee26 = yl/ee22;
ee27 = ee11/xi;
ee28 = ee12/xi;
ee30 = yh * ee7/ee23;
ee32 = yl * ee7/ee24;
ee34 = yh/ee13 - yl/ee14;
ee37 = ee25 - 2 * ee27;
ee38 = ee26 - 2 * ee28;
ee41 = 1/ee19 - 1/ee18;
ee42 = R_pow(ee8, ee20);
ee43 = R_pow(ee9, ee20);
ee44 = 1/ee13;
ee45 = 1/ee14;
ee46 = ee15 * ee16;
ee47 = ee15 * ee17;
ee48 = ee42 * ee1;
ee49 = ee43 * ee1;
ee51 = (ee12/ee19 - ee11/ee18)/xi + ee34/ee1;
ee52 = ee11/(ee15 * ee13);
ee53 = ee12/(ee15 * ee14);
ee57 = ee11/(xi * ee18) - yh/(ee13 * ee1);
ee59 = ee12/(xi * ee19) - yl/(ee14 * ee1);
ee60 = ee11/(xi * ee13);
ee61 = ee12/(xi * ee14);
ee63 = yh * ee10/ee48;
ee65 = yl * ee10/ee49;
ee66 = ee52 - ee30;
ee67 = ee53 - ee32;
ee72 = ee37/ee18;
ee73 = ee38/ee19;
ee74 = ee11 * ee57;
ee75 = ee11/ee46;
ee76 = ee12 * ee59;
ee77 = ee12/ee47;
ee78 = yh * ((ee44 - ee60)/xi + ee30);
ee79 = yl * ((ee45 - ee61)/xi + ee32);
ee80 = ee75 - ee63;
ee81 = ee77 - ee65;
ee82 = 2/xi;
ee83 = ee3 * ee7;
ee84 = ee4 * ee7;
ee85 = 1/ee16;
ee86 = 1/ee17;
ee88 = (ee76 - ee74)/xi + ee73;
ee89 = (ee79 - ee78)/ee1;
ee91 = (ee88 - ee72)/xi + ee89;
ee92 = ee37/ee13;
ee93 = ee38/ee14;
ee96 = yh * (ee44 - ee83/ee23);
ee97 = yh * ee66;
ee98 = yl * (ee45 - ee84/ee24);
ee99 = yl * ee67;
ee100 = ee7 * ee80;
ee101 = ee7 * ee81;
ee102 = ee41 * ee1;
ee103 = ee97 - ee99;
ee104 = ee98 - ee96;
ee109 = 2 * ee25;
ee110 = 2 * ee26;
ee112 = yh * (1/ee46 - ee100)/ee1;
ee114 = yl * (1/ee47 - ee101)/ee1;
ee115 = ee2 + 4;
ee118 = xi * ee41;
ee120 = (ee109 - 6 * ee27)/xi + yh * (ee82 + ee25)/ee21;
ee122 = (ee110 - 6 * ee28)/xi + yl * (ee82 + ee26)/ee22;
ee123 = 2/ee16;
ee124 = 2/ee17;
ee125 = ee11 * ee66;
ee126 = ee12 * ee67;
ee127 = R_pow(ee51, 2);
ee128 = ee92 + ee125;
ee129 = ee93 + ee126;
ee130 = ee15 * ee42;
ee131 = ee15 * ee43;
ee132 = ee128/ee15;
ee134 = (ee72 + ee74/xi)/xi + ee78/ee1;
ee135 = ee129/ee15;
ee137 = (ee73 + ee76/xi)/xi + ee79/ee1;
ee143 = ee134 * ee11 + 2 * (ee57 * ee37);
ee145 = ee137 * ee12 + 2 * (ee59 * ee38);
ee146 = ee120/ee18;
ee147 = ee122/ee19;
ee148 = ee51 * ee34;
ee149 = ee127/xi;
ee150 = ee132 + ee112;
ee151 = ee135 + ee114;
ee153 = ee7 * (ee85 + xi * ee80) + ee85;
ee155 = ee7 * (ee86 + xi * ee81) + ee86;
ee157 = R_pow(ee8, ee115) * ee1;
ee159 = R_pow(ee9, ee115) * ee1;
ee160 = ee85 + ee123;
ee161 = ee86 + ee124;
ee162 = ee3 * ee10;
ee163 = ee4 * ee10;
ee164 = yh * (((((2/ee13 - ee60)/xi + ee30) * ee11 - (ee92 +  2/ee13))/xi - 2 * ee30)/xi - ee112);
ee165 = yl * (((((2/ee14 - ee61)/xi + ee32) * ee12 - (ee93 +  2/ee14))/xi - 2 * ee32)/xi - ee114);
ee167 = (ee145 - ee143)/xi + ee146;
ee168 = ee91 * ee41;
ee169 = ee51 * ee103;
ee170 = (ee165 - ee164)/ee1;
ee173 = 6/xi;
ee176 = ee11/ee130 - yh * ee20/ee157;
ee179 = ee12/ee131 - yl * ee20/ee159;
ee186 = yh * (ee52 - yh * ee153/ee1);
ee187 = yh * (ee83 * (ee160 - ee162/ee48)/ee1 - ee44);
ee188 = yl * (ee53 - yl * ee155/ee1);
ee189 = yl * (ee84 * (ee161 - ee163/ee49)/ee1 - ee45);
ee191 = (ee167 - ee147)/xi + ee170;
ee192 = ee168 + ee149;
ee193 = ee91 * ee34;
ee194 = ee148/xi;
ee195 = ee51 * ee104;
ee197 = (ee37/ee16 + ee11 * ee80)/ee15 + yh * (1/ee130 -  ee10 * ee176)/ee1;
ee199 = (ee38/ee17 + ee12 * ee81)/ee15 + yl * (1/ee131 -  ee10 * ee179)/ee1;
ee200 = ee41 * ee103;
ee201 = R_pow(ee34, 2);
ee202 = 2 * ee63;
ee203 = 2 * ee65;
ee204 = yh * ee150;
ee205 = yl * ee151;
ee206 = ee193 + ee169;
ee207 = ee150 * ee11;
ee209 = ee151 * ee12;
ee215 = ee120/ee13;
ee216 = ee122/ee14;
ee218 = (ee88 - (ee72 + 2 * (ee127/ee41)))/xi + ee89;
ee220 = ee195/xi;
ee221 = ee41 * ee104;
ee222 = ee103 * ee104;
ee223 = ee201/ee1;
ee224 = R_pow(ee104, 2);
ee225 = 2 * ee192;
ee226 = 2 * (ee148/ee118);
ee227 = ee109 + ee173;
ee228 = ee110 + ee173;
ee229 = 2/ee42;
ee230 = 2/ee43;
ee231 = 4 * ee200;
ee232 = 8 * ee194;
ee233 = 8 * ee149;
ee234 = xi * ee10;
ee236 = yh * (((2 * (ee11/(xi * ee16)) - ee123)/xi - ee202)/ee15 -  ee197 * ee7)/ee1;
ee238 = yl * (((2 * (ee12/(xi * ee17)) - ee124)/xi - ee203)/ee15 -  ee199 * ee7)/ee1;
ee239 = ee189 - ee187;

out(j, 0) = ((2 * (ee201/ee102) + 2 * ee104 + ee98 - ee96) * ee34/ee102 +
   ee189 - ee187)/ee102;
out(j, 1) = (((2 * ee103 - ee226) * ee34/ee1 - ee220)/ee41 +
   ee188 - ee186)/ee102;
out(j, 2) = (ee204 - ((ee193 + ee51 * (2 * ee97 - (ee226 + 2 * ee99)))/ee118 +
   ee205))/ee102;
out(j, 3) = -(((ee167 - ((ee218 + 2 * ee91) * ee51/ee41 + ee147))/xi +
   ee170)/ee118);
out(j, 4) = (((((2 * (ee221 - ee223) + 4 * ee221 + 8 * ee223)/ee41 +
   6 * ee104) * ee34/ee102 + 2 * ee189 - 2 * ee187) * ee34 +
   ee224 + 2 * (ee34 * ee239 + ee224))/ee102 + yl * (ee45 -
   ee84 * (ee161 + 4/ee17 - ee163 * (ee230 + 4/ee43 - ee4 * ee20/ee159)/ee1)/ee1) -
   yh * (ee44 - ee83 * (ee160 + 4/ee16 -
   ee162 * (ee229 + 4/ee42 - ee3 * ee20/ee157)/ee1)/ee1))/ee102;
out(j, 5) = ((((((2 * (ee194 + ee200) + ee231 - ee232) * ee34/ee102 -
   6 * ee220)/ee41 + ee188 - ee186) * ee34 + ee222 +
   2 * (ee222 + ee34 * (ee188 - ee186)))/ee1 - ee51 * ee239/xi)/ee41 +
   yl * (yl * (ee155 + ee124 + xi * (2 * ee101 - yl * (ee7 * (ee230 +
   ee234 * ee179) + ee10/ee43)/ee1))/ee1 - ee53) -
   yh * (yh * (ee153 + ee123 + xi * (2 * ee100 - yh * (ee7 * (ee229 +
   ee234 * ee176) + ee10/ee42)/ee1))/ee1 - ee52))/ee102;
out(j, 6) = (((2 * ((ee204 - ee205) * ee34 + R_pow(ee103, 2)) -
   ((ee51 * (ee231 - ee232) + 2 * (ee192 * ee34))/ee41 + 4 * ee169) * ee34/ee118)/ee1 -
   (ee91 * ee104 + ee51 * (2 * ee188 -
   (2 * (ee195/ee118) + 2 * ee186)))/xi)/ee41 + yl * (ee135 +
   yl * ((ee86 - 2 * (ee12/ee17))/ee15 + ee203 - ee7 * (ee77 +
   xi * ee199 - ee65))/ee1) - yh * (ee132 + yh * ((ee85 - 2 * (ee11/ee16))/ee15 +
   ee202 - ee7 * (ee75 + xi * ee197 - ee63))/ee1))/ee102;
out(j, 7) = (yh * ((ee207 + 2 * (ee66 * ee37) - ee215)/ee15 +
   ee236) - ((ee191 * ee34 + ee51 * (3 * ee204 - ((((ee225 -
   ee233) * ee34 + 2 * (ee206 * ee41))/ee41 + 2 * ee206 + 2 * ee169)/ee118 +
   3 * ee205)) + 3 * (ee91 * ee103))/ee118 + yl * ((ee209 +
   2 * (ee67 * ee38) - ee216)/ee15 + ee238)))/ee102;
out(j, 8) =  - ((((((ee145/xi - ee147)/xi + ee165/ee1) * ee12 +
   3 * (ee137 * ee38) + 3 * (ee120 * ee57) - (((ee143/xi -
   ee146)/xi + ee164/ee1) * ee11 + 3 * (ee134 * ee37) + 3 * (ee122 * ee59)))/xi +
   (((24 * ee27 - 6 * ee25)/xi - yh * ee227/ee21)/xi -
   yh * (ee227/xi + yh * (ee109 + ee82)/ee21)/ee21)/ee18 -
   ((ee218 * ee91 + ee51 * (2 * ee191 - ((ee225 + 4 * ee168 -
   ee233)/ee41 + 4 * ee91) * ee51/ee118) + 2 * (ee191 * ee51 +
   R_pow(ee91, 2)))/ee41 + (((24 * ee28 - 6 * ee26)/xi -
   yl * ee228/ee22)/xi - yl * (ee228/xi + yl * (ee110 + ee82)/ee22)/ee22)/ee19))/xi +
   (yl * (((ee216 + (ee93 + 2 * ee129 +
   6/ee14 + ee126)/xi - (ee209 + (2 * ee38 + 6) * ee67))/xi + 3 * ee114)/xi -
   ee238) - yh * (((ee215 + (ee92 + 2 * ee128 +
   6/ee13 + ee125)/xi - (ee207 + (2 * ee37 + 6) * ee66))/xi + 3 * ee112)/xi -
   ee236))/ee1)/ee118);

} else {

ee1 = exp(lpsi);
ee2 = yh/ee1;
ee3 = yl/ee1;
ee6 = exp(-ee2);
ee7 = exp(-ee3);
ee8 = yh * ee6;
ee9 = yl * ee7;
ee10 = ee7 - ee6;
ee11 = ee8 * (ee2 - 1);
ee12 = ee9 - ee8;
ee13 = ee9 * (ee3 - 1);
ee14 = ee10 * ee1;
ee15 = ee13 - ee11;
ee16 = R_pow(ee12, 2);
ee20 = yh * (1 + yh * (ee2 - 3)/ee1) * ee6;
ee22 = yl * (1 + yl * (ee3 - 3)/ee1) * ee7;
ee23 = ee10 * ee15;
ee24 = ee16/ee1;
ee25 = R_pow(ee15, 2);

out(j, 0) = -((ee22 - ((2 * ee15 + ee13 - (2 * (ee16/ee14) +
   ee11)) * ee12/ee14 + ee20))/ee14);
out(j, 1) = 0;
out(j, 2) = 0;
out(j, 3) = 0;
out(j, 4) = -((ee9 * (yl * (7 + yl * (ee3 - 6)/ee1)/ee1 - 1) -
   (((2 * ee22 - (((2 * (ee23 + ee24) + 4 * ee23 - 8 * ee24)/ee10 +
   6 * ee15) * ee12/ee14 + 2 * ee20)) * ee12 + ee25 + 2 * ((ee22 -
   ee20) * ee12 + ee25))/ee14 + ee8 * (yh * (7 + yh * (ee2 -
   6)/ee1)/ee1 - 1)))/ee14);
out(j, 5) = 0;
out(j, 6) = 0;
out(j, 7) = 0;
out(j, 8) = 0;
        
}

}

return out;

}
