// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double xieps = 0.0;

// //' Generalized extreme value (GEV) distribution negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each GEV parameter
// //' @param X1 a design matrix for the GEV location parameter
// //' @param X2 a design matrix for the GEV log scale parameter
// //' @param X3 a design matrix for the GEV shape parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return gevd0 a scalar, the negative log-liklihood
// //' @return gevd12 a matrix, first then second derivatives w.r.t. GEV parameters
// //' @return gevd34 a matrix, third then fourth derivatives w.r.t. GEV parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double gevd0(const Rcpp::List& pars, const arma::mat& X1, const arma::mat& X2, const arma::mat& X3, arma::vec yvec, const arma::uvec& dupid, int dcate)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    xivec = xivec.elem(dupid);
}

double y, mu, lpsi, xi;
double ee1, ee2;
double nllh=0.0;

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

return(nllh);

}

// //' @rdname gevd0
// [[Rcpp::export]]
arma::mat gevd12(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec, const arma::uvec dupid, int dcate)
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

double y, mu, lpsi, xi;
double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9, ee10; 
double ee11, ee12, ee13, ee16, ee17, ee18, ee19, ee20; 
double ee22, ee23; 

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
out(j, 3) = -(ee18 * ee23/(ee5 * ee5 * ee1 * ee1));
out(j, 4) = (xi * ee16 * ee7 - ee10)/ee5/ee1;
out(j, 5) = -(ee19/ee1);
out(j, 6) = -((ee10 + xi * (ee23 * ee2/ee9 - 1) * ee7)/ee5 * ee2/ee1);
out(j, 7) = -(ee19 * ee2/ee1);
out(j, 8) = ((((ee2/ee9 - 2 * (ee11/xi))/R_pow(ee5, ee6 - 1) - 
ee2/ee1)/ee5 + (2 + ee17 - ee22) * ee11/xi)/xi + 
(ee13/(R_pow(ee5, ee6 + 2) * ee1) + (1/ee12 - ee11/(xi * 
ee12))/xi) * ee2/ee1)/xi - (ee20 + 1/(xi * xi)) * ee2/ee9;

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
out(j, 3) = ee5/(ee1 * ee1);
out(j, 4) = ee7/ee1;
out(j, 5) = 0;
out(j, 6) = ee7 * ee2/ee1;
out(j, 7) = 0;
out(j, 8) = 0;

}

}

return out;

}

// //' @rdname gevd0
// [[Rcpp::export]]
arma::mat gevd34(const Rcpp::List& pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec, const arma::uvec dupid, int dcate)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 25);

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    xivec = xivec.elem(dupid);
}

double y, mu, lpsi, xi;

double ee1, ee2, ee3, ee4, ee5, ee6, ee8, ee9; 
double ee10, ee11, ee12, ee13, ee14, ee16, ee17, ee18, ee19; 
double ee20, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29; 
double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee37, ee39; 
double ee40, ee41, ee42, ee43, ee44, ee45, ee46, ee48, ee49; 
double ee50, ee51, ee52, ee54, ee55, ee56, ee57, ee59; 
double ee61, ee62, ee63, ee64, ee65, ee67, ee68, ee69; 
double ee70, ee71, ee72, ee73, ee74, ee75, ee76, ee78; 
double ee80, ee81, ee83, ee84, ee85, ee86, ee87, ee88; 
double ee90, ee91, ee92, ee93, ee94, ee95, ee96, ee97, ee98; 
double ee100, ee101, ee102, ee103, ee104, ee106, ee107, ee109; 
double ee110, ee111, ee112, ee114, ee116, ee117, ee118; 
double ee120, ee124, ee125, ee126, ee127, ee128; 
double ee130, ee134, ee135, ee138, ee139; 
double ee142, ee143, ee144, ee145, ee146, ee147; 
double ee150, ee151, ee152, ee154, ee155, ee158, ee159; 
double ee160, ee161, ee162, ee166, ee167, ee168, ee169; 
double ee170, ee172, ee173, ee174, ee176, ee177, ee178, ee179; 
double ee180, ee181, ee182, ee183, ee184, ee187, ee189; 
double ee192, ee194, ee195, ee196, ee197, ee198, ee199; 
double ee202, ee204, ee205, ee206, ee207, ee208, ee209; 
double ee210, ee211, ee212, ee213, ee214, ee215, ee218; 
double ee220, ee221, ee224, ee225, ee226, ee228, ee229; 
double ee230, ee231, ee232, ee233, ee234, ee235, ee236, ee237, ee238, ee239; 
double ee240, ee241, ee246, ee247, ee248, ee249; 
double ee250, ee251, ee252, ee253, ee254, ee255, ee256, ee257, ee258, ee259; 
double ee261, ee264, ee265, ee266, ee267, ee268, ee269; 
double ee270, ee271, ee272, ee273, ee274, ee275, ee276, ee277, ee278, ee279; 
double ee280, ee281, ee282, ee283, ee284, ee285, ee286, ee287; 
double ee292, ee293, ee294, ee296, ee298, ee299; 
double ee303, ee304, ee305, ee306, ee307, ee308, ee309; 
double ee310, ee319; 
double ee320, ee321, ee322, ee323, ee324, ee325, ee328, ee329; 
double ee330, ee332, ee333, ee334, ee335, ee337, ee339; 
double ee340, ee341, ee342, ee343, ee344, ee345, ee346, ee347, ee348, ee349; 
double ee350, ee351, ee352, ee353, ee354, ee355, ee356, ee357, ee358, ee359; 
double ee360, ee361, ee362, ee363, ee364, ee365, ee366, ee367, ee368, ee369; 
double ee370, ee371, ee372, ee373, ee374, ee375, ee377, ee378, ee379;

for (int j=0; j < nobs; j++) {

y = yvec[j];
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];

if (fabs(xi) >= xieps) {

ee1 = exp(lpsi);
ee2 = (y - mu);
ee3 = xi * ee2;
ee4 = ee3/ee1;
ee5 = (1 + ee4);
ee8 = (-1/xi);
ee9 = ee1 * ee1;
ee11 = (ee8 - 1);
ee12 = xi * xi;
ee13 = ee3 * ee1;
ee14 = ee2/ee1;;
ee16 = (ee11 - 1);
ee17 = 1/ee12;
ee18 = ee13/ee9;
ee19 = (ee14);
ee20 = ee1 * ee1;
ee22 = 2 * (ee20);
ee23 = (ee9);
ee24 = (ee17);
ee25 = log(ee5);
ee26 = xi/ee1;
ee27 = R_pow(ee5, ee16);
ee28 = ee5 * ee5;
ee29 = (ee12);
ee30 = (ee22);
ee31 = (ee18);
ee32 = ee23 * ee23;
ee33 = 2 * xi;
ee34 = (ee25 * ee24);
ee35 = ee29 * ee29;
ee36 = ee13 * ee30;
ee37 = ee2 * ee1;
ee39 = (ee16 - 1);
ee40 = (ee26);
ee41 = ee37/ee9;
ee42 = xi * ee1;
ee43 = ee33/ee35;
ee44 = R_pow(ee5, ee11);
ee45 = R_pow(ee5, ee39);
ee46 = (ee28);
ee48 = ee18 - ee36/ee32;
ee49 = ee42/ee9;
ee50 = 1/ee1;
ee51 = ee46 * ee46;
ee52 = (ee48);
ee54 = (ee11 * ee19);
ee55 = ee14/ee5;
ee56 = (ee41);
ee57 = (ee43);
ee59 = (ee16 * ee19);
ee61 = 2 * (ee14 * ee5);
ee62 = (ee27 * ee54 + ee44 * ee34);
ee63 = (ee45 * ee59 + ee27 * ee34);
ee64 = (ee50);
ee65 = (ee61);
ee67 = (ee11 * ee31);
ee68 = (ee55 * ee24 - ee25 * ee57);
ee69 = (ee17 * ee19);
ee70 = (ee49);
ee71 = ee20 + ee20;
ee72 = (ee71);
ee73 = ee14 * ee19;
ee74 = 2 * ee72;
ee75 = (ee32);
ee76 = (ee74);
ee78 = 2 * (ee22 * ee23);
ee80 = (ee16 * ee31);
ee81 = (ee35);
ee83 = 2 * (ee33 * ee29);
ee84 = ee75 * ee75;
ee85 = (ee78);
ee86 = ee17 * ee31;
ee87 = ee81 * ee81;
ee88 = (ee83);
ee90 = 2 * (ee18 * ee5);
ee91 = ee13 * ee76;
ee92 = (1/xi + 1);
ee93 = ee37 * ee30;
ee94 = ee36 + ee91;
ee95 = ee45 * ee80;
ee96 = (ee90);
ee97 = ee17 * ee40;
ee98 = ee33 * ee88;
ee100 = (ee11 * ee40);
ee101 = (ee94);
ee102 = ee42 * ee30;
ee103 = ee26 * ee40;
ee104 = (2/ee35 - ee98/ee87);
ee106 = (ee8 * ee19);
ee107 = R_pow(ee5, ee8);
ee109 = ee41 - ee93/ee32;
ee110 = ee73/ee28;
ee111 = ee1/ee9;
ee112 = ee36 * ee85;
ee114 = (ee63 * ee54 + ee27 * ee69 + (ee62 * ee34 + ee44 * 
  ee68));
ee116 = (ee109);
ee117 = ee27 * ee67;
ee118 = ee48 - (ee101/ee32 - ee112/ee84);
ee120 = ee49 - ee102/ee32;
ee124 = (ee8 * ee40);
ee125 = (ee11 * ee52);
ee126 = (ee86 + ee11 * ee56);
ee127 = (ee118);
ee128 = (ee120);
ee130 = ee18 * ee31;
ee134 = ee51 * ee51;
ee135 = R_pow(ee5, ee39 - 1);
ee138 = (ee8 * ee70);
ee139 = (ee97 + ee8 * ee64);
ee142 = ee18 * ee19;
ee143 = ee26 * ee31;
ee144 = (ee44 * ee106 + ee107 * ee34);
ee145 = ee26 * ee19;
ee146 = ee41 * ee19;
ee147 = ee50 * ee19;
ee150 = (ee63 * ee67 + ee27 * ee126);
ee151 = (ee27 * ee125 - ee95 * ee67);
ee152 = (ee135 * (ee39 * ee19) + ee45 * ee34);
ee154 = (2 * (ee61 * ee46));
ee155 = (ee111);
ee158 = (ee8 * ee52);
ee159 = (ee110);
ee160 = (ee86 + ee8 * ee56);
ee161 = ee52 * ee31;
ee162 = ee49 * ee31;
ee166 = (ee11 * ee70);
ee167 = (ee16 * ee40);
ee168 = (ee55 * ee57 + ee110 * ee24 + (ee55 * ee57 + ee25 * 
  ee104));
ee169 = (ee97 + ee11 * ee64);
ee170 = (ee43 * ee19);
ee172 = ee146/ee28;
ee173 = ee147/ee28;
ee174 = 2 * (ee26 * ee5);
ee176 = (ee8 * ee31);
ee177 = (ee55);
ee178 = (ee50 * ee40 + ee26 * ee64);
ee179 = (ee17 * ee56 - ee43 * ee31 + ee17 * ee56);
ee180 = (ee17 * ee64 - ee43 * ee40 + ee17 * ee64);
ee181 = (ee174);
ee182 = ee73 * ee65;
ee183 = ee17 * ee52;
ee184 = ee17 * ee70;
ee187 = (ee8 * ee128);
ee189 = (ee152 * ee59 + ee45 * ee69 + (ee63 * ee34 + ee27 * 
  ee68));
ee192 = (ee41 * ee31 + ee18 * ee56);
ee194 = (ee50 * ee31 + ee26 * ee56);
ee195 = (ee184 + ee8 * ee155);
ee196 = (2 * (ee73));
ee197 = (ee26 * ee70 + ee49 * ee40);
ee198 = ee182/ee51;
ee199 = ee43 * ee159;
ee202 = (ee8 * ee127);
ee204 = (ee62 * ee106 + ee44 * ee69 + (ee144 * ee34 + ee107 * 
  ee68));
ee205 = (ee161 + ee18 * ee52);
ee206 = (ee41/ee5 - ee142/ee28);
ee207 = (ee198);
ee208 = ee27 * ee100;
ee209 = ee45 * ee167;
ee210 = (ee50/ee5 - ee145/ee28);
ee211 = (ee183 + ee8 * ee116);
ee212 = (ee26 * ee52 - ee162);
ee213 = ee142 * ee65;
ee214 = ee145 * ee65;
ee215 = ee103 * ee181;
ee218 = (ee63 * ee100 + ee27 * ee169);
ee220 = ee62 * ee69;
ee221 = (ee27 * ee166 + ee95 * ee100);
ee224 = (ee172 + (ee172 - ee213/ee51));
ee225 = ee95 * ee125;
ee226 = ee135 * (ee39 * ee31);
ee228 = (ee173 + (ee173 - ee214/ee51));
ee229 = (2 * (ee90 * ee46));
ee230 = ee104 * ee177;
ee231 = ee101 * ee85;
ee232 = ee52 * ee19;
ee233 = ee17 * ee207;
ee234 = 2 * ee88;
ee235 = ee130 * ee65;
ee236 = ee130 * ee96;
ee237 = ee49 * ee19;
ee238 = ee143 * ee65;
ee239 = ee143 * ee96;
ee240 = ee103 * ee65;
ee241 = ee103 * ee96;
ee246 = (ee189 * ee54 + ee63 * ee69 + (ee63 * ee69 - ee27 * 
  ee170) + (ee114 * ee34 + ee62 * ee68 + (ee62 * ee68 - 
  ee44 * ee168)));
ee247 = ee204 * ee68;
ee248 = ee114 * ee69;
ee249 = ee114 * ee160;
ee250 = ee114 * ee139;
ee251 = (ee152 * ee80 + ee45 * (ee86 + ee16 * ee56));
ee252 = (ee116/ee5 - ee232/ee28 + (ee192/ee28 - ee235/ee51));
ee253 = ee144 * ee168;
ee254 = ee62 * ee179;
ee255 = ee62 * ee180;
ee256 = ee62 * ee170;
ee257 = ee151 * ee158;
ee258 = ee151 * ee138;
ee259 = (ee45 * (ee16 * ee52) - ee226 * ee80);
ee261 = (ee178/ee28 - ee240/ee51);
ee264 = ee117 * ee202;
ee265 = ee117 * ee187;
ee266 = ee209 * ee100;
ee267 = (ee233 + ee199);
ee268 = (2 * (ee52 * ee5 - ee130));
ee269 = (2 * (ee41 * ee5 + ee142));
ee270 = ee230 - ee199;
ee271 = (ee111/ee5 - ee237/ee28 - (ee194/ee28 - ee238/ee51));
ee272 = ee127 * ee31;
ee273 = ee52 * ee52;
ee274 = ee128 * ee31;
ee275 = ee146 * ee65;
ee276 = ee110 * ee57;
ee277 = ee55 * ee104;
ee278 = ee147 * ee65;
ee279 = ee17 * ee224;
ee280 = ee17 * ee159;
ee281 = ee17 * ee228;
ee282 = ee43 * ee56;
ee283 = ee43 * ee206;
ee284 = ee43 * ee177;
ee285 = ee43 * ee64;
ee286 = ee43 * ee210;
ee287 = ee49 * ee52;
ee292 = (ee150 * ee158 + ee117 * ee211);
ee293 = (ee150 * ee138 + ee117 * ee195);
ee294 = (ee189 * ee67 + ee63 * ee126 + (ee63 * ee126 + ee27 * 
  ee179));
ee296 = ee247 - ee253;
ee298 = (ee257 + ee264);
ee299 = (ee258 + ee265);
ee303 = ee114 * ee106 + ee220 + (ee220 - ee44 * ee170) + 
  (ee204 * ee34 + ee144 * ee68 + (ee144 * ee68 - ee107 * 
  ee168));
ee304 = ee248 - ee256;
ee305 = ee249 + ee254;
ee306 = ee250 + ee255;
ee307 = (ee63 * ee125 + ee27 * (ee183 + ee11 * ee116) - 
  (ee251 * ee67 + ee95 * ee126));
ee308 = ee150 * ee160;
ee309 = ee150 * ee139;
ee310 = ee218 * ee139;
ee319 = ee62 * ee211;
ee320 = ee62 * ee160;
ee321 = ee62 * ee195;
ee322 = ee62 * ee139;
ee323 = (ee27 * (ee11 * ee127) - ee225 - (ee259 * ee67 + 
  ee225));
ee324 = ee221 * ee138;
ee325 = ee63 * ee169;
ee328 = (ee127/ee5 + ee161/ee28 + (ee205/ee28 + ee236/ee51));
ee329 = ee205 * ee96;
ee330 = (ee52/ee5 + ee130/ee28);
ee332 = (ee274 + ee287);
ee333 = (ee128/ee5 + ee162/ee28 - (ee212/ee28 + ee239/ee51));
ee334 = (ee197/ee28 - ee241/ee51);
ee335 = (ee109 - ((ee93 + ee37 * ee76)/ee32 - ee93 * ee85/ee84));
ee337 = ee116 * ee19/ee28;
ee339 = ee116 * ee31 + ee52 * ee56;
ee340 = ee192 * ee65;
ee341 = ee117 * ee158;
ee342 = ee117 * ee138;
ee343 = ee95 * ee166;
ee344 = ee194 * ee65;
ee345 = ee178 * ee65;
ee346 = (ee234/ee87 + ((ee234 + ee33 * (2 * (2 * ee29 + 
  ee33 * (ee33))))/ee87 - ee98 * (2 * (ee83 * ee81))/(ee87 * ee87)));
ee347 = ee270 - ee267;
ee348 = (ee111 - ee1 * ee30/ee32);
ee349 = (ee111 * ee31 + ee49 * ee56);
ee350 = (ee118 - ((ee94 + (ee91 + ee13 * (2 * (ee71 + ee72))))/ee32 - 
  ee231/ee84 - ((ee231 + ee36 * (2 * (ee74 * ee23 + ee22 * 
  ee30)))/ee84 - ee112 * (2 * (ee78 * ee75))/(ee84 * ee84))));
ee351 = ee272 + ee273;
ee352 = (ee18/ee5);
ee353 = (ee120 - ((ee102 + ee42 * ee76)/ee32 - ee102 * ee85/ee84));
ee354 = (ee49/ee5 - ee143/ee28);
ee355 = ee212 * ee96;
ee356 = ee197 * ee96;
ee357 = (ee215/ee51);
ee358 = (ee103/ee28);
ee359 = (ee26/ee5);
ee360 = ee41 * ee56;
ee361 = ee275/ee51;
ee362 = ee277 - ee276;
ee363 = ee50 * ee56;
ee364 = ee278/ee51;
ee365 = ee50 * ee64;
ee366 = ee17 * ee252;
ee367 = ee17 * ee261;
ee368 = ee17 * ee116;
ee369 = ee279 + ee283;
ee370 = ee17 * ee206;
ee371 = ee280 + ee284;
ee372 = ee281 + ee286;
ee373 = ee17 * ee210;
ee374 = ee17 * ee155;
ee375 = ee17 * ee271;
ee377 = ee111 * ee19/ee28;
ee378 = log1p(ee4);
ee379 = ee49 * ee70;
    
out(j, 0) = -(ee92 * ee357 + ee266 * ee124);
out(j, 1) = ee92 * 
        ee334 - (ee208 * ee138 + ee221 * ee124);
out(j, 2) = -(ee92 * 
        ee261 - ee17 * ee358 - (ee218 * ee124 + ee208 * ee139));
out(j, 3) = ee92 * ee333 + (ee44 * ee187 - ee342 + 
            (ee151 * ee124 - ee342));
out(j, 4) = ee92 * ee271 - 
            ee17 * ee354 + (ee62 * ee138 + ee44 * ee195 + (ee150 * 
            ee124 + ee117 * ee139));
out(j, 5) = -(ee114 * ee124 + 
            ee322 + (ee322 + ee44 * ee180) - (ee92 * ee228 + 
            ee373 + (ee373 - ee43 * ee359)));
out(j, 6) = -(ee92 * 
            ee328 + (ee44 * ee202 - ee341 - (ee151 * ee176 + 
            ee341)));
out(j, 7) = -(ee92 * ee252 - ee17 * 
            ee330 + (ee62 * ee158 + ee44 * ee211 - (ee150 * ee176 + 
            ee117 * ee160)));
out(j, 8) = -(ee114 * ee176 + 
            ee320 + (ee320 + ee44 * ee179) - (ee92 * ee224 + 
            ee370 + (ee370 - ee43 * ee352)));
out(j, 9) = ee303 + 
            (ee371 + (ee104 * ee378 + ee284) + (ee371 + (ee92 * 
                ee207 + ee280)));
            
out(j, 10) = ee135 * (ee39 * 
            ee40) * ee167 * ee100 * ee124 + ee92 * (ee103 * (2 * 
            (ee103))/ee51 - ee215 * (2 * (ee174 * ee46))/ee134);
out(j, 11) = ee266 * ee138 + (ee209 * ee166 + (ee45 * 
            (ee16 * ee70) + ee226 * ee167) * ee100) * ee124 + 
            ee92 * ((ee103 * (2 * (ee143 + ee49 * ee5)) + ee197 * 
                ee181)/ee51 - ee215 * ee229/ee134);
out(j, 12) = -(ee92 * 
            ((ee178 * ee181 + ee103 * (2 * (ee50 * ee5 + ee145)))/ee51 - 
                ee215 * ee154/ee134) - ee17 * ee357 + (((ee152 * 
            ee167 + ee45 * (ee97 + ee16 * ee64)) * ee100 + ee209 * 
            ee169) * ee124 + ee266 * ee139));
out(j, 13) = ee92 * 
            ((ee26 * ee128 - ee379 + (ee128 * ee40 - ee379))/ee28 + 
                ee356/ee51 - ((ee103 * ee268 - ee356)/ee51 + 
                ee241 * ee229/ee134)) - (ee208 * ee187 - ee324 + 
            ((ee27 * (ee11 * ee128) - ee343 + (ee259 * ee100 - 
                ee343)) * ee124 - ee324));
out(j, 14) = ee92 * 
            ((ee50 * ee70 + ee26 * ee155 + (ee111 * ee40 + ee49 * 
                ee64))/ee28 - ee197 * ee65/ee51 - ((ee178 * ee96 + 
                ee103 * ee269)/ee51 - ee241 * ee154/ee134)) - 
            ee17 * ee334 - (ee218 * ee138 + ee208 * ee195 + ((ee63 * 
            ee166 + ee27 * (ee184 + ee11 * ee155) + (ee251 * 
            ee100 + ee95 * ee169)) * ee124 + ee221 * ee139));
out(j, 15) = -(ee92 * ((ee365 + ee365)/ee28 - ee345/ee51 - 
            ((ee345 + ee103 * ee196)/ee51 - ee240 * ee154/ee134)) - 
            ee367 - (ee367 - ee43 * ee358) - ((ee189 * ee100 + 
            ee325 + (ee325 + ee27 * ee180)) * ee124 + ee310 + 
            (ee310 + ee208 * ee180)));
out(j, 16) = ee92 * 
            (ee353/ee5 + ee274/ee28 + (ee332/ee28 + ee162 * ee96/ee51) - 
                ((ee26 * ee127 - ee287 - ee332)/ee28 + ee355/ee51 + 
                  ((ee355 + ee143 * ee268)/ee51 + ee239 * ee229/ee134))) + 
            (ee44 * (ee8 * ee353) - ee265 - ee299 + (ee323 * 
                ee124 - ee258 - ee299));
out(j, 17) = ee92 * 
            (ee348/ee5 - ee128 * ee19/ee28 + (ee349/ee28 - ee162 * 
                ee65/ee51) - ((ee50 * ee52 + ee26 * ee116 - ee349)/ee28 - 
                ee212 * ee65/ee51 + ((ee194 * ee96 + ee143 * 
                ee269)/ee51 - ee239 * ee154/ee134))) - ee17 * 
            ee333 + (ee62 * ee187 + ee44 * (ee17 * ee128 + ee8 * 
            ee348) - ee293 + (ee307 * ee124 + ee151 * ee139 - 
            ee293));
out(j, 18) = ee114 * ee138 + ee321 + 
            (ee321 + ee44 * (ee374 - ee43 * ee70 + ee374)) + 
            (ee294 * ee124 + ee309 + (ee309 + ee117 * ee180)) - 
            (ee92 * (ee377 + (ee377 - ee237 * ee65/ee51) + ((ee363 + 
                ee363)/ee28 - ee344/ee51 - ((ee344 + ee143 * 
                ee196)/ee51 - ee238 * ee154/ee134))) + ee375 + 
                (ee375 - ee43 * ee354));
out(j, 19) = -(ee246 * 
            ee124 + ee250 + (ee306) + (ee306 + (ee255 - ee44 * 
            (ee285 + (ee285 + (ee104 * ee40 + ee285))))) + (ee372 + 
            (ee104 * ee359 + ee286) + (ee372 + (ee92 * (ee364 + 
            ((ee278 + ee145 * ee196)/ee51 - ee214 * ee154/ee134) + 
            ee364) + ee281))));
out(j, 20) = -(ee92 * 
            (ee350/ee5 + ee272/ee28 + ((ee351)/ee28 + ee161 * 
                ee96/ee51) + ((ee351 + (ee273 + ee18 * ee127))/ee28 + 
                ee329/ee51 + ((ee329 + ee130 * ee268)/ee51 + 
                ee236 * ee229/ee134))) + (ee44 * (ee8 * ee350) - 
            ee264 - ee298 - (ee323 * ee176 + ee257 + ee298)));
out(j, 21) = -(ee92 * (ee335/ee5 - ee127 * ee19/ee28 + 
            ((ee339)/ee28 - ee161 * ee65/ee51) + ((ee339 + (ee41 * 
            ee52 + ee18 * ee116))/ee28 - ee205 * ee65/ee51 + 
            ((ee192 * ee96 + ee130 * ee269)/ee51 - ee236 * ee154/ee134))) - 
            ee17 * ee328 + (ee62 * ee202 + ee44 * (ee17 * ee127 + 
            ee8 * ee335) - ee292 - (ee307 * ee176 + ee151 * ee160 + 
            ee292)));
out(j, 22) = -(ee92 * ((ee360 + ee360)/ee28 - 
            ee340/ee51 - ((ee340 + ee130 * ee196)/ee51 - ee235 * 
            ee154/ee134) - (ee337 + (ee337 - ee232 * ee65/ee51))) - 
            ee366 - (ee366 - ee43 * ee330) + (ee114 * ee158 + 
            ee319 + (ee319 + ee44 * (ee368 - ee43 * ee52 + ee368)) - 
            (ee294 * ee176 + ee308 + (ee308 + ee117 * ee179))));
out(j, 23) = -(ee246 * ee176 + ee249 + (ee305) + (ee305 + 
            (ee254 - ee44 * (ee282 + (ee282 + (ee104 * ee31 + 
                ee282))))) + (ee369 + (ee104 * ee352 + ee283) + 
            (ee369 + (ee92 * (ee361 + ((ee275 + ee142 * ee196)/ee51 - 
                ee213 * ee154/ee134) + ee361) + ee279))));
out(j, 24) = ee246 * 
            ee106 + ee248 + (ee304) + (ee304 - (ee256 + ee44 * 
            (ee104 * ee19))) + ((ee303) * ee34 + ee247 + (ee296) + 
            (ee296 - (ee253 + ee107 * (ee362 - (ee276 + ee198 * 
                ee24) + (ee362 + (ee277 - ee25 * ee346)))))) + 
            (ee347 + (ee230 - ee346 * ee378 + (ee270)) + (ee347 + 
                (ee92 * (ee73 * ee196/ee51 - ee182 * ee154/ee134) - 
                  ee233 - ee267)));

} else {

ee1 = exp(lpsi);
ee2 = y - mu;
ee3 = ee2/ee1;
ee5 = exp(-ee3);
ee6 = ee3 - 3;
ee8 = (((ee3 - 6) * ee2/ee1 + 7) * ee2/ee1 - 1) * ee5 + 1;
ee10 = (ee6 * ee2/ee1 + 1) * ee5 - 1;
ee11 = R_pow(ee1, 2);
ee12 = R_pow(ee1, 3);

// third derivatives
// 1=location, 2=log(scale), 3=shape
// order: 111, 112, 113, 122, 123, 133, 222, 223, 233, 333
out(j, 0) = ee5/ee12;
out(j, 1) = (ee3 - 2) * ee5/ee11;
out(j, 2) = 0;
out(j, 3) = ee10/ee1;
out(j, 4) = 0;
out(j, 5) = 0;
out(j, 6) = ee10 * ee2/ee1;
out(j, 7) = 0;
out(j, 8) = 0;
out(j, 9) = 0;

// fourth derivatives
// 1=location, 2=log(scale), 3=shape
// order: 1111, 1112, 1113, 1122, 1123, 1133, 1222, 1223, 1233, 1333
//          2222, 2223, 2233, 2333, 3333
out(j, 10) = ee5/R_pow(ee1, 4);
out(j, 11) = ee6 * ee5/ee12;
out(j, 12) = 0;
out(j, 13) = ((ee3 - 5) * ee2/ee1 + 4) * ee5/ee11;
out(j, 14) = 0;
out(j, 15) = 0;
out(j, 16) = ee8/ee1;
out(j, 17) = 0;
out(j, 18) = 0;
out(j, 19) = 0;
out(j, 20) = ee8 * ee2/ee1;
out(j, 21) = 0;
out(j, 22) = 0;
out(j, 23) = 0;
out(j, 24) = 0;

}

}

return out;

}
