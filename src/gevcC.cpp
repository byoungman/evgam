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
double gevcd0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = ymat.n_rows;

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    xivec = xivec.elem(dupid);
}

double yl, yh, mu, lpsi, xi;
double ee1, ee2, ee3, ee4;
double nllh=0.0;

for (int j=0; j < nobs; j++) {

yl = ymat(j, 0);
yh = ymat(j, 1);
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];

if (fabs(xi) >= xieps) {

ee1 = -(1/xi);
ee2 = exp(lpsi);
ee3 = xi * (yh - mu)/ee2;

if (ee3 <= -1.0) {
    nllh = 1e20;
    break;
} else {

ee4 = xi * (yl - mu)/ee2;

if (ee4 <= -1.0) {
    nllh = 1e20;
    break;
} else {

nllh -= log(exp(-R_pow(1 + ee3, ee1)) - exp(-R_pow(1 + ee4, ee1)));

}
}

} else {

ee1 = exp(lpsi);
nllh -= log(exp(-exp(-((yh - mu)/ee1))) - exp(-exp(-((yl - mu)/ee1))));

}

}

return(nllh);

}

// //' @rdname gevd0
// [[Rcpp::export]]
arma::mat gevcd12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = ymat.n_rows;
arma::mat out = arma::mat(nobs, 9);

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    xivec = xivec.elem(dupid);
}

double yl, yh, mu, lpsi, xi;

double ee1, ee2, ee3, ee4, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee20, ee21, ee22, ee23, ee24, ee27, ee28, ee29;
double ee31, ee33, ee34, ee35, ee36, ee37, ee38;
double ee40, ee46, ee48;
double ee50, ee52, ee55, ee57, ee58;
double ee60, ee61, ee62, ee63, ee64, ee65, ee66;

for (int j=0; j < nobs; j++) {

yl = ymat(j, 0);
yh = ymat(j, 1);
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];

if (fabs(xi) >= xieps) {

ee1 = exp(lpsi);
ee2 = 1/xi;
ee3 = yh - mu;
ee4 = yl - mu;
ee6 = xi * ee3/ee1;
ee8 = xi * ee4/ee1;
ee9 = 1 + ee6;
ee10 = 1 + ee8;
ee11 = -ee2;
ee12 = 1 + ee2;
ee15 = exp( - R_pow(ee9, ee11));
ee16 = exp( - R_pow(ee10, ee11));
ee17 = R_pow(ee9, ee12);
ee18 = R_pow(ee10, ee12);
ee19 = ee15 - ee16;
ee20 = ee2 + 2;
ee21 = log1p(ee6);
ee22 = log1p(ee8);
ee23 = R_pow(ee9, ee2);
ee24 = R_pow(ee10, ee2);
ee27 = ee19 * ee1;
ee28 = ee3/(ee17 * ee1);
ee29 = ee4/(ee18 * ee1);
ee31 = ee21/(xi * ee23) - ee28;
ee33 = ee22/(xi * ee24) - ee29;
ee34 = R_pow(ee9, ee20);
ee35 = R_pow(ee10, ee20);
ee36 = 2 * ee12;
ee37 = xi * ee12;
ee38 = xi * ee19;
ee40 = ee15 * ee31 - ee16 * ee33;
ee46 = ee16 * ee4/ee18 - ee15 * ee3/ee17;
ee48 = ee16/ee18 - ee15/ee17;
ee50 = ee12 * ee3/(ee34 * ee1);
ee52 = ee12 * ee4/(ee35 * ee1);
ee55 = 1/ee17;
ee57 = 1/R_pow(ee9, ee36) - ee37/ee34;
ee58 = 1/ee18;
ee60 = 1/R_pow(ee10, ee36) - ee37/ee35;
ee61 = ee21/(xi * ee17);
ee62 = ee22/(xi * ee18);
ee63 = ((ee31/ee17 - ee61)/xi + ee50) * ee15;
ee64 = ((ee33/ee18 - ee62)/xi + ee52) * ee16;
ee65 = (ee57 * ee3/ee1 + ee55) * ee15;
ee66 = (ee60 * ee4/ee1 + ee58) * ee16;

out(j, 0) = -(ee48/ee27);
out(j, 1) = -(ee46/ee27);
out(j, 2) = ee40/ee38;
out(j, 3) =  - ((ee57 * ee15 - (ee60 * ee16 + R_pow(ee48, 2)/ee19))/(ee19 * R_pow(ee1, 2)));
out(j, 4) = -((ee65 - (ee66 + ee46 * ee48/ee27))/ee27);
out(j, 5) = -((ee63 + ee40 * ee48/ee38 - ee64)/ee27);
out(j, 6) =  - ((ee65 * ee3 - (ee66 * ee4 + R_pow(ee46, 2)/ee27))/ee27);
out(j, 7) = -((ee63 * ee3 + ee40 * ee46/ee38 - ee64 * ee4)/ee27);
out(j, 8) = (((((1 - 1/ee23) * ee21/xi + ee28) * ee31 + (ee3/(ee9 * ee1) -
   2 * (ee21/xi))/ee23)/xi + (ee50 + (ee55 - ee61)/xi) * ee3/ee1) * ee15 +
   R_pow(ee40, 2)/ee38 - ((((1 - 1/ee24) * ee22/xi +
   ee29) * ee33 + (ee4/(ee10 * ee1) - 2 * (ee22/xi))/ee24)/xi +
   (ee52 + (ee58 - ee62)/xi) * ee4/ee1) * ee16)/ee38;

} else {
    
ee1 = exp(lpsi);
ee2 = yh - mu;
ee3 = yl - mu;
ee6 = exp(-(ee2/ee1));
ee7 = exp(-(ee3/ee1));
ee10 = exp(-ee6);
ee11 = exp(-ee7);
ee12 = ee10 - ee11;
ee13 = ee12 * ee1;
ee14 = ee6 * ee10;
ee15 = ee7 * ee11;
ee16 = 1 - ee6;
ee17 = 1 - ee7;
ee18 = ee14 - ee15;
ee20 = ee14 * ee2 - ee15 * ee3;
ee22 = (ee16 * ee2/ee1 - 1) * ee6 * ee10;
ee24 = (ee17 * ee3/ee1 - 1) * ee7 * ee11;

out(j, 0) = ee18/ee13;
out(j, 1) = ee20/ee13;
out(j, 2) = 0;
out(j, 3) = (ee16 * ee6 * ee10 + R_pow(ee18, 2)/ee12 - ee17 * ee7 * ee11)/(ee12 * R_pow(ee1, 2));
out(j, 4) = (ee22 + ee18 * ee20/ee13 - ee24)/ee13;
out(j, 5) = 0;
out(j, 6) = (ee22 * ee2 + R_pow(ee20, 2)/ee13 - ee24 * ee3)/ee13;
out(j, 7) = 0;
out(j, 8) = 0;

}

}

return out;

}

// //' @rdname gevd0
// [[Rcpp::export]]
arma::mat gevcd34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec xivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = ymat.n_rows;
arma::mat out = arma::mat(nobs, 25);

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    xivec = xivec.elem(dupid);
}

double yl, yh, mu, lpsi, xi;

double ee1, ee2, ee3, ee4, ee5, ee6, ee8, ee9;
double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee20, ee22, ee23, ee24, ee25, ee26, ee27, ee28;
double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee38;
double ee40, ee42, ee43, ee44, ee45, ee46, ee47, ee48, ee49;
double ee50, ee51, ee52, ee53, ee54, ee55, ee56, ee57, ee58, ee59;
double ee60, ee61, ee62, ee64, ee65, ee66, ee67, ee68, ee69;
double ee70, ee71, ee72, ee73, ee74, ee75, ee76, ee77, ee78, ee79;
double ee80, ee81, ee82, ee84, ee85, ee86, ee87, ee88, ee89;
double ee90, ee91, ee92, ee93, ee94, ee95, ee96, ee97, ee98, ee99;
double ee101, ee103, ee104, ee105, ee106, ee107, ee109;
double ee110, ee111, ee112, ee114, ee116, ee117, ee118, ee119;
double ee120, ee122, ee124, ee125, ee126, ee127, ee128, ee129;
double ee130, ee132, ee134, ee135, ee136, ee137, ee138, ee139;
double ee140, ee141, ee142, ee143, ee144, ee145, ee146, ee147, ee148, ee149;
double ee150, ee151, ee152, ee153, ee154, ee155, ee156, ee158;
double ee160, ee161, ee162, ee163, ee164, ee165, ee166, ee167, ee168, ee169;
double ee170, ee171, ee172, ee175, ee176, ee177, ee178, ee179;
double ee180, ee181, ee182, ee183, ee184, ee185, ee186, ee187, ee188, ee189;
double ee190, ee191, ee194, ee197, ee199;
double ee201, ee202, ee203, ee209;
double ee211, ee213, ee215, ee216, ee217, ee218, ee219;
double ee220, ee221, ee222, ee223, ee224, ee226, ee228, ee229;
double ee230, ee231, ee232, ee238;
double ee240, ee241, ee242, ee243, ee244, ee245, ee247, ee249;
double ee250, ee252, ee254, ee255, ee258, ee259;
double ee260, ee261, ee265, ee268, ee269;
double ee270, ee271, ee272, ee278;
double ee280;
double ee291, ee294, ee295, ee296, ee297, ee298;
double ee300, ee302, ee303, ee304, ee309;
double ee311, ee313, ee314, ee315, ee316, ee317, ee318;
double ee320, ee322, ee324, ee326, ee327, ee328;
double ee330, ee331, ee333, ee334, ee336, ee338, ee339;
double ee340, ee341, ee343, ee344, ee345, ee346, ee348;
double ee350, ee351, ee353, ee354, ee355, ee356, ee357, ee358, ee359;
double ee360, ee361, ee362, ee363, ee364, ee366, ee367, ee368, ee369;
double ee370, ee371, ee372, ee373, ee374, ee375, ee376, ee377, ee378;
double ee381, ee382, ee384, ee385, ee386, ee387, ee388;
double ee390, ee391, ee393, ee394, ee395, ee397, ee399;
double ee401, ee403, ee404, ee405, ee406;
double ee411, ee412, ee413, ee414, ee415, ee416, ee417, ee418;
double ee420, ee421, ee422, ee423, ee424, ee425, ee426, ee427;
double ee432, ee433, ee434, ee435, ee436, ee437, ee438, ee439;
double ee440, ee441, ee442, ee443, ee444, ee446, ee448, ee449;
double ee450, ee451, ee452, ee453, ee454, ee455, ee457, ee459;
double ee460, ee461, ee463, ee464, ee465, ee466, ee468;
double ee470, ee471, ee472, ee474, ee476, ee478, ee479;
double ee481, ee482, ee483, ee485, ee486, ee487, ee489;
double ee491, ee492, ee493, ee494, ee495, ee496, ee497, ee498;
double ee501, ee504, ee505, ee506, ee509;
double ee510, ee511, ee512, ee513, ee514, ee517, ee518, ee519;
double ee520, ee521, ee522, ee523, ee524, ee525, ee526, ee527, ee528, ee529;
double ee530, ee531, ee532, ee533, ee534, ee535, ee536;

for (int j=0; j < nobs; j++) {

yl = ymat(j, 0);
yh = ymat(j, 1);
mu = muvec[j];
lpsi = lpsivec[j];
xi = xivec[j];

if (fabs(xi) >= xieps) {

ee1 = exp(lpsi);
ee2 = 1/xi;
ee3 = yh - mu;
ee4 = yl - mu;
ee6 = xi * ee3/ee1;
ee8 = xi * ee4/ee1;
ee9 = 1 + ee2;
ee10 = 1 + ee6;
ee11 = 1 + ee8;
ee12 = ee2 + 2;
ee13 = R_pow(ee10, ee9);
ee14 = R_pow(ee11, ee9);
ee15 = -ee2;
ee16 = log1p(ee6);
ee17 = log1p(ee8);
ee18 = R_pow(ee10, ee12);
ee19 = R_pow(ee11, ee12);
ee22 = exp( - R_pow(ee10, ee15));
ee23 = exp( - R_pow(ee11, ee15));
ee24 = R_pow(xi, 2);
ee25 = R_pow(ee10, ee2);
ee26 = R_pow(ee11, ee2);
ee27 = xi * ee9;
ee28 = ee2 + 3;
ee31 = ee3/(ee13 * ee1);
ee32 = ee4/(ee14 * ee1);
ee33 = ee18 * ee1;
ee34 = ee19 * ee1;
ee36 = ee16/(xi * ee25) - ee31;
ee38 = ee17/(xi * ee26) - ee32;
ee40 = ee9 * ee3/ee33;
ee42 = ee9 * ee4/ee34;
ee43 = 1/ee13;
ee44 = 1/ee14;
ee45 = 2 * ee9;
ee46 = ee10 * ee1;
ee47 = ee11 * ee1;
ee48 = ee16/(xi * ee13);
ee49 = ee17/(xi * ee14);
ee50 = ee3/ee46;
ee51 = ee4/ee47;
ee52 = R_pow(ee10, ee28);
ee53 = R_pow(ee11, ee28);
ee54 = ee22 - ee23;
ee55 = ee16/xi;
ee56 = ee17/xi;
ee57 = ee50 - 2 * ee55;
ee58 = ee51 - 2 * ee56;
ee59 = ee24 * ee18;
ee60 = ee24 * ee19;
ee61 = ee52 * ee1;
ee62 = ee53 * ee1;
ee66 = 1/R_pow(ee10, ee45) - ee27/ee18;
ee68 = 1/R_pow(ee11, ee45) - ee27/ee19;
ee70 = ee12 * ee3/ee61;
ee72 = ee12 * ee4/ee62;
ee73 = xi * ee12;
ee74 = ee16/(ee24 * ee13);
ee75 = ee17/(ee24 * ee14);
ee76 = ee16/ee59;
ee77 = ee17/ee60;
ee78 = 1/ee18;
ee79 = 1/ee19;
ee80 = ee76 - ee70;
ee81 = ee77 - ee72;
ee84 = ee66 * ee3/ee1 + ee43;
ee87 = ee68 * ee4/ee1 + ee44;
ee89 = ee22 * ee36 - ee23 * ee38;
ee90 = ee74 - ee40;
ee91 = ee75 - ee42;
ee93 = (ee36/ee13 - ee48)/xi + ee40;
ee95 = (ee38/ee14 - ee49)/xi + ee42;
ee97 = (ee40 + (ee43 - ee48)/xi) * ee3/ee1;
ee99 = (ee42 + (ee44 - ee49)/xi) * ee4/ee1;
ee104 = ee23 * ee4/ee14 - ee22 * ee3/ee13;
ee105 = 1/ee25;
ee106 = 1/ee26;
ee107 = 2/xi;
ee110 = ee23/ee14 - ee22/ee13;
ee111 = ee57/ee25;
ee112 = ee58/ee26;
ee114 = (((1 - ee105) * ee16/xi + ee31) * ee36 + ee111)/xi +  ee97;
ee116 = (((1 - ee106) * ee17/xi + ee32) * ee38 + ee112)/xi +  ee99;
ee117 = ee93 * ee22;
ee118 = ee95 * ee23;
ee119 = ee84 * ee22;
ee120 = ee87 * ee23;
ee122 = ee27 * ee3/ee33;
ee124 = ee27 * ee4/ee34;
ee125 = 2/ee18;
ee126 = 2/ee19;
ee127 = ee57/ee13;
ee128 = ee58/ee14;
ee129 = ee9 * ee80;
ee130 = ee9 * ee81;
ee132 = (1/ee59 - ee129) * ee3/ee1;
ee134 = (1/ee60 - ee130) * ee4/ee1;
ee135 = ee2 + 4;
ee136 = xi * ee54;
ee137 = ee54 * ee1;
ee138 = ee114 * ee22;
ee139 = ee116 * ee23;
ee140 = ee16 * ee90;
ee141 = ee17 * ee91;
ee142 = xi * ee80;
ee143 = xi * ee81;
ee144 = ee127 + ee140;
ee145 = ee128 + ee141;
ee146 = ee107 + 3;
ee147 = ee119 * ee3;
ee148 = ee120 * ee4;
ee149 = ee122 - ee43;
ee150 = ee124 - ee44;
ee151 = 2/ee13;
ee152 = 2/ee14;
ee153 = ee117 * ee3;
ee154 = ee118 * ee4;
ee155 = ee78 + ee125;
ee156 = ee79 + ee126;
ee158 = ee73 * ee3/ee61;
ee160 = ee73 * ee4/ee62;
ee161 = ee138 - ee139;
ee162 = ee147 - ee148;
ee163 = ee84/ee13;
ee164 = ee87/ee14;
ee165 = ee9 * (ee78 + ee142);
ee166 = ee9 * (ee79 + ee143);
ee167 = ee66 * ee22;
ee168 = ee68 * ee23;
ee169 = ee93/ee13;
ee170 = ee95/ee14;
ee171 = ee90/ee13;
ee172 = ee91/ee14;
ee175 = ee27 * (ee155 - ee158);
ee176 = ee27 * (ee156 - ee160);
ee177 = ee117 - ee118;
ee178 = ee153 - ee154;
ee179 = ee167 - ee168;
ee180 = ee43 - ee122;
ee181 = ee44 - ee124;
ee182 = 2 * ee50;
ee183 = 2 * ee51;
ee184 = ee119 - ee120;
ee185 = R_pow(ee89, 2);
ee186 = ee24 * ee52;
ee187 = ee24 * ee53;
ee188 = ee144/xi;
ee189 = ee145/xi;
ee190 = R_pow(ee10, ee135);
ee191 = R_pow(ee11, ee135);
ee194 = (ee50 + ee107) * ee3/ee46 + (ee182 - 6 * ee55)/xi;
ee197 = (ee51 + ee107) * ee4/ee47 + (ee183 - 6 * ee56)/xi;
ee199 = (ee114/ee13 + 2 * (ee36 * ee90) - ee188)/xi - ee132;
ee201 = (ee116/ee14 + 2 * (ee38 * ee91) - ee189)/xi - ee134;
ee202 = ee190 * ee1;
ee203 = ee191 * ee1;
ee209 = (ee111 + ee16 * ee36/xi)/xi + ee97;
ee211 = (ee112 + ee17 * ee38/xi)/xi + ee99;
ee213 = (ee171 - (ee169 + ee165 + ee78)) * ee3/ee1;
ee215 = (ee172 - (ee170 + ee166 + ee79)) * ee4/ee1;
ee216 = R_pow(ee10, ee146);
ee217 = R_pow(ee11, ee146);
ee218 = ee89 * ee104;
ee219 = ee149/ee13;
ee220 = ee150/ee14;
ee221 = ee73/ee52;
ee222 = ee73/ee53;
ee223 = ee199 * ee22;
ee224 = ee201 * ee23;
ee226 = ((((ee40 + (ee151 - ee48)/xi) * ee16 - (ee127 +  ee151))/xi - 2 * ee40)/xi - ee132) * ee3/ee1;
ee228 = ((((ee42 + (ee152 - ee49)/xi) * ee17 - (ee128 +  ee152))/xi - 2 * ee42)/xi - ee134) * ee4/ee1;
ee229 = (ee209 * ee16 + 2 * (ee57 * ee36))/xi;
ee230 = (ee211 * ee17 + 2 * (ee58 * ee38))/xi;
ee231 = ee194/ee25;
ee232 = ee197/ee26;
ee238 = ee16/ee186 - ee28 * ee3/ee202;
ee240 = ee17/ee187 - ee28 * ee4/ee203;
ee241 = xi * ee28;
ee242 = ee185/xi;
ee243 = R_pow(ee104, 2);
ee244 = ((((3 - ee105) * ee16/xi + ee31) * ee36 + ee57 *  (ee105 + 2/ee25))/xi + 3 * ee97) * ee36;
ee245 = ((((3 - ee106) * ee17/xi + ee32) * ee38 + ee58 *  (ee106 + 2/ee26))/xi + 3 * ee99) * ee38;
ee247 = ee144/ee24 + ee132;
ee249 = ee145/ee24 + ee134;
ee250 = ee89 * ee110;
ee252 = (ee229 - (ee244 + ee231))/xi + ee226;
ee254 = (ee230 - (ee245 + ee232))/xi + ee228;
ee255 = ee161 * ee54;
ee258 = ee213 + (ee48 - ee180 * ee36)/xi;
ee259 = ee215 + (ee49 - ee181 * ee38)/xi;
ee260 = ee165 + ee78;
ee261 = ee166 + ee79;
ee265 = (ee175 - (ee163 + 2 * (ee180/ee13))) * ee3/ee1 -  ee43;
ee268 = (ee176 - (ee164 + 2 * (ee181/ee14))) * ee4/ee1 -  ee44;
ee269 = 2 * (ee185/ee136);
ee270 = 2/ee52;
ee271 = 2/ee53;
ee272 = ee255 - ee242;
ee278 = (ee57/ee18 + ee16 * ee80)/ee24 + (1/ee186 - ee12 *  ee238) * ee3/ee1;
ee280 = (ee58/ee19 + ee17 * ee81)/ee24 + (1/ee187 - ee12 *  ee240) * ee4/ee1;
ee291 = (2 * ee219 + ee175 - ee163) * ee3/ee1 - ee43;
ee294 = (2 * ee220 + ee176 - ee164) * ee4/ee1 - ee44;
ee295 = ee218/xi;
ee296 = ee104 * ee110;
ee297 = ee171 - (ee169 + ee9 * (ee142 - ee36/ee18) + ee78);
ee298 = ee172 - (ee170 + ee9 * (ee143 - ee38/ee19) + ee79);
ee300 = ee219 + xi * ((1/ee216 - ee221) * ee3/ee1 + ee125) *  ee9 - ee163;
ee302 = ee220 + xi * ((1/ee217 - ee222) * ee4/ee1 + ee126) *  ee9 - ee164;
ee303 = 2 * ee70;
ee304 = 2 * ee72;
ee309 = ee107 + 4;
ee311 = ee27 * (2/ee216 - ee221) - ee66/ee13;
ee313 = ee27 * (2/ee217 - ee222) - ee68/ee14;
ee314 = ee252 * ee22;
ee315 = ee254 * ee23;
ee316 = ee223 * ee3;
ee317 = ee224 * ee4;
ee318 = ee178 * ee54;
ee320 = ee258 * ee22 * ee3;
ee322 = ee259 * ee23 * ee4;
ee324 = ee260 * ee3/ee1;
ee326 = ee261 * ee4/ee1;
ee327 = ee291 * ee22;
ee328 = ee294 * ee23;
ee330 = (ee36 * ee149 + ee48)/xi + ee213;
ee331 = ee297 * ee22;
ee333 = (ee38 * ee150 + ee49)/xi + ee215;
ee334 = ee298 * ee23;
ee336 = ee265 * ee22 * ee3;
ee338 = ee268 * ee23 * ee4;
ee339 = ee300 * ee22;
ee340 = ee302 * ee23;
ee341 = ee54 * R_pow(ee1, 2);
ee343 = ee243/ee1;
ee344 = 2 * (ee218/ee136);
ee345 = ee22 * ee311;
ee346 = ee23 * ee313;
ee348 = ee175 * ee3/ee1;
ee350 = ee176 * ee4/ee1;
ee351 = ee314 - ee315;
ee353 = ee138 + ee269 - ee139;
ee354 = ee177 * ee54;
ee355 = ee178 * ee110;
ee356 = ee184 * ee104;
ee357 = ee162 * ee54;
ee358 = ee162 * ee110;
ee359 = ee330 * ee22;
ee360 = ee333 * ee23;
ee361 = ee148 + 2 * (ee243/ee137);
ee362 = ee336 - ee338;
ee363 = ee12/ee52;
ee364 = ee12/ee53;
ee366 = ee250/xi;
ee367 = R_pow(ee110, 2);
ee368 = 2 * (ee250/ee136);
ee369 = 2 * (ee16/(xi * ee18));
ee370 = 2 * (ee17/(xi * ee19));
ee371 = 6/xi;
ee372 = ee345 - ee346;
ee373 = ee223 - ee224;
ee374 = ee316 - ee317;
ee375 = ee199/ee13;
ee376 = ee201/ee14;
ee377 = ee161 * ee104;
ee378 = ee161 * ee110;
ee381 = ee138 + 2 * ee161 + ee269 - ee139;
ee382 = ee177 * ee89;
ee384 = ee177 * ee104 + ee355;
ee385 = ee318 - ee295;
ee386 = ee178 * ee89;
ee387 = ee247 * ee16;
ee388 = ee247/ee13;
ee390 = ee249 * ee17;
ee391 = ee249/ee14;
ee393 = ee184 * ee54;
ee394 = ee356 + ee358;
ee395 = ee357 + ee343;
ee397 = (((ee369 - ee125)/xi - ee303)/ee24 - ee278 * ee9) *  ee3/ee1;
ee399 = (((ee370 - ee126)/xi - ee304)/ee24 - ee280 * ee9) *  ee4/ee1;
ee401 = ee320 - ee322;
ee403 = ee194/ee13;
ee404 = ee197/ee14;
ee405 = ee324 - ee74;
ee406 = ee326 - ee75;
ee411 = ee84 * ee90;
ee412 = ee87 * ee91;
ee413 = ee327 - ee328;
ee414 = ee331 - ee334;
ee415 = ee339 - ee340;
ee416 = R_pow(ee10, ee309);
ee417 = R_pow(ee11, ee309);
ee418 = ee54 * R_pow(ee1, 3);
ee420 = ee296/ee1;
ee421 = ee43 - ee348;
ee422 = ee78 + ee27 * ee80;
ee423 = ee44 - ee350;
ee424 = ee79 + ee27 * ee81;
ee425 = 2 * ee272;
ee426 = 2 * (ee93 * ee90);
ee427 = 2 * (ee95 * ee91);
ee432 = 2 * (ee296/ee137);
ee433 = ee125 - ee158;
ee434 = ee126 - ee160;
ee435 = 4 * ee354;
ee436 = 4/ee18;
ee437 = 4/ee52;
ee438 = 4/ee19;
ee439 = 4/ee53;
ee440 = 8 * ee295;
ee441 = 8 * ee366;
ee442 = 8 * ee242;
ee443 = xi * ee278;
ee444 = xi * ee280;
ee446 = ee241 * ee3/ee202;
ee448 = ee241 * ee4/ee203;
ee449 = xi * (2 * ee129 - (ee9 * (ee270 + ee73 * ee238) +  ee363) * ee3/ee1);
ee450 = xi * (2 * ee130 - (ee9 * (ee271 + ee73 * ee240) +  ee364) * ee4/ee1);
ee451 = ((ee252/ee13 + 3 * (ee114 * ee90) + 3 * (ee247 *  ee36) - (ee387 + 2 * (ee57 * ee90) - ee403)/xi)/xi -  ee397) * ee22;
ee452 = ((ee254/ee14 + 3 * (ee116 * ee91) + 3 * (ee249 *  ee38) - (ee390 + 2 * (ee58 * ee91) - ee404)/xi)/xi -  ee399) * ee23;
ee453 = ee272 * ee104;
ee454 = ee375 + ee9 * (ee76 + ee443 - ee70);
ee455 = ee376 + ee9 * (ee77 + ee444 - ee72);
ee457 = (ee229 - ee231)/xi + ee226;
ee459 = (ee230 - ee232)/xi + ee228;
ee460 = ee377 + ee386;
ee461 = ee378 + ee382;
ee463 = ee388 + ee303;
ee464 = ee391 + ee304;
ee465 = ee405/ee13;
ee466 = ee406/ee14;
ee468 = ee184 * ee89/xi;
ee470 = ee162 * ee89/xi;
ee471 = ee359 - ee360;
ee472 = ee93 * ee149;
ee474 = ee153 + ee344 - ee154;
ee476 = ee117 + ee368 - ee118;
ee478 = ee258/ee13 + ee411;
ee479 = ee95 * ee150;
ee481 = ee259/ee14 + ee412;
ee482 = ee84 * ee149;
ee483 = ee147 - ee361;
ee485 = ee147 + 2 * ee162 - ee361;
ee486 = ee84/ee18;
ee487 = ee179 * ee54;
ee489 = ee179 * ee89/xi;
ee491 = ee179 * ee104/ee1;
ee492 = ee87 * ee150;
ee493 = ee87/ee19;
ee494 = ee265/ee13;
ee495 = ee268/ee14;
ee496 = ee422/ee13;
ee497 = ee424/ee14;
ee498 = ee168 + 2 * (ee367/ee54);
ee501 = ee43 + ee151;
ee504 = 1/ee52;
ee505 = 1/ee416;
ee506 = ee44 + ee152;
ee509 = 1/ee53;
ee510 = 1/ee417;
ee511 = ee425 + ee442;
ee512 = 2 * (ee384 * ee54);
ee513 = 2 * ee385;
ee514 = 2 * ee395;
ee517 = 2 * ee153 + ee344 - 2 * ee154;
ee518 = 2 * (ee433/ee13);
ee519 = 2 * (ee434/ee14);
ee520 = ee182 + ee371;
ee521 = ee183 + ee371;
ee522 = 2 * (ee16/ee18);
ee523 = 2 * (ee17/ee19);
ee524 = ee435 + ee441;
ee525 = 4 * ee318;
ee526 = 4 * ee393;
ee527 = 8 * ee420;
ee528 = 8 * ee343;
ee529 = ee74 - ee324;
ee530 = ee75 - ee326;
ee531 = ee348 - ee43;
ee532 = ee27 * (ee155 + ee436 - ee73 * (ee270 + ee437 -  ee446) * ee3/ee1);
ee533 = ee350 - ee44;
ee534 = ee27 * (ee156 + ee438 - ee73 * (ee271 + ee439 -  ee448) * ee4/ee1);
ee535 = ee241/ee190;
ee536 = ee241/ee191;

out(j, 0) = -((ee345 - ((ee167 + 2 * ee179 - ee498) * ee110/ee54 +
   ee346))/ee418);
out(j, 1) = -((ee339 - ((ee491 + (2 * ee184 - ee432) * ee110)/ee54 +
   ee340))/ee341);
out(j, 2) = -(((ee489 - (2 * ee177 + ee368) * ee110)/ee54 +
   ee331 - ee334)/ee341);
out(j, 3) = -((ee327 - ((ee358 + (2 * ee119 - (2 * ee120 + ee432)) * ee104)/ee137 +
   ee328))/ee137);
out(j, 4) = -(((ee468 - (ee355 + ee476 * ee104)/ee1)/ee54 +
   ee359 - ee360)/ee137);
out(j, 5) = -(((ee378 + (2 * ee117 + ee368 - 2 * ee118) * ee89)/ee136 +
   ee223 - ee224)/ee137);
out(j, 6) = -((ee336 - (ee485 * ee104/ee137 + ee338))/ee137);
out(j, 7) = -(((ee470 - (2 * ee178 + ee344) * ee104/ee1)/ee54 +
   ee320 - ee322)/ee137);
out(j, 8) = -(((ee377 + ee517 * ee89)/ee136 + ee316 - ee317)/ee137);
out(j, 9) = (ee314 + ee381 * ee89/ee136 - ee315)/ee136;
out(j, 10) =  - (((ee313/ee14 + xi * (ee156 * ee68 - xi * (ee510 +
   2/ee417 - ee536) * ee12) * ee9) * ee23 - (((ee167 - ee498) * ee179 +
   (2 * ee372 - ((2 * (ee487 + ee367) + 4 * ee487 -
   8 * ee367)/ee54 + 4 * ee179) * ee110/ee54) * ee110 + 2 * (R_pow(ee179, 2) +
   ee372 * ee110))/ee54 + (ee311/ee13 + xi * (ee155 * ee66 -
   xi * (ee505 + 2/ee416 - ee535) * ee12) * ee9) * ee22))/(ee54 * R_pow(ee1, 4)));
out(j, 11) = -(((ee302/ee14 + ee68 * ee150 + ee27 * (2 * ee493 +
   ee519 - xi * ((ee510 - ee536) * ee4/ee1 + ee509 + ee271) * ee12)) * ee23 -
   (((ee119 - (ee120 + ee432)) * ee179 + (ee339 -
   (((2 * (ee393 + ee420) + ee526 - ee527) * ee110/ee54 +
   4 * ee491)/ee54 + ee340)) * ee110 + ee372 * ee104/ee1 + 2 * (ee184 * ee179 +
   ee415 * ee110))/ee54 + (ee300/ee13 + ee66 * ee149 +
   ee27 * (2 * ee486 + ee518 - xi * ((ee505 - ee535) * ee3/ee1 +
   ee504 + ee270) * ee12)) * ee22))/ee418);
out(j, 12) = -(((ee298/ee14 + ee68 * ee91 + xi * (ee9 * (2 * (ee95/ee19) +
   ee271 - (ee38/ee53 - xi * ee240) * ee12) + ee364) -
   2 * ee497) * ee23 + (ee89 * ee372/xi - (ee476 * ee179 +
   (ee331 + (4 * ee489 - (2 * (ee354 - ee366) + ee435 + ee441) * ee110/ee54)/ee54 -
   ee334) * ee110 + 2 * (ee177 * ee179 +
   ee414 * ee110)))/ee54 - (ee297/ee13 + ee66 * ee90 + xi * (ee9 * (2 * (ee93/ee18) +
   ee270 - (ee36/ee52 - xi * ee238) * ee12) +
   ee363) - 2 * ee496) * ee22)/ee418);
out(j, 13) =  - (((ee294/ee14 + 2 * ee492 + xi * ((ee493 + ee519 -
   xi * (ee509 + ee439 - ee448) * ee12) * ee4/ee1 + ee438) * ee9 -
   ee423/ee14) * ee23 - (((ee483 * ee179 + 2 * (ee415 * ee104) -
   (((ee526 - ee527) * ee104 + 2 * (ee395 * ee110))/ee54 +
   4 * ee356) * ee110/ee54)/ee1 + 2 * (R_pow(ee184, 2) +
   ee413 * ee110))/ee54 + (ee291/ee13 + 2 * ee482 + xi * ((ee486 +
   ee518 - xi * (ee504 + ee437 - ee446) * ee12) * ee3/ee1 +
   ee436) * ee9 - ee421/ee13) * ee22))/ee341);
out(j, 14) = -(((ee333/ee14 + ee479 + ee412 + ee9 * ee434 * ee38 +
   (xi * ee95 * ee9/ee19 - ee497) * ee4/ee1 - (ee466 + ee126 +
   ee450)) * ee23 + (ee415 * ee89/xi - ((ee474 * ee179 +
   ee414 * ee104)/ee1 + (4 * ee468 - (ee524 * ee104 + 2 * (ee385 * ee110))/ee137) * ee110/ee54 +
   2 * (ee471 * ee110 + ee177 * ee184)))/ee54 -
   (ee330/ee13 + ee472 + ee411 + ee9 * ee433 * ee36 +
   (xi * ee93 * ee9/ee18 - ee496) * ee3/ee1 - (ee465 +
   ee125 + ee449)) * ee22)/ee341);
out(j, 15) =  - ((((ee353 * ee179 + 2 * (ee414 * ee89) - ((ee524 * ee89 +
   2 * (ee272 * ee110))/ee54 + 4 * ee382) * ee110/ee54)/xi -
   2 * (ee373 * ee110 + R_pow(ee177, 2)))/ee54 + (ee463 -
   (ee375 + ee9 * (ee443 - ee114/ee18) + (ee369 - 2 * (ee422 * ee36))/xi +
   ee426)) * ee22 - (ee464 - (ee376 + ee9 * (ee444 -
   ee116/ee19) + (ee370 - 2 * (ee424 * ee38))/xi + ee427)) * ee23)/ee341);
out(j, 16) = -((((ee495 + 3 * ee492 + ee534 - ee423 * ee506) * ee4/ee1 -
   ee44) * ee23 - ((ee184 * ee485 + (ee327 + 2 * ee413 -
   ((((ee514 - ee528) * ee110 + 2 * (ee394 * ee54))/ee54 +
   2 * ee394)/ee137 + ee328)) * ee104 + ee362 * ee110)/ee137 +
   ((ee494 + 3 * ee482 + ee532 - ee421 * ee501) * ee3/ee1 -
   ee43) * ee22))/ee137);
out(j, 17) =  - ((((ee481 + 2 * ee479 - (ee261 + 2 * ee466 +
   ee126 + ee450)) * ee4/ee1 + (ee49 - ee423 * ee38)/xi) * ee23 +
   (ee413 * ee89/xi - (ee177 * ee162 + ee401 * ee110 + ee184 * ee517 +
   (2 * ee359 - (((ee513 + ee440) * ee110 + ee512)/(R_pow(ee54, 2) * ee1) +
   2 * ee360)) * ee104 + 2 * (ee394 * ee89/ee136))/ee1)/ee54 -
   ((ee478 + 2 * ee472 - (ee260 + 2 * ee465 +
   ee125 + ee449)) * ee3/ee1 + (ee48 - ee421 * ee36)/xi) * ee22)/ee137);
out(j, 18) = -((((ee353 * ee184 + 2 * (ee471 * ee89))/xi - (ee373 * ee104 +
   ee374 * ee110 + (((ee512 + 8 * (ee218 * ee110/xi)) * ee89 +
   2 * (ee453 * ee110))/ee54 + 2 * (ee384 * ee89))/ee136 +
   2 * (ee177 * ee178))/ee1)/ee54 + ((ee114 * ee149 +
   ee188 + 2 * (ee405 * ee36))/xi + (ee463 - (ee454 + (ee522 -
   ee78)/ee24 + ee426)) * ee3/ee1) * ee22 - ((ee116 * ee150 +
   ee189 + 2 * (ee406 * ee38))/xi + (ee464 - (ee455 + (ee523 -
   ee79)/ee24 + ee427)) *   ee4/ee1) * ee23)/ee137);
out(j, 19) = -(((ee351 * ee110 + (ee223 + ((ee511 * ee110 +
   2 * (ee461 * ee54))/ee54 + 2 * ee461)/ee136 + 2 * ee373 - ee224) * ee89 +
   ee381 * ee177)/ee136 + ee451 - ee452)/ee137);
out(j, 20) =  - ((((3 * (ee84 * ee180) - (ee494 + ee501 * ee531 +
   ee532)) * ee3/ee1 + ee43) * ee22 * ee3 - ((ee483 * ee162 +
   (2 * ee362 - ((ee514 + 4 * ee357 - ee528)/ee54 + 4 * ee162) * ee104/ee137) * ee104 +
   2 * (R_pow(ee162, 2) + ee362 * ee104))/ee137 +
   ((3 * (ee87 * ee181) - (ee495 + ee506 * ee533 +
   ee534)) * ee4/ee1 + ee44) * ee23 * ee4))/ee137);
out(j, 21) = -(((ee362 * ee89/xi - (ee474 * ee162 + (ee320 +
   (4 * ee470 - (ee513 + ee525 + ee440) * ee104/ee137)/ee54 -
   ee322) * ee104 + 2 * (ee178 * ee162 + ee401 * ee104))/ee1)/ee54 +
   ((ee260 + 2 * (ee93 * ee180) + ee125 + ee449 - (ee478 +
   2 * (ee529/ee13))) * ee3/ee1 - (ee36 * ee531 + ee48)/xi) * ee22 * ee3 -
   ((ee261 + 2 * (ee95 * ee181) + ee126 + ee450 -
   (ee481 + 2 * (ee530/ee14))) * ee4/ee1 - (ee38 * ee533 + ee49)/xi) * ee23 * ee4)/ee137);
out(j, 22) =  - ((((ee353 * ee162 + 2 * (ee401 * ee89))/xi -
   ((((ee525 + ee440) * ee89 + 2 * ee453)/ee54 + 4 * ee386) * ee104/ee136 +
   2 * (ee374 * ee104 + R_pow(ee178, 2)))/ee1)/ee54 +
   ((ee388 + (ee78 - ee522)/ee24 + ee303 - (ee454 + ee426)) * ee3/ee1 +
   (ee188 - (ee114 * ee180 + 2 * (ee36 * ee529)))/xi) * ee22 * ee3 -
   ((ee391 + (ee79 - ee523)/ee24 + ee304 -
   (ee455 + ee427)) * ee4/ee1 + (ee189 - (ee116 * ee181 + 2 * (ee38 * ee530)))/xi) * ee23 * ee4)/ee137);
out(j, 23) = -(((ee351 * ee104 + (ee316 + ((ee511 * ee104 +
   2 * (ee460 * ee54))/ee54 + 2 * ee460)/ee136 + 2 * ee374 - ee317) * ee89 +
   ee381 * ee178)/ee136 + ee451 * ee3 - ee452 * ee4)/ee137);
out(j, 24) = ((((ee457 * ee16 + 3 * (ee209 * ee57) - 3 * (ee194 * ee36))/xi -
   ((((24 * ee55 - 6 * ee50)/xi - ee520 * ee3/ee46)/xi -
   ((ee182 + ee107) * ee3/ee46 + ee520/xi) * ee3/ee46)/ee25 +
   ((2 * ee229 - (ee244 + 2 * ee231))/xi + 2 * ee457 +
   2 * ee226) * ee36 + 3 * (ee114 * ee209)))/xi + ((((ee127 +
   2 * ee144 + 6/ee13 + ee140)/xi + ee403 - (ee387 + (2 * ee57 +
   6) * ee90))/xi + 3 * ee132)/xi - ee397) * ee3/ee1) * ee22 +
   (ee161 * ee353 + (((ee425 + 4 * ee255 + ee442)/ee54 + 4 * ee161) * ee89/ee136 +
   2 * ee351) * ee89 + 2 * (ee351 * ee89 +
   R_pow(ee161, 2)))/ee136 - (((ee459 * ee17 + 3 * (ee211 * ee58) -
   3 * (ee197 * ee38))/xi - ((((24 * ee56 - 6 * ee51)/xi -
   ee521 * ee4/ee47)/xi - ((ee183 + ee107) * ee4/ee47 + ee521/xi) * ee4/ee47)/ee26 +
   ((2 * ee230 - (ee245 + 2 * ee232))/xi +
   2 * ee459 + 2 * ee228) * ee38 + 3 * (ee116 * ee211)))/xi +
   ((((ee128 + 2 * ee145 + 6/ee14 + ee141)/xi + ee404 - (ee390 +
   (2 * ee58 + 6) * ee91))/xi + 3 * ee134)/xi - ee399) * ee4/ee1) * ee23)/ee136;
   
} else {

ee1 = exp(lpsi);
ee2 = yh - mu;
ee3 = yl - mu;
ee4 = ee2/ee1;
ee5 = ee3/ee1;
ee8 = exp(-ee4);
ee9 = exp(-ee5);
ee12 = exp(-ee8);
ee13 = exp(-ee9);
ee14 = ee8 * ee12;
ee15 = ee9 * ee13;
ee16 = 1 - ee8;
ee17 = 1 - ee9;
ee18 = ee12 - ee13;
ee20 = ee16 * ee2/ee1;
ee22 = ee17 * ee3/ee1;
ee23 = ee20 - 1;
ee24 = ee22 - 1;
ee26 = ee14 * ee2 - ee15 * ee3;
ee28 = ee23 * ee8 * ee12;
ee30 = ee24 * ee9 * ee13;
ee31 = ee14 - ee15;
ee32 = ee18 * ee1;
ee33 = ee28 * ee2;
ee34 = ee30 * ee3;
ee35 = 3 - ee8;
ee36 = 3 - ee9;
ee38 = ee16 * ee8 * ee12;
ee40 = ee17 * ee9 * ee13;
ee42 = ee35 * ee2/ee1;
ee44 = ee36 * ee3/ee1;
ee45 = ee28 - ee30;
ee46 = ee33 - ee34;
ee47 = (ee42 - 3) * ee8;
ee48 = (ee44 - 3) * ee9;
ee49 = ee38 - ee40;
ee50 = R_pow(ee26, 2);
ee52 = ((ee4 - (ee47 + 3)) * ee2/ee1 + 1) * ee8 * ee12;
ee54 = ((ee5 - (ee48 + 3)) * ee3/ee1 + 1) * ee9 * ee13;
ee56 = (2 - ee8) * ee2/ee1;
ee58 = (2 - ee9) * ee3/ee1;
ee59 = ee50/ee1;
ee60 = ee45 * ee26;
ee61 = (ee56 - 2) * ee8;
ee62 = (ee58 - 2) * ee9;
ee64 = 2 * (ee50/ee32);
ee65 = ee46 * ee18;
ee66 = ee52 * ee2;
ee67 = ee54 * ee3;
ee69 = (ee20 - (ee61 + 2)) * ee8 * ee12;
ee71 = (ee22 - (ee62 + 2)) * ee9 * ee13;
ee73 = (1 - ee35 * ee8) * ee8 * ee12;
ee75 = (1 - ee36 * ee9) * ee9 * ee13;
ee76 = ee31 * ee26;
ee77 = R_pow(ee31, 2);
ee78 = 2 * ee4;
ee79 = 2 * ee5;
ee81 = ee66 - ee67;
ee82 = ee23 * (ee4 - 1);
ee84 = ee33 + ee64 - ee34;
ee85 = ee24 * (ee5 - 1);
ee86 = ee47 + 6;
ee87 = ee48 + 6;
ee88 = ee73 - ee75;
ee89 = 2 * (ee65 - ee59);
ee90 = 8 * ee59;
ee92 = (((ee4 - 6) * ee2/ee1 + 7 - ((2 * (ee4 - 3) + ee78 -  ee86) * ee2/ee1 + 3 * ee82 + 4) * ee8) * ee2/ee1 - 1) *  ee8 * ee12;
ee94 = (((ee5 - 6) * ee3/ee1 + 7 - ((2 * (ee5 - 3) + ee79 -  ee87) * ee3/ee1 + 3 * ee85 + 4) * ee9) * ee3/ee1 - 1) *  ee9 * ee13;
ee95 = ee60 + ee46 * ee31;
ee96 = ee45 * ee18;
ee97 = ee52 - ee54;
ee98 = ee69 - ee71;
ee101 = ee33 + 2 * ee46 + ee64 - ee34;
ee103 = ee49 * ee26/ee1;
ee104 = ee49 * ee18;
ee105 = (ee89 + ee90) * ee31;
ee107 = ee76/ee1;
ee109 = ee18 * R_pow(ee1, 2);
ee110 = ee18 * R_pow(ee1, 3);
ee111 = 2 * (ee76/ee32);
ee112 = 2 * (ee77/ee18);

out(j, 0) = ((ee38 + 2 * ee49 + ee112 - ee40) * ee31/ee18 +
   ee73 - ee75)/ee110;
out(j, 1) = ((ee103 + (2 * ee45 + ee111) * ee31)/ee18 + ee69 -
   ee71)/ee109;
out(j, 2) = 0;
out(j, 3) = ((ee84 * ee31 + 2 * ee60)/ee32 + ee52 - ee54)/ee32;
out(j, 4) = 0;
out(j, 5) = 0;
out(j, 6) = (ee101 * ee26/ee32 + ee66 - ee67)/ee32;
out(j, 7) = 0;
out(j, 8) = 0;
out(j, 9) = 0;
out(j, 10) = (((((2 * (ee104 - ee77) + 4 * ee104 + 8 * ee77)/ee18 +
   4 * ee49) * ee31/ee18 + 2 * ee88) * ee31 + ee49 * (ee38 +
   ee112 - ee40) + 2 * (ee88 * ee31 + R_pow(ee49, 2)))/ee18 +
   (1 - (7 - (6 - ee8) * ee8) * ee8) * ee8 * ee12 - (1 - (7 -
   (6 - ee9) * ee9) * ee9) * ee9 * ee13)/(ee18 * R_pow(ee1, 4));
out(j, 11) = (((ee69 + ((2 * (ee96 - ee107) + 4 * ee96 + 8 * ee107) * ee31/ee18 +
   4 * ee103)/ee18 - ee71) * ee31 + (ee28 +
   ee111 - ee30) * ee49 + ee88 * ee26/ee1 + 2 * (ee98 * ee31 +
   ee45 * ee49))/ee18 + (ee20 - ((ee16 * (ee78 - 1) + 2 * ee56 -
   (ee61 + 8)) * ee8 + 3)) * ee8 * ee12 - (ee22 - ((ee17 * (ee79 -
   1) + 2 * ee58 - (ee62 + 8)) * ee9 + 3)) * ee9 * ee13)/ee110;
out(j, 12) = 0;
out(j, 13) = (((ee84 * ee49 + ((ee105 + 4 * (ee60 * ee18))/ee18 +
   4 * ee60) * ee31/ee18 + 2 * (ee98 * ee26))/ee1 + 2 * (ee97 * ee31 +
   R_pow(ee45, 2)))/ee18 + ((ee4 - ((ee42 - 5) * ee8 +
   5)) * ee2/ee1 + 4 - ((ee78 - ee86) * ee2/ee1 + 2 + 2 * ee82) * ee8) * ee8 * ee12 -
   ((ee5 - ((ee44 - 5) * ee9 + 5)) * ee3/ee1 +
   4 - ((ee79 - ee87) * ee3/ee1 + 2 + 2 * ee85) * ee9) * ee9 * ee13)/ee109;
out(j, 14) = 0;
out(j, 15) = 0;
out(j, 16) = (((((ee105 + 2 * (ee95 * ee18))/ee18 + 2 * ee95)/ee32 +
   ee52 + 2 * ee97 - ee54) * ee26 + ee81 * ee31 + ee45 * ee101)/ee32 +
   ee92 - ee94)/ee32;
out(j, 17) = 0;
out(j, 18) = 0;
out(j, 19) = 0;
out(j, 20) = ((ee46 * ee84 + (((ee89 + 4 * ee65 + ee90)/ee18 +
   4 * ee46) * ee26/ee32 + 2 * ee81) * ee26 + 2 * (ee81 * ee26 +
   R_pow(ee46, 2)))/ee32 + ee92 * ee2 - ee94 * ee3)/ee32;
out(j, 21) = 0;
out(j, 22) = 0;
out(j, 23) = 0;
out(j, 24) = 0;

}

}

return out;

}
