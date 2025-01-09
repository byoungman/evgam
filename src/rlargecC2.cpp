// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

double xi_from_zero_v3(double xi, double eps) 
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

const double xieps3 = 0.0001;

double xi2txi_v3(double xi) 
{
return -log(1.5 / (1.0 + xi) - 1.0);
}

// //' Censored r-largest generalized extreme value (GEV) distribution negative log-likelihood with constrained shape parameter
// //'
// //' @param pars a list of vectors of coefficients for each GEV parameter
// //' @param X1 a design matrix for the GEV location parameter
// //' @param X2 a design matrix for the GEV log scale parameter
// //' @param X3 a design matrix for the GEV transformed shape parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @param drop an integer specifying how many order statistics to consider as censored
// //' @return rlargecd0 a scalar, the negative log-likelihood
// //' @return rlargecd12 a matrix, first then second derivatives w.r.t. GEV parameters
// //' @return rlargecd34 a matrix, third then fourth derivatives w.r.t. GEV parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double rlargecd0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate, int drop)
{
    
arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = ymat.n_rows;
int r = ymat.n_cols - drop;

if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    txivec = txivec.elem(dupid);
}

double mu, lpsi, txi, xi;
arma::rowvec y(r + drop), z(r + drop), w(r + drop), F(r + drop);
double nllh = 0.0;

for (int j=0; j < nobs; j++) {

mu = muvec[j];
lpsi = lpsivec[j];
txi = txivec[j];
xi = 1.5 / (1.0 + exp(-txi)) - 1.0;

xi = xi_from_zero_v3(xi, xieps3);  
txi = xi2txi_v3(xi);

y = ymat.row(j);
z = 1 + xi * (y - mu) / exp(lpsi);

for (int l=0; l < r + drop; l++) {
  w[l] = R_pow(z[l], -1.0/xi);
  if (l < drop) {
    if (z[l] < 0.0) {
      w[l] = 0.0;
    }
  } else {
    if (z[l] < 0.0) {
      nllh = 1e20;
      break;
    }
  }
}

for (int l=0; l < r; l++) {

    if (l == r - 1) {
  
      nllh -= log(exp(-w[l]) - exp(-w[l + drop]));
  
    } else {

      if (l < drop && z[l] <= 0.0) {
      
        nllh += (1 / xi) * log(z[l + drop]);
      
      } else {
        
        nllh -= log(w[l + drop] - w[l]);
      
      }
    
  }

} 

}

return(nllh);

}

// //' @rdname rlargecd0
// [[Rcpp::export]]
arma::mat rlargecd12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate, int drop)
{
    
  arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = ymat.n_rows;
  int r = ymat.n_cols - drop;
  
  if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    txivec = txivec.elem(dupid);
  }
  
  double mu, lpsi, txi, xi;
  arma::rowvec y(r + drop), z(r + drop);
  
  double ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9, ee10;
  double ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
  double ee34, ee35, ee36, ee37, ee38, ee39, ee30, ee31, ee32, ee33;
  double ee43, ee45, ee46, ee47, ee48, ee40, ee41, ee42;
  double ee50, ee56, ee58, ee51, ee53, ee59;
  double ee63, ee64, ee66, ee67, ee68, ee60, ee61, ee65;
  double ee70, ee72, ee73, ee74, ee75, ee76, ee78, ee79;
  double ee80, ee83, ee86;
  
  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    mu = muvec[j];
    lpsi = lpsivec[j];
    txi = txivec[j];
    xi = 1.5 / (1.0 + exp(-txi)) - 1.0;
    
    xi = xi_from_zero_v3(xi, xieps3);  
    txi = xi2txi_v3(xi);
    
    y = ymat.row(j);
    z = 1 + xi * (y - mu) / exp(lpsi);
    
    for (int l=0; l < r; l++) {
      
      if (l == r - 1) {
        
        ee2 = exp(-txi);
        ee3 = 1 + ee2;
        ee5 = 1.5/ee3 - 1;
        ee6 = exp(lpsi);
        ee7 = 1/ee5;
        ee8 = y[l] - mu;
        ee9 = y[l + drop] - mu;
        ee11 = ee5 * ee8/ee6;
        ee13 = ee5 * ee9/ee6;
        ee14 = ee11 + 1;
        ee15 = ee13 + 1;
        ee16 = -ee7;
        ee17 = 1 + ee7;
        ee20 = exp(-R_pow(ee14, ee16));
        ee21 = exp(-R_pow(ee15, ee16));
        ee22 = R_pow(ee14, ee17);
        ee23 = R_pow(ee15, ee17);
        ee24 = ee20 - ee21;
        ee25 = R_pow(ee3, 2);
        ee26 = log1p(ee11);
        ee27 = log1p(ee13);
        ee28 = R_pow(ee14, ee7);
        ee29 = R_pow(ee15, ee7);
        ee34 = 1.5 * (ee8/(ee22 * ee6));
        ee35 = 1.5 * (ee9/(ee23 * ee6));
        ee36 = ee7 + 2;
        ee37 = ee17 * ee5;
        ee38 = ee24 * ee6;
        ee39 = 2 * ee17;
        ee43 = 1.5 * (ee26/(ee28 * ee5)) - ee34;
        ee45 = 1.5 * (ee27/(ee29 * ee5)) - ee35;
        ee46 = R_pow(ee14, ee36);
        ee47 = R_pow(ee15, ee36);
        ee48 = ee3 * ee5;
        ee50 = ee43 * ee20 - ee45 * ee21;
        ee56 = ee21 * ee9/ee23 - ee20 * ee8/ee22;
        ee58 = ee21/ee23 - ee20/ee22;
        ee63 = ee17 * ee8;
        ee64 = ee17 * ee9;
        ee66 = ee25 * ee5 * ee24;
        ee67 = 1.5/ee28;
        ee68 = 1.5/ee29;
        ee70 = 1/R_pow(ee14, ee39) - ee37/ee46;
        ee72 = 1/R_pow(ee15, ee39) - ee37/ee47;
        ee73 = (((ee67 - 1.5) * ee26/ee5 - ee34)/ee5 + 1.5 * (ee63/(ee14 *  ee6))) * ee20;
        ee74 = (((ee68 - 1.5) * ee27/ee5 - ee35)/ee5 + 1.5 * (ee64/(ee15 *  ee6))) * ee21;
        ee75 = (ee70 * ee8/ee6 + 1/ee22) * ee20;
        ee76 = (ee72 * ee9/ee6 + 1/ee23) * ee21;
        ee78 = ee25 * ee24 * ee6;
        ee79 = ee5 * ee24;
        ee80 = R_pow(ee5, 2);
        ee83 = (2.25/ee48 - 3) * ee2/ee3 + 1.5;
        ee86 = (4.5/ee48 - 3) * ee2/ee3 + 1.5;

        out(j, 0) += -(ee58/ee38);
        out(j, 1) += -(ee56/ee38);
        out(j, 2) += ee50 * ee2/ee66;
        out(j, 3) += -((ee70 * ee20 - (ee72 * ee21 + R_pow(ee58, 2)/
          ee24))/(ee24 * R_pow(ee6, 2)));
        out(j, 4) += -((ee75 - (ee76 + ee56 * ee58/ee38))/ee38);
        out(j, 5) += -((ee73/ee22 + ee50 * ee58/ee79 - ee74/ee23) * ee2/
          ee78);
        out(j, 6) += -((ee75 * ee8 - (ee76 * ee9 + R_pow(ee56, 2)/ee38))/
          ee38);
        out(j, 7) += -((ee73 * ee8/ee22 + ee50 * ee56/ee79 - ee74 * ee9/
          ee23) * ee2/ee78);
        out(j, 8) += (((((1.5 - ee67) * ee26/ee5 + ee34) * ee43 * ee2/
          ee25 + (2.25 * (ee2 * ee8/(ee14 * ee25 * ee6)) - ee86 * ee26)/
            ee28)/ee5 + (ee83/ee22 - 1.5 * ((1.5 * (ee26/(ee22 * ee80)) -
              1.5 * (ee63/(ee46 * ee6))) * ee2/ee25)) * ee8/ee6) * ee20 +
              R_pow(ee50, 2) * ee2/ee66 - ((((1.5 - ee68) * ee27/ee5 +
              ee35) * ee45 * ee2/ee25 + (2.25 * (ee2 * ee9/(ee15 * ee25 *
              ee6)) - ee86 * ee27)/ee29)/ee5 + (ee83/ee23 - 1.5 * ((1.5 *
              (ee27/(ee23 * ee80)) - 1.5 * (ee64/(ee47 * ee6))) * ee2/ee25)) *
              ee9/ee6) * ee21) * ee2/ee66;
        
        
      } else {
        
        if (l < drop && z[l] <= 0.0) {
        
          ee2 = exp(-txi);
          ee3 = 1 + ee2;
          ee4 = exp(lpsi);
          ee6 = 1.5/ee3 - 1;
          ee7 = y[l + drop] - mu;
          ee8 = ee6 * ee7;
          ee9 = ee8/ee4;
          ee10 = ee9 + 1;
          ee11 = ee10 * ee4;
          ee12 = R_pow(ee3, 2);
          ee13 = ee8/ee11;
          ee14 = ee10 * ee12;
          ee15 = ee7/ee11;
          ee17 = ee14 * ee6 * ee4;
          ee18 = ee12 * ee6;
          ee20 = 1.5 * ee13 * ee2;
          ee21 = log1p(ee9);
          
          out(j, 0) += -(1/ee11);
          out(j, 1) += -ee15;
          out(j, 2) += (1.5 * ee15 - 1.5 * (ee21/ee6)) * ee2/ee18;
          out(j, 3) += -(ee6/(R_pow(ee10, 2) * R_pow(ee4, 2)));
          out(j, 4) += (1 - ee13)/ee11;
          out(j, 5) += ee20/ee17;
          out(j, 6) += -((ee13 - 1) * ee7/ee11);
          out(j, 7) += ee20 * ee7/ee17;
          out(j, 8) += -(((((2.25 * ee15 + 2.25/ee6)/ee3 - 3) * ee2/ee3 +
            1.5) * ee7/ee11 + (2.25 * (ee2 * ee7/(ee14 * ee4)) - ((4.5/
              (ee3 * ee6) - 3) * ee2/ee3 + 1.5) * ee21)/ee6) * ee2/ee18);
          
        } else {
          
          ee2 = exp(-txi);
          ee3 = 1 + ee2;
          ee5 = 1.5/ee3 - 1;
          ee6 = exp(lpsi);
          ee7 = 1/ee5;
          ee8 = y[l] - mu;
          ee9 = y[l + drop] - mu;
          ee11 = ee5 * ee8/ee6;
          ee13 = ee5 * ee9/ee6;
          ee14 = ee11 + 1;
          ee15 = ee13 + 1;
          ee16 = 1 + ee7;
          ee17 = R_pow(ee14, ee7);
          ee18 = R_pow(ee15, ee7);
          ee19 = R_pow(ee14, ee16);
          ee20 = R_pow(ee15, ee16);
          ee23 = 1/ee18 - 1/ee17;
          ee24 = ee7 + 2;
          ee25 = log1p(ee11);
          ee26 = log1p(ee13);
          ee27 = R_pow(ee3, 2);
          ee28 = ee8/ee19;
          ee29 = ee9/ee20;
          ee30 = R_pow(ee14, ee24);
          ee31 = R_pow(ee15, ee24);
          ee32 = 1/ee19;
          ee33 = 1/ee20;
          ee34 = ee23 * ee6;
          ee35 = ee3 * ee5;
          ee37 = (1.5 * ee28 - 1.5 * ee29)/ee6 + (1.5 * (ee26/ee18) -  1.5 * (ee25/ee17))/ee5;
          ee38 = R_pow(ee5, 2);
          ee39 = ee28 - ee29;
          ee40 = ee32 - ee33;
          ee41 = ee30 * ee6;
          ee42 = ee31 * ee6;
          ee43 = ee16 * ee5;
          ee51 = ee27 * ee5 * ee23;
          ee53 = ee27 * ee23 * ee6;
          ee56 = (2.25/ee35 - 3) * ee2/ee3 + 1.5;
          ee59 = (4.5/ee35 - 3) * ee2/ee3 + 1.5;
          ee60 = ee8/ee30;
          ee61 = ee9/ee31;
          ee65 = 1.5 * (ee25/(ee19 * ee38)) - 1.5 * (ee16 * ee8/ee41);
          ee67 = 1.5 * (ee26/(ee20 * ee38)) - 1.5 * (ee16 * ee9/ee42);
          
          out(j, 0) += ee40/ee34;
          out(j, 1) += ee39/ee34;
          out(j, 2) += -(ee37 * ee2/ee51);
          out(j, 3) += -((ee43 * (1/ee31 - 1/ee30) - R_pow(ee40, 2)/ee23)/
            (ee23 * R_pow(ee6, 2)));
          out(j, 4) += -((((ee61 - ee60) * ee16 * ee5 - ee39 * ee40/ee23)/
            ee6 + ee32 - ee33)/ee34);
          out(j, 5) += (((1.5 * (ee25/ee19) - 1.5 * (ee26/ee20))/ee5 -
            ee37 * ee40/ee23)/ee5 + ee16 * (1.5 * ee61 - 1.5 * ee60)/ee6) *
            ee2/ee53;
          out(j, 6) += (R_pow(ee39, 2)/ee34 + (ee33 - ee43 * ee9/ee42) *
            ee9 - (ee32 - ee43 * ee8/ee41) * ee8)/ee34;
          out(j, 7) += (ee65 * ee8 - (ee37 * ee39/(ee5 * ee23) + ee67 *
            ee9)) * ee2/ee53;
          out(j, 8) += -((((ee56/ee20 - 1.5 * (ee67 * ee2/ee27)) * ee9 -
            (ee56/ee19 - 1.5 * (ee65 * ee2/ee27)) * ee8)/ee6 + (((1.5 *
            ((1.5 * (ee26/(ee18 * ee5)) - 1.5 * (ee9/(ee20 * ee6))) *
            ee26) - 1.5 * ((1.5 * (ee25/(ee17 * ee5)) - 1.5 * (ee8/(ee19 *
            ee6))) * ee25))/ee5 - R_pow(ee37, 2)/ee23) * ee2/ee27 + (2.25 *
            (ee2 * ee9/(ee15 * ee27 * ee6)) - ee59 * ee26)/ee18 -
            (2.25 * (ee2 * ee8/(ee14 * ee27 * ee6)) - ee59 * ee25)/ee17)/
              ee5) * ee2/ee51);
          
        }  
        
      }
      
    }
      
}

return out;

}

// //' @rdname rlargecd0
// [[Rcpp::export]]
arma::mat rlargecd34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat ymat, arma::uvec dupid, int dcate, int drop)
{
    
  arma::vec muvec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lpsivec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = ymat.n_rows;
  int r = ymat.n_cols - drop;
  
  if (dcate == 1) {
    muvec = muvec.elem(dupid);
    lpsivec = lpsivec.elem(dupid);
    txivec = txivec.elem(dupid);
  }
  
double mu, lpsi, txi, xi;
arma::rowvec y(r + drop), z(r + drop);

double ee2, ee3, ee5, ee6, ee7, ee8, ee9, ee10;
double ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
double ee20, ee21, ee22, ee23, ee24, ee27, ee28, ee29, ee25, ee26;
double ee30, ee31, ee32, ee33, ee38, ee39, ee34, ee37;
double ee42, ee43, ee44, ee45, ee46, ee47, ee48, ee49, ee40, ee41;
double ee50, ee51, ee52, ee53, ee54, ee55, ee56, ee57, ee58, ee59;
double ee61, ee63, ee64, ee65, ee60, ee62, ee66, ee67, ee68, ee69;
double ee70, ee71, ee72, ee74, ee76, ee77, ee78, ee79, ee73, ee75;
double ee80, ee82, ee84, ee88, ee81, ee83, ee85, ee86, ee87, ee89;
double ee90, ee91, ee92, ee93, ee94, ee98, ee95, ee96, ee97;
double ee100, ee101, ee102, ee103, ee104, ee105, ee106, ee107, ee108, ee109;
double ee114, ee110, ee111, ee112, ee115, ee117, ee118, ee119;
double ee122, ee124, ee125, ee126, ee127, ee129, ee120, ee121, ee123;
double ee131, ee133, ee134, ee136, ee138, ee130, ee132, ee135, ee137, ee139;
double ee140, ee142, ee143, ee144, ee145, ee141, ee146, ee147, ee148;
double ee150, ee151, ee152, ee155, ee157, ee158, ee159, ee156;
double ee160, ee162, ee164, ee165, ee167, ee168, ee169, ee161, ee163, ee166;
double ee170, ee171, ee172, ee173, ee174, ee175, ee176, ee177, ee178, ee179;
double ee180, ee182, ee184, ee187, ee181, ee185, ee188, ee189;
double ee193, ee195, ee197, ee199, ee190, ee192, ee194;
double ee201, ee202, ee203, ee204, ee205, ee206, ee209, ee200, ee207, ee208;
double ee210, ee211, ee212, ee213, ee215, ee217, ee218, ee219, ee216;
double ee220, ee221, ee222, ee223, ee224, ee225, ee226, ee227, ee228, ee229;
double ee230, ee231, ee232, ee233, ee234, ee235, ee236, ee237, ee238, ee239;
double ee240, ee241, ee242, ee244, ee249, ee243, ee245, ee246, ee247, ee248;
double ee250, ee253, ee256, ee258, ee251, ee252, ee254, ee257, ee259;
double ee260, ee261, ee262, ee267, ee268, ee269, ee263, ee264, ee265, ee266;
double ee270, ee271, ee272, ee274, ee276, ee277, ee278, ee279, ee273, ee275;
double ee280, ee281, ee283, ee286, ee289, ee284, ee285, ee287;
double ee293, ee297, ee298, ee299, ee296;
double ee300, ee301, ee302, ee303, ee304, ee305, ee306, ee307, ee308;
double ee311, ee312, ee310, ee314, ee316, ee317, ee318, ee319;
double ee320, ee322, ee323, ee321, ee325, ee326, ee327, ee328, ee329;
double ee331, ee333, ee337, ee339, ee330, ee334, ee336, ee338;
double ee340, ee341, ee342, ee343, ee344, ee348, ee349;
double ee351, ee352, ee353, ee354, ee355, ee356, ee357, ee358, ee359, ee350;
double ee360, ee361, ee362, ee363, ee365, ee366, ee367, ee368, ee369;
double ee371, ee373, ee375, ee377, ee378, ee379, ee370;
double ee380, ee381, ee382, ee384, ee386, ee389, ee383;
double ee392, ee395, ee390, ee391, ee398, ee399;
double ee401, ee404, ee406, ee409, ee400, ee403, ee405, ee407;
double ee412, ee413, ee414, ee415, ee416, ee417, ee418, ee419, ee411;
double ee420, ee423, ee426, ee427, ee428, ee429, ee421, ee422, ee424, ee425;
double ee431, ee432, ee434, ee435, ee437, ee439, ee433;
double ee440, ee441, ee442, ee443, ee445, ee446, ee447, ee448, ee449;
double ee450, ee451, ee452, ee453, ee454, ee455, ee456, ee457, ee458, ee459;
double ee460, ee461, ee463, ee464, ee465, ee466, ee468, ee469, ee467;
double ee471, ee472, ee477, ee478, ee479, ee470, ee473, ee474, ee476;
double ee480, ee482, ee483, ee485, ee486, ee487, ee488, ee489, ee481, ee484;
double ee490, ee491, ee492, ee493, ee494, ee495, ee496, ee497;
double ee500, ee503, ee505, ee507, ee509;
double ee512, ee513, ee514, ee517, ee518, ee519;
double ee520, ee527, ee528, ee529;
double ee530, ee535, ee536, ee537;
double ee540, ee543, ee545, ee547, ee548, ee549;
double ee550, ee552, ee554, ee556, ee557, ee558, ee559;
double ee560, ee561, ee562, ee563, ee564, ee565, ee566, ee567, ee568, ee569;
double ee570, ee571, ee572, ee577, ee579;
double ee580, ee581, ee583, ee585, ee586, ee587, ee589;
double ee591, ee594, ee595, ee596, ee598;
double ee600, ee601, ee609;
double ee611, ee612, ee613, ee615, ee617, ee618, ee619;
double ee621, ee622, ee624, ee626, ee628, ee629;
double ee632, ee633, ee634, ee635, ee636, ee637, ee638, ee639;
double ee640, ee643, ee644, ee645, ee646, ee647, ee648, ee649;
double ee650, ee651, ee652, ee653, ee654, ee655, ee656, ee657, ee658, ee659;
double ee660;
  
arma::mat out = arma::mat(nobs, 25, arma::fill::zeros);

for (int j=0; j < nobs; j++) {

mu = muvec[j];
lpsi = lpsivec[j];
txi = txivec[j];
xi = 1.5 / (1.0 + exp(-txi)) - 1.0;

xi = xi_from_zero_v3(xi, xieps3);  
txi = xi2txi_v3(xi);

y = ymat.row(j);
z = 1 + xi * (y - mu) / exp(lpsi);

for (int l=0; l < r; l++) {

  if (l == r - 1) {
    
    ee2 = exp(-txi);
    ee3 = 1 + ee2;
    ee5 = 1.5/ee3 - 1;
    ee6 = exp(lpsi);
    ee7 = 1/ee5;
    ee8 = y[l] - mu;
    ee9 = y[l + drop] - mu;
    ee11 = ee5 * ee8/ee6;
    ee13 = ee5 * ee9/ee6;
    ee14 = 1 + ee7;
    ee15 = ee11 + 1;
    ee16 = ee13 + 1;
    ee17 = ee7 + 2;
    ee18 = -ee7;
    ee19 = R_pow(ee3, 2);
    ee20 = R_pow(ee15, ee14);
    ee21 = R_pow(ee16, ee14);
    ee22 = log1p(ee11);
    ee23 = log1p(ee13);
    ee24 = R_pow(ee5, 2);
    ee27 = exp(-R_pow(ee15, ee18));
    ee28 = exp(-R_pow(ee16, ee18));
    ee29 = ee3 * ee5;
    ee30 = R_pow(ee15, ee17);
    ee31 = R_pow(ee16, ee17);
    ee32 = R_pow(ee15, ee7);
    ee33 = R_pow(ee16, ee7);
    ee38 = 1.5 * (ee8/(ee20 * ee6));
    ee39 = 1.5 * (ee9/(ee21 * ee6));
    ee42 = (4.5/ee29 - 3) * ee2/ee3 + 1.5;
    ee43 = ee14 * ee5;
    ee44 = 2 * ee14;
    ee45 = 3 * (ee2/ee3);
    ee46 = 1.5 - ee45;
    ee47 = ee7 + 3;
    ee48 = ee14 * ee8;
    ee49 = ee14 * ee9;
    ee50 = ee32 * ee5;
    ee51 = ee33 * ee5;
    ee52 = 1.5 * (ee22/ee50);
    ee53 = 1.5 * (ee23/ee51);
    ee54 = ee52 - ee38;
    ee55 = ee53 - ee39;
    ee56 = ee30 * ee6;
    ee57 = ee31 * ee6;
    ee58 = 2.25/ee29;
    ee63 = (ee58 - 3) * ee2/ee3 + 1.5;
    ee64 = 1.5 * (ee22/(ee20 * ee24));
    ee65 = 1.5 * (ee23/(ee21 * ee24));
    ee70 = ee64 - 1.5 * (ee48/ee56);
    ee71 = ee65 - 1.5 * (ee49/ee57);
    ee72 = ee27 - ee28;
    ee74 = ee15 * ee19 * ee6;
    ee76 = ee16 * ee19 * ee6;
    ee77 = ee15 * ee6;
    ee78 = ee16 * ee6;
    ee79 = 1.5/ee32;
    ee80 = 1.5/ee33;
    ee82 = ee2 * ee8/ee74;
    ee84 = ee2 * ee9/ee76;
    ee88 = 2.25 * ee82 - ee42 * ee22;
    ee90 = 2.25 * ee84 - ee42 * ee23;
    ee91 = R_pow(ee15, ee44);
    ee92 = R_pow(ee16, ee44);
    ee93 = R_pow(ee15, ee47);
    ee94 = R_pow(ee16, ee47);
    ee98 = 1/ee91 - ee43/ee30;
    ee100 = 1/ee92 - ee43/ee31;
    ee101 = ee70 * ee2;
    ee102 = ee71 * ee2;
    ee103 = ee2/ee19;
    ee104 = ee19 * ee24;
    ee105 = 1 + 2 * ee2;
    ee106 = ee46 * ee5;
    ee107 = 2.25 * (ee2/ee104);
    ee108 = 3 * ee105;
    ee109 = 3 * ee3;
    ee114 = ee14 * ee46;
    ee122 = 1.5 * (ee22/(ee30 * ee24)) - 1.5 * (ee17 * ee8/(ee93 *  ee6));
    ee124 = 1.5 * (ee23/(ee31 * ee24)) - 1.5 * (ee17 * ee9/(ee94 *  ee6));
    ee125 = 1/ee20;
    ee126 = 1/ee21;
    ee127 = 12 * ee2;
    ee129 = ((ee79 - 1.5) * ee22/ee5 - ee38)/ee5 + 1.5 * (ee48/ee77);
    ee131 = ((ee80 - 1.5) * ee23/ee5 - ee39)/ee5 + 1.5 * (ee49/ee78);
    ee133 = ee54 * ee27 - ee55 * ee28;
    ee134 = ee5 * ee17;
    ee136 = ee98 * ee8/ee6;
    ee138 = ee100 * ee9/ee6;
    ee140 = (ee63/ee20 - 1.5 * (ee101/ee19)) * ee8/ee6;
    ee142 = (ee63/ee21 - 1.5 * (ee102/ee19)) * ee9/ee6;
    ee143 = ee88/ee32;
    ee144 = ee90/ee33;
    ee145 = ee114 + ee107;
    ee150 = ee28 * ee9/ee21 - ee27 * ee8/ee20;
    ee151 = 1.5/ee5;
    ee152 = 3 * ee46;
    ee155 = ee28/ee21 - ee27/ee20;
    ee157 = ee151 - 1.5 * ee14;
    ee158 = ee19 * ee5;
    ee159 = ee136 + ee125;
    ee160 = ee138 + ee126;
    ee162 = (((1.5 - ee79) * ee22/ee5 + ee38) * ee54 * ee2/ee19 +  ee143)/ee5 + ee140;
    ee164 = (((1.5 - ee80) * ee23/ee5 + ee39) * ee55 * ee2/ee19 +  ee144)/ee5 + ee142;
    ee165 = 2/ee5;
    ee167 = 2.25 * ee103 - ee106;
    ee168 = ee108 + ee109;
    ee169 = ee129 * ee27;
    ee170 = ee131 * ee28;
    ee171 = ee159 * ee27;
    ee172 = ee160 * ee28;
    ee173 = ee43 * ee8;
    ee174 = ee43 * ee9;
    ee175 = ee165 + 3;
    ee176 = 3 * ee106;
    ee177 = 3 * ee167;
    ee178 = (ee177 - (27 * ee103 + ee176))/ee5;
    ee179 = ee14 * ee122;
    ee180 = ee14 * ee124;
    ee182 = ee134 * ee8/ee77;
    ee184 = ee134 * ee9/ee78;
    ee187 = (((ee178 - ee152)/ee5 + ee108 + ee109 - ee127)/ee3 +  3) * ee2/ee3 - 1.5;
    ee193 = ee88/ee20 + 1.5 * (ee101 * ee22/ee19);
    ee195 = ee90/ee21 + 1.5 * (ee102 * ee23/ee19);
    ee197 = (ee145/ee30 - 1.5 * (ee179 * ee2/ee19)) * ee8/ee6;
    ee199 = (ee145/ee31 - 1.5 * (ee180 * ee2/ee19)) * ee9/ee6;
    ee201 = 1.5 * ee42 + ee152;
    ee202 = 3 * ee42;
    ee203 = ee162 * ee27;
    ee204 = ee164 * ee28;
    ee205 = ee171 * ee8;
    ee206 = ee172 * ee9;
    ee209 = ((ee201/ee5 + ee127 - ee168)/ee3 - 3) * ee2/ee3 +  1.5;
    ee210 = ee14 * (3 - ee182);
    ee211 = ee14 * (3 - ee184);
    ee212 = ee122 * ee5;
    ee213 = ee124 * ee5;
    ee215 = ee169 * ee8/ee20;
    ee217 = ee170 * ee9/ee21;
    ee218 = ee173/ee56;
    ee219 = ee174/ee57;
    ee220 = ee203 - ee204;
    ee221 = ee169/ee20;
    ee222 = ee170/ee21;
    ee223 = ee98 * ee27;
    ee224 = ee100 * ee28;
    ee225 = ee7 + 4;
    ee226 = ee205 - ee206;
    ee227 = ee173/ee77;
    ee228 = ee174/ee78;
    ee229 = ee5 * ee72;
    ee230 = ee72 * ee6;
    ee231 = R_pow(ee15, ee175);
    ee232 = R_pow(ee16, ee175);
    ee233 = ee157/ee30;
    ee234 = ee157/ee31;
    ee235 = ee223 - ee224;
    ee236 = ee215 - ee217;
    ee237 = ee221 - ee222;
    ee238 = ee197 + ee193/ee24;
    ee239 = ee199 + ee195/ee24;
    ee240 = (ee212 + 1.5/ee30) * ee14;
    ee241 = (ee213 + 1.5/ee31) * ee14;
    ee242 = ee171 - ee172;
    ee244 = R_pow(ee133, 2) * ee2;
    ee249 = ee193/ee5;
    ee250 = ee195/ee5;
    ee253 = (2.25 * (ee8/(ee15 * ee3 * ee6)) - 3) * ee2/ee3 +  1.5;
    ee256 = (2.25 * (ee9/(ee16 * ee3 * ee6)) - 3) * ee2/ee3 +  1.5;
    ee258 = ee187 * ee22 + (1.5 * ee253 + ee202) * ee2 * ee8/ee74;
    ee260 = ee187 * ee23 + (1.5 * ee256 + ee202) * ee2 * ee9/ee76;
    ee261 = ee129/ee91;
    ee262 = ee131/ee92;
    ee267 = ee19 * ee72;
    ee268 = ee70/ee20;
    ee269 = ee71/ee21;
    ee270 = 1.5 - ((ee168 - ee127)/ee3 + 3) * ee2/ee3;
    ee271 = ee22/(ee20 * ee5);
    ee272 = ee23/(ee21 * ee5);
    ee274 = ee162/ee20 + 2 * (ee70 * ee54 * ee2/ee19);
    ee276 = ee164/ee21 + 2 * (ee71 * ee55 * ee2/ee19);
    ee277 = 2 * ee227;
    ee278 = 2 * ee228;
    ee279 = ee140 + (ee143 + 1.5 * (ee54 * ee2 * ee22/ee158))/ee5;
    ee280 = ee142 + (ee144 + 1.5 * (ee55 * ee2 * ee23/ee158))/ee5;
    ee281 = ee133 * ee150;
    ee283 = ee46 * ee17 + ee107;
    ee286 = (ee210 - 1) * ee5 * ee8/ee6;
    ee289 = (ee211 - 1) * ee5 * ee9/ee6;
    ee293 = (1 - ee210) * ee5 * ee8/ee6 + 1;
    ee297 = (1 - ee211) * ee5 * ee9/ee6 + 1;
    ee298 = ee218 - ee125;
    ee299 = ee219 - ee126;
    ee300 = ee158 * ee72;
    ee301 = 1.5 * ee271;
    ee302 = 1.5 * ee272;
    ee303 = ((ee274 - ee249)/ee5 - ee197) * ee27;
    ee304 = ((ee276 - ee250)/ee5 - ee199) * ee28;
    ee305 = ee258/ee32;
    ee306 = ee260/ee33;
    ee311 = ee286 - 1;
    ee312 = ee289 - 1;
    ee320 = ((2 * (ee63 * ee70) - 1.5 * ee238) * ee2/ee19 -  ee209/ee20) * ee8/ee6;
    ee322 = ((2 * (ee63 * ee71) - 1.5 * ee239) * ee2/ee19 -  ee209/ee21) * ee9/ee6;
    ee323 = ee267 * ee6;
    ee331 = 1.5 * (ee279 * ee22) + 2 * (ee54 * ee88);
    ee333 = 1.5 * (ee280 * ee23) + 2 * (ee55 * ee90);
    ee337 = 1.5 * (ee22/(ee93 * ee24)) - 1.5 * (ee47 * ee8/(R_pow(ee15, ee225) * ee6));
    ee339 = 1.5 * (ee23/(ee94 * ee24)) - 1.5 * (ee47 * ee9/(R_pow(ee16, ee225) * ee6));
    ee340 = ee125 - ee218;
    ee341 = ee126 - ee219;
    ee342 = 2 - ee182;
    ee343 = 2 - ee184;
    ee344 = ee44 - ee17;
    ee348 = (ee268 + ee233 - (ee261 + ee240)) * ee8/ee6;
    ee349 = ee244/ee158;
    ee351 = (ee269 + ee234 - (ee262 + ee241)) * ee9/ee6;
    ee352 = R_pow(ee150, 2);
    ee353 = ((((4.5 - ee79) * ee22/ee5 + ee38) * ee54 * ee2/ee19 +  3 * ee143)/ee5 + 3 * ee140) * ee54;
    ee354 = ((((4.5 - ee80) * ee23/ee5 + ee39) * ee55 * ee2/ee19 +  3 * ee144)/ee5 + 3 * ee142) * ee55;
    ee355 = ee133 * ee155;
    ee356 = ee331/ee5;
    ee357 = ee333/ee5;
    ee358 = ee122 * ee2;
    ee359 = ee124 * ee2;
    ee360 = ee5 * ee47;
    ee361 = ee157 * ee17;
    ee362 = ee342/ee30;
    ee363 = ee343/ee31;
    ee365 = ee151 - 1.5 * ee17;
    ee366 = 4 * ee2;
    ee367 = 4.5 * ee14;
    ee368 = 4.5/ee5;
    ee369 = ee220 * ee72;
    ee371 = ee311/ee15 - (ee136 + (3 - ee277)/ee20) * ee8/ee6;
    ee373 = ee312/ee16 - (ee138 + (3 - ee278)/ee21) * ee9/ee6;
    ee375 = ((ee356 - ee353) * ee2/ee19 - ee305)/ee5 + ee320;
    ee377 = ((ee357 - ee354) * ee2/ee19 - ee306)/ee5 + ee322;
    ee378 = ee231 * ee6;
    ee379 = ee232 * ee6;
    ee380 = ee54 * ee340;
    ee381 = ee55 * ee341;
    ee382 = ee369 - ee349;
    ee384 = ((ee227 - 2)/ee20 - ee136)/ee20 + (ee362 + ee8/ee378) *  ee14 * ee5;
    ee386 = ((ee228 - 2)/ee21 - ee138)/ee21 + (ee363 + ee9/ee379) *  ee14 * ee5;
    ee389 = (ee283/ee93 - 1.5 * (ee337 * ee17 * ee2/ee19)) *  ee8/ee6 + (ee88/ee30 + 1.5 * (ee358 * ee22/ee19))/ee24;
    ee392 = (ee283/ee94 - 1.5 * (ee339 * ee17 * ee2/ee19)) *  ee9/ee6 + (ee90/ee31 + 1.5 * (ee359 * ee23/ee19))/ee24;
    ee395 = ee281/ee5;
    ee401 = ((ee277 - 3)/ee20 - ee136) * ee8/ee6 - ee293/ee15;
    ee404 = ((ee278 - 3)/ee21 - ee138) * ee9/ee6 - ee297/ee16;
    ee406 = ee14 * ee270 + ee201 * ee2/ee104;
    ee409 = ee43 * (2 - ee50 * ee17)/ee231 - ee98/ee20;
    ee412 = ee43 * (2 - ee51 * ee17)/ee232 - ee100/ee21;
    ee413 = ee54/ee30;
    ee414 = ee55/ee31;
    ee415 = ee150 * ee155;
    ee416 = 2 * (ee244/ee300);
    ee417 = ee303 * ee8;
    ee418 = ee304 * ee9;
    ee419 = ee384 * ee27;
    ee420 = ee386 * ee28;
    ee423 = ee371 * ee27 * ee8/ee20;
    ee426 = ee373 * ee28 * ee9/ee21;
    ee427 = ee375 * ee27;
    ee428 = ee377 * ee28;
    ee429 = ee236 * ee72;
    ee431 = (ee348 + (ee301 - ee380)/ee5) * ee27 * ee8;
    ee432 = ((ee413 - ee212) * ee14 + ee268 + ee233 - ee261) *  ee27;
    ee434 = (ee351 + (ee302 - ee381)/ee5) * ee28 * ee9;
    ee435 = ((ee414 - ee213) * ee14 + ee269 + ee234 - ee262) *  ee28;
    ee437 = ee401 * ee27/ee20;
    ee439 = ee404 * ee28/ee21;
    ee440 = ee409 * ee27;
    ee441 = ee412 * ee28;
    ee442 = ee298 * ee54;
    ee443 = ee299 * ee55;
    ee445 = R_pow(ee15, ee344);
    ee446 = R_pow(ee16, ee344);
    ee447 = ee352/ee6;
    ee448 = 2 * (ee281/ee229);
    ee449 = R_pow(ee6, 2);
    ee450 = ((ee442 + ee301)/ee5 + ee348) * ee27;
    ee451 = ((ee443 + ee302)/ee5 + ee351) * ee28;
    ee452 = ee423 - ee426;
    ee453 = ee427 - ee428;
    ee454 = ee236 * ee155;
    ee455 = ee237 * ee72;
    ee456 = ee129/ee445;
    ee457 = ee131/ee446;
    ee458 = ee242 * ee150;
    ee459 = ee226 * ee72;
    ee460 = ee226 * ee155;
    ee461 = ee440 - ee441;
    ee463 = ee355/ee5;
    ee464 = ee15 * ee70;
    ee465 = ee16 * ee71;
    ee466 = ee206 + 2 * (ee352/ee230);
    ee468 = ((ee58 + 3) * ee2/ee3 - 1.5)/ee5 + 2 * ee114;
    ee469 = ((6 * (2 * ee105 + ee366) + 8 * ee168 - 96 * ee2)/ee3 +  12) * ee2;
    ee471 = ee337 * ee5;
    ee472 = ee339 * ee5;
    ee477 = ee361/ee93;
    ee478 = ee361/ee94;
    ee479 = ee365/ee93;
    ee480 = ee365/ee94;
    ee482 = ee178 + 1.5 * ee46 + 4.5 * ee42;
    ee483 = R_pow(ee155, 2);
    ee485 = 2 * (ee355/ee229);
    ee486 = ee44 - ee175;
    ee487 = 3 - ee360 * ee8/ee77;
    ee488 = 3 - ee360 * ee9/ee78;
    ee489 = 3 * (1 + ee366);
    ee490 = 6 * ee3;
    ee491 = 9 * ee105;
    ee492 = ee303 - ee304;
    ee493 = ee417 - ee418;
    ee494 = ee258/ee20;
    ee495 = ee260/ee21;
    ee496 = ee220 * ee150;
    ee497 = ee220 * ee155;
    ee500 = ee419 - ee420;
    ee503 = ee203 + 2 * ee220 + ee416 - ee204;
    ee505 = ee203 + ee416 - ee204;
    ee507 = ee389 * ee5;
    ee509 = ee392 * ee5;
    ee512 = ee236 * ee133 * ee2/ee19;
    ee513 = ee429 - ee395;
    ee514 = ee454 + ee237 * ee150;
    ee517 = ee237 * ee133 * ee2/ee19;
    ee518 = (((ee479 - ee471 * ee17) * ee14 + ee477) * ee8/ee6 +  2 * ee179) * ee5;
    ee519 = (((ee480 - ee472 * ee17) * ee14 + ee478) * ee9/ee6 +  2 * ee180) * ee5;
    ee520 = ee431 - ee434;
    ee527 = ee242 * ee72;
    ee528 = ee458 + ee460;
    ee529 = ee459 + ee447;
    ee530 = ee437 - ee439;
    ee535 = ee133 * ee235;
    ee536 = ee159 * ee70;
    ee537 = ee160 * ee71;
    ee540 = (ee17 * ee487 + ee7) * ee5 * ee8/ee6;
    ee543 = (ee17 * ee488 + ee7) * ee5 * ee9/ee6;
    ee545 = ((2 * (ee145 * ee122) - 1.5 * (ee389 * ee14)) *  ee2/ee19 - ee406/ee30) * ee8/ee6;
    ee547 = ((2 * (ee145 * ee124) - 1.5 * (ee392 * ee14)) *  ee2/ee19 - ee406/ee31) * ee9/ee6;
    ee548 = ((2 * ee249 - ee274)/ee5 + 2 * ee197 - 2 * (ee129 *  ee70 * ee2/ee19))/ee20;
    ee549 = ((2 * ee250 - ee276)/ee5 + 2 * ee199 - 2 * (ee131 *  ee71 * ee2/ee19))/ee21;
    ee550 = ee267 * ee449;
    ee552 = ee122 * ee157;
    ee554 = ee124 * ee157;
    ee556 = ee415/ee6;
    ee557 = 1.5 * (ee238 * ee22);
    ee558 = 1.5 * (ee239 * ee23);
    ee559 = 2 * ee382;
    ee560 = 2 * (ee70 * ee88);
    ee561 = 2 * (ee71 * ee90);
    ee562 = 2 * (ee415/ee230);
    ee563 = 4 * ee455;
    ee564 = 8 * ee395;
    ee565 = 8 * ee463;
    ee566 = 8 * ee349;
    ee567 = R_pow(ee6, 3);
    ee568 = ee382 * ee150;
    ee569 = ((ee375/ee20 + (3 * (ee238 * ee54) + 3 * (ee162 *  ee70)) * ee2/ee19 - ((ee482/ee20 + ee557 + ee560) * ee2/ee19 -  ee494)/ee5)/ee5 - ee545) * ee27;
    ee570 = ((ee377/ee21 + (3 * (ee239 * ee55) + 3 * (ee164 *  ee71)) * ee2/ee19 - ((ee482/ee21 + ee558 + ee561) * ee2/ee19 -  ee495)/ee5)/ee5 - ee547) * ee28;
    ee571 = ee496 + ee512;
    ee572 = ee497 + ee517;
    ee577 = ((((ee46 * (6 * ee167 - 18 * ee103) + (12 * (ee177 -  ee176) + 9 * (2 * ee167 + 9 * ee103) - 324 * ee103) *  ee2/ee158 + 3 * (4.5 * (ee46 * ee2/ee19) - ee270 * ee5) -  6 * (ee270 * ee5))/ee5 - 3 * ee270)/ee5 + ee489 + ee490 +  ee491 - ee469)/ee3 + 3) * ee2/ee3 - 1.5;
    ee579 = ee450 - ee451;
    ee580 = ee371/R_pow(ee15, ee486);
    ee581 = ee373/R_pow(ee16, ee486);
    ee583 = ((ee464 + ee368 - (ee456 + ee367))/ee30 - 3 * ee240) *  ee8/ee6;
    ee585 = ((ee465 + ee368 - (ee457 + ee367))/ee31 - 3 * ee241) *  ee9/ee6;
    ee586 = ee518 + ee240;
    ee587 = ee519 + ee241;
    ee589 = (ee468/ee30 + ee548 + 2 * (ee552 * ee2/ee19) - (ee507 +  1.5 * (ee358/ee19)) * ee14) * ee8/ee6;
    ee591 = (ee468/ee31 + ee549 + 2 * (ee554 * ee2/ee19) - (ee509 +  1.5 * (ee359/ee19)) * ee14) * ee9/ee6;
    ee594 = (((1.5 * ee187 - (3 * (ee42 * ee46) + 4.5 * ee270))/ee5 +  ee489 + ee490 + ee491 - ee469)/ee3 + 3) * ee2/ee3 - 1.5;
    ee595 = ee432 - ee435;
    ee596 = ee129 * ee298;
    ee598 = ee215 + ee448 - ee217;
    ee600 = ee221 + ee485 - ee222;
    ee601 = ee131 * ee299;
    ee609 = ee242 * ee133/ee5;
    ee611 = ee226 * ee133/ee5;
    ee612 = ee298 * ee159;
    ee613 = ee299 * ee160;
    ee615 = (ee331 * ee2/ee158 - ee305)/ee5 + ee320;
    ee617 = (ee333 * ee2/ee158 - ee306)/ee5 + ee322;
    ee618 = ee535/ee5;
    ee619 = ee205 - ee466;
    ee621 = ee205 + 2 * ee226 - ee466;
    ee622 = ee235 * ee72;
    ee624 = ee235 * ee150/ee6;
    ee626 = ((2 - ee540)/ee93 + ee362) * ee14 * ee5;
    ee628 = ((2 - ee543)/ee94 + ee363) * ee14 * ee5;
    ee629 = ee114 + (ee45 - 1.5)/ee5;
    ee632 = ee233 - ee179 * ee5;
    ee633 = ee234 - ee180 * ee5;
    ee634 = ee224 + 2 * (ee483/ee72);
    ee635 = ee72 * ee449;
    ee636 = ee72 * ee567;
    ee637 = R_pow(ee72, 2);
    ee638 = ee559 + ee566;
    ee639 = 2 * (ee514 * ee72);
    ee640 = 2 * ee513;
    ee643 = 2 * ee215 + ee448 - 2 * ee217;
    ee644 = 2 * (((ee233 - ee240) * ee8/ee6 + ee64) * ee54 *  ee2/ee19);
    ee645 = 2 * (((ee234 - ee241) * ee9/ee6 + ee65) * ee55 *  ee2/ee19);
    ee646 = 2 * ee529;
    ee647 = 2 * ee157;
    ee648 = ee165 + 4;
    ee649 = 3 * ee14;
    ee650 = 3 * ee271;
    ee651 = 3 * ee272;
    ee652 = 3/ee5;
    ee653 = 4 * ee429;
    ee654 = ee563 + ee565;
    ee655 = 4 * ee527;
    ee656 = 4.5 * ee187;
    ee657 = 4.5 * ee46;
    ee658 = ee368 - ee367;
    ee659 = 8 * ee556;
    ee660 = 8 * ee447;
    
    out(j, 0) += -((ee440 - (ee441 + (ee223 + 2 * ee235 - ee634) *
      ee155/ee72))/ee636);
    out(j, 1) += -((ee419 - (ee420 + (ee624 + (2 * ee242 - ee562) *
      ee155)/ee72))/ee635);
    out(j, 2) += -(((ee618 - (2 * ee237 + ee485) * ee155)/ee72 +
      ee432 - ee435) * ee2/ee550);
    out(j, 3) += -((ee437 - ((ee460 + (2 * ee171 - (2 * ee172 + ee562)) *
      ee150)/ee230 + ee439))/ee230);
    out(j, 4) += -((ee450 + (ee609 - (ee454 + ee600 * ee150)/ee6)/
      ee72 - ee451) * ee2/ee323);
    out(j, 5) += -(((ee497 + ee133 * (2 * ee221 + ee485 - 2 * ee222) *
      ee2/ee19)/ee229 + ee303 - ee304) * ee2/ee323);
    out(j, 6) += -((ee423 - (ee426 + ee621 * ee150/ee230))/ee230);
    out(j, 7) += -(((ee611 - (2 * ee236 + ee448) * ee150/ee6)/ee72 +
      ee431 - ee434) * ee2/ee323);
    out(j, 8) += -(((ee496 + ee133 * ee643 * ee2/ee19)/ee229 + ee417 -
      ee418) * ee2/ee323);
    out(j, 9) += (ee503 * ee133 * ee2/ee300 + ee427 - ee428) * ee2/
      ee300;
    out(j, 10) += -(((ee412/ee21 + ((1/ee446 - (ee17 * (2 - ee51 *
      ee47) + ee165 + 3) * ee5)/ee31 + 2 * ee100) * ee14 * ee5/
        ee31) * ee28 - ((ee409/ee20 + ((1/ee445 - (ee17 * (2 - ee50 *
          ee47) + ee165 + 3) * ee5)/ee30 + 2 * ee98) * ee14 * ee5/ee30) *
          ee27 + ((ee223 - ee634) * ee235 + (2 * ee461 - ((2 * (ee622 +
          ee483) + 4 * ee622 - 8 * ee483)/ee72 + 4 * ee235) * ee155/
            ee72) * ee155 + 2 * (ee461 * ee155 + R_pow(ee235, 2)))/
              ee72))/(ee72 * R_pow(ee6, 4)));
    out(j, 11) += -(((ee386/ee21 + ee299 * ee100 + ((2 * ee160 +
      2 * (ee343/ee21))/ee31 - (ee488/ee94 + ee9/(R_pow(ee16, ee648) *
      ee6)) * ee5 * ee17) * ee14 * ee5) * ee28 - (((ee419 - (ee420 +
      ((2 * (ee527 + ee556) + ee655 - ee659) * ee155/ee72 +
      4 * ee624)/ee72)) * ee155 + ee461 * ee150/ee6 + (ee171 - (ee172 +
      ee562)) * ee235 + 2 * (ee500 * ee155 + ee242 * ee235))/
        ee72 + (ee384/ee20 + ee298 * ee98 + ((2 * ee159 + 2 * (ee342/
          ee20))/ee30 - (ee487/ee93 + ee8/(R_pow(ee15, ee648) * ee6)) *
            ee5 * ee17) * ee14 * ee5) * 
            ee27))/ee636);
    out(j, 12) += -((((((ee131 * ee5 - 4.5) * ee14 + ee465 + ee368 -
      ee457)/ee31 + (ee414 - 3 * ee213) * ee14)/ee21 + ((ee131/
        ee232 - ((ee55/ee94 - ee472) * ee17 + ee480)) * ee14 - ee478) *
          ee5 + ee71 * ee100) * ee28 + (ee461 * ee133/ee5 - ((ee432 +
          (4 * ee618 - (2 * (ee455 - ee463) + ee563 + ee565) * ee155/
            ee72)/ee72 - ee435) * ee155 + ee600 * ee235 + 2 * (ee595 *
              ee155 + ee237 * ee235)))/ee72 - ((((ee129 * ee5 - 4.5) *
              ee14 + ee464 + ee368 - ee456)/ee30 + (ee413 - 3 * ee212) * ee14)/
                ee20 + ((ee129/ee231 - ((ee54/ee93 - ee471) * ee17 + ee479)) *
                  ee14 - ee477) * ee5 + ee70 * ee98) * ee27) * ee2/(ee267 *
                  ee567));
    out(j, 13) += -(((((ee21 * ee160 + 4 - 2 * ee184) * ee9/ee379 -
      ((ee543 - 2)/ee94 + (ee184 - 2)/ee31)) * ee14 * ee5 + ee404/
        ee92 + 2 * ee613 - ee297/ee232) * ee28 - ((((ee20 * ee159 +
          4 - 2 * ee182) * ee8/ee378 - ((ee540 - 2)/ee93 + (ee182 -
          2)/ee30)) * ee14 * ee5 + ee401/ee91 + 2 * ee612 - ee293/ee231) *
          ee27 + ((ee619 * ee235 + 2 * (ee500 * ee150) - (((ee655 -
          ee659) * ee150 + 2 * (ee529 * ee155))/ee72 + 4 * ee458) *
          ee155/ee72)/ee6 + 2 * (ee530 * ee155 + R_pow(ee242, 2)))/
            ee72))/ee635);
    out(j, 14) += -(((ee500 * ee133/ee5 - ((ee595 * ee150 + ee598 *
      ee235)/ee6 + (4 * ee609 - (ee654 * ee150 + 2 * (ee513 * ee155))/
        ee230) * ee155/ee72 + 2 * (ee579 * ee155 + ee237 * ee242)))/
          ee72 + ((((ee465 + ee652 - (ee457 + ee649))/ee31 - 2 *
            ee241) * ee9/ee6 + (ee443 + ee651)/ee5 + ee601)/ee21 + (ee131 *
            ee14 * ee5/ee232 + ee633/ee21) * ee9/ee6 + (ee14 * ee55 *
            ee343 + ee647)/ee31 + ee537 - ee519) * ee28 - ((((ee464 +
            ee652 - (ee456 + ee649))/ee30 - 2 * ee240) * ee8/ee6 + (ee442 +
            ee650)/ee5 + ee596)/ee20 + (ee129 * ee14 * ee5/ee231 +
            ee632/ee20) * ee8/ee6 + (ee14 * ee54 * ee342 + ee647)/ee30 +
            ee536 - ee518) * ee27) * ee2/ee550);
    out(j, 15) += -((((ee220 * ee235 + ee133 * (2 * ee432 + 2 * (ee535/
      ee229) - 2 * ee435) * ee2/ee19 - ((ee133 * ee654 * ee2/
        ee19 + 2 * (ee382 * ee155))/ee72 + 4 * ee517) * ee155/ee72)/
          ee5 - 2 * (ee492 * ee155 + R_pow(ee237, 2) * ee2/ee19))/ee72 +
            (ee629/ee30 + ee548 - ((ee507 - ee162/ee30) * ee14 + (2 *
            (ee632 * ee54/ee5) - 2 * ee552) * ee2/ee19)) * ee27 - (ee629/
              ee31 + ee549 - ((ee509 - ee164/ee31) * ee14 + (2 * (ee633 *
                ee55/ee5) - 2 * ee554) * ee2/ee19)) * ee28) * ee2/ee550);
    out(j, 16) += -(((((ee581 - (ee297 + 2 * ee297))/ee232 + ee628 +
      3 * ee613) * ee9/ee6 + ee312/ee31) * ee28 - ((((ee580 -
      (ee293 + 2 * ee293))/ee231 + ee626 + 3 * ee612) * ee8/ee6 +
      ee311/ee30) * ee27 + (ee452 * ee155 + (ee437 + 2 * ee530 - ((((ee646 -
      ee660) * ee155 + 2 * (ee528 * ee72))/ee72 + 2 * ee528)/
        ee230 + ee439)) * ee150 + ee242 * ee621)/ee230))/ee230);
    out(j, 17) = -(((((ee585 + ((ee53 + 4.5)/ee5 - ee367)/ee16 +
      (ee651 - ee381)/ee5 + 2 * ee601)/ee21 + ee537 - ee587) * ee9/
        ee6 + (ee302 - ee297 * ee55/ee31)/ee5) * ee28 + (ee530 * ee133/
          ee5 - (ee520 * ee155 + ee237 * ee226 + ee242 * ee643 +
            (2 * ee450 - (((ee640 + ee564) * ee155 + ee639)/(ee637 *   ee6) +
            2 * ee451)) * ee150 + 2 * (ee528 * ee133/ee229))/ee6)/
              ee72 - (((ee583 + ((ee52 + 4.5)/ee5 - ee367)/ee15 + (ee650 -
                ee380)/ee5 + 2 * ee596)/ee20 + ee536 - ee586) * ee8/ee6 + (ee301 -
                ee293 * ee54/ee30)/ee5) * ee27) * ee2/ee323);
    out(j, 18) += -((((ee505 * ee242 + 2 * (ee579 * ee133 * ee2/
      ee19))/ee5 - (ee492 * ee150 + ee493 * ee155 + (ee133 * (ee639 +
        8 * (ee281 * ee155/ee5)) * ee2/ee19 + 2 * (ee568 * ee155))/
          (ee5 * ee637) + (2 * (ee514 * ee133/ee229) + 2 * (ee236 *
            ee237)) * ee2/ee19)/ee6)/ee72 + (ee589 - (ee644 - (ee162 * ee298 +
            ee249))/ee5) * ee27 - (ee591 - (ee645 - (ee164 * ee299 +
            ee250))/ee5) * ee28) * ee2/ee323);
    out(j, 19) += -(((((ee303 + ((ee638 * ee155 + 2 * (ee572 * ee72))/
      ee72 + 2 * ee572)/ee229 + 2 * ee492 - ee304) * ee133 +
        ee503 * ee237) * ee2/ee19 + ee453 * ee155)/ee229 + ee569 - ee570) *
        ee2/ee323);
    out(j, 20) += -(((ee293/ee30 + (3 * (ee159 * ee340) - ((ee580 +
      ee286 + 2 * ee311 - 1)/ee231 + ee626)) * ee8/ee6) * ee27 *
      ee8 - ((ee619 * ee226 + (2 * ee452 - ((ee646 + 4 * ee459 -
      ee660)/ee72 + 4 * ee226) * ee150/ee230) * ee150 + 2 * (ee452 *
      ee150 + R_pow(ee226, 2)))/ee230 + (ee297/ee31 + (3 * (ee160 *
      ee341) - ((ee581 + ee289 + 2 * ee312 - 1)/ee232 + ee628)) *
      ee9/ee6) * ee28 * ee9))/ee230);
    out(j, 21) += -(((ee452 * ee133/ee5 - ((ee431 + (4 * ee611 -
      (ee640 + ee653 + ee564) * ee150/ee230)/ee72 - ee434) * ee150 +
      ee598 * ee226 + 2 * (ee520 * ee150 + ee236 * ee226))/ee6)/
        ee72 + ((ee586 + (2 * (ee129 * ee340) - (ee583 + (4.5 * ee271 -
          ee380)/ee5 + ee658/ee15))/ee20 - ee536) * ee8/ee6 - (ee311 *
          ee54/ee30 + ee301)/ee5) * ee27 * ee8 - ((ee587 + (2 *
          (ee131 * ee341) - (ee585 + (4.5 * ee272 - ee381)/ee5 + ee658/
            ee16))/ee21 - ee537) * ee9/ee6 - (ee312 * ee55/ee31 + ee302)/
              ee5) * ee28 * ee9) * ee2/ee323);
    out(j, 22) += -((((ee505 * ee226 + 2 * (ee520 * ee133 * ee2/
      ee19))/ee5 - (((ee133 * (ee653 + ee564) * ee2/ee19 + 2 * ee568)/
        ee72 + 4 * ee512) * ee150/ee229 + 2 * (ee493 * ee150 + R_pow(ee236, 2) *
          ee2/ee19))/ee6)/ee72 + (ee589 + (ee249 - (ee162 *
          ee340 + ee644))/ee5) * ee27 * ee8 - (ee591 + (ee250 -
          (ee164 * ee341 + ee645))/ee5) * ee28 * ee9) * ee2/ee323);
    out(j, 23) += -(((((ee417 + ((ee638 * ee150 + 2 * (ee571 * ee72))/
      ee72 + 2 * ee571)/ee229 + 2 * ee493 - ee418) * ee133 +
        ee503 * ee236) * ee2/ee19 + ee453 * ee150)/ee229 + ee569 * ee8 -
        ee570 * ee9) * ee2/ee323);
    out(j, 24) += ((ee220 * ee505 + (ee133 * ((ee559 + 4 * ee369 +
      ee566)/ee72 + 4 * ee220) * ee2/ee300 + 2 * ee453) * ee133 +
      2 * (R_pow(ee220, 2) + ee453 * ee133)) * ee2/ee300 + ((((1.5 *
      (ee615 * ee22) + 3 * (ee279 * ee88) - 3 * (ee258 * ee54))/
        ee5 - ((((2 * ee356 - ee353) * ee2/ee19 - 2 * ee305)/ee5 +
          2 * ee615 + 2 * ee320) * ee54 + 3 * (ee162 * ee279))) * ee2/
            ee19 - ((ee656 - (1.5 * ((((ee657 + 6.75 * ee82) * ee8/ee77 +
              ee127 - ee168)/ee3 - 3) * ee2/ee3 + 1.5) + 3 * (ee253 *
              ee42))) * ee2 * ee8/ee74 - ee577 * 
              ee22)/ee32)/ee5 + ((3 *
              (ee238 * ee63) - (1.5 * (((ee557 + ee560) * ee2/ee19 - ee494)/
                ee24 + ee545) + 3 * (ee209 * ee70))) * ee2/ee19 - ee594/
                  ee20) * ee8/ee6) * ee27 - ((((1.5 * (ee617 * ee23) + 3 * (ee280 *
                    ee90) - 3 * (ee260 * ee55))/ee5 - ((((2 * ee357 - ee354) *
                    ee2/ee19 - 2 * ee306)/ee5 + 2 * ee617 + 2 * ee322) *
                    ee55 + 3 * (ee164 * ee280))) * ee2/ee19 - ((ee656 - (1.5 * ((((ee657 +
                    6.75 * ee84) * ee9/ee78 + ee127 - ee168)/ee3 - 3) *
                    ee2/ee3 + 1.5) + 3 * (ee256 * ee42))) * ee2 * ee9/ee76 -
                    ee577 * ee23)/ee33)/ee5 + ((3 * (ee239 * ee63) - (1.5 *
                    (((ee558 + ee561) * ee2/ee19 - ee495)/ee24 + ee547) + 3 * (ee209 *
                    ee71))) * ee2/ee19 - ee594/ee21) * ee9/ee6) * ee28) *
                    ee2/ee300;
    
    
  } else {
    
    if (l < drop && z[l] <= 0.0) {
      
      ee2 = exp(-txi);
      ee3 = 1 + ee2;
      ee5 = 1.5/ee3 - 1;
      ee6 = exp(lpsi);
      ee7 = y[l + drop] - mu;
      ee8 = ee5 * ee7;
      ee9 = ee8/ee6;
      ee10 = ee9 + 1;
      ee11 = R_pow(ee3, 2);
      ee12 = ee10 * ee6;
      ee13 = 1.5 - 3 * (ee2/ee3);
      ee14 = ee2/ee11;
      ee15 = ee13 * ee5;
      ee16 = ee8/ee12;
      ee17 = 2.25 * ee14;
      ee18 = 1 + 2 * ee2;
      ee19 = ee10 * ee11;
      ee20 = ee2 * ee7;
      ee21 = ee17 - ee15;
      ee22 = 3 * ee18;
      ee23 = 3 * ee3;
      ee26 = ee7/(ee10 * ee3 * ee6);
      ee27 = 12 * ee2;
      ee29 = ee22 + ee23;
      ee30 = ee20/(ee11 * ee6);
      ee32 = ee19 * ee5 * ee6;
      ee33 = ee3 * ee5;
      ee34 = ee19 * ee6;
      ee37 = (4.5/ee33 - 3) * ee2/ee3 + 1.5;
      ee38 = 3 * ee21;
      ee39 = ee10 * ee13;
      ee43 = (4.5 * ee26 - 3) * ee2/ee3;
      ee44 = 1.5 - ((ee29 - ee27)/ee3 + 3) * ee2/ee3;
      ee45 = 3 * ee16;
      ee46 = 4.5 * ee13;
      ee48 = (2.25 * ee26 - 3) * ee2/ee3;
      ee49 = 2 * ee16;
      ee51 = 2.25 * ee30 - ee39;
      ee52 = 3 * ee15;
      ee53 = R_pow(ee10, 2);
      ee54 = (ee38 - (27 * ee14 + ee52))/ee5;
      ee55 = ee43 + 1.5;
      ee57 = (ee46 + 6.75 * (ee20/ee34)) * ee7/ee12;
      ee58 = 3 * ee13;
      ee59 = R_pow(ee6, 2);
      ee61 = ee44 * ee5;
      ee64 = ee5 * (4.5 - ee45) * ee7/ee12;
      ee65 = (2 * (1.5 * ee10 + 1.5 * ee9) + 6 * ee10 - 12 * ee9)/ee10;
      ee66 = ee48 + 1.5;
      ee68 = (ee17 - ee55 * ee5) * ee7/ee12;
      ee71 = 1.5 - 1.5 * ee16;
      ee73 = 2 * (1 + 2 * ee9) + 4 * ee10;
      ee74 = 2 * ee21;
      ee75 = 2.25 * (ee7/ee12);
      ee76 = 4 * ee2;
      ee77 = 8 * ee9;
      ee78 = 9 * ee14;
      ee81 = (((ee54 - ee58)/ee5 + ee22 + ee23 - ee27)/ee3 + 3) *  ee2/ee3 - 1.5;
      ee83 = (ee5 * (2 * ee51 - 18 * ee30) + 9 * (ee10 * ee2/ee11))/ee10 +  ee78;
      ee85 = ee53 * ee11 * ee59;
      ee86 = ee11 * ee5;
      ee87 = (ee77 - ee73)/ee10;
      ee88 = R_pow(ee6, 3);
      ee89 = (((((ee5 * (3 * ee51 - 27 * ee30) + 3 * (ee10 * ee21))/ee10 +  ee38) * ee7/ee12 + ee46) * ee2/ee11 - ee61) * ee7/ee12 +  ((ee54 + 3 * (ee37 * ee71) + 4.5 * (ee48 + ee68 + 1.5) -  ee58)/ee5 + ee57 + 0)/ee3 * ee2/ee3 + 0) * ee2;
      ee91 = (((3 * ee71 - 4.5)/ee5 + ee75)/ee3 * ee2/ee3 + ee68 +  0) * ee2;
      ee94 = ((ee57 + ee27 - ee29)/ee3 - 3) * ee2/ee3 + 1.5;
      ee95 = ee53 * ee59;
      ee96 = R_pow(ee10, 3);
      ee97 = R_pow(ee10, 4);
      ee100 = (ee65 + 9) * ee5 * ee7/ee12;
      ee102 = ((6 * (2 * ee18 + ee76) + 8 * ee29 - 96 * ee2)/ee3 +  12) * ee2;
      ee106 = ee5 * (ee49 - 3) * ee7/ee12 + 1;
      ee110 = ee5 * (3 - ee49) * ee7/ee12 - 1;
      ee111 = ee16 - 1;
      ee112 = R_pow(ee5, 2);
      ee114 = (ee74 + ee17 - (ee83 * ee7/ee12 + ee43 + 1.5) *  ee5) * ee7/ee12;
      ee115 = 1 - ee16;
      ee117 = 2 - ee49;
      ee118 = 3 - ee45;
      ee119 = 3 * (ee66 * ee37);
      ee120 = 3 * (1 + ee76);
      ee121 = 3 * ee44;
      ee122 = 6 * ee3;
      ee123 = 9 * ee18;
      ee124 = log1p(ee9);
      
      out(j, 0) += -(2 * (ee112/(ee96 * ee88)));
      out(j, 1) += ee5 * ee117/ee95;
      out(j, 2) += -((1.5 - ee45) * ee2/ee85);
      out(j, 3) += ee110/ee12;
      out(j, 4) += (1.5 - (ee64 + 1.5 * ee115)) * ee2/ee32;
      out(j, 5) += ee91/ee32;
      out(j, 6) += -(ee106 * ee7/ee12);
      out(j, 7) += -((ee64 - (1.5 + 1.5 * ee111)) * ee2 * ee7/ee32);
      out(j, 8) += ee91 * ee7/ee32;
      out(j, 9) += ((ee81 * ee124 + (1.5 * ee66 + 3 * ee37) * ee2 *
        ee7/ee34)/ee5 + ((((1.5 * ee37 + 3 * ee66)/ee5 + ee57 + ee27 -
        ee29)/ee3 - 3) * ee2/ee3 + 1.5) * ee7/ee12) * ee2/ee86;
      out(j, 10) += -(6 * (R_pow(ee5, 3)/(ee97 * R_pow(ee6, 4))));
      out(j, 11) += ee112 * (ee73 - ee77)/(ee97 * ee88);
      out(j, 12) += -((ee65 - 3) * ee5 * ee2/(ee96 * ee11 * ee88));
      out(j, 13) += (ee5 * (4 - ee87) * ee7/ee12 - 4) * ee5/ee95;
      out(j, 14) += (6 - ((ee65 + 6) * ee5 * ee7/ee12 + 1.5 * ee117)) *
        ee2/ee85;
      out(j, 15) += -(((ee74 - (ee83 * ee5 * ee7/ee12 + 1.5 * (ee118 *
        ee2/ee11)))/ee5 + 1.5 - ((1.5 * ee118 - 4.5)/ee33 + 3) *
        ee2/ee3) * ee2/ee85);
      out(j, 16) += ((ee5 * (6 - ee87) * ee7/ee12 - 7) * ee5 * ee7/
        ee12 + 1)/ee12;
      out(j, 17) += (ee5 * (10.5 - ee100) * ee7/ee12 - (1.5 + 1.5 *
        ee110)) * ee2/ee32;
      out(j, 18) += -((((ee75 + 3 * ((1.5 - ee64)/ee5))/ee3 - 3) *
        ee2/ee3 + ee114 + 1.5 - ee37 * ee115) * ee2/ee32);
      out(j, 19) += -(ee89/ee32);
      out(j, 20) += -((((ee87 - 6) * ee5 * ee7/ee12 + 7) * ee5 * ee7/
        ee12 - 1) * ee7/ee12);
      out(j, 21) += -(((ee100 - 10.5) * ee5 * ee7/ee12 + 1.5 - 1.5 *
        ee106) * ee2 * ee7/ee32);
      out(j, 22) += -((ee111 * ee37 + ((ee75 - 3 * ((ee64 - 1.5)/
        ee5))/ee3 - 3) * ee2/ee3 + ee114 + 1.5) * ee2 * ee7/ee32);
      out(j, 23) += -(ee89 * ee7/ee32);
      out(j, 24) += (((((((4.5 * ee51 - (40.5 * ee30 + 9 * ee39))/
        ee10 - 9 * ee13) * ee2 * ee7/ee34 - (ee55 * ee13 + 2 * (R_pow(ee13, 2) +
          1.5 * ee44) + ee121)) * ee7/ee12 + (1.5 * ee81 -
          (ee119 + 4.5 * ee94))/ee5 + ee120 + ee122 + ee123 - ee102)/
            ee3 + 3) * ee2/ee3 - 1.5) * ee7/ee12 + ((4.5 * ee81 - (1.5 *
              ee94 + ee119)) * ee2 * ee7/ee34 - (((((ee13 * (6 * ee21 -
              18 * ee14) + (12 * (ee38 - ee52) + 9 * (ee74 + ee78) - 324 *
              ee14) * ee2/ee86 + 3 * (4.5 * (ee13 * ee2/ee11) - ee61) - 6 *
              ee61)/ee5 - ee121)/ee5 + ee120 + 
              ee122 + ee123 - ee102)/
                ee3 + 3) * ee2/ee3 - 1.5) * ee124)/ee5) * ee2/ee86;
      
    } else {
      
      ee2 = exp(-txi);
      ee3 = 1 + ee2;
      ee5 = 1.5/ee3 - 1;
      ee6 = exp(lpsi);
      ee7 = 1/ee5;
      ee8 = y[l] - mu;
      ee9 = y[l + drop] - mu;
      ee11 = ee5 * ee8/ee6;
      ee13 = ee5 * ee9/ee6;
      ee14 = 1 + ee7;
      ee15 = ee11 + 1;
      ee16 = ee13 + 1;
      ee17 = ee7 + 2;
      ee18 = R_pow(ee3, 2);
      ee19 = R_pow(ee15, ee14);
      ee20 = R_pow(ee16, ee14);
      ee21 = log1p(ee11);
      ee22 = log1p(ee13);
      ee23 = R_pow(ee5, 2);
      ee24 = R_pow(ee15, ee17);
      ee25 = R_pow(ee16, ee17);
      ee26 = R_pow(ee15, ee7);
      ee27 = R_pow(ee16, ee7);
      ee28 = ee3 * ee5;
      ee29 = ee7 + 3;
      ee32 = (4.5/ee28 - 3) * ee2/ee3 + 1.5;
      ee33 = 3 * (ee2/ee3);
      ee34 = 1.5 - ee33;
      ee37 = 1/ee27 - 1/ee26;
      ee38 = ee8/ee19;
      ee39 = ee9/ee20;
      ee40 = ee24 * ee6;
      ee41 = ee25 * ee6;
      ee42 = 1/ee19;
      ee43 = 1/ee20;
      ee46 = 1.5 * (ee21/(ee19 * ee23));
      ee47 = 1.5 * (ee22/(ee20 * ee23));
      ee48 = R_pow(ee15, ee29);
      ee49 = R_pow(ee16, ee29);
      ee56 = ee46 - 1.5 * (ee14 * ee8/ee40);
      ee57 = ee47 - 1.5 * (ee14 * ee9/ee41);
      ee58 = 2.25/ee28;
      ee60 = ee15 * ee18 * ee6;
      ee62 = ee16 * ee18 * ee6;
      ee65 = (ee58 - 3) * ee2/ee3 + 1.5;
      ee67 = ee2 * ee8/ee60;
      ee69 = ee2 * ee9/ee62;
      ee73 = (1.5 * ee38 - 1.5 * ee39)/ee6 + (1.5 * (ee22/ee27) -  1.5 * (ee21/ee26))/ee5;
      ee75 = 2.25 * ee67 - ee32 * ee21;
      ee77 = 2.25 * ee69 - ee32 * ee22;
      ee78 = 1 + 2 * ee2;
      ee79 = ee2/ee18;
      ee80 = ee42 - ee43;
      ee81 = ee18 * ee23;
      ee82 = ee34 * ee5;
      ee83 = ee38 - ee39;
      ee84 = 3 * ee78;
      ee85 = 3 * ee3;
      ee86 = 2.25 * (ee2/ee81);
      ee87 = ee18 * ee5;
      ee88 = 12 * ee2;
      ee89 = ee14 * ee5;
      ee90 = ee14 * ee34;
      ee95 = ee56 * ee2;
      ee96 = ee57 * ee2;
      ee104 = 1.5 * (ee21/(ee24 * ee23)) - 1.5 * (ee17 * ee8/(ee48 *  ee6));
      ee106 = 1.5 * (ee22/(ee25 * ee23)) - 1.5 * (ee17 * ee9/(ee49 *  ee6));
      ee107 = ee90 + ee86;
      ee108 = 1.5/ee5;
      ee109 = 3 * ee34;
      ee110 = ee5 * ee17;
      ee120 = 1.5 * (ee21/(ee26 * ee5)) - 1.5 * (ee8/(ee19 * ee6));
      ee122 = 1.5 * (ee22/(ee27 * ee5)) - 1.5 * (ee9/(ee20 * ee6));
      ee123 = ee8/ee24;
      ee124 = ee9/ee25;
      ee126 = ee108 - 1.5 * ee14;
      ee127 = ee84 + ee85;
      ee129 = 2.25 * ee79 - ee82;
      ee130 = (ee65/ee19 - 1.5 * (ee95/ee18)) * ee8;
      ee131 = (ee65/ee20 - 1.5 * (ee96/ee18)) * ee9;
      ee132 = ee75/ee26;
      ee133 = ee77/ee27;
      ee134 = 1/ee24;
      ee135 = 1/ee25;
      ee136 = 3 * ee82;
      ee137 = 3 * ee129;
      ee138 = (ee137 - (27 * ee79 + ee136))/ee5;
      ee139 = (ee131 - ee130)/ee6;
      ee140 = ee15 * ee6;
      ee141 = ee16 * ee6;
      ee146 = 1.5 * (ee122 * ee22) - 1.5 * (ee120 * ee21);
      ee147 = 1.5 * (ee21/ee19);
      ee148 = 1.5 * (ee22/ee20);
      ee151 = (((ee138 - ee109)/ee5 + ee84 + ee85 - ee88)/ee3 +  3) * ee2/ee3 - 1.5;
      ee156 = (ee42 - ee89 * ee8/ee40) * ee8;
      ee157 = (ee43 - ee89 * ee9/ee41) * ee9;
      ee158 = ee147 - ee148;
      ee159 = ee135 - ee134;
      ee160 = ee56 * ee8;
      ee161 = ee57 * ee9;
      ee162 = ee37 * ee6;
      ee163 = ee7 + 4;
      ee164 = ee139 + (ee146 * ee2/ee87 + ee133 - ee132)/ee5;
      ee165 = ee14 * ee104;
      ee166 = ee14 * ee106;
      ee167 = ee158/ee23;
      ee168 = ee157 - ee156;
      ee169 = ee75/ee19;
      ee170 = ee77/ee20;
      ee172 = ee14 * (1.5 * ee124 - 1.5 * ee123)/ee6;
      ee174 = ee110 * ee8/ee140;
      ee176 = ee110 * ee9/ee141;
      ee178 = 1.5 * ee32 + ee109;
      ee179 = 3 * ee32;
      ee180 = (ee107/ee24 - 1.5 * (ee165 * ee2/ee18)) * ee8;
      ee181 = (ee107/ee25 - 1.5 * (ee166 * ee2/ee18)) * ee9;
      ee182 = ee160 - ee161;
      ee184 = (ee124 - ee123) * ee14 * ee5;
      ee185 = ee172 + ee167;
      ee188 = ee184/ee6 + ee42 - ee43;
      ee189 = (ee169 + 1.5 * (ee95 * ee21/ee18))/ee23;
      ee190 = (ee170 + 1.5 * (ee96 * ee22/ee18))/ee23;
      ee192 = ee180/ee6 + ee189;
      ee194 = ee181/ee6 + ee190;
      ee197 = ((ee178/ee5 + ee88 - ee127)/ee3 - 3) * ee2/ee3 +  1.5;
      ee200 = ee18 * ee37;
      ee201 = 1.5 - ((ee127 - ee88)/ee3 + 3) * ee2/ee3;
      ee202 = R_pow(ee73, 2);
      ee203 = ee5 * ee37;
      ee204 = ee14 * (3 - ee174);
      ee205 = ee14 * (3 - ee176);
      ee207 = ee34 * ee17 + ee86;
      ee208 = ee73 * ee83;
      ee209 = (ee104 * ee5 + 1.5/ee24) * ee14;
      ee210 = (ee106 * ee5 + 1.5/ee25) * ee14;
      ee213 = R_pow(ee15, ee163);
      ee216 = R_pow(ee16, ee163);
      ee219 = (2.25 * (ee8/(ee15 * ee3 * ee6)) - 3) * ee2/ee3 +  1.5;
      ee222 = (2.25 * (ee9/(ee16 * ee3 * ee6)) - 3) * ee2/ee3 +  1.5;
      ee224 = ee151 * ee21 + (1.5 * ee219 + ee179) * ee2 * ee8/ee60;
      ee226 = ee151 * ee22 + (1.5 * ee222 + ee179) * ee2 * ee9/ee62;
      ee227 = ee200 * ee6;
      ee229 = ee202 * ee2/ee87;
      ee234 = (ee126/ee24 - ee209) * ee8;
      ee235 = (ee126/ee25 - ee210) * ee9;
      ee243 = 1.5 * (ee21/(ee48 * ee23)) - 1.5 * (ee29 * ee8/(ee213 *  ee6));
      ee245 = 1.5 * (ee22/(ee49 * ee23)) - 1.5 * (ee29 * ee9/(ee216 *  ee6));
      ee246 = ee73 * ee80;
      ee247 = R_pow(ee83, 2);
      ee248 = ee5 * ee29;
      ee249 = 1/ee48;
      ee250 = 1/ee49;
      ee251 = 4 * ee2;
      ee252 = ee164 * ee37;
      ee254 = ee130/ee6 + (ee132 + 1.5 * (ee120 * ee2 * ee21/ee87))/ee5;
      ee256 = ee131/ee6 + (ee133 + 1.5 * (ee122 * ee2 * ee22/ee87))/ee5;
      ee257 = ee104 * ee2;
      ee258 = ee106 * ee2;
      ee259 = R_pow(ee80, 2);
      ee260 = (2 - ee174)/ee24;
      ee261 = (2 - ee176)/ee25;
      ee263 = ee108 - 1.5 * ee17;
      ee264 = ee134 - ee135;
      ee265 = ee224/ee26;
      ee266 = ee226/ee27;
      ee267 = ee252 + ee229;
      ee272 = (ee207/ee48 - 1.5 * (ee243 * ee17 * ee2/ee18)) *  ee8;
      ee273 = (ee207/ee49 - 1.5 * (ee245 * ee17 * ee2/ee18)) *  ee9;
      ee274 = ((1 - ee204) * ee5 * ee8/ee6 + 1)/ee24;
      ee275 = ((1 - ee205) * ee5 * ee9/ee6 + 1)/ee25;
      ee279 = (ee204 - 1) * ee5 * ee8/ee6 - 1;
      ee283 = (ee205 - 1) * ee5 * ee9/ee6 - 1;
      ee284 = ee208/ee5;
      ee285 = ((2 * (ee65 * ee56) - 1.5 * ee192) * ee2/ee18 -  ee197/ee19) * ee8;
      ee286 = ((2 * (ee65 * ee57) - 1.5 * ee194) * ee2/ee18 -  ee197/ee20) * ee9;
      ee287 = ee83 * ee80;
      ee289 = ee14 * ee201 + ee178 * ee2/ee81;
      ee293 = (1.5 * (ee56 * ee21) - 1.5 * (ee57 * ee22)) * ee2/ee18 +  ee169 - ee170;
      ee296 = (1.5 * (ee8/ee48) - 1.5 * (ee9/ee49)) * ee17/ee6 +  (1.5 * (ee22/ee25) - 1.5 * (ee21/ee24))/ee23;
      ee299 = ee126 * ee264;
      ee300 = ee126 * ee17;
      ee301 = ee75/ee24;
      ee302 = ee77/ee25;
      ee304 = 1.5 * (ee254 * ee21) + 2 * (ee120 * ee75);
      ee306 = 1.5 * (ee256 * ee22) + 2 * (ee122 * ee77);
      ee307 = ee249 - ee250;
      ee308 = (ee180 - ee181)/ee6;
      ee310 = ee279 * ee8/ee24;
      ee312 = ee283 * ee9/ee25;
      ee314 = ee272/ee6 + (ee301 + 1.5 * (ee257 * ee21/ee18))/ee23;
      ee316 = ee273/ee6 + (ee302 + 1.5 * (ee258 * ee22/ee18))/ee23;
      ee317 = (ee234/ee6 + ee46) * ee8;
      ee318 = (ee235/ee6 + ee47) * ee9;
      ee319 = (ee286 - ee285)/ee6;
      ee320 = ee182 * ee37;
      ee321 = ee247/ee6;
      ee322 = ee306 - ee304;
      ee323 = R_pow(ee6, 2);
      ee325 = (ee265 + ee322 * ee2/ee87 - ee266)/ee5 + ee319;
      ee326 = ee192 * ee8;
      ee327 = ee194 * ee9;
      ee328 = ee139 + ((ee146/ee5 - 2 * (ee202/ee37)) * ee2/ee18 +  ee133 - ee132)/ee5;
      ee329 = ee312 - ee310;
      ee330 = ee188 * ee83;
      ee331 = ee185 * ee37;
      ee334 = ee73 * ee182 * ee2/ee18;
      ee336 = ee246/ee5;
      ee337 = ee182 * ee80;
      ee338 = ee168 * ee80;
      ee339 = ee168 * ee37;
      ee341 = ((ee58 + 3) * ee2/ee3 - 1.5)/ee5 + 2 * ee90;
      ee342 = ((6 * (2 * ee78 + ee251) + 8 * ee127 - 96 * ee2)/ee3 +  12) * ee2;
      ee348 = R_pow(ee37, 2);
      ee350 = 2 * (ee208/ee203);
      ee351 = 3 - ee248 * ee8/ee140;
      ee352 = 3 - ee248 * ee9/ee141;
      ee353 = 3 * (1 + ee251);
      ee356 = 4.5/ee5 - 4.5 * ee14;
      ee357 = 6 * ee3;
      ee358 = 9 * ee78;
      ee359 = ee224/ee19;
      ee360 = ee226/ee20;
      ee361 = ee164 * ee83;
      ee362 = ee164 * ee80;
      ee363 = ee308 + ee293/ee23;
      ee368 = ((ee263/ee48 - ee243 * ee5 * ee17) * ee14 + ee300/ee48) *  ee8;
      ee369 = ((ee263/ee49 - ee245 * ee5 * ee17) * ee14 + ee300/ee49) *  ee9;
      ee370 = ee188 * ee37;
      ee371 = ee275 - ee274;
      ee378 = ee185 * ee73 * ee2/ee18;
      ee380 = ee185 * ee83 + ee337;
      ee381 = ee293/ee5;
      ee382 = ee73 * ee168;
      ee383 = ee284 + ee320;
      ee386 = ee73 * ee159/ee37;
      ee389 = ee296 * ee14 * ee5 + ee299;
      ee390 = ee296 * ee5;
      ee391 = ee338 - ee330;
      ee392 = ee339 - ee321;
      ee395 = (ee17 * ee351 + ee7) * ee5 * ee8/ee6;
      ee398 = (ee17 * ee352 + ee7) * ee5 * ee9/ee6;
      ee399 = ((2 * (ee107 * ee104) - 1.5 * (ee314 * ee14)) *  ee2/ee18 - ee289/ee24) * ee8;
      ee400 = ((2 * (ee107 * ee106) - 1.5 * (ee316 * ee14)) *  ee2/ee18 - ee289/ee25) * ee9;
      ee401 = ee287/ee6;
      ee403 = ee83 * ee159/ee162;
      ee405 = ee89 * ee159;
      ee406 = ee87 * ee37;
      ee407 = ee200 * ee323;
      ee409 = ee158/ee5;
      ee411 = ee80 * ee307;
      ee412 = ee157 + 2 * (ee247/ee162);
      ee414 = ee138 + 1.5 * ee34 + 4.5 * ee32;
      ee415 = 1.5 * (ee192 * ee21);
      ee416 = 1.5 * (ee194 * ee22);
      ee417 = 2 * ee267;
      ee418 = 2 * (ee246/ee37);
      ee419 = 2 * (ee56 * ee75);
      ee420 = 2 * (ee57 * ee77);
      ee421 = 4 * ee331;
      ee422 = 8 * ee284;
      ee423 = 8 * ee336;
      ee424 = 8 * ee229;
      ee425 = R_pow(ee6, 3);
      ee426 = ee267 * ee83;
      ee431 = ((((ee34 * (6 * ee129 - 18 * ee79) + (12 * (ee137 -  ee136) + 9 * (2 * ee129 + 9 * ee79) - 324 * ee79) * ee2/ee87 +  3 * (4.5 * (ee34 * ee2/ee18) - ee201 * ee5) - 6 * (ee201 *  ee5))/ee5 - 3 * ee201)/ee5 + ee353 + ee357 + ee358 -  ee342)/ee3 + 3) * ee2/ee3 - 1.5;
      ee432 = ee361 + ee334;
      ee433 = ee362 + ee378;
      ee434 = ee326 - ee327;
      ee439 = ee328 + 2 * ee164;
      ee440 = (ee341/ee24 + 2 * (ee104 * ee126 * ee2/ee18) - (ee314 *  ee5 + 1.5 * (ee257/ee18)) * ee14) * ee8;
      ee441 = (ee341/ee25 + 2 * (ee106 * ee126 * ee2/ee18) - (ee316 *  ee5 + 1.5 * (ee258/ee18)) * ee14) * ee9;
      ee445 = (((1.5 * ee151 - (3 * (ee32 * ee34) + 4.5 * ee201))/ee5 +  ee353 + ee357 + ee358 - ee342)/ee3 + 3) * ee2/ee3 - 1.5;
      ee447 = (ee390 - 2 * ee386) * ee14 + ee299;
      ee449 = (ee234 - ee235)/ee6 + ee167;
      ee450 = ee318 - ee317;
      ee451 = ee188 * ee73;
      ee452 = ee188 * ee159;
      ee453 = (ee14 * R_pow(ee159, 2) + ee411 * ee17) * ee5;
      ee454 = ee382/ee5;
      ee455 = ee399/ee6;
      ee456 = ee400/ee6;
      ee457 = (ee356/ee24 - ((ee368/ee6 + 2 * ee165) * ee5 + ee209)) *  ee8;
      ee458 = (ee356/ee25 - ((ee369/ee6 + 2 * ee166) * ee5 + ee210)) *  ee9;
      ee459 = ee405 * ee37;
      ee461 = ee259 * ee159/ee37;
      ee463 = ee412 + 2 * ee168 - ee156;
      ee464 = ee37 * ee323;
      ee465 = ee37 * ee425;
      ee467 = (2 - ee395)/ee48 + ee260;
      ee469 = (2 - ee398)/ee49 + ee261;
      ee470 = ee261 - ee260;
      ee471 = ee416 + ee420;
      ee472 = ee250 - ee249;
      ee473 = ee417 - ee424;
      ee474 = 2 * (ee380 * ee37);
      ee476 = 2 * ee383;
      ee477 = 2 * ee392;
      ee478 = 2 * ee403;
      ee480 = 2 * ee160 - (ee350 + 2 * ee161);
      ee481 = 2 * ee104;
      ee482 = 2 * ee106;
      ee483 = 4 * ee370;
      ee484 = ee421 - ee423;
      ee485 = 4 * ee320;
      ee486 = 4.5 * ee151;
      ee487 = 4.5 * ee34;
      ee488 = 8 * ee401;
      ee489 = 8 * ee321;
      
      out(j, 0) += ((ee5 * ee307 * ee17 - ee80 * ee159/ee37) * ee14 *
        ee5 - ee80 * (2 * ee405 - 2 * (ee259/ee37))/ee37)/ee465;
      out(j, 1) += ((ee261 - (ee403 + ee260)) * ee14 * ee5 - ee80 *
        (2 * ee188 - 2 * (ee287/ee162))/ee37)/ee464;
      out(j, 2) += -(((ee390 - ee386) * ee14 + ee299 - ee80 * (2 *
        ee185 - 2 * (ee246/ee203))/ee37) * ee2/ee407);
      out(j, 3) += -(((((2 * ee184 - 2 * (ee287/ee37))/ee6 + 2/ee19 -
        2/ee20) * ee83 - ee338)/ee162 + ee275 - ee274)/ee162);
      out(j, 4) += -(((ee234 - ((((ee409 - ee418)/ee5 + ee172) * ee83 +
        ee337)/ee37 + ee235))/ee6 + (ee409 - ee451/ee37)/ee5) *
        ee2/ee227);
      out(j, 5) += (ee308 + (ee381 - (ee362 + ee73 * ((2 * ee409 -
        ee418)/ee5 + 2 * ee172) * ee2/ee18)/ee37)/ee5) * ee2/ee227;
      out(j, 6) += (ee312 + ee463 * ee83/ee162 - ee310)/ee162;
      out(j, 7) += (ee318 + (ee83 * (2 * ee182 - ee350)/ee6 - ee454)/
        ee37 - ee317) * ee2/ee227;
      out(j, 8) += (ee326 - ((ee361 + ee73 * ee480 * ee2/ee18)/ee203 +
        ee327)) * ee2/ee227;
      out(j, 9) += -(((ee265 + (ee322/ee5 - ee439 * ee73/ee37) * ee2/
        ee18 - ee266)/ee5 + ee319) * ee2/ee406);
      out(j, 10) += -((((ee5 * (1/ee216 - 1/ee213) * ee29 - ee411/
        ee37) * ee5 * ee17 - (ee453 - 2 * ee461)/ee37) * ee14 * ee5 -
          (ee89 * (2 * ee453 - 4 * ee461) - ee259 * (2 * (ee459 + ee259) +
          4 * ee459 - 8 * ee259)/ee348)/ee37)/(ee37 * R_pow(ee6, 4)));
      out(j, 11) += -((((ee351/ee48 - (ee83 * ee307/ee162 + ee352/
        ee49)) * ee5 * ee17 - (ee452 + (ee261 - (ee260 + ee478)) * ee80)/
          ee37) * ee14 * ee5 - (ee89 * (2 * (ee452 + ee470 * ee80) -
            4 * (ee287 * ee159/ee162)) - ee259 * (2 * (ee370 + ee401) +
            ee483 - ee488)/ee348)/ee37)/ee465);
      out(j, 12) += ((((((1.5 * (ee9/ee216) - 1.5 * (ee8/ee213)) *
        ee29/ee6 + (1.5 * (ee21/ee48) - 1.5 * (ee22/ee49))/ee23) * ee5 -
        ee73 * ee307/ee37) * ee17 + ee263 * ee472) * ee14 + ee126 *
        ee472 * ee17) * ee5 - ((ee185 * ee5 - ee418) * ee14 * ee159 +
        (ee389 - (ee80 * (2 * (ee331 + ee336) + ee421 - ee423)/
          ee37 + 4 * (ee73 * ee14 * ee159))/ee37) * ee80 + 2 * (ee389 *
            ee80 + ee185 * ee14 * ee5 * ee159))/ee37) * ee2/(ee200 *
            ee425);
      out(j, 13) += (((ee398 - 2)/ee49 + (ee83 * (2 * ee261 - (ee478 +
        2 * ee260)) - ee168 * ee159)/ee162 + (ee176 - 2)/ee25 -
        ((ee395 - 2)/ee48 + (ee174 - 2)/ee24)) * ee14 * ee5 - (ee80 *
        (4 * ee330 - (ee83 * (ee488 - ee483) + 2 * (ee392 * ee80))/
          ee37)/ee162 + 2 * (ee371 * ee80 - R_pow(ee188, 2)))/ee37)/
            ee464;
      out(j, 14) += (((ee369 - ee368)/ee6 + ee14 * (ee482 - ee481)) *
        ee5 + ee126 * (2/ee24 - 2/ee25) - ((ee447 * ee83 + ee182 *
        ee14 * ee5 * ee159)/ee6 + ee73 * ee470 * ee14 + 2 * (ee449 *
        ee80 + ee188 * ee185) - ((ee83 * ee484 + 2 * (ee383 * ee80))/
          ee162 + 4 * (ee451/ee5)) * ee80/ee37)/ee37) * ee2/ee407;
      out(j, 15) += -((((ee273 - ee272)/ee6 + ((1.5 * (ee106 * ee22) -
        1.5 * (ee104 * ee21)) * ee2/ee18 + ee302 - ee301)/ee23) *
        ee14 * ee5 + (ee90 + (ee33 - 1.5)/ee5) * ee264 + (ee126 *
        (ee481 - ee482) - ee389 * ee73/ee203) * ee2/ee18 - ((ee447 *
        ee73 * ee2/ee18 - ((ee73 * ee484 * ee2/ee18 + 2 * (ee267 *
        ee80))/ee37 + 4 * ee378) * ee80/ee37)/ee5 + ee164 * ee14 * ee159 +
        2 * (ee363 * ee80 + R_pow(ee185, 2) * ee2/ee18))/ee37) *
        ee2/ee407);
      out(j, 16) += -(((((ee275 + 2 * ee371 - (((ee80 * (ee477 + ee489) +
        2 * (ee391 * ee37))/ee37 + 2 * ee391)/ee162 + ee274)) *
        ee83 + ee188 * ee463 - ee329 * ee80)/ee37 + (ee469 * ee9 -
        ee467 * ee8) * ee14 * ee5)/ee6 + ee283/ee25 - ee279/ee24)/
          ee162);
      out(j, 17) += -((((ee188 * ee480 + ((2 * ee234 - ((ee80 * (ee476 -
        ee422) + ee474)/ee348 + 2 * ee235))/ee6 + 2 * ee167) *
        ee83 + 2 * (ee391 * ee73/ee203) - (ee450 * ee80 + ee185 * ee168))/
          ee37 + ee458 - ee457)/ee6 + ((ee148 - ee147)/ee5 - ee371 *
            ee73/ee37)/ee5) * ee2/ee227);
      out(j, 18) += -(((ee440 - ((ee363 * ee83 + ee434 * ee80 + (2 *
        (ee185 * ee182) - 2 * (ee380 * ee73/ee203)) * ee2/ee18 - (ee73 *
        (ee474 - 8 * (ee208 * ee80/ee5)) * ee2/ee18 + 2 * (ee426 *
        ee80))/(ee5 * ee348))/ee37 + ee441))/ee6 + (ee381 - (ee328 *
        ee188 + 2 * (ee449 * ee73 * ee2/ee18))/ee37)/ee5) * ee2/
          ee227);
      out(j, 19) += (((ee360 + (ee414 * ee80 + ee415 + ee419 - ee471) *
        ee2/ee18 - ee359)/ee5 - (ee325 * ee80 + ((ee308 + (ee381 -
        ((ee80 * ee473 + 2 * (ee433 * ee37))/ee37 + 2 * ee433)/
          ee37)/ee5 + 2 * ee363) * ee73 + ee439 * ee185) * ee2/ee18)/ee37)/
            ee5 + (ee399 - ee400)/ee6) * ee2/ee227;
      out(j, 20) += (((((ee477 + 4 * ee339 + ee489)/ee37 + 4 * ee168) *
        ee83/ee162 + 2 * ee329) * ee83 + ee168 * (ee412 - ee156) +
        2 * (ee329 * ee83 + R_pow(ee168, 2)))/ee162 + (ee275 - ee469 *
        ee14 * ee5 * ee9/ee6) * ee9 - (ee274 - ee467 * ee14 *
        ee5 * ee8/ee6) * ee8)/ee162;
      out(j, 21) += ((((ee318 + (ee83 * (ee476 + ee485 - ee422)/ee162 -
        4 * ee454)/ee37 - ee317) * ee83 + (ee160 - (ee161 + ee350)) *
        ee168 + 2 * (ee450 * ee83 + ee182 * ee168))/ee6 - ee329 *
        ee73/ee5)/ee37 + (ee457/ee6 + ee46) * ee8 - (ee458/ee6 +
        ee47) * ee9) * ee2/ee227;
      out(j, 22) += ((ee441/ee6 + ee190) * ee9 + ((2 * (ee434 * ee83 +
        R_pow(ee182, 2) * ee2/ee18) - ((ee73 * (ee485 - ee422) *
        ee2/ee18 + 2 * ee426)/ee37 + 4 * ee334) * ee83/ee203)/ee6 -
        (ee164 * ee168 + ee73 * (2 * ee318 - (2 * ee317 + 2 * (ee382/
          ee203))) * ee2/ee18)/ee5)/ee37 - (ee440/ee6 + ee189) * ee8) *
            ee2/ee227;
      out(j, 23) += ((((ee414/ee19 + ee415 + ee419) * ee2/ee18 - ee359)/
        ee23 + ee455) * ee8 - ((ee325 * ee83 + (ee73 * (3 * ee326 -
          (((ee83 * ee473 + 2 * (ee432 * ee37))/ee37 + 2 * ee432 +
          2 * ee334)/ee203 + 3 * ee327)) + 3 * (ee164 * ee182)) * ee2/
            ee18)/ee203 + (((ee414/ee20 + ee416 + ee420) * ee2/ee18 -
              ee360)/ee23 +   ee456) * ee9)) * ee2/ee227;
      out(j, 24) += -(((((1.5 * (((ee306 * ee2/ee87 - ee266)/ee5 +
        ee286/ee6) * ee22) + 3 * (ee224 * ee120) + 3 * (ee256 * ee77) -
        (1.5 * (((ee304 * ee2/ee87 - ee265)/ee5 + ee285/ee6) * ee21) +
        3 * (ee226 * ee122) + 3 * (ee254 * ee75)))/ee5 - (ee328 *
        ee164 + ee73 * (2 * ee325 - ee73 * ((ee417 + 4 * ee252 -
        ee424)/ee37 + 4 * ee164) * ee2/ee406) + 2 * (ee325 * ee73 +
        R_pow(ee164, 2)))/ee37) * ee2/ee18 + ((ee486 - (1.5 * ((((ee487 +
        6.75 * ee67) * ee8/ee140 + ee88 - ee127)/ee3 - 3) *
        ee2/ee3 + 1.5) + 3 * (ee219 * ee32))) * 
        ee2 * ee8/ee60 -
        ee431 * ee21)/ee26 - ((ee486 - (1.5 * ((((ee487 + 6.75 * ee69) *
        ee9/ee141 + ee88 - ee127)/ee3 - 3) * ee2/ee3 + 1.5) +
        3 * (ee222 * ee32))) * ee2 * ee9/ee62 - ee431 * ee22)/ee27)/
          ee5 + (((3 * (ee194 * ee65) - (1.5 * ((ee471 * ee2/ee18 - ee360)/
            ee23 + ee456) + 3 * (ee197 * ee57))) * ee2/ee18 - ee445/
              ee20) * ee9 - ((3 * (ee192 * ee65) - (1.5 * (((ee415 + ee419) *
                ee2/ee18 - ee359)/ee23 + ee455) + 3 * (ee197 * ee56))) *
                ee2/ee18 - ee445/ee19) * ee8)/ee6) * ee2/ee406);
      
    }
    
  }
  
  }

}

return out;

}

