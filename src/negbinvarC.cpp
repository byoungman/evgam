// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// //' Negative binomial distribution negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for neg bin log mean and log overdispersion
// //' @param X1 a design matrix for the neg bin log mean parameter
// //' @param X2 a design matrix for the neg bin log overdispersion mean parameter
// //' @param ymat a matrix
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return negbind0 a scalar, the negative log-likelihood
// //' @return negbind12 a matrix, first then second derivatives w.r.t. negbinson parameters
// //' @return negbind34 a matrix, third then fourth derivatives w.r.t. negbinson parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double negbind0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, 
                arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat wmat, arma::mat offmat)
{
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, w, pars1, pars2, offset;
  double mu, size, p, sigsq;
  double nllh = 0.0;
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      w = wmat(j, l);
      offset = offmat(j, l);

      mu = exp(log(offset) + pars1);
      sigsq = mu + exp(pars2);
      size = mu * mu / (sigsq - mu);
      p = size / (size + mu);
      nllh -= lgamma(y + size) - lgamma(size) - lgamma(y + 1) + size * log(p) + y * log(1 - p);

    } 
    
  }
  
  return(nllh);
  
}

// //' @rdname negbind0
// [[Rcpp::export]]
arma::mat negbind12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, 
                    arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat wmat, arma::mat offmat)
{
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, w, pars1, pars2, offset;
  
  double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee12, ee13, ee14, ee15, ee16, ee17, ee18;
  double ee21, ee22, ee23, ee24, ee25, ee29;
  double ee31, ee32, ee33, ee34, ee35, ee36, ee39;
  double ee42, ee43, ee44, ee45, ee46, ee47, ee48;
  
  arma::mat out = arma::mat(nobs, 5, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];

    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      w = wmat(j, l);
      offset = offmat(j, l);
      
      ee1 = log(offset);
      ee2 = ee1 + pars1;
      ee3 = exp(ee2);
      ee4 = exp(pars2);
      ee5 = ee3/ee4;
      ee6 = 1 + ee5;
      ee7 = R_pow(ee3, 2);
      ee8 = R_pow((ee6 * ee3), 2);
      ee9 = ee7/ee4;
      ee10 = 1 + 2 * ee5;
      ee12 = ee3/(ee6 * ee4);
      ee13 = ee8 * ee4;
      ee14 = R_pow(ee3, 3);
      ee15 = 1 - ee12;
      ee16 = 1/ee6;
      ee17 = 2/ee6;
      ee18 = ee9 + y;
      ee21 = ee16 - ee14/ee13;
      ee22 = ee17 - ee10 * ee7/ee8;
      ee23 = Rf_digamma(ee18);
      ee24 = Rf_digamma(ee9);
      ee25 = log1p(ee5);
      ee29 = 2 * (ee10 * ee6 * ee7/ee8);
      ee31 = R_pow(ee12, 2) * ee4;
      ee32 = 4/ee6;
      ee33 = ee2 - (ee25 + pars2);
      ee34 = Rf_trigamma(ee18);
      ee35 = Rf_trigamma(ee9);
      ee36 = ee15 * ee4;
      ee39 = ee10 * (4 - ee29) + 1 + 4 * ee5;
      ee42 = ee21 * ee22;
      ee43 = R_pow(ee21, 2);
      ee44 = R_pow(ee22, 2);
      ee45 = 2 * (ee6 * ee14/ee13);
      ee46 = 2 * ee33;
      ee47 = 2 * ee23;
      ee48 = 2 * ee24;
      
      out(j, 0) += w * (-(((ee6 * ee22 + ee46 + ee47 - ee48) * ee3 -
        y * ee22/ee15) * ee3/ee4));
      out(j, 1) += w * (-(ee3 * (y * ee21/ee15 - (ee6 * ee21 + ee23 +
        ee1 + pars1 - (ee24 + ee25 + pars2)) * ee3)/ee4));
      out(j, 2) += w * (-(((ee6 * (4 * ee22 + ee32 - ee39 * ee7/ee8) +
        (4 * ee34 - (ee44/ee31 + 4 * ee35)) * ee7/ee4 + 4 * ee33 +
        4 * ee23 - 4 * ee24) * ee3 - y * ((ee44/ee36 - ee39 * ee3/
          ee8) * ee3 + ee32)/ee15) * ee3/ee4));
      out(j, 3) += w * (((ee6 * (2 * ee21 + ee32 - ((8 - ee29) * ee3/
        ee4 + 2) * ee7/ee8) + (2 * ee34 - (ee42/ee31 + 2 * ee35)) *
          ee7/ee4 + ee46 + ee47 - ee48) * ee3 - y * ((ee42/ee36 - ((6 -
          ee29) * ee3/ee4 + 1) * ee3/ee8) * ee3 + ee17)/ee15) * ee3/
            ee4);
      out(j, 4) += w * (-(ee3 * (y * (((3 - ee45) * ee7/ee8 - ee43/
        ee15) * ee3/ee4 - ee16)/ee15 - ((ee43/ee31 + ee35 - ee34) *
          ee7/ee4 + ((5 - ee45) * ee14/ee13 - 3/ee6) * ee6 + ee24 + ee25 +
          pars2 - (ee23 + ee1 + pars1)) * ee3)/ee4));      
    }
    
  }
  
  return out;
  
}

// //' @rdname negbind0
// [[Rcpp::export]]
arma::mat negbind34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, 
                    arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat wmat, arma::mat offmat)
{
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, w, pars1, pars2, offset;
  
  double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee19;
  double ee21, ee22, ee23, ee25, ee26, ee27, ee28, ee29;
  double ee31, ee32, ee33, ee34, ee35, ee36, ee37, ee38, ee39;
  double ee40, ee42, ee45, ee47;
  double ee50, ee53, ee54, ee57, ee58, ee59;
  double ee60, ee61, ee62, ee63, ee64, ee65, ee67, ee68, ee69;
  double ee70, ee71, ee72, ee73, ee74, ee75, ee76, ee77, ee78, ee79;
  double ee80, ee81, ee82, ee83, ee84, ee85, ee86, ee87, ee89;
  double ee90, ee91, ee92;
  double ee100, ee101, ee102, ee103, ee104, ee106, ee108, ee109;
  double ee112, ee113, ee114, ee115, ee116, ee119;
  double ee124, ee125, ee126, ee127, ee128;
  double ee130, ee134, ee139;
  double ee142, ee144, ee148, ee149;
  double ee150, ee151, ee152, ee154, ee155, ee156, ee157, ee158, ee159;
  double ee161, ee162, ee163, ee164, ee167, ee168, ee169;
  double ee173, ee174, ee175, ee176, ee177, ee178, ee179;
  double ee180, ee181, ee182, ee183, ee184, ee186, ee188;
  double ee190, ee195, ee196, ee198, ee199;
  double ee201, ee202, ee204, ee206, ee208;
  double ee211, ee212, ee213, ee214, ee216, ee218, ee219;
  double ee221, ee223, ee224, ee225, ee226, ee227;
  double ee230, ee235;
  double ee243, ee245, ee247, ee249;
  double ee252, ee254, ee256, ee257, ee258, ee259;
  double ee264, ee267, ee268, ee269;
  double ee272, ee274, ee275, ee276, ee277;
  double ee282, ee283, ee285, ee286, ee287, ee288, ee289;
  double ee290, ee291, ee292, ee293, ee294, ee295, ee296, ee297, ee298, ee299;
  double ee300, ee301, ee303, ee304, ee305, ee306, ee307, ee308, ee309;
  double ee310, ee311, ee312, ee313, ee314, ee315, ee316, ee317, ee318, ee319;
  double ee320, ee321, ee322, ee323;
  
  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      w = wmat(j, l);
      offset = offmat(j, l);

      ee1 = log(offset);
      ee2 = ee1 + pars1;
      ee3 = exp(ee2);
      ee4 = exp(pars2);
      ee5 = ee3/ee4;
      ee6 = 1 + ee5;
      ee7 = R_pow((ee6 * ee3), 2);
      ee8 = R_pow(ee3, 2);
      ee9 = 2 * ee5;
      ee10 = 1 + ee9;
      ee11 = R_pow(ee3, 3);
      ee12 = ee7 * ee4;
      ee13 = 4 * ee5;
      ee14 = 1/ee6;
      ee15 = ee8/ee4;
      ee16 = 2/ee6;
      ee19 = ee10 * ee6 * ee8/ee7;
      ee21 = ee3/(ee6 * ee4);
      ee22 = 1 + ee13;
      ee23 = ee14 - ee11/ee12;
      ee25 = ee10 * ee8/ee7;
      ee26 = 2 * ee19;
      ee27 = 2 * ee10;
      ee28 = ee16 - ee25;
      ee29 = ee15 + y;
      ee31 = ee6 * ee11/ee12;
      ee32 = 1 - ee21;
      ee33 = R_pow(ee21, 2);
      ee34 = 4/ee6;
      ee35 = R_pow(ee6, 2);
      ee36 = 2 * ee31;
      ee37 = R_pow(ee10, 2);
      ee38 = ee22 * ee6;
      ee39 = 3 - ee36;
      ee40 = ee10 * (4 - ee26);
      ee42 = ee40 + 1 + ee13;
      ee45 = (6 - ee26) * ee3/ee4 + 1;
      ee47 = 1 + ee27 + ee13;
      ee50 = ee39 * ee11/ee12 - ee14;
      ee53 = ee34 - ee42 * ee8/ee7;
      ee54 = R_pow(ee4, 2);
      ee57 = ee16 - ee45 * ee8/ee7;
      ee58 = 1 + 2 * ee6;
      ee59 = 8 * ee5;
      ee60 = ee37 + ee38;
      ee61 = ee58 + ee9;
      ee62 = ee35 * ee11;
      ee63 = Rf_trigamma(ee29);
      ee64 = Rf_trigamma(ee15);
      ee65 = 4 - 8 * ee19;
      ee67 = ee61 * ee3/ee4;
      ee68 = Rf_psigamma(ee29, 2);
      ee69 = Rf_psigamma(ee15, 2);
      ee70 = ee33 * ee54;
      ee71 = 4 * ee10;
      ee72 = ee32 * ee4;
      ee73 = R_pow(ee23, 2);
      ee74 = 2 * ee60;
      ee75 = ee47 * ee6;
      ee76 = ee6 * ee33;
      ee77 = 2 * ee67;
      ee78 = 2 * ee22;
      ee79 = 4 * ee22;
      ee80 = ee23 * ee28;
      ee81 = 8/ee6;
      ee82 = ee76 * ee54;
      ee83 = ee62/ee4;
      ee84 = 8 * ee31;
      ee85 = ee62/ee12;
      ee86 = ee10 * ee65;
      ee87 = R_pow(ee28, 2);
      ee89 = ee65 * ee3/ee4;
      ee90 = 2 * ee38;
      ee91 = 4 * ee47;
      ee92 = 6 * ee22;
      ee100 = ee10 * (12 - ((ee86 + ee91) * ee6 + ee74) * ee8/ee7) +  1 + ee92 + ee59;
      ee101 = 2 * (ee75 * ee8/ee7);
      ee102 = 8 * ee85;
      ee103 = Rf_digamma(ee29);
      ee104 = Rf_digamma(ee15);
      ee106 = ee8 * ee68/ee4;
      ee108 = ee8 * ee69/ee4;
      ee109 = log1p(ee5);
      ee112 = (ee6 * (10 - ee84) + ee27) * ee11/ee12 - 7;
      ee113 = ee33 * ee4;
      ee114 = 1 + ee71;
      ee115 = 1 + ee59;
      ee116 = 16/ee6;
      ee119 = ee112 * ee11/ee12 + ee14;
      ee124 = ee10 * (4 - ((ee89 + ee78) * ee6 + ee77) * ee8/ee7) +  (20 - ee101) * ee3/ee4 + 1;
      ee125 = ee6 * (ee79 + ee13);
      ee126 = ee23 * ee53;
      ee127 = ee102 - ee27;
      ee128 = ee2 - (ee109 + pars2);
      ee130 = ee100 * ee8/ee7;
      ee134 = ee50 * ee28;
      ee139 = ee23 * ee57;
      ee142 = (14 - (ee125 - ee10 * ee127) * ee8/ee7) * ee3/ee4 +  1;
      ee144 = ee114 + ee79 + ee59;
      ee148 = 2 * ee75 + ee74;
      ee149 = ee77 + ee90;
      ee150 = 4 * ee6;
      ee151 = 4 * (2 * (ee7 * ee10) + 4 * ee83);
      ee152 = 64 * ee83;
      ee154 = ee73 * ee3/ee4;
      ee155 = 4 * ee63;
      ee156 = 4 * ee64;
      ee157 = ee81 - ee130;
      ee158 = Rf_psigamma(ee29, 3);
      ee159 = Rf_psigamma(ee15, 3);
      ee161 = ee50 * ee32;
      ee162 = ee10 * (ee152 - ee151);
      ee163 = ee80 * ee3;
      ee164 = ee80 * ee8;
      ee167 = ee57 * ee53;
      ee168 = ee87 * ee8;
      ee169 = R_pow(ee53, 2);
      ee173 = ee90 + 2 * ee37;
      ee174 = ee27 + ee150;
      ee175 = ee27 + ee13;
      ee176 = 2 * ee144;
      ee177 = ee78 + ee71;
      ee178 = 2 * ee53;
      ee179 = 2 * ee63;
      ee180 = 2 * ee64;
      ee181 = 4 * (ee7 * ee149);
      ee182 = 4 * (1 + 12 * ee5 + ee27);
      ee183 = 8 - 16 * ee19;
      ee184 = 8 * (ee23 * ee8/ee82);
      ee186 = ee8 * ee158/ee4;
      ee188 = ee8 * ee159/ee4;
      ee190 = ee124 * ee3/ee7;
      ee195 = ee142 * ee8/ee7 - ee16;
      ee196 = ee161 + ee154;
      ee198 = ee50 * ee6 + ee73 * ee8/ee70;
      ee199 = ee134 - ee139;
      ee201 = ee50/ee6 - ee73;
      ee202 = ee32 * ee57;
      ee204 = ee6 * ee57 - ee164/ee70;
      ee206 = ee6 * ee53 - ee168/ee70;
      ee208 = ee163/ee4;
      ee211 = 2 * (ee80/ee82);
      ee212 = 2 * (ee73/(ee76 * ee4));
      ee213 = 2 * ee50;
      ee214 = 2 * ee57;
      ee216 = 2 * ee106 + ee179;
      ee218 = 2 * ee108 + ee180;
      ee219 = 4 * ee67;
      ee221 = 4 * ee106 + ee155;
      ee223 = 4 * ee108 + ee156;
      ee224 = ee34 - ee124 * ee8/ee7;
      ee225 = 6 * ee5;
      ee226 = ee106 + ee63;
      ee227 = ee108 + ee64;
      ee230 = ee119 * ee28;
      ee235 = ((ee6 * (ee182 - (((4 * (2 * (ee7 * ee61) + 4 *  (ee10 * ee35 * ee8)) - 64 * (ee10 * ee35 * ee8)) * ee3/ee4 +  ee181) * ee10/ee7 + 4 * (ee148 * ee3/ee4)) * ee8/ee7) +  2 * ((ee114 + ee150 + ee13) * ee3/ee4)) * ee10 + ee47 *  (4 * ee61 - 8 * (ee10 * ee35 * ee8/ee7)) * ee3/ee4 +  2 * (ee60 * ee22)) * ee3/ee7;
      ee243 = (((4 * ee149 - (ee162 - ee181)/ee7) * ee11/ee12 -  2 * ee115) * ee10 * ee6 - (2 * ((ee58 + ee225) * ee10) +  4 * (ee61 * ee22)) * ee3/ee4) * ee3/ee7;
      ee245 = ee100 * ee3/ee7;
      ee247 = (ee162 - 4 * (ee7 * ee173))/ee7 - 4 * ee173;
      ee249 = (ee6 * (20 - 16 * ee31) + ee71) * ee3/ee4;
      ee252 = ee142 * ee3/ee7;
      ee254 = (ee39 * ee3/ee7 + ee212) * ee8/ee4;
      ee256 = (ee39 * ee8/ee7 - 2 * (ee73/ee32)) * ee3/ee4;
      ee257 = ee50 * ee53;
      ee258 = ee32 * ee53;
      ee259 = ee10 * (16 * ee85 - ee71);
      ee264 = ee10 * (32 - ((ee10 * ee183 + ee176 + 4 * ee177) *  ee6 + 4 * ee60) * ee8/ee7) + (16/ee4 - ((ee6 * (4 * ee144 -  ((4 * (ee7 * ee148) + 4 * (2 * (ee60 * ee7) + 4 * (ee37 *  ee35 * ee8)) - 64 * (ee37 * ee35 * ee8))/ee7 + 4 *  ee148) * ee10 * ee8/ee7) + 2 * (ee115 * ee6 + 3 *  (ee10 * ee22))) * ee10 + ee47 * (6 * ee60 - 8 * (ee37 *  ee35 * ee8/ee7))) * ee3/ee7) * ee3 + 1 + 24 * ee22 +  8 * ee115;
      ee267 = ee10 * (8 - ((ee89 + ee27 + ee78) * ee6 + ee77) *  ee8/ee7) + (24 - ee101) * ee3/ee4 + 2;
      ee268 = ee10/ee7;
      ee269 = ee22 * (ee102 - 8 * ee10);
      ee272 = ee6 * (ee176 + 2 * ee177) * ee8/ee7;
      ee274 = ee6 * (4 - ee84) + ee27;
      ee275 = ee6 * (ee182 + 4 * ee175 + ee59);
      ee276 = ee23 * ee157;
      ee277 = (ee178 + ee34 - (ee42/ee7 + 2 * (ee87/ee82)) * ee8) *  ee28;
      ee282 = ee87 * ee3/ee4;
      ee283 = (ee152 - (4 * (ee7 * ee174) + ee151))/ee7;
      ee285 = ee183 * ee3/ee4;
      ee286 = 12 * ee22;
      ee287 = 12 * ee28;
      ee288 = 2 * (ee119 * ee23 + R_pow(ee50, 2));
      ee289 = 2 * ee119;
      ee290 = 2 * (ee195 * ee28 - R_pow(ee57, 2));
      ee291 = 2 * ee196;
      ee292 = 2 * ee201;
      ee293 = 2 * ee139;
      ee294 = 2 * (ee167 + ee28 * ee224);
      ee295 = 2 * (ee28 * ee157 + ee169);
      ee296 = 2 * ee175;
      ee297 = 2 * ee221;
      ee298 = 2 * ee223;
      ee299 = 2 * ee128;
      ee300 = 2 * ee103;
      ee301 = 2 * ee104;
      ee303 = 4 * ee202 + 8 * ee208;
      ee304 = 4 * ee23;
      ee305 = 4 * ee174;
      ee306 = 4 * ee57;
      ee307 = 4 * ee28;
      ee308 = 4 * ee128;
      ee309 = 4 * ee103;
      ee310 = 4 * ee104;
      ee311 = 4 * ee68;
      ee312 = 4 * ee69;
      ee313 = 6 * ee206;
      ee314 = 6 * ee115;
      ee315 = 6/ee6;
      ee316 = 8 * ee154;
      ee317 = 8 * ee50;
      ee318 = 8 * ee23;
      ee319 = 8 * ee128;
      ee320 = 8 * ee103;
      ee321 = 8 * ee104;
      ee322 = 8 * ee63;
      ee323 = 8 * ee64;
      
      out(j, 0) += -(((ee6 * (ee287 + ee81 - ee130) + (16 * ee63 +
        ee297 - (ee277/ee113 + 16 * ee64 + ee298)) * ee8/ee4 + ee313 +
        ee319 + ee320 - ee321) * ee3 - y * ((((2 * (ee87/ee72) -
        ee42 * ee3/ee7) * ee3 + ee178 + ee34) * ee28/ee72 - ee245) *
        ee3 + ee81)/ee32) * ee3/ee4);
      out(j, 1) += ((ee6 * (ee304 + ee307 + ee81 - ee267 * ee8/ee7) +
        (2 * ee216 + ee322 - ((ee126 + (ee214 + ee16 - (ee268 +
        ee211) * ee8) * ee28)/ee113 + 2 * ee218 + ee323)) * ee8/ee4 +
        4 * ee204 + ee308 + ee309 - ee310) * ee3 - y * (((ee126 +
        (2 * (ee163/ee72) + ee214) * ee28)/ee72 - ee190) * ee3 + ee34)/
          ee32) * ee3/ee4;
      out(j, 2) += ((((ee23 * (ee81 - ((2 + 2 * ee45 + ee13)/ee7 +
        ee211) * ee8) - ee134)/ee113 + 2 * ee227 + ee156 - (2 * ee226 +
        ee155)) * ee8/ee4 + (((28 - (ee274 * ee10 + ee125) * ee8/
          ee7) * ee3/ee4 + 4) * ee8/ee7 - (ee304 + ee81)) * ee6 + 2 *
            ee198 + ee301 - (ee299 + ee300)) * ee3 - y * (((ee134 - ((2 *
            (ee80/ee72) - 2 * (ee45 * ee3/ee7)) * ee3 + ee34) * ee23)/
              ee72 + ee252) * ee3 - ee16)/ee32) * ee3/ee4;
      out(j, 3) += -(ee3 * (y * ((ee112 * ee8/ee7 - (ee256 + ee213 -
        ee14) * ee23/ee32) * ee3/ee4 + ee14)/ee32 - ((((((6 - ee36) *
        ee3/ee7 + ee212) * ee8/ee4 + ee213 - ee34) * ee23/ee33 +
        ee8 * (ee68 - ee69))/ee4 + 3 * ee63 - 3 * ee64) * ee8/ee4 +
        (((ee6 * (16 - ee84) + ee27) * ee11/ee12 - 19) * ee11/ee12 +
        7/ee6) * ee6 + ee103 + ee1 + pars1 - (ee104 + ee109 + pars2)) *
        ee3)/ee4);
      out(j, 4) += -(((ee6 * (ee116 + 32 * ee28 - ee264 * ee8/ee7) +
        (12 * ee221 + 2 * ((16 * ee68 + 2 * (4 * ee186 + ee311)) *
        ee8/ee4 + ee322) + 48 * ee63 - (((ee116 - (((10 * ee53 - 8 *
        (ee168/ee82))/ee6 + 2 * (ee87 + ee53/ee6)) * ee28/ee70 +
        2 * (ee100/ee7)) * ee8) * ee28 + ee169 + ee295)/ee113 + 12 *
        ee223 + 2 * ((16 * ee69 + 2 * (4 * ee188 + ee312)) * ee8/ee4 +
        ee323) + 48 * ee64)) * ee8/ee4 + 16 * ee128 + 16 * ee103 +
        24 * ee206 + 8 * (ee6 * ee157 - ee277 * ee8/ee70) - 16 *
        ee104) * ee3 - y * (((((((2 * (ee258 - ee282) + 4 * ee258 +
        8 * ee282)/ee32 + 6 * ee53) * ee28/ee72 - 2 * ee245) * ee3 +
        ee116) * ee28 + ee169 + ee295)/ee72 - ee264 * ee3/ee7) * ee3 +
        ee116)/ee32) * ee3/ee4);
      out(j, 5) += ((ee6 * (ee287 + ee116 + ee318 - (((80 - ee272)/
        ee4 - ee235) * ee3 + ee10 * (24 - ((ee86 + ee285 + ee296 +
          ee91) * ee6 + ee74 + ee219) * ee8/ee7) + ee286 + 2) * ee8/ee7) +
          (2 * ((2 * (2 * ee186 + 2 * ee68) + 8 * ee68) * ee8/ee4 +
          ee155) + ee297 + 24 * ee63 + 8 * ee216 - ((ee276 + (ee178 +
          ee81 - (((((2 - ee184) * ee28 + ee306)/ee6 + 2 * (ee80 +
          ee57/ee6)) * ee28 + 6 * (ee126/ee6))/ee70 + ee267/ee7) * ee8) *
          ee28 + ee167 + ee294)/ee113 + 2 * ((2 * (2 * ee188 + 2 *
          ee69) + 8 * ee69) * ee8/ee4 + ee156) + ee298 + 24 * ee64 +
          8 * ee218)) * ee8/ee4 + 12 * ee204 + 6 * (ee6 * ee224 - (ee126 +
          (ee214 - 2 * (ee164/ee82)) * ee28) * ee8/ee70) + ee313 +
          ee319 + ee320 - ee321) * ee3 + y * (((((ee190 + ((2 * (ee208 -
          ee202) - ee303) * ee28/ee32 - 6 * ee126)/ee72) * ee3 -
          ee34) * ee28 - (ee276 + ee167 + ee294))/ee72 + (((72 - ee272)/
            ee4 - ee235) * ee3 + ee10 * (12 - ((ee285 + ee296) * ee6 +
              ee219) * ee8/ee7) + 1 + ee92) * ee3/ee7) * ee3 - ee81)/ee32) *
              ee3/ee4;
      out(j, 6) += (((((ee243 + (96 - (ee274 * ee47 + ee275 - ee259) *
        ee8/ee7)/ee4) * ee3 + ee10 * (16 - (ee6 * (2 * ee89 + ee27 +
        ee79) + ee219) * ee8/ee7) + 4) * ee8/ee7 - (ee116 + ee307 +
        ee318)) * ee6 + ((ee23 * (ee116 - ((2 + 2 * ee124 + 2 *
        ee40 + ee59)/ee7 + 2 * (ee126/ee82)) * ee8) + ee28 * (ee16 +
        ee306 - (((ee28 * (4 - ee184) + 8 * ee57) * ee23/ee6 - 2 *
        (ee201 * ee28))/ee70 + ee268) * ee8) - (ee257 + ee290))/ee113 +
        12 * ee64 + 2 * ((2 * (ee188 + ee69) + ee312) * ee8/ee4 +
        ee180) + 4 * ee218 + 4 * ee227 - (12 * ee63 + 2 * ((2 *
        (ee186 + ee68) + ee311) * ee8/ee4 + ee179) + 4 * ee216 + 4 *
        ee226)) * ee8/ee4 + 4 * (ee195 * ee6 + (ee23 * (ee34 - (2 *
        (ee45/ee7) + ee211) * ee8) - ee134) * ee8/ee70) + 4 * ee198 +
        ee310 - (ee308 + ee309 + 8 * ee204)) * ee3 - y * ((((ee243 +
        (52 - (ee275 - (ee47 * ee127 + ee259)) * ee8/ee7)/ee4) *
        ee3 + 1 + ee71) * ee3/ee7 + (((2 * (ee196 * ee28) - ee23 *
        ee303)/ee32 - 4 * ee139) * ee28 * ee3/ee72 + ee257 + ee290 -
        ((2 * (ee126/ee72) - 2 * ee190) * ee3 + ee81) * ee23)/ee72) *
        ee3 - ee34)/ee32) * ee3/ee4;
      out(j, 7) += ((((((((ee247 * ee8/ee7 + 12) * ee3/ee4 + ee286 +
        ee314) * ee6 + ee249 + (ee6 * (6 - 24 * ee31) + 6 * ee10) *
        ee10 - ee269) * ee8/ee7 - 92) * ee3/ee4 - 8) * ee8/ee7 +
        ee116 + 6 * ee23) * ee6 + ((((((ee23 * (6 - ee184)/ee6 - ee292) *
        ee28 + (ee293 - 4 * ee199)/ee6)/ee70 + (3 + 3 * ee142 +
        6 * ee45 + ee225)/ee7) * ee8 - 24/ee6) * ee23 + ee50 * (3 *
        ee57 + ee315 - 3 * ee25) - ee230)/ee113 + 2 * ((3 * ee68 +
        ee186) * ee8/ee4 + ee63) + 6 * ee226 + 6 * ee63 - (2 * ((3 *
        ee69 + ee188) * ee8/ee4 + ee64) + 6 * ee227 + 6 * ee64)) *
        ee8/ee4 + 2 * (ee119 * ee6 + (ee254 + ee213 - ee14) * ee23 *
        ee8/ee70) + ee299 + ee300 - (ee301 + 6 * ee198)) * ee3 - y *
        ((((((ee247 * ee11/ee12 + ee314) * ee6 + ee249 - ee269) *
        ee8/ee7 - 30) * ee3/ee4 - 1) * ee3/ee7 + (ee230 - ((((((ee291 -
        ee316) * ee28 + 2 * (ee199 * ee32))/ee32 + 2 * ee199 -
        ee293)/ee72 + 3 * ee252) * ee3 - ee315) * ee23 + 3 * (ee50 *
        ee57)))/ee72) * ee3 + ee16)/ee32) * ee3/ee4;
      out(j, 8) += -(ee3 * (y * ((((((24 * ee6 + ee305 - ee283) *
        ee11/ee12 - 34) * ee6 - (14 * ee10 + ee78)) * ee11/ee12 + 15) *
        ee8/ee7 - ((ee256 - ee14) * ee50 + ee23 * (ee289 - ((ee291 +
        4 * ee161 - ee316)/ee32 + 4 * ee50) * ee23 * ee3/ee72) +
        ee288)/ee32) * ee3/ee4 - ee14)/ee32 - ((((((((ee23 * (ee184 -
        8) + ee317)/ee6 + ee292) * ee23/ee113 - (18 - ee84) * ee3/
          ee7) * ee8/ee4 + 10/ee6 + ee289 - ee317) * ee23 + (ee254 -
            ee14) * ee50 + ee288)/ee33 + (6 * ee69 + ee8 * (ee159 - ee158)/
              ee4 - 6 * ee68) * ee8)/ee4 + 7 * ee64 - 7 * ee63) * ee8/
                ee4 + (((((ee305 + 56 * ee6 - ee283) * ee11/ee12 - 86) * ee6 -
                  (ee78 + 22 * ee10)) * ee11/ee12 + 65) * ee11/ee12 - 15/ee6) *
                  ee6 + ee104 + ee109 + pars2 - (ee103 + ee1 + pars1)) *
                  ee3)/ee4);
      }

  }

  return out;
  
}

// //' @param pars a list of vectors of coefficients for negbinson log mean
// //' @param X1 a sparse design matrix for the negbinson log mean parameter
// //' @param ymat a matrix
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return negbind0 a scalar, the negative log-likelihood
// //' @return negbind12 a matrix, first then second derivatives w.r.t. negbinson parameters
// //' @return negbind34 a matrix, third then fourth derivatives w.r.t. negbinson parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double negbinspd0(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::mat ymat, 
                  arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat wmat, arma::mat offmat)
{
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, w, pars1, pars2, offset;
  double mu, size, p, sigsq;
  double nllh = 0.0;
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      w = wmat(j, l);
      offset = offmat(j, l);
      
      mu = exp(log(offset) + pars1);
      sigsq = mu + exp(pars2);
      size = mu * mu / (sigsq - mu);
      p = size / (size + mu);
      nllh -= lgamma(y + size) - lgamma(size) - lgamma(y + 1) + size * log(p) + y * log(1 - p);
      
    } 
    
  }
  
  return(nllh);
  
}

// //' @rdname negbinspd0
// [[Rcpp::export]]
arma::mat negbinspd12(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::mat ymat, 
                      arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat wmat, arma::mat offmat)
{
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, w, pars1, pars2, offset;
  
  double ee1, ee2, ee3, ee4, ee5, ee6, ee8, ee9;
  double ee10, ee12, ee13, ee14, ee15, ee16, ee17, ee19;
  double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
  double ee30, ee31, ee32, ee33, ee34, ee35;
  
  arma::mat out = arma::mat(nobs, 5, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      w = wmat(j, l);
      offset = offmat(j, l);
      
      ee1 = exp(pars1);
      ee2 = log(offset);
      ee3 = ee1 + ee2;
      ee4 = exp(pars2);
      ee5 = R_pow(ee3, 2);
      ee6 = ee5/ee4;
      ee8 = ee6 + ee1 + ee2;
      ee9 = ee8 * ee4;
      ee10 = ee5/ee9;
      ee12 = 1 + 2 * (ee3/ee4);
      ee13 = ee12 * ee3;
      ee14 = ee13/ee8;
      ee15 = ee6 + y;
      ee16 = 1 - ee10;
      ee17 = 2 - ee14;
      ee19 = 2 * ee1 + ee2;
      ee20 = 2 * log(ee3);
      ee21 = Rf_digamma(ee15);
      ee22 = Rf_digamma(ee6);
      ee23 = log(ee8);
      ee24 = R_pow(ee10, 2);
      ee25 = 2 * ee14;
      ee26 = 2 * (ee20 - (ee23 + pars2));
      ee27 = 2 * ee21;
      ee28 = 2 * ee22;
      ee29 = Rf_trigamma(ee15);
      ee30 = Rf_trigamma(ee6);
      ee31 = ee8 * ee16;
      ee32 = (1 + 2 * (ee19/ee4)) * ee3;
      ee33 = ee12 * (4 - ee25);
      ee34 = R_pow(ee17, 2);
      ee35 = 2 * ee19;
      
      out(j, 0) += w * (-((2 + ee26 + ee27 - ((ee13 + y * ee17/ee16)/
        ee8 + ee28)) * ee1 * ee3/ee4));
      out(j, 1) += w * (-((ee15/ee8 + ee22 + ee23 + pars2 - (1 + ee20 +
        ee21)) * ee5/ee4));
      out(j, 2) += w * (-((((ee35 - (ee32 + ee33 * ee1) * ee3/ee8)/
        ee5 + (4 * ee29 - (ee34 * ee5/(ee24 * R_pow(ee8, 2) * ee4) +
          4 * ee30)) * ee1/ee4) * ee5 + (ee26 + ee27 - ee28) * ee19 +
          4 * (ee17 * ee1) - y * (((ee34 * ee3/(ee16 * ee4) - ee33) *
          ee1 - ee32) * ee3/ee8 + ee35)/ee31) * ee1/ee4));
      out(j, 3) += w * ((((2 - ((6 - ee25) * ee3/ee4 + 1) * ee3/ee8)/
        ee3 + (2 * ee29 - 2 * ee30) * ee3/ee4 - ((ee16 * ee17 * ee5/
          (ee24 * ee8 * ee4) + 2) * ee3/ee4 + 1)/ee8) * ee3 + 2 + 2 *
            ee16 + ee26 + ee27 - (ee28 + y * (((ee14 - 4) * ee3/ee4 -
            1) * ee3/ee8 + 2)/ee31)) * ee1 * ee3/ee4);
      out(j, 4) += w * (-((ee20 + 3 + ee21 + y * (((3 - 2 * ee10) *
        ee5/ee9 - 1)/ee16 - ee10)/ee8 - ((((R_pow(ee16, 2)/ee24 -
        2) * ee5/ee9 + 5)/ee8 + ee30 - ee29) * ee5/ee4 + ee22 + ee23 +
        pars2)) * ee5/ee4));
      
    }
    
  }
  
  return out;
  
}

// //' @rdname negbinspd0
// [[Rcpp::export]]
arma::mat negbinspd34(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::mat ymat, 
                      arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat wmat, arma::mat offmat)
{
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
  int nobs = nhere.size();
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }
  
  double y, w, pars1, pars2, offset;
  
  double ee1, ee2, ee3, ee4, ee5, ee6, ee8, ee9;
  double ee10, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee29;
  double ee31, ee32, ee33, ee36, ee37, ee38, ee39;
  double ee41, ee42, ee43, ee45, ee46, ee47, ee48;
  double ee50, ee51, ee52, ee54, ee55, ee56, ee58, ee59;
  double ee60, ee62, ee63, ee64, ee65, ee66, ee67, ee68;
  double ee70, ee71, ee72, ee73, ee74, ee75, ee76, ee77, ee78, ee79;
  double ee80, ee81, ee84, ee85, ee87, ee88, ee89;
  double ee90, ee91, ee92, ee93, ee94, ee97, ee98, ee99;
  double ee100, ee101, ee102, ee103, ee104, ee107, ee108, ee109;
  double ee110, ee111, ee112, ee113, ee115, ee116, ee117, ee118, ee119;
  double ee121, ee123, ee124, ee125, ee126, ee127, ee128, ee129;
  double ee130, ee131, ee132, ee133, ee134, ee135, ee138, ee139;
  double ee141, ee142, ee143, ee145, ee146, ee149;
  double ee150, ee151, ee153, ee155, ee156, ee157, ee158, ee159;
  double ee161, ee162, ee163, ee167;
  double ee170, ee171, ee172, ee173, ee175, ee176, ee177, ee179;
  double ee181, ee182, ee184, ee186, ee188;
  double ee190, ee191, ee192, ee193, ee194, ee197;
  double ee201, ee202, ee206, ee209;
  double ee214, ee216, ee218;
  double ee220, ee223, ee226, ee227, ee229;
  double ee230, ee232, ee233, ee234, ee235, ee237, ee239;
  double ee240, ee241, ee242, ee243, ee244, ee245, ee246, ee247, ee248, ee249;
  double ee251, ee252, ee257, ee258;
  double ee261, ee263, ee265;
  double ee271, ee272, ee273, ee275;
  double ee281, ee283, ee284, ee287, ee289;
  double ee293, ee294, ee295, ee296;
  double ee302, ee303, ee304, ee305, ee306, ee307, ee308, ee309;
  double ee310, ee311, ee316, ee317, ee318, ee319;
  double ee321, ee323, ee324, ee325, ee326, ee327, ee328, ee329;
  double ee330, ee331, ee332, ee333, ee334, ee335, ee336, ee337, ee338;
  
  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  
  for (int j=0; j < nobs; j++) {
    
    pars1 = p1vec[j];
    pars2 = p2vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      w = wmat(j, l);
      offset = offmat(j, l);
      
      ee1 = exp(pars1);
      ee2 = log(offset);
      ee3 = ee1 + ee2;
      ee4 = exp(pars2);
      ee5 = R_pow(ee3, 2);
      ee6 = ee5/ee4;
      ee8 = ee6 + ee1 + ee2;
      ee9 = ee3/ee4;
      ee10 = 1 + 2 * ee9;
      ee12 = 2 * ee1 + ee2;
      ee13 = ee8 * ee4;
      ee14 = ee10 * ee3;
      ee15 = ee5/ee13;
      ee16 = ee12/ee4;
      ee17 = 1 + 2 * ee16;
      ee18 = ee14/ee8;
      ee19 = ee17 * ee3;
      ee20 = 1 - ee15;
      ee21 = 2 - ee18;
      ee22 = 4 * ee9;
      ee23 = 2 * ee18;
      ee24 = ee10 * ee1;
      ee25 = ee6 + y;
      ee26 = R_pow(ee15, 2);
      ee27 = ee19 + 2 * ee24;
      ee29 = 4 * ee1 + ee2;
      ee31 = R_pow(ee10, 2) * ee1;
      ee32 = 1 + ee22;
      ee33 = 2 * ee15;
      ee36 = 2 * ee6 + ee1 + ee2;
      ee37 = ee10 * (4 - ee23);
      ee38 = 2 * ee12;
      ee39 = ee19 + ee37 * ee1;
      ee41 = ee39 * ee3/ee8;
      ee42 = R_pow(ee8, 2);
      ee43 = ee38 - ee41;
      ee45 = ((6 - ee23) * ee3/ee4 + 1) * ee3/ee8;
      ee46 = 2 * ee36;
      ee47 = 2 - ee45;
      ee48 = R_pow(ee4, 2);
      ee50 = 2 * (ee8 * ee32);
      ee51 = 2 * ee8;
      ee52 = ee26 * ee42;
      ee54 = (3 - ee33) * ee5/ee13;
      ee55 = Rf_trigamma(ee25);
      ee56 = Rf_trigamma(ee6);
      ee58 = ee8 * ee17 + ee31;
      ee59 = ee14 + ee51;
      ee60 = ee52 * ee48;
      ee62 = 1 + 2 * (ee29/ee4);
      ee63 = ee54 - 1;
      ee64 = 8 * ee6;
      ee65 = Rf_psigamma(ee25, 2);
      ee66 = Rf_psigamma(ee6, 2);
      ee67 = ee62 * ee3;
      ee68 = (ee67 + 2 * (ee17 * ee1)) * ee3;
      ee70 = 2 * (ee27 * ee8);
      ee71 = 8 * ee9;
      ee72 = 2 * ee10;
      ee73 = ee10 * ee12;
      ee74 = R_pow(ee21, 2);
      ee75 = 2 * ee58;
      ee76 = 2 * ee27;
      ee77 = 2 * ee59;
      ee78 = 4 * ee24;
      ee79 = 6 * ee12;
      ee80 = 8 * ee14;
      ee81 = 8 * ee31;
      ee84 = ((ee77 - ee80) * ee3/ee4 + ee50)/ee8;
      ee85 = ee10 * (ee79 - (((ee75 - ee81) * ee3 + ee70)/ee8 +  ee76 + ee78) * ee3/ee8);
      ee87 = 2 * ee19 + 2 * ee73;
      ee88 = ee72 + ee22;
      ee89 = ee46 + 4 * ee8;
      ee90 = ee20 * ee21;
      ee91 = R_pow(ee20, 2);
      ee92 = 2 * ee29;
      ee93 = 4 * ee6;
      ee94 = ee64 - ee46;
      ee97 = ee26 * ee8;
      ee98 = ee10 * (4 - (ee84 + ee22) * ee3/ee8);
      ee99 = ee74 * ee1;
      ee100 = 2 - ee33;
      ee101 = 2 * log(ee3);
      ee102 = Rf_digamma(ee25);
      ee103 = Rf_digamma(ee6);
      ee104 = log(ee8);
      ee107 = ee68 + (ee85 + 4 * ee19) * ee1;
      ee108 = (ee10 * ee94 - ee50)/ee8;
      ee109 = ee98 + ee71;
      ee110 = (ee64 - ee89)/ee8;
      ee111 = 2 * (ee27 * ee3/ee8);
      ee112 = 4 * ee16;
      ee113 = ee1 * ee5;
      ee115 = ee109 * ee1;
      ee116 = ee8 * ee20;
      ee117 = ee20 * ee43;
      ee118 = ee99 * ee5;
      ee119 = ee100 * ee12;
      ee121 = ee5 * ee65/ee4;
      ee123 = ee5 * ee66/ee4;
      ee124 = R_pow(ee3, 3);
      ee125 = R_pow(ee3, 4);
      ee126 = 2 * ee32;
      ee127 = 2 * (ee101 - (ee104 + pars2));
      ee128 = ee92 - ee107/ee8;
      ee129 = 2 * ee102;
      ee130 = 2 * ee103;
      ee131 = 2 * ee55;
      ee132 = 2 * ee56;
      ee133 = 4 * ee12;
      ee134 = ee68 + ee87 * ee1;
      ee135 = ee63 * ee21;
      ee138 = ((6 - ee110) * ee5/ee13 - 7) * ee5/ee13 + 1;
      ee139 = ee91/ee26;
      ee141 = (1 + ee112) * ee3 + ee88 * ee1;
      ee142 = 1 + ee71;
      ee143 = 2 * ee47;
      ee145 = 2 * ee43;
      ee146 = 4 * (ee46 + ee93);
      ee149 = ee97 * ee4;
      ee150 = ee20 * ee47;
      ee151 = ee20 * ee4;
      ee153 = ee119 - (ee115 + ((ee133 - ee111)/ee4 + 1) * ee3) *  ee3/ee8;
      ee155 = ee43/ee5 - ee118/ee60;
      ee156 = ee126 + ee22;
      ee157 = 4 * ee55;
      ee158 = 4 * ee56;
      ee159 = 64 * ee6;
      ee161 = 8 * ee1 + ee2;
      ee162 = Rf_psigamma(ee25, 3);
      ee163 = Rf_psigamma(ee6, 3);
      ee167 = ee141 * ee8;
      ee170 = ((14 - (ee156 - ee108) * ee3/ee8) * ee3/ee4 + 1) *  ee3/ee8 - 2;
      ee171 = ee135 - ee150;
      ee172 = ee97 * ee48;
      ee173 = ee26 * ee4;
      ee175 = ee10 * (ee159 - ee146);
      ee176 = ee47 * ee43;
      ee177 = R_pow(ee43, 2);
      ee179 = 2 * (ee58 * ee3) + ee70;
      ee181 = 2 * (ee59 * ee3/ee4) + ee50;
      ee182 = 2 * (ee8 * ee88);
      ee184 = 2 * ee121 + ee131;
      ee186 = 2 * ee123 + ee132;
      ee188 = 4 * ee15;
      ee190 = ee93 + ee1 + ee2;
      ee191 = 4 * ee17;
      ee192 = 4 * ee10;
      ee193 = 4 * ee47;
      ee194 = 4 * ee36;
      ee197 = ((ee139 - 2) * ee5/ee13 + 3) * ee5/ee13 - 1;
      ee201 = ee138 * ee20 + R_pow(ee63, 2);
      ee202 = ee8 * ee142;
      ee206 = ee90 * ee125/ee60;
      ee209 = ee47/ee3 - ee90 * ee124/ee60;
      ee214 = ee21 * ee128;
      ee216 = ee21 * ee5/ee13;
      ee218 = (2 * ee139 - 2) * ee5/ee13;
      ee220 = ee128/ee5 - ee21 * (ee145 + ee38 - ((ee37 + 2 *  (ee74 * ee124/ee172)) * ee1 + ee19) * ee3/ee8) * ee1 *  ee3/ee60;
      ee223 = ee54 - (ee91 + 1);
      ee226 = (4 - 3 * ee15) * ee5/ee13 - 1;
      ee227 = ee121 + ee55;
      ee229 = ee5 * ee162/ee4;
      ee230 = ee123 + ee56;
      ee232 = ee5 * ee163/ee4;
      ee233 = 16 * ee6;
      ee234 = 2 * ee134;
      ee235 = 2 * ee155;
      ee237 = 2 * (ee12 * ee55) + 4 * (ee113 * ee65/ee4);
      ee239 = 2 * (ee12 * ee56) + 4 * (ee113 * ee66/ee4);
      ee240 = 2 * ee63;
      ee241 = 2 * ee21;
      ee242 = ee127 + ee129;
      ee243 = 2 * ee190;
      ee244 = 4 - ee33;
      ee245 = 6 - ee188;
      ee246 = 6 * ee32;
      ee247 = 6 * ee36;
      ee248 = (((ee27 * (4 * ee59 - ee80) + 2 * (((ee19 + ee78) *  ee3 + 2 * (ee8 * ee12)) * ee10))/ee4 + 2 * (ee58 * ee32)) *  ee3 + (((4 * ee59 - 16 * ee14) * ee3/ee4 + ee182) * ee1 +  4 * ee167 - (((4 * (ee77 + 4 * ee14) - 64 * ee14) * ee3/ee4 +  4 * ee181) * ee10 * ee1 + 4 * (ee179 * ee3/ee4)) * ee3/ee8) *  ee10)/ee42;
      ee249 = ((ee27 * (6 * ee58 - ee81) + 2 * ((ee8 * ee62 +  3 * (ee17 * ee10 * ee1)) * ee10 * ee3)) * ee3 + ee10 *  (4 * (ee134 * ee8) - ((4 * (ee75 + 4 * ee31) - 64 * ee31) *  ee3 + 8 * ee179) * ee10 * ee1 * ee3/ee8))/ee42;
      ee251 = ee201/ee20;
      ee252 = (ee27 * ee94 - 2 * ee167)/ee8;
      ee257 = ee138 * ee21;
      ee258 = (((8 * ee181 - ee175) * ee5/ee13 - 2 * ee202) *  ee10 - (2 * (((1 + 6 * ee9) * ee3 + ee51) * ee10) + 4 *  (ee59 * ee32)) * ee3/ee4)/ee42;
      ee261 = ((1 + 2 * (ee161/ee4)) * ee3 + 2 * (ee62 * ee1)) *  ee3;
      ee263 = (ee10 * ee245 + 12 * ee9) * ee12;
      ee265 = (ee175 - 8 * (ee50 + 2 * (ee10 * ee36))) * ee5/ee13;
      ee271 = (ee218 + 3) * ee5/ee13;
      ee272 = ee220 * ee3;
      ee273 = ee63 * ee43;
      ee275 = ee52 * ee4;
      ee281 = ee10 * (8 - (ee84 + ee72 + ee22) * ee3/ee8);
      ee283 = ee10 * (8 * ee29 - ((((4 * ee58 - 16 * ee31) * ee3 +  2 * (ee8 * ee87))/ee8 + 2 * ee87 + 4 * ee73) * ee1 +  ee234)/ee8) + 12 * (ee17 * ee12);
      ee284 = ee32 * (ee64 - ee247);
      ee287 = ee118/ee13;
      ee289 = (2 * (ee90 * ee3/ee173) - ee192) * ee3/ee8;
      ee293 = (2 * ee88 - (ee10 * (ee233 - ee194) - ee182)/ee8) *  ee3/ee8;
      ee294 = ee12 * ee244;
      ee295 = ee12 * ee65;
      ee296 = ee12 * ee66;
      ee302 = (ee131 - ee132) * ee5/ee4 + 2 * ee20 + ee127 + ee129 -  ee130;
      ee303 = 2 * ee197;
      ee304 = 2 * (ee170 * ee21 - R_pow(ee47, 2));
      ee305 = 2 * (ee153 * ee21 + ee176);
      ee306 = 2 * (ee153/ee5 - (ee117 + ee21 * (ee143 - 2 * ee206) *  ee1) * ee5/ee60);
      ee307 = 2 * (ee226 * ee21);
      ee308 = 2 * ee141;
      ee309 = 2 * (ee214 * ee3 + ee177);
      ee310 = 2 * (ee87 * ee3/ee8);
      ee311 = 2 * ee223;
      ee316 = 2 * ee230 + ee158 - (2 * ee227 + ee157);
      ee317 = 2 * ee142;
      ee318 = 2 * ee237;
      ee319 = 2 * ee239;
      ee321 = 2 * ee184 + ee157;
      ee323 = 2 * ee186 + ee158;
      ee324 = ee242 - ee130;
      ee325 = 2 * ee161;
      ee326 = 3 * ee18;
      ee327 = 4 * ee27;
      ee328 = 4 * ee20;
      ee329 = ee193 + 8 * ee216;
      ee330 = ee146 + 8 * ee89;
      ee331 = 4 * ee65;
      ee332 = 4 * ee66;
      ee333 = 6 * ee43;
      ee334 = 8 * (ee20 * ee125/ee60);
      ee335 = 8 * ee63;
      ee336 = 8 * ee8;
      ee337 = 8 * ee17;
      ee338 = 8 * ee12;
      
      out(j, 0) += -(((((ee12 * (8 * ee55 - 8 * ee56) + ee318 - ee319)/
        ee4 + ee235) * ee1 + ee272) * ee3 + ee324 * ee29 + (4 *
          (ee155 * ee3) + 6 * (ee21 * ee12/ee3)) * ee1 - y * (((((((2 *
          (ee74 * ee3/ee151) - ee37) * ee1 - ee19) * ee3/ee8 + ee145 +
          ee38) * ee21/ee151 - ee191) * ee3 - ee85) * ee1 - ee68)/
            ee8 + ee92)/ee116) * ee1/ee4);
      out(j, 1) += ((((ee321 - ee323) * ee1 - (ee117 + ee21 * (2 +
        ee143 - ((2 + 2 * (ee90 * ee5/ee149)) * ee3/ee4 + 1) * ee3/
          ee8) * ee1) * ee5/ee275)/ee4 + (ee294 - ((ee281 + ee71) * ee1 +
            ((ee79 - ee111)/ee4 + 2) * ee3) * ee3/ee8)/ee5) * ee5 +
            ee302 * ee12 + 2 * ((2 * (ee209 * ee3) + ee241) * ee1) - y *
            ((((ee21 * (2 * ee216 + ee143)/ee20 - 8) * ee3/ee4 - ee98) *
            ee1 + (((ee76 - ee39) * ee3/ee8 - ee38)/ee4 - 1) * ee3) *
            ee3/ee8 + ee119)/ee116) * ee1/ee4;
      out(j, 2) += (((((ee20 * (8 - ((ee289 + 16) * ee3/ee4 + 4) *
        ee3/ee8) - ee135) * ee5/ee149 + 2) * ee3/ee4 + 1)/ee8 + (((26 -
        (ee156 + ee192 - ee108) * ee3/ee8) * ee3/ee4 + 3) * ee3/
          ee8 - 6)/ee3 + ee316 * ee3/ee4) * ee3 + ee303 + ee130 - (2 +
            ee127 + ee129 + ee328 + y * ((((ee108 + (2 - (ee23 + ee241)) *
            ee3/ee4 + 1 - ee126) * ee3/ee8 + ee171/ee20 + 12) * ee3/
              ee4 + 1) * ee3/ee8 - 2)/ee116)) * ee1 * ee3/ee4;
      out(j, 3) += -((ee103 + ee104 + pars2 + y * (((((ee89 - 6 *
        ee6)/ee8 + 3) * ee5/ee13 - 6) * ee5/ee13 + 1)/ee20 - (2 * (ee63/
          ee20) - ee33) * ee5/ee13)/ee8 - ((((((ee218 + 6) * ee5/
            ee13 + ee240 - 4) * ee20/ee26 + 12 - ee110) * ee5/ee13 - 19)/
              ee8 + ee5 * (ee65 - ee66)/ee4 + 3 * ee55 - 3 * ee56) * ee5/
                ee4 + ee101 + 7 + ee102)) * ee5/ee4);
      out(j, 4) += -((((((12 * ee55 - 12 * ee56) * ee29 + 2 * ((2 *
        (2 * ee295 + 4 * (ee113 * ee162/ee4)) + 8 * ee295) * ee1 *
        ee3/ee4 + 2 * (ee29 * ee55)) - 2 * ((2 * (2 * ee296 + 4 * (ee113 *
        ee163/ee4)) + 8 * ee296) * ee1 * ee3/ee4 + 2 * (ee29 *
        ee56))) * ee3 + ee12 * (6 * ee237 - 6 * ee239))/ee4 + 12 *
        (ee155 * ee12) + 6 * ee272 + 8 * (ee21 * ee29/ee3)) * ee1 +
        (((ee325 - (ee261 + (ee283 + 6 * ee67 - ee249) * ee1)/ee8)/
          ee5 - (ee21 * (4 * ee29 - ((((10 * ee43 + 2 * (ee99 + ee38 -
            ee41) - 8 * (ee99 * ee125/ee60)) * ee21 * ee5/ee172 + ee337) *
            ee3 + 2 * ee85) * ee1 + 2 * ee68)/ee8) * ee3 + ee177 +
            ee309) * ee1/ee60) * ee3 + 2 * (ee220 * ee1)) * ee3 + ee324 *
            ee161 - y * ((((((((((2 * (ee117 - ee287) + 4 * ee117 + 8 *
            ee287)/ee20 + ee333) * ee21/ee151 - ee191) * ee3 - ee85) *
            ee1 - ee68)/ee8 + ee92) * ee21 * ee3 + ee177 + ee309)/ee151 +
            ee249 + (ee214/ee151 - 6 * ee62) * ee3 - ee283) * ee1 - ee261)/
              ee8 + ee325)/ee116) * ee1/ee4);
      out(j, 5) += (((((2 * ee321 - 2 * ee323) * ee12 + 2 * (ee184 *
        ee12 + (2 * (2 * ee229 + 2 * ee65) + ee331) * ee1 * ee5/
          ee4) + ee318 - (2 * (ee186 * ee12 + (2 * (2 * ee232 + 2 * ee66) +
            ee332) * ee1 * ee5/ee4) + ee319))/ee4 + ee306 + ee235) *
            ee1 + ((ee244 * ee29 - ((ee263 + ((ee338 - ee310)/ee4 + ee337 -
            ee248) * ee3 + ee85) * ee1 + (((6 * ee29 - 2 * (ee134/
              ee8))/ee4 + 2) * ee3 + (ee112 + ee191) * ee1) * ee3)/ee8)/ee5 -
                (((ee294 + ee145 - (((((ee21 * (2 * (ee90 + 2 - ee45) +
                ee193 - 8 * ee206) * ee1 + 6 * ee117) * ee3/ee173 - ee76) *
                ee3/ee8 + ee79)/ee4 + 2) * ee3 + (ee281 + (2 * (ee74 * ee5/
                  ee149) + 8) * ee3/ee4) * ee1) * ee3/ee8) * ee21 + ee176 + ee305) *
                    ee1 + ee20 * ee128 * ee3) * ee3/ee60) * ee3) * ee3 + ee302 *
                    ee29 + (2 * ((ee306 + ee235) * ee3) + 3 * ((2 * ee209 +
                    2 * (ee21/ee3)) * ee12)) * ee1 + y * ((((((((((ee21 * (2 *
                    (((8 - ee326) * ee3/ee4 + 1) * ee3/ee8 - 2) - ee329) * ee1/
                      ee20 + ee133 - (ee111 + ee333))/ee4 + 1) * ee3 + ee115) * ee3/
                        ee8 - ee119) * ee21 - (ee176 + ee305))/ee20 + ee338 - ee310)/
                          ee4 + ee191 - ee248) * ee3 + ee263) * ee1 - ((((ee234 -
                            ee107)/ee8 - ee92)/ee4 - 1) * ee3 - (2 * ee17 + ee112) * ee1) *
                            ee3)/ee8 - ee100 * ee29)/ee116) * ee1/ee4;
      out(j, 6) += ((((((ee258 + (40 - ee293)/ee4) * ee3 + ee10 *
        (16 - (2 * ee84 + ee72 + ee71) * ee3/ee8)) * ee1 + ((18 * ee12 -
        (ee308 + ee327 - ee252) * ee3/ee8)/ee4 + 4) * ee3) * ee3/
          ee8 + ((10 - ee188) * ee5/ee13 - 8) * ee12)/ee5 + (((ee21 *
            (2 + ee193 - ((((ee21 * (4 - ee334) + 8 * ee47) * ee20 - 2 *
            (ee223 * ee21)) * ee5/ee149 + 2) * ee3/ee4 + 1) * ee3/ee8) -
            ee304) * ee1 + ((2 * ee100 + 4) * ee12 - ((((2 * (ee117 *
            ee3/ee173) - ee327) * ee3/ee8 + 12 * ee12)/ee4 + 4) * ee3 +
            (2 * ee109 + 2 * ee37) * ee1) * ee3/ee8) * ee20 - ee273) *
            ee5/ee275 + (2 * ((2 * (ee232 + ee66) + ee332) * ee5/ee4 +
            ee132) + 4 * ee186 + ee158 - (2 * ((2 * (ee229 + ee65) + ee331) *
            ee5/ee4 + ee131) + 4 * ee184 + ee157)) * ee1)/ee4) * ee5 +
            (ee316 * ee5/ee4 + ee303 + ee130 - (ee242 + ee328)) * ee12 +
            2 * (((2 * (ee170/ee3 + (ee20 * (4 - ((ee289 + 12) * ee3/
              ee4 + 2) * ee3/ee8) - ee135) * ee124/ee60) - 4 * ee209) *
                ee3 - ee241) * ee1) - y * ((((((((ee307 - 4 * ee150)/ee20 -
                ee329) * ee21 * ee5/ee13 + ee304)/ee20 + 24 - ee293)/ee4 +
                ee258) * ee3 + ee192) * ee1 + (((ee273 - ee153 * ee20)/ee20 +
                (ee33 + 6) * ee12 - (((ee111 + ee145 - ee133)/ee4 - 1) * ee3 +
                ee308 - (ee252 + ee115)) * ee3/ee8)/ee4 + 1) * ee3) * ee3/
                  ee8 + (ee245 * ee5/ee13 - 2) * ee12)/ee116) * ee1/ee4;
      out(j, 7) += (((((((((((ee20 * (ee143 + 6 * ee21) - (ee21 *
        (ee311 + 8 * (ee91 * ee125/ee60)) + 4 * ee171))/ee26 - 12) *
        ee3/ee4 + 3 * ee108 - (12 * ee10 + ee246)) * ee3/ee8 + 84) *
        ee3/ee4 + 12) * ee3/ee8 - 24) * ee20 + ee63 * (3 * ee47 +
        6 - ee326) - ee257) * ee5/ee149 - 2) * ee3/ee4 - 1)/ee8 + (((((24 -
        (ee233 - (ee194 + ee336))/ee8) * ee3/ee4 + ee317 + 6 *
        ee10 + ee246 - (ee10 * (24 * ee6 - (ee243 + ee247)) + ee284 -
        (ee265 + ee8 * (4 * ee142 + ee246)))/ee8) * ee3/ee8 - 90) *
        ee3/ee4 - 7) * ee3/ee8 + 14)/ee3 + (2 * ((ee229 + 3 * ee65) *
        ee5/ee4 + ee55) + 6 * ee227 + 6 * ee55 - (2 * ((ee232 +
        3 * ee66) * ee5/ee4 + ee56) + 6 * ee230 + 6 * ee56)) * ee3/
          ee4) * ee3 + 2 + 2 * ((((ee271 + ee240 - 1) * ee20/ee26 +
            6 - ee110) * ee5/ee13 - 7) * ee5/ee13 + 1) + ee127 + ee129 +
            6 * ee20 - (ee130 + 6 * ee197 + y * ((((((((8 * ee21 - 12) *
            ee3/ee4 + ee126 - ee108) * ee3 + ee194 + ee336)/ee8 + ee143 -
            ((2 * ee171 + ee307)/ee20 + 2)) * ee3/ee4 + ee317 - ((ee284 -
            (ee265 + 2 * (ee10 * ee190) + 4 * ee202))/ee8 + 1)) *
            ee3/ee8 + (ee257 - (2 * (ee170 * ee20) + 2 * (ee171 * ee5/ee13) +
            3 * (ee63 * ee47)))/ee20 - 28) * ee3/ee4 - 1) * ee3/ee8 +
            2)/ee116)) * ee1 * ee3/ee4;
      out(j, 8) += -((15 + ee101 + ee102 + y * (((((((ee330 - (ee89 +
        56 * ee6))/ee8 + 18) * ee5/ee4 - (14 * ee36 + ee243 + 20 *
        ee8))/ee8 + ee240 - 7) * ee5/ee13 + 14 - ee251) * ee5/ee13 -
        1)/ee20 - ((2 * ee251 - 4 * (ee63 * ee5/ee13))/ee20 - ((2 *
        ee226 + 4 * ee63)/ee20 - 8 * ee15) * ee5/ee13) * ee5/ee13)/
          ee8 - (((((((((ee20 * (ee334 - 8) + ee311 + ee335) * ee20/
            ee26 + 8) * ee5/ee13 - 18) * ee5/ee13 + 10 + 2 * ee138 - ee335) *
              ee20 + (ee271 - 1) * ee63 + 2 * ee201)/ee26 + ((56 -
              (ee159 - ee330)/ee8) * ee5/ee4 - (ee243 + 22 * ee36 + 36 * ee8))/
                ee8 - 50) * ee5/ee13 + 65)/ee8 + (ee5 * (ee163 - ee162)/
                  ee4 + 6 * ee66 - 6 * ee65) * ee5/ee4 + 7 * ee56 - 7 * ee55) *
                    ee5/ee4 + ee103 + ee104 + pars2)) * ee5/ee4);
    }
    
  }
  
  return out;
  
}
