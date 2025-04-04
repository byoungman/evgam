// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export]]
double ldenom(double left, double alpha, double beta) {
  double bl = beta * left;
  double out = 1.0 / tgamma(alpha + 1.0);
  for (int l=1; l < 13; l++) {
    out += R_pow(bl, l) / tgamma(alpha + l + 1.0);
  }
  out = out * R_pow(bl, alpha) * exp(-bl);
  return log(1 - out);
}

// //' Left-truncated Gamma distribution negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each GEV parameter
// //' @param X1 a design matrix for the shape parameter
// //' @param X2 a design matrix for the rate parameter
// //' @param ymat a matrix
// //' @param leftmat a matrix
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return ltgammad0 a scalar, the negative log-likelihood
// //' @return ltgamma12 a matrix, first then second derivatives w.r.t. parameters
// //' @return ltgamma34 a matrix, third then fourth derivatives w.r.t. parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double ltgammad0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat leftmat)
{

  int nobs = nhere.size();

  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);

  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }

  double y, left, p1, p2, alpha, beta, lnum;

  double nllh = 0.0;

  for (int j=0; j < nobs; j++) {

    p1 = p1vec[j];
    p2 = p2vec[j];
    alpha = exp(p1);
    beta = exp(p2);

    for (int l=0; l < nhere[j]; l++) {

      y = ymat(j, l);
      left = leftmat(j, l);

      lnum = (alpha - 1) * log(y) - beta * y + alpha * p2 - lgamma(alpha);
      nllh += ldenom(left, alpha, beta) - lnum;

    }

  }

return(nllh);

}

// //' @rdname ltgammad0
// [[Rcpp::export]]
arma::mat ltgammad12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat leftmat)
{

int nobs = nhere.size();
arma::mat out = arma::mat(nobs, 5, arma::fill::zeros);

  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);

  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }

  double y, left, p1, p2;

  double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee29;
  double ee30, ee31, ee32, ee35, ee36, ee39;
  double ee40, ee43, ee46, ee47;
  double ee50, ee52, ee53, ee54, ee55, ee56, ee57;
  double ee60, ee61, ee62, ee63, ee66, ee67, ee69;
  double ee70, ee71, ee76, ee79;
  double ee80, ee81, ee82, ee85, ee88, ee89;
  double ee90, ee93, ee96, ee97, ee99;
  double ee100, ee103, ee106, ee107, ee108, ee109;
  double ee110, ee113, ee114, ee115, ee116, ee117, ee118, ee119;
  double ee123, ee125, ee127, ee128, ee129;
  
  for (int j=0; j < nobs; j++) {

    p1 = p1vec[j];
    p2 = p2vec[j];

    for (int l=0; l < nhere[j]; l++) {

      y = ymat(j, l);
      left = leftmat(j, l);

      ee1 = exp(p1);
      ee2 = exp(p2);
      ee3 = left * ee2;
      ee4 = R_pow(ee3, 2);
      ee5 = 13 + ee1;
      ee6 = 12 + ee1;
      ee7 = 11 + ee1;
      ee8 = 10 + ee1;
      ee9 = 9 + ee1;
      ee10 = tgamma(ee5);
      ee11 = 8 + ee1;
      ee12 = tgamma(ee6);
      ee13 = tgamma(ee7);
      ee14 = 7 + ee1;
      ee15 = tgamma(ee8);
      ee16 = 6 + ee1;
      ee17 = tgamma(ee9);
      ee18 = tgamma(ee11);
      ee19 = 5 + ee1;
      ee20 = 4 + ee1;
      ee21 = tgamma(ee14);
      ee22 = tgamma(ee16);
      ee23 = 3 + ee1;
      ee24 = 2 + ee1;
      ee25 = tgamma(ee19);
      ee26 = tgamma(ee20);
      ee27 = ee4/ee10;
      ee29 = ee27 + 1/ee13 + ee3/ee12;
      ee30 = 1 + ee1;
      ee31 = tgamma(ee23);
      ee32 = tgamma(ee24);
      ee35 = ee29 * ee4 + 1/ee17 + ee3/ee15;
      ee36 = R_pow(ee3, ee1);
      ee39 = ee35 * ee4 + 1/ee21 + ee3/ee18;
      ee40 = tgamma(ee30);
      ee43 = ee39 * ee4 + 1/ee25 + ee3/ee22;
      ee46 = ee43 * ee4 + 1/ee31 + ee3/ee26;
      ee47 = 1/ee40;
      ee50 = ee46 * ee4 + ee47 + ee3/ee32;
      ee52 = exp(-ee3);
      ee53 = ee3/ee10;
      ee54 = Rf_digamma(ee6);
      ee55 = Rf_digamma(ee5);
      ee56 = Rf_digamma(ee8);
      ee57 = Rf_digamma(ee7);
      ee60 = ee57/ee13 + ee55 * ee4/ee10 + left * ee54 * ee2/ee12;
      ee61 = Rf_digamma(ee11);
      ee62 = Rf_digamma(ee9);
      ee63 = ee1 - 1;
      ee66 = ee60 * ee4 + ee62/ee17 + left * ee56 * ee2/ee15;
      ee67 = R_pow(ee3, ee63);
      ee69 = 11/ee12 + 12 * ee53;
      ee70 = Rf_digamma(ee16);
      ee71 = Rf_digamma(ee14);
      ee76 = ee66 * ee4 + ee71/ee21 + left * ee61 * ee2/ee18;
      ee79 = ee69 * ee4 + 10 * (ee3/ee13) + 9/ee15;
      ee80 = 1 - ee50 * ee52 * ee36;
      ee81 = Rf_digamma(ee20);
      ee82 = Rf_digamma(ee19);
      ee85 = ee76 * ee4 + ee82/ee25 + left * ee70 * ee2/ee22;
      ee88 = ee79 * ee4 + 7/ee18 + 8 * (ee3/ee17);
      ee89 = Rf_digamma(ee24);
      ee90 = Rf_digamma(ee23);
      ee93 = ee85 * ee4 + ee90/ee31 + left * ee81 * ee2/ee26;
      ee96 = ee88 * ee4 + 5/ee22 + 6 * (ee3/ee21);
      ee97 = Rf_digamma(ee30);
      ee99 = log(left) + p2;
      ee100 = ee50 * ee36;
      ee103 = ee93 * ee4 + ee97/ee40 + left * ee89 * ee2/ee32;
      ee106 = ee96 * ee4 + 3/ee26 + 4 * (ee3/ee25);
      ee107 = 1/ee32;
      ee108 = ee1 * ee67;
      ee109 = ee100 * ee99;
      ee110 = ee103 * ee36;
      ee113 = ee106 * ee4 + ee107 + 2 * (ee3/ee31);
      ee114 = ee108 - ee36;
      ee115 = ee50 * ee114;
      ee116 = ee109 - ee110;
      ee117 = ee113 * ee36;
      ee118 = (((((1/ee12 + 2 * ee53) * ee4 + 1/ee15 + 2 * (left *  ee29 * ee2)) * ee4 + 1/ee18 + 2 * (left * ee35 * ee2)) *  ee4 + 1/ee22 + 2 * (left * ee39 * ee2)) * ee4 + 1/ee26 +  2 * (left * ee43 * ee2)) * ee4;
      ee119 = ee115 + ee117;
      ee123 = ee46 * (2 * ee3 - ee4) + ee118 + (1 - ee3)/ee32 -  ee47;
      ee125 = ee116 * ee52/ee80;
      ee127 = ee123 * ee36 + ee50 * ee1 * ee67;
      ee128 = Rf_digamma(ee1);
      ee129 = log(y);
      
      out(j, 0) += -((ee125 + ee129 + p2 - ee128) * ee1);
      out(j, 1) += ee2 * (y - left * ee119 * ee52/ee80) - ee1;
      out(j, 2) += -((((R_pow(ee116, 2) * ee52/ee80 + (ee109 - 2 *
        ee110) * ee99 - ((((((ee4 * (Rf_trigamma(ee5) - R_pow(ee55, 2))/
          ee10 + (Rf_trigamma(ee7) - R_pow(ee57, 2))/ee13 + ee3 *
            (Rf_trigamma(ee6) - R_pow(ee54, 2))/ee12) * ee4 + (Rf_trigamma(ee9) -
            R_pow(ee62, 2))/ee17 + ee3 * (Rf_trigamma(ee8) -
            R_pow(ee56, 2))/ee15) * ee4 + (Rf_trigamma(ee14) - R_pow(ee71, 2))/
              ee21 + ee3 * (Rf_trigamma(ee11) - R_pow(ee61, 2))/ee18) *
                ee4 + (Rf_trigamma(ee19) - R_pow(ee82, 2))/ee25 + ee3 *
                (Rf_trigamma(ee16) - R_pow(ee70, 
                             2))/ee22) * ee4 + (Rf_trigamma(ee23) -
                               R_pow(ee90, 2))/ee31 + ee3 * (Rf_trigamma(ee20) -
                               R_pow(ee81, 2))/ee26) * ee4 + (Rf_trigamma(ee30) - R_pow(ee97, 2))/
                                 ee40 + ee3 * (Rf_trigamma(ee24) - R_pow(ee89, 2))/
                                   ee32) * ee36) * ee52/ee80 - Rf_trigamma(ee1)) * ee1 + ee125 +
                                     ee129 + p2 - ee128) * ee1);
      out(j, 3) += -(((ee100 + left * (ee127 * ee116 * ee52/ee80 +
        (ee115 + (ee118 + ee107 + 2 * (left * ee46 * ee2)) * ee36) *
        ee99 + ee103 * (ee36 - ee108) - ((((((2 * (left * ee55 * ee2/
          ee10) + ee54/ee12) * ee4 + 2 * (left * ee60 * ee2) + ee56/
            ee15) * ee4 + 2 * (left * ee66 * ee2) + ee61/ee18) * ee4 +
              2 * (left * ee76 * ee2) + ee70/ee22) * ee4 + 2 * (left * ee85 *
              ee2) + ee81/ee26) * ee4 + 2 * (left * ee93 * ee2) + ee89/
                ee32) * ee36) * ee2) * ee52/ee80 + 1) * ee1);
      out(j, 4) += ee2 * (y - left * (ee119 + left * (ee127 * ee119 *
        ee52/ee80 + (ee50 * (ee63 * R_pow(ee3, (ee1 - 2)) - ee67) +
        ee113 * ee67) * ee1 + ee123 * ee114 + (((((10/ee13 + 12 *
        ee27 + 2 * (left * ee69 * ee2)) * ee4 + 2 * (left * ee79 *
        ee2) + 8/ee17) * ee4 + 2 * (left * ee88 * ee2) + 6/ee21) *
        ee4 + 2 * (left * ee96 * ee2) + 4/ee25) * ee4 + 2 * (left *
        ee106 * ee2) + 2/ee31) * ee36 - ee117) * ee2) * ee52/ee80);
      
}

}

return out;

}

// //' @rdname ltgammad0
// [[Rcpp::export]]
arma::mat ltgammad34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat leftmat)
{

  int nobs = nhere.size();
  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);

  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);

  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
    p2vec = p2vec.elem(dupid);
  }

  double y, left, p1, p2;

  double ee1, ee2, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
  double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  double ee20, ee21, ee22, ee23, ee24, ee25, ee27, ee28, ee29;
  double ee30, ee31, ee32, ee34, ee35, ee36, ee37, ee38;
  double ee40, ee41, ee42, ee43, ee44, ee45, ee46, ee48, ee49;
  double ee50, ee51, ee52, ee54, ee56, ee57, ee58, ee59;
  double ee62, ee64, ee65, ee66, ee67, ee68, ee69;
  double ee70, ee72, ee73, ee75, ee76, ee77, ee78, ee79;
  double ee80, ee82, ee83, ee84, ee86, ee87, ee88, ee89;
  double ee90, ee92, ee93, ee95, ee97, ee98, ee99;
  double ee100, ee102, ee103, ee104, ee106, ee108;
  double ee110, ee111, ee112, ee113, ee115, ee118;
  double ee120, ee121, ee122, ee123, ee124, ee127, ee128, ee129;
  double ee131, ee132, ee133, ee134, ee135, ee136, ee137, ee138, ee139;
  double ee140, ee141, ee142, ee143, ee144, ee145, ee146;
  double ee150, ee153, ee155, ee156, ee157, ee158, ee159;
  double ee160, ee161, ee162, ee163, ee164, ee166, ee167;
  double ee171, ee174, ee175, ee176, ee177, ee178, ee179;
  double ee180, ee181, ee183, ee186, ee187, ee188, ee189;
  double ee190, ee191, ee192, ee193, ee194, ee196, ee197;
  double ee200, ee202, ee203, ee204, ee205, ee206, ee207, ee208, ee209;
  double ee210, ee212, ee215, ee217, ee218, ee219;
  double ee220, ee221, ee222, ee223, ee224, ee225, ee228, ee229;
  double ee231, ee232, ee233, ee234, ee235, ee236, ee237;
  double ee240, ee242, ee244, ee245, ee246, ee249;
  double ee250, ee251, ee252, ee253, ee254, ee255, ee256, ee257;
  double ee260, ee261, ee262, ee264, ee265, ee266, ee267;
  double ee270, ee271, ee272, ee273, ee274, ee275, ee276, ee277, ee279;
  double ee283, ee287, ee288, ee289;
  double ee290, ee291, ee292, ee293, ee294, ee295, ee297;
  
  for (int j=0; j < nobs; j++) {

    p1 = p1vec[j];
    p2 = p2vec[j];

    for (int l=0; l < nhere[j]; l++) {

      y = ymat(j, l);
      left = leftmat(j, l);

      ee1 = exp(p1);
      ee2 = exp(p2);
      ee3 = left * ee2;
      ee4 = R_pow(ee3, 2);
      ee5 = 13 + ee1;
      ee6 = 12 + ee1;
      ee7 = 11 + ee1;
      ee8 = tgamma(ee5);
      ee9 = 10 + ee1;
      ee10 = tgamma(ee6);
      ee11 = 9 + ee1;
      ee12 = tgamma(ee7);
      ee13 = 8 + ee1;
      ee14 = tgamma(ee9);
      ee15 = 7 + ee1;
      ee16 = tgamma(ee11);
      ee17 = 6 + ee1;
      ee18 = tgamma(ee13);
      ee19 = ee4/ee8;
      ee20 = ee19 + 1/ee12;
      ee21 = ee20 + ee3/ee10;
      ee22 = 5 + ee1;
      ee23 = tgamma(ee15);
      ee24 = 4 + ee1;
      ee25 = tgamma(ee17);
      ee27 = ee21 * ee4 + 1/ee16;
      ee28 = ee27 + ee3/ee14;
      ee29 = 3 + ee1;
      ee30 = tgamma(ee22);
      ee31 = 2 + ee1;
      ee32 = tgamma(ee24);
      ee34 = ee28 * ee4 + 1/ee23;
      ee35 = ee34 + ee3/ee18;
      ee36 = 1 + ee1;
      ee37 = tgamma(ee29);
      ee38 = tgamma(ee31);
      ee40 = ee35 * ee4 + 1/ee30;
      ee41 = R_pow(ee3, ee1);
      ee42 = ee40 + ee3/ee25;
      ee43 = ee3/ee8;
      ee44 = tgamma(ee36);
      ee45 = Rf_digamma(ee5);
      ee46 = Rf_digamma(ee6);
      ee48 = ee42 * ee4 + 1/ee37;
      ee49 = ee48 + ee3/ee32;
      ee50 = Rf_digamma(ee7);
      ee51 = Rf_digamma(ee9);
      ee52 = 1/ee44;
      ee54 = ee45 * ee4/ee8;
      ee56 = ee50/ee12 + ee54;
      ee57 = ee56 + left * ee46 * ee2/ee10;
      ee58 = Rf_digamma(ee11);
      ee59 = Rf_digamma(ee13);
      ee62 = ee49 * ee4 + ee52 + ee3/ee38;
      ee64 = ee57 * ee4 + ee58/ee16;
      ee65 = ee64 + left * ee51 * ee2/ee14;
      ee66 = Rf_digamma(ee15);
      ee67 = Rf_digamma(ee17);
      ee68 = 11/ee10;
      ee69 = ee68 + 12 * ee43;
      ee70 = 2 * ee43;
      ee72 = 1/ee10 + ee70;
      ee73 = ee1 - 1;
      ee75 = ee65 * ee4 + ee66/ee23;
      ee76 = ee72 * ee4;
      ee77 = 2 * (left * ee21 * ee2);
      ee78 = Rf_digamma(ee22);
      ee79 = ee75 + left * ee59 * ee2/ee18;
      ee80 = Rf_digamma(ee24);
      ee82 = ee76 + 1/ee14 + ee77;
      ee83 = ee69 * ee4;
      ee84 = 9/ee14;
      ee86 = ee83 + 10 * (ee3/ee12) + ee84;
      ee87 = R_pow(ee3, ee73);
      ee88 = ee82 * ee4;
      ee89 = 2 * (left * ee28 * ee2);
      ee90 = Rf_digamma(ee29);
      ee92 = ee79 * ee4 + ee78/ee30;
      ee93 = Rf_digamma(ee31);
      ee95 = ee92 + left * ee67 * ee2/ee25;
      ee97 = ee88 + 1/ee18 + ee89;
      ee98 = exp(-ee3);
      ee99 = ee97 * ee4;
      ee100 = 2 * (left * ee35 * ee2);
      ee102 = ee86 * ee4 + 7/ee18;
      ee103 = ee102 + 8 * (ee3/ee16);
      ee104 = Rf_digamma(ee36);
      ee106 = ee95 * ee4 + ee90/ee37;
      ee108 = ee99 + 1/ee25 + ee100;
      ee110 = log(left) + p2;
      ee111 = ee106 + left * ee80 * ee2/ee32;
      ee112 = ee108 * ee4;
      ee113 = 2 * (left * ee42 * ee2);
      ee115 = ee112 + 1/ee32 + ee113;
      ee118 = ee111 * ee4 + ee104/ee44 + left * ee93 * ee2/ee38;
      ee120 = ee103 * ee4 + 5/ee25;
      ee121 = ee120 + 6 * (ee3/ee23);
      ee122 = ee115 * ee4;
      ee123 = 1/ee38;
      ee124 = ee62 * ee41;
      ee127 = ee118 * ee41;
      ee128 = 1 - ee62 * ee98 * ee41;
      ee129 = 2 * ee3;
      ee131 = ee121 * ee4 + 3/ee32;
      ee132 = ee1 * ee87;
      ee133 = ee124 * ee110;
      ee134 = ee131 + 4 * (ee3/ee30);
      ee135 = 2 * (left * ee45 * ee2/ee8);
      ee136 = ee46/ee10;
      ee137 = ee135 + ee136;
      ee138 = ee137 * ee4;
      ee139 = 2 * (left * ee57 * ee2);
      ee140 = ee129 - ee4;
      ee141 = ee51/ee14;
      ee142 = R_pow(ee46, 2);
      ee143 = R_pow(ee45, 2);
      ee144 = Rf_trigamma(ee6);
      ee145 = Rf_trigamma(ee5);
      ee146 = ee133 - ee127;
      ee150 = ee49 * ee140 + ee122 + (1 - ee3)/ee38 - ee52;
      ee153 = ee134 * ee4 + ee123 + 2 * (ee3/ee37);
      ee155 = ee138 + ee139 + ee141;
      ee156 = 12 * ee19;
      ee157 = 2 * (left * ee69 * ee2);
      ee158 = R_pow(ee51, 2);
      ee159 = R_pow(ee50, 2);
      ee160 = ee132 - ee41;
      ee161 = Rf_trigamma(ee9);
      ee162 = Rf_trigamma(ee7);
      ee163 = ee144 - ee142;
      ee164 = ee145 - ee143;
      ee166 = ee62 * ee1 * ee87;
      ee167 = ee155 * ee4;
      ee171 = ee4 * ee164/ee8 + (ee162 - ee159)/ee12 + ee3 * ee163/ee10;
      ee174 = 10/ee12 + ee156 + ee157;
      ee175 = 2 * (left * ee65 * ee2);
      ee176 = ee59/ee18;
      ee177 = R_pow(ee59, 2);
      ee178 = R_pow(ee58, 2);
      ee179 = ee161 - ee158;
      ee180 = Rf_trigamma(ee13);
      ee181 = Rf_trigamma(ee11);
      ee183 = ee167 + ee175 + ee176;
      ee186 = ee171 * ee4 + (ee181 - ee178)/ee16 + ee3 * ee179/ee14;
      ee187 = ee174 * ee4;
      ee188 = 2 * (left * ee86 * ee2);
      ee189 = R_pow(ee67, 2);
      ee190 = R_pow(ee66, 2);
      ee191 = ee1 - 2;
      ee192 = Rf_trigamma(ee17);
      ee193 = Rf_trigamma(ee15);
      ee194 = ee180 - ee177;
      ee196 = ee150 * ee41 + ee166;
      ee197 = ee183 * ee4;
      ee200 = ee186 * ee4 + (ee193 - ee190)/ee23 + ee3 * ee194/ee18;
      ee202 = ee187 + ee188 + 8/ee16;
      ee203 = R_pow(ee3, ee191);
      ee204 = 2 * (left * ee79 * ee2);
      ee205 = R_pow(ee80, 2);
      ee206 = R_pow(ee78, 2);
      ee207 = ee67/ee25;
      ee208 = Rf_trigamma(ee24);
      ee209 = Rf_trigamma(ee22);
      ee210 = ee192 - ee189;
      ee212 = ee122 + ee123 + 2 * (left * ee49 * ee2);
      ee215 = ee200 * ee4 + (ee209 - ee206)/ee30 + ee3 * ee210/ee25;
      ee217 = ee197 + ee204 + ee207;
      ee218 = ee202 * ee4;
      ee219 = 2 * (left * ee103 * ee2);
      ee220 = R_pow(ee93, 2);
      ee221 = R_pow(ee90, 2);
      ee222 = Rf_trigamma(ee31);
      ee223 = Rf_trigamma(ee29);
      ee224 = ee208 - ee205;
      ee225 = ee62 * ee160;
      ee228 = ee215 * ee4 + (ee223 - ee221)/ee37 + ee3 * ee224/ee32;
      ee229 = ee217 * ee4;
      ee231 = ee218 + ee219 + 6/ee23;
      ee232 = 2 * (left * ee95 * ee2);
      ee233 = R_pow(ee104, 2);
      ee234 = ee80/ee32;
      ee235 = Rf_trigamma(ee36);
      ee236 = ee222 - ee220;
      ee237 = ee153 * ee41;
      ee240 = ee228 * ee4 + (ee235 - ee233)/ee44 + ee3 * ee236/ee38;
      ee242 = ee229 + ee232 + ee234;
      ee244 = ee231 * ee4 + 2 * (left * ee121 * ee2);
      ee245 = ee240 * ee41;
      ee246 = ee212 * ee41;
      ee249 = ee242 * ee4 + 2 * (left * ee111 * ee2) + ee93/ee38;
      ee250 = ee244 + 4/ee30;
      ee251 = ee73 * ee203;
      ee252 = 2 * ee127;
      ee253 = 2 * ee41;
      ee254 = (ee133 - ee252) * ee110;
      ee255 = ee249 * ee41;
      ee256 = ee250 * ee4;
      ee257 = 4 * ee3;
      ee260 = R_pow(ee146, 2) * ee98/ee128;
      ee261 = ee225 + ee237;
      ee262 = ee153 * ee87;
      ee264 = ee256 + 2 * (left * ee134 * ee2) + 2/ee37;
      ee265 = ((((2 * (ee20 + left * (ee70 + 2/ee10) * ee2) +  2 * ee19 + 2 * (left * ee72 * ee2)) * ee4 + 2 * (ee27 +  left * (ee76 + ee77 + 2/ee14) * ee2) + 2 * (left * ee82 *  ee2)) * ee4 + 2 * (ee34 + left * (ee88 + ee89 + 2/ee18) *  ee2) + 2 * (left * ee97 * ee2)) * ee4 + 2 * (ee40 + left *  (ee99 + ee100 + 2/ee25) * ee2) + 2 * (left * ee108 *  ee2)) * ee4;
      ee266 = ee251 - ee87;
      ee267 = 2 - ee129;
      ee270 = ee196 * ee146 * ee98/ee128;
      ee271 = (ee62 * ee266 + ee262) * ee1;
      ee272 = (ee225 + ee246) * ee110;
      ee273 = ee260 + ee254;
      ee274 = ee118 * (ee41 - ee132);
      ee275 = ee264 * ee41;
      ee276 = 2 * ee132;
      ee277 = 2 * ee4;
      ee279 = (ee150 * ee87 + ee62 * ee73 * ee203 + ee212 * ee87) *  ee1 + (ee49 * ee267 + ee115 * (ee257 - ee4) + ee265 -  ee123) * ee41;
      ee283 = ee196 * ee261 * ee98/ee128 + ee271;
      ee287 = ee273 - ee245;
      ee288 = (ee166 + ee246) * ee110;
      ee289 = ee255 + ee118 * ee1 * ee87;
      ee290 = ee41 + ee253;
      ee291 = ee276 - ee253;
      ee292 = 2 * ee87;
      ee293 = 2/ee44;
      ee294 = 24 * ee43;
      ee295 = ee257 - ee277;
      ee297 = left * (ee270 + ee272 + ee274 - ee255) * ee2;
      
      out(j, 0) += -((((((ee146 * ee110 - (ee245 + 2 * (ee245 + ee127 *
        ee110))) * ee110 + (ee254 + 2 * ee260 + 2 * (ee254 - ee245) -
        ee245) * ee146 * ee98/ee128 - ((((((ee4 * (Rf_psigamma(ee5, 2) -
        (3 * ee145 - ee143) * ee45)/ee8 + (Rf_psigamma(ee7, 2) -
        (3 * ee162 - ee159) * ee50)/ee12 + ee3 * (Rf_psigamma(ee6, 2) -
        (3 * ee144 - ee142) * ee46)/ee10) * ee4 + (Rf_psigamma(ee11, 2) -
        (3 * ee181 - ee178) * ee58)/ee16 + ee3 *
        (Rf_psigamma(ee9, 2) - (3 * ee161 - ee158) * ee51)/ee14) *
        ee4 + (Rf_psigamma(ee15, 2) - (3 * ee193 - ee190) * ee66)/ee23 +
        ee3 * (Rf_psigamma(ee13, 2) - (3 * ee180 - ee177) * ee59)/
          ee18) * ee4 + (Rf_psigamma(ee22, 2) - (3 * ee209 - ee206) *
            ee78)/ee30 + ee3 * (Rf_psigamma(ee17, 2) - (3 * ee192 - ee189) *
            ee67)/ee25) * ee4 + (Rf_psigamma(ee29, 2) - (3 * ee223 -
            ee221) * ee90)/ee37 + ee3 * (Rf_psigamma(ee24, 2) - (3 *
            ee208 - ee205) * ee80)/ee32) * ee4 + (Rf_psigamma(ee36, 2) -
            (3 * ee235 - ee233) * ee104)/ee44 + ee3 * (Rf_psigamma(ee31, 2) -
            (3 * ee222 - ee220) * ee93)/ee38) * ee41) * ee98/ee128 -
            Rf_psigamma(ee1, 2)) * ee1 + 3 * (ee287 * ee98/ee128) -
            3 * Rf_trigamma(ee1)) * ee1 + ee146 * ee98/ee128 + log(y) +
            p2 - Rf_digamma(ee1)) * ee1);
      out(j, 1) += -(((((ee146 * (2 * (ee124 + left * (ee288 - ee289) *
        ee2) + left * (ee196 * ee98/ee128 - 1) * ee146 * ee2) +
        left * ee287 * ee196 * ee2) * ee98/ee128 + (2 * ee124 + left *
        (ee288 - 2 * ee289) * ee2) * ee110 - (ee252 + left * (ee273 +
        ee240 * ee160 + ((((((ee163/ee10 + 2 * (ee3 * ee164/
          ee8)) * ee4 + ee179/ee14 + 2 * (left * ee171 * ee2)) * ee4 +
            ee194/ee18 + 2 * (left * ee186 * ee2)) * ee4 + ee210/ee25 +
            2 * (left * ee200 * ee2)) * ee4 + ee224/ee32 + 2 * (left * ee215 *
            ee2)) * ee4 + ee236/ee38 + 2 * (left * ee228 * ee2)) *
            ee41) * ee2)) * ee1 + ee124 + ee297) * ee98/ee128 + 1) * ee1);
      out(j, 2) += -(left * ((ee196 * (ee62 * (ee41 * ee110 + ee253) +
        2 * ee297 - ee127) + left * ee279 * ee146 * ee2) * ee98/
          ee128 + ee272 + ee62 * ee291 + ee274 + 2 * ee246 + left * ((ee62 *
            ((ee251 - ee292) * ee1 + ee41) + ee212 * ee291 + (ee265 +
            2 * (ee48 + left * (ee112 + ee113 + 2/ee32) * ee2) + 2 *
            (left * ee115 * ee2)) * ee41) * ee110 + ee249 * (ee253 -
            ee276) + ee118 * ((ee292 - ee251) * ee1 - ee41) - (ee270 + (((((2 *
            (ee56 + left * (2 * ee136 + ee135) * ee2) + 2 * ee54 +
            2 * (left * ee137 * ee2)) * ee4 + 2 * (ee64 + left * (ee138 +
            2 * ee141 + ee139) * ee2) + 2 * (left * ee155 * ee2)) *
            ee4 + 2 * (ee75 + left * (ee167 + 2 * ee176 + ee175) * ee2) +
            2 * (left * ee183 * ee2)) * ee4 + 2 * (ee92 + left * (ee197 +
            2 * ee207 + ee204) * ee2) + 2 * (left * ee217 * ee2)) *
            ee4 + 2 * (ee106 + left * (ee229 + 2 * ee234 + ee232) * ee2) +
            2 * (left * ee242 * ee2)) * ee41)) * ee2 - ee255) * ee98 *
            ee1 * ee2/ee128);
      out(j, 3) += ee2 * (y - left * (ee261 + left * (ee196 * (ee153 *
        ee290 + 3 * ee225 + left * (ee283 + ee150 * ee160 + ee275 -
        ee237) * ee2) * ee98/ee128 + (ee49 * (6 * ee3 - 3 * ee4) +
        (3 - 3 * ee3)/ee38 + 3 * ee122 - 3/ee44) * ee160 + (ee134 *
        ee140 + ee256 + ee267/ee37 - ee123) * ee290 + 3 * ee271 +
        left * ((ee279 * ee261 + (ee283 + (ee49 * (ee129 - ee277) +
        ee122 + (1 - ee129)/ee38 - ee293) * ee160 + ee275 - 2 * ee237) *
        ee196) * ee98/ee128 + ((ee62 * (ee191 * R_pow(ee3, (ee1 -
        3)) - ee203) + ee153 * ee203) * 
        ee73 + (ee49 * ee295 +
        ee267/ee38 + 2 * ee122 - ee293) * ee266 + 2 * (ee264 *
        ee87) - 2 * ee262) * ee1 + (ee49 * (ee4 + 2 - ee257) + ee115 *
        ee295 + ee265 + ee52 - (2 - ee3)/ee38) * ee160 + ee237 +
        (((((2 * (ee68 + ee294) + ee294) * ee4 + 2 * (ee83 + ee84 +
        left * (ee156 + ee157 + 20/ee12) * ee2) + 2 * (left * ee174 *
        ee2)) * ee4 + 2 * (ee102 + left * (ee187 + 16/ee16 + ee188) *
        ee2) + 2 * (left * ee202 * ee2)) * ee4 + 2 * (ee120 + left *
        (ee218 + 12/ee23 + ee219) * ee2) + 2 * (left * ee231 *
        
        ee2)) * ee4 + 2 * (ee131 + left * (ee244 + 8/ee30) * ee2) +
        2 * (left * ee250 * ee2)) * ee41 - 2 * ee275) * ee2) *
        ee2) * ee98/ee128);
      
    }

  }

  return out;

}

// //' Left-truncated fixed-alpha Gamma distribution negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each GEV parameter
// //' @param X1 a design matrix for the shape parameter
// //' @param ymat a matrix
// //' @param leftmat a matrix
// //' @param alphamat a matrix
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return ltgammabd0 a scalar, the negative log-likelihood
// //' @return ltgammab12 a matrix, first then second derivatives w.r.t. parameters
// //' @return ltgammab34 a matrix, third then fourth derivatives w.r.t. parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double ltgammabd0(Rcpp::List pars, arma::mat X1, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat leftmat, arma::mat alphamat)
{
  
  int nobs = nhere.size();
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
  }
  
  double y, left, p1, alpha, beta, lnum;
  
  double nllh = 0.0;
  
  for (int j=0; j < nobs; j++) {
    
    p1 = p1vec[j];
    beta = exp(p1);
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      left = leftmat(j, l);
      alpha = alphamat(j, l);
      
      lnum = (alpha - 1) * log(y) - beta * y + alpha * p1 - lgamma(alpha);
      nllh += ldenom(left, alpha, beta) - lnum;      
    }
    
  }
  
  return(nllh);
  
}

// //' @rdname ltgammabd0
// [[Rcpp::export]]
arma::mat ltgammabd12(Rcpp::List pars, arma::mat X1, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat leftmat, arma::mat alphamat)
{
  
  int nobs = nhere.size();
  arma::mat out = arma::mat(nobs, 2, arma::fill::zeros);
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
  }
  
  double y, left, alpha, p1;
  
  double ee1, ee2, ee3, ee5, ee7, ee9;
  double ee11, ee13, ee15, ee17, ee19;
  double ee20, ee22, ee24, ee26, ee29;
  double ee31, ee33, ee36, ee37;
  double ee40, ee41, ee44, ee46, ee48;
  double ee51, ee54, ee57, ee58;
  double ee61, ee62, ee66, ee67;
  double ee70, ee71, ee73, ee75, ee78;
  double ee82;
  
  for (int j=0; j < nobs; j++) {
    
    p1 = p1vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      left = leftmat(j, l);
      alpha = alphamat(j, l);
      
      ee1 = exp(p1);
      ee2 = left * ee1;
      ee3 = R_pow(ee2, 2);
      ee5 = tgamma(13 + alpha);
      ee7 = tgamma(12 + alpha);
      ee9 = tgamma(11 + alpha);
      ee11 = tgamma(10 + alpha);
      ee13 = tgamma(9 + alpha);
      ee15 = tgamma(8 + alpha);
      ee17 = tgamma(7 + alpha);
      ee19 = tgamma(6 + alpha);
      ee20 = ee3/ee5;
      ee22 = tgamma(5 + alpha);
      ee24 = ee20 + 1/ee9 + ee2/ee7;
      ee26 = tgamma(4 + alpha);
      ee29 = ee24 * ee3 + 1/ee13 + ee2/ee11;
      ee31 = tgamma(3 + alpha);
      ee33 = tgamma(2 + alpha);
      ee36 = ee29 * ee3 + 1/ee17 + ee2/ee15;
      ee37 = R_pow(ee2, alpha);
      ee40 = ee36 * ee3 + 1/ee22 + ee2/ee19;
      ee41 = ee2/ee5;
      ee44 = ee40 * ee3 + 1/ee31 + ee2/ee26;
      ee46 = 1/tgamma(1 + alpha);
      ee48 = 11/ee7 + 12 * ee41;
      ee51 = ee48 * ee3 + 10 * (ee2/ee9) + 9/ee11;
      ee54 = ee44 * ee3 + ee46 + ee2/ee33;
      ee57 = ee51 * ee3 + 7/ee15 + 8 * (ee2/ee13);
      ee58 = alpha - 1;
      ee61 = ee57 * ee3 + 5/ee19 + 6 * (ee2/ee17);
      ee62 = R_pow(ee2, ee58);
      ee66 = ee61 * ee3 + 3/ee26 + 4 * (ee2/ee22);
      ee67 = exp(-ee2);
      ee70 = ee66 * ee3 + 1/ee33 + 2 * (ee2/ee31);
      ee71 = ee70 * ee37;
      ee73 = alpha * ee62 - ee37;
      ee75 = ee54 * ee73 + ee71;
      ee78 = 1 - ee54 * ee67 * ee37;
      ee82 = ee44 * (2 * ee2 - ee3) + (((((1/ee7 + 2 * ee41) *  ee3 + 1/ee11 + 2 * (left * ee24 * ee1)) * ee3 + 1/ee15 +  2 * (left * ee29 * ee1)) * ee3 + 1/ee19 + 2 * (left *  ee36 * ee1)) * ee3 + 1/ee26 + 2 * (left * ee40 * ee1)) *  ee3 + (1 - ee2)/ee33 - ee46;
      
      out(j, 0) += ee1 * (y - left * ee75 * ee67/ee78) - alpha;
      out(j, 1) += ee1 * (y - left * (ee75 + left * ((ee82 * ee37 +
        alpha * ee54 * ee62) * ee75 * ee67/ee78 + ee82 * ee73 + (((((10/
          ee9 + 12 * ee20 + 2 * (left * ee48 * ee1)) * ee3 + 2 *
            (left * ee51 * ee1) + 8/ee13) * ee3 + 2 * (left * ee57 * ee1) +
            6/ee17) * ee3 + 2 * (left * ee61 * ee1) + 4/ee22) * ee3 +
            2 * (left * ee66 * ee1) + 2/ee31) * ee37 + alpha * (ee54 *
            (ee58 * R_pow(ee2, (alpha - 2)) - ee62) + ee70 * ee62) - ee71) *
            ee1) * ee67/ee78);
      
    }
    
  }
  
  return out;
  
}

// //' @rdname ltgammabd0
// [[Rcpp::export]]
arma::mat ltgammabd34(Rcpp::List pars, arma::mat X1, arma::mat ymat, arma::uvec dupid, int dcate, arma::ivec nhere, arma::mat leftmat, arma::mat alphamat)
{
  
  int nobs = nhere.size();
  arma::mat out = arma::mat(nobs, 2, arma::fill::zeros);
  
  arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
  
  if (dcate == 1) {
    p1vec = p1vec.elem(dupid);
  }
  
  double y, left, alpha, p1;

  double ee1, ee2, ee3, ee5, ee7, ee9;
  double ee11, ee12, ee14, ee15, ee16, ee17, ee19;
  double ee21, ee23, ee24, ee25, ee26, ee28;
  double ee30, ee31, ee32, ee34, ee36, ee37, ee39;
  double ee40, ee42, ee43, ee44, ee46, ee48, ee49;
  double ee50, ee52, ee53, ee54, ee55, ee56, ee57, ee59;
  double ee60, ee61, ee63, ee64, ee65, ee66, ee67, ee68;
  double ee70, ee71, ee72, ee73, ee75, ee76, ee77, ee78;
  double ee80, ee81, ee82, ee83, ee84, ee85, ee86, ee88;
  double ee91, ee92, ee93, ee94, ee95, ee96, ee97, ee99;
  double ee102, ee103, ee105, ee106, ee107, ee108, ee109;
  double ee110, ee111, ee112, ee114, ee115, ee117, ee118, ee119;
  double ee121, ee122, ee123, ee124, ee126, ee127, ee128, ee129;
  double ee133, ee136, ee137, ee139;
  double ee140, ee141, ee142, ee145, ee146, ee147, ee148;
  double ee150, ee153, ee155, ee156, ee158, ee159;
  double ee160, ee161, ee165, ee168, ee169;
  double ee170, ee171, ee172, ee174, ee175, ee176, ee177, ee179;
  double ee180, ee181, ee182, ee184, ee186, ee187, ee188, ee189;
  double ee191, ee192, ee193, ee194, ee195, ee197, ee198, ee199;
  double ee200, ee201, ee203, ee205, ee207, ee208;
  double ee210, ee211, ee212, ee213, ee217, ee218;
  double ee220, ee221, ee223, ee224, ee225, ee226, ee227, ee228, ee229;
  double ee230, ee231, ee232, ee233, ee234, ee235;
  double ee243, ee244, ee245, ee248;
  double ee250, ee251, ee252, ee253, ee254, ee255, ee256, ee257;
  double ee265, ee267, ee268;
  double ee272, ee276, ee277, ee278, ee279;
  double ee280, ee282, ee283, ee284, ee286, ee287, ee288, ee289;
  double ee290, ee291, ee292, ee294;
  
  for (int j=0; j < nobs; j++) {
    
    p1 = p1vec[j];
    
    for (int l=0; l < nhere[j]; l++) {
      
      y = ymat(j, l);
      left = leftmat(j, l);
      alpha = alphamat(j, l);
      
      ee1 = exp(p1);
      ee2 = left * ee1;
      ee3 = R_pow(ee2, 2);
      ee5 = tgamma(13 + alpha);
      ee7 = tgamma(12 + alpha);
      ee9 = tgamma(11 + alpha);
      ee11 = tgamma(10 + alpha);
      ee12 = ee3/ee5;
      ee14 = tgamma(9 + alpha);
      ee15 = ee2/ee5;
      ee16 = ee12 + 1/ee9;
      ee17 = ee16 + ee2/ee7;
      ee19 = tgamma(8 + alpha);
      ee21 = tgamma(7 + alpha);
      ee23 = ee17 * ee3 + 1/ee14;
      ee24 = 11/ee7;
      ee25 = ee23 + ee2/ee11;
      ee26 = ee24 + 12 * ee15;
      ee28 = tgamma(6 + alpha);
      ee30 = tgamma(5 + alpha);
      ee31 = ee26 * ee3;
      ee32 = 9/ee11;
      ee34 = ee25 * ee3 + 1/ee21;
      ee36 = ee31 + 10 * (ee2/ee9) + ee32;
      ee37 = ee34 + ee2/ee19;
      ee39 = tgamma(4 + alpha);
      ee40 = 2 * ee15;
      ee42 = ee36 * ee3 + 7/ee19;
      ee43 = 1/ee7;
      ee44 = ee43 + ee40;
      ee46 = tgamma(3 + alpha);
      ee48 = ee37 * ee3 + 1/ee30;
      ee49 = ee42 + 8 * (ee2/ee14);
      ee50 = ee48 + ee2/ee28;
      ee52 = tgamma(2 + alpha);
      ee53 = ee44 * ee3;
      ee54 = 2 * (left * ee17 * ee1);
      ee55 = ee53 + 1/ee11;
      ee56 = ee55 + ee54;
      ee57 = R_pow(ee2, alpha);
      ee59 = ee49 * ee3 + 5/ee28;
      ee60 = alpha - 1;
      ee61 = ee56 * ee3;
      ee63 = ee50 * ee3 + 1/ee46;
      ee64 = ee59 + 6 * (ee2/ee21);
      ee65 = ee63 + ee2/ee39;
      ee66 = 2 * (left * ee25 * ee1);
      ee67 = 12 * ee12;
      ee68 = 2 * (left * ee26 * ee1);
      ee70 = tgamma(1 + alpha);
      ee71 = ee61 + 1/ee19;
      ee72 = ee71 + ee66;
      ee73 = R_pow(ee2, ee60);
      ee75 = 10/ee9 + ee67;
      ee76 = ee75 + ee68;
      ee77 = 1/ee70;
      ee78 = ee72 * ee3;
      ee80 = ee64 * ee3 + 3/ee39;
      ee81 = ee76 * ee3;
      ee82 = 2 * (left * ee37 * ee1);
      ee83 = 2 * (left * ee36 * ee1);
      ee84 = ee80 + 4 * (ee2/ee30);
      ee85 = 8/ee14;
      ee86 = ee78 + 1/ee28;
      ee88 = ee81 + ee83 + ee85;
      ee91 = ee65 * ee3 + ee77 + ee2/ee52;
      ee92 = ee86 + ee82;
      ee93 = 1/ee52;
      ee94 = ee88 * ee3;
      ee95 = 2 * (left * ee49 * ee1);
      ee96 = ee92 * ee3;
      ee97 = 2 * (left * ee50 * ee1);
      ee99 = ee96 + 1/ee39 + ee97;
      ee102 = ee84 * ee3 + ee93 + 2 * (ee2/ee46);
      ee103 = 6/ee21;
      ee105 = ee94 + ee95 + ee103;
      ee106 = 2 * ee2;
      ee107 = alpha - 2;
      ee108 = ee99 * ee3;
      ee109 = ee105 * ee3;
      ee110 = ee109 + 2 * (left * ee64 * ee1);
      ee111 = R_pow(ee2, ee107);
      ee112 = 4/ee30;
      ee114 = ee110 + ee112;
      ee115 = 24 * ee15;
      ee117 = alpha * ee73 - ee57;
      ee118 = exp(-ee2);
      ee119 = 2/ee7;
      ee121 = 2 * (ee16 + left * (ee40 + ee119) * ee1);
      ee122 = 2 * ee12;
      ee123 = 2 * (left * ee44 * ee1);
      ee124 = ee106 - ee3;
      ee126 = ee121 + ee122 + ee123;
      ee127 = ee114 * ee3;
      ee128 = ee102 * ee57;
      ee129 = 2/ee11;
      ee133 = ee65 * ee124 + ee108 + (1 - ee2)/ee52 - ee77;
      ee136 = ee126 * ee3 + 2 * (ee23 + left * (ee53 + ee54 +  ee129) * ee1);
      ee137 = 2 * (left * ee56 * ee1);
      ee139 = ee127 + 2 * (left * ee84 * ee1) + 2/ee46;
      ee140 = ee136 + ee137;
      ee141 = 2/ee19;
      ee142 = 4 * ee2;
      ee145 = ee140 * ee3 + 2 * (ee34 + left * (ee61 + ee66 +  ee141) * ee1);
      ee146 = 2 * (left * ee72 * ee1);
      ee147 = ee91 * ee117;
      ee148 = ee145 + ee146;
      ee150 = ee60 * ee111 - ee73;
      ee153 = 1 - ee91 * ee118 * ee57;
      ee155 = 2 * (ee24 + ee115);
      ee156 = alpha * ee91;
      ee158 = ee133 * ee57 + ee156 * ee73;
      ee159 = ee155 + ee115;
      ee160 = 2/ee28;
      ee161 = 20/ee9;
      ee165 = ee148 * ee3 + 2 * (ee48 + left * (ee78 + ee82 +  ee160) * ee1) + 2 * (left * ee92 * ee1);
      ee168 = ee159 * ee3;
      ee169 = 2 * (ee31 + ee32 + left * (ee67 + ee68 + ee161) *  ee1);
      ee170 = 2 * (left * ee76 * ee1);
      ee171 = ee147 + ee128;
      ee172 = ee102 * ee73;
      ee174 = ee168 + ee169 + ee170;
      ee175 = ee165 * ee3;
      ee176 = 16/ee14;
      ee177 = ee139 * ee57;
      ee179 = ee174 * ee3;
      ee180 = 2 * (ee42 + left * (ee81 + ee176 + ee83) * ee1);
      ee181 = 2 * (left * ee88 * ee1);
      ee182 = 2 * ee57;
      ee184 = ee91 * ee150 + ee172;
      ee186 = ee179 + ee180 + ee181;
      ee187 = 12/ee21;
      ee188 = 2 - ee106;
      ee189 = alpha * ee184;
      ee191 = ee186 * ee3;
      ee192 = ee57 + ee182;
      ee193 = 2 * (ee59 + left * (ee94 + ee187 + ee95) * ee1);
      ee194 = 2 * (left * ee105 * ee1);
      ee195 = 2 * ee3;
      ee197 = ee191 + ee193 + ee194;
      ee198 = alpha - 3;
      ee199 = ee158 * ee171;
      ee200 = R_pow(ee2, ee198);
      ee201 = 8/ee30;
      ee203 = ee199 * ee118/ee153;
      ee205 = ee108 + ee93 + 2 * (left * ee65 * ee1);
      ee207 = ee197 * ee3;
      ee208 = 2 * (ee80 + left * (ee110 + ee201) * ee1);
      ee210 = ee207 + ee208 + 2 * (left * ee114 * ee1);
      ee211 = 3 * ee3;
      ee212 = 6 * ee2;
      ee213 = ee102 * ee111;
      ee217 = ee65 * ee188 + ee99 * (ee142 - ee3) + ee175 - ee93;
      ee218 = ee139 * ee73;
      ee220 = ee107 * ee200 - ee111;
      ee221 = (ee91 * ee220 + ee213) * ee60;
      ee223 = ee217 * ee57 + alpha * (ee133 * ee73 + ee91 * ee60 *  ee111 + ee205 * ee73);
      ee224 = ee210 * ee57;
      ee225 = 2/ee70;
      ee226 = 4 * ee57;
      ee227 = 2 * ee128;
      ee228 = 2 * ee108;
      ee229 = 2 * ee218;
      ee230 = 2 * ee177;
      ee231 = 2 * ee73;
      ee232 = 3 * ee2;
      ee233 = 3/ee70;
      ee234 = ee142 - ee195;
      ee235 = ee223 * ee171;
      ee243 = ee84 * ee124 + ee127 + ee188/ee46 - ee93;
      ee244 = ee73 + ee231;
      ee245 = 4 * ee3;
      ee248 = left * (ee203 + ee133 * ee117 + ee177 + ee189 -  ee128) * ee1;
      ee250 = (ee203 + (ee65 * (ee106 - ee195) + ee108 + (1 -  ee106)/ee52 - ee225) * ee117 + ee177 + ee189 - ee227) *  ee158 + ee235;
      ee251 = ee243 * ee192;
      ee252 = ee192 + ee226;
      ee253 = 2 * ee172;
      ee254 = 2 * (left * ee99 * ee1);
      ee255 = ee142 - ee245;
      ee256 = 4/ee70;
      ee257 = ee212 - ee211;
      ee265 = ee102 * ee244;
      ee267 = ee102 * ee252 + 7 * ee147;
      ee268 = ee102 * ee192;
      ee272 = ee65 * ee255 + (2 - ee142)/ee52 + ee228 - ee256;
      ee276 = ee65 * ee257 + (3 - ee232)/ee52 + 3 * ee108 - ee233;
      ee277 = ee139 * ee192;
      ee278 = ee210 * ee73;
      ee279 = ee175 + ee254;
      ee280 = ((((2 * (ee43 + 4 * ee15) + 2 * (ee119 + 6 * ee15) +  4 * ee15) * ee3 + 2 * (ee55 + left * (ee121 + 2 * ee17 +  ee122 + ee123) * ee1) + 2 * (2 * ee53 + ee129 + left *  (ee126 + 4 * ee17) * ee1) + 2 * (left * ee126 * ee1)) *  ee3 + 2 * (ee71 + left * (ee136 + 2 * ee25 + ee137) *  ee1) + 2 * (2 * ee61 + ee141 + left * (ee140 + 4 * ee25) *  ee1) + 2 * (left * ee140 * ee1)) * ee3 + 2 * (ee86 +  left * (ee145 + 2 * ee37 + ee146) * ee1) + 2 * (2 * ee78 +  ee160 + left * (ee148 + 4 * ee37) * ee1) + 2 * (left *  ee148 * ee1)) * ee3;
      ee282 = 18 * ee2;
      ee283 = ee182 + ee226;
      ee284 = 2 + ee195;
      ee286 = ee211 + 6 - 12 * ee2;
      ee287 = 4 * ee128;
      ee288 = ee142 - ee211;
      ee289 = 7 * ee2;
      ee290 = 7 * ee3;
      ee291 = 7/ee70;
      ee292 = 8 * ee2;
      ee294 = left * (ee250 * ee118/ee153 + (ee65 * (ee3 + 2 -  ee142) + ee99 * ee234 + ee175 + ee77 - (2 - ee2)/ee52) *  ee117 + ee128 + ee224 + alpha * (ee221 + (ee65 * ee234 +  ee188/ee52 + ee228 - ee225) * ee150 + ee229 - ee253) -  ee230) * ee1;
      
      out(j, 0) += ee1 * (y - left * (ee171 + left * (ee158 * (ee268 +
        3 * ee147 + ee248) * ee118/ee153 + ee276 * ee117 + ee251 +
        3 * ee189 + ee294) * ee1) * ee118/ee153);
      out(j, 1) += ee1 * (y - left * (ee171 + left * (ee158 * (ee267 +
        left * (ee158 * (ee102 * (ee57 + ee226) + 5 * ee147 + ee248) *
        ee118/ee153 + (ee65 * (10 * ee2 - 5 * ee3) + (5 - 5 *
        ee2)/ee52 + 5 * ee108 - 5/ee70) * ee117 + ee251 + ee230 +
        5 * ee189 + ee294 - ee227) * ee1) * ee118/ee153 + (ee65 * (14 *
        ee2 - ee290) + (7 - ee289)/ee52 + 7 * ee108 - ee291) * ee117 +
        ee177 + 2 * ee251 + 7 * ee189 + left * (((ee158 * (ee267 +
        ee248) * ee118/ee153 + (ee65 * (ee292 - ee290) + (4 -
        ee289)/ee52 + 4 * ee108 - (6 * ee91 + 
        ee291)) * ee117 +
        ee102 * (alpha * ee244 - (ee192 + ee182 + 8 * ee57)) + ee139 *
        ee252 + 3 * (ee205 * ee117 + ee156 * ee150) + 4 * ee189 +
        left * (((ee203 + (ee65 * (ee106 - ee211) + ee108 + (1 - ee232)/
          ee52 - ee233) * ee117 + ee177 + ee189 - ee268) * ee158 +
            ee235) * ee118/ee153 + (ee65 * (ee284 - ee212) + ee99 * ee288 +
            ee175 + ee225 - (3 - ee106)/ee52) * ee117 + ee224 + ee227 +
            alpha * (ee221 + (ee65 * ee288 + (2 - ee232)/ee52 + ee228 -
            ee233) * ee150 + ee229 - ee265) - ee277) * ee1) * ee158 +
            
            ee223 * (ee102 * ee283 + 6 * ee147 + ee248)) * ee118/
              ee153 + (ee65 * (12 + 6 * ee3 - 24 * ee2) + ee99 * (ee282 -
                12 * ee3) + 3 * ee279 + 3 * ee175 + 6/ee70 - (12 - ee212)/
                  ee52) * ee117 + (ee84 * (ee284 - ee142) + ee114 * (ee212 -
                    ee195) + (ee142 - 4)/ee46 + ee208 + 2 * ee207 + 2/ee52) * ee192 +
                    alpha * ((ee65 * (ee282 - 9 * ee3) + (9 - 9 * ee2)/ee52 +
                    9 * ee108 - 9/ee70) * ee150 + ee243 * (ee244 + ee231 + 4 *
                    ee73) + 3 * (ee221 + ee205 * ee150 + ee218) + 3 * ee221 -
                    3 * ee184) + left * (((ee203 + 
                    ee272 * ee117 + ee230 +
                    2 * ee189 - ee287) * ee223 + ee158 * ((ee65 * (2 + ee245 -
                    ee292) + ee99 * ee255 + ee175 + ee256 - (4 - ee142)/ee52) *
                    ee117 + ee224 + (2 * ee250 - ee199) * ee118/ee153 + ee287 +
                    alpha * (ee221 + ee272 * ee150 + ee229 - ee253) - (2 * (ee177 +
                    alpha * ee102 * ee73) + ee230)) + ee171 * ((ee99 * (6 -
                    ee142) + ee165 * (ee212 - ee3) + ee280 - 2 * ee65) * ee57 +
                    alpha * ((ee133 * ee111 + ee91 * ee107 * ee200 + 2 * (ee205 *
                    ee111)) * ee60 + (ee175 + 2 * (ee63 + left * (ee96 + 
                    ee97 +
                    2/ee39) * ee1) + ee254) * ee73 + 2 * (ee217 * ee73)))) *
                    ee118/ee153 + (ee65 * (ee212 - (ee3 + 6)) + ee99 * ee286 +
                    ee165 * ee257 + ee280 + (3 - ee2)/ee52 - ee77) * ee117 +
                    ee277 + ((((2 * (ee75 + left * (2 * ee26 + ee155 + ee115) *
                    ee1) + 2 * (ee161 + 24 * ee12 + left * (ee159 + 4 * ee26) *
                    ee1) + 2 * (left * ee159 * ee1) + 72 * ee12) * ee3 + 2 * (ee81 +
                    ee85 + left * (ee168 + 2 * ee36 + ee169 + ee170) * ee1) +
                    2 * (ee176 + 2 * ee81 + left * (ee174 + 4 * ee36) * ee1) +
                    2 * (left * ee174 * 
                    ee1)) * ee3 + 2 * (ee94 + ee103 +
                    left * (ee179 + 2 * ee49 + ee180 + ee181) * ee1) + 2 * (ee187 +
                    2 * ee94 + left * (ee186 + 4 * ee49) * ee1) + 2 * (left *
                    ee186 * ee1)) * ee3 + 2 * (ee109 + ee112 + left * (ee191 +
                    2 * ee64 + ee193 + ee194) * ee1) + 2 * (2 * ee109 + ee201 +
                    left * (ee197 + 4 * ee64) * ee1) + 2 * (left * ee197 * ee1)) *
                    ee57 + alpha * (((ee91 * (ee198 * R_pow(ee2, (alpha - 4)) -
                    ee200) + ee102 * ee200) * ee107 + ee276 * ee220 + ee139 *
                    ee111 - ee213) * ee60 + (ee65 * ee286 + ee99 * 
                    (2 *
                    ee234 - ee195) + ee175 + 2 * ee279 + ee233 - (6 - ee232)/ee52) *
                    ee150 + ee265 + ee278 + 2 * (ee139 * ee60 * ee111 + ee278) -
                    (2 * (ee102 * ee60 * ee111 + ee218) + ee229)) - (ee128 +
                    ee224 + 2 * (ee224 + alpha * ee139 * ee73))) * ee1 - ee139 *
                    ee283) * ee1 - ee128) * ee1) * ee118/ee153);    }
    
  }
  
  return out;
  
}
