// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

double ltind2(double x, double u) {
  double out = 0.0;
  if (x < u) {
    out = 1.0;
  }
  return out;
}

 //' Asymmetric generalised Gaussian log-likelihood
 //'
 //' @param pars a list of vectors of coefficients for each conditional EVD parameter
 //' @param X1 a design matrix for nu
 //' @param X2 a design matrix for (transformed) kappa1
 //' @param X3 a design matrix for (transformed) kappa2
 //' @param X4 a design matrix for (transformed) delta
 //' @param ymat a matrix
 //' @param xmat a matrix
 //' @param wmat a matrix
 //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
 //' @return aggaussd0 a scalar, the negative log-likelihood
 //' @return aggaussd12 a matrix, first then second derivatives w.r.t. parameters
 //' @return aggaussd34 a matrix, third then fourth derivatives w.r.t. parameters
 //' @examples
 //' ## to follow
 //' @export
 // [[Rcpp::export]]
 double aggaussd0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, 
                    arma::mat X4, 
                 arma::mat ymat, arma::mat wmat, 
                 arma::uvec dupid, int dcate, arma::uvec nhere, 
                 double C, double epsilon = 0.1)
 {
   
   arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
   arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
   arma::vec p3vec = X3 * Rcpp::as<arma::vec>(pars[2]);
   arma::vec p4vec = X4 * Rcpp::as<arma::vec>(pars[3]);

   int nobs = ymat.n_rows;
   int ncol = ymat.n_cols;
   
   if (dcate == 1) {
     p1vec = p1vec.elem(dupid);
     p2vec = p2vec.elem(dupid);
     p3vec = p3vec.elem(dupid);
     p4vec = p4vec.elem(dupid);
   }
   
   double y, res, w, p1, p2, p3, p4;
   double nu, kappa1, kappa2, delta, kappa, S;
   double nllh = 0.0;

   for (int j=0; j < nobs; j++) {
     
     p1 = p1vec[j];
     p2 = p2vec[j];
     p3 = p3vec[j];
     p4 = p4vec[j];

     for (int l=0; l < ncol; l++) {
       
       y = ymat(j, l);

       if (std::isfinite(y)) {
         
         w = wmat(j, l);
         nu = p1;
         kappa1 = exp(p2);
         kappa2 = exp(p3);
         delta = exp(p4);
         
         nllh += w * (log(kappa1 + kappa2) + lgamma(1 / delta) - p4);
         
         res = y - nu;
         S = 1.0 / (1.0 + exp(-res / C));
         kappa = (1.0 - S) * kappa1 + S * kappa2;

         // if (y <= p1) {
         // 
         //   nllh += R_pow((nu - y) / kappa, delta);
         // 
         // } else {
         // 
         //   nllh += R_pow((y - nu) / kappa, delta);
         // 
         // }
         
         res = R_pow(res * res + epsilon * epsilon, 0.5 * delta);
         
         nllh += w * (res / R_pow(kappa, delta));

         // res = sqrt(res * res);
         // nllh += R_pow(res / kappa, delta);
         
       }
     }
   }
   
   return(nllh);
   
 }
 
 //' @rdname aggaussd0
 // [[Rcpp::export]]
 arma::mat aggaussd12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, 
                      arma::mat X4, 
                      arma::mat ymat, arma::mat wmat, 
                      arma::uvec dupid, int dcate, arma::uvec nhere, 
                      double C, double epsilon = 0.1)
 {
   
   arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
   arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
   arma::vec p3vec = X3 * Rcpp::as<arma::vec>(pars[2]);
   arma::vec p4vec = X4 * Rcpp::as<arma::vec>(pars[3]);
   
   int nobs = ymat.n_rows;
   int ncol = ymat.n_cols;
   
   arma::mat out = arma::mat(nobs, 14, arma::fill::zeros);

   if (dcate == 1) {
     p1vec = p1vec.elem(dupid);
     p2vec = p2vec.elem(dupid);
     p3vec = p3vec.elem(dupid);
     p4vec = p4vec.elem(dupid);
   }
   
   double y, w, p1, p2, p3, p4;

   for (int j=0; j < nobs; j++) {
     
     p1 = p1vec[j];
     p2 = p2vec[j];
     p3 = p3vec[j];
     p4 = p4vec[j];
     
     for (int l=0; l < ncol; l++) {
       
       y = ymat(j, l);

       if (std::isfinite(y)) {
         
         w = wmat(j, l);
         
         double ee1, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
         double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
         double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
         double ee31, ee32, ee33, ee34, ee36, ee37, ee38;
         double ee40, ee41, ee42, ee44, ee45, ee46, ee47;
         
         ee1 = y - p1;
         ee3 = exp(-(ee1/C));
         ee4 = 1 + ee3;
         ee5 = exp(p4);
         ee6 = exp(p2);
         ee7 = exp(p3);
         ee8 = 1 - 1/ee4;
         ee9 = ee8 * ee6;
         ee10 = ee9 + ee7/ee4;
         ee11 = R_pow(ee1, 2);
         ee12 = ee11 + R_pow(epsilon, 2);
         ee13 = ee5/2;
         ee14 = R_pow(ee12, ee13);
         ee15 = ee5 - 1;
         ee16 = 1 + ee5;
         ee17 = R_pow(ee10, ee16);
         ee18 = ee6 - ee7;
         ee19 = ee13 - 1;
         ee20 = R_pow(ee12, ee19);
         ee21 = ee6 + ee7;
         ee22 = log(ee10);
         ee23 = log(ee12);
         ee24 = R_pow(ee10, (2 + ee5));
         ee25 = R_pow(ee10, ee15);
         ee26 = R_pow(ee4, 2);
         ee27 = R_pow(ee10, (2 * ee5));
         ee28 = R_pow(ee10, (ee5 - 2));
         ee29 = R_pow(ee10, ee5);
         ee31 = ee23/2 - ee22;
         ee32 = 1 + ee5 * ee31;
         ee33 = ee24 * ee4;
         ee34 = ee20 * ee1;
         ee36 = ee14 * ee3 * ee18;
         ee37 = 1/ee5;
         ee38 = ee28 * ee8;
         ee40 = ee34/ee29 + ee36/(C * ee17 * ee26);
         ee41 = ee14 * ee8;
         ee42 = ee14 * ee5;
         ee44 = ee32/ee17 * ee14;
         ee45 = 1/ee21;
         ee46 = C * ee26;
         ee47 = R::digamma(ee37);
         
         out(j, 0) += w * (-(ee40 * ee5));
         out(j, 1) += w * ((ee45 - ee41 * ee5/ee17) * ee6);
         out(j, 2) += w * ((ee45 - ee42/(ee17 * ee4)) * ee7);
         out(j, 3) += w * (ee42 * ee31/ee29 - (1 + ee47/ee5));
         out(j, 4) += w * (-(((((ee25 * (1 - 2 * (ee3/ee4)) + ee28 *
           ee3 * ee18 * ee15/ee26) * ee14/C - ee25 * ee20 * ee5 * ee1)/
             ee27 - (ee34/ee17 + 2 * (ee36/(C * ee24 * ee26))) * ee5) * ee3 *
               ee18/ee46 - (ee20 + 2 * (R_pow(ee12, (ee13 - 2)) * ee19 *
               ee11))/ee29) * ee5));
         out(j, 5) += w * (-((((ee25 + ee38 * ee18 * ee15)/ee27 - 2 *
           (ee8 * ee18 * ee5/ee24)) * ee14 * ee3/ee46 - ee20 * ee8 * ee5 *
           ee1/ee17) * ee6 * ee5));
         out(j, 6) += w * (-((((ee28 * ee18 * ee15/ee4 - ee25)/ee27 -
           2 * (ee18 * ee5/ee33)) * ee14 * ee3/(C * ee4) - ee20 * ee5 *
           ee1/ee17) * ee7 * ee5/ee4));
         out(j, 7) += w * (-(ee40 * ee32 * ee5));
         out(j, 8) += w * (((1 - ee6/ee21)/ee21 - ((ee25 + ee38 * ee6 *
           ee15)/ee27 - 2 * (ee9 * ee5/ee24)) * ee14 * ee8 * ee5) *
           ee6);
         out(j, 9) += w * (-((1/R_pow(ee21, 2) - ee41 * ee16 * ee5/ee33) *
           ee6 * ee7));
         out(j, 10) += w * (-(ee44 * ee8 * ee6 * ee5));
         out(j, 11) += w * (((1 - ee7/ee21)/ee21 - ((ee25 + ee28 * ee7 *
           ee15/ee4)/ee27 - 2 * (ee7 * ee5/ee33)) * ee14 * ee5/ee4) *
           ee7);
         out(j, 12) += w * (-(ee44 * ee7 * ee5/ee4));
         out(j, 13) += w * (((0.5 + ee5 * (ee23/4 - ee22/2)) * ee23 -
           ee32 * ee22) * ee14 * ee5/ee29 + (ee47 + R::trigamma(ee37)/ee5)/
             ee5);
         }}}

   return out;

 }

 //' @rdname aggaussd0
 // [[Rcpp::export]]
 arma::mat aggaussd34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, 
                       arma::mat X4, 
                      arma::mat ymat, arma::mat wmat, 
                      arma::uvec dupid, int dcate, arma::uvec nhere, 
                      double C, double epsilon = 0.1)
 {
   
   arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
   arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
   arma::vec p3vec = X3 * Rcpp::as<arma::vec>(pars[2]);
   arma::vec p4vec = X4 * Rcpp::as<arma::vec>(pars[3]);
   
   int nobs = ymat.n_rows;
   int ncol = ymat.n_cols;
   
   arma::mat out = arma::mat(nobs, 20, arma::fill::zeros);
   
   if (dcate == 1) {
     p1vec = p1vec.elem(dupid);
     p2vec = p2vec.elem(dupid);
     p3vec = p3vec.elem(dupid);
     p4vec = p4vec.elem(dupid);
   }
   
   double y, w, p1, p2, p3, p4;
   
   for (int j=0; j < nobs; j++) {
     
     p1 = p1vec[j];
     p2 = p2vec[j];
     p3 = p3vec[j];
     p4 = p4vec[j];
     
     for (int l=0; l < ncol; l++) {
       
       y = ymat(j, l);
       w = wmat(j, l);
       
       if (std::isfinite(y)) {
         
         double ee1, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
         double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
         double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
         double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee37, ee38, ee39;
         double ee41, ee42, ee43, ee44, ee45, ee46, ee47;
         double ee51, ee52, ee53, ee54, ee55, ee56, ee57, ee58, ee59;
         double ee60, ee61, ee62, ee63, ee64, ee65, ee68, ee69;
         double ee70, ee71, ee72, ee73, ee75, ee76, ee77, ee79;
         double ee80, ee81, ee82, ee83, ee84, ee85, ee86, ee87, ee89;
         double ee90, ee91, ee92, ee93, ee95, ee96, ee99;
         double ee100, ee103, ee104, ee105, ee106, ee107, ee108, ee109;
         double ee110, ee111, ee113, ee115, ee116, ee117, ee118, ee119;
         double ee120, ee121, ee122, ee123, ee124, ee125;
         
         ee1 = y - p1;
         ee3 = exp(-(ee1/C));
         ee4 = 1 + ee3;
         ee5 = exp(p4);
         ee6 = exp(p3);
         ee7 = exp(p2);
         ee8 = 1 - 1/ee4;
         ee9 = ee8 * ee7;
         ee10 = ee6/ee4;
         ee11 = ee9 + ee10;
         ee12 = ee5 - 1;
         ee13 = R_pow(ee1, 2);
         ee14 = ee13 + R_pow(epsilon, 2);
         ee15 = ee5/2;
         ee16 = ee5 - 2;
         ee17 = R_pow(ee11, ee12);
         ee18 = log(ee11);
         ee19 = ee7 - ee6;
         ee20 = R_pow(ee11, ee16);
         ee21 = 2 * ee5;
         ee22 = ee15 - 1;
         ee23 = log(ee14);
         ee24 = R_pow(ee14, ee15);
         ee25 = R_pow(ee14, ee22);
         ee26 = R_pow(ee4, 2);
         ee27 = ee5 * ee18;
         ee28 = ee23/2;
         ee29 = R_pow(ee11, ee21);
         ee30 = ee20 * ee8;
         ee31 = R_pow(ee11, (2 + ee5));
         ee32 = 1 - 2 * (ee3/ee4);
         ee33 = R_pow(ee11, ee5);
         ee34 = R_pow(ee11, (3 + ee5));
         ee35 = ee17 * ee32;
         ee36 = ee5 * (ee18 + ee28);
         ee37 = R_pow(ee11, (2 * ee12));
         ee38 = R_pow(ee11, (4 * ee5));
         ee39 = ee35 + ee20 * ee3 * ee19 * ee12/ee26;
         ee41 = 1 + ee5;
         ee42 = 1 + ee36;
         ee43 = 2 * (1 + 2 * ee27);
         ee44 = ee7 + ee6;
         ee45 = R_pow(ee11, (1 + ee21));
         ee46 = R_pow(ee11, ee41);
         ee47 = ee17 + ee30 * ee19 * ee12;
         ee51 = ee20 * ee19 * ee12/ee4 - ee17;
         ee52 = ee15 - 2;
         ee53 = R_pow(ee11, (ee21 - 1));
         ee54 = ee20 * ee6;
         ee55 = R_pow(ee14, ee52);
         ee56 = C * ee26;
         ee57 = ee17 + ee30 * ee7 * ee12;
         ee58 = ee17 + ee54 * ee12/ee4;
         ee59 = ee34 * ee4;
         ee60 = ee17 * ee25;
         ee61 = ee25 + 2 * (ee55 * ee22 * ee13);
         ee62 = 1 + ee27;
         ee63 = ee43 + 2 * ee42;
         ee64 = 8 * ee27;
         ee65 = C * ee4;
         ee68 = ee39 * ee24/C - ee60 * ee5 * ee1;
         ee69 = ((ee12 * ee18 + 2) * ee5 - 1) * ee20;
         ee70 = R_pow(ee11, (ee5 - 3));
         ee71 = ee63 - ee64;
         ee72 = ee31 * ee4;
         ee73 = ee17 * ee62;
         ee75 = ((1.5 + ee5 * (ee18 + ee23/4)) * ee23 + (3 + 9 *  ee27 - (ee43 + 4 * ee42)) * ee18) * ee5 + 1;
         ee76 = ee8 * ee19;
         ee77 = ee19 * ee5;
         ee79 = ee6 * ee5;
         ee80 = ee28 - 2 * ee18;
         ee81 = ee28 - ee18;
         ee82 = (ee58 * ee33 + ee37 * ee6 * ee5/ee4) * ee17;
         ee83 = ee39 * ee25;
         ee84 = ee51 * ee53;
         ee85 = ee11 * ee4;
         ee86 = ee37 * ee8;
         ee87 = ee76 * ee5;
         ee89 = 1/ee5;
         ee90 = 2 + ee5 * ee81;
         ee91 = 8 * (ee77/ee59);
         ee92 = ee79/ee59;
         ee93 = ee5 * ee23;
         ee95 = (ee57 * ee33 + ee86 * ee7 * ee5) * ee17;
         ee96 = (ee37 * ee5 + R_pow(ee11, (ee21 - 2)) * ee12) * ee17;
         ee99 = ee83 * ee1/C - ee61 * ee17;
         ee100 = ee47 * ee53;
         ee103 = ee75/ee46 * ee24;
         ee104 = ee69 * ee8;
         ee105 = ee61/ee46;
         ee106 = R_pow(ee11, (3 * ee5 - 3));
         ee107 = ee20 + ee70 * ee8 * ee7 * ee16;
         ee108 = ee20 + ee70 * ee6 * ee16/ee4;
         ee109 = ee25 * ee8;
         ee110 = ee25 * ee90;
         ee111 = ee8 * ee71;
         ee113 = ee9 * ee5;
         ee115 = ee19 * ee16/ee85;
         ee116 = R_pow(ee44, 2);
         ee117 = 2 * (ee68/ee45);
         ee118 = 2 * ee82;
         ee119 = 2 * (ee57/ee45);
         ee120 = 2 * ee84;
         ee121 = 2 * (ee7/ee44);
         ee122 = 2 * (ee6/ee44);
         ee123 = 8 * (ee87/ee34);
         ee124 = 8 * ee92;
         ee125 = ee93/4;
         
         out(j, 0) += w * (-((((((((ee115 - 6) * ee3/ee4 + 3) * ee3 * ee19 *
           ee12/ee26 + ee11 * (1 - ((2 * (1 + 2 * ee3) + 2 * ee4 - 8 *
           ee3)/ee4 + 2) * ee3/ee4)) * ee20 * ee24/C - ee83 * ee5 *
           ee1)/C - ee99 * ee5)/ee29 - (ee99/ee29 + ((2 * ((ee39 * ee33 +
           ee37 * ee3 * ee19 * ee5/ee26) * ee17 * ee24/C) + 2 * (ee68 *
           ee53))/ee38 + ee117 - (2 * (ee25 * ee1/ee31) + 8 * (ee24 *
           ee3 * ee19/(C * ee34 * ee26))) * ee5) * ee3 * ee19/ee56 -
           ee105) * ee5) * ee3 * ee19/ee56 + ee55 * (4 * (ee52 * ee13/
             ee14) + 6) * ee22 * ee1/ee33) * ee5));
         out(j, 1) += w * (-(((((ee47 * ee32 + (ee76 * ee16/ee11 + 2) * ee20 *
           ee3 * ee19 * ee12/ee26) * ee24/C - ee47 * ee25 * ee5 *
           ee1)/ee29 - (((2 * ((ee47 * ee33 + ee86 * ee19 * ee5) * ee17) +
           2 * ee100)/ee38 - ee123) * ee24 * ee3/ee56 + 2 * (ee68 *
           ee8/ee45)) * ee19 * ee5) * ee3/ee56 - ((ee47/ee29 - 2 * (ee87/
             ee31)) * ee25 * ee3 * ee1/ee56 - ee61 * ee8/ee46) * ee5) *
               ee7 * ee5));
         out(j, 2) += w * (-(((((ee51 * ee32 + ee20 * (ee115 - 2) * ee3 *
           ee19 * ee12/ee26) * ee24/C - ee51 * ee25 * ee5 * ee1)/ee29 -
           (((2 * ((ee51 * ee33 + ee37 * ee19 * ee5/ee4) * ee17) + ee120)/
             ee38 - ee91) * ee24 * ee3/ee65 + ee117) * ee19 * ee5/ee4) *
               ee3/ee65 - ((ee51/ee29 - 2 * (ee77/ee72)) * ee25 * ee3 *
               ee1/ee65 - ee105) * ee5) * ee6 * ee5/ee4));
         out(j, 3) += w * (-(((((ee39 * ee5 * ee23/2 + ee69 * ee3 * ee19/
           ee26 + ee35 * ee62) * ee24/C - (ee60 * (2 + ee36) * ee1 + 2 *
             (ee68 * ee18)) * ee5)/ee29 - (ee110 * ee1/ee46 + ee24 * ee71 *
             ee3 * ee19/(C * ee31 * ee26)) * ee5) * ee3 * ee19/ee56 -
             (((ee22 * ee23 + 1) * ee5 + 2 * ee22) * ee55 * ee13 + ee25 *
             (1 + ee93/2) - ee61 * ee5 * ee18)/ee33) * ee5));
         out(j, 4) += w * (-(((((ee107 * ee19 + 2 * (ee20 * ee7)) * ee8 *
           ee12 + ee17)/ee29 - ((2 * (ee95 * ee19) + 2 * (ee100 * ee7))/
             ee38 + (2 * (ee47/ee45) - ee123) * ee7) * ee8 * ee5) * ee24 *
               ee3/ee56 - (ee57/ee29 - 2 * (ee113/ee31)) * ee25 * ee8 *
               ee5 * ee1) * ee7 * ee5));
         out(j, 5) += w * (-(((((ee70 * ee19 * ee16/ee4 - ee20) * ee8 + ee20/
           ee4) * ee12/ee29 - (((2 * (ee96 * ee19/ee4) + ee120)/ee38 -
             ee91) * ee8 + 2 * (ee47/(ee45 * ee4))) * ee5) * ee24 * ee3/
               ee65 + ee109 * ee41 * ee5 * ee1/ee31) * ee7 * ee6 * ee5/ee4));
         out(j, 6) += w * (-((((ee47 * ee5 * ee80 + ee104 * ee19 + ee73)/
           ee29 - ee111 * ee19 * ee5/ee31) * ee24 * ee3/ee56 - ee109 * ee90 *
             ee5 * ee1/ee46) * ee7 * ee5));
         out(j, 7) += w * (-(((((ee108 * ee19 - 2 * ee54) * ee12/ee4 - ee17)/
           ee29 - ((2 * (ee82 * ee19) + 2 * (ee84 * ee6))/ee38 + (2 *
             (ee51/ee45) - ee91) * ee6) * ee5/ee4) * ee24 * ee3/ee65 -
             (ee58/ee29 - 2 * (ee79/ee72)) * ee25 * ee5 * ee1) * ee6 * ee5/
               ee4));
         out(j, 8) += w * (-((((ee51 * ee5 * ee80 + ee69 * ee19/ee4 - ee73)/
           ee29 - ee71 * ee19 * ee5/ee72) * ee24 * ee3/ee65 - ee110 *
             ee5 * ee1/ee46) * ee6 * ee5/ee4));
         out(j, 9) += w * (-((ee75 * ee24 * ee3 * ee19/(C * ee46 * ee26) +
           (((1.5 + ee125) * ee23 - (3 + ee5 * (ee23 - ee18)) * ee18) *
           ee5 + 1) * ee25 * ee1/ee33) * ee5));
         out(j, 10) += w * (((1 - (3 - ee121) * ee7/ee44)/ee44 - ((((ee9 *
           ee16/ee11 + 3) * ee12 + 1) * ee8 * ee7 + ee10)/ee31 - ((2 *
           ee95 + 2 * (ee57 * ee53))/ee38 + ee119 - 8 * (ee113/ee34)) *
           ee8 * ee7 * ee5) * ee24 * ee8 * ee5) * ee7);
         out(j, 11) += w * (-(((ee107 * ee12/ee29 - (((2 * ee96 + 2 * (ee106 *
           ee12))/ee38 - 8 * (ee5/ee34)) * ee8 * ee7 + ee119) * ee5) *
           ee24 * ee8 * ee5/ee4 + (1 - ee121)/ee116) * ee7 * ee6));
         out(j, 12) += w * (-(((ee57 * ee5 * ee80 + ee104 * ee7 + ee73)/ee29 -
           ee111 * ee7 * ee5/ee31) * ee24 * ee8 * ee7 * ee5));
         out(j, 13) += w * (-((((ee108/ee29 - 2 * ee92) * ee12 - ((ee118 +
           2 * (ee106 * ee6 * ee12/ee4))/ee38 - ee124) * ee5) * ee24 *
           ee8 * ee5/ee4 + (1 - ee122)/ee116) * ee7 * ee6));
         out(j, 14) += w * (-(((ee12 * ee81 + 2 + ee64 - ee63) * ee5 - 1)/
           ee31 * ee24 * ee8 * ee7 * ee6 * ee5/ee4));
         out(j, 15) += w * (-(ee103 * ee8 * ee7 * ee5));
         out(j, 16) += w * (((1 - (3 - ee122) * ee6/ee44)/ee44 - ((((2 + ee6 *
           ee16/ee85) * ee12 + ee5) * ee6/ee4 + ee9)/ee31 - ((ee118 +
           2 * (ee58 * ee53))/ee38 + 2 * (ee58/ee45) - ee124) * ee6 *
           ee5/ee4) * ee24 * ee5/ee4) * ee6);
         out(j, 17) += w * (-(((ee58 * ee5 * ee80 + ee69 * ee6/ee4 + ee73)/
           ee29 - ee71 * ee6 * ee5/ee72) * ee24 * ee6 * ee5/ee4));
         out(j, 18) += w * (-(ee103 * ee6 * ee5/ee4));
         out(j, 19) += w * ((((((0.5 + ee125)/2 + 0.5) * ee23 - ((0.5 * ee23 -
           ee18) * ee5 + ee62/2 + 1) * ee18) * ee5 + 0.5) * ee23 -
           ee75 * ee18) * ee24 * ee5/ee33 - ((3 * R::trigamma(ee89) + R::psigamma(ee89, 2)/
             ee5)/ee5 + R::digamma(ee89))/ee5);
         
         }}}    
   
   return out;
   
 }
