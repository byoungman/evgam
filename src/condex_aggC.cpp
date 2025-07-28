// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

double ltind(double x, double u) {
  double out = 0.0;
  if (x < u) {
    out = 1.0;
  }
  return out;
}

//' Conditional extreme value model negative log-likelihood for asymmetric generalised Gaussian
 //'
 //' @param pars a list of vectors of coefficients for each conditional EVD parameter
 //' @param X1 a design matrix for (transformed) alpha
 //' @param X2 a design matrix for (transformed) beta
 //' @param X3 a design matrix for nu
 //' @param X4 a design matrix for (transformed) kappa1
 //' @param X5 a design matrix for (transformed) kappa2
 //' @param X6 a design matrix for (transformed) delta
 //' @param ymat a matrix
 //' @param xmat a matrix
 //' @param wmat a matrix
 //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
 //' @return condexaggd0 a scalar, the negative log-likelihood
 //' @return condexaggd12 a matrix, first then second derivatives w.r.t. parameters
 //' @return condexaggd34 a matrix, third then fourth derivatives w.r.t. parameters
 //' @examples
 //' ## to follow
 //' @export
 // [[Rcpp::export]]
 double condexaggd0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, 
                    arma::mat X4, arma::mat X5, arma::mat X6, 
                 arma::mat ymat, arma::mat xmat, arma::mat wmat, 
                 arma::uvec dupid, int dcate, arma::uvec nhere,
                 double C, double epsilon = 0.1)
 {
   
   arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
   arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
   arma::vec p3vec = X3 * Rcpp::as<arma::vec>(pars[2]);
   arma::vec p4vec = X4 * Rcpp::as<arma::vec>(pars[3]);
   arma::vec p5vec = X5 * Rcpp::as<arma::vec>(pars[4]);
   arma::vec p6vec = X6 * Rcpp::as<arma::vec>(pars[5]);

   int nobs = ymat.n_rows;
   int ncol = ymat.n_cols;
   
   if (dcate == 1) {
     p1vec = p1vec.elem(dupid);
     p2vec = p2vec.elem(dupid);
     p3vec = p3vec.elem(dupid);
     p4vec = p4vec.elem(dupid);
     p5vec = p5vec.elem(dupid);
     p6vec = p6vec.elem(dupid);
   }
   
   double y, w, yi, z, res, p1, p2, p3, p4, p5, p6;
   double alpha, beta, nu, kappa, kappa1, kappa2, delta, S;
   double nllh = 0.0;
   
   for (int j=0; j < nobs; j++) {
     
     p1 = p1vec[j];
     p2 = p2vec[j];
     p3 = p3vec[j];
     p4 = p4vec[j];
     p5 = p5vec[j];
     p6 = p6vec[j];

     for (int l=0; l < ncol; l++) {
       
       y = ymat(j, l);
       yi = xmat(j, l);
       
       if (arma::is_finite(y) & arma::is_finite(yi)) {
         
         w = wmat(j, l);
         
         alpha = 2.0 / (1.0 + exp(-p1)) - 1.0;
         beta = 1.0 / (1.0 + exp(-p2));
         z = (y - alpha * yi) / R_pow(yi, beta);
         
         w = wmat(j, l);
         nu = p3;
         kappa1 = exp(p4);
         kappa2 = exp(p5);
         delta = exp(p6);
         
         nllh += w * (log(kappa1 + kappa2) + lgamma(1 / delta) - p6);
         
         res = z - nu;
         S = 1.0 / (1.0 + exp(-res / C));
         kappa = (1.0 - S) * kappa1 + S * kappa2;
         
         res = R_pow(res * res + epsilon * epsilon, 0.5 * delta);
         
         nllh += w * (res / R_pow(kappa, delta));
         
       }
     }
   }
   
   return(nllh);
   
 }
 
 //' @rdname condexaggd0
 // [[Rcpp::export]]
 arma::mat condexaggd12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, 
                        arma::mat X4, arma::mat X5, arma::mat X6, 
                        arma::mat ymat, arma::mat xmat, arma::mat wmat, 
                        arma::uvec dupid, int dcate, arma::uvec nhere,
                        double C, double epsilon = 0.1)
 {
   
   arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
   arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
   arma::vec p3vec = X3 * Rcpp::as<arma::vec>(pars[2]);
   arma::vec p4vec = X4 * Rcpp::as<arma::vec>(pars[3]);
   arma::vec p5vec = X5 * Rcpp::as<arma::vec>(pars[4]);
   arma::vec p6vec = X6 * Rcpp::as<arma::vec>(pars[5]);
   
   int nobs = ymat.n_rows;
   int ncol = ymat.n_cols;
   
   if (dcate == 1) {
     p1vec = p1vec.elem(dupid);
     p2vec = p2vec.elem(dupid);
     p3vec = p3vec.elem(dupid);
     p4vec = p4vec.elem(dupid);
     p5vec = p5vec.elem(dupid);
     p6vec = p6vec.elem(dupid);
   }
   
   double y, w, yi, z, res, p1, p2, p3, p4, p5, p6;
   double alpha, beta, nu, kappa, kappa1, kappa2, delta, S;
   
   arma::mat out = arma::mat(nobs, 27, arma::fill::zeros);
   
   for (int j=0; j < nobs; j++) {
     
     p1 = p1vec[j];
     p2 = p2vec[j];
     p3 = p3vec[j];
     p4 = p4vec[j];
     p5 = p5vec[j];
     p6 = p6vec[j];
     
     for (int l=0; l < ncol; l++) {
       
       y = ymat(j, l);
       yi = xmat(j, l);
       
       if (arma::is_finite(y) & arma::is_finite(yi)) {
         
         w = wmat(j, l);
         
         double ee2, ee4, ee5, ee6, ee7;
         double ee10, ee11, ee13, ee15, ee16, ee17, ee18, ee19;
         double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
         double ee30, ee31, ee32, ee33, ee34, ee35, ee37, ee38, ee39;
         double ee40, ee41, ee42, ee43, ee44, ee45, ee46, ee47, ee48;
         double ee50, ee51, ee53, ee54, ee56, ee57, ee58, ee59;
         double ee60, ee62, ee64, ee65, ee66, ee67, ee68, ee69;
         double ee71, ee72, ee73, ee74, ee75, ee78;
         double ee80, ee81, ee82, ee84, ee85, ee86, ee88;
         double ee90, ee91, ee92, ee93, ee94, ee98;
         double ee102, ee106, ee107, ee108, ee109;
         double ee110, ee111, ee113, ee116, ee118, ee119;
         double ee120, ee121, ee122, ee123, ee124, ee125;
         
         ee2 = exp(-p2);
         ee4 = exp(-p1);
         ee5 = 1 + ee2;
         ee6 = 1 + ee4;
         ee7 = 1/ee5;
         ee10 = y - yi * (2/ee6 - 1);
         ee11 = R_pow(yi, ee7);
         ee13 = ee10/ee11 - p3;
         ee15 = exp(-(ee13/C));
         ee16 = 1 + ee15;
         ee17 = exp(p6);
         ee18 = exp(p4);
         ee19 = exp(p5);
         ee20 = 1 - 1/ee16;
         ee21 = ee20 * ee18;
         ee22 = ee21 + ee19/ee16;
         ee23 = R_pow(ee13, 2);
         ee24 = ee23 + R_pow(epsilon, 2);
         ee25 = ee17/2;
         ee26 = ee25 - 1;
         ee27 = ee17 - 1;
         ee28 = R_pow(ee24, ee25);
         ee29 = R_pow(ee24, ee26);
         ee30 = R_pow(ee16, 2);
         ee31 = ee18 - ee19;
         ee32 = 1 + ee17;
         ee33 = R_pow(ee22, ee32);
         ee34 = R_pow(ee22, ee27);
         ee35 = ee29 * ee13;
         ee37 = 2 * ee18 - 2 * ee19;
         ee38 = R_pow(ee22, (2 + ee17));
         ee39 = log(yi);
         ee40 = R_pow(ee22, (2 * ee17));
         ee41 = R_pow(ee22, (ee17 - 2));
         ee42 = R_pow(ee22, ee17);
         ee43 = R_pow(ee6, 2);
         ee44 = R_pow(ee5, 2);
         ee45 = C * ee30;
         ee46 = 1 - ee7;
         ee47 = log(ee24);
         ee48 = log(ee22);
         ee50 = ee28 * ee15 * ee31;
         ee51 = ee18 + ee19;
         ee53 = ee47/2 - ee48;
         ee54 = R_pow(yi, ee46);
         ee56 = ee29 * ee34 * ee13;
         ee57 = ee35/ee33;
         ee58 = R_pow(ee24, (ee25 - 2));
         ee59 = 1 + ee17 * ee53;
         ee60 = 2/ee5;
         ee62 = C * ee33 * ee30;
         ee64 = C * ee38 * ee30;
         ee65 = ee15/ee16;
         ee66 = ee35/ee42;
         ee67 = ee38 * ee16;
         ee68 = ee66 + ee50/ee62;
         ee69 = ee58 * ee23;
         ee71 = ee28 * ee37 * ee15;
         ee72 = ee41 * ee20;
         ee73 = 2 * ee46;
         ee74 = ee4 * ee17;
         ee75 = ee56 * ee37;
         ee78 = ee35 * ee20 * ee17/ee33;
         ee80 = ee35 * ee17/ee33;
         ee81 = ee57 + 2 * (ee50/ee64);
         ee82 = ee69 * ee26;
         ee84 = ee41 * ee37;
         ee85 = ee41 * ee15;
         ee86 = 1 - 2 * ee65;
         ee88 = 1/ee17;
         ee90 = 2 * ee57 + 2 * (ee71/ee64);
         ee91 = C * ee16;
         ee92 = R_pow(yi, (1 - ee60));
         ee93 = R_pow(yi, ee73);
         ee94 = R_pow(yi, ee60);
         ee98 = (((ee34 * ee86 + ee85 * ee31 * ee27/ee30) * ee28/C -  ee56 * ee17)/ee40 - ee81 * ee17) * ee15 * ee31/ee45 -  (ee29 + 2 * ee82)/ee42;
         ee102 = ((ee34 + ee72 * ee31 * ee27)/ee40 - 2 * (ee20 *  ee31 * ee17/ee38)) * ee28 * ee15/ee45 - ee78;
         ee106 = ((ee41 * ee31 * ee27/ee16 - ee34)/ee40 - 2 * (ee31 *  ee17/ee67)) * ee28 * ee15/ee91 - ee80;
         ee107 = ee68/ee11;
         ee108 = ee75 * ee17;
         ee109 = ee28 * (ee59/ee33);
         ee110 = ee28 * ee20;
         ee111 = ee28 * ee17;
         ee113 = ee84 * ee15 * ee27;
         ee116 = (2 - ee39/ee5) * ee2/ee5 - 1;
         ee118 = ee90 * ee31 * ee17;
         ee119 = 1/ee51;
         ee120 = 2 - 4 * ee65;
         ee121 = 2 - 4 * (ee4/ee6);
         ee122 = 2 * ee34;
         ee123 = R::digamma(ee88);
         ee124 = ee54 * ee13;
         ee125 = ee94 * ee44;
         
         out(j, 0) += w * (-(ee74 * (ee54 * (ee71/ee62 + 2 * ee66))/
           ee43));
         out(j, 1) += w * (-(ee107 * ee2 * ee17 * ee39 * ee10/ee44));
         out(j, 2) += w * (-(ee68 * ee17));
         out(j, 3) += w * ((ee119 - ee110 * ee17/ee33) * ee18);
         out(j, 4) += w * ((ee119 - ee111/(ee33 * ee16)) * ee19);
         out(j, 5) += w * (ee111 * ee53/ee42 - (1 + ee123/ee17));
         out(j, 6) += w * (-(ee74 * (ee93 * (((ee34 * ((4 - 8 * ee65) *
           ee4/(C * ee43) - R_pow(yi, (1 - (ee7 + ee73))) * ee121) *
           ee31 + ee41 * R_pow(ee37, 2) * ee15 * ee4 * ee27/(ee45 * ee43)) *
           ee28 - 2 * (ee75 * ee4 * ee17/ee43))/ee40 - ee90 * ee37 *
           ee4 * ee17/ee43) * ee15/ee45 - (ee29 * (4 * (ee93 * ee4/
             ee43) + ee124 * ee121) + 8 * (ee93 * ee58 * ee23 * ee4 * ee26/
               ee43))/ee42)/ee43));
         out(j, 7) += w * (-(ee4 * ee2 * ee17 * ee39 * (ee92 * (((ee34 *
           (ee120 * ee10/C - 2 * ee11) + ee113 * ee10/ee45) * ee28 *
           ee31 - ee108 * ee10)/ee40 - ee118 * ee10) * ee15/ee45 - (ee29 *
           (2 * ee124 + 2 * (ee92 * ee10)) + 4 * (ee92 * ee58 * ee23 *
           ee26 * ee10))/ee42)/(ee43 * ee44)));
         out(j, 8) += w * (-(ee74 * (ee54 * ((((ee34 * ee120 + ee113/
           ee30) * ee28 * ee31/C - ee108)/ee40 - ee118) * ee15/ee45 - (2 *
             ee29 + 4 * ee82)/ee42))/ee43));
         out(j, 9) += w * (-(ee4 * ee18 * ee17 * (ee54 * (((ee72 * ee37 *
           ee27 + ee122)/ee40 - 2 * (ee20 * ee37 * ee17/ee38)) * ee28 *
           ee15/ee45 - 2 * ee78))/ee43));
         out(j, 10) += w * (-(ee4 * ee19 * ee17 * (ee54 * (((ee84 * ee27/
           ee16 - ee122)/ee40 - 2 * (ee37 * ee17/ee67)) * ee28 * ee15/
             ee91 - 2 * ee80))/(ee16 * ee43)));
         out(j, 11) += w * (-(ee74 * (ee54 * (ee35 * (2 + ee17 * (ee47 -
           2 * ee48))/ee42 + ee28 * ee59 * ee37 * ee15/ee62))/ee43));
         out(j, 12) += w * (-(((((ee34 * (ee86 * ee2 * ee39 * ee10/(C *
           ee44) + ee11 * ee116) + ee85 * ee2 * ee31 * ee27 * ee39 *
           ee10/(ee45 * ee44)) * ee28 - ee56 * ee2 * ee17 * ee39 * ee10/
             ee44)/ee40 - ee81 * ee2 * ee17 * ee39 * ee10/ee44)/ee94 *
               ee15 * ee31/ee45 + ((ee116 * ee13/ee11 - ee2 * ee39 * ee10/
                 ee125) * ee29 - 2 * (ee69 * ee2 * ee26 * ee39 * ee10/ee125))/
                   ee42) * ee2 * ee17 * ee39 * ee10/ee44));
         out(j, 13) += w * (-(ee98/ee11 * ee2 * ee17 * ee39 * ee10/ee44));
         out(j, 14) += w * (-(ee102/ee11 * ee2 * ee18 * ee17 * ee39 *
           ee10/ee44));
         out(j, 15) += w * (-(ee106/ee11 * ee2 * ee19 * ee17 * ee39 *
           ee10/(ee16 * ee44)));
         out(j, 16) += w * (-(ee107 * ee59 * ee2 * ee17 * ee39 * ee10/
           ee44));
         out(j, 17) += w * (-(ee98 * ee17));
         out(j, 18) += w * (-(ee102 * ee18 * ee17));
         out(j, 19) += w * (-(ee106 * ee19 * ee17/ee16));
         out(j, 20) += w * (-(ee68 * ee59 * ee17));
         out(j, 21) += w * (((1 - ee18/ee51)/ee51 - ((ee34 + ee72 * ee18 *
           ee27)/ee40 - 2 * (ee21 * ee17/ee38)) * ee28 * ee20 * ee17) *
           ee18);
         out(j, 22) += w * (-((1/R_pow(ee51, 2) - ee110 * ee32 * ee17/
           ee67) * ee18 * ee19));
         out(j, 23) += w * (-(ee109 * ee20 * ee18 * ee17));
         out(j, 24) += w * (((1 - ee19/ee51)/ee51 - ((ee34 + ee41 * ee19 *
           ee27/ee16)/ee40 - 2 * (ee19 * ee17/ee67)) * ee28 * ee17/
             ee16) * ee19);
         out(j, 25) += w * (-(ee109 * ee19 * ee17/ee16));
         out(j, 26) += w * (ee28 * ((0.5 + ee17 * (ee47/4 - ee48/2)) *
           ee47 - ee59 * ee48) * ee17/ee42 + (ee123 + R::trigamma(ee88)/
             ee17)/ee17);
         }}}    
   
   return out;
   
 }
 
 //' @rdname condexaggd0
 // [[Rcpp::export]]
 arma::mat condexaggd34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, 
                        arma::mat X4, arma::mat X5, arma::mat X6, 
                        arma::mat ymat, arma::mat xmat, arma::mat wmat, 
                        arma::uvec dupid, int dcate, arma::uvec nhere,
                        double C, double epsilon = 0.1)
 {
   
   arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
   arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
   arma::vec p3vec = X3 * Rcpp::as<arma::vec>(pars[2]);
   arma::vec p4vec = X4 * Rcpp::as<arma::vec>(pars[3]);
   arma::vec p5vec = X5 * Rcpp::as<arma::vec>(pars[4]);
   arma::vec p6vec = X6 * Rcpp::as<arma::vec>(pars[5]);
   
   int nobs = ymat.n_rows;
   int ncol = ymat.n_cols;
   
   if (dcate == 1) {
     p1vec = p1vec.elem(dupid);
     p2vec = p2vec.elem(dupid);
     p3vec = p3vec.elem(dupid);
     p4vec = p4vec.elem(dupid);
     p5vec = p5vec.elem(dupid);
     p6vec = p6vec.elem(dupid);
   }
   
   double y, w, yi, z, res, p1, p2, p3, p4, p5, p6;
   double alpha, beta, nu, kappa, kappa1, kappa2, delta, S;
   
   arma::mat out = arma::mat(nobs, 56, arma::fill::zeros);
   
   for (int j=0; j < nobs; j++) {
     
     p1 = p1vec[j];
     p2 = p2vec[j];
     p3 = p3vec[j];
     p4 = p4vec[j];
     p5 = p5vec[j];
     p6 = p6vec[j];
     
     for (int l=0; l < ncol; l++) {
       
       y = ymat(j, l);
       yi = xmat(j, l);
       
       if (arma::is_finite(y) & arma::is_finite(yi)) {
         
         w = wmat(j, l);
         
       //   double ee2, ee4, ee5, ee6, ee7;
       //   double ee10, ee11, ee13, ee15, ee16, ee17, ee18, ee19;
       //   double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
       //   double ee30, ee31, ee32, ee33, ee34, ee35, ee37, ee38, ee39;
       //   double ee40, ee41, ee42, ee43, ee44, ee45, ee46, ee47, ee48;
       //   double ee50, ee51, ee53, ee54, ee56, ee57, ee58, ee59;
       //   double ee60, ee62, ee64, ee65, ee66, ee67, ee68, ee69;
       //   double ee71, ee72, ee73, ee74, ee75, ee78;
       //   double ee80, ee81, ee82, ee84, ee85, ee86, ee88;
       //   double ee90, ee91, ee92, ee93, ee94, ee98;
       //   double ee102, ee106, ee107, ee108, ee109;
       //   double ee110, ee111, ee113, ee116, ee118, ee119;
       //   double ee120, ee121, ee122, ee123, ee124, ee125;
       //   
       //   ee2 = exp(-p2);
       //   ee4 = exp(-p1);
       //   ee5 = 1 + ee2;
       //   ee6 = 1 + ee4;
       //   ee7 = 1/ee5;
       //   ee10 = y - yi * (2/ee6 - 1);
       //   ee11 = R_pow(yi, ee7);
       //   ee13 = ee10/ee11 - p3;
       //   ee15 = exp(-(ee13/C));
       //   ee16 = 1 + ee15;
       //   ee17 = exp(p6);
       //   ee18 = exp(p4);
       //   ee19 = exp(p5);
       //   ee20 = 1 - 1/ee16;
       //   ee21 = ee20 * ee18;
       //   ee22 = ee21 + ee19/ee16;
       //   ee23 = R_pow(ee13, 2);
       //   ee24 = ee23 + R_pow(epsilon, 2);
       //   ee25 = ee17/2;
       //   ee26 = ee25 - 1;
       //   ee27 = ee17 - 1;
       //   ee28 = R_pow(ee24, ee25);
       //   ee29 = R_pow(ee24, ee26);
       //   ee30 = R_pow(ee16, 2);
       //   ee31 = ee18 - ee19;
       //   ee32 = 1 + ee17;
       //   ee33 = R_pow(ee22, ee32);
       //   ee34 = R_pow(ee22, ee27);
       //   ee35 = ee29 * ee13;
       //   ee37 = 2 * ee18 - 2 * ee19;
       //   ee38 = R_pow(ee22, (2 + ee17));
       //   ee39 = log(yi);
       //   ee40 = R_pow(ee22, (2 * ee17));
       //   ee41 = R_pow(ee22, (ee17 - 2));
       //   ee42 = R_pow(ee22, ee17);
       //   ee43 = R_pow(ee6, 2);
       //   ee44 = R_pow(ee5, 2);
       //   ee45 = C * ee30;
       //   ee46 = 1 - ee7;
       //   ee47 = log(ee24);
       //   ee48 = log(ee22);
       //   ee50 = ee28 * ee15 * ee31;
       //   ee51 = ee18 + ee19;
       //   ee53 = ee47/2 - ee48;
       //   ee54 = R_pow(yi, ee46);
       //   ee56 = ee29 * ee34 * ee13;
       //   ee57 = ee35/ee33;
       //   ee58 = R_pow(ee24, (ee25 - 2));
       //   ee59 = 1 + ee17 * ee53;
       //   ee60 = 2/ee5;
       //   ee62 = C * ee33 * ee30;
       //   ee64 = C * ee38 * ee30;
       //   ee65 = ee15/ee16;
       //   ee66 = ee35/ee42;
       //   ee67 = ee38 * ee16;
       //   ee68 = ee66 + ee50/ee62;
       //   ee69 = ee58 * ee23;
       //   ee71 = ee28 * ee37 * ee15;
       //   ee72 = ee41 * ee20;
       //   ee73 = 2 * ee46;
       //   ee74 = ee4 * ee17;
       //   ee75 = ee56 * ee37;
       //   ee78 = ee35 * ee20 * ee17/ee33;
       //   ee80 = ee35 * ee17/ee33;
       //   ee81 = ee57 + 2 * (ee50/ee64);
       //   ee82 = ee69 * ee26;
       //   ee84 = ee41 * ee37;
       //   ee85 = ee41 * ee15;
       //   ee86 = 1 - 2 * ee65;
       //   ee88 = 1/ee17;
       //   ee90 = 2 * ee57 + 2 * (ee71/ee64);
       //   ee91 = C * ee16;
       //   ee92 = R_pow(yi, (1 - ee60));
       //   ee93 = R_pow(yi, ee73);
       //   ee94 = R_pow(yi, ee60);
       //   ee98 = (((ee34 * ee86 + ee85 * ee31 * ee27/ee30) * ee28/C -  ee56 * ee17)/ee40 - ee81 * ee17) * ee15 * ee31/ee45 -  (ee29 + 2 * ee82)/ee42;
       //   ee102 = ((ee34 + ee72 * ee31 * ee27)/ee40 - 2 * (ee20 *  ee31 * ee17/ee38)) * ee28 * ee15/ee45 - ee78;
       //   ee106 = ((ee41 * ee31 * ee27/ee16 - ee34)/ee40 - 2 * (ee31 *  ee17/ee67)) * ee28 * ee15/ee91 - ee80;
       //   ee107 = ee68/ee11;
       //   ee108 = ee75 * ee17;
       //   ee109 = ee28 * (ee59/ee33);
       //   ee110 = ee28 * ee20;
       //   ee111 = ee28 * ee17;
       //   ee113 = ee84 * ee15 * ee27;
       //   ee116 = (2 - ee39/ee5) * ee2/ee5 - 1;
       //   ee118 = ee90 * ee31 * ee17;
       //   ee119 = 1/ee51;
       //   ee120 = 2 - 4 * ee65;
       //   ee121 = 2 - 4 * (ee4/ee6);
       //   ee122 = 2 * ee34;
       //   ee123 = R::digamma(ee88);
       //   ee124 = ee54 * ee13;
       //   ee125 = ee94 * ee44;
       //   
       //   out(j, 0) += w * (-(ee74 * (ee54 * (ee71/ee62 + 2 * ee66))/
       //     ee43));
       //   out(j, 1) += w * (-(ee107 * ee2 * ee17 * ee39 * ee10/ee44));
       //   out(j, 2) += w * (-(ee68 * ee17));
       //   out(j, 3) += w * ((ee119 - ee110 * ee17/ee33) * ee18 - 1);
       //   out(j, 4) += w * ((ee119 - ee111/(ee33 * ee16)) * ee19);
       //   out(j, 5) += w * (ee111 * ee53/ee42 - ee123/ee17);
       //   out(j, 6) += w * (-(ee74 * (ee93 * (((ee34 * ((4 - 8 * ee65) *
       //     ee4/(C * ee43) - R_pow(yi, (1 - (ee7 + ee73))) * ee121) *
       //     ee31 + ee41 * R_pow(ee37, 2) * ee15 * ee4 * ee27/(ee45 * ee43)) *
       //     ee28 - 2 * (ee75 * ee4 * ee17/ee43))/ee40 - ee90 * ee37 *
       //     ee4 * ee17/ee43) * ee15/ee45 - (ee29 * (4 * (ee93 * ee4/
       //       ee43) + ee124 * ee121) + 8 * (ee93 * ee58 * ee23 * ee4 * ee26/
       //         ee43))/ee42)/ee43));
       //   out(j, 7) += w * (-(ee4 * ee2 * ee17 * ee39 * (ee92 * (((ee34 *
       //     (ee120 * ee10/C - 2 * ee11) + ee113 * ee10/ee45) * ee28 *
       //     ee31 - ee108 * ee10)/ee40 - ee118 * ee10) * ee15/ee45 - (ee29 *
       //     (2 * ee124 + 2 * (ee92 * ee10)) + 4 * (ee92 * ee58 * ee23 *
       //     ee26 * ee10))/ee42)/(ee43 * ee44)));
       //   out(j, 8) += w * (-(ee74 * (ee54 * ((((ee34 * ee120 + ee113/
       //     ee30) * ee28 * ee31/C - ee108)/ee40 - ee118) * ee15/ee45 - (2 *
       //       ee29 + 4 * ee82)/ee42))/ee43));
       //   out(j, 9) += w * (-(ee4 * ee18 * ee17 * (ee54 * (((ee72 * ee37 *
       //     ee27 + ee122)/ee40 - 2 * (ee20 * ee37 * ee17/ee38)) * ee28 *
       //     ee15/ee45 - 2 * ee78))/ee43));
       //   out(j, 10) += w * (-(ee4 * ee19 * ee17 * (ee54 * (((ee84 * ee27/
       //     ee16 - ee122)/ee40 - 2 * (ee37 * ee17/ee67)) * ee28 * ee15/
       //       ee91 - 2 * ee80))/(ee16 * ee43)));
       //   out(j, 11) += w * (-(ee74 * (ee54 * (ee35 * (2 + ee17 * (ee47 -
       //     2 * ee48))/ee42 + ee28 * ee59 * ee37 * ee15/ee62))/ee43));
       //   out(j, 12) += w * (-(((((ee34 * (ee86 * ee2 * ee39 * ee10/(C *
       //     ee44) + ee11 * ee116) + ee85 * ee2 * ee31 * ee27 * ee39 *
       //     ee10/(ee45 * ee44)) * ee28 - ee56 * ee2 * ee17 * ee39 * ee10/
       //       ee44)/ee40 - ee81 * ee2 * ee17 * ee39 * ee10/ee44)/ee94 *
       //         ee15 * ee31/ee45 + ((ee116 * ee13/ee11 - ee2 * ee39 * ee10/
       //           ee125) * ee29 - 2 * (ee69 * ee2 * ee26 * ee39 * ee10/ee125))/
       //             ee42) * ee2 * ee17 * ee39 * ee10/ee44));
       //   out(j, 13) += w * (-(ee98/ee11 * ee2 * ee17 * ee39 * ee10/ee44));
       //   out(j, 14) += w * (-(ee102/ee11 * ee2 * ee18 * ee17 * ee39 *
       //     ee10/ee44));
       //   out(j, 15) += w * (-(ee106/ee11 * ee2 * ee19 * ee17 * ee39 *
       //     ee10/(ee16 * ee44)));
       //   out(j, 16) += w * (-(ee107 * ee59 * ee2 * ee17 * ee39 * ee10/
       //     ee44));
       //   out(j, 17) += w * (-(ee98 * ee17));
       //   out(j, 18) += w * (-(ee102 * ee18 * ee17));
       //   out(j, 19) += w * (-(ee106 * ee19 * ee17/ee16));
       //   out(j, 20) += w * (-(ee68 * ee59 * ee17));
       //   out(j, 21) += w * (((1 - ee18/ee51)/ee51 - ((ee34 + ee72 * ee18 *
       //     ee27)/ee40 - 2 * (ee21 * ee17/ee38)) * ee28 * ee20 * ee17) *
       //     ee18);
       //   out(j, 22) += w * (-((1/R_pow(ee51, 2) - ee110 * ee32 * ee17/
       //     ee67) * ee18 * ee19));
       //   out(j, 23) += w * (-(ee109 * ee20 * ee18 * ee17));
       //   out(j, 24) += w * (((1 - ee19/ee51)/ee51 - ((ee34 + ee41 * ee19 *
       //     ee27/ee16)/ee40 - 2 * (ee19 * ee17/ee67)) * ee28 * ee17/
       //       ee16) * ee19);
       //   out(j, 25) += w * (-(ee109 * ee19 * ee17/ee16));
       //   out(j, 26) += w * (ee28 * ((0.5 + ee17 * (ee47/4 - ee48/2)) *
       //     ee47 - ee59 * ee48) * ee17/ee42 + (ee123 + R::trigamma(ee88)/
       //       ee17)/ee17);
       
       double ee2, ee3, ee5, ee6, ee7;
       double ee10, ee11, ee13, ee15, ee16, ee17, ee18, ee19;
       double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
       double ee30, ee31, ee32, ee33, ee34, ee35, ee37, ee38, ee39;
       double ee40, ee41, ee42, ee43, ee44, ee45, ee46, ee47, ee48, ee49;
       double ee50, ee51, ee52, ee53, ee54, ee55, ee56, ee57, ee58;
       double ee60, ee61, ee63, ee64, ee65, ee67, ee68, ee69;
       double ee70, ee71, ee72, ee73, ee74, ee75, ee76, ee77, ee78, ee79;
       double ee80, ee81, ee82, ee83, ee84, ee85, ee88, ee89;
       double ee90, ee92, ee93, ee94, ee95, ee96, ee97, ee98, ee99;
       double ee100, ee101, ee102, ee103, ee104, ee106, ee108, ee109;
       double ee110, ee111, ee114, ee115, ee116, ee119;
       double ee120, ee121, ee122, ee123, ee128, ee129;
       double ee130, ee131, ee132, ee133, ee134, ee135, ee136, ee137, ee138, ee139;
       double ee140, ee142, ee143, ee144, ee145, ee146, ee147, ee148;
       double ee150, ee151, ee152, ee154, ee156, ee159;
       double ee162, ee165, ee168, ee169;
       double ee170, ee171, ee173, ee174, ee176, ee177, ee178;
       double ee180, ee181, ee182, ee183, ee185, ee187, ee188, ee189;
       double ee190, ee192, ee194, ee196, ee198, ee199;
       double ee202, ee204, ee208, ee209;
       double ee210, ee211, ee212, ee214, ee215, ee217;
       double ee220, ee221, ee222, ee223, ee224, ee226, ee227, ee228, ee229;
       double ee230, ee232, ee233, ee234, ee239;
       double ee240, ee241, ee243, ee244, ee245, ee246, ee247, ee249;
       double ee250, ee251, ee252, ee254, ee255, ee256, ee257, ee258, ee259;
       double ee260, ee262, ee263, ee264, ee266, ee267, ee268, ee269;
       double ee270, ee271, ee272, ee273, ee274, ee275, ee276, ee278, ee279;
       double ee280, ee282, ee283, ee284, ee285, ee286, ee287, ee288, ee289;
       double ee291, ee292, ee293, ee294, ee298;
       double ee301, ee302, ee303, ee304, ee307, ee309;
       double ee310, ee311, ee312, ee314, ee315, ee317, ee318, ee319;
       double ee320, ee321, ee324;
       double ee330, ee333, ee335, ee337, ee338, ee339;
       double ee341, ee343, ee345, ee346, ee347, ee348, ee349;
       double ee351, ee352, ee354, ee356, ee357, ee359;
       double ee361, ee363, ee366, ee367, ee368;
       double ee370, ee371, ee372, ee373, ee376, ee377, ee378, ee379;
       double ee381, ee383, ee384, ee385, ee386, ee388;
       double ee391, ee392, ee393, ee394, ee397, ee398, ee399;
       double ee400, ee401, ee405, ee409;
       double ee412, ee415, ee419;
       double ee423, ee427;
       double ee431, ee435, ee437, ee438;
       double ee446, ee448, ee449;
       double ee451, ee452, ee453, ee456, ee457, ee459;
       double ee461, ee462, ee467, ee469;
       double ee470, ee474, ee475, ee476;
       double ee480, ee482, ee484, ee488, ee489;
       double ee491, ee492, ee494, ee497, ee498;
       double ee500, ee501, ee502, ee503, ee508, ee509;
       double ee510, ee515, ee516, ee517, ee518;
       double ee520, ee522, ee523, ee524, ee526, ee527, ee528;
       double ee530, ee535, ee537, ee539;
       double ee540, ee541, ee542, ee544, ee546, ee547;
       double ee550, ee552, ee553, ee554, ee555, ee556, ee557, ee558, ee559;
       double ee560, ee561, ee562, ee563, ee567, ee568, ee569;
       double ee570, ee571, ee573, ee574;
       
       ee2 = exp(-p2);
       ee3 = 1 + ee2;
       ee5 = exp(-p1);
       ee6 = 1 + ee5;
       ee7 = 1/ee3;
       ee10 = y - yi * (2/ee6 - 1);
       ee11 = R_pow(yi, ee7);
       ee13 = ee10/ee11 - p3;
       ee15 = exp(-(ee13/C));
       ee16 = 1 + ee15;
       ee17 = exp(p6);
       ee18 = exp(p5);
       ee19 = exp(p4);
       ee20 = 1 - 1/ee16;
       ee21 = ee20 * ee19;
       ee22 = ee18/ee16;
       ee23 = ee21 + ee22;
       ee24 = ee17 - 1;
       ee25 = R_pow(ee13, 2);
       ee26 = ee25 + R_pow(epsilon, 2);
       ee27 = ee17/2;
       ee28 = R_pow(ee23, ee24);
       ee29 = ee19 - ee18;
       ee30 = ee17 - 2;
       ee31 = ee27 - 1;
       ee32 = R_pow(ee23, ee30);
       ee33 = R_pow(ee16, 2);
       ee34 = R_pow(ee26, ee31);
       ee35 = log(yi);
       ee37 = 2 * ee19 - 2 * ee18;
       ee38 = log(ee23);
       ee39 = 2 * ee17;
       ee40 = 1 - ee7;
       ee41 = R_pow(ee26, ee27);
       ee42 = R_pow(ee3, 2);
       ee43 = R_pow(ee6, 2);
       ee44 = ee15/ee16;
       ee45 = C * ee33;
       ee46 = log(ee26);
       ee47 = 2/ee3;
       ee48 = ee27 - 2;
       ee49 = 2 * ee40;
       ee50 = R_pow(ee23, ee39);
       ee51 = R_pow(ee26, ee48);
       ee52 = R_pow(yi, ee40);
       ee53 = ee17 * ee38;
       ee54 = 1 - 2 * ee44;
       ee55 = ee32 * ee20;
       ee56 = R_pow(ee23, ee17);
       ee57 = ee35/ee3;
       ee58 = R_pow(ee23, (2 + ee17));
       ee60 = R_pow(yi, (1 - ee47));
       ee61 = R_pow(yi, ee47);
       ee63 = ee34 * ee28 * ee13;
       ee64 = 2 - 4 * (ee5/ee6);
       ee65 = R_pow(yi, ee49);
       ee67 = (2 - ee57) * ee2/ee3;
       ee68 = ee46/2;
       ee69 = ee32 * ee15;
       ee70 = ee67 - 1;
       ee71 = 1 + ee17;
       ee72 = ee32 * ee37;
       ee73 = R_pow(ee23, ee71);
       ee74 = 2 - 4 * ee44;
       ee75 = R_pow(ee23, (3 + ee17));
       ee76 = 2 * ee24;
       ee77 = 4 * ee17;
       ee78 = R_pow(ee23, ee76);
       ee79 = R_pow(ee23, ee77);
       ee80 = R_pow(ee23, (1 + ee39));
       ee81 = ee32 * ee29;
       ee82 = ee51 * ee25;
       ee83 = R_pow(ee23, (ee39 - 1));
       ee84 = ee28 * ee54;
       ee85 = ee28 + ee55 * ee29 * ee24;
       ee88 = ee81 * ee24/ee16 - ee28;
       ee89 = C * ee43;
       ee90 = ee84 + ee69 * ee29 * ee24/ee33;
       ee92 = ee72 * ee15 * ee24;
       ee93 = ee34 * ee13;
       ee94 = ee52 * ee13;
       ee95 = 2 * ee28;
       ee96 = ee17 * (ee68 + ee38);
       ee97 = C * ee16;
       ee98 = ee63 * ee37;
       ee99 = ee82 * ee31;
       ee100 = ee7 + ee49;
       ee101 = 2 * ee53;
       ee102 = ee45 * ee43;
       ee103 = ee2 * ee35;
       ee104 = ee61 * ee42;
       ee106 = R_pow(yi, (1 - ee100));
       ee108 = 1 + ee96;
       ee109 = 2 * (1 + ee101);
       ee110 = ee45 * ee42;
       ee111 = ee60 * ee10;
       ee114 = (4 - 8 * ee44) * ee5/ee89 - ee106 * ee64;
       ee115 = 1 + ee53;
       ee116 = ee65 * ee5;
       ee119 = ee74 * ee10/C - 2 * ee11;
       ee120 = C * ee42;
       ee121 = ee103 * ee10;
       ee122 = ee11 * ee70;
       ee123 = R_pow(ee23, (ee17 - 3));
       ee128 = ee54 * ee2 * ee35 * ee10/ee120 + ee122;
       ee129 = ee34 + 2 * ee99;
       ee130 = R_pow(ee37, 2);
       ee131 = 3/ee3;
       ee132 = ee75 * ee16;
       ee133 = ee28 * ee119;
       ee134 = ee109 + 2 * ee108;
       ee135 = 8 * ee53;
       ee136 = ((ee24 * ee38 + 2) * ee17 - 1) * ee32;
       ee137 = ee98 * ee17;
       ee138 = ee28 * ee128;
       ee139 = ee133 + ee92 * ee10/ee45;
       ee140 = ee28 * ee114;
       ee142 = ee28 * ee74 + ee92/ee33;
       ee143 = ee32 * ee18;
       ee144 = ee134 - ee135;
       ee145 = C * ee61;
       ee146 = ee116/ee43;
       ee147 = ee58 * ee16;
       ee148 = ee138 + ee69 * ee2 * ee29 * ee24 * ee35 * ee10/ee110;
       ee150 = ee140 * ee29 + ee32 * ee130 * ee15 * ee5 * ee24/ee102;
       ee151 = ee20 * ee29;
       ee152 = ee29 * ee17;
       ee154 = 2 * ee94 + 2 * ee111;
       ee156 = 4 * ee146 + ee94 * ee64;
       ee159 = ee55 * ee37 * ee24 + ee95;
       ee162 = ee72 * ee24/ee16 - ee95;
       ee165 = ee70 * ee13/ee11 - ee121/ee104;
       ee168 = ee90 * ee41/C - ee63 * ee17;
       ee169 = ee23 * ee16;
       ee170 = ee28 + ee55 * ee19 * ee24;
       ee171 = ee28 + ee143 * ee24/ee16;
       ee173 = C * ee58 * ee33;
       ee174 = ee17 * ee46;
       ee176 = ee90 * ee34 * ee13;
       ee177 = ee78 * ee20;
       ee178 = ee151 * ee17;
       ee180 = C * ee75 * ee33;
       ee181 = ee19 + ee18;
       ee182 = ee165 * ee34;
       ee183 = ee145 * ee42;
       ee185 = ee60 * ee51 * ee25;
       ee187 = ee65 * ee51 * ee25;
       ee188 = ee129 * ee28;
       ee189 = ee176/C;
       ee190 = ee182 - 2 * (ee82 * ee2 * ee31 * ee35 * ee10/ee104);
       ee192 = ee34 * ee154 + 4 * (ee185 * ee31 * ee10);
       ee194 = ee34 * ee156 + 8 * (ee187 * ee5 * ee31/ee43);
       ee196 = 2 * ee34 + 4 * ee99;
       ee198 = ee68 - 2 * ee38;
       ee199 = ee68 - ee38;
       ee202 = ee139 * ee41 * ee29 - ee137 * ee10;
       ee204 = ee150 * ee41 - 2 * (ee98 * ee5 * ee17/ee43);
       ee208 = ee142 * ee41 * ee29/C - ee137;
       ee209 = ee88 * ee83;
       ee210 = ee93 * ee20;
       ee211 = ee78 * ee15;
       ee212 = ee28 * ee115;
       ee214 = ee70/ee11 + ee121/ee183;
       ee215 = ee16 * ee43;
       ee217 = ee29 * ee30/ee169;
       ee220 = 2 * ee115;
       ee221 = 2 + ee17 * ee199;
       ee222 = 8 * (ee152/ee132);
       ee223 = ee111/C;
       ee224 = R_pow(yi, (1 - ee131));
       ee226 = ee148 * ee41 - ee63 * ee2 * ee17 * ee35 * ee10/ee42;
       ee227 = ee189 - ee188;
       ee228 = R_pow(ee26, (ee27 - 3));
       ee229 = ee41 * ee144;
       ee230 = ee123 * ee20;
       ee232 = ((1.5 + ee17 * (ee46/4 + ee38)) * ee46 + (3 + 9 *  ee53 - (ee109 + 4 * ee108)) * ee38) * ee17 + 1;
       ee233 = ee20 * ee37;
       ee234 = ee31 * ee46;
       ee239 = 2 * ee223 - 2 * ee52;
       ee240 = 2 * ee15;
       ee241 = 2 + ee220;
       ee243 = 4 * (ee116/ee89) - ee52 * ee64;
       ee244 = 4 * ee38;
       ee245 = ee18 * ee17;
       ee246 = R_pow(yi, (1 + ee49 - ee7));
       ee247 = R_pow(yi, (ee49 - ee7));
       ee249 = (ee90 * ee56 + ee211 * ee29 * ee17/ee33) * ee41 *  ee28;
       ee250 = (ee85 * ee56 + ee177 * ee29 * ee17) * ee28;
       ee251 = (ee171 * ee56 + ee78 * ee18 * ee17/ee16) * ee28;
       ee252 = (ee88 * ee56 + ee78 * ee29 * ee17/ee16) * ee28;
       ee254 = ee142 * ee34 * ee13;
       ee255 = ee85 * ee83;
       ee256 = ee85/ee50;
       ee257 = ee162 * ee83;
       ee258 = ee88/ee50;
       ee259 = ee93 * ee221;
       ee260 = ee93/ee58;
       ee262 = ee41 * ee37 * ee15;
       ee263 = (ee234 + 1) * ee17;
       ee264 = ee178/ee58;
       ee266 = ee16 * ee42;
       ee267 = ee43 * ee42;
       ee268 = ee152/ee147;
       ee269 = 1 + ee240;
       ee270 = 1 + ee174/2;
       ee271 = 2 * ee209;
       ee272 = ee241 + ee17 * (ee46 - ee244);
       ee273 = 2 + ee96;
       ee274 = 8 * (ee178/ee75);
       ee275 = ee5 * ee19;
       ee276 = R_pow(yi, ee131);
       ee278 = (ee170 * ee56 + ee177 * ee19 * ee17) * ee28;
       ee279 = ee129/ee73;
       ee280 = (ee78 * ee17 + R_pow(ee23, (ee39 - 2)) * ee24) * ee28;
       ee282 = ee254 * ee29/C;
       ee283 = ee159 * ee83;
       ee284 = ee136 * ee20;
       ee285 = ee93 * ee272;
       ee286 = ee34 * ee270;
       ee287 = ee51 * ee154;
       ee288 = ee51 * ee156;
       ee289 = ee228 * ee25;
       ee291 = ee41 * ee15 * ee29;
       ee292 = ee28 * ee196;
       ee293 = ee32 + ee230 * ee19 * ee30;
       ee294 = ee32 + ee123 * ee18 * ee30/ee16;
       ee298 = ee123 * ee29 * ee30/ee16 - ee32;
       ee301 = (ee217 - 6) * ee15/ee16 + 3;
       ee302 = ee20 * ee144;
       ee303 = ee21 * ee17;
       ee304 = ee37 * ee17;
       ee307 = (2 + ee57) * ee2/ee3 - 1;
       ee309 = 2 * (ee168/ee80);
       ee310 = 2 * ee81;
       ee311 = 2 * ee16;
       ee312 = 4 * (ee103/ee42);
       ee314 = C * ee73 * ee33;
       ee315 = ee275 * ee17;
       ee317 = ee5 * ee18 * ee17;
       ee318 = ee5 * ee17;
       ee319 = ee174/4;
       ee320 = ee11 * ee42;
       ee321 = R_pow(yi, (2 - ee131));
       ee324 = ee227/ee50;
       ee330 = (ee256 - 2 * ee264) * ee34 * ee13 * ee15;
       ee333 = (ee258 - 2 * ee268) * ee34 * ee13 * ee15;
       ee335 = ee148 * ee34 * ee13;
       ee337 = ee85 * ee34 * ee13;
       ee338 = ee170/ee50;
       ee339 = ee171/ee50;
       ee341 = ee159 * ee34 * ee13;
       ee343 = ee162 * ee34 * ee13;
       ee345 = ee88 * ee34 * ee13;
       ee346 = ee232 * ee41;
       ee347 = ee214 * ee16;
       ee348 = ee136 * ee37;
       ee349 = ee136 * ee15;
       ee351 = (ee263 + 2 * ee31) * ee51 * ee25;
       ee352 = ee63 * ee273;
       ee354 = ee285/ee73 + ee229 * ee37 * ee15/ee173;
       ee356 = ee259/ee73 + ee229 * ee15 * ee29/ee173;
       ee357 = ee289 * ee48;
       ee359 = ee41 * ((2 * ee250 + 2 * ee255)/ee79 - ee274) *  ee15;
       ee361 = ee41 * ((2 * ee252 + ee271)/ee79 - ee222) * ee15;
       ee363 = ee262 * ee29 * ee17;
       ee366 = (ee151 * ee30/ee23 + 2) * ee32 * ee15;
       ee367 = ee80 * ee16;
       ee368 = ee78 * ee37;
       ee370 = ee32 * (ee217 - 2) * ee15;
       ee371 = ee32 * ee19;
       ee372 = ee32/ee16;
       ee373 = ee123 * ee15;
       ee376 = ee233 * ee17;
       ee377 = ee303/ee58;
       ee378 = ee16 * ee239;
       ee379 = (2 * (ee168 * ee83) + 2 * (ee249/C))/ee79;
       ee381 = ee37 * ee30/ee169;
       ee383 = 1/ee17;
       ee384 = 2 * (ee226/ee80);
       ee385 = 2 * (ee202/ee80);
       ee386 = 2 * (ee208 * ee83);
       ee388 = 2 * ee260 + 8 * (ee291/ee180);
       ee391 = 2 * ee269 + ee311 - 8 * ee15;
       ee392 = 4 * (ee25 * ee48/ee26);
       ee393 = 8 * (ee304/ee132);
       ee394 = ee145 * ee33;
       ee397 = ee15 * ee2 * ee35 * ee10;
       ee398 = ee5 * ee2;
       ee399 = ee245/ee147;
       ee400 = ee245/ee132;
       ee401 = ee60 * ee34;
       ee405 = ((((ee301 * ee15 * ee29 * ee24/ee33 + ee23 * (1 -  (ee391/ee16 + 2) * ee15/ee16)) * ee41 * ee32/C - ee176 *  ee17)/C - ee227 * ee17)/ee50 - (ee324 + (ee379 + ee309 -  ee388 * ee17) * ee15 * ee29/ee45 - ee279) * ee17) * ee15 *  ee29/ee45 + ee51 * ee13 * (ee392 + 6) * ee31/ee56;
       ee409 = (((ee90 * ee17 * ee46/2 + ee349 * ee29/ee33 + ee84 *  ee115) * ee41/C - (ee352 + 2 * (ee168 * ee38)) * ee17)/ee50 -  ee356 * ee17) * ee15 * ee29/ee45 - (ee351 + ee286 - ee129 *  ee17 * ee38)/ee56;
       ee412 = (((ee85 * ee54 + ee366 * ee29 * ee24/ee33) * ee41/C -  ee337 * ee17)/ee50 - (ee359/ee45 + 2 * (ee168 * ee20/ee80)) *  ee29 * ee17) * ee15/ee45 - (ee330/ee45 - ee129 * ee20/ee73) *  ee17;
       ee415 = (((ee88 * ee54 + ee370 * ee29 * ee24/ee33) * ee41/C -  ee345 * ee17)/ee50 - (ee361/ee97 + ee309) * ee29 * ee17/ee16) *  ee15/ee97 - (ee333/ee97 - ee279) * ee17;
       ee419 = (((ee293 * ee29 + 2 * ee371) * ee20 * ee24 + ee28)/ee50 -  ((2 * (ee278 * ee29) + 2 * (ee255 * ee19))/ee79 + (2 *  (ee85/ee80) - ee274) * ee19) * ee20 * ee17) * ee41 *  ee15/ee45 - (ee338 - 2 * ee377) * ee34 * ee13 * ee20 *  ee17;
       ee423 = (((ee294 * ee29 - 2 * ee143) * ee24/ee16 - ee28)/ee50 -  ((2 * (ee251 * ee29) + 2 * (ee209 * ee18))/ee79 + (2 *  (ee88/ee80) - ee222) * ee18) * ee17/ee16) * ee41 *  ee15/ee97 - (ee339 - 2 * ee399) * ee34 * ee13 * ee17;
       ee427 = ((ee85 * ee17 * ee198 + ee284 * ee29 + ee212)/ee50 -  ee302 * ee29 * ee17/ee58) * ee41 * ee15/ee45 - ee210 *  ee221 * ee17/ee73;
       ee431 = ((ee88 * ee17 * ee198 + ee136 * ee29/ee16 - ee212)/ee50 -  ee144 * ee29 * ee17/ee147) * ee41 * ee15/ee97 - ee259 *  ee17/ee73;
       ee435 = ((ee298 * ee20 + ee372) * ee24/ee50 - (((2 * (ee280 *  ee29/ee16) + ee271)/ee79 - ee222) * ee20 + 2 * (ee85/ee367)) *  ee17) * ee41 * ee15/ee97 + ee210 * ee71 * ee17/ee58;
       ee437 = (ee148 * ee56 + ee211 * ee2 * ee29 * ee17 * ee35 *  ee10/ee110) * ee41 * ee28;
       ee438 = ee202 * ee83;
       ee446 = ((ee32 + ee230 * ee29 * ee30) * ee37 + ee310) *  ee15 * ee24;
       ee448 = (ee298 * ee37 - ee310) * ee15 * ee24;
       ee449 = ee190/ee73;
       ee451 = ee354 * ee29 * ee17;
       ee452 = ee192/ee73;
       ee453 = ee194 * ee29;
       ee456 = ee150 * ee34 * ee13;
       ee457 = ee282 - ee188 * ee37;
       ee459 = ee341 * ee17;
       ee461 = (ee32 * ee54 + ee373 * ee29 * ee30/ee33) * ee37;
       ee462 = ee343 * ee17;
       ee467 = ee346 * ee15 * ee29/ee314 + (((1.5 + ee319) * ee46 -  (3 + ee17 * (ee46 - ee38)) * ee38) * ee17 + 1) * ee34 *  ee13/ee56;
       ee469 = ee232/ee73 * ee41;
       ee470 = ee347 + ee397/ee183;
       ee474 = ee301 * ee2 * ee35 * ee10/ee120;
       ee475 = ee348 * ee15;
       ee476 = ee98 * ee273;
       ee480 = ee93 * ee15 * ee29 * ee17/ee173;
       ee482 = ee41 * ((2 * (ee250 * ee37) + 2 * (ee283 * ee29))/ee79 -  8 * (ee233 * ee29 * ee17/ee75)) * ee15;
       ee484 = ee41 * ((2 * (ee252 * ee37) + 2 * (ee257 * ee29))/ee79 -  8 * (ee37 * ee29 * ee17/ee132)) * ee15;
       ee488 = ee368 * ee15 * ee17;
       ee489 = R_pow(ee23, (3 * ee17 - 3));
       ee491 = ee28 * (2 + ee101);
       ee492 = ee32 * ee119;
       ee494 = ee32 * ee74 * ee29;
       ee497 = (ee381 - 8) * ee15/ee16 + 4;
       ee498 = ee114 * ee29;
       ee500 = ee16 * ee243;
       ee501 = (2 * (ee249 * ee37/C) + ee386)/ee79;
       ee502 = (2 * ee189 - ee292)/ee50;
       ee503 = ee196/ee73;
       ee508 = (6 + 9 * ee57) * ee2/ee3;
       ee509 = R_pow(ee181, 2);
       ee510 = 0.5 + ee319;
       ee515 = 16 * ee5;
       ee516 = 2 * (ee204/ee80);
       ee517 = 2 * (ee208/ee80);
       ee518 = 2 * ee251;
       ee520 = 2 * ee282 - ee292 * ee37;
       ee522 = 2 * ee256 - 4 * ee264;
       ee523 = 2 * (ee170/ee80);
       ee524 = 2 * ee257;
       ee526 = 2 * ee258 - 4 * ee268;
       ee527 = 2 * ee287;
       ee528 = 2 * ee288;
       ee530 = 2 * ee263 + 4 * ee31;
       ee535 = 2 * ((2 + 2 * ee57) * ee2/ee3 - 1) + 3 + 4 * ee307;
       ee537 = 2 * ee307 - ee312;
       ee539 = 2 * (1 + 2 * ee2) + 2 * ee3;
       ee540 = 2 * (2 * ee51 + 4 * ee357);
       ee541 = 2 * (ee19/ee181);
       ee542 = 2 * (ee18/ee181);
       ee544 = 2 * (ee321 * ee139 * ee34 * ee13 * ee29) - ee52 *  ee192 * ee28 * ee37;
       ee546 = 4 * ee260 + 8 * (ee262/ee180);
       ee547 = 4 * ee51;
       ee550 = 4 * ee269 + 4 * ee16 - 16 * ee15;
       ee552 = 4 * (1 + 2 * ee5) + 4 * ee6;
       ee553 = 8 * (ee363/ee180);
       ee554 = 8 * ee51;
       ee555 = 8 * (ee376/ee75);
       ee556 = 8 * ee400;
       ee557 = 8 * ee2;
       ee558 = ee97 * ee43;
       ee559 = ee145 * ee16;
       ee560 = ee394 * ee42;
       ee561 = C * ee276;
       ee562 = ee52 * ee34;
       ee563 = ee401 * ee13;
       ee567 = ee224 * ee139 * ee34 * ee13 * ee2;
       ee568 = ee246 * ee32;
       ee569 = R_pow(yi, (ee100 - 1));
       ee570 = ee11 * ee214;
       ee571 = ee11 * ee73;
       ee573 = ee65 * ee34 * ee13;
       ee574 = ee61 * ee190;
       
       out(j, 0) += -((((((ee28 * (ee52 * ((((ee569 * (2 * (ee106 *
         ee243) - 4 * ee64) - 4 * (ee243 * ee15/ee16))/C + ee515 - ee552)/
           ee6 - 4) * ee5/ee6 + 2) - ee15 * ee5 * (ee52 * (4 * (ee500 +
             4 * (ee65 * ee15 * ee5/ee89)) + 4 * ee500) - 64 * (ee246 *
             ee15 * ee5/ee89))/ee102) + ee568 * ee114 * ee37 * ee15 *
             ee5 * ee24/ee102) * ee29 + ee568 * (ee130 * ee15 * ee5 *
             ee30/(C * ee23 * ee33 * ee43) + 2 * ee498) * ee37 * ee15 * ee5 *
             ee24/ee102) * ee41 - ee246 * (4 * ee456 - ee194 * ee28 *
             ee37/ee65) * ee5 * ee17/ee43)/ee50 - ((2 * (ee246 * ee150 *
             ee34 * ee13) - ee52 * ee194 * ee28 * ee37)/ee50 + ee246 *
             ((2 * (ee204 * ee83) + 2 * ((ee150 * ee56 + ee78 * ee130 * ee15 *
             ee5 * ee17/ee102) * ee41 * ee28))/ee79 + ee516 - ee37 *
             ee546 * ee5 * ee17/ee43) * ee37 * ee15/ee45) * ee5 * ee17/
               ee43) * ee15/ee45 + ee52 * (ee194 * ee37 * ee15 * ee5 * ee17/
                 (ee314 * ee43) + (ee13 * (2 * (16 * (ee65 * ee228 * ee25 *
                   ee5 * ee48/ee43) + ee528) + 8 * ee288) * ee5 * ee31/ee43 + ee569 *
                   ee34 * (6 * (ee64 * ee5/ee43) + ee106 * ee13 * (2 - ((ee552 -
                   ee515)/ee6 + 4) * ee5/ee6)))/ee56)) * ee5 * ee17/ee43);
       out(j, 1) += -((((((((2 * (ee52 * ee239) - 4 * ee65) * ee5/
         ee43 - (2 * (ee243 * ee10/ee11) + ee5 * (ee52 * (4 * (ee378 +
           2 * (ee60 * ee15 * ee10/C)) + 4 * ee378) - 32 * (ee247 * ee15 *
           ee10/C))/ee215) * ee15/ee16)/C - ee64 * (ee223 - ee52)) *
           ee23 + ee15 * ee24 * (ee321 * (ee497 * ee10/C - 4 * ee11) *
           ee37 * ee5/ee43 + ee247 * ee114 * ee29 * ee10)/ee45) * ee41 *
           ee32 * ee29 - (ee544 * ee5/ee43 + ee247 * ee150 * ee34 *
           ee13 * ee10) * ee17)/ee50 - ((ee544/ee50 - 4 * (ee247 * ee34 *
           ee13 * ee37 * ee15 * ee29 * ee17 * ee10/ee173)) * ee5/
             ee43 + (ee37 * ee5 * (ee321 * (2 * ee438 + 2 * ((ee139 * ee56 +
               ee488 * ee10/ee45) * ee41 * ee28 * ee29))/ee79 - 8 * (ee247 *
               ee41 * ee37 * ee15 * ee29 * ee17 * ee10/ee180))/ee43 +
               2 * (ee247 * ee204 * ee29 * ee10/ee80)) * ee15/ee45 - ee453 *
               ee10/ee571) * ee17) * ee15/ee45 + (ee34 * (ee64 * (ee94 +
               ee111) + 8 * ee146) + ee13 * (2 * (ee288 * ee10/ee11) + ee52 *
               (2 * (ee527 + 8 * (ee60 * ee228 * ee25 * ee48 * ee10)) +
               4 * ee287) * ee5/ee43) * ee31)/ee56) * ee5 * ee2 * ee17 * ee35/
                 ee267);
       out(j, 2) += -(((((ee23 * (ee54 * ee243 - ee65 * (4 * (ee311 +
         ee240) + 8 * ee16 - 32 * ee15) * ee15 * ee5/ee102) + ee65 *
         (ee497 * ee37 * ee5/ee89 + ee498) * ee15 * ee24/ee33) * ee41 *
         ee32 * ee29/C - ee65 * (ee456 + ee520 * ee5/ee43) * ee17)/
           ee50 - ee17 * (ee65 * ((((ee386 + 2 * ((ee142 * ee56 + ee488/
             ee33) * ee41 * ee28 * ee29/C))/ee79 - ee553) * ee37 * ee5/
               ee43 + 2 * (ee204 * ee29/ee80)) * ee15/ee45 + (ee520/ee50 -
                 4 * (ee93 * ee37 * ee15 * ee29 * ee17/ee173)) * ee5/ee43) -
                 ee453/ee73)) * ee15/ee45 + (ee13 * (ee528 + ee65 * (2 * (ee547 +
                 8 * ee357) + ee554) * ee5/ee43) * ee31 + ee562 * ee64)/
                   ee56) * ee5 * ee17/ee43);
       out(j, 3) += -(ee315 * (ee65 * (((ee85 * ee114 + (ee233 * ee30/
         ee23 + 4) * ee32 * ee37 * ee15 * ee5 * ee24/ee102) * ee41 -
           2 * (ee341 * ee5 * ee17/ee43))/ee50 - (ee41 * ((2 * ((ee159 *
           ee56 + ee177 * ee37 * ee17) * ee28) + 2 * ee283)/ee79 -
           ee555) * ee37 * ee15 * ee5/ee102 + 2 * (ee204 * ee20/ee80)) *
           ee17) * ee15/ee45 - ee17 * (ee573 * (2 * (ee159/ee50) - 4 *
           (ee376/ee58)) * ee15 * ee5/ee102 - ee194 * ee20/ee73))/ee43);
       out(j, 4) += -(ee317 * (ee65 * (((ee88 * ee114 + ee32 * (ee381 -
         4) * ee37 * ee15 * ee5 * ee24/ee102) * ee41 - 2 * (ee343 *
         ee5 * ee17/ee43))/ee50 - (ee41 * ((2 * ((ee162 * ee56 +
         ee368 * ee17/ee16) * ee28) + ee524)/ee79 - ee393) * ee37 * ee15 *
         ee5/ee558 + ee516) * ee17/ee16) * ee15/ee97 - ee17 * (ee573 *
         (2 * (ee162/ee50) - 4 * (ee304/ee147)) * ee15 * ee5/
           ee558 - ee194/ee73))/ee215);
       out(j, 5) += -(ee318 * (ee65 * (((ee150 * ee17 * ee46/2 + ee136 *
         ee130 * ee15 * ee5/ee102 + ee140 * ee115 * ee29) * ee41 -
         (ee98 * (ee241 + ee174) * ee5/ee43 + 2 * (ee204 * ee38)) *
         ee17)/ee50 - ee354 * ee37 * ee5 * ee17/ee43) * ee15/ee45 -
         (ee286 * ee156 + ee187 * (2 * ((2 + 2 * ee234) * ee17) + 8 *
         ee31) * ee5/ee43 - ee194 * ee17 * ee38)/ee56)/ee43);
       out(j, 6) += -((((((ee28 * (ee10 * (ee60 * (2 * ee570 - ee312) -
         (((2 * (ee378/ee11) - 16 * (ee224 * ee15 * ee10/C)) * ee2 *
         ee35/ee42 + 4 * (ee52 * ee470))/ee16 + 2 * (ee239 * ee2 *
         ee35/ee320)) * ee15/ee16)/C - ee52 * ee537) + ee224 * ((ee32 *
         ee128 + ee373 * ee2 * ee29 * ee30 * ee35 * ee10/ee110) *
         ee37 + 2 * (ee492 * ee2 * ee29 * ee35/ee42)) * ee15 * ee24 *
         ee10/ee45) * ee41 - ee567 * ee17 * ee35 * ee10/ee42) * ee29 -
         ee17 * ee10 * (ee52 * ee190 * ee28 * ee37 + ee567 * ee29 *
         ee35/ee42))/ee50 - ((2 * (ee224 * ee148 * ee34 * ee13) -
         ee192 * ee28 * ee2 * ee35/ee320)/ee50 + ee224 * ((2 * (ee437 *
         ee37) + 2 * (ee438 * ee2 * ee35/ee42))/ee79 + (ee385 -
         8 * (ee363 * ee10/ee180)) * ee2 * ee35/ee42) * ee15/ee45 - (ee192/
           ee571 + 4 * (ee224 * ee34 * ee13 * ee15 * ee29 * ee17 *
             ee10/ee173)) * ee2 * ee35/ee42) * ee29 * ee17 * ee10) * ee15/
               ee45 - (ee13 * (2 * (ee52 * (2 * (ee165 * ee51) - 4 * (ee289 *
                 ee2 * ee48 * ee35 * ee10/ee104))) - 4 * (ee287 * ee2 *
                 ee35/ee320)) * ee31 * ee10 + ee401 * ((2 * ee70 - ee312) *
                 ee10 + ee11 * ee13 * ee537))/ee56) * ee5 * ee2 * ee17 * ee35/
                   ee267);
       out(j, 7) += -((((((ee28 * (ee54 * ee239 - ee60 * ee550 * ee15 *
         ee10/ee45) + ee60 * ((ee461 + ee494) * ee10/C + ee492 *
         ee29) * ee15 * ee24/ee33) * ee41/C - ee60 * ee139 * ee34 *
         ee13 * ee17) * ee29 - ee60 * ee457 * ee17 * ee10)/ee50 - ee152 *
         (ee60 * (((ee501 - ee553) * ee10 + ee385) * ee15/ee45 +
         (ee502 - 4 * ee480) * ee10) - ee452)) * ee15/ee45 + (ee13 *
         (ee527 + ee60 * (ee540 + ee547) * ee10) * ee31 + 2 * ee562)/
           ee56) * ee5 * ee2 * ee17 * ee35/ee267);
       out(j, 8) += -(ee398 * ee19 * ee17 * ee35 * (ee60 * (((ee446 *
         ee10/ee45 + ee85 * ee119) * ee41 - ee459 * ee10)/ee50 - (ee482 *
         ee10/ee45 + 2 * (ee202 * ee20/ee80)) * ee17) * ee15/
           ee45 - ee17 * (ee563 * ee522 * ee15 * ee10/ee45 - ee192 * ee20/
             ee73))/ee267);
       out(j, 9) += -(ee398 * ee18 * ee17 * ee35 * (ee60 * (((ee448 *
         ee10/ee45 + ee88 * ee119) * ee41 - ee462 * ee10)/ee50 - (ee484 *
         ee10/ee97 + ee385) * ee17/ee16) * ee15/ee97 - ee17 *
         (ee563 * ee526 * ee15 * ee10/ee97 - ee452))/(ee215 * ee42));
       out(j, 10) += -(ee398 * ee17 * ee35 * (ee60 * (((ee139 * ee17 *
         ee46/2 + ee475 * ee10/ee45 + ee133 * ee115) * ee41 * ee29 -
         (ee476 * ee10 + 2 * (ee202 * ee38)) * ee17)/ee50 - ee451 *
         ee10) * ee15/ee45 - (ee286 * ee154 + ee185 * ee530 * ee10 -
         ee192 * ee17 * ee38)/ee56)/ee267);
       out(j, 11) += -(ee318 * (ee52 * ((((((ee461 + 2 * ee494) * ee15 *
         ee24/ee33 + ee28 * (2 - (ee550/ee16 + 4) * ee15/ee16)) *
         ee41/C - ee254 * ee17) * ee29/C - ee457 * ee17)/ee50 - ((ee501 +
         ee517 - ee546 * ee29 * ee17) * ee15/ee45 +   ee502 -
         ee503) * ee29 * ee17) * ee15/ee45 + ee13 * (ee540 + ee554) *
         ee31/ee56))/ee43);
       out(j, 12) += -(ee315 * (ee52 * ((((ee446/ee33 + ee85 * ee74) *
         ee41/C - ee459)/ee50 - (ee482/ee45 + 2 * (ee208 * ee20/
           ee80)) * ee17) * ee15/ee45 - (ee93 * ee522 * ee15/ee45 - ee20 *
             ee196/ee73) * ee17))/ee43);
       out(j, 13) += -(ee317 * (ee52 * ((((ee448/ee33 + ee88 * ee74) *
         ee41/C - ee462)/ee50 - (ee484/ee97 + ee517) * ee17/ee16) *
         ee15/ee97 - (ee93 * ee526 * ee15/ee97 - ee503) * ee17))/
           ee215);
       out(j, 14) += -(ee318 * (ee52 * ((((ee142 * ee17 * ee46/2 +
         ee475/ee33 + ee212 * ee74) * ee41 * ee29/C - (ee476 + 2 * (ee208 *
         ee38)) * ee17)/ee50 - ee451) * ee15/ee45 - (ee34 * (2 +
         ee174) + ee82 * ee530 - ee196 * ee17 * ee38)/ee56))/ee43);
       out(j, 15) += -(ee315 * (ee52 * ((((ee293 * ee37 + 4 * ee371) *
         ee20 * ee24 + ee95)/ee50 - ((2 * (ee278 * ee37) + 2 * (ee283 *
         ee19))/ee79 + (2 * (ee159/ee80) - ee555) * ee19) * ee20 *
         ee17) * ee41 * ee15/ee45 - ee210 * (2 * ee338 - 4 * ee377) *
         ee17))/ee43);
       out(j, 16) += -(ee275 * ee18 * ee17 * (ee52 * ((((ee123 * ee37 *
         ee30/ee16 - 2 * ee32) * ee20 + 2 * ee372) * ee24/ee50 -
         (((2 * (ee280 * ee37/ee16) + ee524)/ee79 - ee393) * ee20 +
         2 * (ee159/ee367)) * ee17) * ee41 * ee15/ee97 - ee210 * (ee76 -
         ee77) * ee17/ee58))/ee215);
       out(j, 17) += -(ee315 * (ee52 * (((ee159 * ee17 * ee198 + ee284 *
         ee37 + ee491)/ee50 - ee302 * ee37 * ee17/ee58) * ee41 *
         ee15/ee45 - ee210 * ee272 * ee17/ee73))/ee43);
       out(j, 18) += -(ee317 * (ee52 * ((((ee294 * ee37 - 4 * ee143) *
         ee24/ee16 - ee95)/ee50 - ((2 * (ee251 * ee37) + 2 * (ee257 *
         ee18))/ee79 + (2 * (ee162/ee80) - ee393) * ee18) * ee17/
           ee16) * ee41 * ee15/ee97 - ee93 * (2 * ee339 - 4 * ee399) *
             ee17))/ee215);
       out(j, 19) += -(ee317 * (ee52 * (((ee162 * ee17 * ee198 + ee348/
         ee16 - ee491)/ee50 - ee144 * ee37 * ee17/ee147) * ee41 *
           ee15/ee97 - ee285 * ee17/ee73))/ee215);
       out(j, 20) += -(ee318 * (ee52 * (ee346 * ee37 * ee15/ee314 +
         (((2 + 2 * ee510) * ee46 - ((2 * ee46 - ee244) * ee17 + ee220 +
         4) * ee38) * ee17 + 2) * ee34 * ee13/ee56))/ee43);
       out(j, 21) += -((((((((((((2 * ee67 + ee570 - 2)/ee11 - 2 *
         (ee214 * ee15/ee16)) * ee10/C + ee508 - ee535) * ee35 + ee557 -
         ee539)/ee3 - 2) * ee2/ee3 + 1)/ee11 - ((2 * ee470 + 2 *
         ee347)/ee11 - 8 * (ee397/(ee561 * ee42))) * ee15 * ee2 * ee35 *
         ee10/ee110) * ee23 + (ee474 + 3 * ee122) * ee15 * ee2 *
         ee29 * ee24 * ee35 * ee10/(ee561 * ee33 * ee42)) * ee41 * ee32 -
         (2 * ee335 + ee574 * ee28) * ee2 * ee17 * ee35 * ee10/
           (ee276 * ee42))/ee50 - ((ee190 * ee28/ee11 + ee335/ee276)/ee50 +
             (((2 * (ee226 * ee83) + 2 * ee437)/ee79 + ee384 - ee388 *
             ee2 * ee17 * ee35 * ee10/ee42) * ee15 * ee29/ee45 + ee574/
               ee73)/ee276) * ee2 * ee17 * ee35 * ee10/ee42) * ee15 * ee29/
                 ee45 + (ee26 * (ee11 * ((((ee508 - ee535) * ee35 + ee557 - ee539)/
                   ee3 - 2) * ee2/ee3 + 1) * ee13 - 3 * (ee70 * ee2 * ee35 *
                     ee10/ee42))/ee61 - ee13 * (6 * ee165 - 4 * (ee25 * ee2 *
                     ee48 * ee35 * ee10/(ee61 * ee26 * ee42))) * ee2 * ee31 * ee35 *
                     ee10/ee320) * ee51/ee56) * ee2 * ee17 * ee35 * ee10/ee42);
       out(j, 22) += -((((((ee214 * ee54 - ee391 * ee15 * ee2 * ee35 *
         ee10/ee560) * ee23 + (ee474 + ee122) * ee15 * ee29 * ee24/
           (ee61 * ee33)) * ee41 * ee32/C - (ee227 * ee2 * ee35 * ee10/
             ee42 + ee335) * ee17/ee61)/ee50 - (((ee324 - 2 * ee480) *
               ee2 * ee35 * ee10/ee42 + ((ee379 - 8 * (ee291 * ee17/ee180)) *
               ee2 * ee35 * ee10/ee42 + ee384) * ee15 * ee29/ee45)/ee61 +
               ee449) * ee17) * ee15 * ee29/ee45 - (ee26 * ee70/ee11 + ee13 *
               (2 * ee165 - (4 + ee392) * ee2 * ee35 * ee10/ee104) * ee31) *
               ee51/ee56) * ee2 * ee17 * ee35 * ee10/ee42);
       out(j, 23) += -(((((ee85 * ee128 + ee366 * ee2 * ee29 * ee24 *
         ee35 * ee10/ee110) * ee41 - ee337 * ee2 * ee17 * ee35 * ee10/
           ee42)/ee50 - (ee359 * ee2 * ee35 * ee10/ee110 + 2 * (ee226 *
             ee20/ee80)) * ee29 * ee17) * ee15/ee394 - (ee330 * ee2 *
             ee35 * ee10/ee560 + ee190 * ee20/ee73) * ee17) * ee2 * ee19 *
             ee17 * ee35 * ee10/ee42);
       out(j, 24) += -(((((ee88 * ee128 + ee370 * ee2 * ee29 * ee24 *
         ee35 * ee10/ee110) * ee41 - ee345 * ee2 * ee17 * ee35 * ee10/
           ee42)/ee50 - (ee361 * ee2 * ee35 * ee10/(ee97 * ee42) +
             ee384) * ee29 * ee17/ee16) * ee15/ee559 - (ee333 * ee2 * ee35 *
             ee10/(ee559 * ee42) + ee449) * ee17) * ee2 * ee18 * ee17 *
             ee35 * ee10/ee266);
       out(j, 25) += -(((((ee148 * ee17 * ee46/2 + ee349 * ee2 * ee29 *
         ee35 * ee10/ee110 + ee138 * ee115) * ee41 - (ee352 * ee2 *
         ee35 * ee10/ee42 + 2 * (ee226 * ee38)) * ee17)/ee50 - ee356 *
         ee2 * ee17 * ee35 * ee10/ee42)/ee61 * ee15 * ee29/ee45 +
         (ee182 * ee270 - (ee190 * ee17 * ee38 + ee351 * ee2 * ee35 *
         ee10/ee104))/ee56) * ee2 * ee17 * ee35 * ee10/ee42);
       out(j, 26) += -(ee405/ee11 * ee2 * ee17 * ee35 * ee10/ee42);
       out(j, 27) += -(ee412/ee11 * ee2 * ee19 * ee17 * ee35 * ee10/
         ee42);
       out(j, 28) += -(ee415/ee11 * ee2 * ee18 * ee17 * ee35 * ee10/
         ee266);
       out(j, 29) += -(ee409/ee11 * ee2 * ee17 * ee35 * ee10/ee42);
       out(j, 30) += -(ee419/ee11 * ee2 * ee19 * ee17 * ee35 * ee10/
         ee42);
       out(j, 31) += -(ee435/ee11 * ee2 * ee19 * ee18 * ee17 * ee35 *
         ee10/ee266);
       out(j, 32) += -(ee427/ee11 * ee2 * ee19 * ee17 * ee35 * ee10/
         ee42);
       out(j, 33) += -(ee423/ee11 * ee2 * ee18 * ee17 * ee35 * ee10/
         ee266);
       out(j, 34) += -(ee431/ee11 * ee2 * ee18 * ee17 * ee35 * ee10/
         ee266);
       out(j, 35) += -(ee467/ee11 * ee2 * ee17 * ee35 * ee10/ee42);
       out(j, 36) += -(ee405 * ee17);
       out(j, 37) += -(ee412 * ee19 * ee17);
       out(j, 38) += -(ee415 * ee18 * ee17/ee16);
       out(j, 39) += -(ee409 * ee17);
       out(j, 40) += -(ee419 * ee19 * ee17);
       out(j, 41) += -(ee435 * ee19 * ee18 * ee17/ee16);
       out(j, 42) += -(ee427 * ee19 * ee17);
       out(j, 43) += -(ee423 * ee18 * ee17/ee16);
       out(j, 44) += -(ee431 * ee18 * ee17/ee16);
       out(j, 45) += -(ee467 * ee17);
       out(j, 46) += ((1 - (3 - ee541) * ee19/ee181)/ee181 - ((((ee21 *
         ee30/ee23 + 3) * ee24 + 1) * ee20 * ee19 + ee22)/ee58 -
         ((2 * ee278 + 2 * (ee170 * ee83))/ee79 + ee523 - 8 * (ee303/
           ee75)) * ee20 * ee19 * ee17) * ee41 * ee20 * ee17) * ee19;
       out(j, 47) += -(((ee293 * ee24/ee50 - (((2 * ee280 + 2 * (ee489 *
         ee24))/ee79 - 8 * (ee17/ee75)) * ee20 * ee19 + ee523) *
         ee17) * ee41 * ee20 * ee17/ee16 + (1 - ee541)/ee509) * ee19 *
         ee18);
       out(j, 48) += -(((ee170 * ee17 * ee198 + ee284 * ee19 + ee212)/
         ee50 - ee302 * ee19 * ee17/ee58) * ee41 * ee20 * ee19 * ee17);
       out(j, 49) += -((((ee294/ee50 - 2 * ee400) * ee24 - ((ee518 +
         2 * (ee489 * ee18 * ee24/ee16))/ee79 - ee556) * ee17) * ee41 *
         ee20 * ee17/ee16 + (1 - ee542)/ee509) * ee19 * ee18);
       out(j, 50) += -(((ee24 * ee199 + 2 + ee135 - ee134) * ee17 -
         1)/ee58 * ee41 * ee20 * ee19 * ee18 * ee17/ee16);
       out(j, 51) += -(ee469 * ee20 * ee19 * ee17);
       out(j, 52) += ((1 - (3 - ee542) * ee18/ee181)/ee181 - ((((2 +
         ee18 * ee30/ee169) * ee24 + ee17) * ee18/ee16 + ee21)/ee58 -
         ((ee518 + 2 * (ee171 * ee83))/ee79 + 2 * (ee171/ee80) - ee556) *
         ee18 * ee17/ee16) * ee41 * ee17/ee16) * ee18;
       out(j, 53) += -(((ee171 * ee17 * ee198 + ee136 * ee18/ee16 +
         ee212)/ee50 - ee144 * ee18 * ee17/ee147) * ee41 * ee18 * ee17/
           ee16);
       out(j, 54) += -(ee469 * ee18 * ee17/ee16);
       out(j, 55) += ((((ee510/2 + 0.5) * ee46 - ((0.5 * ee46 - ee38) *
         ee17 + ee115/2 + 1) * ee38) * ee17 + 0.5) * ee46 - ee232 *
         ee38) * ee41 * ee17/ee56 - ((3 * R::trigamma(ee383) + R::psigamma(ee383, 2)/
           ee17)/ee17 + R::digamma(ee383))/ee17;
       }}}
   
   return out;
   
 }

// //' Conditional extreme value model negative log-likelihood for asymmetric generalised Gaussian
 // //'
 // //' @param pars a list of vectors of coefficients for each conditional EVD parameter
 // //' @param X1 a design matrix for (transformed) alpha
 // //' @param X2 a design matrix for (transformed) beta
 // //' @param X3 a design matrix for nu
 // //' @param X4 a design matrix for (transformed) kappa1
 // //' @param X5 a design matrix for (transformed) kappa2
 // //' @param X6 a design matrix for (transformed) delta
 // //' @param ymat a matrix
 // //' @param xmat a matrix
 // //' @param wmat a matrix
 // //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
 // //' @return condexd0 a scalar, the negative log-likelihood
 // //' @return condexd12 a matrix, first then second derivatives w.r.t. parameters
 // //' @return condexd34 a matrix, third then fourth derivatives w.r.t. parameters
 // //' @examples
 // //' ## to follow
 // //' @export
 // // [[Rcpp::export]]
 // double condexaggd0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, 
 //                    arma::mat X4, arma::mat X5, arma::mat X6, 
 //                    arma::mat ymat, arma::mat xmat, arma::mat wmat, 
 //                    arma::uvec dupid, int dcate, arma::uvec nhere)
 // {
 //   
 //   arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
 //   arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
 //   arma::vec p3vec = X3 * Rcpp::as<arma::vec>(pars[2]);
 //   arma::vec p4vec = X4 * Rcpp::as<arma::vec>(pars[3]);
 //   arma::vec p5vec = X5 * Rcpp::as<arma::vec>(pars[4]);
 //   arma::vec p6vec = X6 * Rcpp::as<arma::vec>(pars[5]);
 //   
 //   int nobs = ymat.n_rows;
 //   int ncol = ymat.n_cols;
 //   
 //   if (dcate == 1) {
 //     p1vec = p1vec.elem(dupid);
 //     p2vec = p2vec.elem(dupid);
 //     p3vec = p3vec.elem(dupid);
 //     p4vec = p4vec.elem(dupid);
 //     p5vec = p5vec.elem(dupid);
 //     p6vec = p6vec.elem(dupid);
 //   }
 //   
 //   double y, w, yi, z, p1, p2, p3, p4, p5, p6;
 //   double alpha, beta, nu, kappa1, kappa2, delta;
 //   double nllh = 0.0;
 //   
 //   for (int j=0; j < nobs; j++) {
 //     
 //     p1 = p1vec[j];
 //     p2 = p2vec[j];
 //     p3 = p3vec[j];
 //     p4 = p4vec[j];
 //     p5 = p5vec[j];
 //     p6 = p6vec[j];
 //     
 //     for (int l=0; l < ncol; l++) {
 //       
 //       y = ymat(j, l);
 //       yi = xmat(j, l);
 //       
 //       if (arma::is_finite(y) & arma::is_finite(yi)) {
 //         
 //         w = wmat(j, l);
 //         
 //         alpha = 2.0 / (1.0 + exp(-p1)) - 1.0;
 //         beta = 1.0 / (1.0 + exp(-p2));
 //         nu = p3;
 //         kappa1 = exp(p4);
 //         kappa2 = exp(p5);
 //         delta = exp(p6);
 //         z = (y - alpha * yi) / R_pow(yi, beta);
 //         nllh += w * log(kappa1 + kappa2) + lgamma(1 / delta) - log(delta);
 //         if (z < nu) {
 //           nllh += w * R_pow((nu - z) / kappa1, delta);
 //         } else {
 //           nllh += w * R_pow((z - nu) / kappa2, delta);
 //         }
 //         
 //       }
 //     }
 //   }
 //   
 //   return(nllh);
 //   
 // }
 // 
 // //' @rdname condexaggd0
 // // [[Rcpp::export]]
 // arma::mat condexaggd12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, 
 //                        arma::mat X4, arma::mat X5, arma::mat X6, 
 //                        arma::mat ymat, arma::mat xmat, arma::mat wmat, 
 //                        arma::uvec dupid, int dcate, arma::uvec nhere)
 // {
 //   
 //   arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
 //   arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
 //   arma::vec p3vec = X3 * Rcpp::as<arma::vec>(pars[2]);
 //   arma::vec p4vec = X4 * Rcpp::as<arma::vec>(pars[3]);
 //   arma::vec p5vec = X5 * Rcpp::as<arma::vec>(pars[4]);
 //   arma::vec p6vec = X6 * Rcpp::as<arma::vec>(pars[5]);
 //   
 //   int nobs = ymat.n_rows;
 //   int ncol = ymat.n_cols;
 //   
 //   if (dcate == 1) {
 //     p1vec = p1vec.elem(dupid);
 //     p2vec = p2vec.elem(dupid);
 //     p3vec = p3vec.elem(dupid);
 //     p4vec = p4vec.elem(dupid);
 //     p5vec = p5vec.elem(dupid);
 //     p6vec = p6vec.elem(dupid);
 //   }
 //   
 //   double y, w, yi, z, p1, p2, p3, p4, p5, p6;
 //   double alpha, beta, nu, kappa1, kappa2, delta;
 //   
 //   arma::mat out = arma::mat(nobs, 27, arma::fill::zeros);
 //   
 //   // double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
 //   // double ee12, ee13, ee15, ee16, ee17, ee18, ee19;
 //   // double ee20, ee21, ee22, ee23, ee24, ee25, ee27, ee28, ee29;
 //   // double ee30, ee33, ee34, ee35, ee36, ee37, ee38, ee39;
 //   // double ee40, ee41, ee42, ee43, ee44, ee45, ee46, ee47, ee48, ee49, ee50;
 //   
 //   for (int j=0; j < nobs; j++) {
 //     
 //     p1 = p1vec[j];
 //     p2 = p2vec[j];
 //     p3 = p3vec[j];
 //     p4 = p4vec[j];
 //     p5 = p5vec[j];
 //     p6 = p6vec[j];
 //     
 //     for (int l=0; l < ncol; l++) {
 //       
 //       y = ymat(j, l);
 //       yi = xmat(j, l);
 //       
 //       if (arma::is_finite(y) & arma::is_finite(yi)) {
 //         
 //         w = wmat(j, l);
 //         
 //         alpha = 2.0 / (1.0 + exp(-p1)) - 1.0;
 //         beta = 1.0 / (1.0 + exp(-p2));
 //         nu = p3;
 //         kappa1 = exp(p4);
 //         kappa2 = exp(p5);
 //         delta = exp(p6);
 //         z = (y - alpha * yi) / R_pow(yi, beta);
 //         
 //         if (z < nu) {
 //           
 //           // ee2 = exp(-p2);
 //           // ee3 = exp(p6);
 //           // ee5 = 1 + ee2;
 //           // ee6 = exp(-p1);
 //           // ee7 = exp(p4);
 //           // ee8 = 1 + ee6;
 //           // ee9 = 1/ee5;
 //           // ee12 = y - yi * (2/ee8 - 1);
 //           // ee13 = R_pow(yi, ee9);
 //           // ee15 = p3 - ee12/ee13;
 //           // ee16 = ee15/ee7;
 //           // ee17 = ee3 - 1;
 //           // ee18 = R_pow(ee16, ee17);
 //           // ee19 = exp(p5);
 //           // ee20 = R_pow(ee16, (ee3 - 2));
 //           // ee21 = log(yi);
 //           // ee22 = R_pow(ee8, 2);
 //           // ee23 = R_pow(ee5, 2);
 //           // ee24 = 1 - ee9;
 //           // ee25 = ee7 + ee19;
 //           // ee27 = log(ee15) - p4;
 //           // ee28 = R_pow(yi, ee24);
 //           // ee29 = ee22 * ee7;
 //           // ee30 = ee3 * ee27;
 //           // ee33 = ee20 * ee17 * ee15/ee7;
 //           // ee34 = 1 + ee30;
 //           // ee35 = ee18 + ee33;
 //           // ee36 = ee23 * ee7;
 //           // ee37 = 1/ee3;
 //           // ee38 = R_pow(ee7, 2);
 //           // ee39 = ee28 * ee18;
 //           // ee40 = ee35 * ee3;
 //           // ee42 = ee18 * ee34 * ee3;
 //           // ee43 = ee18 * ee3;
 //           // ee44 = ee20 * ee2;
 //           // ee45 = R_pow(ee16, ee3);
 //           // ee46 = 2/ee5;
 //           // ee47 = R::digamma(ee37);
 //           // ee48 = ee7/ee25;
 //           // ee49 = ee19/ee25;
 //           // ee50 = ee13 * ee23;
 //           // 
 //           // out(j, 0) += w * (2 * (ee39 * ee6 * ee3/ee29));
 //           // out(j, 1) += w * (ee18 * ee2 * ee3 * ee21 * ee12/(ee50 * ee7));
 //           // out(j, 2) += w * (ee43/ee7);
 //           // out(j, 3) += w * (ee48 - ee43 * ee15/ee7);
 //           // out(j, 4) += w * (ee49);
 //           // out(j, 5) += w * (ee45 * ee3 * ee27 - (1 + ee47/ee3));
 //           // out(j, 6) += w * ((4 * (R_pow(yi, (2 * ee24)) * ee20 * ee6 * ee17/
 //           //   ee29) - ee39 * (2 - 4 * (ee6/ee8))) * ee6 * ee3/ee29);
 //           // out(j, 7) += w * ((2 * (R_pow(yi, (1 - ee46)) * ee20 * ee17 * ee12/
 //           //   ee7) - 2 * ee39) * ee6 * ee2 * ee3 * ee21/(ee22 * ee23 * ee7));
 //           // out(j, 8) += w * (2 * (ee28 * ee20 * ee6 * ee3 * ee17/(ee22 * ee38)));
 //           // out(j, 9) += w * (-(ee6 * ee3 * (ee28 * (2 * ee33 + 2 * ee18))/ee29));
 //           // out(j, 10) += w * (0);
 //           // out(j, 11) += w * (ee18 * ee6 * ee3 * (ee28 * (2 + 2 * ee30))/ee29);
 //           // out(j, 12) += w * ((((2 - ee21/ee5) * ee2/ee5 - 1) * ee18/ee13 +
 //           //   ee44 * ee17 * ee21 * ee12/(R_pow(yi, ee46) * ee23 * ee7)) *
 //           //   ee2 * ee3 * ee21 * ee12/ee36);
 //           // out(j, 13) += w * (ee44 * ee3 * ee17 * ee21 * ee12/(ee50 * ee38));
 //           // out(j, 14) += w * (-(ee35/ee13 * ee2 * ee3 * ee21 * ee12/ee36));
 //           // out(j, 15) += w * (0);
 //           // out(j, 16) += w * (ee34/ee13 * ee18 * ee2 * ee3 * ee21 * ee12/ee36);
 //           // out(j, 17) += w * (ee20 * ee3 * ee17/ee38);
 //           // out(j, 18) += w * (-(ee40/ee7));
 //           // out(j, 19) += w * (0);
 //           // out(j, 20) += w * (ee42/ee7);
 //           // out(j, 21) += w * (ee40 * ee15/ee7 + (1 - ee48) * ee7/ee25);
 //           // out(j, 22) += w * (-(ee7 * ee19/R_pow(ee25, 2)));
 //           // out(j, 23) += w * (-(ee42 * ee15/ee7));
 //           // out(j, 24) += w * ((1 - ee49) * ee19/ee25);
 //           // out(j, 25) += w * (0);
 //           // out(j, 26) += w * (ee45 * ee34 * ee3 * ee27 + (ee47 + R::trigamma(ee37)/
 //           //   ee3)/ee3);
 //           
 //           double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
 //           double ee12, ee13, ee15, ee16, ee17, ee18, ee19;
 //           double ee20, ee21, ee22, ee23, ee24, ee25, ee27, ee28, ee29;
 //           double ee30, ee33, ee34, ee35, ee36, ee37, ee38, ee39;
 //           double ee40, ee42, ee43, ee44, ee45, ee46, ee47, ee48, ee49;
 //           double ee50;
 //           
 //           ee2 = exp(-p2);
 //           ee3 = exp(p6);
 //           ee5 = 1 + ee2;
 //           ee6 = exp(-p1);
 //           ee7 = exp(p4);
 //           ee8 = 1 + ee6;
 //           ee9 = 1/ee5;
 //           ee12 = y - yi * (2/ee8 - 1);
 //           ee13 = R_pow(yi, ee9);
 //           ee15 = p3 - ee12/ee13;
 //           ee16 = ee15/ee7;
 //           ee17 = ee3 - 1;
 //           ee18 = R_pow(ee16, ee17);
 //           ee19 = exp(p5);
 //           ee20 = R_pow(ee16, (ee3 - 2));
 //           ee21 = log(yi);
 //           ee22 = R_pow(ee8, 2);
 //           ee23 = R_pow(ee5, 2);
 //           ee24 = 1 - ee9;
 //           ee25 = ee7 + ee19;
 //           ee27 = log(ee15) - p4;
 //           ee28 = R_pow(yi, ee24);
 //           ee29 = ee22 * ee7;
 //           ee30 = ee3 * ee27;
 //           ee33 = ee20 * ee17 * ee15/ee7;
 //           ee34 = 1 + ee30;
 //           ee35 = ee18 + ee33;
 //           ee36 = ee23 * ee7;
 //           ee37 = 1/ee3;
 //           ee38 = R_pow(ee7, 2);
 //           ee39 = ee28 * ee18;
 //           ee40 = ee35 * ee3;
 //           ee42 = ee18 * ee34 * ee3;
 //           ee43 = ee18 * ee3;
 //           ee44 = ee20 * ee2;
 //           ee45 = R_pow(ee16, ee3);
 //           ee46 = 2/ee5;
 //           ee47 = R::digamma(ee37);
 //           ee48 = ee7/ee25;
 //           ee49 = ee19/ee25;
 //           ee50 = ee13 * ee23;
 //           
 //           out(j, 0) += w * (2 * (ee39 * ee6 * ee3/ee29));
 //           out(j, 1) += w * (ee18 * ee2 * ee3 * ee21 * ee12/(ee50 * ee7));
 //           out(j, 2) += w * (ee43/ee7);
 //           out(j, 3) += w * (ee48 - ee43 * ee15/ee7);
 //           out(j, 4) += w * (ee49);
 //           out(j, 5) += w * (ee45 * ee3 * ee27 - (1 + ee47/ee3));
 //           out(j, 6) += w * ((4 * (R_pow(yi, (2 * ee24)) * ee20 * ee6 *
 //             ee17/ee29) - ee39 * (2 - 4 * (ee6/ee8))) * ee6 * ee3/ee29);
 //           out(j, 7) += w * ((2 * (R_pow(yi, (1 - ee46)) * ee20 * ee17 *
 //             ee12/ee7) - 2 * ee39) * ee6 * ee2 * ee3 * ee21/(ee22 * ee23 *
 //             ee7));
 //           out(j, 8) += w * (2 * (ee28 * ee20 * ee6 * ee3 * ee17/(ee22 *
 //             ee38)));
 //           out(j, 9) += w * (-(ee6 * ee3 * (ee28 * (2 * ee33 + 2 * ee18))/
 //             ee29));
 //           out(j, 10) += w * (0);
 //           out(j, 11) += w * (ee18 * ee6 * ee3 * (ee28 * (2 + 2 * ee30))/
 //             ee29);
 //           out(j, 12) += w * ((((2 - ee21/ee5) * ee2/ee5 - 1) * ee18/ee13 +
 //             ee44 * ee17 * ee21 * ee12/(R_pow(yi, ee46) * ee23 * ee7)) *
 //             ee2 * ee3 * ee21 * ee12/ee36);
 //           out(j, 13) += w * (ee44 * ee3 * ee17 * ee21 * ee12/(ee50 * ee38));
 //           out(j, 14) += w * (-(ee35/ee13 * ee2 * ee3 * ee21 * ee12/ee36));
 //           out(j, 15) += w * (0);
 //           out(j, 16) += w * (ee34/ee13 * ee18 * ee2 * ee3 * ee21 * ee12/
 //             ee36);
 //           out(j, 17) += w * (ee20 * ee3 * ee17/ee38);
 //           out(j, 18) += w * (-(ee40/ee7));
 //           out(j, 19) += w * (0);
 //           out(j, 20) += w * (ee42/ee7);
 //           out(j, 21) += w * (ee40 * ee15/ee7 + (1 - ee48) * ee7/ee25);
 //           out(j, 22) += w * (-(ee7 * ee19/R_pow(ee25, 2)));
 //           out(j, 23) += w * (-(ee42 * ee15/ee7));
 //           out(j, 24) += w * ((1 - ee49) * ee19/ee25);
 //           out(j, 25) += w * (0);
 //           out(j, 26) += w * (ee45 * ee34 * ee3 * ee27 + (ee47 + R::trigamma(ee37)/
 //             ee3)/ee3);
 //           
 //         } else {
 //           
 //           // ee2 = exp(-p2);
 //           // ee3 = exp(p6);
 //           // ee5 = 1 + ee2;
 //           // ee6 = exp(-p1);
 //           // ee7 = exp(p5);
 //           // ee8 = 1 + ee6;
 //           // ee9 = 1/ee5;
 //           // ee12 = y - yi * (2/ee8 - 1);
 //           // ee13 = R_pow(yi, ee9);
 //           // ee15 = ee12/ee13 - p3;
 //           // ee16 = ee15/ee7;
 //           // ee17 = ee3 - 1;
 //           // ee18 = R_pow(ee16, ee17);
 //           // ee19 = exp(p4);
 //           // ee20 = R_pow(ee16, (ee3 - 2));
 //           // ee21 = log(yi);
 //           // ee22 = R_pow(ee8, 2);
 //           // ee23 = R_pow(ee5, 2);
 //           // ee24 = 1 - ee9;
 //           // ee25 = ee19 + ee7;
 //           // ee27 = log(ee15) - p5;
 //           // ee28 = R_pow(yi, ee24);
 //           // ee29 = ee22 * ee7;
 //           // ee30 = ee3 * ee27;
 //           // ee33 = ee20 * ee15 * ee17/ee7;
 //           // ee34 = 1 + ee30;
 //           // ee35 = ee18 + ee33;
 //           // ee36 = ee23 * ee7;
 //           // ee37 = 1/ee3;
 //           // ee38 = R_pow(ee7, 2);
 //           // ee39 = ee28 * ee18;
 //           // ee40 = ee18 * ee15;
 //           // ee41 = ee20 * ee2;
 //           // ee42 = R_pow(ee16, ee3);
 //           // ee43 = 2/ee5;
 //           // ee44 = R::digamma(ee37);
 //           // ee45 = ee19/ee25;
 //           // ee46 = ee7/ee25;
 //           // ee47 = ee13 * ee23;
 //           // 
 //           // out(j, 0) += w * (-(2 * (ee39 * ee6 * ee3/ee29)));
 //           // out(j, 1) += w * (-(ee18 * ee2 * ee3 * ee21 * ee12/(ee47 * ee7)));
 //           // out(j, 2) += w * (-(ee18 * ee3/ee7));
 //           // out(j, 3) += w * (ee45);
 //           // out(j, 4) += w * (ee46 - ee40 * ee3/ee7);
 //           // out(j, 5) += w * (ee42 * ee3 * ee27 - (1 + ee44/ee3));
 //           // out(j, 6) += w * ((4 * (R_pow(yi, (2 * ee24)) * ee20 * ee6 * ee17/
 //           //   ee29) + ee39 * (2 - 4 * (ee6/ee8))) * ee6 * ee3/ee29);
 //           // out(j, 7) += w * ((2 * ee39 + 2 * (R_pow(yi, (1 - ee43)) * ee20 *
 //           //   ee17 * ee12/ee7)) * ee6 * ee2 * ee3 * ee21/(ee22 * ee23 * ee7));
 //           // out(j, 8) += w * (2 * (ee28 * ee20 * ee6 * ee3 * ee17/(ee22 * ee38)));
 //           // out(j, 9) += w * (0);
 //           // out(j, 10) += w * (ee6 * ee3 * (ee28 * (2 * ee33 + 2 * ee18))/ee29);
 //           // out(j, 11) += w * (-(ee18 * ee6 * ee3 * (ee28 * (2 + 2 * ee30))/
 //           //   ee29));
 //           // out(j, 12) += w * (-((ee18 * ((2 - ee21/ee5) * ee2/ee5 - 1)/ee13 -
 //           //   ee41 * ee17 * ee21 * ee12/(R_pow(yi, ee43) * ee23 * ee7)) *
 //           //   ee2 * ee3 * ee21 * ee12/ee36));
 //           // out(j, 13) += w * (ee41 * ee3 * ee17 * ee21 * ee12/(ee47 * ee38));
 //           // out(j, 14) += w * (0);
 //           // out(j, 15) += w * (ee35/ee13 * ee2 * ee3 * ee21 * ee12/ee36);
 //           // out(j, 16) += w * (-(ee18 * (ee34/ee13) * ee2 * ee3 * ee21 * ee12/
 //           //   ee36));
 //           // out(j, 17) += w * (ee20 * ee3 * ee17/ee38);
 //           // out(j, 18) += w * (0);
 //           // out(j, 19) += w * (ee35 * ee3/ee7);
 //           // out(j, 20) += w * (-(ee18 * ee34 * ee3/ee7));
 //           // out(j, 21) += w * ((1 - ee45) * ee19/ee25);
 //           // out(j, 22) += w * (-(ee19 * ee7/R_pow(ee25, 2)));
 //           // out(j, 23) += w * (0);
 //           // out(j, 24) += w * (ee35 * ee15 * ee3/ee7 + (1 - ee46) * ee7/ee25);
 //           // out(j, 25) += w * (-(ee40 * ee34 * ee3/ee7));
 //           // out(j, 26) += w * (ee42 * ee34 * ee3 * ee27 + (ee44 + R::trigamma(ee37)/
 //           //   ee3)/ee3);
 //           
 //           double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
 //           double ee12, ee13, ee15, ee16, ee17, ee18, ee19;
 //           double ee20, ee21, ee22, ee23, ee24, ee25, ee27, ee28, ee29;
 //           double ee30, ee33, ee34, ee35, ee36, ee37, ee38, ee39;
 //           double ee40, ee41, ee42, ee43, ee44, ee45, ee46, ee47;
 //           
 //           ee2 = exp(-p2);
 //           ee3 = exp(p6);
 //           ee5 = 1 + ee2;
 //           ee6 = exp(-p1);
 //           ee7 = exp(p5);
 //           ee8 = 1 + ee6;
 //           ee9 = 1/ee5;
 //           ee12 = y - yi * (2/ee8 - 1);
 //           ee13 = R_pow(yi, ee9);
 //           ee15 = ee12/ee13 - p3;
 //           ee16 = ee15/ee7;
 //           ee17 = ee3 - 1;
 //           ee18 = R_pow(ee16, ee17);
 //           ee19 = exp(p4);
 //           ee20 = R_pow(ee16, (ee3 - 2));
 //           ee21 = log(yi);
 //           ee22 = R_pow(ee8, 2);
 //           ee23 = R_pow(ee5, 2);
 //           ee24 = 1 - ee9;
 //           ee25 = ee19 + ee7;
 //           ee27 = log(ee15) - p5;
 //           ee28 = R_pow(yi, ee24);
 //           ee29 = ee22 * ee7;
 //           ee30 = ee3 * ee27;
 //           ee33 = ee20 * ee15 * ee17/ee7;
 //           ee34 = 1 + ee30;
 //           ee35 = ee18 + ee33;
 //           ee36 = ee23 * ee7;
 //           ee37 = 1/ee3;
 //           ee38 = R_pow(ee7, 2);
 //           ee39 = ee28 * ee18;
 //           ee40 = ee18 * ee15;
 //           ee41 = ee20 * ee2;
 //           ee42 = R_pow(ee16, ee3);
 //           ee43 = 2/ee5;
 //           ee44 = R::digamma(ee37);
 //           ee45 = ee19/ee25;
 //           ee46 = ee7/ee25;
 //           ee47 = ee13 * ee23;
 //           
 //           out(j, 0) += w * (-(2 * (ee39 * ee6 * ee3/ee29)));
 //           out(j, 1) += w * (-(ee18 * ee2 * ee3 * ee21 * ee12/(ee47 * ee7)));
 //           out(j, 2) += w * (-(ee18 * ee3/ee7));
 //           out(j, 3) += w * (ee45);
 //           out(j, 4) += w * (ee46 - ee40 * ee3/ee7);
 //           out(j, 5) += w * (ee42 * ee3 * ee27 - (1 + ee44/ee3));
 //           out(j, 6) += w * ((4 * (R_pow(yi, (2 * ee24)) * ee20 * ee6 *
 //             ee17/ee29) + ee39 * (2 - 4 * (ee6/ee8))) * ee6 * ee3/ee29);
 //           out(j, 7) += w * ((2 * ee39 + 2 * (R_pow(yi, (1 - ee43)) * ee20 *
 //             ee17 * ee12/ee7)) * ee6 * ee2 * ee3 * ee21/(ee22 * ee23 *
 //             ee7));
 //           out(j, 8) += w * (2 * (ee28 * ee20 * ee6 * ee3 * ee17/(ee22 *
 //             ee38)));
 //           out(j, 9) += w * (0);
 //           out(j, 10) += w * (ee6 * ee3 * (ee28 * (2 * ee33 + 2 * ee18))/
 //             ee29);
 //           out(j, 11) += w * (-(ee18 * ee6 * ee3 * (ee28 * (2 + 2 * ee30))/
 //             ee29));
 //           out(j, 12) += w * (-((ee18 * ((2 - ee21/ee5) * ee2/ee5 - 1)/
 //             ee13 - ee41 * ee17 * ee21 * ee12/(R_pow(yi, ee43) * ee23 * ee7)) *
 //               ee2 * ee3 * ee21 * ee12/ee36));
 //           out(j, 13) += w * (ee41 * ee3 * ee17 * ee21 * ee12/(ee47 * ee38));
 //           out(j, 14) += w * (0);
 //           out(j, 15) += w * (ee35/ee13 * ee2 * ee3 * ee21 * ee12/ee36);
 //           out(j, 16) += w * (-(ee18 * (ee34/ee13) * ee2 * ee3 * ee21 *
 //             ee12/ee36));
 //           out(j, 17) += w * (ee20 * ee3 * ee17/ee38);
 //           out(j, 18) += w * (0);
 //           out(j, 19) += w * (ee35 * ee3/ee7);
 //           out(j, 20) += w * (-(ee18 * ee34 * ee3/ee7));
 //           out(j, 21) += w * ((1 - ee45) * ee19/ee25);
 //           out(j, 22) += w * (-(ee19 * ee7/R_pow(ee25, 2)));
 //           out(j, 23) += w * (0);
 //           out(j, 24) += w * (ee35 * ee15 * ee3/ee7 + (1 - ee46) * ee7/
 //             ee25);
 //           out(j, 25) += w * (-(ee40 * ee34 * ee3/ee7));
 //           out(j, 26) += w * (ee42 * ee34 * ee3 * ee27 + (ee44 + R::trigamma(ee37)/
 //             ee3)/ee3);
 //           
 //         }
 //         
 //       }}}    
 //   
 //   return out;
 //   
 // }
 // 
 //  //' @rdname condexaggd0
 // // [[Rcpp::export]]
 // arma::mat condexaggd34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3,
 //                        arma::mat X4, arma::mat X5, arma::mat X6, 
 //                        arma::mat ymat, arma::mat xmat, arma::mat wmat, 
 //                        arma::uvec dupid, int dcate, arma::uvec nhere)
 // {
 //   
 //   arma::vec p1vec = X1 * Rcpp::as<arma::vec>(pars[0]);
 //   arma::vec p2vec = X2 * Rcpp::as<arma::vec>(pars[1]);
 //   arma::vec p3vec = X3 * Rcpp::as<arma::vec>(pars[2]);
 //   arma::vec p4vec = X4 * Rcpp::as<arma::vec>(pars[3]);
 //   arma::vec p5vec = X5 * Rcpp::as<arma::vec>(pars[4]);
 //   arma::vec p6vec = X6 * Rcpp::as<arma::vec>(pars[5]);
 //   
 //   int nobs = ymat.n_rows;
 //   int ncol = ymat.n_cols;
 //   
 //   if (dcate == 1) {
 //     p1vec = p1vec.elem(dupid);
 //     p2vec = p2vec.elem(dupid);
 //     p3vec = p3vec.elem(dupid);
 //     p4vec = p4vec.elem(dupid);
 //     p5vec = p5vec.elem(dupid);
 //     p6vec = p6vec.elem(dupid);
 //   }
 //   
 //   double y, w, yi, z, p1, p2, p3, p4, p5, p6;
 //   double alpha, beta, nu, kappa1, kappa2, delta;
 //   
 //   arma::mat out = arma::mat(nobs, 56, arma::fill::zeros);
 //   
 //   double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
 //   double ee12, ee13, ee15, ee16, ee17, ee18, ee19;
 //   double ee20, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
 //   double ee30, ee31, ee32, ee33, ee34, ee35, ee36, ee37, ee39;
 //   double ee42, ee44, ee45, ee46, ee47;
 //   double ee50, ee51, ee52, ee53, ee54, ee55, ee56, ee57, ee58;
 //   double ee61, ee62, ee63, ee66, ee67;
 //   double ee70, ee72, ee74, ee77, ee78, ee79;
 //   double ee80, ee82, ee83, ee84, ee85, ee86, ee87, ee88;
 //   double ee91, ee92, ee94, ee95, ee96, ee97, ee98;
 //   double ee101, ee102, ee103, ee105, ee108, ee109;
 //   double ee110, ee111;
 //   
 //   double ee38, ee40, ee43, ee48, ee49, ee59, ee60, ee64, ee69;
 //   double ee76, ee81, ee89, ee90, ee93, ee100, ee104;
 //   double ee112, ee113, ee114, ee115;
 //   
 //   for (int j=0; j < nobs; j++) {
 //     
 //     p1 = p1vec[j];
 //     p2 = p2vec[j];
 //     p3 = p3vec[j];
 //     p4 = p4vec[j];
 //     p5 = p5vec[j];
 //     p6 = p6vec[j];
 //     
 //     for (int l=0; l < ncol; l++) {
 //       
 //       y = ymat(j, l);
 //       yi = xmat(j, l);
 //       
 //       if (arma::is_finite(y) & arma::is_finite(yi)) {
 //         
 //         w = wmat(j, l);
 //         
 //         alpha = 2.0 / (1.0 + exp(-p1)) - 1.0;
 //         beta = 1.0 / (1.0 + exp(-p2));
 //         nu = p3;
 //         kappa1 = exp(p4);
 //         kappa2 = exp(p5);
 //         delta = exp(p6);
 //         z = (y - alpha * yi) / R_pow(yi, beta);
 //         
 //         if (z < nu) {
 //           
 //           ee2 = exp(-p2);
 //           ee3 = 1 + ee2;
 //           ee5 = exp(-p1);
 //           ee6 = exp(p6);
 //           ee7 = 1 + ee5;
 //           ee8 = 1/ee3;
 //           ee9 = exp(p4);
 //           ee12 = y - yi * (2/ee7 - 1);
 //           ee13 = R_pow(yi, ee8);
 //           ee15 = p3 - ee12/ee13;
 //           ee16 = ee15/ee9;
 //           ee17 = ee6 - 1;
 //           ee18 = ee6 - 2;
 //           ee19 = R_pow(ee16, ee18);
 //           ee20 = log(yi);
 //           ee22 = log(ee15) - p4;
 //           ee23 = 1 - ee8;
 //           ee24 = R_pow(ee3, 2);
 //           ee25 = R_pow(ee7, 2);
 //           ee26 = R_pow(ee16, ee17);
 //           ee27 = R_pow(yi, ee23);
 //           ee28 = exp(p5);
 //           ee29 = R_pow(ee16, (ee6 - 3));
 //           ee30 = ee25 * ee9;
 //           ee31 = 2/ee3;
 //           ee32 = ee6 * ee22;
 //           ee33 = ee17 * ee22;
 //           ee34 = ee9 + ee28;
 //           ee35 = ee20/ee3;
 //           ee36 = R_pow(ee9, 2);
 //           ee37 = ee24 * ee9;
 //           ee38 = 2 * ee23;
 //           ee40 = (ee33 + 2) * ee6 - 1;
 //           ee43 = (2 - ee35) * ee2/ee3 - 1;
 //           ee45 = 1 + ee32;
 //           ee46 = 2 - 4 * (ee5/ee7);
 //           ee47 = 2 * ee17;
 //           ee48 = ee5 * ee6;
 //           ee49 = R_pow(yi, (1 - ee31));
 //           ee50 = R_pow(yi, ee38);
 //           ee51 = R_pow(yi, ee31);
 //           ee52 = ee40 * ee19;
 //           ee53 = ee25 * ee24;
 //           ee54 = 2 * ee19;
 //           ee55 = ee27 * ee26;
 //           ee58 = ee29 * ee18 * ee15/ee9;
 //           ee59 = ee29 * ee2;
 //           ee60 = ee53 * ee9;
 //           ee63 = (3 + ee32) * ee6 * ee22 + 1;
 //           ee64 = ee27 * ee19;
 //           ee66 = ee51 * ee24 * ee9;
 //           ee69 = ee52 * ee15/ee9 + ee26 * ee45;
 //           ee72 = ee43 * ee19/ee13 + ee59 * ee18 * ee20 * ee12/ee66;
 //           ee76 = ee19 * ee17 * ee15/ee9;
 //           ee78 = ee25 * ee36;
 //           ee79 = ee24 * ee36;
 //           ee81 = 1/ee6;
 //           ee83 = 2 * ((ee33 + 1) * ee6) + ee47;
 //           ee84 = 2 * (ee19 + ee58);
 //           ee85 = 2 * (1 + ee47);
 //           ee88 = ee5 * ee2 * ee6 * ee20;
 //           ee89 = R_pow(ee9, 3);
 //           ee90 = ee69 * ee6;
 //           ee92 = ee63 * ee26 * ee6;
 //           ee93 = ee26 + ee76;
 //           ee94 = ee19 * ee6;
 //           ee97 = (2 + ee35) * ee2/ee3 - 1;
 //           ee98 = R_pow(ee34, 2);
 //           ee100 = ee17/ee9 + ee85;
 //           ee101 = ee84 + ee54;
 //           ee102 = 2 * ee26;
 //           ee103 = 2 * (ee9/ee34);
 //           ee104 = 2 * (ee28/ee34);
 //           ee108 = 2 * (ee49 * ee29 * ee18 * ee12/ee9) - 2 * ee64;
 //           ee109 = 2 + 2 * ee32;
 //           ee110 = 4 * ee19;
 //           ee112 = 4 * (ee50 * ee29 * ee5 * ee18/ee30) - ee64 * ee46;
 //           ee113 = ee49 * ee19;
 //           ee114 = ee13 * ee72;
 //           ee115 = ee50 * ee19;
 //           
 //           out(j, 0) += w * (ee48 * (ee55 * (2 - ((4 * (1 + 2 * ee5) + 4 * ee7 -
 //             16 * ee5)/ee7 + 4) * ee5/ee7) + ee50 * (2 * (R_pow(yi, (1 -
 //             (ee8 + ee38))) * ee112) - 4 * (ee19 * ee46)) * ee5 * ee17/
 //               ee30)/ee30);
 //           out(j, 1) += w * (((2 * (ee27 * ee108) - 4 * ee115) * ee5 * ee17/
 //             ee30 - ee46 * (ee113 * ee17 * ee12/ee9 - ee55)) * ee5 * ee2 *
 //               ee6 * ee20/ee60);
 //           out(j, 2) += w * (ee112 * ee5 * ee6 * ee17/ee78);
 //           out(j, 3) += w * (-(ee48 * (ee50 * (2 * (2 * ee58 + ee54) + ee110) *
 //             ee5 * ee17/ee30 - ee27 * ee93 * ee46)/ee30));
 //           out(j, 4) += w * (0);
 //           out(j, 5) += w * (ee48 * (ee115 * (2 * ((2 + 2 * ee33) * ee6) + 4 *
 //             ee17) * ee5/ee30 - ee55 * ee45 * ee46)/ee30);
 //           out(j, 6) += w * (ee88 * (ee49 * (2 * ee114 - 4 * (ee19 * ee2 * ee20/
 //             ee24)) * ee17 * ee12/ee9 - ee55 * (2 * ee97 - 4 * (ee2 *
 //               ee20/ee24)))/ee60);
 //           out(j, 7) += w * (ee108 * ee5 * ee2 * ee6 * ee17 * ee20/(ee53 * ee36));
 //           out(j, 8) += w * (-(ee88 * (ee49 * ee101 * ee17 * ee12/ee9 - ee27 *
 //             (2 * ee76 + ee102))/ee60));
 //           out(j, 9) += w * (0);
 //           out(j, 10) += w * (ee88 * (ee113 * ee83 * ee12/ee9 - ee55 * ee109)/
 //             ee60);
 //           out(j, 11) += w * (2 * (ee27 * ee29 * ee5 * ee6 * ee17 * ee18/(ee25 *
 //             ee89)));
 //           out(j, 12) += w * (-(ee48 * ee17 * (ee27 * ee101)/ee78));
 //           out(j, 13) += w * (0);
 //           out(j, 14) += w * (ee19 * ee5 * ee6 * (ee27 * ee83)/ee78);
 //           out(j, 15) += w * (ee48 * (ee27 * ((ee84 + ee110) * ee17 * ee15/
 //             ee9 + ee102))/ee30);
 //           out(j, 16) += w * (0);
 //           out(j, 17) += w * (-(ee48 * (ee27 * (ee26 * ee109 + ee19 * ee83 *
 //             ee15/ee9))/ee30));
 //           out(j, 18) += w * (0);
 //           out(j, 19) += w * (0);
 //           out(j, 20) += w * (ee26 * ee5 * ee6 * (ee27 * ((2 * ee45 + 4) * ee6 *
 //             ee22 + 2))/ee30);
 //           out(j, 21) += w * (((((((6 + 9 * ee35) * ee2/ee3 - (2 * ((2 + 2 *
 //             ee35) * ee2/ee3 - 1) + 3 + 4 * ee97)) * ee20 + 8 * ee2 - (2 *
 //             (1 + 2 * ee2) + 2 * ee3))/ee3 - 2) * ee2/ee3 + 1) * ee26/
 //               ee13 + (ee43 * ee54 + ee114)/ee51 * ee2 * ee17 * ee20 * ee12/
 //                 ee37) * ee2 * ee6 * ee20 * ee12/ee37);
 //           out(j, 22) += w * (ee72 * ee2 * ee6 * ee17 * ee20 * ee12/ee79);
 //           out(j, 23) += w * (-((ee54 * ee2 * ee17 * ee20 * ee12/ee37 + ee13 *
 //             ee93 * ee43)/ee51 * ee2 * ee20 * ee12/ee37 * ee6));
 //           out(j, 24) += w * (0);
 //           out(j, 25) += w * ((ee52 * ee2 * ee20 * ee12/ee66 + ee43 * ee26 *
 //             ee45/ee13) * ee2 * ee6 * ee20 * ee12/ee37);
 //           out(j, 26) += w * (ee59 * ee6 * ee17 * ee18 * ee20 * ee12/(ee13 *
 //             ee24 * ee89));
 //           out(j, 27) += w * (-(ee94/ee13 * ee2 * ee6 * ee17 * ee20 * ee12/
 //             ee79));
 //           out(j, 28) += w * (0);
 //           out(j, 29) += w * (ee40/ee13 * ee19 * ee2 * ee6 * ee20 * ee12/ee79);
 //           out(j, 30) += w * (ee100/ee13 * ee19 * ee2 * ee20 * ee15 * ee12/
 //             ee24/ee9 * ee6);
 //           out(j, 31) += w * (0);
 //           out(j, 32) += w * (-(ee69/ee13 * ee2 * ee6 * ee20 * ee12/ee37));
 //           out(j, 33) += w * (0);
 //           out(j, 34) += w * (0);
 //           out(j, 35) += w * (ee63/ee13 * ee26 * ee2 * ee6 * ee20 * ee12/ee37);
 //           out(j, 36) += w * (ee29 * ee6 * ee17 * ee18/ee89);
 //           out(j, 37) += w * (-(ee94 * ee6 * ee17/ee36));
 //           out(j, 38) += w * (0);
 //           out(j, 39) += w * (ee52 * ee6/ee36);
 //           out(j, 40) += w * (ee100 * ee19 * ee15/ee9 * ee6);
 //           out(j, 41) += w * (0);
 //           out(j, 42) += w * (-(ee90/ee9));
 //           out(j, 43) += w * (0);
 //           out(j, 44) += w * (0);
 //           out(j, 45) += w * (ee92/ee9);
 //           out(j, 46) += w * ((1 - (3 - ee103) * ee9/ee34) * ee9/ee34 - (ee17 *
 //             ee15/ee9 + ee85) * ee19 * ee6 * ee15/ee9);
 //           out(j, 47) += w * (-((1 - ee103) * ee9 * ee28/ee98));
 //           out(j, 48) += w * (ee90 * ee15/ee9);
 //           out(j, 49) += w * (-((1 - ee104) * ee9 * ee28/ee98));
 //           out(j, 50) += w * (0);
 //           out(j, 51) += w * (-(ee92 * ee15/ee9));
 //           out(j, 52) += w * ((1 - (3 - ee104) * ee28/ee34) * ee28/ee34);
 //           out(j, 53) += w * (0);
 //           out(j, 54) += w * (0);
 //           out(j, 55) += w * (ee63 * R_pow(ee16, ee6) * ee6 * ee22 - ((3 * R::trigamma(ee81) +
 //             R::psigamma(ee81, 2)/ee6)/ee6 + R::digamma(ee81))/ee6);
 //           
 //        } else {
 //           
 //           ee2 = exp(-p2);
 //           ee3 = 1 + ee2;
 //           ee5 = exp(-p1);
 //           ee6 = exp(p6);
 //           ee7 = 1 + ee5;
 //           ee8 = 1/ee3;
 //           ee9 = exp(p5);
 //           ee12 = y - yi * (2/ee7 - 1);
 //           ee13 = R_pow(yi, ee8);
 //           ee15 = ee12/ee13 - p3;
 //           ee16 = ee15/ee9;
 //           ee17 = ee6 - 1;
 //           ee18 = ee6 - 2;
 //           ee19 = R_pow(ee16, ee18);
 //           ee20 = log(yi);
 //           ee22 = log(ee15) - p5;
 //           ee23 = 1 - ee8;
 //           ee24 = R_pow(ee3, 2);
 //           ee25 = R_pow(ee7, 2);
 //           ee26 = R_pow(ee16, ee17);
 //           ee27 = R_pow(yi, ee23);
 //           ee28 = exp(p4);
 //           ee29 = R_pow(ee16, (ee6 - 3));
 //           ee30 = ee25 * ee9;
 //           ee31 = 2/ee3;
 //           ee32 = ee6 * ee22;
 //           ee33 = ee17 * ee22;
 //           ee34 = ee28 + ee9;
 //           ee35 = ee20/ee3;
 //           ee36 = R_pow(ee9, 2);
 //           ee37 = 2 * ee23;
 //           ee39 = (ee33 + 2) * ee6 - 1;
 //           ee42 = (2 - ee35) * ee2/ee3 - 1;
 //           ee44 = 1 + ee32;
 //           ee45 = 2 - 4 * (ee5/ee7);
 //           ee46 = R_pow(yi, (1 - ee31));
 //           ee47 = R_pow(yi, ee37);
 //           ee50 = R_pow(yi, ee31) * ee24 * ee9;
 //           ee51 = ee39 * ee19;
 //           ee52 = ee25 * ee24;
 //           ee53 = ee24 * ee9;
 //           ee54 = 2 * ee17;
 //           ee55 = ee5 * ee6;
 //           ee56 = ee27 * ee26;
 //           ee57 = ee19 * ee42;
 //           ee58 = ee19 * ee15;
 //           ee61 = ee29 * ee15 * ee18/ee9;
 //           ee62 = ee29 * ee2;
 //           ee63 = ee52 * ee9;
 //           ee66 = (3 + ee32) * ee6 * ee22 + 1;
 //           ee67 = ee27 * ee19;
 //           ee70 = ee51 * ee15/ee9 + ee26 * ee44;
 //           ee72 = ee57/ee13 - ee62 * ee18 * ee20 * ee12/ee50;
 //           ee74 = ee58 * ee17/ee9;
 //           ee77 = ee25 * ee36;
 //           ee78 = ee24 * ee36;
 //           ee79 = 1/ee6;
 //           ee80 = 2 * (ee19 + ee61);
 //           ee82 = 2 * ((ee33 + 1) * ee6) + ee54;
 //           ee83 = 2 * ee19;
 //           ee84 = R_pow(ee9, 3);
 //           ee85 = ee26 * ee66;
 //           ee86 = ee26 + ee74;
 //           ee87 = ee19 * ee2;
 //           ee88 = ee19 * ee6;
 //           ee91 = (2 + ee35) * ee2/ee3 - 1;
 //           ee92 = R_pow(ee34, 2);
 //           ee94 = ee17/ee9 + 2 * (1 + ee54);
 //           ee95 = ee80 + ee83;
 //           ee96 = 2 * ee26;
 //           ee97 = 2 * (ee28/ee34);
 //           ee98 = 2 * (ee9/ee34);
 //           ee101 = 2 * ee67 + 2 * (ee46 * ee29 * ee18 * ee12/ee9);
 //           ee102 = 2 + 2 * ee32;
 //           ee103 = 4 * ee19;
 //           ee105 = 4 * (ee47 * ee29 * ee5 * ee18/ee30) + ee67 * ee45;
 //           ee108 = ee5 * ee2 * ee6 * ee20;
 //           ee109 = ee46 * ee19;
 //           ee110 = ee13 * ee72;
 //           ee111 = ee47 * ee19;
 //           
 //           out(j, 0) += w * (-(ee5 * (ee5 * ee17 * (ee47 * (2 * (R_pow(yi, (1 -
 //             (ee8 + ee37))) * ee105) + 4 * (ee19 * ee45)))/ee30 + ee56 *
 //             (2 - ((4 * (1 + 2 * ee5) + 4 * ee7 - 16 * ee5)/ee7 + 4) *
 //             ee5/ee7)) * ee6/ee30));
 //           out(j, 1) += w * (-((ee45 * (ee56 + ee109 * ee17 * ee12/ee9) + (2 *
 //             (ee27 * ee101) + 4 * ee111) * ee5 * ee17/ee30) * ee5 * ee2 *
 //             ee6 * ee20/ee63));
 //           out(j, 2) += w * (-(ee105 * ee5 * ee6 * ee17/ee77));
 //           out(j, 3) += w * (0);
 //           out(j, 4) += w * (-(ee55 * (ee27 * ee86 * ee45 + ee47 * (2 * (2 *
 //             ee61 + ee83) + ee103) * ee5 * ee17/ee30)/ee30));
 //           out(j, 5) += w * (ee55 * (ee56 * ee44 * ee45 + ee111 * (2 * ((2 +
 //             2 * ee33) * ee6) + 4 * ee17) * ee5/ee30)/ee30);
 //           out(j, 6) += w * ((ee17 * ee12 * (ee46 * (2 * ee110 - 4 * (ee87 *
 //             ee20/ee24)))/ee9 + ee56 * (2 * ee91 - 4 * (ee2 * ee20/ee24))) *
 //             ee5 * ee2 * ee6 * ee20/ee63);
 //           out(j, 7) += w * (-(ee101 * ee5 * ee2 * ee6 * ee17 * ee20/(ee52 *
 //             ee36)));
 //           out(j, 8) += w * (0);
 //           out(j, 9) += w * (-(ee108 * (ee27 * (2 * ee74 + ee96) + ee46 * ee95 *
 //             ee17 * ee12/ee9)/ee63));
 //           out(j, 10) += w * (ee108 * (ee56 * ee102 + ee109 * ee82 * ee12/ee9)/
 //             ee63);
 //           out(j, 11) += w * (-(2 * (ee27 * ee29 * ee5 * ee6 * ee17 * ee18/
 //             (ee25 * ee84))));
 //           out(j, 12) += w * (0);
 //           out(j, 13) += w * (-(ee55 * ee17 * (ee27 * ee95)/ee77));
 //           out(j, 14) += w * (ee19 * ee5 * ee6 * (ee27 * ee82)/ee77);
 //           out(j, 15) += w * (0);
 //           out(j, 16) += w * (0);
 //           out(j, 17) += w * (0);
 //           out(j, 18) += w * (-(ee55 * (ee27 * (ee15 * (ee80 + ee103) * ee17/
 //             ee9 + ee96))/ee30));
 //           out(j, 19) += w * (ee55 * (ee27 * (ee26 * ee102 + ee58 * ee82/ee9))/
 //             ee30);
 //           out(j, 20) += w * (-(ee26 * ee5 * ee6 * (ee27 * ((2 * ee44 + 4) *
 //             ee6 * ee22 + 2))/ee30));
 //           out(j, 21) += w * (-(((((((6 + 9 * ee35) * ee2/ee3 - (2 * ((2 + 2 *
 //             ee35) * ee2/ee3 - 1) + 3 + 4 * ee91)) * ee20 + 8 * ee2 -
 //             (2 * (1 + 2 * ee2) + 2 * ee3))/ee3 - 2) * ee2/ee3 + 1) * ee26/
 //               ee13 - (2 * ee57 + ee110) * ee2 * ee17 * ee20 * ee12/ee50) *
 //                 ee2 * ee6 * ee20 * ee12/ee53));
 //           out(j, 22) += w * (ee72 * ee2 * ee6 * ee17 * ee20 * ee12/ee78);
 //           out(j, 23) += w * (0);
 //           out(j, 24) += w * ((ee86 * ee42/ee13 - 2 * (ee87 * ee17 * ee20 *
 //             ee12/ee50)) * ee2 * ee20 * ee12/ee53 * ee6);
 //           out(j, 25) += w * (-((ee26 * ee42 * ee44/ee13 - ee51 * ee2 * ee20 *
 //             ee12/ee50) * ee2 * ee6 * ee20 * ee12/ee53));
 //           out(j, 26) += w * (-(ee62 * ee6 * ee17 * ee18 * ee20 * ee12/(ee13 *
 //             ee24 * ee84)));
 //           out(j, 27) += w * (0);
 //           out(j, 28) += w * (-(ee88/ee13 * ee2 * ee6 * ee17 * ee20 * ee12/
 //             ee78));
 //           out(j, 29) += w * (ee39/ee13 * ee19 * ee2 * ee6 * ee20 * ee12/ee78);
 //           out(j, 30) += w * (0);
 //           out(j, 31) += w * (0);
 //           out(j, 32) += w * (0);
 //           out(j, 33) += w * (-(ee94/ee13 * ee19 * ee15 * ee2 * ee20 * ee12/
 //             ee24/ee9 * ee6));
 //           out(j, 34) += w * (ee70/ee13 * ee2 * ee6 * ee20 * ee12/ee53);
 //           out(j, 35) += w * (-(ee66/ee13 * ee26 * ee2 * ee6 * ee20 * ee12/
 //             ee53));
 //           out(j, 36) += w * (-(ee29 * ee6 * ee17 * ee18/ee84));
 //           out(j, 37) += w * (0);
 //           out(j, 38) += w * (-(ee88 * ee6 * ee17/ee36));
 //           out(j, 39) += w * (ee51 * ee6/ee36);
 //           out(j, 40) += w * (0);
 //           out(j, 41) += w * (0);
 //           out(j, 42) += w * (0);
 //           out(j, 43) += w * (-(ee19 * ee94 * ee15/ee9 * ee6));
 //           out(j, 44) += w * (ee70 * ee6/ee9);
 //           out(j, 45) += w * (-(ee85 * ee6/ee9));
 //           out(j, 46) += w * ((1 - (3 - ee97) * ee28/ee34) * ee28/ee34);
 //           out(j, 47) += w * (-((1 - ee97) * ee28 * ee9/ee92));
 //           out(j, 48) += w * (0);
 //           out(j, 49) += w * (-((1 - ee98) * ee28 * ee9/ee92));
 //           out(j, 50) += w * (0);
 //           out(j, 51) += w * (0);
 //           out(j, 52) += w * ((1 - (3 - ee98) * ee9/ee34) * ee9/ee34 - ((ee16 +
 //             4) * ee17 + 2) * ee19 * ee15 * ee6/ee9);
 //           out(j, 53) += w * (ee70 * ee15 * ee6/ee9);
 //           out(j, 54) += w * (-(ee85 * ee15 * ee6/ee9));
 //           out(j, 55) += w * (R_pow(ee16, ee6) * ee66 * ee6 * ee22 - ((3 * R::trigamma(ee79) +
 //             R::psigamma(ee79, 2)/ee6)/ee6 + R::digamma(ee79))/ee6);
 //           
 //        }
 //           
 //       }}}    
 //   
 //   return out;
 //   
 // }
 // 
