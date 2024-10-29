// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double logtwopi = log(2.0 * M_PI);

//' Conditional extreme value model negative log-likelihood
 //'
 //' @param pars a list of vectors of coefficients for each conditional EVD parameter
 //' @param X1 a design matrix for (transformed) alpha
 //' @param X2 a design matrix for (transformed) beta
 //' @param X3 a design matrix for mu
 //' @param X4 a design matrix for (transformed) sigma
 //' @param ymat a matrix
 //' @param xmat a matrix
 //' @param wmat a matrix
 //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
 //' @return condexd0 a scalar, the negative log-likelihood
 //' @return condexd12 a matrix, first then second derivatives w.r.t. parameters
 //' @return condexd34 a matrix, third then fourth derivatives w.r.t. parameters
 //' @examples
 //' ## to follow
 //' @export
 // [[Rcpp::export]]
 double condexd0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, 
                 arma::mat ymat, arma::mat xmat, arma::mat wmat, 
                 arma::uvec dupid, int dcate, arma::uvec nhere)
 {
   
   arma::vec tavec = X1 * Rcpp::as<arma::vec>(pars[0]);
   arma::vec tbvec = X2 * Rcpp::as<arma::vec>(pars[1]);
   arma::vec muvec = X3 * Rcpp::as<arma::vec>(pars[2]);
   arma::vec tsigvec = X4 * Rcpp::as<arma::vec>(pars[3]);
   
   int nobs = ymat.n_rows;
   int ncol = ymat.n_cols;
   
   if (dcate == 1) {
     tavec = tavec.elem(dupid);
     tbvec = tbvec.elem(dupid);
     muvec = muvec.elem(dupid);
     tsigvec = tsigvec.elem(dupid);
   }
   
   double y, w, yi, talpha, tbeta, mu, tsigma;
   double alpha, beta, sigma, mu2, sigma2, yib, res;
   double nllh=0.0;
   
   for (int j=0; j < nobs; j++) {
     
     talpha = tavec[j];
     tbeta = tbvec[j];
     mu = muvec[j];
     tsigma = tsigvec[j];
     
     for (int l=0; l < ncol; l++) {
       
       y = ymat(j, l);
       yi = xmat(j, l);
       
       if (arma::is_finite(y) & arma::is_finite(yi)) {
         
         w = wmat(j, l);
         
         alpha = (2.0 / (1.0 + exp(-talpha)) - 1.0);
         beta = 1.0 - exp(tbeta);
         sigma = exp(tsigma);
         yib = R_pow(yi, beta);
         mu2 = alpha * yi + yib * mu;
         sigma2 = sigma * yib;
         res = (y - mu2) / sigma2;
         nllh += w * (0.5 * logtwopi + log(sigma2) + 0.5 * res * res); 
         
       }
     }
   }
   
   return(nllh);
   
 }
 
 //' @rdname condexd0
 // [[Rcpp::export]]
 arma::mat condexd12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, 
                     arma::mat ymat, arma::mat xmat, arma::mat wmat, 
                     arma::uvec dupid, int dcate, arma::uvec nhere)
 {
   
   arma::vec tavec = X1 * Rcpp::as<arma::vec>(pars[0]);
   arma::vec tbvec = X2 * Rcpp::as<arma::vec>(pars[1]);
   arma::vec muvec = X3 * Rcpp::as<arma::vec>(pars[2]);
   arma::vec tsigvec = X4 * Rcpp::as<arma::vec>(pars[3]);
   
   int nobs = ymat.n_rows;
   int ncol = ymat.n_cols;
   
   if (dcate == 1) {
     tavec = tavec.elem(dupid);
     tbvec = tbvec.elem(dupid);
     muvec = muvec.elem(dupid);
     tsigvec = tsigvec.elem(dupid);
   }
   
   double y, w, yi, talpha, tbeta, mu, tsigma;
   
   arma::mat out = arma::mat(nobs, 14, arma::fill::zeros);
   
   double ee1, ee2, ee3, ee4, ee6, ee7, ee8, ee9;
   double ee10, ee14, ee15, ee16, ee17, ee18, ee19;
   double ee20, ee21, ee22, ee23, ee26, ee27, ee29;
   double ee30, ee33, ee34, ee35, ee36;
   
   for (int j=0; j < nobs; j++) {
     
     talpha = tavec[j];
     tbeta = tbvec[j];
     mu = muvec[j];
     tsigma = tsigvec[j];
     
     for (int l=0; l < ncol; l++) {
       
       y = ymat(j, l);
       yi = xmat(j, l);
       
       if (arma::is_finite(y) & arma::is_finite(yi)) {
         
         w = wmat(j, l);
         
         ee1 = exp(tbeta);
         ee2 = 1 - ee1;
         ee3 = exp(tsigma);
         ee4 = R_pow(yi, ee2);
         ee6 = exp(-talpha);
         ee7 = ee4 * ee3;
         ee8 = 1 + ee6;
         ee9 = R_pow(yi, ee1);
         ee10 = R_pow(ee7, 2);
         ee14 = y - yi * (2/ee8 + mu/ee9 - 1);
         ee15 = 2 * ee2;
         ee16 = R_pow(ee3, 2);
         ee17 = log(yi);
         ee18 = R_pow(yi, ee15);
         ee19 = R_pow(ee8, 2);
         ee20 = mu/ee3;
         ee21 = ee20 + ee7 * ee14/ee10;
         ee22 = ee18 * ee16;
         ee23 = ee19 * ee16;
         ee26 = 2 * ee1;
         ee27 = ee4 - ee4 * ee1 * ee17;
         ee29 = R_pow(yi, (1 + ee15 - ee1)) * ee16;
         ee30 = R_pow(ee14, 2);
         ee33 = mu * ee18;
         ee34 = ee4 - 2 * (ee29/ee10);
         ee35 = ee22/ee10;
         ee36 = R_pow(yi, (ee26 - 1));
         
         out(j, 0) += w * (-(2 * (ee36 * ee6 * ee14/ee23)));
         out(j, 1) += w * ((ee21 * ee14/ee7 - 1) * ee1 * ee17);
         out(j, 2) += w * (-(ee14/(ee4 * ee16)));
         out(j, 3) += w * (1 - ee30/ee10);
         out(j, 4) += w * ((4 * (R_pow(yi, ee26) * ee6/ee19) + ee36 * (2 -
           4 * (ee6/ee8)) * ee14) * ee6/ee23);
         out(j, 5) += w * (-(yi * (2 * (ee14/ee10) + 2 * (R_pow(yi, (ee1 -
           1)) * ee21/ee3)) * ee6 * ee1 * ee17/ee19));
         out(j, 6) += w * (2 * (ee9 * ee6/ee23));
         out(j, 7) += w * (4 * (yi * ee6 * ee14/(ee19 * ee10)));
         out(j, 8) += w * (((((2 * (ee29 * ee14/ee10) + ee33) * ee1 * ee17 +
           ee14 * ee27) * ee3/ee10 + mu * (ee27/ee7 + ee18 * ee1 * ee3 *
           ee17/ee10)) * ee14/ee7 + (R_pow(ee21, 2) - ee35) * ee1 *
           ee17 - ee27/ee4) * ee1 * ee17);
         out(j, 9) += w * ((((1/ee3 - ee18 * ee3/ee10)/ee4 - ee7/ee10) * ee14 -
           ee20) * ee1 * ee17/ee3);
         out(j, 10) += w * (((((ee14 * ee34 - ee33)/ee4 - ee7 * ee21) * ee14 +
           ee22)/ee10 - 1) * ee1 * ee17);
         out(j, 11) += w * (1/ee16);
         out(j, 12) += w * (ee14 * (ee4 + R_pow(yi, (ee15 + ee1 - 1)))/ee10);
         out(j, 13) += w * (1 - ((ee34/ee4 - ee35) * ee30 + ee22)/ee10);
         
       }}}    
   
   return out;
   
 }
 
 //' @rdname condexd0
 // [[Rcpp::export]]
 arma::mat condexd34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, 
                     arma::mat ymat, arma::mat xmat, arma::mat wmat, 
                     arma::uvec dupid, int dcate, arma::uvec nhere)
 {
   
   arma::vec tavec = X1 * Rcpp::as<arma::vec>(pars[0]);
   arma::vec tbvec = X2 * Rcpp::as<arma::vec>(pars[1]);
   arma::vec muvec = X3 * Rcpp::as<arma::vec>(pars[2]);
   arma::vec tsigvec = X4 * Rcpp::as<arma::vec>(pars[3]);
   
   int nobs = ymat.n_rows;
   int ncol = ymat.n_cols;
   
   if (dcate == 1) {
     tavec = tavec.elem(dupid);
     tbvec = tbvec.elem(dupid);
     muvec = muvec.elem(dupid);
     tsigvec = tsigvec.elem(dupid);
   }
   
   double y, w, yi, talpha, tbeta, mu, tsigma;
   
   arma::mat out = arma::mat(nobs, 20, arma::fill::zeros);
   
   double ee1, ee2, ee3, ee4, ee5, ee6, ee8, ee9;
   double ee10, ee11, ee12, ee13, ee17, ee18;
   double ee20, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
   double ee30, ee32, ee34, ee35, ee38, ee39;
   double ee40, ee41, ee42, ee43, ee44, ee45, ee48;
   double ee52, ee53, ee54, ee58, ee59;
   double ee60, ee61, ee62, ee63, ee64, ee65, ee67, ee69;
   double ee70, ee71, ee74, ee75, ee77, ee78, ee79;
   double ee82;
   
   for (int j=0; j < nobs; j++) {
     
     talpha = tavec[j];
     tbeta = tbvec[j];
     mu = muvec[j];
     tsigma = tsigvec[j];
     
     for (int l=0; l < ncol; l++) {
       
       y = ymat(j, l);
       yi = xmat(j, l);
       
       if (arma::is_finite(y) & arma::is_finite(yi)) {
         
         w = wmat(j, l);
         
         ee1 = exp(tbeta);
         ee2 = 1 - ee1;
         ee3 = R_pow(yi, ee2);
         ee4 = exp(tsigma);
         ee5 = ee3 * ee4;
         ee6 = R_pow(ee5, 2);
         ee8 = exp(-talpha);
         ee9 = log(yi);
         ee10 = 1 + ee8;
         ee11 = 2 * ee2;
         ee12 = R_pow(yi, ee1);
         ee13 = R_pow(ee4, 2);
         ee17 = y - yi * (2/ee10 + mu/ee12 - 1);
         ee18 = R_pow(yi, ee11);
         ee20 = ee3 * ee1 * ee9;
         ee23 = ee3 - ee20;
         ee24 = R_pow(yi, (1 + ee11 - ee1));
         ee25 = 4 * ee2;
         ee26 = R_pow(ee10, 2);
         ee27 = mu * ee18;
         ee28 = ee24 * ee13;
         ee29 = R_pow(yi, ee25);
         ee30 = 2 * (ee29 * ee13/ee6);
         ee32 = mu/ee4 + ee5 * ee17/ee6;
         ee34 = 2 * ee3 + ee3;
         ee35 = ee17 * ee23;
         ee38 = ee3 - 2 * (ee28/ee6);
         ee39 = R_pow(yi, (2 - ee1));
         ee40 = ee26 * ee6;
         ee41 = ee23/ee5;
         ee42 = 2 * ee1;
         ee43 = 2 * ee18;
         ee44 = ee18 * ee1;
         ee45 = ee18 * ee4;
         ee48 = ((2 * (ee28 * ee17/ee6) + ee27) * ee1 * ee9 + ee35) *  ee4/ee6 + mu * (ee41 + ee44 * ee4 * ee9/ee6);
         ee52 = 1/ee4 - ee45/ee6;
         ee53 = 2 - 4 * (ee8/ee10);
         ee54 = R_pow(yi, (1 + ee25 - ee1));
         ee58 = ee17 * ee38 - ee27;
         ee59 = 2 * ee39;
         ee60 = ee43 - ee30;
         ee61 = 8 * (ee54 * ee13/ee6);
         ee62 = ee3 - (ee34 - ee20) * ee1 * ee9;
         ee63 = ee3 * ee23;
         ee64 = ee18 - ee30;
         ee65 = ee26 * ee13;
         ee67 = ee52 * ee1 * ee9;
         ee69 = ee34 * ee23 + 2 * (ee29 * ee1 * ee13 * ee9/ee6);
         ee70 = (ee59 - 4 * (R_pow(yi, (2 + ee11 - ee1)) * ee13/ee6))/ee3;
         ee71 = (ee43 + ee18 - ee30) * ee13;
         ee74 = ee35 + ee27 * ee1 * ee9;
         ee75 = 2 * ee24;
         ee77 = yi * ee8;
         ee78 = ee3 - (8 * ee24 - ee61) * ee13/ee6;
         ee79 = R_pow(yi, ee42);
         ee82 = R_pow(yi, (3 - ee42)) * ee13/ee6;
         
         out(j, 0) += w * (-((ee53 * (2 * ee79 + 4 * ee79) * ee8/ee26 + R_pow(yi, (ee42 -
           1)) * (2 - ((4 * (1 + 2 * ee8) + 4 * ee10 - 16 *
           ee8)/ee10 + 4) * ee8/ee10) * ee17) * ee8/ee65));
         out(j, 1) += w * (yi * ((ee17/ee6 + R_pow(yi, (ee1 - 1)) * ee32/ee4) *
           ee53 + 8 * (ee77/ee40)) * ee8 * ee1 * ee9/ee26);
         out(j, 2) += w * (-(ee12 * ee53 * ee8/ee65));
         out(j, 3) += w * (-(yi * (2 * (ee53 * ee17) + 8 * (ee77/ee26)) * ee8/
           ee40));
         out(j, 4) += w * (-(((4 * ee39 * ee1 * ee4 * ee9 * ee32 + ee12 * (2 *
           ee23 + 4 * (ee24 * ee1 * ee13 * ee9/ee6)) * ee17)/ee6 +
           2 * (ee12 * ee48/ee4)) * ee8 * ee1 * ee9/ee26));
         out(j, 5) += w * (-((2 * (ee12 * ee52/ee4) - 2 * (ee39/ee6)) * ee8 *
           ee1 * ee9/ee26));
         out(j, 6) += w * (-(((ee70 - 2 * ee82) * ee17 + 2 * (ee12 * ee58) -
           2 * (ee39 * ee4 * ee32)) * ee8 * ee1 * ee9/ee40));
         out(j, 7) += w * 0;
         out(j, 8) += w * (-((ee59 + 2 * R_pow(yi, (ee11 + ee1))) * ee8/ee40));
         out(j, 9) += w * ((ee70 + 2 * (ee12 * ee38) - 4 * ee82) * ee8 * ee17/
           ee40);
         out(j, 10) += w * (((((((2 * (ee3 * (ee63 - ee44 * ee9)) + 8 * (ee54 *
           ee1 * ee13 * ee9/ee6)) * ee17 + 4 * (ee18 * ee74)) * ee13/
             ee6 + mu * ee34 * ee23) * ee1 * ee9 + ee17 * ee62) * ee4/
               ee6 + mu * (ee69 * ee1 * ee4 * ee9/ee6 + ee62/ee5)) * ee17/
                 ee5 + (3 * (ee48 * ee32) - ee69 * ee13/ee6) * ee1 * ee9 - ee62/
                   ee3) * ee1 * ee9);
         out(j, 11) += w * (((((ee60 * ee1 * ee9 - ee63) * ee4/ee6 + ee41) *
           ee17/ee3 - ee48)/ee4 + 2 * (ee67 * ee32)) * ee1 * ee9);
         out(j, 12) += w * (((((ee74 * (1 - 2 * (ee18 * ee13/ee6)) + (ee75 +
           4 * ee24 - ee61) * ee1 * ee13 * ee9 * ee17/ee6 + mu * (ee1 *
           ee9 * ee64 - ee63))/ee3 - ee3 * ee48 * ee4) * ee17 + (2 *
           (ee58 * ee32) - ee60 * ee4) * ee1 * ee4 * ee9)/ee6 - (1/ee3 -
           ee3 * ee13/ee6) * ee23) * ee1 * ee9);
         out(j, 13) += w * (-(2 * (ee67/ee4)));
         out(j, 14) += w * (-(((ee60/ee3 + ee4 * (ee3 * ee52 - 2 * (ee24 *
           ee4/ee6)) + ee3) * ee17 - (ee27 + ee45 * ee32)) * ee1 * ee9/
             ee6));
         out(j, 15) += w * (((((ee17 * ee78 - mu * ee64)/ee3 - (ee32 * ee38 +
           2 * (ee3 * ee58 * ee4/ee6)) * ee4) * ee17 + ee71)/ee6 - 1) *
           ee1 * ee9);
         out(j, 16) += w * 0;
         out(j, 17) += w * (-(2 * (ee18/ee6)));
         out(j, 18) += w * ((ee64/ee3 + ee3 - (ee75 + ee24 + ee24) * ee13/
           ee6) * ee17/ee6);
         out(j, 19) += w * (1 - ((ee78/ee3 - ee34 * ee13 * ee38/ee6) * R_pow(ee17, 2) +
           ee71)/ee6);
         
       }}}    
   
   return out;
   
 }
 
 //' Conditional extreme value model negative log-likelihood
 //' using sparse design matrices
 //'
 //' @param pars a list of vectors of coefficients for each conditional EVD parameter
 //' @param X1 a sparse design matrix for (transformed) alpha
 //' @param X2 a sparse design matrix for (transformed) beta
 //' @param X3 a sparse design matrix for mu
 //' @param X4 a sparse design matrix for (transformed) sigma
 //' @param ymat a matrix
 //' @param xmat a matrix
 //' @param wmat a matrix
 //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
 //' @return condexspd0 a scalar, the negative log-likelihood
 //' @return condexspd12 a matrix, first then second derivatives w.r.t. parameters
 //' @return condexspd34 a matrix, third then fourth derivatives w.r.t. parameters
 //' @examples
 //' ## to follow
 //' @export
 // [[Rcpp::export]]
 double condexspd0(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::sp_mat X3, arma::sp_mat X4, 
                 arma::mat ymat, arma::mat xmat, arma::mat wmat, 
                 arma::uvec dupid, int dcate, arma::uvec nhere)
 {
   
   arma::vec tavec = X1 * Rcpp::as<arma::vec>(pars[0]);
   arma::vec tbvec = X2 * Rcpp::as<arma::vec>(pars[1]);
   arma::vec muvec = X3 * Rcpp::as<arma::vec>(pars[2]);
   arma::vec tsigvec = X4 * Rcpp::as<arma::vec>(pars[3]);
   
   int nobs = ymat.n_rows;
   int ncol = ymat.n_cols;
   
   if (dcate == 1) {
     tavec = tavec.elem(dupid);
     tbvec = tbvec.elem(dupid);
     muvec = muvec.elem(dupid);
     tsigvec = tsigvec.elem(dupid);
   }
   
   double y, w, yi, talpha, tbeta, mu, tsigma;
   double alpha, beta, sigma, mu2, sigma2, yib, res;
   double nllh=0.0;
   
   for (int j=0; j < nobs; j++) {
     
     talpha = tavec[j];
     tbeta = tbvec[j];
     mu = muvec[j];
     tsigma = tsigvec[j];
     
     for (int l=0; l < ncol; l++) {
       
       y = ymat(j, l);
       yi = xmat(j, l);
       
       if (arma::is_finite(y) & arma::is_finite(yi)) {
         
         w = wmat(j, l);
         
         alpha = (2.0 / (1.0 + exp(-talpha)) - 1.0);
         beta = 1.0 - exp(tbeta);
         sigma = exp(tsigma);
         yib = R_pow(yi, beta);
         mu2 = alpha * yi + yib * mu;
         sigma2 = sigma * yib;
         res = (y - mu2) / sigma2;
         nllh += w * (0.5 * logtwopi + log(sigma2) + 0.5 * res * res); 
         
       }
     }
   }
   
   return(nllh);
   
 }
 
 //' @rdname condexspd0
 // [[Rcpp::export]]
 arma::mat condexspd12(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::sp_mat X3, arma::sp_mat X4, 
                     arma::mat ymat, arma::mat xmat, arma::mat wmat, 
                     arma::uvec dupid, int dcate, arma::uvec nhere)
 {
   
   arma::vec tavec = X1 * Rcpp::as<arma::vec>(pars[0]);
   arma::vec tbvec = X2 * Rcpp::as<arma::vec>(pars[1]);
   arma::vec muvec = X3 * Rcpp::as<arma::vec>(pars[2]);
   arma::vec tsigvec = X4 * Rcpp::as<arma::vec>(pars[3]);
   
   int nobs = ymat.n_rows;
   int ncol = ymat.n_cols;
   
   if (dcate == 1) {
     tavec = tavec.elem(dupid);
     tbvec = tbvec.elem(dupid);
     muvec = muvec.elem(dupid);
     tsigvec = tsigvec.elem(dupid);
   }
   
   double y, w, yi, talpha, tbeta, mu, tsigma;
   
   arma::mat out = arma::mat(nobs, 14, arma::fill::zeros);
   
   double ee1, ee2, ee3, ee4, ee6, ee7, ee8, ee9;
   double ee10, ee14, ee15, ee16, ee17, ee18, ee19;
   double ee20, ee21, ee22, ee23, ee26, ee27, ee29;
   double ee30, ee33, ee34, ee35, ee36;
   
   for (int j=0; j < nobs; j++) {
     
     talpha = tavec[j];
     tbeta = tbvec[j];
     mu = muvec[j];
     tsigma = tsigvec[j];
     
     for (int l=0; l < ncol; l++) {
       
       y = ymat(j, l);
       yi = xmat(j, l);
       
       if (arma::is_finite(y) & arma::is_finite(yi)) {
         
         w = wmat(j, l);
         
         ee1 = exp(tbeta);
         ee2 = 1 - ee1;
         ee3 = exp(tsigma);
         ee4 = R_pow(yi, ee2);
         ee6 = exp(-talpha);
         ee7 = ee4 * ee3;
         ee8 = 1 + ee6;
         ee9 = R_pow(yi, ee1);
         ee10 = R_pow(ee7, 2);
         ee14 = y - yi * (2/ee8 + mu/ee9 - 1);
         ee15 = 2 * ee2;
         ee16 = R_pow(ee3, 2);
         ee17 = log(yi);
         ee18 = R_pow(yi, ee15);
         ee19 = R_pow(ee8, 2);
         ee20 = mu/ee3;
         ee21 = ee20 + ee7 * ee14/ee10;
         ee22 = ee18 * ee16;
         ee23 = ee19 * ee16;
         ee26 = 2 * ee1;
         ee27 = ee4 - ee4 * ee1 * ee17;
         ee29 = R_pow(yi, (1 + ee15 - ee1)) * ee16;
         ee30 = R_pow(ee14, 2);
         ee33 = mu * ee18;
         ee34 = ee4 - 2 * (ee29/ee10);
         ee35 = ee22/ee10;
         ee36 = R_pow(yi, (ee26 - 1));
         
         out(j, 0) += w * (-(2 * (ee36 * ee6 * ee14/ee23)));
         out(j, 1) += w * ((ee21 * ee14/ee7 - 1) * ee1 * ee17);
         out(j, 2) += w * (-(ee14/(ee4 * ee16)));
         out(j, 3) += w * (1 - ee30/ee10);
         out(j, 4) += w * ((4 * (R_pow(yi, ee26) * ee6/ee19) + ee36 * (2 -
           4 * (ee6/ee8)) * ee14) * ee6/ee23);
         out(j, 5) += w * (-(yi * (2 * (ee14/ee10) + 2 * (R_pow(yi, (ee1 -
           1)) * ee21/ee3)) * ee6 * ee1 * ee17/ee19));
         out(j, 6) += w * (2 * (ee9 * ee6/ee23));
         out(j, 7) += w * (4 * (yi * ee6 * ee14/(ee19 * ee10)));
         out(j, 8) += w * (((((2 * (ee29 * ee14/ee10) + ee33) * ee1 * ee17 +
           ee14 * ee27) * ee3/ee10 + mu * (ee27/ee7 + ee18 * ee1 * ee3 *
           ee17/ee10)) * ee14/ee7 + (R_pow(ee21, 2) - ee35) * ee1 *
           ee17 - ee27/ee4) * ee1 * ee17);
         out(j, 9) += w * ((((1/ee3 - ee18 * ee3/ee10)/ee4 - ee7/ee10) * ee14 -
           ee20) * ee1 * ee17/ee3);
         out(j, 10) += w * (((((ee14 * ee34 - ee33)/ee4 - ee7 * ee21) * ee14 +
           ee22)/ee10 - 1) * ee1 * ee17);
         out(j, 11) += w * (1/ee16);
         out(j, 12) += w * (ee14 * (ee4 + R_pow(yi, (ee15 + ee1 - 1)))/ee10);
         out(j, 13) += w * (1 - ((ee34/ee4 - ee35) * ee30 + ee22)/ee10);
         
       }}}    
   
   return out;
   
 }
 
 //' @rdname condexspd0
 // [[Rcpp::export]]
 arma::mat condexspd34(Rcpp::List pars, arma::sp_mat X1, arma::sp_mat X2, arma::sp_mat X3, arma::sp_mat X4, 
                     arma::mat ymat, arma::mat xmat, arma::mat wmat, 
                     arma::uvec dupid, int dcate, arma::uvec nhere)
 {
   
   arma::vec tavec = X1 * Rcpp::as<arma::vec>(pars[0]);
   arma::vec tbvec = X2 * Rcpp::as<arma::vec>(pars[1]);
   arma::vec muvec = X3 * Rcpp::as<arma::vec>(pars[2]);
   arma::vec tsigvec = X4 * Rcpp::as<arma::vec>(pars[3]);
   
   int nobs = ymat.n_rows;
   int ncol = ymat.n_cols;
   
   if (dcate == 1) {
     tavec = tavec.elem(dupid);
     tbvec = tbvec.elem(dupid);
     muvec = muvec.elem(dupid);
     tsigvec = tsigvec.elem(dupid);
   }
   
   double y, w, yi, talpha, tbeta, mu, tsigma;
   
   arma::mat out = arma::mat(nobs, 20, arma::fill::zeros);
   
   double ee1, ee2, ee3, ee4, ee5, ee6, ee8, ee9;
   double ee10, ee11, ee12, ee13, ee17, ee18;
   double ee20, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
   double ee30, ee32, ee34, ee35, ee38, ee39;
   double ee40, ee41, ee42, ee43, ee44, ee45, ee48;
   double ee52, ee53, ee54, ee58, ee59;
   double ee60, ee61, ee62, ee63, ee64, ee65, ee67, ee69;
   double ee70, ee71, ee74, ee75, ee77, ee78, ee79;
   double ee82;
   
   for (int j=0; j < nobs; j++) {
     
     talpha = tavec[j];
     tbeta = tbvec[j];
     mu = muvec[j];
     tsigma = tsigvec[j];
     
     for (int l=0; l < ncol; l++) {
       
       y = ymat(j, l);
       yi = xmat(j, l);
       
       if (arma::is_finite(y) & arma::is_finite(yi)) {
         
         w = wmat(j, l);
         
         ee1 = exp(tbeta);
         ee2 = 1 - ee1;
         ee3 = R_pow(yi, ee2);
         ee4 = exp(tsigma);
         ee5 = ee3 * ee4;
         ee6 = R_pow(ee5, 2);
         ee8 = exp(-talpha);
         ee9 = log(yi);
         ee10 = 1 + ee8;
         ee11 = 2 * ee2;
         ee12 = R_pow(yi, ee1);
         ee13 = R_pow(ee4, 2);
         ee17 = y - yi * (2/ee10 + mu/ee12 - 1);
         ee18 = R_pow(yi, ee11);
         ee20 = ee3 * ee1 * ee9;
         ee23 = ee3 - ee20;
         ee24 = R_pow(yi, (1 + ee11 - ee1));
         ee25 = 4 * ee2;
         ee26 = R_pow(ee10, 2);
         ee27 = mu * ee18;
         ee28 = ee24 * ee13;
         ee29 = R_pow(yi, ee25);
         ee30 = 2 * (ee29 * ee13/ee6);
         ee32 = mu/ee4 + ee5 * ee17/ee6;
         ee34 = 2 * ee3 + ee3;
         ee35 = ee17 * ee23;
         ee38 = ee3 - 2 * (ee28/ee6);
         ee39 = R_pow(yi, (2 - ee1));
         ee40 = ee26 * ee6;
         ee41 = ee23/ee5;
         ee42 = 2 * ee1;
         ee43 = 2 * ee18;
         ee44 = ee18 * ee1;
         ee45 = ee18 * ee4;
         ee48 = ((2 * (ee28 * ee17/ee6) + ee27) * ee1 * ee9 + ee35) *  ee4/ee6 + mu * (ee41 + ee44 * ee4 * ee9/ee6);
         ee52 = 1/ee4 - ee45/ee6;
         ee53 = 2 - 4 * (ee8/ee10);
         ee54 = R_pow(yi, (1 + ee25 - ee1));
         ee58 = ee17 * ee38 - ee27;
         ee59 = 2 * ee39;
         ee60 = ee43 - ee30;
         ee61 = 8 * (ee54 * ee13/ee6);
         ee62 = ee3 - (ee34 - ee20) * ee1 * ee9;
         ee63 = ee3 * ee23;
         ee64 = ee18 - ee30;
         ee65 = ee26 * ee13;
         ee67 = ee52 * ee1 * ee9;
         ee69 = ee34 * ee23 + 2 * (ee29 * ee1 * ee13 * ee9/ee6);
         ee70 = (ee59 - 4 * (R_pow(yi, (2 + ee11 - ee1)) * ee13/ee6))/ee3;
         ee71 = (ee43 + ee18 - ee30) * ee13;
         ee74 = ee35 + ee27 * ee1 * ee9;
         ee75 = 2 * ee24;
         ee77 = yi * ee8;
         ee78 = ee3 - (8 * ee24 - ee61) * ee13/ee6;
         ee79 = R_pow(yi, ee42);
         ee82 = R_pow(yi, (3 - ee42)) * ee13/ee6;
         
         out(j, 0) += w * (-((ee53 * (2 * ee79 + 4 * ee79) * ee8/ee26 + R_pow(yi, (ee42 -
           1)) * (2 - ((4 * (1 + 2 * ee8) + 4 * ee10 - 16 *
           ee8)/ee10 + 4) * ee8/ee10) * ee17) * ee8/ee65));
         out(j, 1) += w * (yi * ((ee17/ee6 + R_pow(yi, (ee1 - 1)) * ee32/ee4) *
           ee53 + 8 * (ee77/ee40)) * ee8 * ee1 * ee9/ee26);
         out(j, 2) += w * (-(ee12 * ee53 * ee8/ee65));
         out(j, 3) += w * (-(yi * (2 * (ee53 * ee17) + 8 * (ee77/ee26)) * ee8/
           ee40));
         out(j, 4) += w * (-(((4 * ee39 * ee1 * ee4 * ee9 * ee32 + ee12 * (2 *
           ee23 + 4 * (ee24 * ee1 * ee13 * ee9/ee6)) * ee17)/ee6 +
           2 * (ee12 * ee48/ee4)) * ee8 * ee1 * ee9/ee26));
         out(j, 5) += w * (-((2 * (ee12 * ee52/ee4) - 2 * (ee39/ee6)) * ee8 *
           ee1 * ee9/ee26));
         out(j, 6) += w * (-(((ee70 - 2 * ee82) * ee17 + 2 * (ee12 * ee58) -
           2 * (ee39 * ee4 * ee32)) * ee8 * ee1 * ee9/ee40));
         out(j, 7) += w * 0;
         out(j, 8) += w * (-((ee59 + 2 * R_pow(yi, (ee11 + ee1))) * ee8/ee40));
         out(j, 9) += w * ((ee70 + 2 * (ee12 * ee38) - 4 * ee82) * ee8 * ee17/
           ee40);
         out(j, 10) += w * (((((((2 * (ee3 * (ee63 - ee44 * ee9)) + 8 * (ee54 *
           ee1 * ee13 * ee9/ee6)) * ee17 + 4 * (ee18 * ee74)) * ee13/
             ee6 + mu * ee34 * ee23) * ee1 * ee9 + ee17 * ee62) * ee4/
               ee6 + mu * (ee69 * ee1 * ee4 * ee9/ee6 + ee62/ee5)) * ee17/
                 ee5 + (3 * (ee48 * ee32) - ee69 * ee13/ee6) * ee1 * ee9 - ee62/
                   ee3) * ee1 * ee9);
         out(j, 11) += w * (((((ee60 * ee1 * ee9 - ee63) * ee4/ee6 + ee41) *
           ee17/ee3 - ee48)/ee4 + 2 * (ee67 * ee32)) * ee1 * ee9);
         out(j, 12) += w * (((((ee74 * (1 - 2 * (ee18 * ee13/ee6)) + (ee75 +
           4 * ee24 - ee61) * ee1 * ee13 * ee9 * ee17/ee6 + mu * (ee1 *
           ee9 * ee64 - ee63))/ee3 - ee3 * ee48 * ee4) * ee17 + (2 *
           (ee58 * ee32) - ee60 * ee4) * ee1 * ee4 * ee9)/ee6 - (1/ee3 -
           ee3 * ee13/ee6) * ee23) * ee1 * ee9);
         out(j, 13) += w * (-(2 * (ee67/ee4)));
         out(j, 14) += w * (-(((ee60/ee3 + ee4 * (ee3 * ee52 - 2 * (ee24 *
           ee4/ee6)) + ee3) * ee17 - (ee27 + ee45 * ee32)) * ee1 * ee9/
             ee6));
         out(j, 15) += w * (((((ee17 * ee78 - mu * ee64)/ee3 - (ee32 * ee38 +
           2 * (ee3 * ee58 * ee4/ee6)) * ee4) * ee17 + ee71)/ee6 - 1) *
           ee1 * ee9);
         out(j, 16) += w * 0;
         out(j, 17) += w * (-(2 * (ee18/ee6)));
         out(j, 18) += w * ((ee64/ee3 + ee3 - (ee75 + ee24 + ee24) * ee13/
           ee6) * ee17/ee6);
         out(j, 19) += w * (1 - ((ee78/ee3 - ee34 * ee13 * ee38/ee6) * R_pow(ee17, 2) +
           ee71)/ee6);
         
       }}}    
   
   return out;
   
 }
 