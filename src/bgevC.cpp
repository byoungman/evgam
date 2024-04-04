// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

const double xieps = 0.0;

double ell1(double a, double xi)
{
return R_pow(-log(a), -xi);
}

double ell2(double a)
{
  return log(-log(a));
}

double iF(double p, double qalpha, double sbeta, double xi, double alpha, double beta)
{
  double out = qalpha;
  out += (ell1(p, xi) - ell1(alpha, xi)) * sbeta / (ell1(1 - .5 * beta, xi) - ell1(.5 * beta, xi));
  return out;
}

// //' Blended generalized extreme value (bGEV) distribution negative log-likelihood
// //'
// //' @param pars a list of vectors of coefficients for each bGEV parameter
// //' @param X1 a design matrix for the bGEV location parameter
// //' @param X2 a design matrix for the bGEV log scale parameter
// //' @param X3 a design matrix for the bGEV shape parameter
// //' @param yvec a vector
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @return bgevd0 a scalar, the negative log-likelihood
// //' @return bgevd12 a matrix, first then second derivatives w.r.t. bGEV parameters
// //' @return bgevd34 a matrix, third then fourth derivatives w.r.t. bGEV parameters
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double bgevd0(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec, arma::uvec dupid, int dcate, arma::vec other)
{
    
arma::vec qavec = X1 * Rcpp::as<arma::vec>(pars[0]);
arma::vec lsbvec = X2 * Rcpp::as<arma::vec>(pars[1]);
arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
int nobs = yvec.size();

double psuba = other[0];
double psubb = other[1];
double alpha = other[2];
double beta = other[3];
double hbeta = 0.5 * beta;

if (dcate == 1) {
    qavec = qavec.elem(dupid);
    lsbvec = lsbvec.elem(dupid);
    txivec = txivec.elem(dupid);
}

double y, qalpha, lsbeta, txi, xi, tqalpha, sbeta, tsbeta, iFa, iFb;
double z1d, z1, t1, z2d, z2, t2;
double y2, dx, px, lH, fiF, giG;
double nllh=0.0;


for (int j=0; j < nobs; j++) {

y = yvec[j];
qalpha = qavec[j];
lsbeta = lsbvec[j];
txi = txivec[j];
sbeta = exp(lsbeta);
xi = 1.5 / (1.0 + exp(-txi)) - 0.5;

iFa = iF(psuba, qalpha, sbeta, xi, alpha, beta);
iFb = iF(psubb, qalpha, sbeta, xi, alpha, beta);

tqalpha = iFa - (iFb - iFa) * (ell2(alpha) - ell2(psuba)) / (ell2(psuba) - ell2(psubb));
tsbeta = (iFb - iFa) * (ell2(hbeta) - ell2(1.0 - hbeta)) / (ell2(psuba) - ell2(psubb));

if (y < iFa) { // Gumbel
  
  z2d = (ell2(hbeta) - ell2(1.0 - hbeta)) / tsbeta;
  z2 = (y - tqalpha) * (ell2(hbeta) - ell2(1.0 - hbeta)) / tsbeta - ell2(alpha);
  t2 = exp(-z2);
  // G = exp(-t2);
  // nllh += log(G * t2 * z2d);
  // nllh += log(G) + log(t2) + log(z2d);
  nllh += t2 + z2 - log(z2d);
  
} else {
  
  if (y > iFb) { // GEV
    
    z1d = (ell1(1.0 - hbeta, xi) - ell1(hbeta, xi)) / sbeta;
    z1 = (y - qalpha) * z1d + ell1(alpha, xi);
    t1 = R_pow(z1, -1/xi);
    // F = exp(-t1);
    // nllh += log(F * R_pow(z1, -(1 + 1/xi)) * z1d / xi);
    // nllh += log(F) + log(R_pow(z1, -(1 + 1/xi))) + log(ell1(1.0 - hbeta, xi) - ell1(hbeta, xi)) - log(sb)- log(xi);
    nllh += t1 + (1 + 1 / xi) * log(z1) - log(z1d) + log(xi);
    
  } else { // bGEV
    
    z2d = (ell2(hbeta) - ell2(1.0 - hbeta)) / tsbeta;
    z2 = (y - tqalpha) * (ell2(hbeta) - ell2(1.0 - hbeta)) / tsbeta - ell2(alpha);
    t2 = exp(-z2);
    z1d = (ell1(1.0 - hbeta, xi) - ell1(hbeta, xi)) / sbeta;
    z1 = (y - qalpha) * z1d + ell1(alpha, xi);
    t1 = R_pow(z1, -1/xi);
    
    y2 = (y - iFa) / (iFb - iFa);
    px = (R_pow(y2, 5)*(70*R_pow(y2, 4)-315*R_pow(y2, 3)+540*R_pow(y2, 2)-420*(y2)+126));
    dx = 630 * (R_pow(y2, 4)) * (R_pow(1 - y2, 4)) / (iFb - iFa);
    // lF = -(t1 + (1 + 1 / xi) * log(z1) - log(z1d) + log(xi));
    // lG = -(t2 + z2 - log(z2d));
    // H = exp(px * lF + (1 - px) * lG);
    // lH = px * log(F) + (1 - px) * log(G);
    lH = px * t1 + (1 - px) * t2;
    nllh += lH;
    // nllh += log(- dx * t1 + px * f / F + dx * t2 + (1 - px) * g / G);
    fiF = R_pow(z1, -(1 + 1/xi)) * z1d / xi;
    giG = t2 * z2d;
    nllh += -log(- dx * t1 + px * fiF + dx * t2 + (1 - px) * giG);
    
    // ee4 = 1.5/(1 + exp(-txi)) - 0.5;
    // ee5 = -log(alpha);
    // ee7 = 1/R_pow(ee5, ee4);
    // ee8 = -log(psuba);
    // ee14 = 1/R_pow((-log(1 - hbeta)), ee4) - 1/R_pow((-log(hbeta)), ee4);
    // ee16 = 1/R_pow(ee8, ee4);
    // ee17 = exp(lsbeta);
    // ee18 = -log(psubb);
    // ee22 = ee7 + 1/R_pow(ee18, ee4) - (ee7 + ee16);
    // ee23 = ee22 * ee17;
    // ee27 = y - ((ee16 - ee7) * ee17/ee14 + qalpha);
    // ee29 = ee14 * ee27/ee23;
    // ee30 = log(ee8);
    // ee31 = ee30 - log(ee18);
    // ee35 = (ee14 * (70 * ee29 - 315) * ee27/ee23 + 540) * R_pow(ee29, 2) + 
    // 126 - 420 * ee29;
    //             ee37 = ee35 * R_pow(ee29, 5);
    //             ee40 = ee14 * (y - qalpha)/ee17 + ee7;
    //             ee41 = 1/ee4;
    //             ee42 = exp(-(ee14 * ee31 * (y - ((ee16 - (ee22 * (log(ee5) - 
    //             ee30)/ee31 + ee7)) * ee17/ee14 + qalpha))/ee23 - log(ee5)));
    //             ee43 = R_pow(ee40, ee41);
    //             ee44 = (1 - ee37) * ee42;
    //             nllh += ee37/ee43 + ee44 + log(ee22) + lsbeta - (log((ee35 * ee14 * 
    //             ee27/(R_pow(ee40, (1 + ee41)) * ee4 * ee17) + R_pow((1 - ee29), 4) * 
    // (630 * ee42 - 630/ee43)) * R_pow(ee29, 4) + ee44 * ee31) + log(ee14));

  }
}

}

return(nllh);

}

// // [[Rcpp::export]]
// double bgevd0_test(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec, arma::uvec dupid, int dcate, arma::vec other)
// {
//   
//   arma::vec qavec = X1 * Rcpp::as<arma::vec>(pars[0]);
//   arma::vec lsbvec = X2 * Rcpp::as<arma::vec>(pars[1]);
//   arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
//   int nobs = yvec.size();
//   
//   double psuba = other[0];
//   double psubb = other[1];
//   double alpha = other[2];
//   double beta = other[3];
//   double hbeta = 0.5 * beta;
//   
//   if (dcate == 1) {
//     qavec = qavec.elem(dupid);
//     lsbvec = lsbvec.elem(dupid);
//     txivec = txivec.elem(dupid);
//   }
//   
//   double y, qalpha, lsbeta, txi, xi, tqalpha, sbeta, tsbeta, iFa, iFb;
//   double z1d, z1, t1, z2d, z2, t2;
//   double ee1, ee2;
//   double nllh=0.0;
//   
//   for (int j=0; j < nobs; j++) {
//     
//     y = yvec[j];
//     qalpha = qavec[j];
//     lsbeta = lsbvec[j];
//     txi = txivec[j];
//     sbeta = exp(lsbeta);
//     xi = 1.5 / (1.0 + exp(-txi)) - 0.5;
//     
//     iFa = iF(psuba, qalpha, sbeta, xi, alpha, beta);
//     iFb = iF(psubb, qalpha, sbeta, xi, alpha, beta);
//     
//     tqalpha = iFa - (iFb - iFa) * (ell2(alpha) - ell2(psuba)) / (ell2(psuba) - ell2(psubb));
//     tsbeta = (iFb - iFa) * (ell2(hbeta) - ell2(1.0 - hbeta)) / (ell2(psuba) - ell2(psubb));
//     
//     if (y < iFa) { // Gumbel
//       
//       // z2d = (ell2(hbeta) - ell2(1.0 - hbeta)) / tsbeta;
//       // z2 = (y - tqalpha) * (ell2(hbeta) - ell2(1.0 - hbeta)) / tsbeta - ell2(alpha);
//       // t2 = exp(-z2);
//       // // G = exp(-t2);
//       // // nllh += log(G * t2 * z2d);
//       // // nllh += log(G) + log(t2) + log(z2d);
//       // nllh += -log(exp(-exp(-((y - (iFa - (iFb - iFa) * (log(-log(alpha)) - 
//       //   log(-log(psuba)))/(log(-log(psuba)) - log(-log(psubb))))) * 
//       //   (log(-log(hbeta)) - log(-log(1 - hbeta)))/((iFb - iFa) * 
//       //   (log(-log(hbeta)) - log(-log(1 - hbeta)))/(log(-log(psuba)) - 
//       //   log(-log(psubb)))) - log(-log(alpha)))))) - log(exp(-((y - 
//       //   (iFa - (iFb - iFa) * (log(-log(alpha)) - log(-log(psuba)))/(log(-log(psuba)) - 
//       //   log(-log(psubb))))) * (log(-log(hbeta)) - log(-log(1 - 
//       //   hbeta)))/((iFb - iFa) * (log(-log(hbeta)) - log(-log(1 - 
//       //   hbeta)))/(log(-log(psuba)) - log(-log(psubb)))) - log(-log(alpha))))) - 
//       //   log(((log(-log(hbeta)) - log(-log(1 - hbeta)))/((iFb - iFa) * 
//       //   (log(-log(hbeta)) - log(-log(1 - hbeta)))/(log(-log(psuba)) - 
//       //   log(-log(psubb))))));
//       nllh += -log(exp(-exp(-((y-((qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))*(log(-log(alpha))-log(-log(psuba)))/(log(-log(psuba))-log(-log(psubb)))))*(log(-log(hbeta))-log(-log(1-hbeta)))/(((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))*(log(-log(hbeta))-log(-log(1-hbeta)))/(log(-log(psuba))-log(-log(psubb))))-log(-log(alpha))))))-log(exp(-((y-((qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))*(log(-log(alpha))-log(-log(psuba)))/(log(-log(psuba))-log(-log(psubb)))))*(log(-log(hbeta))-log(-log(1-hbeta)))/(((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))*(log(-log(hbeta))-log(-log(1-hbeta)))/(log(-log(psuba))-log(-log(psubb))))-log(-log(alpha)))))-log(((log(-log(hbeta))-log(-log(1-hbeta)))/(((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))*(log(-log(hbeta))-log(-log(1-hbeta)))/(log(-log(psuba))-log(-log(psubb))))));
//       
//     } else {
//       
//       if (y > iFb) { // GEV
//         
//         z1d = (ell1(1.0 - hbeta, xi) - ell1(hbeta, xi)) / sbeta;
//         z1 = (y - qalpha) * z1d + ell1(alpha, xi);
//         t1 = R_pow(z1, -1/xi);
//         // F = exp(-t1);
//         // nllh += log(F * R_pow(z1, -(1 + 1/xi)) * z1d / xi);
//         // nllh += log(F) + log(R_pow(z1, -(1 + 1/xi))) + log(ell1(1.0 - hbeta, xi) - ell1(hbeta, xi)) - log(sb)- log(xi);
//         // nllh += log((1.5/(1+exp(-txi))-.5))-log(exp(-(R_pow(((y-qalpha)*((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))/exp(lsbeta)+(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5)))),-1/(1.5/(1+exp(-txi))-.5)))))+(1/1+(1.5/(1+exp(-txi))-.5))*log(((y-qalpha)*((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))/exp(lsbeta)+(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5)))))-log((((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))/exp(lsbeta)));
//         nllh += log((1.5/(1+exp(-txi))-.5))+(R_pow(((y-qalpha)*((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))/exp(lsbeta)+(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5)))),-1/(1.5/(1+exp(-txi))-.5)))+(1+1/(1.5/(1+exp(-txi))-.5))*log(((y-qalpha)*((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))/exp(lsbeta)+(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5)))))-log((((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))/exp(lsbeta)));
//         
//       } else { // bGEV
//         
//         // z2d = (ell2(hbeta) - ell2(1.0 - hbeta)) / tsbeta;
//         // z2 = (y - qa2) * (ell2(hbeta) - ell2(1.0 - hbeta)) / tsbeta - ell2(alpha);
//         // t2 = exp(-z2);
//         // z1d = (ell1(1.0 - hbeta, xi) - ell1(hbeta, xi)) / sbeta;
//         // z1 = (y - qa) * z1d + ell1(alpha, xi);
//         // t1 = R_pow(z1, -1/xi);
//         // 
//         // px = ...;
//         // dx = ...;
//         // H = exp(px * log(F) + (1 - px) * log(G));
//         // // lH = px * log(F) + (1 - px) * log(G);
//         // lH = px * t1 + (1 - px) * t2;
//         // nllh += lH;
//         // // nllh += log(- dx * t1 + px * f / F + dx * t2 + (1 - px) * g / G);
//         // fiF = R_pow(z1, -(1 + 1/xi)) * z1d / xi;
//         // giG = t2 * z2d;
//         // nllh += log(- dx * t1 + px * fiF + dx * t2 + (1 - px) * g / G);
//         nllh += ((R_pow(((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))),5)*(70*R_pow(((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))),4)-315*R_pow(((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))),3)+540*R_pow(((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))),2)-420*((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))))+126)))*(R_pow(((y-qalpha)*((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))/exp(lsbeta)+(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5)))),-1/(1.5/(1+exp(-txi))-.5)))+(1-((R_pow(((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))),5)*(70*R_pow(((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))),4)-315*R_pow(((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))),3)+540*R_pow(((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))),2)-420*((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))))+126))))*exp(-((y-((qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))*(log(-log(alpha))-log(-log(psuba)))/(log(-log(psuba))-log(-log(psubb)))))*(log(-log(hbeta))-log(-log(1-hbeta)))/(((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))*(log(-log(hbeta))-log(-log(1-hbeta)))/(log(-log(psuba))-log(-log(psubb))))-log(-log(alpha))))-log(-630*(R_pow(((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))),4))*(R_pow((1-((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))))),4))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))*(R_pow(((y-qalpha)*((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))/exp(lsbeta)+(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5)))),-1/(1.5/(1+exp(-txi))-.5)))+((R_pow(((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))),5)*(70*R_pow(((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))),4)-315*R_pow(((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))),3)+540*R_pow(((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))),2)-420*((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))))+126)))*(R_pow(((y-qalpha)*((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))/exp(lsbeta)+(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5)))),-(1+1/(1.5/(1+exp(-txi))-.5))))*(((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))/exp(lsbeta))/(1.5/(1+exp(-txi))-.5)+630*(R_pow(((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))),4))*(R_pow((1-((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))))),4))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))*exp(-((y-((qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))*(log(-log(alpha))-log(-log(psuba)))/(log(-log(psuba))-log(-log(psubb)))))*(log(-log(hbeta))-log(-log(1-hbeta)))/(((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))*(log(-log(hbeta))-log(-log(1-hbeta)))/(log(-log(psuba))-log(-log(psubb))))-log(-log(alpha))))+(1-((R_pow(((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))),5)*(70*R_pow(((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))),4)-315*R_pow(((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))),3)+540*R_pow(((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))),2)-420*((y-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))/((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))))+126))))*exp(-((y-((qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))*(log(-log(alpha))-log(-log(psuba)))/(log(-log(psuba))-log(-log(psubb)))))*(log(-log(hbeta))-log(-log(1-hbeta)))/(((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))*(log(-log(hbeta))-log(-log(1-hbeta)))/(log(-log(psuba))-log(-log(psubb))))-log(-log(alpha))))*((log(-log(hbeta))-log(-log(1-hbeta)))/(((qalpha+((R_pow(-log(psubb),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5)))))-(qalpha+((R_pow(-log(psuba),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(alpha),-(1.5/(1+exp(-txi))-.5))))*exp(lsbeta)/((R_pow(-log(1-hbeta),-(1.5/(1+exp(-txi))-.5)))-(R_pow(-log(hbeta),-(1.5/(1+exp(-txi))-.5))))))*(log(-log(hbeta))-log(-log(1-hbeta)))/(log(-log(psuba))-log(-log(psubb))))));
//       }
//     }
//     
//   }
//   
//   return(nllh);
//   
// }

// //' @rdname bgevd0
// [[Rcpp::export]]
arma::mat bgevd12(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec, arma::uvec dupid, int dcate, arma::vec other)
{
  
  arma::vec qavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lsbvec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = yvec.size();
  arma::mat out = arma::mat(nobs, 9, arma::fill::zeros);
  
  double psuba = other[0];
  double psubb = other[1];
  double alpha = other[2];
  double beta = other[3];
  double hbeta = 0.5 * beta;
  
  if (dcate == 1) {
    qavec = qavec.elem(dupid);
    lsbvec = lsbvec.elem(dupid);
    txivec = txivec.elem(dupid);
  }
  
  double y, qalpha, tqalpha, lsbeta, txi, xi, sbeta, tsbeta, iFa, iFb;
  
  // double ee16, ee27, ee29, ee51, ee52, ee53, ee65, ee66;
  // double ee94;
  // double ee106;
  // 
  // double ee34, ee60, ee77, ee79, ee86, ee91, ee97;
  // double ee100, ee105, ee118, ee119, ee121, ee128, ee129;
  // 
  // double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
  // double ee10, ee11, ee12, ee13, ee14, ee15, ee17, ee18, ee19;
  // double ee20, ee21, ee22, ee23, ee24, ee25, ee26, ee28;
  // double ee30, ee31, ee32, ee33, ee35, ee36, ee37, ee38, ee39;
  // double ee40, ee41, ee42, ee43, ee44, ee45, ee46, ee47, ee48, ee49;
  // double ee50, ee54, ee55, ee56, ee57, ee58, ee59;
  // double ee61, ee62, ee64, ee67, ee68, ee69;
  // double ee70, ee71, ee72, ee74, ee75, ee78;
  // double ee80, ee81, ee82, ee83, ee84, ee85, ee88, ee89;
  // double ee90, ee92, ee93, ee95, ee99;
  // double ee101, ee102, ee103, ee104, ee108, ee109;
  // double ee110, ee111, ee112, ee113, ee114, ee115, ee116, ee117;
  // double ee122, ee123, ee124, ee125, ee126, ee127;
  // double ee130, ee131, ee137;
  // double ee140, ee141, ee144, ee145, ee147;
  // double ee150, ee151, ee152, ee153, ee154, ee157, ee159;
  // double ee160, ee162, ee163, ee165, ee168;
  // double ee171, ee173, ee175, ee177, ee178, ee179;
  // double ee181, ee184, ee185, ee186, ee187, ee188, ee189;
  // double ee191, ee192, ee193, ee194, ee198, ee199;
  // double ee200, ee201;
  // double ee212, ee213, ee214;
  // double ee225, ee227;
  // double ee234;
  // double ee242, ee243, ee247;
  // double ee250, ee251, ee252, ee255, ee258;
  // double ee261, ee265, ee268, ee269;
  // double ee270, ee271, ee272, ee273, ee274, ee275, ee276;
  // 
  double ee3, ee5, ee6, ee8, ee9;
  double ee12, ee13, ee15, ee17, ee18, ee19;
  double ee21, ee22, ee23, ee25, ee26, ee27, ee29;
  double ee30, ee32, ee33, ee37, ee38, ee39;
  double ee41, ee42, ee43, ee45, ee46;
  double ee50, ee51, ee52, ee53, ee54, ee56, ee57, ee58;
  double ee60, ee61, ee62, ee63, ee64, ee65, ee66;
  double ee71, ee72, ee74, ee75, ee79;
  double ee80, ee85, ee86, ee87, ee88;
  double ee90, ee92, ee93, ee96, ee98, ee99;
  double ee102, ee107, ee109;
  double ee112, ee114, ee115, ee117, ee118;
  double ee122, ee124, ee126, ee127, ee129;
  double ee131, ee132, ee133, ee134, ee135, ee136, ee137;
  double ee141, ee143, ee146, ee148;
  double ee151, ee153, ee155, ee156, ee157, ee158;
  double ee160, ee161, ee162, ee165, ee167, ee169;
  double ee173, ee176, ee177, ee178;
  double ee187;
  double ee194, ee197;
  double ee200, ee202, ee206, ee209;
  double ee214, ee216;
  double ee221, ee227;
  double ee244, ee248;
  double ee260, ee266;
  double ee279;
  double ee283, ee284, ee285, ee286, ee288, ee289;
  double ee291, ee292, ee293, ee295, ee297;
  double ee300, ee304, ee309;
  double ee310;
  double ee323, ee325, ee327, ee329;
  double ee334, ee338;
  double ee348;
  double ee351, ee356;
  double ee362, ee365, ee367, ee368, ee369;
  double ee370, ee371, ee376;
  double ee380, ee382, ee388;
  double ee395;
  double ee402, ee404, ee408;
  double ee424, ee427, ee428;
  double ee430, ee431, ee438;
  double ee445, ee447;
  double ee452;
  double ee465, ee468;
  double ee471;
  double ee487, ee489;
  double ee494, ee495, ee496, ee497, ee498, ee499;
  double ee500, ee501, ee502, ee503, ee504, ee506, ee508, ee509;
  double ee510, ee511, ee513, ee514, ee515, ee516, ee518, ee519;
  double ee520, ee522, ee523, ee525, ee527, ee528, ee529;
  double ee533, ee534, ee535, ee537, ee538;
  double ee541, ee543, ee544, ee546, ee548;
  double ee551, ee555, ee558;
  double ee563;
  double ee593, ee596, ee599;
  double ee605, ee607, ee609;
  double ee610, ee611, ee613, ee615, ee617, ee618;
  double ee627;
  double ee631, ee634, ee636, ee637, ee638, ee639;
  double ee641, ee643, ee644, ee646;
  double ee651, ee653, ee659;
  double ee661, ee666;
  double ee682, ee686, ee688;
  double ee691;
  double ee701;
  double ee712, ee714, ee717;
  double ee722, ee723;
  double ee732;
  double ee742;
  double ee751, ee757;
  double ee760;
  double ee776, ee778;
  double ee784;
  double ee791, ee797, ee799;
  double ee802, ee808;
  double ee815, ee817;
  double ee825;
  double ee838;
  double ee841, ee843, ee846;
  double ee856, ee859;
  double ee861, ee869;
  double ee873;
  double ee887;
  double ee899;
  double ee906;
  double ee926;
  double ee947;
  double ee957;
  double ee997;
  double ee1001;
  double ee1020, ee1024;
  double ee1031, ee1033, ee1037;
  double ee1115;
  double ee1131;
  double ee1141;
  double ee1152, ee1156;
  double ee1181, ee1186, ee1188;
  double ee1191, ee1193, ee1199;
  double ee1215, ee1216;
  double ee1223, ee1226;
  double ee1232, ee1238;
  double ee1249;
  double ee1259;
  double ee1261, ee1269;
  double ee1271, ee1276, ee1278;
  double ee1293, ee1297;
  double ee1311;
  double ee1322, ee1329;
  double ee1341, ee1343, ee1347;
  double ee1356;
  double ee1362, ee1368;
  double ee1372;
  double ee1402;
  double ee1413;
  
  double ee2, ee7, ee10, ee11, ee14, ee16, ee20, ee24, ee28, ee31, ee34, ee35, ee36;
  double ee44, ee47, ee48, ee59, ee77, ee78, ee82, ee84, ee89, ee91, ee94, ee97;
  double ee100, ee101, ee104, ee105, ee106, ee116, ee119, ee121, ee125, ee128, ee130;
  
  // 
  // double ee2, ee3, ee5, ee6, ee7, ee8, ee9;
  // double ee10, ee11, ee12, ee13, ee14, ee15, ee16, ee17, ee18, ee19;
  // double ee21, ee22, ee23, ee24, ee25, ee26, ee27, ee28, ee29;
  // double ee30, ee31, ee32, ee34, ee35, ee36, ee37, ee38;
  // double ee42, ee43, ee44, ee45, ee46, ee47, ee48;
  // double ee50, ee52, ee53, ee57, ee58, ee59;
  // double ee61, ee62, ee64, ee65, ee66;
  // double ee71, ee72, ee77, ee78, ee79;
  // double ee80, ee82, ee84, ee86, ee88, ee89;
  // double ee90, ee91, ee92, ee94, ee97;
  // double ee100, ee101, ee102, ee104, ee105, ee106;
  // double ee112, ee114, ee116, ee117, ee118, ee119;
  // double ee121, ee124, ee125, ee128, ee129;
  // double ee130, ee131;
  // 
  // double ee20, ee33, ee51, ee56, ee60;
  
  double ee4;

  for (int j=0; j < nobs; j++) {
    
    y = yvec[j];
    qalpha = qavec[j];
    lsbeta = lsbvec[j];
    txi = txivec[j];
    sbeta = exp(lsbeta);
    xi = 1.5 / (1.0 + exp(-txi)) - 0.5;
    
    iFa = iF(psuba, qalpha, sbeta, xi, alpha, beta);
    iFb = iF(psubb, qalpha, sbeta, xi, alpha, beta);
    
    tqalpha = iFa - (iFb - iFa) * (ell2(alpha) - ell2(psuba)) / (ell2(psuba) - ell2(psubb));
    tsbeta = (iFb - iFa) * (ell2(hbeta) - ell2(1.0 - hbeta)) / (ell2(psuba) - ell2(psubb));
    
    if (y < iFa) { // Gumbel
      
      ee2 = exp(-txi);
      ee3 = 1 + ee2;
      ee5 = 1.5/ee3 - 0.5;
      ee6 = -log(psuba);
      ee7 = -log(alpha);
      ee8 = R_pow(ee7, ee5);
      ee9 = -log(psubb);
      ee10 = 1/ee8;
      ee11 = -log(1 - hbeta);
      ee12 = -log(hbeta);
      ee13 = R_pow(ee6, ee5);
      ee14 = log(ee6);
      ee15 = R_pow(ee11, ee5);
      ee16 = R_pow(ee12, ee5);
      ee17 = 1/ee13;
      ee18 = R_pow(ee9, ee5);
      ee19 = log(ee9);
      ee21 = 1/ee15 - 1/ee16;
      ee22 = 1/ee18;
      ee23 = ee10 + ee17;
      ee24 = ee10 + ee22;
      ee25 = ee14 - ee19;
      ee26 = ee24 - ee23;
      ee27 = log(ee11);
      ee28 = log(ee12);
      ee29 = exp(lsbeta);
      ee30 = log(ee7);
      ee31 = ee28 - ee27;
      ee32 = ee26 * ee29;
      ee34 = 1.5 * (ee27/ee15) - 1.5 * (ee28/ee16);
      ee35 = ee21 * ee25;
      ee36 = ee30 - ee14;
      ee37 = R_pow(ee3, 2);
      ee38 = ee14/ee13;
      ee42 = ee17 - (ee26 * ee36/ee25 + ee10);
      ee43 = R_pow(ee32 * ee31/ee35, 2);
      ee44 = R_pow(ee31, 2);
      ee45 = 1.5 * ee38;
      ee46 = 1.5 - 3 * (ee2/ee3);
      ee47 = ee23 - ee24;
      ee48 = ee19/ee18;
      ee50 = ee42 * ee29/ee21;
      ee52 = y - (ee50 + qalpha);
      ee53 = ee43 * ee21;
      ee57 = ee34 * ee47/ee21 + 1.5 * ee48 - ee45;
      ee58 = ee53 * ee25;
      ee59 = ee17 - ee10;
      ee61 = ee34 * ee59/ee21;
      ee62 = ee30/ee8;
      ee64 = ee32 * ee44;
      ee65 = 3 * ee38;
      ee66 = exp(-(ee35 * ee52/ee32 - ee30));
      ee71 = ee45 - (ee57 * ee36/ee25 + ee61 + 1.5 * ee62);
      ee72 = R_pow(ee29, 2);
      ee77 = (ee46/ee15 + 2.25 * (ee2 * ee27/(ee15 * ee37))) *  ee27 - (ee46/ee16 + 2.25 * (ee2 * ee28/(ee16 * ee37))) *  ee28;
      ee78 = (ee46/ee13 + 2.25 * (ee2 * ee14/(ee13 * ee37))) *  ee14;
      ee79 = R_pow(ee21, 2);
      ee80 = R_pow(ee25, 2);
      ee82 = ee57 * ee29 * ee44;
      ee84 = ee71 * ee25/ee26;
      ee86 = ee64 * ee52/ee58;
      ee88 = ee77 * ee47;
      ee89 = (ee46/ee18 + 2.25 * (ee2 * ee19/(ee18 * ee37))) *  ee19;
      ee90 = (ee34 * (2 * (ee22 - ee10) - 2 * ee59)/ee21 + ee65 -  3 * ee48) * ee34;
      ee91 = ee43 * ee79;
      ee92 = R_pow(ee26, 2);
      ee94 = ee42 * ee25/ee26;
      ee97 = (ee88 + ee90 * ee2/ee37)/ee21 + ee89 - ee78;
      ee100 = ee82 * ee52/ee58 + ee84;
      ee101 = ee86 + ee94;
      ee102 = ee92 * ee72;
      ee104 = 2 * ee61 + 3 * ee62;
      ee105 = R_pow(ee35/ee32, 2);
      ee106 = ee91 * ee80;
      ee112 = ee97 * ee52;
      ee114 = (ee78 - (ee97 * ee36/ee25 + (ee77 * ee59 + ee34 *  (ee65 - ee104) * ee2/ee37)/ee21 + (ee46/ee8 + 2.25 *  (ee2 * ee30/(ee8 * ee37))) * ee30)) * ee25/ee26;
      ee116 = (ee57 * (2 * ee86 + 2 * ee14 - 2 * ee30)/ee25 +  ee65 - ee104) * ee57 * ee2;
      ee117 = ee82/ee58;
      ee118 = R_pow(ee57, 2);
      ee119 = ee53 * ee80;
      ee121 = (ee25/ee26 - ee26 * ee72 * ee44/(ee91 * ee25)) *  ee42;
      ee124 = ee71 * ee26 * ee29/ee21;
      ee125 = ee64/ee58;
      ee128 = (2 * (ee102 * ee44/ee106) - 1) * ee52;
      ee129 = 1/ee105;
      ee130 = 2 * ee42;
      ee131 = ee66 - 1;
      
      out(j, 0) = ee21 * ee131 * ee25/ee32;
      out(j, 1) = (((1 + ee14 - ee30) * ee26/ee25 + ee17 - ee10) * ee29/ee21 +
        qalpha - y) * ee26 * ee29 * ee44/ee58 + ee101 * ee66 -
        ee94;
      out(j, 2) = -((ee100 * ee131 + ee57 * ee26 * ee72 * ee44/ee106) * ee2/ee37);
      out(j, 3) = 1 * ee79 * ee66 * ee80/ee102;
      out(j, 4) = (ee101 * ee21 * ee25/ee32 - ee125) * ee66 + ee125;
      out(j, 5) = ((ee117 - ee100 * ee21 * ee25/ee32) * ee66 - ee117) * ee2/ee37;
      out(j, 6) = ((((ee26 * (2 * ee52 - 2 * (ee32/ee35)) * ee29/ee35 +
        ee129) * ee44/ee43 + 1) * ee26/ee25 + ee130) * ee29/ee21 +
        qalpha - y) * ee26 * ee29 * ee44/ee58 + (R_pow(ee101, 2) +
        ee121 + ee64 * (y - ((2 * (ee92 * ee29 * ee44 * ee52/ee119) +
        ee130) * ee29/ee21 + qalpha))/ee58) * ee66 - ee121;
      out(j, 7) = (((ee57 * (ee50 + ee128) + ee124) * ee29 * ee44/ee58 -
        (ee100 * ee101 + ee84)) * ee66 + ee84 - (((((ee129 -
        2 * (ee102/(ee79 * ee80))) * ee44/ee43 + 1 + ee14 - ee30) * ee26/ee25 +
        ee17 - ee10) * ee29/ee21 + ee128) * ee57 + ee124) * ee29 * ee44/ee58) * ee2/ee37;
      out(j, 8) =  - ((((ee116 * ee29/(ee37 * ee21) - ee112) * ee29 * ee44/ee58 -
        (ee114 + R_pow(ee100, 2) * ee2/ee37)) * ee66 +
        ee114 - (((((ee88 + (ee90 - 2 * (ee118 * ee26 * ee72 * ee44/ee119)) * ee2/ee37)/ee21 +
        ee89 - ee78) * ee26 + ee118 * ee2 * ee44/(ee105 * ee43 * ee37))/ee25 +
        ee116/ee37) * ee29/ee21 -
        ee112) * ee29 * ee44/ee58) * ee2/ee37);
      
    } else {
      
      if (y > iFb) { // GEV
        
        ee2 = exp(-txi);
        ee3 = 1 + ee2;
        ee5 = 1.5/ee3 - 0.5;
        ee6 = -log(1 - hbeta);
        ee7 = -log(hbeta);
        ee8 = R_pow(ee6, ee5);
        ee9 = R_pow(ee7, ee5);
        ee10 = -log(alpha);
        ee11 = exp(lsbeta);
        ee13 = 1/ee8 - 1/ee9;
        ee14 = y - qalpha;
        ee15 = R_pow(ee10, ee5);
        ee16 = ee13 * ee14;
        ee18 = ee16/ee11 + 1/ee15;
        ee19 = 1/ee5;
        ee20 = log(ee6);
        ee21 = log(ee7);
        ee22 = 1 + ee19;
        ee23 = 1.5 * (ee20/ee8);
        ee24 = 1.5 * (ee21/ee9);
        ee25 = R_pow(ee3, 2);
        ee26 = log(ee10);
        ee27 = ee23 - ee24;
        ee30 = ee27 * ee14/ee11 + 1.5 * (ee26/ee15);
        ee31 = R_pow(ee18, ee22);
        ee32 = 1.5 - 3 * (ee2/ee3);
        ee33 = log(ee18);
        ee34 = R_pow(ee18, ee19);
        ee35 = R_pow(ee5, 2);
        ee36 = R_pow(ee18, ee19 + 2);
        ee37 = ee30 * ee22;
        ee38 = ee18 * ee11;
        ee42 = (ee32/ee8 + 2.25 * (ee2 * ee20/(ee8 * ee25))) * ee20 -  (ee32/ee9 + 2.25 * (ee2 * ee21/(ee9 * ee25))) * ee21;
        ee44 = ee37/ee36 + 1.5 * (ee33/(ee31 * ee35));
        ee45 = ee30 * ee2;
        ee46 = R_pow(ee13/ee11, 2);
        ee47 = R_pow(ee11, 2);
        ee51 = ee42 * ee14/ee11 + (ee32/ee15 + 2.25 * (ee2 * ee26/(ee15 *  ee25))) * ee26;
        ee53 = ((1 - 1/ee34)/ee5 + 1) * ee13;
        ee56 = ((1.5 * (ee13/ee5) + ee23 - ee24)/ee31 - ee44 * ee13)/ee5 -  (ee22 * (ee23 - (ee30 * ee13/ee18 + ee24)) + 1.5 * (ee13/ee35))/ee18;
        ee57 = ee30/ee31;
        ee58 = ee18 * ee25;
        ee60 = ee34 * ee5;
        ee61 = ((4.5/(ee3 * ee5) - 3) * ee2/ee3 + 1.5) * ee33;
        ee65 = ee16/ee38;
        ee66 = (1/ee31 - ee22 * ee13 * ee14/(ee36 * ee11))/ee5;
        
        out(j, 0) = -(ee53/ee38);
        out(j, 1) = -(ee53 * ee14/ee38 - 1);
        out(j, 2) = ((ee57 + (1.5/ee34 - 1.5) * ee33/ee5 + 1.5)/ee5 +
          ee27/ee13 - ee37/ee18) * ee2/ee25;
        out(j, 3) =  - ((1 - 1/ee60) * ee22 * R_pow(ee13, 2)/(R_pow(ee18, 2) * ee47));
        out(j, 4) = -((ee66 - (1 - ee65) * ee22/ee18) * ee13/ee11);
        out(j, 5) = -(ee56 * ee2/(ee25 * ee11));
        out(j, 6) = -((((ee65 - 1) * ee22/ee18 + ee66) * ee14 - ee13/(ee46 * ee11)) * ee13/ee11 + 1);
        out(j, 7) = -((ee56 * ee14/ee11 + ee27 * (1/ee13 - ee13/(ee46 * ee47))) * ee2/ee25);
        out(j, 8) = ((((ee44 * ee30 - 2.25/ee5)/ee3 + 3) * ee2/ee3 +
          (ee61 + (1.5 * ((ee57 + 1.5 * (ee33/ee60)) * ee33/ee5) + 1.5 * (ee30/ee18)) * ee2/ee25 -
          (ee61 + 1.5 * (ee45/ee58))/ee34)/ee5 -
          ((ee51 + 1.5 * (ee45/(ee25 * ee5)))/ee31 + 1.5))/ee5 +
          ((ee51 - R_pow(ee30, 2) * ee2/ee58) * ee22 + 1.5 * (ee45/(ee25 * ee35)))/ee18 +
          R_pow(ee27, 2) * ee2/(ee46 * ee25 * ee47) -
          ee42/ee13) * ee2/ee25;
        
      } else {
        
        ee3 = -log(1 - hbeta);
        ee5 = exp(-txi);
        ee6 = 1 + ee5;
        ee8 = 1.5/ee6 - 0.5;
        ee9 = R_pow(ee3, ee8);
        ee12 = -log(hbeta);
        ee13 = R_pow(ee12, ee8);
        ee15 = 1/ee9 - 1/ee13;
        ee17 = -log(psuba);
        ee18 = R_pow(ee17, ee8);
        ee19 = 1/ee18;
        ee21 = -log(alpha);
        ee22 = R_pow(ee21, ee8);
        ee23 = 1/ee22;
        ee25 = exp(lsbeta);
        ee26 = (ee19 - ee23) * ee25;
        ee27 = ee26/ee15;
        ee29 = y - (ee27 + qalpha);
        ee30 = ee15 * ee29;
        ee32 = -log(psubb);
        ee33 = R_pow(ee32, ee8);
        ee37 = ee23 + 1/ee33 - (ee23 + ee19);
        ee38 = ee37 * ee25;
        ee39 = ee30/ee38;
        ee41 = 70 * ee39 - 315;
        ee42 = ee15 * ee41;
        ee43 = ee42 * ee29;
        ee45 = ee43/ee38 + 540;
        ee46 = R_pow(ee39, 2);
        ee50 = ee45 * ee46 + 126 - 420 * ee39;
        ee51 = R_pow(ee39, 5);
        ee52 = ee50 * ee51;
        ee53 = y - qalpha;
        ee54 = ee15 * ee53;
        ee56 = ee54/ee25 + ee23;
        ee57 = 1/ee8;
        ee58 = R_pow(ee56, ee57);
        ee60 = 1 - ee52;
        ee61 = log(ee17);
        ee62 = log(ee32);
        ee63 = ee61 - ee62;
        ee64 = ee15 * ee63;
        ee65 = log(ee21);
        ee66 = ee65 - ee61;
        ee71 = (ee19 - (ee37 * ee66/ee63 + ee23)) * ee25;
        ee72 = ee71/ee15;
        ee74 = y - (ee72 + qalpha);
        ee75 = ee64 * ee74;
        ee79 = exp(-(ee75/ee38 - ee65));
        ee80 = ee60 * ee79;
        ee85 = ee50 * ee15;
        ee86 = ee85 * ee29;
        ee87 = 1 + ee57;
        ee88 = R_pow(ee56, ee87);
        ee90 = ee88 * ee8 * ee25;
        ee92 = 1 - ee39;
        ee93 = R_pow(ee92, 4);
        ee96 = 630 * ee79 - 630/ee58;
        ee98 = ee86/ee90 + ee93 * ee96;
        ee99 = R_pow(ee39, 4);
        ee102 = ee98 * ee99 + ee80 * ee63;
        ee107 = ee15/ee38;
        ee109 = 5 * (ee107 * ee99);
        ee112 = 2 * (ee107 * ee39);
        ee114 = 70 * ee107;
        ee115 = ee15 * ee114;
        ee117 = ee42 + ee115 * ee29;
        ee118 = ee117/ee38;
        ee122 = ee45 * ee112 + ee118 * ee46 - 420 * ee107;
        ee124 = ee50 * ee109 + ee122 * ee51;
        ee126 = ee64/ee38;
        ee127 = ee79 * ee126;
        ee129 = ee124 * ee79 + ee60 * ee127;
        ee131 = ee57 - 1;
        ee132 = R_pow(ee56, ee131);
        ee133 = ee15/ee25;
        ee134 = ee57 * ee133;
        ee135 = ee132 * ee134;
        ee136 = ee52 * ee135;
        ee137 = R_pow(ee58, 2);
        ee141 = R_pow(ee92, 3);
        ee143 = 4 * (ee107 * ee141);
        ee146 = 630 * ee135;
        ee148 = 630 * ee127 - ee146/ee137;
        ee151 = ee122 * ee15;
        ee153 = ee85 + ee151 * ee29;
        ee155 = ee87 - 1;
        ee156 = R_pow(ee56, ee155);
        ee157 = ee87 * ee133;
        ee158 = ee156 * ee157;
        ee160 = ee158 * ee8 * ee25;
        ee161 = ee86 * ee160;
        ee162 = R_pow(ee90, 2);
        ee165 = ee143 * ee96 + ee93 * ee148 - (ee153/ee90 - ee161/ee162);
        ee167 = R_pow(ee39, 3);
        ee169 = 4 * (ee107 * ee167);
        ee173 = ee165 * ee99 - ee98 * ee169 + ee129 * ee63;
        ee176 = ee124 * ee127;
        ee177 = ee122 * ee109;
        ee178 = ee118 * ee112;
        ee187 = ee178 + (ee115 + ee115)/ee38 * ee46 + (ee45 * (2 *  (ee107 * ee107)) + ee178);
        ee194 = ee177 + ee187 * ee51 + (ee50 * (5 * (ee107 * ee169)) +  ee177);
        ee197 = ee127 * ee126;
        ee200 = ee176 - ee194 * ee79 + (ee176 + ee60 * ee197);
        ee202 = ee124 * ee135;
        ee206 = R_pow(ee56, ee131 - 1);
        ee209 = ee206 * (ee131 * ee133) * ee134;
        ee214 = 2 * (ee135 * ee58);
        ee216 = R_pow(ee137, 2);
        ee221 = R_pow(ee92, 2);
        ee227 = ee143 * ee148;
        ee244 = ee153 * ee160;
        ee248 = R_pow(ee56, ee155 - 1);
        ee260 = R_pow(ee162, 2);
        ee266 = ee165 * ee169;
        ee279 = R_pow(ee102, 2);
        ee283 = ee64 * ee72;
        ee284 = ee283/ee38;
        ee285 = ee75 * ee38;
        ee286 = R_pow(ee38, 2);
        ee288 = ee284 + ee285/ee286;
        ee289 = ee79 * ee288;
        ee291 = ee15 * ee27;
        ee292 = ee291/ee38;
        ee293 = ee30 * ee38;
        ee295 = ee292 + ee293/ee286;
        ee297 = 5 * (ee295 * ee99);
        ee300 = 2 * (ee295 * ee39);
        ee304 = ee15 * ee38/ee286;
        ee309 = 70 * ee295;
        ee310 = ee15 * ee309;
        ee323 = ee42 * ee27;
        ee325 = ee323 + ee310 * ee29;
        ee327 = ee43 * ee38;
        ee329 = ee325/ee38 + ee327/ee286;
        ee334 = ee118 * ee300 + ((ee115 * ee27 + ee15 * (70 * ee304) *  ee29 + ee310)/ee38 + ee117 * ee38/ee286) * ee46 + (ee45 *  (2 * (ee107 * ee295 + ee304 * ee39)) + ee329 * ee112) -  420 * ee304;
        ee338 = 4 * (ee295 * ee167);
        ee348 = ee45 * ee300 + ee329 * ee46 - 420 * ee295;
        ee351 = ee122 * ee297 + ee334 * ee51 + (ee50 * (5 * (ee107 *  ee338 + ee304 * ee99)) + ee348 * ee109);
        ee356 = ee50 * ee297 + ee348 * ee51;
        ee362 = ee289 * ee126 - ee79 * (ee64 * ee38/ee286);
        ee365 = ee124 * ee289 - ee351 * ee79 + (ee356 * ee127 +  ee60 * ee362);
        ee367 = ee54 * ee25;
        ee368 = R_pow(ee25, 2);
        ee369 = ee367/ee368;
        ee370 = ee57 * ee369;
        ee371 = ee132 * ee370;
        ee376 = ee15 * ee25/ee368;
        ee380 = ee206 * (ee131 * ee369);
        ee382 = ee132 * (ee57 * ee376) + ee380 * ee134;
        ee388 = 2 * (ee371 * ee58);
        ee395 = 3 * (ee295 * ee221);
        ee402 = 630 * ee371;
        ee404 = 630 * ee289 - ee402/ee137;
        ee408 = 4 * (ee295 * ee141);
        ee424 = ee348 * ee15;
        ee427 = ee87 * ee369;
        ee428 = ee156 * ee427;
        ee430 = ee428 * ee8 * ee25;
        ee431 = ee90 - ee430;
        ee438 = ee248 * (ee155 * ee369);
        ee445 = ee85 * ee27;
        ee447 = ee445 + ee424 * ee29;
        ee452 = 2 * (ee431 * ee90);
        ee465 = ee86 * ee431;
        ee468 = ee408 * ee96 + ee93 * ee404 - (ee447/ee90 + ee465/ee162);
        ee471 = 3 * (ee295 * ee46);
        ee487 = ee356 * ee79 + ee60 * ee289;
        ee489 = ee468 * ee99 - ee98 * ee338 + ee487 * ee63;
        ee494 = 1.5 * ee5;
        ee495 = R_pow(ee6, 2);
        ee496 = ee494/ee495;
        ee497 = ee61 * ee496;
        ee498 = ee18 * ee497;
        ee499 = R_pow(ee18, 2);
        ee500 = ee498/ee499;
        ee501 = ee65 * ee496;
        ee502 = ee22 * ee501;
        ee503 = R_pow(ee22, 2);
        ee504 = ee502/ee503;
        ee506 = (ee500 - ee504) * ee25;
        ee508 = log(ee3);
        ee509 = ee508 * ee496;
        ee510 = ee9 * ee509;
        ee511 = R_pow(ee9, 2);
        ee513 = log(ee12);
        ee514 = ee513 * ee496;
        ee515 = ee13 * ee514;
        ee516 = R_pow(ee13, 2);
        ee518 = ee510/ee511 - ee515/ee516;
        ee519 = ee26 * ee518;
        ee520 = R_pow(ee15, 2);
        ee522 = ee506/ee15 - ee519/ee520;
        ee523 = ee15 * ee522;
        ee525 = ee523 - ee518 * ee29;
        ee527 = ee62 * ee496;
        ee528 = ee33 * ee527;
        ee529 = R_pow(ee33, 2);
        ee533 = ee528/ee529 + ee504 - (ee500 + ee504);
        ee534 = ee533 * ee25;
        ee535 = ee30 * ee534;
        ee537 = ee525/ee38 + ee535/ee286;
        ee538 = 70 * ee537;
        ee541 = ee15 * ee538 - ee518 * ee41;
        ee543 = ee42 * ee522;
        ee544 = ee541 * ee29 + ee543;
        ee546 = ee43 * ee534;
        ee548 = ee544/ee38 + ee546/ee286;
        ee551 = 2 * (ee537 * ee39);
        ee555 = ee548 * ee46 + ee45 * ee551 - 420 * ee537;
        ee558 = 4 * (ee537 * ee167);
        ee563 = ee518/ee38 - ee15 * ee534/ee286;
        ee593 = ee548 * ee112 + ee45 * (2 * (ee107 * ee537 - ee563 *  ee39)) + (((ee541 + (ee115 * ee522 - (ee15 * (70 * ee563) +  ee518 * ee114) * ee29))/ee38 + ee117 * ee534/ee286) *  ee46 + ee118 * ee551) + 420 * ee563;
        ee596 = 5 * (ee537 * ee99);
        ee599 = ee555 * ee109 + ee50 * (5 * (ee107 * ee558 - ee563 *  ee99)) + (ee593 * ee51 + ee122 * ee596);
        ee605 = (ee500 - (ee504 + ee533 * ee66/ee63)) * ee25;
        ee607 = ee71 * ee518;
        ee609 = ee605/ee15 - ee607/ee520;
        ee610 = ee64 * ee609;
        ee611 = ee518 * ee63;
        ee613 = ee610 - ee611 * ee74;
        ee615 = ee75 * ee534;
        ee617 = ee613/ee38 + ee615/ee286;
        ee618 = ee79 * ee617;
        ee627 = ee79 * (ee611/ee38 - ee64 * ee534/ee286) + ee618 *  ee126;
        ee631 = ee555 * ee51 + ee50 * ee596;
        ee634 = ee599 * ee79 - ee124 * ee618 - (ee60 * ee627 + ee631 *  ee127);
        ee636 = log(ee56);
        ee637 = R_pow(ee8, 2);
        ee638 = ee496/ee637;
        ee639 = ee636 * ee638;
        ee641 = ee518 * ee53;
        ee643 = ee504 + ee641/ee25;
        ee644 = ee57 * ee643;
        ee646 = ee58 * ee639 + ee132 * ee644;
        ee651 = ee518/ee25;
        ee653 = ee638 * ee133;
        ee659 = ee132 * ee639 + ee206 * (ee131 * ee643);
        ee661 = ee132 * (ee57 * ee651 + ee653) + ee659 * ee134;
        ee666 = 2 * (ee646 * ee58);
        ee682 = 4 * (ee537 * ee141);
        ee686 = 630 * ee646;
        ee688 = 630 * ee618 + ee686/ee137;
        ee691 = 3 * (ee537 * ee221);
        ee701 = ee555 * ee15 - ee50 * ee518;
        ee712 = ee87 * ee643;
        ee714 = ee88 * ee639 + ee156 * ee712;
        ee717 = (ee88 * ee496 - ee714 * ee8) * ee25;
        ee722 = ee85 * ee522;
        ee723 = ee701 * ee29 + ee722;
        ee732 = ee156 * ee639 + ee248 * (ee155 * ee643);
        ee742 = 2 * (ee717 * ee90);
        ee751 = ee86 * ee717;
        ee757 = ee723/ee90 - ee751/ee162 - (ee93 * ee688 + ee682 *  ee96);
        ee760 = 3 * (ee537 * ee46);
        ee776 = ee60 * ee618 + ee631 * ee79;
        ee778 = ee757 * ee99 + ee98 * ee558 - ee776 * ee63;
        ee784 = ee52 * ee371;
        ee791 = ee291 * ee38;
        ee797 = 2 * (ee38 * ee38);
        ee799 = R_pow(ee286, 2);
        ee802 = ee292 - ee791/ee286 + ((ee293 - ee791)/ee286 - ee293 *  ee797/ee799);
        ee808 = ee348 * ee297;
        ee815 = ee329 * ee300;
        ee817 = ee310 * ee27;
        ee825 = ee325 * ee38;
        ee838 = ee45 * (2 * (ee802 * ee39 - ee295 * ee295)) - ee815 +  (((ee323 - ee817 + (ee15 * (70 * ee802) * ee29 - ee817))/ee38 -  ee825/ee286 + ((ee327 - ee825)/ee286 - ee327 * ee797/ee799)) *  ee46 - ee815) - 420 * ee802;
        ee841 = ee50 * (5 * (ee802 * ee99 - ee295 * ee338)) - ee808 +  (ee838 * ee51 - ee808);
        ee843 = ee356 * ee289;
        ee846 = ee283 * ee38;
        ee856 = ee289 * ee288 + ee79 * (ee284 - ee846/ee286 + ((ee285 -  ee846)/ee286 - ee285 * ee797/ee799));
        ee859 = ee841 * ee79 + ee843 + (ee843 + ee60 * ee856);
        ee861 = ee356 * ee371;
        ee869 = ee369 - ee367 * (2 * (ee25 * ee25))/R_pow(ee368, 2);
        ee873 = ee132 * (ee57 * ee869) - ee380 * ee370;
        ee887 = ee408 * ee404;
        ee899 = ee424 * ee27;
        ee906 = ee447 * ee431;
        ee926 = ee468 * ee338;
        ee947 = 2 * (ee534 * ee38);
        ee957 = (ee525 * ee38 - ee535)/ee286 + ee293 * ee947/ee799 -  ((ee523 + ee518 * ee27)/ee38 - ee291 * ee534/ee286);
        ee997 = ee548 * ee300 + ee45 * (2 * (ee957 * ee39 + ee295 *  ee537)) + (((ee541 * ee27 - ee543 + ((ee15 * (70 * ee957) -  ee518 * ee309) * ee29 + ee310 * ee522))/ee38 + ee325 *  ee534/ee286 + ((ee544 * ee38 - ee546)/ee286 + ee327 *  ee947/ee799)) * ee46 + ee329 * ee551) - 420 * ee957;
        ee1001 = ee555 * ee297 + ee50 * (5 * (ee957 * ee99 + ee295 *  ee558)) + (ee997 * ee51 + ee348 * ee596);
        ee1020 = ee79 * ((ee613 * ee38 - ee615)/ee286 + ee285 *  ee947/ee799 - ((ee610 + ee611 * ee72)/ee38 - ee283 *  ee534/ee286)) - ee618 * ee288;
        ee1024 = ee1001 * ee79 - ee356 * ee618 + (ee60 * ee1020 -  ee631 * ee289);
        ee1031 = ee641 * ee25/ee368;
        ee1033 = ee638 * ee369;
        ee1037 = ee132 * (ee57 * ee1031 + ee1033) + ee659 * ee370;
        ee1115 = ee52 * ee646;
        ee1131 = ee496 - ee494 * (2 * (ee5 * ee6))/R_pow(ee495, 2);
        ee1141 = (ee498 * ee497 - ee18 * (ee61 * ee1131))/ee499 - ee498 * (2 * (ee498 * ee18))/R_pow(ee499, 2);
        ee1152 = (ee502 * ee501 - ee22 * (ee65 * ee1131))/ee503 - ee502 * (2 * (ee502 * ee22))/R_pow(ee503, 2);
        ee1156 = ee506 * ee518;
        ee1181 = (ee510 * ee509 - ee9 * (ee508 * ee1131))/ee511 - ee510 * (2 * (ee510 * ee9))/R_pow(ee511, 2) - ((ee515 * ee514 - ee13 * (ee513 * ee1131))/ee516 - ee515 * (2 * (ee515 * ee13))/R_pow(ee516, 2));
        ee1186 = 2 * (ee518 * ee15);
        ee1188 = R_pow(ee520, 2);
        ee1191 = (ee1141 - ee1152) * ee25/ee15 + ee1156/ee520 -  ((ee26 * ee1181 - ee1156)/ee520 + ee519 * ee1186/ee1188);
        ee1193 = ee518 * ee522;
        ee1199 = ee525 * ee534;
        ee1215 = (ee528 * ee527 - ee33 * (ee62 * ee1131))/ee529 - ee528 * (2 * (ee528 * ee33))/R_pow(ee529, 2) + ee1152 - (ee1141 + ee1152);
        ee1216 = ee1215 * ee25;
        ee1223 = (ee15 * ee1191 - ee1193 - (ee1181 * ee29 + ee1193))/ee38 +  ee1199/ee286 + ((ee1199 + ee30 * ee1216)/ee286 + ee535 *  ee947/ee799);
        ee1226 = ee518 * ee538;
        ee1232 = ee541 * ee522;
        ee1238 = ee544 * ee534;
        ee1249 = ee548 * ee551;
        ee1259 = (((ee15 * (70 * ee1223) - ee1226 - (ee1181 * ee41 +  ee1226)) * ee29 + ee1232 + (ee1232 + ee42 * ee1191))/ee38 +  ee1238/ee286 + ((ee1238 + ee43 * ee1216)/ee286 + ee546 *  ee947/ee799)) * ee46 + ee1249 + (ee1249 + ee45 * (2 *  (ee1223 * ee39 + ee537 * ee537))) - 420 * ee1223;
        ee1261 = ee555 * ee596;
        ee1269 = ee1259 * ee51 + ee1261 + (ee1261 + ee50 * (5 *  (ee1223 * ee99 + ee537 * ee558)));
        ee1271 = ee631 * ee646;
        ee1276 = ee1152 + ee1181 * ee53/ee25;
        ee1278 = ee638 * ee643;
        ee1293 = ee636 * (ee1131/ee637 + ee496 * (2 * (ee496 * ee8))/R_pow(ee637, 2)) + ee643/ee56 * ee638;
        ee1297 = ee132 * (ee57 * ee1276 - ee1278) - ee659 * ee644 -  (ee58 * ee1293 + ee646 * ee639);
        ee1311 = ee605 * ee518;
        ee1322 = ee611 * ee609;
        ee1329 = ee613 * ee534;
        ee1341 = ee79 * ((ee64 * ((ee1141 - (ee1152 + ee1215 * ee66/ee63)) *  ee25/ee15 + ee1311/ee520 - ((ee71 * ee1181 - ee1311)/ee520 +  ee607 * ee1186/ee1188)) - ee1322 - (ee1181 * ee63 * ee74 +  ee1322))/ee38 + ee1329/ee286 + ((ee1329 + ee75 * ee1216)/ee286 +  ee615 * ee947/ee799)) - ee618 * ee617;
        ee1343 = ee631 * ee618;
        ee1347 = ee60 * ee1341 - ee1343 + (ee1269 * ee79 - ee1343);
        ee1356 = ee555 * ee518;
        ee1362 = ee701 * ee522;
        ee1368 = ee723 * ee717;
        ee1372 = ee714 * ee496;
        ee1402 = ee682 * ee688;
        ee1413 = ee757 * ee558;
        
        out(j, 0) =  ee129 - (ee124/ee58 - ee136/ee137) - ee173/ee102;
        out(j, 1) =  ee487 - (ee356/ee58 - ee784/ee137) + 1 - ee489/ee102;
        out(j, 2) =  ee631/ee58 + ee1115/ee137 - ee776 - ee533/ee37 - (ee778/ee102 - ee518/ee15);
          out(j, 3) =  ee200 + (ee194/ee58 - ee202/ee137 - ((ee52 * ee209 + ee202)/ee137 - ee136 * ee214/ee216)) - (((4 * (ee107 * (3 * (ee107 * ee221))) * ee96 + ee227 + (ee227 + ee93 * (630 * ee197 + (630 * ee209/ee137 - ee146 * ee214/ee216))) + ((ee151 + ee187 * ee15 * ee29 + ee151)/ee90 - ee244/ee162 - ((ee86 * (ee248 * (ee155 * ee133) * ee157 * ee8 * ee25) + ee244)/ee162 - ee161 * (2 * (ee160 * ee90))/ee260))) * ee99 - ee266 - (ee266 - ee98 * (4 * (ee107 * (3 * (ee107 * ee46))))) + ee200 * ee63)/ee102 - ee173 * ee173/ee279);
          out(j, 4) =  ee365 + (ee351/ee58 - ee124 * ee371/ee137 - ((ee52 * ee382 + ee356 * ee135)/ee137 - ee136 * ee388/ee216)) - (((4 * (ee107 * ee395 - ee304 * ee141) * ee96 + ee143 * ee404 + (ee408 * ee148 + ee93 * (630 * ee362 + (630 * ee382/ee137 - ee146 * ee388/ee216))) + ((ee151 * ee27 + ee334 * ee15 * ee29 + ee424)/ee90 + ee153 * ee431/ee162 + ((ee86 * (ee160 - (ee156 * (ee87 * ee376) + ee438 * ee157) * ee8 * ee25) - ee447 * ee160)/ee162 - ee161 * ee452/ee260))) * ee99 - ee165 * ee338 - (ee468 * ee169 - ee98 * (4 * (ee107 * ee471 + ee304 * ee167))) + ee365 * ee63)/ee102 - ee173 * ee489/ee279);
          out(j, 5) =  ee634 - (ee599/ee58 + ee124 * ee646/ee137 - ((ee631 * ee135 - ee52 * ee661)/ee137 + ee136 * ee666/ee216)) - ((ee165 * ee558 - (ee93 * (630 * ee627 - (630 * ee661/ee137 - ee146 * ee666/ee216)) + ee682 * ee148 + (ee143 * ee688 + 4 * (ee107 * ee691 + ee563 * ee141) * ee96) + ((ee701 + ((ee593 * ee15 - ee122 * ee518) * ee29 + ee151 * ee522))/ee90 - ee153 * ee717/ee162 - ((ee723 * ee160 + ee86 * ((ee158 * ee496 - (ee156 * (ee87 * ee651 + ee653) + ee732 * ee157) * ee8) * ee25))/ee162 - ee161 * ee742/ee260))) * ee99 - (ee757 * ee169 + ee98 * (4 * (ee107 * ee760 - ee563 * ee167))) + ee634 * ee63)/ee102 - ee173 * ee778/ee279);
          out(j, 6) =  ee859 - (ee841/ee58 + ee861/ee137 - ((ee52 * ee873 - ee861)/ee137 + ee784 * ee388/ee216)) - (((4 * (ee802 * ee141 + ee295 * ee395) * ee96 + ee887 + (ee887 + ee93 * (630 * ee856 - (630 * ee873/ee137 + ee402 * ee388/ee216))) - ((ee445 - ee899 + (ee838 * ee15 * ee29 - ee899))/ee90 - ee906/ee162 + ((ee86 * (ee431 - ((ee156 * (ee87 * ee869) - ee438 * ee427) * ee8 * ee25 + ee430)) - ee906)/ee162 - ee465 * ee452/ee260))) * ee99 - ee926 - (ee926 + ee98 * (4 * (ee802 * ee167 - ee295 * ee471))) + ee859 * ee63)/ee102 - ee489 * ee489/ee279);
          out(j, 7) =  ee1024 - (ee1001/ee58 + ee356 * ee646/ee137 - ((ee631 * ee371 - ee52 * ee1037)/ee137 + ee784 * ee666/ee216)) - (((4 * (ee957 * ee141 - ee295 * ee691) * ee96 - ee408 * ee688 + (ee93 * (630 * ee1020 + (630 * ee1037/ee137 - ee402 * ee666/ee216)) - ee682 * ee404) - ((ee701 * ee27 - ee722 + ((ee997 * ee15 - ee348 * ee518) * ee29 + ee424 * ee522))/ee90 - ee447 * ee717/ee162 + ((ee723 * ee431 + ee86 * (ee717 - (ee428 * ee496 - (ee156 * (ee87 * ee1031 + ee1033) + ee732 * ee427) * ee8) * ee25))/ee162 - ee465 * ee742/ee260))) * ee99 + ee468 * ee558 - (ee757 * ee338 + ee98 * (4 * (ee957 * ee167 + ee295 * ee760))) + ee1024 * ee63)/ee102 - ee489 * ee778/ee279);
          out(j, 8) = ee1269/ee58 + ee1271/ee137 + ((ee1271 + ee52 * ee1297)/ee137 + ee1115 * ee666/ee216) - ee1347 - (ee1215/ee37 + ee533 * ee533/R_pow(ee37, 2)) - (((((ee1259 * ee15 - ee1356 - (ee1356 + ee50 * ee1181)) * ee29 + ee1362 + (ee1362 + ee85 * ee1191))/ee90 - ee1368/ee162 - ((ee1368 - ee86 * ((ee88 * ee1131 + ee1372 + ((ee156 * (ee87 * ee1276 - ee1278) - ee732 * ee712 - (ee88 * ee1293 + ee714 * ee639)) * ee8 + ee1372)) * ee25))/ee162 - ee751 * ee742/ee260) - (ee93 * (630 * ee1341 + (630 * ee1297/ee137 + ee686 * ee666/ee216)) - ee1402 + (4 * (ee1223 * ee141 - ee537 * ee691) * ee96 - ee1402))) * ee99 + ee1413 + (ee1413 + ee98 * (4 * (ee1223 * ee167 + ee537 * ee760))) - ee1347 * ee63)/ee102 - ee778 * ee778/ee279 - (ee1181/ee15 + ee518 * ee518/ee520));
        
        }
      }}
  
  return out;
  
}

// //' @rdname bgevd0
// [[Rcpp::export]]
arma::mat bgevd34(Rcpp::List pars, arma::mat X1, arma::mat X2, arma::mat X3, arma::vec yvec, arma::uvec dupid, int dcate, arma::vec other)
{

  arma::vec qavec = X1 * Rcpp::as<arma::vec>(pars[0]);
  arma::vec lsbvec = X2 * Rcpp::as<arma::vec>(pars[1]);
  arma::vec txivec = X3 * Rcpp::as<arma::vec>(pars[2]);
  int nobs = yvec.size();
  arma::mat out = arma::mat(nobs, 25, arma::fill::zeros);
  
  double psuba = other[0];
  double psubb = other[1];
  double alpha = other[2];
  double beta = other[3];
  double hbeta = 0.5 * beta;
  
  if (dcate == 1) {
    qavec = qavec.elem(dupid);
    lsbvec = lsbvec.elem(dupid);
    txivec = txivec.elem(dupid);
  }
  
  double y, qalpha, tqalpha, lsbeta, txi, xi, sbeta, tsbeta, iFa, iFb;
  
  double ee10, ee100, ee1000, ee1001, ee1002, ee1003, ee1004, ee1006, ee1007;
  double ee101, ee1010, ee1011, ee1012, ee1014, ee1015, ee1017, ee102, ee1021, ee1023;
  double ee1024, ee1027, ee103, ee1034, ee1036, ee1038, ee104, ee1040, ee1041, ee1044;
  double ee1047, ee1049, ee105, ee1050, ee1052, ee106, ee1065, ee107, ee1071, ee1074;
  double ee1077, ee1078, ee1079, ee108, ee1080, ee1081, ee1082, ee1083, ee1085, ee1087;
  double ee1088, ee1089, ee109, ee1090, ee1091, ee1092, ee1093, ee1094, ee1095, ee1096;
  double ee1097, ee1098, ee11, ee110, ee1100, ee1102, ee1105, ee1106, ee1108, ee111;
  double ee1110, ee1112, ee1113, ee1114, ee1115, ee1116, ee1117, ee112, ee1120, ee1122;
  double ee1124, ee1125, ee1128, ee113, ee1131, ee1134, ee1135, ee1138, ee114, ee1143;
  double ee1145, ee1146, ee1148, ee115, ee1151, ee1152, ee1159, ee116, ee1162, ee1165;
  double ee1166, ee1167, ee117, ee1171, ee1172, ee1173, ee1174, ee1175, ee1176, ee1178;
  double ee1179, ee118, ee1181, ee1182, ee1184, ee1185, ee1186, ee1187, ee1188, ee1189;
  double ee119, ee1192, ee1197, ee12, ee120, ee1200, ee1205, ee121, ee1210, ee1211;
  double ee1214, ee1217, ee122, ee1221, ee1228, ee1229, ee123, ee1230, ee1239, ee124;
  double ee1241, ee1243, ee1245, ee125, ee1255, ee1258, ee126, ee1264, ee1267, ee127;
  double ee1270, ee1272, ee128, ee1281, ee1283, ee1285, ee1289, ee129, ee1291, ee1292;
  double ee1293, ee1294, ee1296, ee1298, ee1299, ee13, ee130, ee1300, ee1301, ee1302;
  double ee1303, ee1304, ee1306, ee1308, ee131, ee1310, ee1313, ee1315, ee132, ee1321;
  double ee1323, ee1327, ee133, ee1331, ee1337, ee134, ee1341, ee1342, ee1344, ee1346;
  double ee1350, ee1352, ee1353, ee1357, ee136, ee1360, ee1364, ee1367, ee1369, ee137;
  double ee1370, ee1375, ee138, ee1380, ee1382, ee1384, ee1389, ee139, ee1393, ee1394;
  double ee1395, ee14, ee140, ee1402, ee1403, ee1405, ee1408, ee141, ee1412, ee1413;
  double ee1415, ee1417, ee1418, ee1419, ee142, ee1421, ee1422, ee1424, ee1425, ee143;
  double ee1433, ee144, ee1441, ee1446, ee1447, ee1448, ee145, ee1450, ee1453, ee1455;
  double ee1456, ee1458, ee146, ee1461, ee1465, ee1468, ee147, ee1473, ee148, ee1484;
  double ee1487, ee1489, ee149, ee1492, ee1493, ee1498, ee15, ee150, ee1506, ee1507;
  double ee151, ee1510, ee1511, ee1512, ee1514, ee1516, ee1517, ee1518, ee152, ee1520;
  double ee1522, ee1526, ee1529, ee153, ee1530, ee1531, ee1533, ee1537, ee154, ee1540;
  double ee1541, ee1542, ee1543, ee1545, ee155, ee1555, ee156, ee1564, ee1567, ee1568;
  double ee157, ee1570, ee1572, ee1573, ee1577, ee158, ee1581, ee1583, ee1588, ee159;
  double ee1590, ee1592, ee1596, ee16, ee160, ee1600, ee1601, ee1603, ee1605, ee161;
  double ee1610, ee1612, ee1614, ee1616, ee162, ee1622, ee1624, ee1627, ee163, ee1630;
  double ee1637, ee164, ee1645, ee165, ee1650, ee1653, ee1655, ee166, ee1660, ee1662;
  double ee1664, ee1666, ee1668, ee167, ee1670, ee1672, ee1674, ee1676, ee168, ee1681;
  double ee1685, ee1686, ee1688, ee1689, ee169, ee1691, ee1693, ee1694, ee1696, ee1697;
  double ee1698, ee17, ee170, ee1703, ee1709, ee171, ee1714, ee1715, ee172, ee1723;
  double ee1724, ee1728, ee173, ee1730, ee1731, ee174, ee1740, ee1742, ee175, ee1751;
  double ee1754, ee1756, ee176, ee1760, ee1764, ee1765, ee1769, ee1770, ee1771, ee1772;
  double ee1774, ee1776, ee1777, ee178, ee1781, ee1783, ee179, ee1792, ee1795, ee1799;
  double ee18, ee180, ee1801, ee1803, ee1805, ee1808, ee1809, ee181, ee1811, ee1818;
  double ee182, ee1820, ee1821, ee1822, ee1825, ee1828, ee1829, ee183, ee1831, ee184;
  double ee1840, ee1841, ee1843, ee1848, ee1849, ee185, ee1851, ee1852, ee1854, ee1858;
  double ee186, ee1862, ee1863, ee1866, ee1868, ee1869, ee187, ee1873, ee1877, ee1879;
  double ee188, ee1880, ee1883, ee1886, ee189, ee1894, ee1898, ee19, ee190, ee1904;
  double ee191, ee1910, ee1912, ee1914, ee192, ee1922, ee1925, ee1929, ee193, ee1932;
  double ee194, ee1940, ee1942, ee1949, ee195, ee196, ee1961, ee197, ee1972, ee1975;
  double ee198, ee1987, ee199, ee1993, ee2, ee20, ee200, ee2001, ee2003, ee201;
  double ee2011, ee2015, ee202, ee2020, ee2024, ee2026, ee203, ee2039, ee204, ee2043;
  double ee2049, ee205, ee2055, ee206, ee2063, ee207, ee2070, ee209, ee2097, ee21;
  double ee210, ee211, ee2110, ee2114, ee2119, ee212, ee2124, ee2129, ee213, ee2138;
  double ee214, ee215, ee2152, ee216, ee217, ee2170, ee218, ee2199, ee22, ee220;
  double ee2204, ee221, ee2215, ee222, ee2225, ee223, ee224, ee2243, ee2248, ee225;
  double ee2251, ee2253, ee2257, ee2266, ee2268, ee227, ee228, ee2281, ee2294, ee23;
  double ee230, ee2301, ee2302, ee2305, ee231, ee2317, ee2319, ee232, ee233, ee2336;
  double ee2338, ee234, ee2340, ee235, ee236, ee237, ee2372, ee238, ee239, ee24;
  double ee240, ee2407, ee242, ee243, ee244, ee2444, ee245, ee246, ee247, ee2476;
  double ee248, ee249, ee2492, ee25, ee250, ee2502, ee251, ee252, ee253, ee254;
  double ee255, ee256, ee257, ee258, ee2581, ee2585, ee259, ee2591, ee2599, ee26;
  double ee260, ee2603, ee261, ee2614, ee2619, ee262, ee2621, ee263, ee2632, ee264;
  double ee2641, ee2642, ee265, ee2652, ee2656, ee2658, ee266, ee2660, ee2664, ee268;
  double ee2680, ee2682, ee269, ee2690, ee2695, ee2699, ee27, ee270, ee2702, ee2707;
  double ee271, ee272, ee2721, ee2726, ee273, ee2732, ee274, ee2740, ee2744, ee2746;
  double ee2748, ee275, ee2750, ee276, ee2760, ee2767, ee2771, ee2775, ee278, ee279;
  double ee2794, ee2798, ee28, ee280, ee2800, ee2803, ee2807, ee281, ee2811, ee2813;
  double ee2818, ee282, ee2820, ee2822, ee2826, ee283, ee2830, ee2833, ee284, ee2840;
  double ee2842, ee2844, ee2845, ee2851, ee2856, ee286, ee2866, ee287, ee2874, ee2879;
  double ee288, ee2882, ee2884, ee289, ee2893, ee2895, ee2897, ee2899, ee29, ee290;
  double ee2903, ee291, ee2914, ee2917, ee2918, ee292, ee2921, ee293, ee2933, ee2939;
  double ee294, ee2942, ee295, ee2951, ee2957, ee296, ee2969, ee297, ee298, ee2981;
  double ee299, ee2991, ee2997, ee3, ee30, ee300, ee3003, ee301, ee3019, ee302;
  double ee303, ee3031, ee3034, ee304, ee3046, ee305, ee3052, ee306, ee3061, ee3063;
  double ee307, ee3071, ee3077, ee308, ee3087, ee3090, ee3094, ee31, ee310, ee3101;
  double ee3105, ee311, ee3111, ee3117, ee312, ee3123, ee313, ee3135, ee314, ee3142;
  double ee315, ee316, ee317, ee319, ee32, ee320, ee321, ee322, ee323, ee325;
  double ee327, ee329, ee33, ee330, ee331, ee333, ee334, ee335, ee336, ee338;
  double ee339, ee34, ee340, ee341, ee344, ee346, ee348, ee349, ee35, ee350;
  double ee351, ee353, ee357, ee358, ee359, ee36, ee360, ee362, ee363, ee366;
  double ee367, ee368, ee369, ee37, ee371, ee373, ee374, ee377, ee378, ee38;
  double ee380, ee381, ee382, ee383, ee386, ee389, ee39, ee391, ee393, ee395;
  double ee396, ee397, ee398, ee399, ee4, ee40, ee401, ee403, ee404, ee405;
  double ee406, ee407, ee408, ee41, ee410, ee411, ee412, ee415, ee416, ee419;
  double ee42, ee420, ee421, ee423, ee424, ee426, ee427, ee428, ee429, ee43;
  double ee430, ee432, ee433, ee435, ee436, ee439, ee44, ee441, ee442, ee443;
  double ee444, ee446, ee447, ee448, ee45, ee450, ee451, ee452, ee454, ee456;
  double ee459, ee46, ee460, ee461, ee462, ee463, ee465, ee468, ee469, ee47;
  double ee470, ee471, ee472, ee473, ee477, ee479, ee48, ee483, ee484, ee485;
  double ee486, ee49, ee490, ee492, ee494, ee495, ee498, ee5, ee50, ee500;
  double ee501, ee502, ee503, ee506, ee51, ee510, ee512, ee514, ee515, ee517;
  double ee518, ee52, ee520, ee521, ee522, ee523, ee524, ee526, ee527, ee528;
  double ee53, ee530, ee533, ee535, ee536, ee538, ee54, ee541, ee543, ee544;
  double ee546, ee547, ee549, ee55, ee550, ee553, ee556, ee557, ee558, ee56;
  double ee560, ee563, ee569, ee57, ee570, ee572, ee573, ee574, ee577, ee578;
  double ee579, ee58, ee581, ee582, ee584, ee585, ee587, ee588, ee589, ee59;
  double ee590, ee591, ee593, ee594, ee595, ee596, ee598, ee599, ee6, ee60;
  double ee600, ee601, ee603, ee605, ee606, ee608, ee609, ee61, ee610, ee611;
  double ee613, ee614, ee616, ee619, ee62, ee620, ee621, ee622, ee623, ee624;
  double ee626, ee627, ee629, ee63, ee630, ee632, ee634, ee635, ee636, ee64;
  double ee640, ee641, ee647, ee65, ee651, ee652, ee653, ee654, ee655, ee657;
  double ee659, ee66, ee660, ee661, ee662, ee663, ee665, ee666, ee667, ee67;
  double ee671, ee672, ee674, ee68, ee681, ee683, ee688, ee69, ee690, ee695;
  double ee696, ee697, ee698, ee699, ee7, ee70, ee700, ee701, ee703, ee704;
  double ee705, ee706, ee707, ee708, ee709, ee71, ee710, ee711, ee713, ee714;
  double ee715, ee716, ee717, ee718, ee72, ee720, ee723, ee724, ee725, ee726;
  double ee728, ee73, ee730, ee731, ee732, ee733, ee734, ee736, ee737, ee739;
  double ee74, ee744, ee745, ee749, ee75, ee750, ee754, ee755, ee756, ee758;
  double ee759, ee76, ee761, ee762, ee763, ee764, ee765, ee766, ee768, ee769;
  double ee77, ee770, ee772, ee773, ee775, ee777, ee778, ee78, ee782, ee784;
  double ee787, ee788, ee789, ee79, ee791, ee796, ee8, ee80, ee800, ee802;
  double ee804, ee808, ee81, ee811, ee812, ee813, ee814, ee815, ee817, ee819;
  double ee82, ee821, ee822, ee823, ee827, ee829, ee83, ee830, ee831, ee832;
  double ee833, ee836, ee838, ee839, ee84, ee842, ee847, ee848, ee849, ee85;
  double ee851, ee852, ee857, ee858, ee86, ee860, ee861, ee863, ee864, ee866;
  double ee87, ee870, ee873, ee874, ee875, ee876, ee878, ee879, ee88, ee880;
  double ee881, ee883, ee885, ee886, ee888, ee889, ee89, ee890, ee891, ee892;
  double ee893, ee894, ee896, ee9, ee90, ee900, ee901, ee904, ee905, ee907;
  double ee910, ee911, ee912, ee92, ee921, ee924, ee928, ee93, ee931, ee932;
  double ee934, ee935, ee936, ee937, ee939, ee94, ee940, ee941, ee942, ee943;
  double ee944, ee946, ee947, ee948, ee95, ee951, ee954, ee957, ee958, ee959;
  double ee96, ee961, ee964, ee967, ee968, ee969, ee97, ee971, ee974, ee977;
  double ee978, ee979, ee98, ee983, ee985, ee987, ee988, ee989, ee99, ee990;
  double ee992, ee994, ee996, ee997;
  
  for (int j=0; j < nobs; j++) {
    
    y = yvec[j];
    qalpha = qavec[j];
    lsbeta = lsbvec[j];
    txi = txivec[j];
    sbeta = exp(lsbeta);
    xi = 1.5 / (1.0 + exp(-txi)) - 0.5;
    
    iFa = iF(psuba, qalpha, sbeta, xi, alpha, beta);
    iFb = iF(psubb, qalpha, sbeta, xi, alpha, beta);
    
    tqalpha = iFa - (iFb - iFa) * (ell2(alpha) - ell2(psuba)) / (ell2(psuba) - ell2(psubb));
    tsbeta = (iFb - iFa) * (ell2(hbeta) - ell2(1.0 - hbeta)) / (ell2(psuba) - ell2(psubb));
    
    if (y < iFa) { // Gumbel
      
      ee2 = exp(-txi);
      ee3 = 1 + ee2;
      ee5 = 1.5/ee3 - 0.5;
      ee6 = -log(psuba);
      ee7 = -log(alpha);
      ee8 = -log(1 - hbeta);
      ee9 = -log(hbeta);
      ee10 = -log(psubb);
      ee11 = R_pow(ee7, ee5);
      ee12 = 1/ee11;
      ee13 = R_pow(ee6, ee5);
      ee14 = R_pow(ee8, ee5);
      ee15 = R_pow(ee9, ee5);
      ee16 = log(ee6);
      ee17 = R_pow(ee10, ee5);
      ee18 = 1/ee13;
      ee19 = log(ee10);
      ee21 = 1/ee14 - 1/ee15;
      ee22 = 1/ee17;
      ee23 = log(ee8);
      ee24 = log(ee9);
      ee25 = ee12 + ee18;
      ee26 = ee12 + ee22;
      ee27 = ee16 - ee19;
      ee28 = ee26 - ee25;
      ee29 = exp(lsbeta);
      ee30 = log(ee7);
      ee31 = R_pow(ee3, 2);
      ee33 = 1.5 * (ee23/ee14) - 1.5 * (ee24/ee15);
      ee34 = ee24 - ee23;
      ee35 = 1.5 - 3 * (ee2/ee3);
      ee36 = ee16/ee13;
      ee37 = ee28 * ee29;
      ee38 = ee30 - ee16;
      ee39 = 1.5 * ee36;
      ee40 = ee21 * ee27;
      ee41 = ee19/ee17;
      ee42 = ee25 - ee26;
      ee46 = ee18 - (ee28 * ee38/ee27 + ee12);
      ee47 = 1.5 * ee41;
      ee48 = R_pow(ee37 * ee34/ee40, 2);
      ee49 = R_pow(ee34, 2);
      ee53 = ee33 * ee42/ee21 + ee47 - ee39;
      ee54 = ee18 - ee12;
      ee56 = ee46 * ee29/ee21;
      ee58 = y - (ee56 + qalpha);
      ee59 = ee48 * ee21;
      ee60 = ee30/ee11;
      ee61 = ee59 * ee27;
      ee62 = 1.5 * ee60;
      ee63 = ee35/ee14;
      ee64 = ee35/ee15;
      ee67 = ee63 + 2.25 * (ee2 * ee23/(ee14 * ee31));
      ee68 = ee64 + 2.25 * (ee2 * ee24/(ee15 * ee31));
      ee70 = ee67 * ee23 - ee68 * ee24;
      ee71 = ee33 * ee54;
      ee72 = ee71/ee21;
      ee73 = ee35/ee13;
      ee75 = ee53 * ee38/ee27;
      ee77 = ee73 + 2.25 * (ee2 * ee16/(ee13 * ee31));
      ee80 = ee39 - (ee75 + ee72 + ee62);
      ee81 = 3 * ee36;
      ee82 = ee77 * ee16;
      ee83 = ee35/ee17;
      ee84 = R_pow(ee27, 2);
      ee85 = ee37 * ee49;
      ee87 = ee83 + 2.25 * (ee2 * ee19/(ee17 * ee31));
      ee88 = ee22 - ee12;
      ee89 = R_pow(ee29, 2);
      ee90 = ee87 * ee19;
      ee92 = ee80 * ee27/ee28;
      ee93 = ee70 * ee42;
      ee94 = (ee33 * (2 * ee88 - 2 * ee54)/ee21 + ee81 - 3 * ee41) *  ee33;
      ee95 = R_pow(ee21, 2);
      ee97 = ee53 * ee29 * ee49;
      ee99 = ee85 * ee58/ee61;
      ee102 = (ee93 + ee94 * ee2/ee31)/ee21 + ee90 - ee82;
      ee103 = ee48 * ee95;
      ee104 = R_pow(ee28, 2);
      ee106 = exp(-(ee40 * ee58/ee37 - ee30));
      ee108 = ee46 * ee27/ee28;
      ee111 = ee97 * ee58/ee61 + ee92;
      ee112 = ee99 + ee108;
      ee115 = 1.5 - ((3 * (1 + 2 * ee2) + 3 * ee3 - 12 * ee2)/ee3 +  3) * ee2/ee3;
      ee116 = ee39 - ee62;
      ee117 = ee70 * ee54;
      ee118 = ee104 * ee89;
      ee119 = ee103 * ee84;
      ee121 = ee118 * ee49/ee119;
      ee122 = 3 * ee60;
      ee123 = ee47 - ee62;
      ee125 = 2 * ee72 + ee122;
      ee126 = ee35/ee11;
      ee127 = R_pow(ee53, 2);
      ee128 = ee59 * ee84;
      ee130 = ee126 + 2.25 * (ee2 * ee30/(ee11 * ee31));
      ee133 = ee33 * ee116 * ee2/ee31;
      ee134 = ee130 * ee30;
      ee137 = ee80 * ee28 * ee29/ee21;
      ee138 = 2 * ee121;
      ee139 = ee117 + ee133;
      ee142 = R_pow(ee33, 2) * ee2/ee31;
      ee144 = ee102 * ee38/ee27;
      ee145 = ee33 * (ee81 - ee125);
      ee146 = 2 * ee139;
      ee147 = ee102 * ee58;
      ee148 = ee103 * ee27;
      ee150 = ee145 * ee2/ee31;
      ee151 = 2 * ee30;
      ee152 = 2 * ee16;
      ee155 = ee31 * ee21;
      ee160 = (ee82 - (ee144 + (ee117 + ee150)/ee21 + ee134)) *  ee27/ee28;
      ee163 = (ee53 * (ee56 + (ee138 - 1) * ee58) + ee137) * ee29 *  ee49/ee61;
      ee164 = ee127 * ee28;
      ee165 = ee27/ee28;
      ee166 = 3 * ee82;
      ee167 = ee111 * ee21;
      ee168 = R_pow(ee40/ee37, 2);
      ee169 = ee48 * ee84;
      ee170 = ee85/ee61;
      ee173 = ((ee53 * (2 * ee99 + ee152 - ee151)/ee27 + ee81 -  ee125) * ee53 * ee2 * ee29/ee155 - ee147) * ee29 * ee49/ee61;
      ee176 = R_pow(ee111, 2);
      ee178 = (ee115/ee14 + (1.5 * ee67 + 3 * ee63) * ee2 * ee23/ee31) *  ee23 - (ee115/ee15 + (1.5 * ee68 + 3 * ee64) * ee2 *  ee24/ee31) * ee24;
      ee179 = (ee115/ee13 + (1.5 * ee77 + 3 * ee73) * ee2 * ee16/ee31) *  ee16;
      ee180 = ee97/ee61;
      ee181 = ee127 * ee2;
      ee182 = R_pow(ee112, 2);
      ee183 = (ee165 - ee28 * ee89 * ee49/ee148) * ee46;
      ee185 = ee85 * (y - ((2 * (ee104 * ee29 * ee49 * ee58/ee128) +  2 * ee46) * ee29/ee21 + qalpha))/ee61;
      ee187 = 2 * (ee70 * ee21 + ee142) - 8 * ee142;
      ee188 = 8 * ee121;
      ee190 = ee102 * ee28 + ee181/ee31;
      ee191 = ee93 + (ee94 - 2 * (ee164 * ee89 * ee49/ee128)) *  ee2/ee31;
      ee192 = ee111 * ee112;
      ee194 = ee112 * ee21;
      ee195 = ee137 - ee53 * ee58;
      ee197 = ee191/ee21 + ee90;
      ee199 = (((ee21 * (2 * (ee70 * ee88 + ee33 * ee123 * ee2/ee31) -  ee146) + ee28 * ee187)/ee21 + ee33 * (2 * ee123 - 2 *  ee116) * ee2/ee31 + 2 * (ee70 * ee88 + ee33 * ee123 *  ee2/ee31) - ee146)/ee21 + ee166 - 3 * ee90) * ee33 +  ee70 * (3 * ee116 - 3 * ee123);
      ee200 = ee178 * ee42;
      ee202 = ee167 * ee27/ee37;
      ee204 = ee176 * ee2/ee31;
      ee205 = (ee115/ee17 + (1.5 * ee87 + 3 * ee83) * ee2 * ee19/ee31) *  ee19;
      ee207 = ee194 * ee27/ee37;
      ee211 = (ee199 * ee2/ee31 + ee200)/ee21 + ee205 - ee179;
      ee213 = 1 - ee106;
      ee214 = 2 * (ee190 * ee58);
      ee215 = 2 * ee163;
      ee216 = 2 * ee195;
      ee217 = 2 * ee92;
      ee220 = 2 * ee80;
      ee221 = 3 - ee138;
      ee222 = 3 * ee134;
      ee223 = y - (2 * ee56 + qalpha);
      ee224 = ee173 - (ee160 + ee204);
      ee227 = ee164 * ee2 * ee89 * ee49;
      ee228 = ee168 * ee48;
      ee231 = ee48 * ee31 * ee95 * ee84;
      ee232 = ee207 - ee170;
      ee234 = ee118/(ee95 * ee84);
      ee235 = 2 * ee106;
      ee236 = 2/ee168;
      ee238 = 3 * ee106 - 2;
      ee239 = ee106 - 1;
      ee244 = ee211 * ee58;
      ee245 = (ee144 + (ee117 + (ee53 * (ee216 + 8 * (ee53 * ee104 *  ee89 * ee49 * ee58/ee119)) * ee2/ee31 - ee214) * ee29 *  ee49/ee169 + ee150)/ee21 + ee134 - ee82) * ee28;
      ee247 = ee102 * ee46;
      ee248 = ((((ee28 * (4 * (ee53 * ee80 * ee2 * ee29/ee155 -  ee147) + 8 * (ee227 * ee58/ee231)) - ee214) * ee29 *  ee49/ee169 + (2 * (ee53 * ee80 * ee28 * ee89 * ee49/ee128) +  3 * ee145) * ee2/ee31 + 3 * ee117)/ee21 + 3 * ee144 +  ee222 - ee166) * ee53 - 3 * (ee102 * ee80)) * ee2;
      ee249 = ee197 - ee82;
      ee255 = (ee111 * ee28 + ee53 * (ee85 * (y - ((ee28 * (2 +  ee16 - ee30)/ee27 + ee18 - ee12) * ee29/ee21 + qalpha))/ee61 +  ee108 + 1)) * ee29 * ee49/ee61;
      ee257 = (ee179 - (ee211 * ee38/ee27 + (ee178 * ee54 + (ee33 *  (ee166 - (((ee54 * ee187 + 2 * (ee139 * ee21))/ee21 +  ee146 + 2 * ee133)/ee21 + ee222)) + 3 * (ee70 * ee116)) *  ee2/ee31)/ee21 + (ee115/ee11 + (1.5 * ee130 + 3 * ee126) *  ee2 * ee30/ee31) * ee30)) * ee27/ee28;
      ee258 = ee163 - (ee192 + ee92);
      ee259 = ee163 - ee92;
      ee261 = (ee53 * ee46 * ee221 + (ee39 + ee220 - (((ee53 *  (ee188 - 4) * ee58 + ee216) * ee28 * ee29 * ee49/ee169 +  ee71)/ee21 + ee75 + ee62)) * ee28) * ee29/ee21;
      ee268 = ee53 * (ee220 + ee81 - (ee53 * (2 * (ee28 * ee46 *  ee89 * ee49/ee148) + ee151 - ee152)/ee27 + (2 * (ee195 *  ee28 * ee29 * ee49/ee169) + 2 * ee71)/ee21 + ee122)) *  ee2/ee31;
      ee270 = ee53 * ee49/ee48;
      ee272 = ee228 * ee31;
      ee276 = ee112 * ee213;
      ee278 = ee182 + ee183 + ee185;
      ee279 = ((ee236 - 2 * ee234) * ee49/ee48 + 1 + 2 * (1 -  ee138))/ee168;
      ee280 = (ee165 - ee28 * ee221 * ee89 * ee49/ee148) * ee46;
      ee286 = ee28 * (2 * y - ((ee28 * (2 + ee152 - ee151)/ee27 +  ee18 + ee18 - 2/ee11) * ee29/ee21 + 2 * qalpha)) * ee29 *  ee49/ee61 + 1 + 2 * ee108;
      ee289 = ee28 * (8 - ee188) * ee29/ee40;
      ee291 = (4 - ee188) * ee58 + 4 * ee223;
      ee292 = ee49/ee48;
      ee293 = 2 * ee190;
      ee294 = 2 * ee173;
      ee295 = 2 * ee160;
      ee297 = 2 * ee202 + 2 * (ee180 - ee202);
      ee298 = 2 * ee180;
      ee299 = 2 * ee232;
      ee300 = 2 * ee207;
      ee301 = 2 * ee183;
      ee302 = 2 * ee185;
      ee303 = 2 * ee170;
      ee304 = ee235 - 2;
      ee305 = 4 * ee46;
      ee306 = R_pow(ee106, 2);
      
      out(j, 0) = ((2 * ee213 + ee235 - 2) * ee106 + 1) * R_pow(ee21, 3) * ee106 * R_pow(ee27, 3)/(R_pow(ee28, 3) * R_pow(ee29, 3));
      out(j, 1) = (((ee194 * ee238 * ee27/ee37 + 2 * (ee276 * ee21 * ee27/ee37 -
        ee170) + ee303) * ee106 + ee276 * (1 + ee106) * ee21 * ee27/ee37 -
        ee170) * ee21 * ee27/ee37 - ee292) * ee106 +
        (ee21 * (ee299 - ee300) * ee27/ee37 + (ee170 * ee21 * ee27/ee37 +
        ee292));
      out(j, 2) = (((ee167 * (ee306 - 1) * ee27/ee37 + ee180 + (2 * (ee167 * ee239 * ee27/ee37 +
        ee180) - (ee167 * ee238 * ee27/ee37 +
        ee298)) * ee106) * ee21 * ee27/ee29 + ee270) * ee106 +
        ee21 * ee297 * ee27/ee29 - 2 * ee270) * ee2/(ee31 * ee28);
      out(j, 3) = (ee278 * ee21 * ee27/ee37 + ee112 * (ee300 - (ee299 +
        ee303)) * ee106 - ee286 * ee28 * ee29 * ee49/ee61) * ee106 -
        ee28 * (2 * ee112 - ee286) * ee29 * ee49/ee61;
      out(j, 4) = ((ee258 * ee21 * ee27/ee37 + ee255) * ee106 + (ee202 +
        ee180) * ee112 - (ee255 + ee111 * ee232)) * ee2/ee31;
      out(j, 5) = ((ee111 * (ee297 - ee298) * ee106 * ee2/ee31 - ((ee197 +
        2 * (ee111 * ee53 * ee2/ee31) - ee82) * ee29 * ee49/ee61 +
        ee224 * ee21 * ee27/ee37)) * ee106 - (ee82 - ee197) * ee29 * ee49/ee61) * ee2/ee31;
      out(j, 6) = ((((ee279 + (ee291 - ee289) * ee28 * ee29/ee40) * ee49/ee48 +
        1) * ee28/ee27 + ee305) * ee29/ee21 + qalpha -
        y) * ee28 * ee29 * ee49/ee61 + (((ee182 * ee304 + 2 * (ee182 * ee213 +
        ee183 + ee185) - (ee301 + ee302)) * ee106 + ee182 +
        3 * ee183 + 3 * ee185) * ee112 + ee280 + ee85 * (y - ((ee291 * ee104 * ee29 * ee49/ee128 +
        ee305) * ee29/ee21 + qalpha))/ee61) * ee106 -
        (ee112 * (ee301 + 2 * ee182 + ee302 - 2 * ee278) +
        ee280);
      out(j, 7) = (((ee261 + ((ee28 * (2 * (ee85 * ee223/ee61) + ee16 -
        ee30)/ee27 + ee18 - ee12) * ee29/ee21 + qalpha - y) * ee53) * ee29 * ee49/ee61 +
        ee111 * (ee182 * ee306 - (ee183 +
        ee185)) + ee112 * ((2 * (ee163 + ee192 * ee239 - ee92) + ee217 -
        (ee192 * ee238 + ee215)) * ee106 + ee215 - (ee192 + ee217)) -
        ee92) * ee106 + ee92 - ((((((ee279 + ee28 * (2 * ee223 -
        ee289) * ee29/ee40) * ee49/ee48 + 1 + ee16 - ee30) * ee28/ee27 +
        ee18 - ee12) * ee29/ee21 + qalpha - y) * ee53 + ee261) * ee29 * ee49/ee61 +
        ee112 * (ee215 - (2 * ee258 + 2 * ee192 +
        ee217)))) * ee2/ee31;
      out(j, 8) = (((((ee249 * ee28 + ee127 * ((ee236 - 4 * ee234) * ee49/ee48 +
        2) * ee2/ee31) * ee49/ee228 + ((ee191 - (ee127 * (2 -
        ee188) * ee2/ee31 + ee293) * ee28 * ee89 * ee49/ee128)/ee21 +
        ee90 - ee82) * ee28)/ee27 + ee247 + ee268 - ee245) * ee29/ee21 -
        ee147) * ee29 * ee49/ee61 + (((ee245 - (ee247 +
        ee268)) * ee29/ee21 + ee147) * ee29 * ee49/ee61 + ee160 -
        (ee224 * ee112 + 2 * (ee259 * ee111 * ee2/ee31))) * ee106 -
        (ee160 + ee111 * (ee215 - (2 * ee259 + ee217)) * ee2/ee31)) * ee2/ee31;
      out(j, 9) = -((((((((ee199 - ee53 * (ee28 * (4 * ee102 - 8 * (ee227/ee231)) +
        ee293) * ee89 * ee49/ee128) * ee2/ee31 + ee200)/ee21 +
        ee205 - ee179) * ee28 + (ee197 + 2 * ee249 + 2 * (ee181 * ee49/(ee272 * ee28)) -
        ee82) * ee53 * ee2 * ee49/ee272)/ee27 -
        ee248/ee31) * ee29/ee21 - ee244) * ee29 * ee49/ee61 +
        ((ee244 + ee248 * ee29/ee155) * ee29 * ee49/ee61 + ((ee176 * ee304 * ee2/ee31 +
        ee294 - (2 * (ee173 + ee176 * ee239 * ee2/ee31 -
        ee160) + ee295)) * ee106 + ee204 + 3 * ee160 -
        3 * ee173) * ee111 * ee2/ee31 + ee257) * ee106 - (ee257 +
        ee111 * (2 * ee224 + ee295 + 2 * ee204 - ee294) * ee2/ee31)) * ee2/ee31);
      
    } else {
      
      if (y > iFb) { // GEV

        ee2 = exp(-txi);
        ee3 = 1 + ee2;
        ee5 = 1.5/ee3 - 0.5;
        ee6 = -log(1 - hbeta);
        ee7 = -log(hbeta);
        ee8 = R_pow(ee6, ee5);
        ee9 = R_pow(ee7, ee5);
        ee10 = -log(alpha);
        ee11 = exp(lsbeta);
        ee12 = y - qalpha;
        ee13 = R_pow(ee10, ee5);
        ee15 = 1/ee8 - 1/ee9;
        ee16 = ee15 * ee12;
        ee18 = ee16/ee11 + 1/ee13;
        ee19 = 1/ee5;
        ee20 = log(ee6);
        ee21 = log(ee7);
        ee22 = 1 + ee19;
        ee23 = R_pow(ee3, 2);
        ee24 = log(ee10);
        ee25 = ee20/ee8;
        ee26 = ee21/ee9;
        ee27 = 1.5 * ee25;
        ee28 = 1.5 * ee26;
        ee29 = ee27 - ee28;
        ee30 = R_pow(ee18, ee22);
        ee31 = 1.5 - 3 * (ee2/ee3);
        ee32 = 1.5 * (ee24/ee13);
        ee35 = ee29 * ee12/ee11 + ee32;
        ee36 = ee19 + 2;
        ee37 = log(ee18);
        ee38 = R_pow(ee18, ee36);
        ee39 = R_pow(ee5, 2);
        ee40 = ee31/ee8;
        ee41 = ee31/ee9;
        ee44 = ee40 + 2.25 * (ee2 * ee20/(ee8 * ee23));
        ee45 = ee41 + 2.25 * (ee2 * ee21/(ee9 * ee23));
        ee46 = ee44 * ee20;
        ee47 = ee45 * ee21;
        ee48 = R_pow(ee18, ee19);
        ee49 = ee31/ee13;
        ee50 = ee46 - ee47;
        ee51 = 1/ee30;
        ee54 = ee35 * ee22/ee38;
        ee55 = ee49 + 2.25 * (ee2 * ee24/(ee13 * ee23));
        ee57 = ee35/ee30 + 1.5 * (ee37/(ee48 * ee5));
        ee60 = ee50 * ee12/ee11 + ee55 * ee24;
        ee64 = (4.5/(ee3 * ee5) - 3) * ee2/ee3 + 1.5;
        ee65 = 2 * ee22;
        ee66 = ee35 * ee2;
        ee67 = ee18 * ee23;
        ee71 = 1.5 * (ee15/ee5) + ee27 - ee28;
        ee72 = ee54 + 1.5 * (ee37/(ee30 * ee39));
        ee73 = R_pow(ee18, ee65);
        ee75 = 12 * ee2;
        ee76 = 3 * (1 + 2 * ee2);
        ee77 = 3 * ee3;
        ee78 = R_pow(ee18, ee19 + 3);
        ee79 = ee73 * ee5;
        ee80 = ee23 * ee5;
        ee81 = R_pow(ee15/ee11, 2);
        ee84 = ee71/ee30;
        ee85 = 1.5 - ((ee76 + ee77 - ee75)/ee3 + 3) * ee2/ee3;
        ee86 = R_pow(ee11, 2);
        ee88 = ee35 * ee15/ee18;
        ee90 = ee18 * ee11;
        ee92 = ee64 * ee37 + 1.5 * (ee66/ee67);
        ee94 = ee35 * ee36/ee78;
        ee95 = ee22/ee38;
        ee96 = ee94 + 1.5 * (ee37/(ee38 * ee39));
        ee98 = ee16/ee90;
        ee99 = 1.5 * (ee37/(ee30 * ee5));
        ee101 = 1/ee79 - ee95;
        ee102 = 2 * ee88;
        ee103 = 2/ee30;
        ee104 = ee60 + 1.5 * (ee66/ee80);
        ee105 = ee35 * ee64;
        ee108 = R_pow(ee35, 2) * ee2/ee67;
        ee110 = R_pow(ee18, ee22 + ee65) * ee5;
        ee111 = ee81 * ee86;
        ee112 = ee64 * ee15;
        ee113 = ee23 * ee39;
        ee116 = R_pow(ee15, 2);
        ee117 = 1.5 * (ee15/ee39);
        ee121 = ((ee57/ee30 - ee99)/ee5 - ee54) * ee15 + ee84;
        ee122 = ee60 - ee108;
        ee123 = ee72 * ee35;
        ee133 = ee96 * ee35 * ee22 * ee2/ee23 + (1.5 * (ee72 * ee2 *  ee37/ee23) - ee92/ee30)/ee39 - (ee60 * ee22 + 1.5 * (ee66/ee113))/ee38;
        ee138 = (ee22 * ee29 + ee117)/ee38;
        ee140 = (ee85/ee8 + (1.5 * ee44 + 3 * ee40) * ee2 * ee20/ee23) *  ee20 - (ee85/ee9 + (1.5 * ee45 + 3 * ee41) * ee2 * ee21/ee23) *  ee21;
        ee141 = (ee84 - ee72 * ee15)/ee30;
        ee145 = ee22 * ee15 * ee12/(ee38 * ee11);
        ee146 = (1.5 * (ee57 * ee2 * ee37/ee80) - ee92/ee48)/ee5;
        ee147 = ee71/ee38;
        ee151 = ee15 * ee101 * ee12/ee11 + ee51;
        ee153 = 2 * ee98;
        ee154 = 2 * (ee116/ee111);
        ee155 = 2 * ee29;
        ee156 = 2/ee38;
        ee157 = (ee121 * (ee51 - ee103) + ee141 + 2 * (ee57 * ee15/ee79))/ee5;
        ee160 = (((((3 * (2.25 * (ee2/ee23) - ee31 * ee5) - (27 *  (ee2/ee23) + 3 * (ee31 * ee5)))/ee5 - 3 * ee31)/ee5 +  ee76 + ee77 - ee75)/ee3 + 3) * ee2/ee3 - 1.5) * ee37;
        ee162 = ((ee46 - ((ee60 * ee15 + ee35 * (3 * ee25 - (ee102 +  3 * ee26)) * ee2/ee23)/ee18 + ee47)) * ee22 + (ee112 +  3 * ((ee27 - (ee88 + ee28)) * ee2/ee23))/ee39)/ee18 +  (2 * (ee72 * ee71 * ee2/ee23) - (ee133 * ee15 + ((ee112 +  3 * (ee29 * ee2/ee23))/ee5 + ee46 - ee47)/ee30))/ee5;
        ee163 = ee104/ee30;
        ee164 = (ee57 * (ee16/(ee79 * ee11) + ee51 - ee51) + ee99)/ee5;
        ee166 = ee140 * ee12/ee11;
        ee169 = R_pow(ee57, 2)/ee5;
        ee170 = ee138 + (ee147 - ee96 * ee15) * ee22;
        ee171 = (ee145 - ee51)/ee30;
        ee172 = (ee85/ee13 + (1.5 * ee55 + 3 * ee49) * ee2 * ee24/ee23) *  ee24;
        ee175 = ((ee27 + ee155 - (ee28 + ee102)) * ee12/ee11 + ee32) *  ee15/ee18 + ee28;
        ee180 = R_pow(ee18, 2/ee5 + 3) * ee5;
        ee181 = R_pow(ee18, 2);
        ee183 = ee81 * ee23 * ee86;
        ee184 = ee81 * ee11;
        ee185 = ee22 * (1/ee38 + ee156 - ee15 * ee36 * ee12/(ee78 *  ee11));
        ee186 = ee23 * ee11;
        ee187 = R_pow(ee29, 2);
        ee189 = ee101/ee30;
        ee190 = ee36/ee78;
        ee192 = 1.5 * ee122 + 2 * ee105;
        ee193 = 2 * (ee151/ee30);
        ee194 = 2 * (ee16/(ee110 * ee11));
        ee195 = 2/ee110;
        ee196 = ee103 - ee51;
        
        out(j, 0) =  - (((ee22 * (2/ee180 - ee190) - (ee195 - 2 * ee189)/ee5)/ee5 +
          2 * (ee22/R_pow(ee18, 3))) * R_pow(ee15, 3)/R_pow(ee11, 3));
        out(j, 1) = -((((ee171 + ee151 * ee196 + (ee189 - ee195) * ee15 * ee12/ee11)/ee5 +
          (ee15 * (1/ee180 - ee190) * ee12/ee11 +
          ee156) * ee22)/ee5 - ee22 * (2 - ee153)/ee181) * ee116/ee86);
        out(j, 2) = -((((ee121 * ee196 + ee57 * ((1/ee73 - 2/ee73)/ee5 -
          ee95) * ee15 - ee141)/ee5 + (((ee57/ee38 - 1.5 * (ee37/(ee38 * ee5)))/ee5 -
          ee94) * ee15 + ee147) * ee22 + ee138)/ee5 -
          (ee22 * (ee155 - ee102) + ee117)/ee181) * ee15 * ee2/(ee23 * ee86));
        out(j, 3) = -((((ee185 + (2 * ee171 + ee193 - ee194)/ee5) * ee15 * ee12/ee11 -
          ee51)/ee5 - (ee15 * (3 - ee153) * ee12/ee90 -
          1) * ee22/ee18) * ee15/ee11);
        out(j, 4) = -((((ee164 + ee54 - (ee157 - ee170) * ee12/ee11) * ee15 -
          ee84)/ee5 + (ee22 * (ee27 - ee175) + 1.5 * ((1 - ee98) * ee15/ee39))/ee18) * ee2/ee186);
        out(j, 5) = -(ee162 * ee2/ee186);
        out(j, 6) = -(((((ee185 - (ee194 + 2 * ((ee51 - ee145)/ee30) -
          ee193)/ee5) * ee15 * ee12/ee11 - ee51)/ee5 + (ee15 * (ee153 -
          3) * ee12/ee90 + 1) * ee22/ee18) * ee12 - ee15 * (ee154 -
          3)/ee184) * ee15/ee11 - 1);
        out(j, 7) = -(((((ee164 + (ee170 - ee157) * ee12/ee11 + ee54) * ee15 -
          ee84)/ee5 - ((ee175 - ee27) * ee22 + 1.5 * ((ee98 -
          1) * ee15/ee39))/ee18) * ee12/ee11 + (ee15 * (3 - ee154)/ee111 -
          1/ee15) * ee29) * ee2/ee23);
        out(j, 8) = -(((ee162 * ee12 + (ee50 * ee15 + ee187 * (2 - ee154) * ee2/ee23)/ee184)/ee11 -
          ee50/ee15) * ee2/ee23);
        out(j, 9) = ((((((ee104 * (ee51 + ee103 - ee51) + (2 * ee169 -
          2 * ee123) * ee2/ee23 + 2 * ((ee123 - ee169) * ee2/ee23 +
          ee146 - ee163) - 2 * ee146) * ee57 + 4.5 * ee31 + 6.75 * (ee2/ee80))/ee5 +
          ee133 * ee35 + ee75 - (2 * (ee104 * ee72) +
          ee76 + ee77))/ee3 - 3) * ee2/ee3 + (ee160 + ((1.5 * ((ee123 * ee2/ee23 +
          ee146 - ee163) * ee37) - 2 * (ee57 * ee92))/ee5 -
          ee192/ee18) * ee2/ee23 - (ee160 - ee192 * ee2/ee67)/ee48)/ee5 +
          (ee166 + (ee105 + 3 * ee60) * ee2/ee80 + ee172)/ee30 +
          1.5)/ee5 + ee140/ee15 - (((ee166 + ee172 - (ee60 + 2 * ee60 -
          2 * ee108) * ee35 * ee2/ee67) * ee22 + (ee105 + 3 * ee122) * ee2/ee113)/ee18 +
          (ee46 + 2 * ee50 - (ee47 + 2 * (ee187 * ee15 * ee2/ee183))) * ee29 * ee2/ee183)) * ee2/ee23;
        
      } else {
        
        ee2 = -log(psuba);
        ee4 = exp(-txi);
        ee5 = 1 + ee4;
        ee7 = 1.5/ee5 - 0.5;
        ee8 = -ee7;
        ee9 = R_pow(ee2, ee8);
        ee11 = -log(alpha);
        ee12 = R_pow(ee11, ee8);
        ee14 = exp(lsbeta);
        ee15 = (ee9 - ee12) * ee14;
        ee18 = -log(1 - hbeta);
        ee19 = R_pow(ee18, ee8);
        ee21 = -log(hbeta);
        ee22 = R_pow(ee21, ee8);
        ee23 = ee19 - ee22;
        ee24 = ee15/ee23;
        ee25 = qalpha + ee24;
        ee26 = y - ee25;
        ee28 = -log(psubb);
        ee29 = R_pow(ee28, ee8);
        ee31 = (ee29 - ee12) * ee14;
        ee32 = ee31/ee23;
        ee34 = qalpha + ee32 - ee25;
        ee35 = ee26/ee34;
        ee36 = R_pow(ee35, 5);
        ee38 = 1 - 1;
        ee39 = ee26 * ee38;
        ee40 = R_pow(ee34, 2);
        ee42 = 1/ee34 + ee39/ee40;
        ee44 = 2 * (ee42 * ee35);
        ee46 = R_pow(ee35, 3);
        ee48 = 4 * (ee42 * ee46);
        ee50 = R_pow(ee35, 2);
        ee52 = 3 * (ee42 * ee50);
        ee57 = 540 * ee44 + (70 * ee48 - 315 * ee52) - 420 * ee42;
        ee59 = R_pow(ee35, 4);
        ee61 = 5 * (ee42 * ee59);
        ee69 = 70 * ee59 - 315 * ee46 + 540 * ee50 - 420 * ee35 +  126;
        ee71 = ee36 * ee57 + ee61 * ee69;
        ee72 = log(ee11);
        ee73 = log(ee2);
        ee74 = ee72 - ee73;
        ee76 = log(ee28);
        ee77 = ee73 - ee76;
        ee81 = log(ee21);
        ee82 = log(ee18);
        ee83 = ee81 - ee82;
        ee84 = (y - (ee25 - ee34 * ee74/ee77)) * ee83;
        ee86 = ee34 * ee83/ee77;
        ee90 = exp(-(ee84/ee86 - ee72));
        ee92 = ee36 * ee69;
        ee93 = 1 - ee92;
        ee97 = (1 - ee38 * ee74/ee77) * ee83;
        ee100 = ee38 * ee83/ee77;
        ee101 = ee84 * ee100;
        ee102 = R_pow(ee86, 2);
        ee104 = ee97/ee86 + ee101/ee102;
        ee105 = ee90 * ee104;
        ee107 = ee71 * ee90 + ee93 * ee105;
        ee108 = y - qalpha;
        ee109 = ee108 * ee23;
        ee111 = ee109/ee14 + ee12;
        ee113 = -1/ee7;
        ee114 = ee113 - 1;
        ee115 = R_pow(ee111, ee114);
        ee116 = ee23/ee14;
        ee117 = ee113 * ee116;
        ee118 = ee115 * ee117;
        ee120 = R_pow(ee111, ee113);
        ee124 = 630 * ee48;
        ee125 = 1 - ee35;
        ee126 = R_pow(ee125, 4);
        ee127 = ee124 * ee126;
        ee129 = -630 * ee59;
        ee130 = R_pow(ee125, 3);
        ee132 = 4 * (ee42 * ee130);
        ee134 = ee127 + ee129 * ee132;
        ee136 = ee129 * ee126;
        ee137 = ee136 * ee38;
        ee139 = ee134/ee34 - ee137/ee40;
        ee141 = ee136/ee34;
        ee145 = 1 + 1/ee7;
        ee146 = ee145 + 1;
        ee148 = R_pow(ee111,  - ee146);
        ee149 = ee145 * ee116;
        ee150 = ee148 * ee149;
        ee153 = R_pow(ee111,  - ee145);
        ee155 = ee92 * ee150 - ee71 * ee153;
        ee156 = ee155 * ee116;
        ee159 = 630 * ee59;
        ee161 = ee159 * ee132 - ee127;
        ee163 = ee159 * ee126;
        ee164 = ee163 * ee38;
        ee166 = ee161/ee34 - ee164/ee40;
        ee168 = ee163/ee34;
        ee172 = ee83/ee86;
        ee174 = ee93 * ee90;
        ee175 = ee83 * ee100;
        ee176 = ee175/ee102;
        ee179 = ee139 * ee120 - ee141 * ee118 + ee156/ee7 + (ee166 *  ee90 + ee168 * ee105) + (ee107 * ee172 - ee174 * ee176);
        ee181 = ee92 * ee153;
        ee182 = ee181 * ee116;
        ee188 = ee141 * ee120 + ee182/ee7 + ee168 * ee90 + ee174 *  ee172;
        ee191 = ee71 * ee105;
        ee192 = ee61 * ee57;
        ee194 = ee38/ee40;
        ee196 = 2 * (ee38 * ee34);
        ee197 = ee39 * ee196;
        ee198 = R_pow(ee40, 2);
        ee201 = ee194 + ee197/ee198 + ee194;
        ee204 = 5 * (ee42 * ee48 + ee201 * ee59);
        ee210 = 4 * (ee42 * ee52 + ee201 * ee46);
        ee215 = 3 * (ee42 * ee44 + ee201 * ee50);
        ee221 = 2 * (ee42 * ee42 + ee201 * ee35);
        ee225 = 70 * ee210 - 315 * ee215 + 540 * ee221 - 420 * ee201;
        ee228 = ee192 + ee204 * ee69 + (ee36 * ee225 + ee192);
        ee232 = ee97 * ee100;
        ee233 = ee232/ee102;
        ee235 = 2 * (ee100 * ee86);
        ee236 = ee101 * ee235;
        ee237 = R_pow(ee102, 2);
        ee240 = ee233 + ee236/ee237 + ee233;
        ee242 = ee105 * ee104 - ee90 * ee240;
        ee245 = ee191 - ee228 * ee90 + (ee191 + ee93 * ee242);
        ee246 = ee71 * ee118;
        ee249 = ee114 - 1;
        ee250 = R_pow(ee111, ee249);
        ee251 = ee114 * ee116;
        ee252 = ee250 * ee251;
        ee253 = ee252 * ee117;
        ee258 = ee124 * ee132;
        ee259 = 630 * ee210;
        ee261 = ee258 - ee259 * ee126;
        ee262 = R_pow(ee125, 2);
        ee264 = 3 * (ee42 * ee262);
        ee268 = 4 * (ee42 * ee264 - ee201 * ee130);
        ee271 = ee261 + (ee258 + ee129 * ee268);
        ee273 = ee134 * ee38;
        ee274 = ee273/ee40;
        ee276 = ee137 * ee196;
        ee279 = ee271/ee34 - ee274 - (ee274 - ee276/ee198);
        ee281 = ee139 * ee118;
        ee286 = ee146 + 1;
        ee288 = R_pow(ee111,  - ee286);
        ee289 = ee146 * ee116;
        ee290 = ee288 * ee289;
        ee291 = ee290 * ee149;
        ee293 = ee71 * ee150;
        ee297 = ee92 * ee291 - ee293 - (ee293 - ee228 * ee153);
        ee298 = ee297 * ee116;
        ee303 = ee159 * ee268 - ee258 - ee261;
        ee305 = ee161 * ee38;
        ee306 = ee305/ee40;
        ee308 = ee164 * ee196;
        ee311 = ee303/ee34 - ee306 - (ee306 - ee308/ee198);
        ee313 = ee166 * ee105;
        ee320 = ee107 * ee176;
        ee322 = ee175 * ee235;
        ee323 = ee322/ee237;
        ee327 = ee279 * ee120 - ee281 - (ee281 - ee141 * ee253) +  ee298/ee7 + (ee311 * ee90 + ee313 + (ee313 + ee168 *  ee242)) + (ee245 * ee172 - ee320 - (ee320 - ee174 * ee323));
        ee329 = ee179 * ee179;
        ee330 = R_pow(ee188, 2);
        ee334 = ee71 * ee242;
        ee335 = ee228 * ee105;
        ee336 = ee334 - ee335;
        ee338 = 2 * (ee38 * ee38);
        ee340 = ee38 * ee196;
        ee344 = 2 * (ee196 * ee40);
        ee346 = R_pow(ee198, 2);
        ee349 = ee340/ee198;
        ee351 = (ee39 * ee338 - ee340)/ee198 - ee197 * ee344/ee346 -  ee349 - ee349;
        ee353 = ee201 * ee48;
        ee360 = ee204 * ee57;
        ee362 = ee61 * ee225;
        ee363 = ee362 + ee360;
        ee366 = ee201 * ee52;
        ee371 = 4 * (ee351 * ee46 - ee366 - (ee42 * ee215 + ee366));
        ee374 = ee201 * ee44;
        ee383 = ee201 * ee42;
        ee396 = 5 * (ee351 * ee59 - ee353 - (ee42 * ee210 + ee353)) *  ee69 - ee360 - ee363 + (ee36 * (70 * ee371 - 315 * (3 *  (ee351 * ee50 - ee374 - (ee42 * ee221 + ee374))) + 540 *  (2 * (ee351 * ee35 - ee383 - (ee42 * ee201 + ee383))) -  420 * ee351) - ee362 - ee363);
        ee401 = ee105 * ee240;
        ee404 = 2 * (ee100 * ee100);
        ee406 = ee232 * ee235;
        ee410 = 2 * (ee235 * ee102);
        ee412 = R_pow(ee237, 2);
        ee415 = ee406/ee237;
        ee420 = ee242 * ee104 - ee401 - (ee401 + ee90 * ((ee101 *  ee404 - ee406)/ee237 - ee236 * ee410/ee412 - ee415 -  ee415));
        ee424 = ee336 - (ee396 * ee90 + ee335) + (ee336 + (ee334 +  ee93 * ee420));
        ee426 = ee228 * ee118;
        ee428 = ee71 * ee253;
        ee429 = ee428 + ee426;
        ee432 = R_pow(ee111, ee249 - 1);
        ee436 = ee432 * (ee249 * ee116) * ee251 * ee117;
        ee442 = ee124 * ee268;
        ee443 = ee259 * ee132;
        ee444 = ee442 - ee443;
        ee448 = ee444 - (630 * ee371 * ee126 + ee443);
        ee456 = ee201 * ee264;
        ee461 = 4 * (ee42 * (3 * (ee42 * (2 * (ee42 * ee125)) -  ee201 * ee262)) - ee456 - (ee351 * ee130 + ee456));
        ee468 = ee271 * ee38/ee40;
        ee470 = ee273 * ee196;
        ee472 = ee468 - ee470/ee198;
        ee483 = ee279 * ee118;
        ee485 = ee139 * ee253;
        ee486 = ee483 - ee485;
        ee494 = R_pow(ee111,  -(ee286 + 1));
        ee500 = ee71 * ee291;
        ee502 = ee228 * ee150;
        ee503 = ee500 - ee502;
        ee518 = ee303 * ee38/ee40;
        ee520 = ee305 * ee196;
        ee522 = ee518 - ee520/ee198;
        ee533 = ee311 * ee105;
        ee535 = ee166 * ee242;
        ee536 = ee533 + ee535;
        ee544 = ee245 * ee176;
        ee546 = ee107 * ee323;
        ee547 = ee544 - ee546;
        ee560 = ee327 * ee179;
        ee569 = R_pow(ee330, 2);
        ee574 = ee32 - ee24;
        ee578 = (ee24 - ee574 * ee74/ee77) * ee83;
        ee579 = ee578/ee86;
        ee581 = ee574 * ee83/ee77;
        ee582 = ee84 * ee581;
        ee584 = ee579 + ee582/ee102;
        ee585 = ee90 * ee584;
        ee587 = ee578 * ee100;
        ee588 = ee587/ee102;
        ee589 = ee581 * ee86;
        ee590 = 2 * ee589;
        ee591 = ee101 * ee590;
        ee594 = ee97 * ee581;
        ee595 = ee594/ee102;
        ee596 = ee588 + ee591/ee237 + ee595;
        ee598 = ee585 * ee104 - ee90 * ee596;
        ee600 = ee24/ee34;
        ee601 = ee26 * ee574;
        ee603 = ee600 + ee601/ee40;
        ee605 = 2 * (ee603 * ee35);
        ee608 = 4 * (ee603 * ee46);
        ee611 = 3 * (ee603 * ee50);
        ee616 = 540 * ee605 + (70 * ee608 - 315 * ee611) - 420 *  ee603;
        ee619 = ee24 * ee38;
        ee620 = ee619/ee40;
        ee621 = ee574 * ee34;
        ee622 = 2 * ee621;
        ee623 = ee39 * ee622;
        ee626 = ee574/ee40;
        ee627 = ee620 + ee623/ee198 + ee626;
        ee630 = 5 * (ee42 * ee608 + ee627 * ee59);
        ee636 = 4 * (ee42 * ee611 + ee627 * ee46);
        ee641 = 3 * (ee42 * ee605 + ee627 * ee50);
        ee647 = 2 * (ee42 * ee603 + ee627 * ee35);
        ee651 = 70 * ee636 - 315 * ee641 + 540 * ee647 - 420 * ee627;
        ee654 = 5 * (ee603 * ee59);
        ee657 = ee61 * ee616 + ee630 * ee69 + (ee36 * ee651 + ee654 *  ee57);
        ee659 = ee71 * ee598 - ee657 * ee105;
        ee661 = 2 * (ee38 * ee574);
        ee667 = 2 * (ee622 * ee40);
        ee672 = ee38 * ee622/ee198;
        ee674 = (ee39 * ee661 - ee619 * ee196)/ee198 - ee197 * ee667/ee346 -  ee672 - ee672;
        ee688 = ee61 * ee651 + ee630 * ee57;
        ee697 = 4 * (ee674 * ee46 - ee201 * ee611 - (ee42 * ee641 +  ee627 * ee52));
        ee725 = 5 * (ee674 * ee59 - ee201 * ee608 - (ee42 * ee636 +  ee627 * ee48)) * ee69 - ee204 * ee616 - ee688 + (ee36 *  (70 * ee697 - 315 * (3 * (ee674 * ee50 - ee201 * ee605 -  (ee42 * ee647 + ee627 * ee44))) + 540 * (2 * (ee674 *  ee35 - ee201 * ee603 - (ee42 * ee627 + ee627 * ee42))) -  420 * ee674) - ee654 * ee225 - ee688);
        ee732 = ee36 * ee616 + ee654 * ee69;
        ee739 = 2 * (ee100 * ee581);
        ee745 = 2 * (ee590 * ee102);
        ee750 = ee232 * ee590/ee237;
        ee755 = ee598 * ee104 - ee105 * ee596 - (ee585 * ee240 +  ee90 * ((ee101 * ee739 - ee587 * ee235)/ee237 - ee236 *  ee745/ee412 - ee750 - ee750));
        ee759 = ee659 - (ee725 * ee90 + ee228 * ee585) + (ee659 +  (ee732 * ee242 + ee93 * ee755));
        ee761 = ee109 * ee14;
        ee762 = R_pow(ee14, 2);
        ee763 = ee761/ee762;
        ee764 = ee113 * ee763;
        ee765 = ee115 * ee764;
        ee768 = ee23 * ee14;
        ee769 = ee768/ee762;
        ee770 = ee113 * ee769;
        ee772 = ee114 * ee763;
        ee773 = ee250 * ee772;
        ee775 = ee115 * ee770 + ee773 * ee117;
        ee778 = ee71 * ee775 + ee657 * ee118;
        ee784 = ee432 * (ee249 * ee763);
        ee788 = ee252 * ee770 + (ee250 * (ee114 * ee769) + ee784 *  ee251) * ee117;
        ee796 = 3 * (ee603 * ee262);
        ee800 = 4 * (ee42 * ee796 - ee627 * ee130);
        ee802 = 630 * ee636;
        ee804 = ee124 * ee800 - ee802 * ee132;
        ee808 = 4 * (ee603 * ee130);
        ee811 = ee804 - (630 * ee697 * ee126 + ee259 * ee808);
        ee812 = 630 * ee608;
        ee813 = ee812 * ee268;
        ee815 = 2 * (ee603 * ee125);
        ee827 = 4 * (ee42 * (3 * (ee42 * ee815 - ee627 * ee262)) -  ee627 * ee264 - (ee674 * ee130 + ee201 * ee796));
        ee838 = ee124 * ee808 - ee802 * ee126;
        ee839 = ee812 * ee132;
        ee842 = ee838 + (ee839 + ee129 * ee800);
        ee847 = ee842 * ee38/ee40 - ee273 * ee622/ee198;
        ee849 = ee812 * ee126;
        ee851 = ee849 + ee129 * ee808;
        ee852 = ee851 * ee38;
        ee866 = ee134 * ee574;
        ee870 = ee137 * ee622;
        ee873 = ee842/ee34 - ee866/ee40 - (ee852/ee40 - ee870/ee198);
        ee876 = ee873 * ee118 - ee139 * ee775;
        ee879 = ee136 * ee574;
        ee881 = ee851/ee34 - ee879/ee40;
        ee888 = ee494 * (ee286 * ee763);
        ee894 = ee145 * ee769;
        ee900 = ee146 * ee763;
        ee901 = ee288 * ee900;
        ee904 = ee901 * ee149 - ee148 * ee894;
        ee907 = ee71 * ee904 - ee657 * ee150;
        ee910 = ee145 * ee763;
        ee911 = ee148 * ee910;
        ee931 = ee159 * ee800 - ee839 - ee838;
        ee936 = ee931 * ee38/ee40 - ee305 * ee622/ee198;
        ee939 = ee159 * ee808 - ee849;
        ee940 = ee939 * ee38;
        ee954 = ee161 * ee574;
        ee958 = ee164 * ee622;
        ee961 = ee931/ee34 - ee954/ee40 - (ee940/ee40 - ee958/ee198);
        ee964 = ee961 * ee105 + ee166 * ee598;
        ee967 = ee163 * ee574;
        ee969 = ee939/ee34 - ee967/ee40;
        ee977 = ee83 * ee581;
        ee978 = ee977/ee102;
        ee987 = ee71 * ee585 - ee657 * ee90 + (ee732 * ee105 + ee93 *  ee598);
        ee989 = ee175 * ee590;
        ee990 = ee989/ee237;
        ee992 = ee987 * ee176 - ee107 * ee990;
        ee996 = ee732 * ee90 + ee93 * ee585;
        ee1014 = ee92 * ee911 - ee732 * ee153;
        ee1017 = ee1014 * ee116 - ee181 * ee769;
        ee1027 = ee881 * ee120 - ee141 * ee765 + ee1017/ee7 + (ee969 *  ee90 + ee168 * ee585) + (ee996 * ee172 - ee174 * ee978);
        ee1044 = ee92 * ee904 - ee732 * ee150 - (ee71 * ee911 -  ee657 * ee153);
        ee1047 = ee1044 * ee116 - ee155 * ee769;
        ee1065 = ee873 * ee120 - ee139 * ee765 - (ee881 * ee118 -  ee141 * ee775) + ee1047/ee7 + (ee961 * ee90 + ee166 *  ee585 + (ee969 * ee105 + ee168 * ee598)) + (ee987 * ee172 -  ee107 * ee978 - (ee996 * ee176 - ee174 * ee990));
        ee1071 = 2 * (ee1027 * ee188);
        ee1077 = 1.5 * ee4;
        ee1078 = R_pow(ee5, 2);
        ee1079 = ee1077/ee1078;
        ee1080 = ee73 * ee1079;
        ee1081 = ee9 * ee1080;
        ee1082 = ee72 * ee1079;
        ee1083 = ee12 * ee1082;
        ee1085 = (ee1081 - ee1083) * ee14;
        ee1087 = ee82 * ee1079;
        ee1088 = ee19 * ee1087;
        ee1089 = ee81 * ee1079;
        ee1090 = ee22 * ee1089;
        ee1091 = ee1088 - ee1090;
        ee1092 = ee15 * ee1091;
        ee1093 = R_pow(ee23, 2);
        ee1095 = ee1085/ee23 - ee1092/ee1093;
        ee1096 = ee1095/ee34;
        ee1097 = ee76 * ee1079;
        ee1098 = ee29 * ee1097;
        ee1100 = (ee1098 - ee1083) * ee14;
        ee1102 = ee31 * ee1091;
        ee1105 = ee1100/ee23 - ee1102/ee1093 - ee1095;
        ee1106 = ee26 * ee1105;
        ee1108 = ee1096 + ee1106/ee40;
        ee1110 = 5 * (ee1108 * ee59);
        ee1112 = ee1105/ee40;
        ee1113 = ee1095 * ee38;
        ee1114 = ee1113/ee40;
        ee1115 = ee1105 * ee34;
        ee1116 = 2 * ee1115;
        ee1117 = ee39 * ee1116;
        ee1120 = ee1112 + (ee1114 + ee1117/ee198);
        ee1124 = 2 * (ee1120 * ee35 + ee42 * ee1108);
        ee1128 = 3 * (ee1108 * ee50);
        ee1131 = 4 * (ee1120 * ee46 + ee42 * ee1128);
        ee1135 = 2 * (ee1108 * ee35);
        ee1138 = 3 * (ee1120 * ee50 + ee42 * ee1135);
        ee1143 = 540 * ee1124 + (70 * ee1131 - 315 * ee1138) - 420 *  ee1120;
        ee1148 = 4 * (ee1108 * ee46);
        ee1151 = 5 * (ee1120 * ee59 + ee42 * ee1148);
        ee1159 = 70 * ee1148 - 315 * ee1128 + 540 * ee1135 - 420 *  ee1108;
        ee1162 = ee1110 * ee57 + ee36 * ee1143 + (ee1151 * ee69 +  ee61 * ee1159);
        ee1165 = ee1105 * ee83/ee77;
        ee1166 = ee97 * ee1165;
        ee1167 = ee1166/ee102;
        ee1171 = (ee1095 - ee1105 * ee74/ee77) * ee83;
        ee1172 = ee1171 * ee100;
        ee1173 = ee1172/ee102;
        ee1174 = ee1165 * ee86;
        ee1175 = 2 * ee1174;
        ee1176 = ee101 * ee1175;
        ee1179 = ee1167 + (ee1173 + ee1176/ee237);
        ee1181 = ee1171/ee86;
        ee1182 = ee84 * ee1165;
        ee1184 = ee1181 + ee1182/ee102;
        ee1185 = ee90 * ee1184;
        ee1187 = ee90 * ee1179 - ee1185 * ee104;
        ee1189 = ee1162 * ee105 + ee71 * ee1187;
        ee1192 = ee1151 * ee57 + ee61 * ee1143;
        ee1197 = ee38 * ee1116/ee198;
        ee1200 = 2 * (ee38 * ee1105);
        ee1205 = 2 * (ee1116 * ee40);
        ee1210 = ee1197 + ((ee1113 * ee196 - ee39 * ee1200)/ee198 +  ee197 * ee1205/ee346) + ee1197;
        ee1228 = 4 * (ee1120 * ee52 + ee42 * ee1138 + (ee1210 *  ee46 + ee201 * ee1128));
        ee1255 = ee1192 + (5 * (ee1120 * ee48 + ee42 * ee1131 +  (ee1210 * ee59 + ee201 * ee1148)) * ee69 + ee204 * ee1159) +  (ee1110 * ee225 + ee36 * (70 * ee1228 - 315 * (3 * (ee1120 *  ee44 + ee42 * ee1124 + (ee1210 * ee50 + ee201 * ee1135))) +  540 * (2 * (ee1120 * ee42 + ee42 * ee1120 + (ee1210 *  ee35 + ee201 * ee1108))) - 420 * ee1210) + ee1192);
        ee1264 = ee232 * ee1175/ee237;
        ee1267 = 2 * (ee100 * ee1165);
        ee1272 = 2 * (ee1175 * ee102);
        ee1281 = ee1187 * ee104 + ee105 * ee1179 - (ee90 * (ee1264 +  ((ee1172 * ee235 - ee101 * ee1267)/ee237 + ee236 * ee1272/ee412) +  ee1264) - ee1185 * ee240);
        ee1285 = ee1110 * ee69 + ee36 * ee1159;
        ee1289 = ee1189 - (ee1255 * ee90 - ee228 * ee1185) + (ee1189 +  (ee93 * ee1281 - ee1285 * ee242));
        ee1291 = log(ee111);
        ee1292 = R_pow(ee7, 2);
        ee1293 = ee1079/ee1292;
        ee1294 = ee1291 * ee1293;
        ee1296 = ee108 * ee1091;
        ee1298 = ee1083 + ee1296/ee14;
        ee1299 = ee114 * ee1298;
        ee1301 = ee115 * ee1294 - ee250 * ee1299;
        ee1303 = ee1293 * ee116;
        ee1304 = ee1091/ee14;
        ee1306 = ee1303 - ee113 * ee1304;
        ee1308 = ee1301 * ee117 + ee115 * ee1306;
        ee1310 = ee1162 * ee118 + ee71 * ee1308;
        ee1313 = ee113 * ee1298;
        ee1315 = ee120 * ee1294 - ee115 * ee1313;
        ee1323 = ee250 * ee1294 - ee432 * (ee249 * ee1298);
        ee1331 = (ee1323 * ee251 + ee250 * (ee1303 - ee114 * ee1304)) *  ee117 + ee252 * ee1306;
        ee1337 = 630 * ee1131;
        ee1341 = 3 * (ee1108 * ee262);
        ee1344 = 4 * (ee1120 * ee130 - ee42 * ee1341);
        ee1346 = ee1337 * ee132 + ee124 * ee1344;
        ee1350 = 4 * (ee1108 * ee130);
        ee1353 = ee1346 - (630 * ee1228 * ee126 - ee259 * ee1350);
        ee1357 = 2 * (ee1108 * ee125);
        ee1367 = 4 * (ee1120 * ee264 + ee42 * (3 * (ee1120 * ee262 -  ee42 * ee1357)) - (ee1210 * ee130 - ee201 * ee1341));
        ee1369 = 630 * ee1148;
        ee1370 = ee1369 * ee268;
        ee1380 = ee1337 * ee126 - ee124 * ee1350;
        ee1382 = ee1369 * ee132;
        ee1384 = ee1380 + (ee129 * ee1344 - ee1382);
        ee1389 = ee1384 * ee38/ee40 + ee273 * ee1116/ee198;
        ee1393 = ee1369 * ee126;
        ee1394 = ee129 * ee1350 + ee1393;
        ee1395 = ee1394 * ee38;
        ee1408 = ee134 * ee1105;
        ee1412 = ee137 * ee1116;
        ee1415 = ee1384/ee34 + ee1408/ee40 + (ee1395/ee40 - ee1412/ee198);
        ee1418 = ee1415 * ee118 + ee139 * ee1308;
        ee1422 = ee136 * ee1105;
        ee1424 = ee1394/ee34 - ee1422/ee40;
        ee1433 = ee494 * (ee286 * ee1298) + ee288 * ee1294;
        ee1441 = ee145 * ee1304 + ee1303;
        ee1447 = ee146 * ee1298;
        ee1450 = ee288 * ee1447 + ee148 * ee1294;
        ee1453 = ee1450 * ee149 - ee148 * ee1441;
        ee1455 = ee1162 * ee150 + ee71 * ee1453;
        ee1458 = ee145 * ee1298;
        ee1461 = ee148 * ee1458 + ee153 * ee1294;
        ee1484 = ee1382 + ee159 * ee1344 - ee1380;
        ee1489 = ee1484 * ee38/ee40 + ee305 * ee1116/ee198;
        ee1492 = ee1393 - ee159 * ee1350;
        ee1493 = ee1492 * ee38;
        ee1507 = ee161 * ee1105;
        ee1511 = ee164 * ee1116;
        ee1514 = ee1484/ee34 + ee1507/ee40 - (ee1493/ee40 + ee1511/ee198);
        ee1517 = ee1514 * ee105 + ee166 * ee1187;
        ee1520 = ee163 * ee1105;
        ee1522 = ee1492/ee34 + ee1520/ee40;
        ee1530 = ee83 * ee1165;
        ee1531 = ee1530/ee102;
        ee1540 = ee1162 * ee90 - ee71 * ee1185 + (ee93 * ee1187 -  ee1285 * ee105);
        ee1542 = ee175 * ee1175;
        ee1543 = ee1542/ee237;
        ee1545 = ee1540 * ee176 + ee107 * ee1543;
        ee1555 = ee93 * ee1185 + ee1285 * ee90;
        ee1567 = ee1285 * ee153 + ee92 * ee1461;
        ee1570 = ee1567 * ee116 - ee181 * ee1304;
        ee1572 = ee182 * ee1079;
        ee1583 = ee141 * ee1315 - ee1424 * ee120 + (ee1570/ee7 -  ee1572/ee1292) + (ee1522 * ee90 - ee168 * ee1185) + (ee174 *  ee1531 - ee1555 * ee172);
        ee1600 = ee1285 * ee150 + ee92 * ee1453 - (ee1162 * ee153 +  ee71 * ee1461);
        ee1603 = ee1600 * ee116 - ee155 * ee1304;
        ee1605 = ee156 * ee1079;
        ee1624 = ee1415 * ee120 + ee139 * ee1315 - (ee141 * ee1308 -  ee1424 * ee118) + (ee1603/ee7 - ee1605/ee1292) + (ee1514 *  ee90 - ee166 * ee1185 + (ee1522 * ee105 + ee168 * ee1187)) +  (ee1540 * ee172 + ee107 * ee1531 - (ee174 * ee1543 -  ee1555 * ee176));
        ee1630 = 2 * (ee1583 * ee188);
        ee1645 = ee179 * ee1027;
        ee1650 = ee578 * ee581;
        ee1660 = ee585 * ee584 + ee90 * (ee579 - ee1650/ee102 +  ((ee582 - ee1650)/ee102 - ee582 * ee590/ee237));
        ee1662 = ee657 * ee585;
        ee1664 = ee24 * ee574;
        ee1672 = ee600 - ee1664/ee40 + ((ee601 - ee1664)/ee40 -  ee601 * ee622/ee198);
        ee1676 = 2 * (ee1672 * ee35 - ee603 * ee603);
        ee1681 = 4 * (ee1672 * ee46 - ee603 * ee611);
        ee1686 = 3 * (ee1672 * ee50 - ee603 * ee605);
        ee1691 = 540 * ee1676 + (70 * ee1681 - 315 * ee1686) - 420 *  ee1672;
        ee1693 = ee630 * ee616;
        ee1696 = ee627 * ee608;
        ee1698 = ee619 * ee622;
        ee1703 = 2 * (ee621 + ee574 * ee574);
        ee1714 = ee620 - ee1698/ee198 + ((ee39 * ee1703 - ee1698)/ee198 -  ee623 * ee667/ee346) + (ee626 - ee574 * ee622/ee198);
        ee1723 = ee627 * ee611;
        ee1728 = 4 * (ee42 * ee1686 - ee1723 + (ee1714 * ee46 -  ee1723));
        ee1731 = ee627 * ee605;
        ee1740 = ee627 * ee603;
        ee1751 = ee654 * ee651;
        ee1756 = 5 * (ee1672 * ee59 - ee603 * ee608);
        ee1760 = ee61 * ee1691 - ee1693 + (5 * (ee42 * ee1681 -  ee1696 + (ee1714 * ee59 - ee1696)) * ee69 - ee1693) +  (ee36 * (70 * ee1728 - 315 * (3 * (ee42 * ee1676 - ee1731 +  (ee1714 * ee50 - ee1731))) + 540 * (2 * (ee42 * ee1672 -  ee1740 + (ee1714 * ee35 - ee1740))) - 420 * ee1714) -  ee1751 + (ee1756 * ee57 - ee1751));
        ee1765 = ee654 * ee616;
        ee1769 = ee36 * ee1691 - ee1765 + (ee1756 * ee69 - ee1765);
        ee1771 = ee732 * ee598;
        ee1774 = ee585 * ee596;
        ee1776 = ee587 * ee590;
        ee1781 = 2 * (ee589 + ee581 * ee581);
        ee1795 = ee1660 * ee104 - ee1774 - (ee1774 + ee90 * (ee588 -  ee1776/ee237 + ((ee101 * ee1781 - ee1776)/ee237 - ee591 *  ee745/ee412) + (ee595 - ee594 * ee590/ee237)));
        ee1799 = ee71 * ee1660 - ee1662 - (ee1760 * ee90 + ee1662) +  (ee1769 * ee105 + ee1771 + (ee1771 + ee93 * ee1795));
        ee1801 = 2 * (ee14 * ee14);
        ee1803 = R_pow(ee762, 2);
        ee1805 = ee763 - ee761 * ee1801/ee1803;
        ee1809 = ee115 * (ee113 * ee1805) - ee773 * ee764;
        ee1811 = ee657 * ee765;
        ee1818 = ee769 - ee768 * ee1801/ee1803;
        ee1821 = ee773 * ee770;
        ee1829 = ee115 * (ee113 * ee1818) - ee1821 + ((ee250 * (ee114 *  ee1805) - ee784 * ee772) * ee117 - ee1821);
        ee1831 = ee732 * ee775;
        ee1841 = 4 * (ee1672 * ee130 + ee603 * ee796);
        ee1843 = ee802 * ee808;
        ee1848 = ee124 * ee1841 - ee1843 - (630 * ee1728 * ee126 +  ee1843);
        ee1849 = 630 * ee1681;
        ee1851 = ee812 * ee800;
        ee1852 = ee1849 * ee132 + ee1851;
        ee1858 = ee627 * ee796;
        ee1863 = 4 * (ee42 * (3 * (ee1672 * ee262 + ee603 * ee815)) -  ee1858 - (ee1714 * ee130 + ee1858));
        ee1869 = ee842 * ee574;
        ee1879 = ee812 * ee808;
        ee1880 = ee1849 * ee126 + ee1879;
        ee1883 = ee1880 + (ee1879 + ee129 * ee1841);
        ee1886 = ee852 * ee622;
        ee1898 = ee873 * ee765;
        ee1904 = ee851 * ee574;
        ee1912 = ee1883/ee34 - ee1904/ee40 - ((ee1904 + ee879)/ee40 -  ee879 * ee622/ee198);
        ee1914 = ee881 * ee775;
        ee1925 = ee901 * ee894;
        ee1932 = ee732 * ee904;
        ee1940 = ee901 * ee910 + ee148 * (ee145 * ee1805);
        ee1942 = ee657 * ee911;
        ee1949 = ee1044 * ee769;
        ee1961 = ee931 * ee574;
        ee1972 = ee159 * ee1841 - ee1879 - ee1880;
        ee1975 = ee940 * ee622;
        ee1987 = ee961 * ee585;
        ee1993 = ee939 * ee574;
        ee2001 = ee1972/ee34 - ee1993/ee40 - ((ee1993 + ee967)/ee40 -  ee967 * ee622/ee198);
        ee2003 = ee969 * ee598;
        ee2011 = ee987 * ee978;
        ee2015 = ee978 - ee977 * ee590/ee237;
        ee2020 = ee732 * ee585;
        ee2024 = ee1769 * ee90 + ee2020 + (ee2020 + ee93 * ee1660);
        ee2026 = ee996 * ee990;
        ee2039 = ee1065 * ee1027;
        ee2043 = ee881 * ee765;
        ee2049 = ee732 * ee911;
        ee2055 = ee1014 * ee769;
        ee2063 = ee969 * ee585;
        ee2070 = ee996 * ee978;
        ee2097 = ee90 * ((ee1171 * ee581 - ee1182)/ee102 + ee582 *  ee1175/ee237 - (ee1181 - ee578 * ee1165/ee102)) - ee1185 *  ee584;
        ee2110 = (ee1095 * ee574 - ee1106)/ee40 + ee601 * ee1116/ee198 -  (ee1096 - ee24 * ee1105/ee40);
        ee2114 = 2 * (ee2110 * ee35 + ee603 * ee1108);
        ee2119 = 4 * (ee2110 * ee46 + ee603 * ee1128);
        ee2124 = 3 * (ee2110 * ee50 + ee603 * ee1135);
        ee2129 = 540 * ee2114 + (70 * ee2119 - 315 * ee2124) - 420 *  ee2110;
        ee2138 = 2 * (ee574 * ee1105 + ee1115);
        ee2152 = (ee1113 * ee622 - ee39 * ee2138)/ee198 + ee623 *  ee1205/ee346 - (ee1114 - ee619 * ee1116/ee198) - (ee1112 -  ee574 * ee1116/ee198);
        ee2170 = 4 * (ee1120 * ee611 + ee42 * ee2124 + (ee2152 *  ee46 + ee627 * ee1128));
        ee2199 = 5 * (ee2110 * ee59 + ee603 * ee1148);
        ee2204 = ee1151 * ee616 + ee61 * ee2129 + (5 * (ee1120 *  ee608 + ee42 * ee2119 + (ee2152 * ee59 + ee627 * ee1148)) *  ee69 + ee630 * ee1159) + (ee1110 * ee651 + ee36 * (70 *  ee2170 - 315 * (3 * (ee1120 * ee605 + ee42 * ee2114 +  (ee2152 * ee50 + ee627 * ee1135))) + 540 * (2 * (ee1120 *  ee603 + ee42 * ee2110 + (ee2152 * ee35 + ee627 * ee1108))) -  420 * ee2152) + (ee2199 * ee57 + ee654 * ee1143));
        ee2215 = ee1110 * ee616 + ee36 * ee2129 + (ee2199 * ee69 +  ee654 * ee1159);
        ee2225 = 2 * (ee581 * ee1165 + ee1174);
        ee2243 = ee2097 * ee104 + ee585 * ee1179 - (ee90 * ((ee1172 *  ee590 - ee101 * ee2225)/ee237 + ee591 * ee1272/ee412 -  (ee1173 - ee587 * ee1175/ee237) - (ee1167 - ee594 * ee1175/ee237)) -  ee1185 * ee596);
        ee2248 = ee1162 * ee585 + ee71 * ee2097 - (ee2204 * ee90 -  ee657 * ee1185) + (ee2215 * ee105 + ee732 * ee1187 +  (ee93 * ee2243 - ee1285 * ee598));
        ee2251 = ee1293 * ee763;
        ee2253 = ee1296 * ee14/ee762;
        ee2257 = ee1301 * ee764 + ee115 * (ee2251 - ee113 * ee2253);
        ee2266 = ee1293 * ee769;
        ee2268 = ee1091 * ee14/ee762;
        ee2281 = ee1301 * ee770 + ee115 * (ee2266 - ee113 * ee2268) +  ((ee1323 * ee772 + ee250 * (ee2251 - ee114 * ee2253)) *  ee117 + ee773 * ee1306);
        ee2294 = 4 * (ee2110 * ee130 - ee603 * ee1341);
        ee2301 = ee1337 * ee808 + ee124 * ee2294 - (630 * ee2170 *  ee126 - ee802 * ee1350);
        ee2302 = 630 * ee2119;
        ee2305 = ee2302 * ee132 + ee812 * ee1344;
        ee2317 = 4 * (ee1120 * ee796 + ee42 * (3 * (ee2110 * ee262 -  ee603 * ee1357)) - (ee2152 * ee130 - ee627 * ee1341));
        ee2319 = ee1369 * ee800;
        ee2336 = ee2302 * ee126 - ee812 * ee1350;
        ee2338 = ee1369 * ee808;
        ee2340 = ee2336 + (ee129 * ee2294 - ee2338);
        ee2372 = ee2340/ee34 + ee851 * ee1105/ee40 + ((ee1422 +  ee1394 * ee574)/ee40 - ee879 * ee1116/ee198);
        ee2407 = ee1450 * ee910 - ee148 * (ee145 * ee2253 + ee2251);
        ee2444 = ee2338 + ee159 * ee2294 - ee2336;
        ee2476 = ee2444/ee34 + ee939 * ee1105/ee40 - ((ee1492 *  ee574 - ee1520)/ee40 + ee967 * ee1116/ee198);
        ee2492 = ee1531 - ee977 * ee1175/ee237;
        ee2502 = ee2215 * ee90 - ee732 * ee1185 + (ee93 * ee2097 -  ee1285 * ee585);
        ee2581 = ee179 * ee1583;
        ee2585 = ee1095 * ee1105;
        ee2591 = ee1079 - ee1077 * (2 * (ee4 * ee5))/R_pow(ee1078, 2);
        ee2599 = ee12 * (ee72 * ee2591) + ee1083 * ee1082;
        ee2603 = ee1100 * ee1091;
        ee2614 = ee19 * (ee82 * ee2591) + ee1088 * ee1087 - (ee22 *  (ee81 * ee2591) + ee1090 * ee1089);
        ee2619 = 2 * (ee1091 * ee23);
        ee2621 = R_pow(ee1093, 2);
        ee2632 = ee1085 * ee1091;
        ee2641 = (ee9 * (ee73 * ee2591) + ee1081 * ee1080 - ee2599) *  ee14/ee23 - ee2632/ee1093 - ((ee15 * ee2614 + ee2632)/ee1093 -  ee1092 * ee2619/ee2621);
        ee2642 = (ee29 * (ee76 * ee2591) + ee1098 * ee1097 - ee2599) *  ee14/ee23 - ee2603/ee1093 - ((ee31 * ee2614 + ee2603)/ee1093 -  ee1102 * ee2619/ee2621) - ee2641;
        ee2652 = (ee2585 - ee26 * ee2642)/ee40 + ee1106 * ee1116/ee198 -  (ee2641/ee34 - ee2585/ee40);
        ee2656 = 5 * (ee2652 * ee59 + ee1108 * ee1148);
        ee2658 = ee1110 * ee1143;
        ee2660 = ee1113 * ee1116;
        ee2664 = 2 * (ee1105 * ee1105 + ee2642 * ee34);
        ee2680 = (ee2660 - ee39 * ee2664)/ee198 + ee1117 * ee1205/ee346 -  (ee2641 * ee38/ee40 - ee2660/ee198) - (ee2642/ee40 -  ee1105 * ee1116/ee198);
        ee2682 = ee1120 * ee1108;
        ee2690 = ee1120 * ee1128;
        ee2695 = 3 * (ee2652 * ee50 + ee1108 * ee1135);
        ee2699 = 4 * (ee2680 * ee46 + ee2690 + (ee2690 + ee42 *  ee2695));
        ee2702 = ee1120 * ee1135;
        ee2707 = 2 * (ee2652 * ee35 + ee1108 * ee1108);
        ee2721 = ee1120 * ee1148;
        ee2726 = 4 * (ee2652 * ee46 + ee1108 * ee1128);
        ee2732 = ee1151 * ee1159;
        ee2740 = 70 * ee2726 - 315 * ee2695 + 540 * ee2707 - 420 *  ee2652;
        ee2744 = ee2656 * ee57 + ee2658 + (ee2658 + ee36 * (540 *  (2 * (ee2680 * ee35 + ee2682 + (ee2682 + ee42 * ee2652))) +  (70 * ee2699 - 315 * (3 * (ee2680 * ee50 + ee2702 + (ee2702 +  ee42 * ee2707)))) - 420 * ee2680)) + (5 * (ee2680 *  ee59 + ee2721 + (ee2721 + ee42 * ee2726)) * ee69 + ee2732 +  (ee2732 + ee61 * ee2740));
        ee2746 = ee1162 * ee1185;
        ee2748 = ee1171 * ee1165;
        ee2750 = ee2642 * ee83/ee77;
        ee2760 = (ee2641 - ee2642 * ee74/ee77) * ee83;
        ee2767 = ee90 * ((ee2748 - ee84 * ee2750)/ee102 + ee1182 *  ee1175/ee237 - (ee2760/ee86 - ee2748/ee102)) - ee1185 *  ee1184;
        ee2771 = ee1172 * ee1175;
        ee2775 = 2 * (ee1165 * ee1165 + ee2750 * ee86);
        ee2794 = ee1185 * ee1179;
        ee2798 = ee90 * ((ee2771 - ee101 * ee2775)/ee237 + ee1176 *  ee1272/ee412 - (ee2760 * ee100/ee102 - ee2771/ee237) -  (ee97 * ee2750/ee102 - ee1166 * ee1175/ee237)) - ee2794 -  (ee2767 * ee104 + ee2794);
        ee2800 = ee1285 * ee1187;
        ee2803 = ee1110 * ee1159;
        ee2807 = ee2656 * ee69 + ee2803 + (ee2803 + ee36 * ee2740);
        ee2811 = ee2744 * ee90 - ee2746 - (ee2746 + ee71 * ee2767) +  (ee93 * ee2798 - ee2800 - (ee2807 * ee105 + ee2800));
        ee2813 = ee1285 * ee1308;
        ee2818 = 2 * (ee1079 * ee7);
        ee2820 = R_pow(ee1292, 2);
        ee2822 = ee2591/ee1292 + ee1079 * ee2818/ee2820;
        ee2826 = ee1291 * ee2822 + ee1298/ee111 * ee1293;
        ee2830 = ee1293 * ee1298;
        ee2833 = ee108 * ee2614/ee14 + ee2599;
        ee2840 = ee1301 * ee1306;
        ee2842 = ee1293 * ee1304;
        ee2844 = ee2842 + ee2822 * ee116;
        ee2845 = ee2614/ee14;
        ee2851 = (ee1301 * ee1294 - ee115 * ee2826 - (ee1323 * ee1299 +  ee250 * (ee2830 - ee114 * ee2833))) * ee117 + ee2840 +  (ee2840 - ee115 * (ee2844 + (ee2842 - ee113 * ee2845)));
        ee2856 = ee1162 * ee1315;
        ee2866 = ee1315 * ee1294 - ee120 * ee2826 - (ee1301 * ee1313 +  ee115 * (ee2830 - ee113 * ee2833));
        ee2874 = ee1337 * ee1350;
        ee2879 = 4 * (ee2652 * ee130 - ee1108 * ee1341);
        ee2882 = 630 * ee2699 * ee126 - ee2874 - (ee2874 + ee124 *  ee2879);
        ee2884 = ee1120 * ee1341;
        ee2893 = 4 * (ee2680 * ee130 - ee2884 - (ee2884 + ee42 *  (3 * (ee2652 * ee262 - ee1108 * ee1357))));
        ee2895 = ee1369 * ee1344;
        ee2897 = 630 * ee2726;
        ee2899 = ee2897 * ee132 + ee2895;
        ee2903 = ee1384 * ee1105;
        ee2914 = ee1369 * ee1350;
        ee2917 = ee2897 * ee126 - ee2914;
        ee2918 = ee129 * ee2879 - ee2914 + ee2917;
        ee2921 = ee1395 * ee1116;
        ee2933 = ee1415 * ee1315;
        ee2939 = ee1424 * ee1308;
        ee2942 = ee1394 * ee1105;
        ee2951 = ee2918/ee34 + ee2942/ee40 + ((ee136 * ee2642 +  ee2942)/ee40 - ee1422 * ee1116/ee198);
        ee2957 = ee1285 * ee1453;
        ee2969 = ee1450 * ee1441;
        ee2981 = ee1162 * ee1461;
        ee2991 = ee1450 * ee1458 - ee148 * (ee145 * ee2833 + ee2830) +  (ee1461 * ee1294 - ee153 * ee2826);
        ee2997 = ee1600 * ee1304;
        ee3003 = ee1603 * ee1079;
        ee3019 = ee1484 * ee1105;
        ee3031 = ee2917 - (ee2914 + ee159 * ee2879);
        ee3034 = ee1493 * ee1116;
        ee3046 = ee1514 * ee1185;
        ee3052 = ee1492 * ee1105;
        ee3061 = ee3031/ee34 + ee3052/ee40 + ((ee3052 - ee163 *  ee2642)/ee40 + ee1520 * ee1116/ee198);
        ee3063 = ee1522 * ee1187;
        ee3071 = ee1540 * ee1531;
        ee3077 = ee83 * ee2750/ee102 - ee1530 * ee1175/ee237;
        ee3087 = ee1555 * ee1543;
        ee3090 = ee1285 * ee1185;
        ee3094 = ee93 * ee2767 - ee3090 + (ee2807 * ee90 - ee3090);
        ee3101 = ee1624 * ee1583;
        ee3105 = ee1424 * ee1315;
        ee3111 = ee1285 * ee1461;
        ee3117 = ee1567 * ee1304;
        ee3123 = ee1570 * ee1079;
        ee3135 = ee1522 * ee1185;
        ee3142 = ee1555 * ee1531;
        
        out(j, 0) =  ee424 + (ee396 * ee120 - ee426 - ee429 - (ee429 + (ee92 * ee436 + ee428))) - ((((ee448 + (ee444 + (ee442 + ee129 * ee461)))/ee34 - ee468 - ee472 - (ee472 - ((ee470 + ee137 * ee338)/ee198 - ee276 * ee344/ee346))) * ee120 - ee483 - ee486 - (ee486 - (ee485 - ee141 * ee436)) + (ee92 * (ee494 * (ee286 * ee116) * ee289 * ee149) - ee500 - ee503 - (ee503 - (ee396 * ee153 + ee502))) * ee116/ee7 + (((ee159 * ee461 - ee442 - ee444 - ee448)/ee34 - ee518 - ee522 - (ee522 - ((ee520 + ee164 * ee338)/ee198 - ee308 * ee344/ee346))) * ee90 + ee533 + ee536 + (ee536 + (ee535 + ee168 * ee420))) + (ee424 * ee172 - ee544 - ee547 - (ee547 - (ee546 + ee174 * (ee175 * ee404/ee237 - ee322 * ee410/ee412)))))/ee188 - ee560/ee330 - ((ee560 + ee179 * ee327)/ee330 - ee329 * (2 * (ee179 * ee188))/ee569));
          out(j, 1) =  ee759 + (ee725 * ee120 - ee228 * ee765 - ee778 - (ee778 + (ee92 * ee788 + ee732 * ee253))) - ((((ee811 + (ee804 + (ee813 + ee129 * ee827)))/ee34 - ee271 * ee574/ee40 - ee847 - (ee847 - ((ee852 * ee196 + ee137 * ee661)/ee198 - ee276 * ee667/ee346))) * ee120 - ee279 * ee765 - ee876 - (ee876 - (ee881 * ee253 - ee141 * ee788)) + ((ee92 * ((ee888 * ee289 - ee288 * (ee146 * ee769)) * ee149 - ee290 * ee894) - ee732 * ee291 - ee907 - (ee907 - (ee725 * ee153 + ee228 * ee911))) * ee116 - ee297 * ee769)/ee7 + (((ee159 * ee827 - ee813 - ee804 - ee811)/ee34 - ee303 * ee574/ee40 - ee936 - (ee936 - ((ee940 * ee196 + ee164 * ee661)/ee198 - ee308 * ee667/ee346))) * ee90 + ee311 * ee585 + ee964 + (ee964 + (ee969 * ee242 + ee168 * ee755))) + (ee759 * ee172 - ee245 * ee978 - ee992 - (ee992 - (ee996 * ee323 + ee174 * (ee175 * ee739/ee237 - ee322 * ee745/ee412)))))/ee188 - ee327 * ee1027/ee330 - ((ee1065 * ee179 + ee179 * ee1065)/ee330 - ee329 * ee1071/ee569));
          out(j, 2) =  ee1289 + (ee1310 + (ee1255 * ee120 + ee228 * ee1315) + (ee1285 * ee253 + ee92 * ee1331 + ee1310)) - ((((ee1353 + (ee1346 + (ee129 * ee1367 - ee1370)))/ee34 + ee271 * ee1105/ee40 - ee1389 - (ee1389 + ((ee137 * ee1200 + ee1395 * ee196)/ee198 - ee276 * ee1205/ee346))) * ee120 + ee279 * ee1315 - ee1418 - (ee1418 - (ee141 * ee1331 - ee1424 * ee253)) + (((ee1285 * ee291 + ee92 * ((ee1433 * ee289 - ee288 * (ee146 * ee1304 + ee1303)) * ee149 - ee290 * ee1441) - ee1455 - (ee1455 - (ee1255 * ee153 + ee228 * ee1461))) * ee116 - ee297 * ee1304)/ee7 - ee298 * ee1079/ee1292) + (((ee1370 + ee159 * ee1367 - ee1346 - ee1353)/ee34 + ee303 * ee1105/ee40 - ee1489 - (ee1489 - ((ee1493 * ee196 - ee164 * ee1200)/ee198 + ee308 * ee1205/ee346))) * ee90 - ee311 * ee1185 + ee1517 + (ee1517 + (ee1522 * ee242 + ee168 * ee1281))) + (ee1289 * ee172 + ee245 * ee1531 - ee1545 - (ee1545 + (ee174 * (ee175 * ee1267/ee237 - ee322 * ee1272/ee412) + ee1555 * ee323))))/ee188 - ee327 * ee1583/ee330 - ((ee1624 * ee179 + ee179 * ee1624)/ee330 - ee329 * ee1630/ee569));
          out(j, 3) =  ee1799 + (ee71 * ee1809 - ee1811 + (ee1760 * ee120 - ee1811) + (ee92 * ee1829 - ee1831 + (ee1769 * ee118 - ee1831))) - ((((ee1848 + (ee1852 + (ee1851 + ee129 * ee1863)))/ee34 - ee1869/ee40 - ((ee1869 + ee866)/ee40 - ee866 * ee622/ee198) - (ee1883 * ee38/ee40 - ee1886/ee198 - ((ee1886 + ee137 * ee1703)/ee198 - ee870 * ee667/ee346))) * ee120 - ee1898 - (ee1898 + ee139 * ee1809) - (ee1912 * ee118 - ee1914 - (ee1914 + ee141 * ee1829)) + ((ee92 * ((ee888 * ee900 + ee288 * (ee146 * ee1805)) * ee149 - ee1925 - (ee1925 + ee148 * (ee145 * ee1818))) - ee1932 - (ee1769 * ee150 + ee1932) - (ee71 * ee1940 - ee1942 - (ee1760 * ee153 + ee1942))) * ee116 - ee1949 - (ee1949 + ee155 * ee1818))/ee7 + (((ee159 * ee1863 - ee1851 - ee1852 - ee1848)/ee34 - ee1961/ee40 - ((ee1961 + ee954)/ee40 - ee954 * ee622/ee198) - (ee1972 * ee38/ee40 - ee1975/ee198 - ((ee1975 + ee164 * ee1703)/ee198 - ee958 * ee667/ee346))) * ee90 + ee1987 + (ee1987 + ee166 * ee1660) + (ee2001 * ee105 + ee2003 + (ee2003 + ee168 * ee1795))) + (ee1799 * ee172 - ee2011 - (ee2011 + ee107 * ee2015) - (ee2024 * ee176 - ee2026 - (ee2026 + ee174 * (ee175 * ee1781/ee237 - ee989 * ee745/ee412)))))/ee188 - ee2039/ee330 - ((ee2039 + ee179 * (ee1912 * ee120 - ee2043 - (ee2043 + ee141 * ee1809) + ((ee92 * ee1940 - ee2049 - (ee1769 * ee153 + ee2049)) * ee116 - ee2055 - (ee2055 + ee181 * ee1818))/ee7 + (ee2001 * ee90 + ee2063 + (ee2063 + ee168 * ee1660)) + (ee2024 * ee172 - ee2070 - (ee2070 + ee174 * ee2015))))/ee330 - ee1645 * ee1071/ee569));
          out(j, 4) =  ee2248 + (ee1162 * ee765 + ee71 * ee2257 + (ee2204 * ee120 + ee657 * ee1315) + (ee1285 * ee775 + ee92 * ee2281 + (ee2215 * ee118 + ee732 * ee1308))) - ((((ee2301 + (ee2305 + (ee129 * ee2317 - ee2319)))/ee34 + ee842 * ee1105/ee40 - ((ee1384 * ee574 - ee1408)/ee40 + ee866 * ee1116/ee198) - (ee2340 * ee38/ee40 + ee852 * ee1116/ee198 + ((ee137 * ee2138 + ee1395 * ee622)/ee198 - ee870 * ee1205/ee346))) * ee120 + ee873 * ee1315 - (ee1415 * ee765 + ee139 * ee2257) - (ee2372 * ee118 + ee881 * ee1308 - (ee141 * ee2281 - ee1424 * ee775)) + (((ee1285 * ee904 + ee92 * ((ee1433 * ee900 - ee288 * (ee146 * ee2253 + ee2251)) * ee149 - ee901 * ee1441 - (ee1450 * ee894 - ee148 * (ee145 * ee2268 + ee2266))) - (ee2215 * ee150 + ee732 * ee1453) - (ee1162 * ee911 + ee71 * ee2407 - (ee2204 * ee153 + ee657 * ee1461))) * ee116 - ee1044 * ee1304 - (ee1600 * ee769 - ee155 * ee2268))/ee7 - ee1047 * ee1079/ee1292) + (((ee2319 + ee159 * ee2317 - ee2305 - ee2301)/ee34 + ee931 * ee1105/ee40 - ((ee1484 * ee574 - ee1507)/ee40 + ee954 * ee1116/ee198) - (ee2444 * ee38/ee40 + ee940 * ee1116/ee198 - ((ee1493 * ee622 - ee164 * ee2138)/ee198 + ee958 * ee1205/ee346))) * ee90 - ee961 * ee1185 + (ee1514 * ee585 + ee166 * ee2097) + (ee2476 * ee105 + ee969 * ee1187 + (ee1522 * ee598 + ee168 * ee2243))) + (ee2248 * ee172 + ee987 * ee1531 - (ee1540 * ee978 - ee107 * ee2492) - (ee2502 * ee176 + ee996 * ee1543 + (ee174 * (ee175 * ee2225/ee237 - ee989 * ee1272/ee412) + ee1555 * ee990))))/ee188 - ee1065 * ee1583/ee330 - ((ee1624 * ee1027 + ee179 * (ee2372 * ee120 + ee881 * ee1315 - (ee141 * ee2257 - ee1424 * ee765) + (((ee1285 * ee911 + ee92 * ee2407 - (ee2215 * ee153 + ee732 * ee1461)) * ee116 - ee1014 * ee1304 - (ee1567 * ee769 - ee181 * ee2268))/ee7 - ee1017 * ee1079/ee1292) + (ee2476 * ee90 - ee969 * ee1185 + (ee1522 * ee585 + ee168 * ee2097)) + (ee2502 * ee172 + ee996 * ee1531 + (ee174 * ee2492 + ee1555 * ee978))))/ee330 - ee1645 * ee1630/ee569));
          out(j, 5) =  ee2811 - (ee2807 * ee118 + ee2813 + (ee2813 + ee92 * ee2851) + (ee2744 * ee120 + ee2856 + (ee2856 + ee71 * ee2866))) - ((((ee2882 + (ee129 * ee2893 - ee2895 - ee2899))/ee34 + ee2903/ee40 + ((ee2903 - ee134 * ee2642)/ee40 + ee1408 * ee1116/ee198) + (ee2918 * ee38/ee40 + ee2921/ee198 + ((ee137 * ee2664 + ee2921)/ee198 - ee1412 * ee1205/ee346))) * ee120 + ee2933 + (ee2933 + ee139 * ee2866) - (ee141 * ee2851 - ee2939 - (ee2951 * ee118 + ee2939)) + (((ee2807 * ee150 + ee2957 + (ee2957 + ee92 * ((ee1433 * ee1447 - ee288 * (ee146 * ee2833 + ee2830) + (ee1450 * ee1294 - ee148 * ee2826)) * ee149 - ee2969 - (ee2969 - ee148 * (ee2844 + (ee145 * ee2845 + ee2842))))) - (ee2744 * ee153 + ee2981 + (ee2981 + ee71 * ee2991))) * ee116 - ee2997 - (ee2997 - ee155 * ee2845))/ee7 - ee3003/ee1292 - ((ee3003 - ee156 * ee2591)/ee1292 - ee1605 * ee2818/ee2820)) + (((ee2899 + (ee2895 + ee159 * ee2893) - ee2882)/ee34 + ee3019/ee40 + ((ee3019 - ee161 * ee2642)/ee40 + ee1507 * ee1116/ee198) - (ee3031 * ee38/ee40 + ee3034/ee198 + ((ee3034 - ee164 * ee2664)/ee198 + ee1511 * ee1205/ee346))) * ee90 - ee3046 - (ee3046 + ee166 * ee2767) + (ee3061 * ee105 + ee3063 + (ee3063 + ee168 * ee2798))) + (ee2811 * ee172 + ee3071 + (ee3071 - ee107 * ee3077) + (ee174 * (ee175 * ee2775/ee237 - ee1542 * ee1272/ee412) + ee3087 + (ee3094 * ee176 + ee3087))))/ee188 - ee3101/ee330 - ((ee3101 + ee179 * (ee141 * ee2866 - ee3105 - (ee2951 * ee120 + ee3105) + (((ee2807 * ee153 + ee3111 + (ee3111 + ee92 * ee2991)) * ee116 - ee3117 - (ee3117 - ee181 * ee2845))/ee7 - ee3123/ee1292 - ((ee3123 - ee182 * ee2591)/ee1292 - ee1572 * ee2818/ee2820)) + (ee3061 * ee90 - ee3135 - (ee3135 + ee168 * ee2767)) - (ee174 * ee3077 + ee3142 + (ee3094 * ee172 + ee3142))))/ee330 - ee2581 * ee1630/ee569));

          ee2 = -log(psuba);
          ee4 = exp(-txi);
          ee5 = 1 + ee4;
          ee7 = 1.5/ee5 - 0.5;
          ee8 = -ee7;
          ee9 = R_pow(ee2, ee8);
          ee11 = -log(alpha);
          ee12 = R_pow(ee11, ee8);
          ee14 = exp(lsbeta);
          ee15 = (ee9 - ee12) * ee14;
          ee18 = -log(1 - hbeta);
          ee19 = R_pow(ee18, ee8);
          ee21 = -log(hbeta);
          ee22 = R_pow(ee21, ee8);
          ee23 = ee19 - ee22;
          ee24 = ee15/ee23;
          ee25 = qalpha + ee24;
          ee26 = y - ee25;
          ee28 = -log(psubb);
          ee29 = R_pow(ee28, ee8);
          ee31 = (ee29 - ee12) * ee14;
          ee32 = ee31/ee23;
          ee34 = qalpha + ee32 - ee25;
          ee35 = ee26/ee34;
          ee36 = R_pow(ee35, 5);
          ee37 = ee24/ee34;
          ee38 = ee32 - ee24;
          ee39 = ee26 * ee38;
          ee40 = R_pow(ee34, 2);
          ee42 = ee37 + ee39/ee40;
          ee44 = 2 * (ee42 * ee35);
          ee46 = R_pow(ee35, 3);
          ee48 = 4 * (ee42 * ee46);
          ee50 = R_pow(ee35, 2);
          ee52 = 3 * (ee42 * ee50);
          ee57 = 540 * ee44 + (70 * ee48 - 315 * ee52) - 420 * ee42;
          ee59 = R_pow(ee35, 4);
          ee61 = 5 * (ee42 * ee59);
          ee69 = 70 * ee59 - 315 * ee46 + 540 * ee50 - 420 * ee35 +  126;
          ee71 = ee36 * ee57 + ee61 * ee69;
          ee72 = log(ee11);
          ee73 = log(ee2);
          ee74 = ee72 - ee73;
          ee76 = log(ee28);
          ee77 = ee73 - ee76;
          ee81 = log(ee21);
          ee82 = log(ee18);
          ee83 = ee81 - ee82;
          ee84 = (y - (ee25 - ee34 * ee74/ee77)) * ee83;
          ee86 = ee34 * ee83/ee77;
          ee90 = exp(-(ee84/ee86 - ee72));
          ee92 = ee36 * ee69;
          ee93 = 1 - ee92;
          ee97 = (ee24 - ee38 * ee74/ee77) * ee83;
          ee98 = ee97/ee86;
          ee100 = ee38 * ee83/ee77;
          ee101 = ee84 * ee100;
          ee102 = R_pow(ee86, 2);
          ee104 = ee98 + ee101/ee102;
          ee105 = ee90 * ee104;
          ee107 = ee71 * ee90 + ee93 * ee105;
          ee108 = y - qalpha;
          ee109 = ee108 * ee23;
          ee111 = ee109/ee14 + ee12;
          ee113 = -1/ee7;
          ee114 = ee113 - 1;
          ee115 = R_pow(ee111, ee114);
          ee116 = ee109 * ee14;
          ee117 = R_pow(ee14, 2);
          ee118 = ee116/ee117;
          ee119 = ee113 * ee118;
          ee120 = ee115 * ee119;
          ee122 = R_pow(ee111, ee113);
          ee126 = 630 * ee48;
          ee127 = 1 - ee35;
          ee128 = R_pow(ee127, 4);
          ee129 = ee126 * ee128;
          ee131 = -630 * ee59;
          ee132 = R_pow(ee127, 3);
          ee134 = 4 * (ee42 * ee132);
          ee136 = ee129 + ee131 * ee134;
          ee138 = ee131 * ee128;
          ee139 = ee138 * ee38;
          ee141 = ee136/ee34 - ee139/ee40;
          ee143 = ee138/ee34;
          ee147 = 1 + 1/ee7;
          ee148 = ee147 + 1;
          ee150 = R_pow(ee111,  - ee148);
          ee151 = ee147 * ee118;
          ee152 = ee150 * ee151;
          ee155 = R_pow(ee111,  - ee147);
          ee157 = ee92 * ee152 - ee71 * ee155;
          ee158 = ee23/ee14;
          ee160 = ee92 * ee155;
          ee161 = ee23 * ee14;
          ee162 = ee161/ee117;
          ee164 = ee157 * ee158 - ee160 * ee162;
          ee167 = 630 * ee59;
          ee169 = ee167 * ee134 - ee129;
          ee171 = ee167 * ee128;
          ee172 = ee171 * ee38;
          ee174 = ee169/ee34 - ee172/ee40;
          ee176 = ee171/ee34;
          ee180 = ee83/ee86;
          ee182 = ee93 * ee90;
          ee183 = ee83 * ee100;
          ee184 = ee183/ee102;
          ee187 = ee141 * ee122 - ee143 * ee120 + ee164/ee7 + (ee174 *  ee90 + ee176 * ee105) + (ee107 * ee180 - ee182 * ee184);
          ee189 = ee160 * ee158;
          ee195 = ee143 * ee122 + ee189/ee7 + ee176 * ee90 + ee182 *  ee180;
          ee198 = ee24 * ee38;
          ee200 = ee37 - ee198/ee40;
          ee201 = ee39 - ee198;
          ee203 = ee38 * ee34;
          ee204 = 2 * ee203;
          ee205 = ee39 * ee204;
          ee206 = R_pow(ee40, 2);
          ee209 = ee200 + (ee201/ee40 - ee205/ee206);
          ee213 = 2 * (ee209 * ee35 - ee42 * ee42);
          ee218 = 4 * (ee209 * ee46 - ee42 * ee52);
          ee223 = 3 * (ee209 * ee50 - ee42 * ee44);
          ee228 = 540 * ee213 + (70 * ee218 - 315 * ee223) - 420 *  ee209;
          ee230 = ee61 * ee57;
          ee235 = 5 * (ee209 * ee59 - ee42 * ee48);
          ee238 = ee36 * ee228 - ee230 + (ee235 * ee69 - ee230);
          ee240 = ee71 * ee105;
          ee243 = ee97 * ee100;
          ee245 = ee98 - ee243/ee102;
          ee246 = ee101 - ee243;
          ee248 = ee100 * ee86;
          ee249 = 2 * ee248;
          ee250 = ee101 * ee249;
          ee251 = R_pow(ee102, 2);
          ee254 = ee245 + (ee246/ee102 - ee250/ee251);
          ee256 = ee105 * ee104 + ee90 * ee254;
          ee259 = ee238 * ee90 + ee240 + (ee240 + ee93 * ee256);
          ee260 = ee14 * ee14;
          ee261 = 2 * ee260;
          ee262 = ee116 * ee261;
          ee263 = R_pow(ee117, 2);
          ee265 = ee118 - ee262/ee263;
          ee266 = ee113 * ee265;
          ee268 = ee114 - 1;
          ee269 = R_pow(ee111, ee268);
          ee270 = ee114 * ee118;
          ee271 = ee269 * ee270;
          ee273 = ee115 * ee266 - ee271 * ee119;
          ee275 = ee71 * ee120;
          ee281 = 630 * ee218;
          ee283 = ee126 * ee134;
          ee284 = ee281 * ee128 + ee283;
          ee286 = R_pow(ee127, 2);
          ee288 = 3 * (ee42 * ee286);
          ee291 = 4 * (ee209 * ee132 + ee42 * ee288);
          ee294 = ee284 + (ee283 + ee131 * ee291);
          ee296 = ee136 * ee38;
          ee299 = ee296 + ee139;
          ee301 = ee139 * ee204;
          ee304 = ee294/ee34 - ee296/ee40 - (ee299/ee40 - ee301/ee206);
          ee306 = ee141 * ee120;
          ee311 = ee148 + 1;
          ee313 = R_pow(ee111,  - ee311);
          ee314 = ee148 * ee118;
          ee315 = ee313 * ee314;
          ee317 = ee147 * ee265;
          ee319 = ee315 * ee151 + ee150 * ee317;
          ee321 = ee71 * ee152;
          ee325 = ee92 * ee319 - ee321 - (ee238 * ee155 + ee321);
          ee327 = ee157 * ee162;
          ee329 = ee161 * ee261;
          ee331 = ee162 - ee329/ee263;
          ee334 = ee325 * ee158 - ee327 - (ee327 + ee160 * ee331);
          ee339 = ee167 * ee291 - ee283 - ee284;
          ee341 = ee169 * ee38;
          ee344 = ee341 + ee172;
          ee346 = ee172 * ee204;
          ee349 = ee339/ee34 - ee341/ee40 - (ee344/ee40 - ee346/ee206);
          ee351 = ee174 * ee105;
          ee358 = ee107 * ee184;
          ee360 = ee183 * ee249;
          ee362 = ee184 - ee360/ee251;
          ee366 = ee304 * ee122 - ee306 - (ee306 + ee143 * ee273) +  ee334/ee7 + (ee349 * ee90 + ee351 + (ee351 + ee176 *  ee256)) + (ee259 * ee180 - ee358 - (ee358 + ee182 * ee362));
          ee368 = ee187 * ee187;
          ee369 = R_pow(ee195, 2);
          ee373 = ee198 + ee198;
          ee381 = ee201 * ee204;
          ee386 = 2 * (ee203 + ee38 * ee38);
          ee391 = 2 * (ee204 * ee40);
          ee393 = R_pow(ee206, 2);
          ee397 = ee200 - (ee373/ee40 - ee198 * ee204/ee206) + ((ee201 -  ee373)/ee40 - ee381/ee206 - ((ee381 + ee39 * ee386)/ee206 -  ee205 * ee391/ee393));
          ee399 = ee209 * ee42;
          ee407 = ee209 * ee52;
          ee412 = 4 * (ee397 * ee46 - ee407 - (ee407 + ee42 * ee223));
          ee415 = ee209 * ee44;
          ee427 = ee61 * ee228;
          ee429 = ee235 * ee57;
          ee430 = ee429 + ee427;
          ee433 = ee209 * ee48;
          ee442 = ee36 * (540 * (2 * (ee397 * ee35 - ee399 - (ee399 +  ee42 * ee209))) + (70 * ee412 - 315 * (3 * (ee397 * ee50 -  ee415 - (ee415 + ee42 * ee213)))) - 420 * ee397) - ee427 -  ee430 + (5 * (ee397 * ee59 - ee433 - (ee433 + ee42 *  ee218)) * ee69 - ee429 - ee430);
          ee444 = ee238 * ee105;
          ee446 = ee71 * ee256;
          ee447 = ee444 + ee446;
          ee450 = ee105 * ee254;
          ee452 = ee243 + ee243;
          ee460 = ee246 * ee249;
          ee465 = 2 * (ee248 + ee100 * ee100);
          ee470 = 2 * (ee249 * ee102);
          ee472 = R_pow(ee251, 2);
          ee479 = ee256 * ee104 + ee450 + (ee450 + ee90 * (ee245 -  (ee452/ee102 - ee243 * ee249/ee251) + ((ee246 - ee452)/ee102 -  ee460/ee251 - ((ee460 + ee101 * ee465)/ee251 - ee250 *  ee470/ee472))));
          ee483 = ee442 * ee90 + ee444 + ee447 + (ee447 + (ee446 +  ee93 * ee479));
          ee485 = 2 * (ee260 + ee260);
          ee490 = 2 * (ee261 * ee117);
          ee492 = R_pow(ee263, 2);
          ee495 = ee265 - ((ee262 + ee116 * ee485)/ee263 - ee262 *  ee490/ee492);
          ee498 = ee271 * ee266;
          ee503 = R_pow(ee111, ee268 - 1);
          ee510 = ee115 * (ee113 * ee495) - ee498 - ((ee269 * (ee114 *  ee265) - ee503 * (ee268 * ee118) * ee270) * ee119 + ee498);
          ee512 = ee71 * ee273;
          ee514 = ee238 * ee120;
          ee515 = ee514 + ee512;
          ee524 = ee281 * ee134;
          ee526 = ee126 * ee291;
          ee527 = ee524 + ee526;
          ee528 = 630 * ee412 * ee128 + ee524 + ee527;
          ee530 = ee209 * ee288;
          ee541 = 4 * (ee397 * ee132 + ee530 + (ee530 + ee42 * (3 *  (ee209 * ee286 + ee42 * (2 * (ee42 * ee127))))));
          ee547 = ee294 * ee38;
          ee550 = ee547 + ee296;
          ee558 = ee299 * ee204;
          ee570 = ee304 * ee120;
          ee572 = ee141 * ee273;
          ee573 = ee570 + ee572;
          ee581 = R_pow(ee111,  -(ee311 + 1));
          ee589 = ee315 * ee317;
          ee596 = ee71 * ee319;
          ee598 = ee238 * ee152;
          ee599 = ee598 + ee596;
          ee606 = ee325 * ee162;
          ee608 = ee157 * ee331;
          ee609 = ee606 + ee608;
          ee629 = ee339 * ee38;
          ee632 = ee629 + ee341;
          ee640 = ee344 * ee204;
          ee652 = ee349 * ee105;
          ee654 = ee174 * ee256;
          ee655 = ee652 + ee654;
          ee663 = ee259 * ee184;
          ee665 = ee107 * ee362;
          ee666 = ee663 + ee665;
          ee681 = ee366 * ee187;
          ee690 = R_pow(ee369, 2);
          ee695 = 1.5 * ee4;
          ee696 = R_pow(ee5, 2);
          ee697 = ee695/ee696;
          ee698 = ee73 * ee697;
          ee699 = ee9 * ee698;
          ee700 = ee72 * ee697;
          ee701 = ee12 * ee700;
          ee703 = (ee699 - ee701) * ee14;
          ee705 = ee82 * ee697;
          ee706 = ee19 * ee705;
          ee707 = ee81 * ee697;
          ee708 = ee22 * ee707;
          ee709 = ee706 - ee708;
          ee710 = ee15 * ee709;
          ee711 = R_pow(ee23, 2);
          ee713 = ee703/ee23 - ee710/ee711;
          ee714 = ee713/ee34;
          ee715 = ee76 * ee697;
          ee716 = ee29 * ee715;
          ee718 = (ee716 - ee701) * ee14;
          ee720 = ee31 * ee709;
          ee723 = ee718/ee23 - ee720/ee711 - ee713;
          ee724 = ee26 * ee723;
          ee726 = ee714 + ee724/ee40;
          ee728 = 5 * (ee726 * ee59);
          ee730 = ee713 * ee38;
          ee731 = ee730 - ee724;
          ee732 = ee24 * ee723;
          ee733 = ee732 + ee730;
          ee736 = ee723 * ee34;
          ee737 = 2 * ee736;
          ee744 = 2 * (ee38 * ee723 + ee736);
          ee749 = 2 * (ee737 * ee40);
          ee755 = ee714 - ee732/ee40;
          ee761 = (ee731 + ee733)/ee40 + ee201 * ee737/ee206 - ((ee731 *  ee204 - ee39 * ee744)/ee206 + ee205 * ee749/ee393) -  (ee755 - (ee733/ee40 - ee198 * ee737/ee206));
          ee766 = ee39 * ee737;
          ee769 = ee731/ee40 + ee766/ee206 - ee755;
          ee778 = 3 * (ee726 * ee50);
          ee784 = 2 * (ee726 * ee35);
          ee787 = 3 * (ee769 * ee50 + ee42 * ee784);
          ee791 = 4 * (ee761 * ee46 + ee209 * ee778 - (ee769 * ee52 +  ee42 * ee787));
          ee800 = 2 * (ee769 * ee35 + ee42 * ee726);
          ee814 = 4 * (ee726 * ee46);
          ee817 = 5 * (ee769 * ee59 + ee42 * ee814);
          ee823 = 4 * (ee769 * ee46 + ee42 * ee778);
          ee829 = 540 * ee800 + (70 * ee823 - 315 * ee787) - 420 *  ee769;
          ee831 = ee817 * ee57 + ee61 * ee829;
          ee848 = 70 * ee814 - 315 * ee778 + 540 * ee784 - 420 * ee726;
          ee852 = ee728 * ee228 + ee36 * (540 * (2 * (ee761 * ee35 +  ee209 * ee726 - (ee769 * ee42 + ee42 * ee769))) + (70 *  ee791 - 315 * (3 * (ee761 * ee50 + ee209 * ee784 - (ee769 *  ee44 + ee42 * ee800)))) - 420 * ee761) - ee831 + (5 *  (ee761 * ee59 + ee209 * ee814 - (ee769 * ee48 + ee42 *  ee823)) * ee69 + ee235 * ee848 - ee831);
          ee857 = (ee713 - ee723 * ee74/ee77) * ee83;
          ee858 = ee857/ee86;
          ee860 = ee723 * ee83/ee77;
          ee861 = ee84 * ee860;
          ee863 = ee858 + ee861/ee102;
          ee864 = ee90 * ee863;
          ee873 = ee728 * ee57 + ee36 * ee829 + (ee817 * ee69 + ee61 *  ee848);
          ee875 = ee857 * ee100;
          ee876 = ee875 - ee861;
          ee878 = ee860 * ee86;
          ee879 = 2 * ee878;
          ee880 = ee101 * ee879;
          ee883 = ee97 * ee860;
          ee885 = ee858 - ee883/ee102;
          ee886 = ee876/ee102 + ee880/ee251 - ee885;
          ee889 = ee90 * ee886 - ee864 * ee104;
          ee891 = ee873 * ee105 + ee71 * ee889;
          ee896 = ee883 + ee875;
          ee905 = 2 * (ee100 * ee860 + ee878);
          ee910 = 2 * (ee879 * ee102);
          ee924 = ee889 * ee104 + ee105 * ee886 + (ee90 * ((ee876 +  ee896)/ee102 + ee246 * ee879/ee251 - ((ee876 * ee249 -  ee101 * ee905)/ee251 + ee250 * ee910/ee472) - (ee885 -  (ee896/ee102 - ee243 * ee879/ee251))) - ee864 * ee254);
          ee928 = ee728 * ee69 + ee36 * ee848;
          ee932 = ee852 * ee90 - ee238 * ee864 + ee891 + (ee891 +  (ee93 * ee924 - ee928 * ee256));
          ee934 = log(ee111);
          ee935 = R_pow(ee7, 2);
          ee936 = ee697/ee935;
          ee937 = ee934 * ee936;
          ee939 = ee108 * ee709;
          ee941 = ee701 + ee939/ee14;
          ee942 = ee114 * ee941;
          ee944 = ee115 * ee937 - ee269 * ee942;
          ee946 = ee936 * ee265;
          ee947 = ee939 * ee14;
          ee948 = ee947/ee117;
          ee951 = ee948 - ee947 * ee261/ee263;
          ee959 = ee269 * ee937 - ee503 * (ee268 * ee941);
          ee961 = ee936 * ee118;
          ee968 = ee961 - ee113 * ee948;
          ee971 = ee944 * ee266 + ee115 * (ee946 - ee113 * ee951) -  ((ee959 * ee270 + ee269 * (ee961 - ee114 * ee948)) *  ee119 + ee271 * ee968);
          ee977 = ee944 * ee119 + ee115 * ee968;
          ee979 = ee873 * ee120 + ee71 * ee977;
          ee983 = ee113 * ee941;
          ee985 = ee122 * ee937 - ee115 * ee983;
          ee994 = 4 * (ee726 * ee132);
          ee997 = 630 * ee823;
          ee1001 = 3 * (ee726 * ee286);
          ee1004 = 4 * (ee769 * ee132 - ee42 * ee1001);
          ee1006 = ee997 * ee134 + ee126 * ee1004;
          ee1007 = 630 * ee791 * ee128 - ee281 * ee994 + ee1006;
          ee1014 = 2 * (ee726 * ee127);
          ee1021 = 4 * (ee761 * ee132 - ee209 * ee1001 + (ee769 *  ee288 + ee42 * (3 * (ee769 * ee286 - ee42 * ee1014))));
          ee1023 = 630 * ee814;
          ee1024 = ee1023 * ee291;
          ee1034 = ee997 * ee128 - ee126 * ee994;
          ee1036 = ee1023 * ee134;
          ee1038 = ee1034 + (ee131 * ee1004 - ee1036);
          ee1040 = ee136 * ee723;
          ee1041 = ee1038 * ee38 - ee1040;
          ee1047 = ee138 * ee723;
          ee1049 = ee1023 * ee128;
          ee1050 = ee131 * ee994 + ee1049;
          ee1052 = ee1047 + ee1050 * ee38;
          ee1074 = ee139 * ee737;
          ee1077 = ee1038/ee34 + ee1040/ee40 + (ee1052/ee40 - ee1074/ee206);
          ee1080 = ee1077 * ee120 + ee141 * ee977;
          ee1085 = ee1050/ee34 - ee1047/ee40;
          ee1094 = ee581 * (ee311 * ee941) + ee313 * ee937;
          ee1102 = ee147 * ee948 + ee961;
          ee1105 = ee148 * ee941;
          ee1108 = ee313 * ee1105 + ee150 * ee937;
          ee1120 = ee1108 * ee151 - ee150 * ee1102;
          ee1122 = ee873 * ee152 + ee71 * ee1120;
          ee1125 = ee147 * ee941;
          ee1128 = ee150 * ee1125 + ee155 * ee937;
          ee1134 = ee709/ee14;
          ee1143 = ee928 * ee152 + ee92 * ee1120 - (ee873 * ee155 +  ee71 * ee1128);
          ee1145 = ee709 * ee14;
          ee1146 = ee1145/ee117;
          ee1148 = ee1143 * ee162 - ee157 * ee1146;
          ee1152 = ee928 * ee155 + ee92 * ee1128;
          ee1176 = ee1036 + ee167 * ee1004 - ee1034;
          ee1178 = ee169 * ee723;
          ee1179 = ee1176 * ee38 - ee1178;
          ee1186 = ee1049 - ee167 * ee994;
          ee1188 = ee171 * ee723;
          ee1189 = ee1186 * ee38 - ee1188;
          ee1211 = ee172 * ee737;
          ee1214 = ee1176/ee34 + ee1178/ee40 - (ee1189/ee40 + ee1211/ee206);
          ee1217 = ee1214 * ee105 + ee174 * ee889;
          ee1221 = ee1186/ee34 + ee1188/ee40;
          ee1229 = ee83 * ee860;
          ee1230 = ee1229/ee102;
          ee1239 = ee873 * ee90 - ee71 * ee864 + (ee93 * ee889 - ee928 *  ee105);
          ee1241 = ee183 * ee879;
          ee1243 = ee1230 - ee1241/ee251;
          ee1245 = ee1239 * ee184 - ee107 * ee1243;
          ee1258 = ee93 * ee864 + ee928 * ee90;
          ee1270 = ee1152 * ee158 - ee160 * ee1134;
          ee1272 = ee189 * ee697;
          ee1283 = ee143 * ee985 - ee1085 * ee122 + (ee1270/ee7 -  ee1272/ee935) + (ee1221 * ee90 - ee176 * ee864) + (ee182 *  ee1230 - ee1258 * ee180);
          ee1300 = ee1143 * ee158 - ee157 * ee1134 - (ee1152 * ee162 -  ee160 * ee1146);
          ee1302 = ee164 * ee697;
          ee1321 = ee1077 * ee122 + ee141 * ee985 - (ee143 * ee977 -  ee1085 * ee120) + (ee1300/ee7 - ee1302/ee935) + (ee1214 *  ee90 - ee174 * ee864 + (ee1221 * ee105 + ee176 * ee889)) +  (ee1239 * ee180 + ee107 * ee1230 + (ee182 * ee1243 +  ee1258 * ee184));
          ee1327 = 2 * (ee1283 * ee195);
          ee1342 = ee187 * ee1283;
          ee1346 = ee713 * ee723;
          ee1352 = ee697 - ee695 * (2 * (ee4 * ee5))/R_pow(ee696, 2);
          ee1360 = ee12 * (ee72 * ee1352) + ee701 * ee700;
          ee1364 = ee718 * ee709;
          ee1375 = ee19 * (ee82 * ee1352) + ee706 * ee705 - (ee22 *  (ee81 * ee1352) + ee708 * ee707);
          ee1380 = 2 * (ee709 * ee23);
          ee1382 = R_pow(ee711, 2);
          ee1393 = ee703 * ee709;
          ee1402 = (ee9 * (ee73 * ee1352) + ee699 * ee698 - ee1360) *  ee14/ee23 - ee1393/ee711 - ((ee15 * ee1375 + ee1393)/ee711 -  ee710 * ee1380/ee1382);
          ee1403 = (ee29 * (ee76 * ee1352) + ee716 * ee715 - ee1360) *  ee14/ee23 - ee1364/ee711 - ((ee31 * ee1375 + ee1364)/ee711 -  ee720 * ee1380/ee1382) - ee1402;
          ee1405 = ee1346 - ee26 * ee1403;
          ee1412 = ee1402/ee34 - ee1346/ee40;
          ee1413 = ee1405/ee40 + ee724 * ee737/ee206 - ee1412;
          ee1417 = 5 * (ee1413 * ee59 + ee726 * ee814);
          ee1419 = ee728 * ee829;
          ee1421 = ee731 * ee737;
          ee1425 = 2 * (ee723 * ee723 + ee1403 * ee34);
          ee1446 = (ee1421 - ee39 * ee1425)/ee206 + ee766 * ee749/ee393 -  ((ee1346 + ee1402 * ee38 + ee1405)/ee40 - ee1421/ee206) +  (ee1412 - ((ee24 * ee1403 + ee1346)/ee40 - ee732 * ee737/ee206));
          ee1448 = ee769 * ee726;
          ee1456 = ee769 * ee778;
          ee1461 = 3 * (ee1413 * ee50 + ee726 * ee784);
          ee1465 = 4 * (ee1446 * ee46 + ee1456 + (ee1456 + ee42 *  ee1461));
          ee1468 = ee769 * ee784;
          ee1473 = 2 * (ee1413 * ee35 + ee726 * ee726);
          ee1487 = ee769 * ee814;
          ee1492 = 4 * (ee1413 * ee46 + ee726 * ee778);
          ee1498 = ee817 * ee848;
          ee1506 = 70 * ee1492 - 315 * ee1461 + 540 * ee1473 - 420 *  ee1413;
          ee1510 = ee1417 * ee57 + ee1419 + (ee1419 + ee36 * (540 *  (2 * (ee1446 * ee35 + ee1448 + (ee1448 + ee42 * ee1413))) +  (70 * ee1465 - 315 * (3 * (ee1446 * ee50 + ee1468 + (ee1468 +  ee42 * ee1473)))) - 420 * ee1446)) + (5 * (ee1446 *  ee59 + ee1487 + (ee1487 + ee42 * ee1492)) * ee69 + ee1498 +  (ee1498 + ee61 * ee1506));
          ee1512 = ee873 * ee864;
          ee1514 = ee857 * ee860;
          ee1516 = ee1403 * ee83/ee77;
          ee1518 = ee1514 - ee84 * ee1516;
          ee1526 = (ee1402 - ee1403 * ee74/ee77) * ee83;
          ee1529 = ee1526/ee86 - ee1514/ee102;
          ee1533 = ee90 * (ee1518/ee102 + ee861 * ee879/ee251 - ee1529) -  ee864 * ee863;
          ee1537 = ee876 * ee879;
          ee1541 = 2 * (ee860 * ee860 + ee1516 * ee86);
          ee1564 = ee864 * ee886;
          ee1568 = ee90 * ((ee1537 - ee101 * ee1541)/ee251 + ee880 *  ee910/ee472 - ((ee1514 + ee1526 * ee100 + ee1518)/ee102 -  ee1537/ee251) + (ee1529 - ((ee97 * ee1516 + ee1514)/ee102 -  ee883 * ee879/ee251))) - ee1564 - (ee1533 * ee104 + ee1564);
          ee1570 = ee928 * ee889;
          ee1573 = ee728 * ee848;
          ee1577 = ee1417 * ee69 + ee1573 + (ee1573 + ee36 * ee1506);
          ee1581 = ee1510 * ee90 - ee1512 - (ee1512 + ee71 * ee1533) +  (ee93 * ee1568 - ee1570 - (ee1577 * ee105 + ee1570));
          ee1583 = ee928 * ee977;
          ee1588 = 2 * (ee697 * ee7);
          ee1590 = R_pow(ee935, 2);
          ee1592 = ee1352/ee935 + ee697 * ee1588/ee1590;
          ee1596 = ee934 * ee1592 + ee941/ee111 * ee936;
          ee1600 = ee936 * ee941;
          ee1601 = ee108 * ee1375;
          ee1603 = ee1601/ee14 + ee1360;
          ee1610 = ee944 * ee968;
          ee1612 = ee936 * ee948;
          ee1614 = ee1612 + ee1592 * ee118;
          ee1616 = ee1601 * ee14/ee117;
          ee1622 = (ee944 * ee937 - ee115 * ee1596 - (ee959 * ee942 +  ee269 * (ee1600 - ee114 * ee1603))) * ee119 + ee1610 +  (ee1610 - ee115 * (ee1614 + (ee1612 - ee113 * ee1616)));
          ee1627 = ee873 * ee985;
          ee1637 = ee985 * ee937 - ee122 * ee1596 - (ee944 * ee983 +  ee115 * (ee1600 - ee113 * ee1603));
          ee1645 = ee997 * ee994;
          ee1650 = 4 * (ee1413 * ee132 - ee726 * ee1001);
          ee1653 = 630 * ee1465 * ee128 - ee1645 - (ee1645 + ee126 *  ee1650);
          ee1655 = ee769 * ee1001;
          ee1664 = 4 * (ee1446 * ee132 - ee1655 - (ee1655 + ee42 *  (3 * (ee1413 * ee286 - ee726 * ee1014))));
          ee1666 = ee1023 * ee1004;
          ee1668 = 630 * ee1492;
          ee1670 = ee1668 * ee134 + ee1666;
          ee1674 = ee1038 * ee723;
          ee1685 = ee1023 * ee994;
          ee1688 = ee1668 * ee128 - ee1685;
          ee1689 = ee131 * ee1650 - ee1685 + ee1688;
          ee1691 = ee1050 * ee723;
          ee1694 = ee138 * ee1403 + ee1691;
          ee1697 = ee1052 * ee737;
          ee1709 = ee1077 * ee985;
          ee1715 = ee1085 * ee977;
          ee1724 = ee1689/ee34 + ee1691/ee40 + (ee1694/ee40 - ee1047 *  ee737/ee206);
          ee1730 = ee928 * ee1120;
          ee1742 = ee1108 * ee1102;
          ee1754 = ee873 * ee1128;
          ee1764 = ee1108 * ee1125 - ee150 * (ee147 * ee1603 + ee1600) +  (ee1128 * ee937 - ee155 * ee1596);
          ee1770 = ee1143 * ee1134;
          ee1772 = ee1375/ee14;
          ee1777 = ee928 * ee1128;
          ee1781 = ee1577 * ee155 + ee1777 + (ee1777 + ee92 * ee1764);
          ee1783 = ee1152 * ee1146;
          ee1792 = ee1300 * ee697;
          ee1808 = ee1176 * ee723;
          ee1820 = ee1688 - (ee1685 + ee167 * ee1650);
          ee1822 = ee1186 * ee723;
          ee1825 = ee1822 - ee171 * ee1403;
          ee1828 = ee1189 * ee737;
          ee1840 = ee1214 * ee864;
          ee1852 = ee1820/ee34 + ee1822/ee40 + (ee1825/ee40 + ee1188 *  ee737/ee206);
          ee1854 = ee1221 * ee889;
          ee1862 = ee1239 * ee1230;
          ee1866 = ee1229 * ee879;
          ee1868 = ee83 * ee1516/ee102 - ee1866/ee251;
          ee1873 = ee928 * ee864;
          ee1877 = ee93 * ee1533 - ee1873 + (ee1577 * ee90 - ee1873);
          ee1879 = ee1258 * ee1243;
          ee1894 = ee1321 * ee1283;
          ee1898 = ee1085 * ee985;
          ee1904 = ee1152 * ee1134;
          ee1910 = ee1270 * ee697;
          ee1922 = ee1221 * ee864;
          ee1929 = ee1258 * ee1230;
          
          out(j, 6) =  ee483 - (ee92 * ee510 - ee512 - ee515 + (ee442 * ee122 - ee514 - ee515)) - ((((ee528 + (ee527 + (ee526 + ee131 * ee541)))/ee34 - ee547/ee40 - (ee550/ee40 - ee296 * ee204/ee206) - ((ee550 + ee299)/ee40 - ee558/ee206 - ((ee558 + ee139 * ee386)/ee206 - ee301 * ee391/ee393))) * ee122 - ee570 - ee573 - (ee573 + (ee572 + ee143 * ee510)) + ((ee92 * ((ee581 * (ee311 * ee118) * ee314 + ee313 * (ee148 * ee265)) * ee151 + ee589 + (ee589 + ee150 * (ee147 * ee495))) - ee596 - ee599 - (ee442 * ee155 + ee598 + ee599)) * ee158 - ee606 - ee609 - (ee609 + (ee608 + ee160 * (ee331 - ((ee329 + ee161 * ee485)/ee263 - ee329 * ee490/ee492)))))/ee7 + (((ee167 * ee541 - ee526 - ee527 - ee528)/ee34 - ee629/ee40 - (ee632/ee40 - ee341 * ee204/ee206) - ((ee632 + ee344)/ee40 - ee640/ee206 - ((ee640 + ee172 * ee386)/ee206 - ee346 * ee391/ee393))) * ee90 + ee652 + ee655 + (ee655 + (ee654 + ee176 * ee479))) + (ee483 * ee180 - ee663 - ee666 - (ee666 + (ee665 + ee182 * (ee362 - ((ee360 + ee183 * ee465)/ee251 - ee360 * ee470/ee472))))))/ee195 - ee681/ee369 - ((ee681 + ee187 * ee366)/ee369 - ee368 * (2 * (ee187 * ee195))/ee690));
          out(j, 7) =  ee932 - (ee928 * ee273 + ee92 * ee971 - ee979 + (ee852 * ee122 + ee238 * ee985 - ee979)) - ((((ee1007 + (ee1006 + (ee131 * ee1021 - ee1024)))/ee34 + ee294 * ee723/ee40 - (ee1041/ee40 + ee296 * ee737/ee206) - ((ee1041 - ee1052)/ee40 + ee299 * ee737/ee206 + ((ee139 * ee744 + ee1052 * ee204)/ee206 - ee301 * ee749/ee393))) * ee122 + ee304 * ee985 - ee1080 - (ee1080 + (ee143 * ee971 - ee1085 * ee273)) + (((ee928 * ee319 + ee92 * ((ee1094 * ee314 - ee313 * (ee148 * ee948 + ee961)) * ee151 - ee315 * ee1102 + (ee1108 * ee317 - ee150 * (ee147 * ee951 + ee946))) - ee1122 - (ee852 * ee155 + ee238 * ee1128 + ee1122)) * ee158 - ee325 * ee1134 - ee1148 - (ee1148 + (ee1152 * ee331 - ee160 * (ee1146 - ee1145 * ee261/ee263))))/ee7 - ee334 * ee697/ee935) + (((ee1024 + ee167 * ee1021 - ee1006 - ee1007)/ee34 + ee339 * ee723/ee40 - (ee1179/ee40 + ee341 * ee737/ee206) - ((ee1179 + ee1189)/ee40 + ee344 * ee737/ee206 - ((ee1189 * ee204 - ee172 * ee744)/ee206 + ee346 * ee749/ee393))) * ee90 - ee349 * ee864 + ee1217 + (ee1217 + (ee1221 * ee256 + ee176 * ee924))) + (ee932 * ee180 + ee259 * ee1230 - ee1245 - (ee1245 - (ee182 * (ee1243 - ((ee183 * ee905 + ee1229 * ee249)/ee251 - ee360 * ee910/ee472)) + ee1258 * ee362))))/ee195 - ee366 * ee1283/ee369 - ((ee1321 * ee187 + ee187 * ee1321)/ee369 - ee368 * ee1327/ee690));
          out(j, 8) =  ee1581 - (ee1577 * ee120 + ee1583 + (ee1583 + ee92 * ee1622) + (ee1510 * ee122 + ee1627 + (ee1627 + ee71 * ee1637))) - ((((ee1653 + (ee131 * ee1664 - ee1666 - ee1670))/ee34 + ee1674/ee40 + ((ee1674 - ee136 * ee1403)/ee40 + ee1040 * ee737/ee206) + ((ee1689 * ee38 - ee1691 - ee1694)/ee40 + ee1697/ee206 + ((ee139 * ee1425 + ee1697)/ee206 - ee1074 * ee749/ee393))) * ee122 + ee1709 + (ee1709 + ee141 * ee1637) - (ee143 * ee1622 - ee1715 - (ee1724 * ee120 + ee1715)) + (((ee1577 * ee152 + ee1730 + (ee1730 + ee92 * ((ee1094 * ee1105 - ee313 * (ee148 * ee1603 + ee1600) + (ee1108 * ee937 - ee150 * ee1596)) * ee151 - ee1742 - (ee1742 - ee150 * (ee1614 + (ee147 * ee1616 + ee1612))))) - (ee1510 * ee155 + ee1754 + (ee1754 + ee71 * ee1764))) * ee158 - ee1770 - (ee1770 - ee157 * ee1772) - (ee1781 * ee162 - ee1783 - (ee1783 - ee160 * (ee1375 * ee14/ee117))))/ee7 - ee1792/ee935 - ((ee1792 - ee164 * ee1352)/ee935 - ee1302 * ee1588/ee1590)) + (((ee1670 + (ee1666 + ee167 * ee1664) - ee1653)/ee34 + ee1808/ee40 + ((ee1808 - ee169 * ee1403)/ee40 + ee1178 * ee737/ee206) - ((ee1820 * ee38 - ee1822 - ee1825)/ee40 + ee1828/ee206 + ((ee1828 - ee172 * ee1425)/ee206 + ee1211 * ee749/ee393))) * ee90 - ee1840 - (ee1840 + ee174 * ee1533) + (ee1852 * ee105 + ee1854 + (ee1854 + ee176 * ee1568))) + (ee1581 * ee180 + ee1862 + (ee1862 - ee107 * ee1868) + (ee1877 * ee184 - ee1879 - (ee182 * (ee1868 - ((ee183 * ee1541 + ee1866)/ee251 - ee1241 * ee910/ee472)) + ee1879))))/ee195 - ee1894/ee369 - ((ee1894 + ee187 * (ee143 * ee1637 - ee1898 - (ee1724 * ee122 + ee1898) + ((ee1781 * ee158 - ee1904 - (ee1904 - ee160 * ee1772))/ee7 - ee1910/ee935 - ((ee1910 - ee189 * ee1352)/ee935 - ee1272 * ee1588/ee1590)) + (ee1852 * ee90 - ee1922 - (ee1922 + ee176 * ee1533)) - (ee182 * ee1868 + ee1929 + (ee1877 * ee180 + ee1929))))/ee369 - ee1342 * ee1327/ee690));      
          
          ee2 = -log(psuba);
            ee4 = exp(-txi);
            ee5 = 1 + ee4;
          ee7 = 1.5/ee5 - 0.5;
          ee8 = -ee7;
          ee9 = R_pow(ee2, ee8);
          ee10 = log(ee2);
            ee11 = 1.5 * ee4;
          ee12 = R_pow(ee5, 2);
          ee13 = ee11/ee12;
          ee14 = ee10 * ee13;
          ee15 = ee9 * ee14;
          ee17 = -log(alpha);
            ee18 = R_pow(ee17, ee8);
          ee19 = log(ee17);
            ee20 = ee19 * ee13;
          ee21 = ee18 * ee20;
          ee23 = exp(lsbeta);
            ee24 = (ee15 - ee21) * ee23;
          ee27 = -log(1 - hbeta);
            ee28 = R_pow(ee27, ee8);
          ee30 = -log(hbeta);
            ee31 = R_pow(ee30, ee8);
          ee32 = ee28 - ee31;
          ee35 = (ee9 - ee18) * ee23;
          ee36 = log(ee27);
            ee37 = ee36 * ee13;
          ee38 = ee28 * ee37;
          ee39 = log(ee30);
            ee40 = ee39 * ee13;
          ee41 = ee31 * ee40;
          ee42 = ee38 - ee41;
          ee43 = ee35 * ee42;
          ee44 = R_pow(ee32, 2);
          ee46 = ee24/ee32 - ee43/ee44;
          ee48 = -log(psubb);
            ee49 = R_pow(ee48, ee8);
          ee50 = log(ee48);
            ee51 = ee50 * ee13;
          ee52 = ee49 * ee51;
          ee54 = (ee52 - ee21) * ee23;
          ee57 = (ee49 - ee18) * ee23;
          ee58 = ee57 * ee42;
          ee61 = ee54/ee32 - ee58/ee44 - ee46;
          ee62 = ee46 * ee61;
          ee64 = qalpha + ee35/ee32;
          ee65 = y - ee64;
          ee66 = ee4 * ee5;
          ee67 = 2 * ee66;
          ee68 = ee11 * ee67;
          ee69 = R_pow(ee12, 2);
          ee71 = ee13 - ee68/ee69;
          ee72 = ee50 * ee71;
          ee75 = ee49 * ee72 + ee52 * ee51;
          ee76 = ee19 * ee71;
          ee79 = ee18 * ee76 + ee21 * ee20;
          ee81 = (ee75 - ee79) * ee23;
          ee83 = ee54 * ee42;
          ee86 = ee36 * ee71;
          ee89 = ee28 * ee86 + ee38 * ee37;
          ee90 = ee39 * ee71;
          ee93 = ee31 * ee90 + ee41 * ee40;
          ee94 = ee89 - ee93;
          ee96 = ee57 * ee94 + ee83;
          ee99 = 2 * (ee42 * ee32);
            ee100 = ee58 * ee99;
          ee101 = R_pow(ee44, 2);
          ee105 = ee10 * ee71;
          ee108 = ee9 * ee105 + ee15 * ee14;
          ee110 = (ee108 - ee79) * ee23;
          ee112 = ee24 * ee42;
          ee116 = ee35 * ee94 + ee112;
          ee118 = ee43 * ee99;
          ee121 = ee110/ee32 - ee112/ee44 - (ee116/ee44 - ee118/ee101);
            ee122 = ee81/ee32 - ee83/ee44 - (ee96/ee44 - ee100/ee101) - 
              ee121;
              ee124 = ee62 - ee65 * ee122;
          ee127 = qalpha + ee57/ee32 - ee64;
          ee128 = R_pow(ee127, 2);
          ee130 = ee65 * ee61;
          ee132 = 2 * (ee61 * ee127);
            ee133 = ee130 * ee132;
          ee134 = R_pow(ee128, 2);
          ee140 = ee124/ee128 + ee133/ee134 - (ee121/ee127 - ee62/ee128);
            ee141 = ee65/ee127;
          ee142 = R_pow(ee141, 4);
          ee146 = ee46/ee127 + ee130/ee128;
          ee147 = R_pow(ee141, 3);
          ee149 = 4 * (ee146 * ee147);
            ee152 = 5 * (ee140 * ee142 + ee146 * ee149);
            ee156 = R_pow(ee141, 2);
          ee161 = 70 * ee142 - 315 * ee147 + 540 * ee156 - 420 * ee141 + 
            126;
            ee164 = 5 * (ee146 * ee142);
            ee167 = 3 * (ee146 * ee156);
            ee171 = 2 * (ee146 * ee141);
            ee175 = 70 * ee149 - 315 * ee167 + 540 * ee171 - 420 * ee146;
          ee176 = ee164 * ee175;
          ee178 = R_pow(ee141, 5);
          ee182 = 4 * (ee140 * ee147 + ee146 * ee167);
            ee187 = 3 * (ee140 * ee156 + ee146 * ee171);
            ee193 = 2 * (ee140 * ee141 + ee146 * ee146);
            ee197 = 70 * ee182 - 315 * ee187 + 540 * ee193 - 420 * ee140;
          ee200 = ee152 * ee161 + ee176 + (ee176 + ee178 * ee197);
            ee201 = y - qalpha;
          ee204 = ee201 * ee32/ee23 + ee18;
          ee206 = -1/ee7;
          ee207 = R_pow(ee204, ee206);
          ee211 = ee164 * ee161 + ee178 * ee175;
          ee212 = log(ee204);
            ee213 = R_pow(ee7, 2);
          ee214 = ee13/ee213;
          ee215 = ee212 * ee214;
          ee217 = ee206 - 1;
          ee218 = R_pow(ee204, ee217);
          ee221 = ee21 + ee201 * ee42/ee23;
          ee222 = ee206 * ee221;
          ee224 = ee207 * ee215 - ee218 * ee222;
          ee225 = ee211 * ee224;
          ee227 = ee178 * ee161;
          ee231 = 2 * (ee13 * ee7);
            ee232 = ee13 * ee231;
          ee233 = R_pow(ee213, 2);
          ee235 = ee71/ee213 + ee232/ee233;
          ee237 = ee221/ee204;
          ee239 = ee212 * ee235 + ee237 * ee214;
          ee243 = ee217 - 1;
          ee244 = R_pow(ee204, ee243);
          ee245 = ee217 * ee221;
          ee247 = ee218 * ee215 - ee244 * ee245;
          ee249 = ee214 * ee221;
          ee252 = ee201 * ee94/ee23 + ee79;
          ee254 = ee249 - ee206 * ee252;
          ee257 = ee224 * ee215 - ee207 * ee239 - (ee247 * ee222 + 
            ee218 * ee254);
            ee261 = 1 - ee227;
          ee262 = ee19 - ee10;
          ee264 = ee10 - ee50;
          ee268 = ee39 - ee36;
          ee269 = (y - (ee64 - ee127 * ee262/ee264)) * ee268;
          ee271 = ee127 * ee268/ee264;
          ee275 = exp(-(ee269/ee271 - ee19));
            ee279 = (ee46 - ee61 * ee262/ee264) * ee268;
          ee281 = ee61 * ee268/ee264;
          ee282 = ee279 * ee281;
          ee284 = ee122 * ee268/ee264;
          ee286 = ee282 - ee269 * ee284;
          ee287 = R_pow(ee271, 2);
          ee289 = ee269 * ee281;
          ee291 = 2 * (ee281 * ee271);
            ee292 = ee289 * ee291;
          ee293 = R_pow(ee287, 2);
          ee299 = (ee121 - ee122 * ee262/ee264) * ee268;
          ee303 = ee286/ee287 + ee292/ee293 - (ee299/ee271 - ee282/ee287);
            ee307 = ee279/ee271 + ee289/ee287;
          ee308 = ee275 * ee307;
          ee310 = ee275 * ee303 - ee308 * ee307;
          ee312 = ee211 * ee308;
          ee316 = ee261 * ee310 - ee312 + (ee200 * ee275 - ee312);
            ee319 = -630 * ee142;
          ee320 = 1 - ee141;
          ee321 = R_pow(ee320, 4);
          ee322 = ee319 * ee321;
          ee323 = ee322/ee127;
          ee325 = R_pow(ee320, 3);
          ee327 = 4 * (ee146 * ee325);
            ee329 = 630 * ee149;
          ee330 = ee329 * ee321;
          ee331 = ee319 * ee327 + ee330;
          ee333 = ee322 * ee61;
          ee335 = ee331/ee127 - ee333/ee128;
          ee336 = ee335 * ee224;
          ee339 = R_pow(ee320, 2);
          ee341 = 3 * (ee146 * ee339);
            ee344 = 4 * (ee140 * ee325 - ee146 * ee341);
            ee346 = ee329 * ee327;
          ee348 = 630 * ee182;
          ee350 = ee348 * ee321 - ee346;
          ee351 = ee319 * ee344 - ee346 + ee350;
          ee353 = ee331 * ee61;
          ee357 = ee322 * ee122 + ee353;
          ee359 = ee333 * ee132;
          ee362 = ee351/ee127 + ee353/ee128 + (ee357/ee128 - ee359/ee134);
            ee367 = 1 + 1/ee7;
          ee369 = R_pow(ee204, -ee367);
          ee371 = ee367 + 1;
          ee373 = R_pow(ee204, -ee371);
          ee374 = ee367 * ee221;
          ee377 = ee373 * ee374 + ee369 * ee215;
          ee378 = ee211 * ee377;
          ee380 = ee371 + 1;
          ee382 = R_pow(ee204, -ee380);
          ee383 = ee371 * ee221;
          ee386 = ee382 * ee383 + ee373 * ee215;
          ee389 = ee367 * ee252 + ee249;
          ee395 = ee386 * ee374 - ee373 * ee389 + (ee377 * ee215 - 
            ee369 * ee239);
            ee398 = ee200 * ee369 + ee378 + (ee378 + ee227 * ee395);
            ee399 = ee32/ee23;
          ee403 = ee211 * ee369 + ee227 * ee377;
          ee404 = ee42/ee23;
          ee405 = ee403 * ee404;
          ee407 = ee227 * ee369;
          ee408 = ee94/ee23;
          ee411 = ee398 * ee399 - ee405 - (ee405 - ee407 * ee408);
            ee415 = ee403 * ee399 - ee407 * ee404;
          ee416 = ee415 * ee13;
          ee419 = ee407 * ee399;
          ee421 = ee416 - ee419 * ee71;
          ee423 = ee419 * ee13;
          ee424 = ee423 * ee231;
          ee429 = 630 * ee142;
          ee432 = ee350 - (ee346 + ee429 * ee344);
            ee435 = ee330 - ee429 * ee327;
          ee436 = ee435 * ee61;
          ee439 = ee429 * ee321;
          ee441 = ee436 - ee439 * ee122;
          ee443 = ee439 * ee61;
          ee444 = ee443 * ee132;
          ee447 = ee432/ee127 + ee436/ee128 + (ee441/ee128 + ee444/ee134);
            ee451 = ee435/ee127 + ee443/ee128;
          ee452 = ee451 * ee308;
          ee454 = ee439/ee127;
          ee459 = ee261 * ee275;
          ee460 = ee268 * ee284;
          ee462 = ee268 * ee281;
          ee463 = ee462 * ee291;
          ee465 = ee460/ee287 - ee463/ee293;
          ee469 = ee261 * ee308 + ee211 * ee275;
          ee470 = ee462/ee287;
          ee471 = ee469 * ee470;
          ee473 = ee268/ee271;
          ee477 = ee323 * ee257 - ee336 - (ee362 * ee207 + ee336) + 
            (ee411/ee7 - ee416/ee213 - (ee421/ee213 - ee424/ee233)) + 
            (ee447 * ee275 - ee452 - (ee452 + ee454 * ee310)) - (ee459 * 
            ee465 + ee471 + (ee316 * ee473 + ee471));
            ee484 = ee323 * ee207 + ee419/ee7 + ee454 * ee275 + ee459 * 
              ee473;
              ee500 = ee323 * ee224 - ee335 * ee207 + (ee415/ee7 - ee423/ee213) + 
                (ee451 * ee275 - ee454 * ee308) + (ee459 * ee470 - ee469 * 
                ee473);
            ee501 = ee500 * ee500;
          ee502 = R_pow(ee484, 2);
          ee506 = ee124 * ee132;
          ee510 = 2 * (ee61 * ee61 + ee122 * ee127);
            ee515 = 2 * (ee132 * ee128);
            ee517 = R_pow(ee134, 2);
          ee520 = ee46 * ee122;
          ee521 = ee121 * ee61;
          ee522 = ee520 + ee521;
          ee523 = ee52 * ee72;
          ee538 = ee71 - ((ee11 * (2 * (ee4 * ee4 + ee66)) + ee68)/ee69 - 
            ee68 * (2 * (ee67 * ee12))/R_pow(ee69, 2));
            ee543 = ee21 * ee76;
          ee549 = ee543 + ee79 * ee20 + (ee18 * (ee19 * ee538) + ee543);
            ee553 = ee81 * ee42;
          ee556 = ee54 * ee94;
          ee557 = ee556 + ee553;
          ee563 = ee38 * ee86;
          ee570 = ee41 * ee90;
          ee577 = ee563 + ee89 * ee37 + (ee28 * (ee36 * ee538) + ee563) - 
            (ee570 + ee93 * ee40 + (ee31 * (ee39 * ee538) + ee570));
            ee582 = ee96 * ee99;
          ee588 = 2 * (ee42 * ee42 + ee94 * ee32);
            ee593 = 2 * (ee99 * ee44);
            ee595 = R_pow(ee101, 2);
          ee600 = ee15 * ee105;
          ee610 = ee110 * ee42;
          ee613 = ee24 * ee94;
          ee614 = ee613 + ee610;
          ee624 = ee116 * ee99;
          ee634 = (ee600 + ee108 * ee14 + (ee9 * (ee10 * ee538) + 
            ee600) - ee549) * ee23/ee32 - ee610/ee44 - (ee614/ee44 - 
            ee112 * ee99/ee101) - ((ee614 + (ee35 * ee577 + ee613))/ee44 - 
            ee624/ee101 - ((ee43 * ee588 + ee624)/ee101 - ee118 * 
            ee593/ee595));
            ee635 = (ee523 + ee75 * ee51 + (ee49 * (ee50 * ee538) + 
              ee523) - ee549) * ee23/ee32 - ee553/ee44 - (ee557/ee44 - 
              ee83 * ee99/ee101) - ((ee557 + (ee57 * ee577 + ee556))/ee44 - 
              ee582/ee101 - ((ee58 * ee588 + ee582)/ee101 - ee100 * 
              ee593/ee595)) - ee634;
          ee651 = (ee506 - ee130 * ee510)/ee134 + ee133 * ee515/ee517 - 
            ((ee522 + (ee520 - ee65 * ee635))/ee128 - ee506/ee134) + 
            (ee634/ee127 - ee521/ee128 - (ee522/ee128 - ee62 * ee132/ee134));
            ee653 = ee140 * ee149;
          ee660 = ee152 * ee175;
          ee662 = ee164 * ee197;
          ee663 = ee660 + ee662;
          ee666 = ee140 * ee167;
          ee671 = 4 * (ee651 * ee147 + ee666 + (ee666 + ee146 * ee187));
            ee674 = ee140 * ee171;
          ee683 = ee140 * ee146;
          ee696 = 5 * (ee651 * ee142 + ee653 + (ee653 + ee146 * ee182)) * 
            ee161 + ee660 + ee663 + (ee663 + (ee662 + ee178 * (70 * 
            ee671 - 315 * (3 * (ee651 * ee156 + ee674 + (ee674 + 
            ee146 * ee193))) + 540 * (2 * (ee651 * ee141 + ee683 + 
            (ee683 + ee146 * ee140))) - 420 * ee651)));
            ee698 = ee200 * ee224;
          ee700 = ee211 * ee257;
          ee701 = ee698 + ee700;
          ee704 = ee224 * ee239;
          ee709 = 2 * (ee13 * ee13 - ee71 * ee7);
            ee711 = ee71 * ee231;
          ee715 = 2 * (ee231 * ee213);
            ee717 = R_pow(ee233, 2);
          ee725 = ee237 * ee235;
          ee734 = ee212 * ((ee13 * ee709 - ee711)/ee233 - ee232 * 
            ee715/ee717 - (ee538/ee213 + ee711/ee233)) - ee725 - 
            (ee725 + (ee252/ee204 - ee221 * ee221/R_pow(ee204, 2)) * ee214);
            ee754 = ee247 * ee254;
          ee756 = ee214 * ee252;
          ee758 = ee756 + ee235 * ee221;
          ee761 = ee549 + ee201 * ee577/ee23;
          ee768 = ee257 * ee215 - ee704 - (ee704 + ee207 * ee734) - 
            ((ee247 * ee215 - ee218 * ee239 - ((ee244 * ee215 - R_pow(ee204, (ee243 - 
            1)) * (ee243 * ee221)) * ee245 + ee244 * (ee249 - 
            ee217 * ee252))) * ee222 + ee754 + (ee754 - ee218 * 
            (ee758 + (ee756 - ee206 * ee761))));
            ee773 = ee286 * ee291;
          ee777 = 2 * (ee281 * ee281 + ee284 * ee271);
            ee782 = 2 * (ee291 * ee287);
            ee784 = R_pow(ee293, 2);
          ee787 = ee279 * ee284;
          ee788 = ee299 * ee281;
          ee789 = ee787 + ee788;
          ee791 = ee635 * ee268/ee264;
          ee813 = ee308 * ee303;
          ee817 = ee275 * ((ee773 - ee289 * ee777)/ee293 + ee292 * 
            ee782/ee784 - ((ee789 + (ee787 - ee269 * ee791))/ee287 - 
            ee773/ee293) + ((ee634 - ee635 * ee262/ee264) * ee268/ee271 - 
            ee788/ee287 - (ee789/ee287 - ee282 * ee291/ee293))) - 
            ee813 - (ee310 * ee307 + ee813);
            ee819 = ee211 * ee310;
          ee821 = ee200 * ee308;
          ee822 = ee821 + ee819;
          ee827 = ee261 * ee817 - ee819 - ee822 + (ee696 * ee275 - 
            ee821 - ee822);
            ee830 = ee335 * ee257;
          ee832 = ee362 * ee224;
          ee833 = ee832 + ee830;
          ee836 = ee140 * ee341;
          ee847 = 4 * (ee651 * ee325 - ee836 - (ee836 + ee146 * (3 * 
            (ee140 * ee339 - ee146 * (2 * (ee146 * ee320))))));
            ee849 = ee329 * ee344;
          ee851 = ee348 * ee327;
          ee852 = ee851 + ee849;
          ee857 = 630 * ee671 * ee321 - ee851 - ee852;
          ee860 = ee351 * ee61;
          ee863 = ee331 * ee122;
          ee864 = ee860 - ee863;
          ee874 = ee357 * ee132;
          ee890 = ee200 * ee377;
          ee892 = ee211 * ee395;
          ee893 = ee890 + ee892;
          ee912 = ee386 * ee389;
          ee921 = ee377 * ee239;
          ee932 = ee398 * ee404;
          ee934 = ee403 * ee408;
          ee935 = ee932 - ee934;
          ee943 = ee411 * ee13;
          ee946 = ee415 * ee71;
          ee947 = ee943 - ee946;
          ee957 = ee421 * ee231;
          ee974 = ee432 * ee61;
          ee977 = ee435 * ee122;
          ee978 = ee974 - ee977;
          ee988 = ee441 * ee132;
          ee1000 = ee447 * ee308;
          ee1002 = ee451 * ee310;
          ee1003 = ee1000 + ee1002;
          ee1010 = ee316 * ee470;
          ee1011 = ee469 * ee465;
          ee1012 = ee1010 - ee1011;
          ee1015 = ee460 * ee291;
          ee1034 = ee477 * ee500;

          out(j, 9) = ee696 * ee207 + ee698 + ee701 + (ee701 + 
          (ee700 + ee227 * ee768)) - ee827 - ((ee323 * ee768 - 
          ee830 - ee833 - (((ee319 * ee847 - ee849 - ee852 + ee857)/ee127 + 
          ee860/ee128 + (ee864/ee128 + ee353 * ee132/ee134) + ((ee864 - 
          (ee322 * ee635 + ee863))/ee128 + ee874/ee134 + ((ee333 * 
          ee510 + ee874)/ee134 - ee359 * ee515/ee517))) * ee207 + 
          ee832 + ee833) + (((ee696 * ee369 + ee890 + ee893 + (ee893 + 
          (ee892 + ee227 * (((R_pow(ee204, -(ee380 + 1)) * (ee380 * ee221) + 
          ee382 * ee215) * ee383 - ee382 * (ee371 * ee252 + 
          ee249) + (ee386 * ee215 - ee373 * ee239)) * ee374 - 
          ee912 - (ee912 - ee373 * (ee758 + (ee367 * ee761 + 
          ee756))) + (ee395 * ee215 - ee921 - (ee921 + ee369 * 
          ee734)))))) * ee399 - ee932 - ee935 - (ee935 - (ee934 - 
          ee407 * (ee577/ee23))))/ee7 - ee943/ee213 - (ee947/ee213 - 
          ee416 * ee231/ee233) - ((ee947 - (ee946 - ee419 * ee538))/ee213 - 
          ee957/ee233 - ((ee957 + ee423 * ee709)/ee233 - ee424 * 
          ee715/ee717))) + (((ee857 - (ee852 + (ee849 + ee429 * 
          ee847)))/ee127 + ee974/ee128 + (ee978/ee128 + ee436 * 
          ee132/ee134) + ((ee978 - (ee977 - ee439 * ee635))/ee128 + 
          ee988/ee134 + ((ee988 - ee443 * ee510)/ee134 + ee444 * 
          ee515/ee517))) * ee275 - ee1000 - ee1003 - (ee1003 + 
          (ee1002 + ee454 * ee817))) - (ee1012 - (ee459 * (ee268 * 
          ee791/ee287 - ee1015/ee293 - ((ee462 * ee777 + ee1015)/ee293 - 
          ee463 * ee782/ee784)) + ee1011) + (ee827 * ee473 + ee1010 + 
          ee1012)))/ee484 - ee1034/ee502 - ((ee1034 + ee500 * ee477)/ee502 - 
          ee501 * (2 * (ee500 * ee484))/R_pow(ee502, 2)));
          }}}
    
return out;

}

