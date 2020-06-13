// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// Probit

double exipd0(arma::vec yvec, arma::uvec zvec, arma::vec pars, arma::vec nmax, arma::mat X, arma::uvec dupid, int dcate)
{   

arma::vec thetavec = X * pars;

if (dcate == 1) thetavec = thetavec.elem(dupid);    

int nobs = yvec.size();
    
double theta;

double nllh=0.0;

for (int j=0; j < nobs; j++) {

theta = thetavec[j];
theta = R::pnorm(theta, 0.0, 1.0, 1, 0);

nllh += nmax[0] * theta / yvec[j];
if (zvec[j] == 1) nllh -= log(theta);

}

return(nllh);

}

arma::mat exipd12(arma::vec yvec, arma::uvec zvec, arma::vec pars, arma::vec nmax, arma::mat X, arma::uvec dupid, int dcate)
{   

arma::vec thetavec = X * pars;

if (dcate == 1) thetavec = thetavec.elem(dupid);    

int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 2);
double theta, y, ee1, ee2;
    
for (int j=0; j < nobs; j++) {

theta = thetavec[j];
y = yvec[j];
ee1 = R::dnorm(theta, 0.0, 1.0, 0);
ee2 = R::pnorm(theta, 0.0, 1.0, 1, 0);

if (zvec[j] == 1) {

out(j, 0) = ee1 * (nmax[0]/y - 1/ee2);
out(j, 1) = -(ee1 * (nmax[0] * theta/y - (ee1/ee2 + theta)/ee2));
    
} else {
    
out(j, 0) = ee1 * (nmax[0]/y);
out(j, 1) = -(ee1 * (nmax[0] * theta/y));

}
}

return(out);

}

arma::mat exipd34(arma::vec yvec, arma::uvec zvec, arma::vec pars, arma::vec nmax, arma::mat X, arma::uvec dupid, int dcate)
{   

arma::vec thetavec = X * pars;

if (dcate == 1) thetavec = thetavec.elem(dupid);    

int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 2);
double theta, y, ee1, ee2, ee3, ee4, ee5, ee6, ee7;
    
for (int j=0; j < nobs; j++) {

theta = thetavec[j];
y = yvec[j];
ee1 = R::dnorm(theta, 0.0, 1.0, 0);
ee2 = R::pnorm(theta, 0.0, 1.0, 1, 0);
ee3 = theta * theta;
ee4 = ee2/ee1;
ee5 = 2 * ee4;
ee6 = 3 - ee3;
ee7 = theta * ee1;

if (zvec[j] == 1) {

out(j, 0) = -(ee2 * (nmax[0] * (1 - ee3)/y - 
            (1 - ((ee5 + 2 * theta) * ee2/ee1 + theta * (ee4 + theta)))/ee1));
out(j, 1) = ee2 * (nmax[0] * theta * ee6/y - 
            ((2 + 2 * (1 - 2 * ee3) + ((2 * (ee2 - ee7) - 
                (4 * ee7 + 8 * ee2))/ee1 - 4 * theta) * ee2/ee1 - 
                theta * (ee5 + 3 * theta)) * ee2/ee1 + theta * 
                ee6)/ee1);

} else {
    
out(j, 0) = -(ee2 * (nmax[0] * (1 - ee3)/y));
out(j, 1) = ee2 * (nmax[0] * theta * ee6/y);

}
}

return(out);

}

// Logistic

double exild0(arma::vec yvec, arma::uvec zvec, arma::vec pars, arma::vec nmax, arma::mat X, arma::uvec dupid, int dcate)
{   

arma::vec thetavec = X * pars;

if (dcate == 1) thetavec = thetavec.elem(dupid);    

int nobs = yvec.size();
    
double theta;

double nllh=0.0;

for (int j=0; j < nobs; j++) {

theta = 1 / (1 + exp(-thetavec[j]));

nllh += nmax[0] * theta / yvec[j];
if (zvec[j] == 1) nllh -= log(theta);

}

return(nllh);

}

arma::mat exild12(arma::vec yvec, arma::uvec zvec, arma::vec pars, arma::vec nmax, arma::mat X, arma::uvec dupid, int dcate)
{   

arma::vec thetavec = X * pars;

if (dcate == 1) thetavec = thetavec.elem(dupid);    

int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 2);
double theta, y, ee1, ee2, ee3;
    
for (int j=0; j < nobs; j++) {

theta = thetavec[j];
y = yvec[j];
ee1 = exp(-theta);
ee2 = 1 + ee1;
ee3 = y * ee2;

if (zvec[j] == 1) {

out(j, 0) = ee1 * (nmax[0]/ee3 - 1)/ee2;
out(j, 1) = -(ee1 * (nmax[0] * (1 - 2 * (ee1/ee2))/ee3 - 
    (1 + (1/(R_pow(1/ee2, 2) * ee2 * ee2) - 2) * ee1/ee2))/ee2);
    
} else {
    
out(j, 0) = ee1 * (nmax[0]/ee3)/ee2;
out(j, 1) = -(ee1 * (nmax[0] * (1 - 2 * (ee1/ee2))/ee3)/ee2);

}
}

return(out);

}

arma::mat exild34(arma::vec yvec, arma::uvec zvec, arma::vec pars, arma::vec nmax, arma::mat X, arma::uvec dupid, int dcate)
{   

arma::vec thetavec = X * pars;

if (dcate == 1) thetavec = thetavec.elem(dupid);    

int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 2);
double theta, y, ee1, ee2, ee3, ee4, ee6, ee7, ee8, ee9;
double ee10, ee11, ee12, ee15, ee16, ee17;
double ee20, ee22, ee23, ee29, ee30, ee31;
    
for (int j=0; j < nobs; j++) {

theta = thetavec[j];
y = yvec[j];
ee1 = exp(-theta);
ee2 = 1 + ee1;
ee3 = 1 + 2 * ee1;
ee4 = 2 * ee3;
ee6 = R_pow(1/ee2, 2);
ee7 = ee4 + 2 * ee2;
ee8 = ee1/ee2;
ee9 = 1 - 2 * ee8;
ee10 = 4 * ee1;
ee11 = 8 * ee1;
ee12 = ee6 * ee2 * ee2;
ee15 = 2 + (ee7 - ee11)/ee2;
ee16 = ee6 * ee2;
ee17 = 1 + (2/ee12 - 2) * ee1/ee2;
ee20 = 2 * (1 + ee10);
ee22 = ee15 * ee1/ee2;
ee23 = 4 * ee2;
ee29 = ((4 * (ee4 + ee10) + 8 * ee7 - 64 * ee1)/ee2 + 8) * ee1;
ee30 = 6 * ee3;
ee31 = y * ee2;

if (zvec[j] == 1) {

out(j, 0) = ee1 * (nmax[0] * (1 - ee22)/ee31 - 
    (1 + (((ee17 + 2 * ee9)/ee16 + ee11 - ee7)/ee2 - 2) * 
    ee1/ee2))/ee2;
out(j, 1) = -(ee1 * 
    (nmax[0] * (1 - (2 + (ee20 + ee23 + ee30 - ee29)/ee2) * 
    ee1/ee2)/ee31 + ((2 + (ee20 + ((2 * ee15 + (2 * 
    (3 * ee8 - 1) - (8 * ee9 + 8 * (ee1/(ee6 * R_pow(ee2, 3)))))/ee12) * 
    ee1/ee2 - (ee17 * ee9 + 2 + 2 * (1 + R_pow(ee9, 2) - ee22)))/ee16 + 
    ee23 + ee30 - ee29)/ee2) * ee1/ee2 - 1))/ee2);

} else {
    
out(j, 0) = ee1 * (nmax[0] * (1 - ee22)/ee31)/ee2;
out(j, 1) = -(ee1 * 
    (nmax[0] * (1 - (2 + (ee20 + ee23 + ee30 - ee29)/ee2) * 
    ee1/ee2)/ee31)/ee2);


}
}

return(out);

}

// cloglog

double exicd0(arma::vec yvec, arma::uvec zvec, arma::vec pars, arma::vec nmax, arma::mat X, arma::uvec dupid, int dcate)
{   

arma::vec thetavec = X * pars;

if (dcate == 1) thetavec = thetavec.elem(dupid);    

int nobs = yvec.size();
    
double theta;

double nllh=0.0;

for (int j=0; j < nobs; j++) {

theta = 1 - exp(-exp(thetavec[j]));

nllh += nmax[0] * theta / yvec[j];
if (zvec[j] == 1) nllh -= log(theta);

}

return(nllh);

}

arma::mat exicd12(arma::vec yvec, arma::uvec zvec, arma::vec pars, arma::vec nmax, arma::mat X, arma::uvec dupid, int dcate)
{   

arma::vec thetavec = X * pars;

if (dcate == 1) thetavec = thetavec.elem(dupid);    

int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 2);
double theta, y, ee1, ee3, ee4, ee5;
    
for (int j=0; j < nobs; j++) {

theta = thetavec[j];
y = yvec[j];

ee1 = exp(theta);
ee3 = exp(-ee1);
ee4 = 1 - ee3;
ee5 = ee3 * ee1;

if (zvec[j] == 1) {

out(j, 0) = ee5 * (nmax[0]/y - 1/ee4);
out(j, 1) = ee5 * (nmax[0] * (1 - ee1)/y - 1 * (1 - (1 + ee3/ee4) * ee1)/ee4);

} else {
    
out(j, 0) = ee5 * nmax[0]/y;
out(j, 1) = ee5 * (nmax[0] * (1 - ee1)/y);

}
}

return(out);

}

arma::mat exicd34(arma::vec yvec, arma::uvec zvec, arma::vec pars, arma::vec nmax, arma::mat X, arma::uvec dupid, int dcate)
{   

arma::vec thetavec = X * pars;

if (dcate == 1) thetavec = thetavec.elem(dupid);    

int nobs = yvec.size();
arma::mat out = arma::mat(nobs, 2);
double theta, y, ee1, ee3, ee4, ee5, ee6, ee7, ee8, ee9;
double ee10, ee11;
    
for (int j=0; j < nobs; j++) {

theta = thetavec[j];
y = yvec[j];

ee1 = exp(theta);
ee3 = exp(-ee1);
ee4 = 1 - ee3;
ee5 = 1 - ee1;
ee6 = ee3 * ee1;
ee7 = 3 - ee1;
ee8 = ee4 * ee5;
ee9 = (1 + 2 * (ee3/ee4)) * ee1;
ee10 = ee7 * ee1;
ee11 = (6 - ee1) * ee1;

if (zvec[j] == 1) {

out(j, 0) = ee6 * (nmax[0] * (1 - ee10)/y - (1 - ((1 + 2 * ee5 - 
              ee9) * ee3/ee4 + 3 - ee1) * ee1)/ee4);

out(j, 1) = ee6 * (nmax[0] * (1 - (7 - ee11) * 
              ee1)/y - (1 - (((1 - ee9) * ee5 + 2 + 2 * (R_pow(ee5, 2) + 
              1 - ee10) - (((2 * (ee8 + ee6) + 4 * ee8 - 8 * ee6)/ee4 + 
              4 * ee5) * ee3/ee4 + 2 * ee7) * ee1) * ee3/ee4 + 
              7 - ee11) * ee1)/ee4);

} else {
    
out(j, 0) = ee6 * (nmax[0] * (1 - ee10)/y);

out(j, 1) = ee6 * (nmax[0] * (1 - (7 - ee11) * ee1)/y);

}
}

return(out);

}

// //' Extremal index (i.e. modified Frechet distribution) negative log-likelihood
// //'
// //' @param yvec a vector
// //' @param zvec a vector
// //' @param pars a list of vectors of coefficients for extremal index
// //' @param nmax an integer, the number of running max used
// //' @param X a design matrix for the extremal index, pre-link
// //' @param dupid a scalar or vector, identifying duplicates in Xs; -1 corresponds to no duplicates
// //' @param link an integer identifying the link function: 1=probit, 2=logit, 3=cloglog
// //' @return exid0 a scalar, the negative log-liklihood
// //' @return exid12 a matrix, first then second derivatives w.r.t. extremal index
// //' @return exid34 a matrix, third then fourth derivatives w.r.t. extremal index
// //' @examples
// //' ## to follow
// //' @export
// [[Rcpp::export]]
double exid0(arma::vec yvec, arma::uvec zvec, arma::vec pars, arma::vec nmax, arma::mat X, arma::uvec dupid, int dcate, int link)
{   

double out;
    
if (link == 0) { // probit
  out = exipd0(yvec, zvec, pars, nmax, X, dupid, dcate);
} else {
  if (link == 1) { // logistic
    out = exild0(yvec, zvec, pars, nmax, X, dupid, dcate);
  } else { // cloglog
    out = exicd0(yvec, zvec, pars, nmax, X, dupid, dcate);
  }
}

return out;
}

// //' @rdname exid0
// [[Rcpp::export]]
arma::mat exid12(arma::vec yvec, arma::uvec zvec, arma::vec pars, arma::vec nmax, arma::mat X, arma::uvec dupid, int dcate, int link)
{   

arma::mat out;
    
if (link == 0) { // probit
  out = exipd12(yvec, zvec, pars, nmax, X, dupid, dcate);
} else {
  if (link == 1) { // logistic
    out = exild12(yvec, zvec, pars, nmax, X, dupid, dcate);
  } else { // cloglog
    out = exicd12(yvec, zvec, pars, nmax, X, dupid, dcate);
  }
}

return out;

}

// //' @rdname exid0
// [[Rcpp::export]]
arma::mat exid34(arma::vec yvec, arma::uvec zvec, arma::vec pars, arma::vec nmax, arma::mat X, arma::uvec dupid, int dcate, int link)
{   

arma::mat out;
    
if (link == 0) { // probit
  out = exipd34(yvec, zvec, pars, nmax, X, dupid, dcate);
} else {
  if (link == 1) { // logistic
    out = exild34(yvec, zvec, pars, nmax, X, dupid, dcate);
  } else { // cloglog
    out = exicd34(yvec, zvec, pars, nmax, X, dupid, dcate);
  }
}

return out;

}

//' Running maximum
//' 
//' Running \eqn{n}-value maximum and data frame with variable swapped for running maximum
//'
//' @param y a vector
//' @param n an integer giving the number of observations to calculate running maxmimum over; defaults to 2
//'
//' @return \code{runmax} returns a vector of the same dimension as \code{y}
//'
//' @examples
//' runmax(runif(10), 5)
//'
//' @name runmax
//'
//' @export
// [[Rcpp::export]]
arma::vec runmax(arma::vec y, int n)
{   

int nobs = y.size();
arma::vec ymax = arma::vec(nobs - n + 1);
    
for (int j=0; j < nobs - n + 1; j++) {

ymax[j] = y[j];

for (int k=1; k < n; k++) {
    
    if (y[j + k] > ymax[j]) ymax[j] = y[j + k];

}

}
        
return(ymax);

}
