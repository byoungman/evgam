// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::export(.gH1)]]
Rcpp::List gH1(arma::mat gh, arma::mat X1, const arma::uvec dupid, int dcate, int sand, int deriv)
{

Rcpp::List out(2);
arma::mat g;

if (dcate == 1) {
    X1 = X1.rows(dupid);
}

if (deriv > 1) {

arma::mat H;
H = X1.t() * (X1.each_col() % gh.col(1));
out(1) = H;
}

X1.each_col() %= gh.col(0);

if (sand == 0) {
    g = sum(X1, 0);
} else {
    g = X1;
}

out(0) = g;

return(out);

}

// [[Rcpp::export(.gH2)]]
Rcpp::List gH2(arma::mat gh, arma::mat X1, arma::mat X2, const arma::uvec dupid, int dcate, int sand, int deriv)
{

Rcpp::List out(2);

arma::mat g;

if (dcate == 1) {
    X1 = X1.rows(dupid);
    X2 = X2.rows(dupid);
}

if (deriv > 1) {
int s1 = 0;
int e1 = X1.n_cols - 1;
int s2 = e1 + 1;
int e2 = e1 + X2.n_cols;
arma::mat H =  arma::mat(e2 + 1, e2 + 1);
H.submat(s1, s1, e1, e1) = X1.t() * (X1.each_col() % gh.col(2));
H.submat(s2, s1, e2, e1) = X2.t() * (X1.each_col() % gh.col(3));
H.submat(s1, s2, e1, e2) = H.submat(s2, s1, e2, e1).t();
H.submat(s2, s2, e2, e2) = X2.t() * (X2.each_col() % gh.col(4));
out(1) = H;
}

X1.each_col() %= gh.col(0);
X2.each_col() %= gh.col(1);

if (sand == 0) {
    g = join_rows(sum(X1, 0), sum(X2, 0));
} else {
    g = join_rows(X1, X2);
}

out(0) = g;

return(out);

}

// [[Rcpp::export(.gH3)]]
Rcpp::List gH3(arma::mat gh, arma::mat X1, arma::mat X2, arma::mat X3, const arma::uvec dupid, int dcate, int sand, int deriv)
{

Rcpp::List out(2);
arma::mat g;

if (dcate == 1) {
    X1 = X1.rows(dupid);
    X2 = X2.rows(dupid);
    X3 = X3.rows(dupid);
}

if (deriv > 1) {
int s1 = 0;
int e1 = X1.n_cols - 1;
int s2 = e1 + 1;
int e2 = e1 + X2.n_cols;
int s3 = e2 + 1;
int e3 = e2 + X3.n_cols;
arma::mat H =  arma::mat(e3 + 1, e3 + 1);
H.submat(s1, s1, e1, e1) = X1.t() * (X1.each_col() % gh.col(3));
H.submat(s2, s1, e2, e1) = X2.t() * (X1.each_col() % gh.col(4));
H.submat(s3, s1, e3, e1) = X3.t() * (X1.each_col() % gh.col(5));
H.submat(s1, s2, e1, e2) = H.submat(s2, s1, e2, e1).t();
H.submat(s1, s3, e1, e3) = H.submat(s3, s1, e3, e1).t();
H.submat(s2, s2, e2, e2) = X2.t() * (X2.each_col() % gh.col(6));
H.submat(s3, s2, e3, e2) = X3.t() * (X2.each_col() % gh.col(7));
H.submat(s2, s3, e2, e3) = H.submat(s3, s2, e3, e2).t();
H.submat(s3, s3, e3, e3) = X3.t() * (X3.each_col() % gh.col(8));
out(1) = H;
}

X1.each_col() %= gh.col(0);
X2.each_col() %= gh.col(1);
X3.each_col() %= gh.col(2);

if (sand == 0) {
    g = join_rows(join_rows(sum(X1, 0), sum(X2, 0)), sum(X3, 0));
} else {
    g = join_rows(join_rows(X1, X2), X3);
}

out(0) = g;

return(out);

}

// added with evgam_0.1.2 (05/04/2020)

// [[Rcpp::export(.gH4)]]
Rcpp::List gH4(arma::mat gh, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, const arma::uvec dupid, int dcate, int sand, int deriv)
{

Rcpp::List out(2);
arma::mat g;

if (dcate == 1) {
    X1 = X1.rows(dupid);
    X2 = X2.rows(dupid);
    X3 = X3.rows(dupid);
    X4 = X4.rows(dupid);
}

if (deriv > 1) {

int s1 = 0;
int e1 = X1.n_cols - 1;
int s2 = e1 + 1;
int e2 = e1 + X2.n_cols;
int s3 = e2 + 1;
int e3 = e2 + X3.n_cols;
int s4 = e3 + 1;
int e4 = e3 + X4.n_cols;

arma::mat H =  arma::mat(e4 + 1, e4 + 1);

H.submat(s1, s1, e1, e1) = X1.t() * (X1.each_col() % gh.col(4));
H.submat(s2, s1, e2, e1) = X2.t() * (X1.each_col() % gh.col(5));
H.submat(s3, s1, e3, e1) = X3.t() * (X1.each_col() % gh.col(6));
H.submat(s4, s1, e4, e1) = X4.t() * (X1.each_col() % gh.col(7));
H.submat(s1, s2, e1, e2) = H.submat(s2, s1, e2, e1).t();
H.submat(s1, s3, e1, e3) = H.submat(s3, s1, e3, e1).t();
H.submat(s1, s4, e1, e4) = H.submat(s4, s1, e4, e1).t();

H.submat(s2, s2, e2, e2) = X2.t() * (X2.each_col() % gh.col(8));
H.submat(s3, s2, e3, e2) = X3.t() * (X2.each_col() % gh.col(9));
H.submat(s4, s2, e4, e2) = X4.t() * (X2.each_col() % gh.col(10));
H.submat(s2, s3, e2, e3) = H.submat(s3, s2, e3, e2).t();
H.submat(s2, s4, e2, e4) = H.submat(s4, s2, e4, e2).t();

H.submat(s3, s3, e3, e3) = X3.t() * (X3.each_col() % gh.col(11));
H.submat(s4, s3, e4, e3) = X4.t() * (X3.each_col() % gh.col(12));
H.submat(s3, s4, e3, e4) = H.submat(s4, s3, e4, e3).t();

H.submat(s4, s4, e4, e4) = X4.t() * (X4.each_col() % gh.col(13));

out(1) = H;
    
}

X1.each_col() %= gh.col(0);
X2.each_col() %= gh.col(1);
X3.each_col() %= gh.col(2);
X4.each_col() %= gh.col(3);

if (sand == 0) {
    g = join_rows(join_rows(join_rows(sum(X1, 0), sum(X2, 0)), sum(X3, 0)), sum(X4, 0));
} else {
    g = join_rows(join_rows(join_rows(X1, X2), X3), X4);
}

out(0) = g;

return(out);

}

// [[Rcpp::export(.gH5)]]
Rcpp::List gH5(arma::mat gh, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, arma::mat X5, const arma::uvec dupid, int dcate, int sand, int deriv)
{

Rcpp::List out(2);
arma::mat g;

if (dcate == 1) {
  X1 = X1.rows(dupid);
  X2 = X2.rows(dupid);
  X3 = X3.rows(dupid);
  X4 = X4.rows(dupid);
  X5 = X5.rows(dupid);
}

if (deriv > 1) {

int s1 = 0;
int e1 = X1.n_cols - 1;
int s2 = e1 + 1;
int e2 = e1 + X2.n_cols;
int s3 = e2 + 1;
int e3 = e2 + X3.n_cols;
int s4 = e3 + 1;
int e4 = e3 + X4.n_cols;
int s5 = e4 + 1;
int e5 = e4 + X5.n_cols;

arma::mat H =  arma::mat(e5 + 1, e5 + 1);

H.submat(s1, s1, e1, e1) = X1.t() * (X1.each_col() % gh.col(5));
H.submat(s2, s1, e2, e1) = X2.t() * (X1.each_col() % gh.col(6));
H.submat(s3, s1, e3, e1) = X3.t() * (X1.each_col() % gh.col(7));
H.submat(s4, s1, e4, e1) = X4.t() * (X1.each_col() % gh.col(8));
H.submat(s5, s1, e5, e1) = X5.t() * (X1.each_col() % gh.col(9));
H.submat(s1, s2, e1, e2) = H.submat(s2, s1, e2, e1).t();
H.submat(s1, s3, e1, e3) = H.submat(s3, s1, e3, e1).t();
H.submat(s1, s4, e1, e4) = H.submat(s4, s1, e4, e1).t();
H.submat(s1, s5, e1, e5) = H.submat(s5, s1, e5, e1).t();

H.submat(s2, s2, e2, e2) = X2.t() * (X2.each_col() % gh.col(10));
H.submat(s3, s2, e3, e2) = X3.t() * (X2.each_col() % gh.col(11));
H.submat(s4, s2, e4, e2) = X4.t() * (X2.each_col() % gh.col(12));
H.submat(s5, s2, e5, e2) = X5.t() * (X2.each_col() % gh.col(13));
H.submat(s2, s3, e2, e3) = H.submat(s3, s2, e3, e2).t();
H.submat(s2, s4, e2, e4) = H.submat(s4, s2, e4, e2).t();
H.submat(s2, s5, e2, e5) = H.submat(s5, s2, e5, e2).t();

H.submat(s3, s3, e3, e3) = X3.t() * (X3.each_col() % gh.col(14));
H.submat(s4, s3, e4, e3) = X4.t() * (X3.each_col() % gh.col(15));
H.submat(s5, s3, e5, e3) = X5.t() * (X3.each_col() % gh.col(16));
H.submat(s3, s4, e3, e4) = H.submat(s4, s3, e4, e3).t();
H.submat(s3, s5, e3, e5) = H.submat(s5, s3, e5, e3).t();

H.submat(s4, s4, e4, e4) = X4.t() * (X4.each_col() % gh.col(17));
H.submat(s5, s4, e5, e4) = X5.t() * (X4.each_col() % gh.col(18));
H.submat(s4, s5, e4, e5) = H.submat(s5, s4, e5, e4).t();

H.submat(s5, s5, e5, e5) = X5.t() * (X5.each_col() % gh.col(19));

out(1) = H;
    
}

X1.each_col() %= gh.col(0);
X2.each_col() %= gh.col(1);
X3.each_col() %= gh.col(2);
X4.each_col() %= gh.col(3);
X5.each_col() %= gh.col(4);

if (sand == 0) {
    g = join_rows(join_rows(join_rows(join_rows(sum(X1, 0), sum(X2, 0)), sum(X3, 0)), sum(X4, 0)), sum(X5, 0));
} else {
    g = join_rows(join_rows(join_rows(join_rows(X1, X2), X3), X4), X5);
}

out(0) = g;

return(out);

}

// [[Rcpp::export(.gH6)]]
Rcpp::List gH6(arma::mat gh, arma::mat X1, arma::mat X2, arma::mat X3, arma::mat X4, arma::mat X5, arma::mat X6, const arma::uvec dupid, int dcate, int sand, int deriv)
{

Rcpp::List out(2);
arma::mat g;

if (dcate == 1) {
  X1 = X1.rows(dupid);
  X2 = X2.rows(dupid);
  X3 = X3.rows(dupid);
  X4 = X4.rows(dupid);
  X5 = X5.rows(dupid);
  X6 = X6.rows(dupid);
}

if (deriv > 1) {

int s1 = 0;
int e1 = X1.n_cols - 1;
int s2 = e1 + 1;
int e2 = e1 + X2.n_cols;
int s3 = e2 + 1;
int e3 = e2 + X3.n_cols;
int s4 = e3 + 1;
int e4 = e3 + X4.n_cols;
int s5 = e4 + 1;
int e5 = e4 + X5.n_cols;
int s6 = e5 + 1;
int e6 = e5 + X6.n_cols;

arma::mat H =  arma::mat(e6 + 1, e6 + 1);

H.submat(s1, s1, e1, e1) = X1.t() * (X1.each_col() % gh.col(6));
H.submat(s2, s1, e2, e1) = X2.t() * (X1.each_col() % gh.col(7));
H.submat(s3, s1, e3, e1) = X3.t() * (X1.each_col() % gh.col(8));
H.submat(s4, s1, e4, e1) = X4.t() * (X1.each_col() % gh.col(9));
H.submat(s5, s1, e5, e1) = X5.t() * (X1.each_col() % gh.col(10));
H.submat(s6, s1, e6, e1) = X6.t() * (X1.each_col() % gh.col(11));
H.submat(s1, s2, e1, e2) = H.submat(s2, s1, e2, e1).t();
H.submat(s1, s3, e1, e3) = H.submat(s3, s1, e3, e1).t();
H.submat(s1, s4, e1, e4) = H.submat(s4, s1, e4, e1).t();
H.submat(s1, s5, e1, e5) = H.submat(s5, s1, e5, e1).t();
H.submat(s1, s6, e1, e6) = H.submat(s6, s1, e6, e1).t();

H.submat(s2, s2, e2, e2) = X2.t() * (X2.each_col() % gh.col(12));
H.submat(s3, s2, e3, e2) = X3.t() * (X2.each_col() % gh.col(13));
H.submat(s4, s2, e4, e2) = X4.t() * (X2.each_col() % gh.col(14));
H.submat(s5, s2, e5, e2) = X5.t() * (X2.each_col() % gh.col(15));
H.submat(s6, s2, e6, e2) = X6.t() * (X2.each_col() % gh.col(16));
H.submat(s2, s3, e2, e3) = H.submat(s3, s2, e3, e2).t();
H.submat(s2, s4, e2, e4) = H.submat(s4, s2, e4, e2).t();
H.submat(s2, s5, e2, e5) = H.submat(s5, s2, e5, e2).t();
H.submat(s2, s6, e2, e6) = H.submat(s6, s2, e6, e2).t();

H.submat(s3, s3, e3, e3) = X3.t() * (X3.each_col() % gh.col(17));
H.submat(s4, s3, e4, e3) = X4.t() * (X3.each_col() % gh.col(18));
H.submat(s5, s3, e5, e3) = X5.t() * (X3.each_col() % gh.col(19));
H.submat(s6, s3, e6, e3) = X6.t() * (X3.each_col() % gh.col(20));
H.submat(s3, s4, e3, e4) = H.submat(s4, s3, e4, e3).t();
H.submat(s3, s5, e3, e5) = H.submat(s5, s3, e5, e3).t();
H.submat(s3, s6, e3, e6) = H.submat(s6, s3, e6, e3).t();

H.submat(s4, s4, e4, e4) = X4.t() * (X4.each_col() % gh.col(21));
H.submat(s5, s4, e5, e4) = X5.t() * (X4.each_col() % gh.col(22));
H.submat(s6, s4, e6, e4) = X6.t() * (X4.each_col() % gh.col(23));
H.submat(s4, s5, e4, e5) = H.submat(s5, s4, e5, e4).t();
H.submat(s4, s6, e4, e6) = H.submat(s6, s4, e6, e4).t();

H.submat(s5, s5, e5, e5) = X5.t() * (X5.each_col() % gh.col(24));
H.submat(s6, s5, e6, e5) = X6.t() * (X5.each_col() % gh.col(25));
H.submat(s5, s6, e5, e6) = H.submat(s6, s5, e6, e5).t();

H.submat(s6, s6, e6, e6) = X6.t() * (X6.each_col() % gh.col(26));

out(1) = H;
    
}

X1.each_col() %= gh.col(0);
X2.each_col() %= gh.col(1);
X3.each_col() %= gh.col(2);
X4.each_col() %= gh.col(3);
X5.each_col() %= gh.col(4);
X6.each_col() %= gh.col(5);

if (sand == 0) {
    g = join_rows(join_rows(join_rows(join_rows(join_rows(sum(X1, 0), sum(X2, 0)), sum(X3, 0)), sum(X4, 0)), sum(X5, 0)), sum(X6, 0));
} else {
    g = join_rows(join_rows(join_rows(join_rows(join_rows(X1, X2), X3), X4), X5), X6);
}

out(0) = g;

return(out);

}

arma::sp_mat schur_spmat(arma::sp_mat x, arma::vec y) {
  int nc = x.n_cols;
  for (int i = 0; i < nc; i++) {
    x.col(i) %= y;
  }
  return x;
}

// [[Rcpp::export(.gH1sp)]]
Rcpp::List gHsp1(arma::mat gh, arma::sp_mat X1, const arma::uvec dupid, int dcate, int sand, int deriv)
{

Rcpp::List out(2);
arma::mat g;

// if (dcate == 1) {
//     X1 = X1.rows(dupid);
// }

if (deriv > 1) {

arma::mat H;
// H = X1.t() * (X1.each_col() % gh.col(1));
H = X1.t() * schur_spmat(X1, gh.col(1));
out(1) = H;
}

// X1.each_col() %= gh.col(0);
X1 = schur_spmat(X1, gh.col(0));

if (sand == 0) {
    g = sum(X1, 0);
} else {
    g = X1;
}

out(0) = g;

return(out);

}

// [[Rcpp::export(.gHsp2)]]
Rcpp::List gHsp2(arma::mat gh, arma::sp_mat X1, arma::sp_mat X2, const arma::uvec dupid, int dcate, int sand, int deriv)
{

Rcpp::List out(2);

arma::mat g;

// if (dcate == 1) {
//     X1 = X1.rows(dupid);
//     X2 = X2.rows(dupid);
// }

if (deriv > 1) {
int s1 = 0;
int e1 = X1.n_cols - 1;
int s2 = e1 + 1;
int e2 = e1 + X2.n_cols;
arma::mat H =  arma::mat(e2 + 1, e2 + 1);
// H.submat(s1, s1, e1, e1) = X1.t() * (X1.each_col() % gh.col(2));
H.submat(s1, s1, e1, e1) = X1.t() * schur_spmat(X1, gh.col(2));
// H.submat(s2, s1, e2, e1) = X2.t() * (X1.each_col() % gh.col(3));
H.submat(s2, s1, e2, e1) = X2.t() * schur_spmat(X1, gh.col(3));
H.submat(s1, s2, e1, e2) = H.submat(s2, s1, e2, e1).t();
// H.submat(s2, s2, e2, e2) = X2.t() * (X2.each_col() % gh.col(4));
H.submat(s2, s2, e2, e2) = X2.t() * schur_spmat(X2, gh.col(4));
out(1) = H;
}

// X1.each_col() %= gh.col(0);
// X2.each_col() %= gh.col(1);
X1 = schur_spmat(X1, gh.col(0));
X2 = schur_spmat(X2, gh.col(1));

if (sand == 0) {
    g = join_rows(sum(X1, 0), sum(X2, 0));
} else {
    g = join_rows(X1, X2);
}

out(0) = g;

return(out);

}

// [[Rcpp::export(.gHsp3)]]
Rcpp::List gHsp3(arma::mat gh, arma::sp_mat X1, arma::sp_mat X2, arma::sp_mat X3, const arma::uvec dupid, int dcate, int sand, int deriv)
{

Rcpp::List out(2);
arma::mat g;

// if (dcate == 1) {
//     X1 = X1.rows(dupid);
//     X2 = X2.rows(dupid);
//     X3 = X3.rows(dupid);
// }

if (deriv > 1) {
int s1 = 0;
int e1 = X1.n_cols - 1;
int s2 = e1 + 1;
int e2 = e1 + X2.n_cols;
int s3 = e2 + 1;
int e3 = e2 + X3.n_cols;
arma::mat H =  arma::mat(e3 + 1, e3 + 1);
// H.submat(s1, s1, e1, e1) = X1.t() * (X1.each_col() % gh.col(3));
H.submat(s1, s1, e1, e1) = X1.t() * schur_spmat(X1, gh.col(3));
// H.submat(s2, s1, e2, e1) = X2.t() * (X1.each_col() % gh.col(4));
H.submat(s2, s1, e2, e1) = X2.t() * schur_spmat(X1, gh.col(4));
// H.submat(s3, s1, e3, e1) = X3.t() * (X1.each_col() % gh.col(5));
H.submat(s3, s1, e3, e1) = X3.t() * schur_spmat(X1, gh.col(5));
H.submat(s1, s2, e1, e2) = H.submat(s2, s1, e2, e1).t();
H.submat(s1, s3, e1, e3) = H.submat(s3, s1, e3, e1).t();
// H.submat(s2, s2, e2, e2) = X2.t() * (X2.each_col() % gh.col(6));
H.submat(s2, s2, e2, e2) = X2.t() * schur_spmat(X2, gh.col(6));
// H.submat(s3, s2, e3, e2) = X3.t() * (X2.each_col() % gh.col(7));
H.submat(s3, s2, e3, e2) = X3.t() * schur_spmat(X2, gh.col(7));
H.submat(s2, s3, e2, e3) = H.submat(s3, s2, e3, e2).t();
// H.submat(s3, s3, e3, e3) = X3.t() * (X3.each_col() % gh.col(8));
H.submat(s3, s3, e3, e3) = X3.t() * schur_spmat(X3, gh.col(8));
out(1) = H;
}

// X1.each_col() %= gh.col(0);
// X2.each_col() %= gh.col(1);
// X3.each_col() %= gh.col(2);
X1 = schur_spmat(X1, gh.col(0));
X2 = schur_spmat(X2, gh.col(1));
X3 = schur_spmat(X3, gh.col(2));

if (sand == 0) {
    g = join_rows(join_rows(sum(X1, 0), sum(X2, 0)), sum(X3, 0));
} else {
    g = join_rows(join_rows(X1, X2), X3);
}

out(0) = g;

return(out);

}

// added with evgam_0.1.2 (05/04/2020)

// [[Rcpp::export(.gHsp4)]]
Rcpp::List gHsp4(arma::mat gh, arma::sp_mat X1, arma::sp_mat X2, arma::sp_mat X3, arma::sp_mat X4, const arma::uvec dupid, int dcate, int sand, int deriv)
{

Rcpp::List out(2);
arma::mat g;

// if (dcate == 1) {
//     X1 = X1.rows(dupid);
//     X2 = X2.rows(dupid);
//     X3 = X3.rows(dupid);
//     X4 = X4.rows(dupid);
// }

if (deriv > 1) {

int s1 = 0;
int e1 = X1.n_cols - 1;
int s2 = e1 + 1;
int e2 = e1 + X2.n_cols;
int s3 = e2 + 1;
int e3 = e2 + X3.n_cols;
int s4 = e3 + 1;
int e4 = e3 + X4.n_cols;

arma::mat H =  arma::mat(e4 + 1, e4 + 1);

// H.submat(s1, s1, e1, e1) = X1.t() * (X1.each_col() % gh.col(4));
H.submat(s1, s1, e1, e1) = X1.t() * schur_spmat(X1, gh.col(4));
// H.submat(s2, s1, e2, e1) = X2.t() * (X1.each_col() % gh.col(5));
H.submat(s2, s1, e2, e1) = X2.t() * schur_spmat(X1, gh.col(5));
// H.submat(s3, s1, e3, e1) = X3.t() * (X1.each_col() % gh.col(6));
H.submat(s3, s1, e3, e1) = X3.t() * schur_spmat(X1, gh.col(6));
// H.submat(s4, s1, e4, e1) = X4.t() * (X1.each_col() % gh.col(7));
H.submat(s4, s1, e4, e1) = X4.t() * schur_spmat(X1, gh.col(7));
H.submat(s1, s2, e1, e2) = H.submat(s2, s1, e2, e1).t();
H.submat(s1, s3, e1, e3) = H.submat(s3, s1, e3, e1).t();
H.submat(s1, s4, e1, e4) = H.submat(s4, s1, e4, e1).t();

// H.submat(s2, s2, e2, e2) = X2.t() * (X2.each_col() % gh.col(8));
H.submat(s2, s2, e2, e2) = X2.t() * schur_spmat(X2, gh.col(8));
// H.submat(s3, s2, e3, e2) = X3.t() * (X2.each_col() % gh.col(9));
H.submat(s3, s2, e3, e2) = X3.t() * schur_spmat(X2, gh.col(9));
// H.submat(s4, s2, e4, e2) = X4.t() * (X2.each_col() % gh.col(10));
H.submat(s4, s2, e4, e2) = X4.t() * schur_spmat(X2, gh.col(10));
H.submat(s2, s3, e2, e3) = H.submat(s3, s2, e3, e2).t();
H.submat(s2, s4, e2, e4) = H.submat(s4, s2, e4, e2).t();

// H.submat(s3, s3, e3, e3) = X3.t() * (X3.each_col() % gh.col(11));
H.submat(s3, s3, e3, e3) = X3.t() * schur_spmat(X3, gh.col(11));
// H.submat(s4, s3, e4, e3) = X4.t() * (X3.each_col() % gh.col(12));
H.submat(s4, s3, e4, e3) = X4.t() * schur_spmat(X3, gh.col(12));
H.submat(s3, s4, e3, e4) = H.submat(s4, s3, e4, e3).t();

// H.submat(s4, s4, e4, e4) = X4.t() * (X4.each_col() % gh.col(13));
H.submat(s4, s4, e4, e4) = X4.t() * schur_spmat(X4, gh.col(13));

out(1) = H;
    
}

// X1.each_col() %= gh.col(0);
// X2.each_col() %= gh.col(1);
// X3.each_col() %= gh.col(2);
// X4.each_col() %= gh.col(3);
X1 = schur_spmat(X1, gh.col(0));
X2 = schur_spmat(X2, gh.col(1));
X3 = schur_spmat(X3, gh.col(2));
X4 = schur_spmat(X4, gh.col(3));

if (sand == 0) {
    g = join_rows(join_rows(join_rows(sum(X1, 0), sum(X2, 0)), sum(X3, 0)), sum(X4, 0));
} else {
    g = join_rows(join_rows(join_rows(X1, X2), X3), X4);
}

out(0) = g;

return(out);

}

// [[Rcpp::export(.gHsp5)]]
Rcpp::List gHsp5(arma::mat gh, arma::sp_mat X1, arma::sp_mat X2, arma::sp_mat X3, arma::sp_mat X4, arma::sp_mat X5, const arma::uvec dupid, int dcate, int sand, int deriv)
{

Rcpp::List out(2);
arma::mat g;

// if (dcate == 1) {
//   X1 = X1.rows(dupid);
//   X2 = X2.rows(dupid);
//   X3 = X3.rows(dupid);
//   X4 = X4.rows(dupid);
//   X5 = X5.rows(dupid);
// }

if (deriv > 1) {

int s1 = 0;
int e1 = X1.n_cols - 1;
int s2 = e1 + 1;
int e2 = e1 + X2.n_cols;
int s3 = e2 + 1;
int e3 = e2 + X3.n_cols;
int s4 = e3 + 1;
int e4 = e3 + X4.n_cols;
int s5 = e4 + 1;
int e5 = e4 + X5.n_cols;

arma::mat H =  arma::mat(e5 + 1, e5 + 1);

// H.submat(s1, s1, e1, e1) = X1.t() * (X1.each_col() % gh.col(5));
H.submat(s1, s1, e1, e1) = X1.t() * schur_spmat(X1, gh.col(5));
// H.submat(s2, s1, e2, e1) = X2.t() * (X1.each_col() % gh.col(6));
H.submat(s2, s1, e2, e1) = X2.t() * schur_spmat(X1, gh.col(6));
// H.submat(s3, s1, e3, e1) = X3.t() * (X1.each_col() % gh.col(7));
H.submat(s3, s1, e3, e1) = X3.t() * schur_spmat(X1, gh.col(7));
// H.submat(s4, s1, e4, e1) = X4.t() * (X1.each_col() % gh.col(8));
H.submat(s4, s1, e4, e1) = X4.t() * schur_spmat(X1, gh.col(8));
// H.submat(s5, s1, e5, e1) = X5.t() * (X1.each_col() % gh.col(9));
H.submat(s5, s1, e5, e1) = X5.t() * schur_spmat(X1, gh.col(9));
H.submat(s1, s2, e1, e2) = H.submat(s2, s1, e2, e1).t();
H.submat(s1, s3, e1, e3) = H.submat(s3, s1, e3, e1).t();
H.submat(s1, s4, e1, e4) = H.submat(s4, s1, e4, e1).t();
H.submat(s1, s5, e1, e5) = H.submat(s5, s1, e5, e1).t();

// H.submat(s2, s2, e2, e2) = X2.t() * (X2.each_col() % gh.col(10));
H.submat(s2, s2, e2, e2) = X2.t() * schur_spmat(X2, gh.col(10));
// H.submat(s3, s2, e3, e2) = X3.t() * (X2.each_col() % gh.col(11));
H.submat(s3, s2, e3, e2) = X3.t() * schur_spmat(X2, gh.col(11));
// H.submat(s4, s2, e4, e2) = X4.t() * (X2.each_col() % gh.col(12));
H.submat(s4, s2, e4, e2) = X4.t() * schur_spmat(X2, gh.col(12));
// H.submat(s5, s2, e5, e2) = X5.t() * (X2.each_col() % gh.col(13));
H.submat(s5, s2, e5, e2) = X5.t() * schur_spmat(X2, gh.col(13));
H.submat(s2, s3, e2, e3) = H.submat(s3, s2, e3, e2).t();
H.submat(s2, s4, e2, e4) = H.submat(s4, s2, e4, e2).t();
H.submat(s2, s5, e2, e5) = H.submat(s5, s2, e5, e2).t();

// H.submat(s3, s3, e3, e3) = X3.t() * (X3.each_col() % gh.col(14));
H.submat(s3, s3, e3, e3) = X3.t() * schur_spmat(X3, gh.col(14));
// H.submat(s4, s3, e4, e3) = X4.t() * (X3.each_col() % gh.col(15));
H.submat(s4, s3, e4, e3) = X4.t() * schur_spmat(X3, gh.col(15));
// H.submat(s5, s3, e5, e3) = X5.t() * (X3.each_col() % gh.col(16));
H.submat(s5, s3, e5, e3) = X5.t() * schur_spmat(X3, gh.col(16));
H.submat(s3, s4, e3, e4) = H.submat(s4, s3, e4, e3).t();
H.submat(s3, s5, e3, e5) = H.submat(s5, s3, e5, e3).t();

// H.submat(s4, s4, e4, e4) = X4.t() * (X4.each_col() % gh.col(17));
H.submat(s4, s4, e4, e4) = X4.t() * schur_spmat(X4, gh.col(17));
// H.submat(s5, s4, e5, e4) = X5.t() * (X4.each_col() % gh.col(18));
H.submat(s5, s4, e5, e4) = X5.t() * schur_spmat(X4, gh.col(18));
H.submat(s4, s5, e4, e5) = H.submat(s5, s4, e5, e4).t();

// H.submat(s5, s5, e5, e5) = X5.t() * (X5.each_col() % gh.col(19));
H.submat(s5, s5, e5, e5) = X5.t() * schur_spmat(X5, gh.col(19));

out(1) = H;
    
}

// X1.each_col() %= gh.col(0);
// X2.each_col() %= gh.col(1);
// X3.each_col() %= gh.col(2);
// X4.each_col() %= gh.col(3);
// X5.each_col() %= gh.col(4);
X1 = schur_spmat(X1, gh.col(0));
X2 = schur_spmat(X2, gh.col(1));
X3 = schur_spmat(X3, gh.col(2));
X4 = schur_spmat(X4, gh.col(3));
X5 = schur_spmat(X5, gh.col(4));

if (sand == 0) {
    g = join_rows(join_rows(join_rows(join_rows(sum(X1, 0), sum(X2, 0)), sum(X3, 0)), sum(X4, 0)), sum(X5, 0));
} else {
    g = join_rows(join_rows(join_rows(join_rows(X1, X2), X3), X4), X5);
}

out(0) = g;

return(out);

}

// [[Rcpp::export(.gHsp6)]]
Rcpp::List gHsp6(arma::mat gh, arma::sp_mat X1, arma::sp_mat X2, arma::sp_mat X3, arma::sp_mat X4, arma::sp_mat X5, arma::sp_mat X6, const arma::uvec dupid, int dcate, int sand, int deriv)
{

Rcpp::List out(2);
arma::mat g;

// if (dcate == 1) {
//   X1 = X1.rows(dupid);
//   X2 = X2.rows(dupid);
//   X3 = X3.rows(dupid);
//   X4 = X4.rows(dupid);
//   X5 = X5.rows(dupid);
//   X6 = X6.rows(dupid);
// }

if (deriv > 1) {

int s1 = 0;
int e1 = X1.n_cols - 1;
int s2 = e1 + 1;
int e2 = e1 + X2.n_cols;
int s3 = e2 + 1;
int e3 = e2 + X3.n_cols;
int s4 = e3 + 1;
int e4 = e3 + X4.n_cols;
int s5 = e4 + 1;
int e5 = e4 + X5.n_cols;
int s6 = e5 + 1;
int e6 = e5 + X6.n_cols;

arma::mat H =  arma::mat(e6 + 1, e6 + 1);

// H.submat(s1, s1, e1, e1) = X1.t() * (X1.each_col() % gh.col(6));
H.submat(s1, s1, e1, e1) = X1.t() * schur_spmat(X1, gh.col(6));
// H.submat(s2, s1, e2, e1) = X2.t() * (X1.each_col() % gh.col(7));
H.submat(s2, s1, e2, e1) = X2.t() * schur_spmat(X1, gh.col(7));
// H.submat(s3, s1, e3, e1) = X3.t() * (X1.each_col() % gh.col(8));
H.submat(s3, s1, e3, e1) = X3.t() * schur_spmat(X1, gh.col(8));
// H.submat(s4, s1, e4, e1) = X4.t() * (X1.each_col() % gh.col(9));
H.submat(s4, s1, e4, e1) = X4.t() * schur_spmat(X1, gh.col(9));
// H.submat(s5, s1, e5, e1) = X5.t() * (X1.each_col() % gh.col(10));
H.submat(s5, s1, e5, e1) = X5.t() * schur_spmat(X1, gh.col(10));
// H.submat(s6, s1, e6, e1) = X6.t() * (X1.each_col() % gh.col(11));
H.submat(s6, s1, e6, e1) = X6.t() * schur_spmat(X1, gh.col(11));
H.submat(s1, s2, e1, e2) = H.submat(s2, s1, e2, e1).t();
H.submat(s1, s3, e1, e3) = H.submat(s3, s1, e3, e1).t();
H.submat(s1, s4, e1, e4) = H.submat(s4, s1, e4, e1).t();
H.submat(s1, s5, e1, e5) = H.submat(s5, s1, e5, e1).t();
H.submat(s1, s6, e1, e6) = H.submat(s6, s1, e6, e1).t();

// H.submat(s2, s2, e2, e2) = X2.t() * (X2.each_col() % gh.col(12));
H.submat(s2, s2, e2, e2) = X2.t() * schur_spmat(X2, gh.col(12));
// H.submat(s3, s2, e3, e2) = X3.t() * (X2.each_col() % gh.col(13));
H.submat(s3, s2, e3, e2) = X3.t() * schur_spmat(X2, gh.col(13));
// H.submat(s4, s2, e4, e2) = X4.t() * (X2.each_col() % gh.col(14));
H.submat(s4, s2, e4, e2) = X4.t() * schur_spmat(X2, gh.col(14));
// H.submat(s5, s2, e5, e2) = X5.t() * (X2.each_col() % gh.col(15));
H.submat(s5, s2, e5, e2) = X5.t() * schur_spmat(X2, gh.col(15));
// H.submat(s6, s2, e6, e2) = X6.t() * (X2.each_col() % gh.col(16));
H.submat(s6, s2, e6, e2) = X6.t() * schur_spmat(X2, gh.col(16));
H.submat(s2, s3, e2, e3) = H.submat(s3, s2, e3, e2).t();
H.submat(s2, s4, e2, e4) = H.submat(s4, s2, e4, e2).t();
H.submat(s2, s5, e2, e5) = H.submat(s5, s2, e5, e2).t();
H.submat(s2, s6, e2, e6) = H.submat(s6, s2, e6, e2).t();

// H.submat(s3, s3, e3, e3) = X3.t() * (X3.each_col() % gh.col(17));
H.submat(s3, s3, e3, e3) = X3.t() * schur_spmat(X3, gh.col(17));
// H.submat(s4, s3, e4, e3) = X4.t() * (X3.each_col() % gh.col(18));
H.submat(s4, s3, e4, e3) = X4.t() * schur_spmat(X3, gh.col(18));
// H.submat(s5, s3, e5, e3) = X5.t() * (X3.each_col() % gh.col(19));
H.submat(s5, s3, e5, e3) = X5.t() * schur_spmat(X3, gh.col(19));
// H.submat(s6, s3, e6, e3) = X6.t() * (X3.each_col() % gh.col(20));
H.submat(s6, s3, e6, e3) = X6.t() * schur_spmat(X3, gh.col(20));
H.submat(s3, s4, e3, e4) = H.submat(s4, s3, e4, e3).t();
H.submat(s3, s5, e3, e5) = H.submat(s5, s3, e5, e3).t();
H.submat(s3, s6, e3, e6) = H.submat(s6, s3, e6, e3).t();

// H.submat(s4, s4, e4, e4) = X4.t() * (X4.each_col() % gh.col(21));
H.submat(s4, s4, e4, e4) = X4.t() * schur_spmat(X4, gh.col(21));
// H.submat(s5, s4, e5, e4) = X5.t() * (X4.each_col() % gh.col(22));
H.submat(s5, s4, e5, e4) = X5.t() * schur_spmat(X4, gh.col(22));
// H.submat(s6, s4, e6, e4) = X6.t() * (X4.each_col() % gh.col(23));
H.submat(s6, s4, e6, e4) = X6.t() * schur_spmat(X4, gh.col(23));
H.submat(s4, s5, e4, e5) = H.submat(s5, s4, e5, e4).t();
H.submat(s4, s6, e4, e6) = H.submat(s6, s4, e6, e4).t();

// H.submat(s5, s5, e5, e5) = X5.t() * (X5.each_col() % gh.col(24));
H.submat(s5, s5, e5, e5) = X5.t() * schur_spmat(X5, gh.col(24));
// H.submat(s6, s5, e6, e5) = X6.t() * (X5.each_col() % gh.col(25));
H.submat(s6, s5, e6, e5) = X6.t() * schur_spmat(X5, gh.col(25));
H.submat(s5, s6, e5, e6) = H.submat(s6, s5, e6, e5).t();

// H.submat(s6, s6, e6, e6) = X6.t() * (X6.each_col() % gh.col(26));
H.submat(s6, s6, e6, e6) = X6.t() * schur_spmat(X6, gh.col(26));

out(1) = H;
    
}

// X1.each_col() %= gh.col(0);
// X2.each_col() %= gh.col(1);
// X3.each_col() %= gh.col(2);
// X4.each_col() %= gh.col(3);
// X5.each_col() %= gh.col(4);
// X6.each_col() %= gh.col(5);
X1 = schur_spmat(X1, gh.col(0));
X2 = schur_spmat(X2, gh.col(1));
X3 = schur_spmat(X3, gh.col(2));
X4 = schur_spmat(X4, gh.col(3));
X5 = schur_spmat(X5, gh.col(4));
X6 = schur_spmat(X6, gh.col(5));

if (sand == 0) {
    g = join_rows(join_rows(join_rows(join_rows(join_rows(sum(X1, 0), sum(X2, 0)), sum(X3, 0)), sum(X4, 0)), sum(X5, 0)), sum(X6, 0));
} else {
    g = join_rows(join_rows(join_rows(join_rows(join_rows(X1, X2), X3), X4), X5), X6);
}

out(0) = g;

return(out);

}
