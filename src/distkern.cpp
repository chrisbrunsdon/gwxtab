#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

double square(double d) {
  return d*d;
}

double hav(double theta) {
  return square(sin(theta/2.0));
}


// [[Rcpp::export]]
NumericMatrix dist2 (NumericMatrix x1, NumericMatrix x2){
  int n1 = x1.nrow();
  int n2 = x2.nrow();
  NumericMatrix dist(n1,n2);
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++) {
      dist(i,j) = square(x1(i,0) - x2(j,0));
      dist(i,j) += square(x1(i,1) - x2(j,1));
    }
  return dist;
}


// [[Rcpp::export]]
NumericMatrix bisq2(NumericMatrix x1, NumericMatrix x2, double h){
  int n1 = x1.nrow();
  int n2 = x2.nrow();
  double dist;
  double h2 = square(h);
  NumericMatrix w(n1,n2);
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++) {
      dist = square(x1(i,0) - x2(j,0));
      dist += square(x1(i,1) - x2(j,1));
      if (dist > h2) {
        w(i,j) = 0.0;
      }
      else {
        w(i,j) = square(1 - dist/h2);
      }
    }
    return w;
}

// [[Rcpp::export]]
NumericMatrix bisq2ad(NumericMatrix x1, NumericMatrix x2, NumericVector h){
  int n1 = x1.nrow();
  int n2 = x2.nrow();
  double dist;
  NumericVector h2 = h;
  for (int j = 0; j < h.length(); j++) h2(j) = square(h(j));
  NumericMatrix w(n1,n2);
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++) {
      dist = square(x1(i,0) - x2(j,0));
      dist += square(x1(i,1) - x2(j,1));
      if (dist > h2(j)) {
        w(i,j) = 0.0;
      }
      else {
        w(i,j) = square(1 - dist/h2(j));
      }
    }
    return w;
}


// [[Rcpp::export]]
NumericVector bisq(NumericMatrix x1, NumericVector x2, double h){
  int n = x1.nrow();
  double dist;
  double h2 = square(h);
  NumericVector w(n);
  for (int i = 0; i < n; i++) {
      dist = square(x1(i,0) - x2(0));
      dist += square(x1(i,1) - x2(1));
      if (dist > h2) {
        w(i) = 0.0;
      }
      else {
        w(i) = square(1 - dist/h2);
      }
    }
    return w;
}


// [[Rcpp::export]]
NumericMatrix gcdist2 (NumericMatrix x1, NumericMatrix x2){
  double d2r = M_PI/180.0;
  int n1 = x1.nrow();
  int n2 = x2.nrow();
  double d0, d1, lat1, lat2, lon1, lon2, a;
  NumericMatrix dist(n1,n2);
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++) {
      lat1 = x1(i,1) * d2r;
      lat2 = x2(j,1) * d2r;
      lon1 = x1(i,0) * d2r;
      lon2 = x2(j,0) * d2r;
      dist(i,j)=acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon1-lon2));
      d0 = std::abs(lon1 - lon2);
      d1 = std::abs(lat1 - lat2);
      a = hav(d1) + cos(lat1) * cos(lat2)*hav(d0);
      dist(i,j) = 2*atan2(sqrt(a),sqrt(1.0-a));
      // dist(i,j) = 2.0*asin(hav(d1) + cos(lat1)*cos(lat2)*hav(d0));
    }
    return dist;
}




