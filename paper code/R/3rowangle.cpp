#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Calculates the dot product of two vectors
double dot_product(NumericVector a, NumericVector b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Calculates the magnitude of a vector
double magnitude(NumericVector a) {
  return sqrt(dot_product(a, a));
}

// Calculates the angle between two vectors
double angle(NumericVector a, NumericVector b) {
  return acos(dot_product(a, b) / (magnitude(a) * magnitude(b)));
}

// [[Rcpp::export]]
NumericVector sliding_window_angles(NumericMatrix points) {
  int n = points.nrow();
  NumericVector angles(n - 2);
  double dot_prod;
  
  for (int i = 0; i < n - 2; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    NumericVector p3 = points.row(i + 2);
    
    NumericVector a = p2 - p1;
    NumericVector b = p3 - p2;
    dot_prod = dot_product(a, b);
    angles[i] = angle(a, b);
  }
  
  return angles;
}
