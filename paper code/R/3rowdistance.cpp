#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector by3rowdistance(NumericMatrix points) {
  int n = points.nrow();
  NumericVector distances(n - 2);
  int k = 0;
  int j; 
  
  for (int i = 0; i < n - 3; ++i) {
    j = i + 2;
    double dx = points(i, 0) - points(j, 0);
    double dy = points(i, 1) - points(j, 1);
    double dz = points(i, 2) - points(j, 2);
    distances[k++] = sqrt(dx * dx + dy * dy + dz * dz);
  }

return distances;
}



/*
 NumericVector pairwise_distances(DataFrame df) {
 int n = df.nrow();
 NumericVector x = df["x"], y = df["y"], z = df["z"];
 NumericVector distances(n * (n - 1) / 2);
 int k = 0;
 
 for (int i = 0; i < n; ++i) {
 for (int j = i + 1; j < n; ++j) {
 double dx = x[i] - x[j];
 double dy = y[i] - y[j];
 double dz = z[i] - z[j];
 distances[k++] = sqrt(dx * dx + dy * dy + dz * dz);
 }
 }
 
 return distances;
 }
*/