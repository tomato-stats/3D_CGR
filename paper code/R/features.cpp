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
  double x = (dot_product(a, b) / (magnitude(a) * magnitude(b)));
  if(x > 1){
    x = 1.0;
  } else if(x < -1){
    x = -1.0;
  }
  return acos(x);
}

// Function to calculate Euclidean distance between two vectors
double euclidean_distance(NumericVector x, NumericVector y) {
  double dist = 0;
  for(int i = 0; i < x.length(); i++) {
    dist += pow(x[i] - y[i], 2);
  }
  return sqrt(dist);
}

// [[Rcpp::export]]
NumericVector by3rowangle(NumericMatrix points) {
  int n = points.nrow();
  NumericVector angles(n - 2);
  
  for (int i = 0; i < n - 2; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    NumericVector p3 = points.row(i + 2);
    
    NumericVector a = p1 - p2;
    NumericVector b = p3 - p2;
    angles[i] = angle(a, b);
  }
  
  return angles;
}

// [[Rcpp::export]]
NumericVector by3rowdistance(NumericMatrix points) {
  int n = points.nrow();
  int m = points.ncol(); 
  NumericVector distances(n - 2);
  int k = 0;
  int j; 
  
  if(m == 3){
    for (int i = 0; i < n - 3; ++i) {
      j = i + 2;
      double di = points(i, 0) - points(j, 0);
      double dj = points(i, 1) - points(j, 1);
      double dk = points(i, 2) - points(j, 2);
      distances[k++] = sqrt(di * di + dj * dj + dk * dk);
    } 
  } else if(m == 4) {
    for (int i = 0; i < n - 3; ++i) {
      j = i + 2;      
      double dr = points(i, 0) - points(j, 0);
      double di = points(i, 0) - points(j, 0);
      double dj = points(i, 1) - points(j, 1);
      double dk = points(i, 2) - points(j, 2);
      distances[k++] = sqrt(dr * dr + di * di + dj * dj + dk * dk);
    }
  } else{
    for (int i = 0; i < n - 3; ++i) {
      j = i + 2;
      double dist = 0;
      for (int p = 0; p < m; ++p) {
        double d = points(i, p) - points(j, p);
        dist += d * d;
      }
      distances[k++] = sqrt(dist);
    }
  }
return distances;
}


// [[Rcpp::export]]
NumericVector vertices_dist(NumericMatrix points, const NumericMatrix &baseCoords) {
  int n = points.nrow();
  int m = points.ncol(); 
  NumericMatrix distances(n, baseCoords.ncol());
  int k = 0;
  int j; 
  
  if(m == 3){
    for (int i = 0; i < n - 3; ++i) {
      j = i + 2;
      double di = points(i, 0) - points(j, 0);
      double dj = points(i, 1) - points(j, 1);
      double dk = points(i, 2) - points(j, 2);
      distances[k++] = sqrt(di * di + dj * dj + dk * dk);
      
      
    } 
  } else if(m == 4) {
    for (int i = 0; i < n - 3; ++i) {
      j = i + 2;      
      double dr = points(i, 0) - points(j, 0);
      double di = points(i, 0) - points(j, 0);
      double dj = points(i, 1) - points(j, 1);
      double dk = points(i, 2) - points(j, 2);
      distances[k++] = sqrt(dr * dr + di * di + dj * dj + dk * dk);
    }
  } else{
    for (int i = 0; i < n - 3; ++i) {
      j = i + 2;
      double dist = 0;
      for (int p = 0; p < m; ++p) {
        double d = points(i, p) - points(j, p);
        dist += d * d;
      }
      distances[k++] = sqrt(dist);
    }
  }
  return distances;
}


// [[Rcpp::export]]
NumericMatrix distance_matrix(NumericMatrix mat1, NumericMatrix mat2) {
  int nrow1 = mat1.nrow();
  int nrow2 = mat2.nrow();
  
  NumericMatrix out(nrow1, nrow2);
  
  for(int i = 0; i < nrow1; i++) {
    for(int j = 0; j < nrow2; j++) {
      NumericVector vec1 = mat1.row(i);
      NumericVector vec2 = mat2.row(j);
      out(i, j) = euclidean_distance(vec1, vec2);
    }
  }
  
  return out;
}





