#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;


// Calculates the cross product of two vectors 
NumericVector cross_product(NumericVector a, NumericVector b) {
  NumericVector cross(3);
  cross[0] = a[1]*b[2] - a[2]*b[1];
  cross[1] = a[2]*b[0] - a[0]*b[2];
  cross[2] = a[0]*b[1] - a[1]*b[0];
  return cross;
}

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

// Calculates the oriented angle between two vectors
double oriented_angle(NumericVector a, NumericVector b, NumericVector n) {
  double x =  (dot_product(a, b) / (magnitude(a) * magnitude(b)));
  if(x > 1){
    x = 1.0;
  } else if(x < -1){
    x = -1.0;
  }
  double theta = std::acos(x);
  double sign = dot_product(n, cross_product(a, b)) < 0 ? -1.0 : 1.0;
  return sign * theta;
}

double oriented_distance(NumericVector a, NumericVector b, NumericVector n) {
  NumericVector ba = b - a;
  NumericVector normalized_n = n / magnitude(n);  
  return dot_product(ba, normalized_n);
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
  // double dot_prod;
  
  for (int i = 0; i < n - 2; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    NumericVector p3 = points.row(i + 2);
    
    NumericVector a = p1 - p2;
    NumericVector b = p3 - p2;
    // dot_prod = dot_product(a, b);
    angles[i] = angle(a, b);
    // if(Rcpp::NumericVector::is_na(dot_prod)) angles[i] = angle(a, b) ;
    // else if(dot_prod < 0) angles[i] = angle(a, b) * -1;
  }
  
  return angles;
}

// [[Rcpp::export]]
NumericVector by3roworientedangle(NumericMatrix points, NumericVector v = NumericVector::create(1.0, 0.0, 0.0)) {
  int n = points.nrow();
  NumericVector angles(n - 2);
  // double dot_prod;
  
  for (int i = 0; i < n - 2; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    NumericVector p3 = points.row(i + 2);
    
    NumericVector a = p1 - p2;
    NumericVector b = p3 - p2;
    //angles[i] = oriented_angle(a, b, cross_product(a, b));
    angles[i] = oriented_angle(a, b, p2);
  }
  
  return angles;
}

// [[Rcpp::export]]
NumericVector by3roworienteddistance(NumericMatrix points, NumericVector v = NumericVector::create(1.0, 0.0, 0.0)) {
  int n = points.nrow();
  NumericVector distances(n - 2);
  // double dot_prod;
  
  for (int i = 0; i < n - 2; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    NumericVector p3 = points.row(i + 2);
    
    NumericVector a = p1 - p2;
    NumericVector b = p3 - p2;
    distances[i] = oriented_distance(a, b, v);
  }
  
  return distances;
}

// [[Rcpp::export]]
NumericVector orienteddistance(NumericMatrix points) {
  int n = points.nrow();
  NumericVector distances(n - 1);
  // double dot_prod;
  
  for (int i = 1; i < n - 1; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);

    distances[i] = oriented_distance(p1, p2, cross_product(p1, p2));
  }
  
  return distances;
}

// [[Rcpp::export]]
NumericVector by3rowdotprod(NumericMatrix points) {
  int n = points.nrow();
  NumericVector dotprod(n - 2);
  
  for (int i = 0; i < n - 2; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    NumericVector p3 = points.row(i + 2);
    
    NumericVector a = p1 - p2;
    NumericVector b = p3 - p2;
    dotprod[i] = dot_product(a, b);
  }
  
  return dotprod;
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
NumericVector distances_from_origin(NumericMatrix points) {
  int n = points.nrow();
  int m = points.ncol(); 
  NumericVector distances(n);
  int k = 0;
  
  if(m == 3){
    for (int i = 0; i < n; ++i) {
      double di = points(i, 0);
      double dj = points(i, 1);
      double dk = points(i, 2);
      distances[k++] = sqrt(di * di + dj * dj + dk * dk);
    } 
  } else if(m == 4) {
    for (int i = 0; i < n; ++i) {
      double dr = points(i, 0);
      double di = points(i, 0);
      double dj = points(i, 1);
      double dk = points(i, 2);
      distances[k++] = sqrt(dr * dr + di * di + dj * dj + dk * dk);
    }
  } else{
    for (int i = 0; i < n; ++i) {
      double dist = 0;
      for (int p = 0; p < m; ++p) {
        double d = points(i, p);
        dist += d * d;
      }
      distances[k++] = sqrt(dist);
    }
  }
  return distances;
}

// [[Rcpp::export]]
NumericMatrix xy_from_3rowangle(NumericMatrix points){
  int n = points.nrow();
  NumericMatrix output(n-2, 2);
  
  for (int i = 0; i < n - 2; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    NumericVector p3 = points.row(i + 2);
    
    NumericVector a = p1 - p2;
    NumericVector b = p3 - p2;
    double angle_val = angle(a, b);
    output(i, 0) = cos(angle_val);
    output(i, 1) = sin(angle_val);
  }
  return output;
}




// [[Rcpp::export]]
NumericVector by3rowscalartripleproduct(NumericMatrix points) {
  int n = points.nrow();
  int veclen = points.ncol();
  NumericVector dotprod(n - 2);
  NumericVector c(veclen, 0.0);
  
  c[0] = 1;
  for (int i = 0; i < n - 2; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    NumericVector p3 = points.row(i + 2);
    
    NumericVector a = p1 - p2;
    NumericVector b = p3 - p2;
    dotprod[i] = dot_product(cross_product(a, b), c);
  }
  
  return dotprod;
}

// [[Rcpp::export]]
NumericVector scalartripleproduct1(NumericMatrix points) {
  int n = points.nrow();
  NumericVector dotprod(n - 2);

  for (int i = 0; i < n - 2; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    NumericVector p3 = points.row(i + 2);
    
    dotprod[i] = dot_product(cross_product(p1, p2), p3);
  }
  
  return dotprod;
}

// [[Rcpp::export]]
NumericVector scalartripleproduct2(NumericMatrix points) {
  int n = points.nrow();
  int veclen = points.ncol();
  NumericVector dotprod(n - 1);
  NumericVector c(veclen, 0.0);
  
  c[0] = 1;
  
  for (int i = 0; i < n - 1; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    
    dotprod[i] = dot_product(cross_product(p1, p2), c);
  }
  
  return dotprod;
}

// [[Rcpp::export]]
NumericVector scalartripleproduct3(NumericMatrix points) {
  int n = points.nrow();
  int veclen = points.ncol();
  NumericVector dotprod(n - 1);
  NumericVector c(veclen, 0.0);
  
  c[1] = 1;
  
  for (int i = 0; i < n - 1; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    
    dotprod[i] = dot_product(cross_product(p1, p2), c);
  }
  
  return dotprod;
}

// [[Rcpp::export]]
NumericVector scalartripleproduct4(NumericMatrix points) {
  int n = points.nrow();
  int veclen = points.ncol();
  NumericVector dotprod(n - 1);
  NumericVector c(veclen, 0.0);
  
  c[2] = 1;
  
  for (int i = 0; i < n - 1; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    
    dotprod[i] = dot_product(cross_product(p1, p2), c);
  }
  
  return dotprod;
}

// [[Rcpp::export]]
NumericVector scalartripleproduct5(NumericMatrix points) {
  int n = points.nrow();
  int veclen = points.ncol();
  NumericVector dotprod(n - 1);
  NumericVector c(veclen, 0.0);
  
  c[0] = -0.5 * 2 * std::sqrt(1.0/3.0);
  c[1] = -0.5 * 2 * std::sqrt(1.0/3.0);
  c[2] = -0.5 * 2 * std::sqrt(1.0/3.0);
  
  for (int i = 0; i < n - 1; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    
    dotprod[i] = dot_product(cross_product(p1, p2), c);
  }
  
  return dotprod;
}
// [[Rcpp::export]]
NumericVector scalartripleproduct6(NumericMatrix points) {
  int n = points.nrow();
  int veclen = points.ncol();
  NumericVector dotprod(n - 1);
  NumericVector c(veclen, 0.0);
  
  c[0] = (1.0-0.5) * 2 * std::sqrt(1.0/3.0);
  c[1] = -0.5 * 2 * std::sqrt(1.0/3.0);
  c[2] = (1.0-0.5) * 2 * std::sqrt(1.0/3.0);
  
  for (int i = 0; i < n - 1; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    
    dotprod[i] = dot_product(cross_product(p1, p2), c);
  }
  
  return dotprod;
}
// [[Rcpp::export]]
NumericVector scalartripleproduct7(NumericMatrix points) {
  int n = points.nrow();
  int veclen = points.ncol();
  NumericVector dotprod(n - 1);
  NumericVector c(veclen, 0.0);
  
  c[0] = -0.5 * 2 * std::sqrt(1.0/3.0);
  c[1] = (1.0-0.5) * 2 * std::sqrt(1.0/3.0);
  c[2] = (1.0-0.5) * 2 * std::sqrt(1.0/3.0);
  
  for (int i = 0; i < n - 1; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    
    dotprod[i] = dot_product(cross_product(p1, p2), c);
  }
  
  return dotprod;
}
// [[Rcpp::export]]
NumericVector scalartripleproduct8(NumericMatrix points) {
  int n = points.nrow();
  int veclen = points.ncol();
  NumericVector dotprod(n - 1);
  NumericVector c(veclen, 0.0);
  
  c[0] = (1.0-0.5) * 2 * std::sqrt(1.0/3.0);
  c[1] = (1.0-0.5) * 2 * std::sqrt(1.0/3.0);
  c[2] = -0.5 * 2 * std::sqrt(1.0/3.0);
  
  for (int i = 0; i < n - 1; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    
    dotprod[i] = dot_product(cross_product(p1, p2), c);
  }
  
  return dotprod;
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





