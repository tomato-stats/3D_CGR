
#include <cmath>
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]


//-----------------------------------------
// Calculations
//-----------------------------------------

// [[Rcpp::export]]
// Calculates the angle between two vectors
double angle(arma::rowvec a, arma::rowvec b) {
  double x = (arma::dot(a, b) / (arma::norm(a) * arma::norm(b)));
  if(x > 1){
    x = 1.0;
  } else if(x < -1){
    x = -1.0;
  }
  return acos(x);
}

// [[Rcpp::export]]
// Calculates the oriented angle between two vectors
// For our particular application, when the norm is zero, we set the angle between points to be zero, This can literally only happen when there are many consecutive identical nucleotides in the sequence. 
double oriented_angle(arma::rowvec  a, arma::rowvec  b, arma::rowvec  n) {
  if (arma::norm(a) * arma::norm(b) == 0.0) {
    return 0.0;
  }
  
  double x =  (arma::dot(a, b) / (arma::norm(a) *  arma::norm(b)));
  if(x > 1){
    x = 1.0;
  } else if(x < -1){
    x = -1.0;
  }
  double theta = std::acos(x);
  double sign = arma::dot(n, arma::cross(a, b)) < 0 ? -1.0 : 1.0;
  return sign * theta;
}

// [[Rcpp::export]]
double oriented_distance(arma::rowvec a, arma::rowvec b, arma::rowvec n) {
  double dist = sqrt(sum(pow(b - a, 2.0)));
  double sign = arma::dot(b-a, n) < 0 ? -1.0 : 1.0;
  return sign * dist; 
}


//-----------------------------------------
// Shape features 
//-----------------------------------------

// angles and edge lengths in a 3-point sliding window
// These functions are still here, but have essentially been superceded by 
// oriented measurements. Likely to be removed soon. 

// [[Rcpp::export]]
NumericVector by3rowangle(arma::mat points) {
  int n = points.n_rows;
  NumericVector angles(n - 2);
  
  for (int i = 0; i < n - 2; ++i) {
    arma::rowvec p1 = points.row(i);
    arma::rowvec p2 = points.row(i + 1);
    arma::rowvec p3 = points.row(i + 2);
    
    arma::rowvec a = p1 - p2;
    arma::rowvec b = p3 - p2;
    angles[i] = angle(a, b);
  }
  
  return angles;
}

// [[Rcpp::export]]
NumericVector by3rowdistance(NumericMatrix points) {
  int n = points.nrow();
  NumericVector distances(n - 2);
  int k = 0;
  
  for(int i = 0; i < n-2; ++i) {
    NumericVector point1 = points(i, _); // get coordinates of point1
    NumericVector point3 = points(i+2, _); // get coordinates of point3
    
    double dist = sqrt(sum(pow(point1 - point3, 2.0))); // calculate Euclidean distance
    distances[k++] = dist; // store distance in vector
  }
  return distances;
}

// oriented angles and edge lengths in a 3-point sliding window


// [[Rcpp::export]]
NumericVector by3roworientedangle1(arma::mat points) {
  // oriented angle in a 3-point sliding window relative to the vector pointing to the vertex
  int n = points.n_rows;
  NumericVector angles(n - 2);
  
  for (int i = 0; i < n - 2; ++i) {
    arma::rowvec p1 = points.row(i);
    arma::rowvec p2 = points.row(i + 1);
    arma::rowvec p3 = points.row(i + 2);
    
    arma::rowvec a = p1 - p2;
    arma::rowvec b = p3 - p2;
    //angles[i] = oriented_angle(a, b, arma::cross(a, b)); this removes the direction component
    angles[i] = oriented_angle(a, b, p2);
  }
  
  return angles;
}

// [[Rcpp::export]]
NumericVector by3roworientedangle2(arma::mat points, arma::rowvec v) {
  // oriented angles between every point and two points away through the origin 
  // relative to vector v
  int n = points.n_rows;
  NumericVector angles(n - 2);
  
  for (int i = 0; i < n - 2; ++i) {
    arma::rowvec p1 = points.row(i);
    arma::rowvec p3 = points.row(i + 2);

    angles[i] = oriented_angle(p1, p3, v);
  }
  
  return angles;
}

// [[Rcpp::export]]
NumericVector by3roworientedangle3(arma::mat points) {
  int n = points.n_rows;
  NumericVector angles(n - 2);
  
  for (int i = 0; i < n - 2; ++i) {
    arma::rowvec p1 = points.row(i);
    arma::rowvec p2 = points.row(i + 1);
    arma::rowvec p3 = points.row(i + 2);

    angles[i] = oriented_angle(p1, p3, p2);
  }
  
  return angles;
}

// [[Rcpp::export]]
NumericVector by3roworientedangle4(arma::mat points, arma::rowvec v) {
  int n = points.n_rows;
  NumericVector angles(n - 2);
  
  for (int i = 0; i < n - 2; ++i) {
    arma::rowvec p1 = points.row(i);
    arma::rowvec p2 = points.row(i + 1);
    arma::rowvec p3 = points.row(i + 2);
    
    arma::rowvec a = p1 - p2;
    arma::rowvec b = p3 - p2;
    angles[i] = oriented_angle(a, b, v);
  }
  
  return angles;
}

// [[Rcpp::export]]
NumericVector by3roworienteddistance1(NumericMatrix points) {
  int n = points.nrow();
  NumericVector distances(n - 2);
  
  for (int i = 0; i < n - 2; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    NumericVector p3 = points.row(i + 2);
    
    NumericVector a = p1 - p2;
    NumericVector b = p3 - p2;
    distances[i] = oriented_distance(a, b, p2);
  }
  
  return distances;
}

// [[Rcpp::export]]
NumericVector by3roworienteddistance2(NumericMatrix points, NumericVector v = NumericVector::create(1.0, 0.0, 0.0)) {
  int n = points.nrow();
  NumericVector distances(n - 2);
  
  for (int i = 0; i < n - 2; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p3 = points.row(i + 2);

    distances[i] = oriented_distance(p1, p3, v);
  }
  
  return distances;
}


// [[Rcpp::export]]
NumericVector by3roworienteddistance3(NumericMatrix points) {
  int n = points.nrow();
  NumericVector angles(n - 2);

  for (int i = 0; i < n - 2; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    NumericVector p3 = points.row(i + 2);

    angles[i] = oriented_distance(p1, p3, p2);
  }

  return angles;
}


// [[Rcpp::export]]
NumericVector by3roworienteddistance4(NumericMatrix points, NumericVector v = NumericVector::create(1.0, 0.0, 0.0)) {
  int n = points.nrow();
  NumericVector distances(n - 2);
  
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
    arma::rowvec p1 = points.row(i);
    arma::rowvec p2 = points.row(i + 1);

    distances[i] = oriented_distance(p1, p2, arma::cross(p1, p2));
  }
  
  return distances;
}

// [[Rcpp::export]]
NumericVector orienteddistance1(NumericMatrix points, NumericVector v = NumericVector::create(1.0, 0.0, 0.0)) {
  int n = points.nrow();
  NumericVector distances(n - 1);
  
  for (int i = 1; i < n - 1; ++i) {
    NumericVector p1 = points.row(i);
    NumericVector p2 = points.row(i + 1);
    
    distances[i] = oriented_distance(p1, p2, v);
  }
  
  return distances;
}

// [[Rcpp::export]]
NumericVector orientedangle1(arma::mat points, arma::rowvec v) {
  int n = points.n_rows;
  NumericVector angles(n - 1);
  
  for (int i = 0; i < n - 1; ++i) {
    arma::rowvec p1 = points.row(i);
    arma::rowvec p2 = points.row(i + 1);
    
    angles[i] = oriented_angle(p1, p2, v);
  }
  
  return angles;
}


// [[Rcpp::export]]
NumericVector by2rowdotprod(NumericMatrix points) {
  int n = points.nrow();
  NumericVector dotprod(n - 1);
  
  for (int i = 0; i < n - 1; ++i) {
    arma::rowvec p1 = points.row(i);
    arma::rowvec p2 = points.row(i + 1);

    dotprod[i] = arma::dot(p1, p2);
  }
  
  return dotprod;
}

// [[Rcpp::export]]
NumericVector by3rowdotprod(NumericMatrix points) {
  int n = points.nrow();
  NumericVector dotprod(n - 2);
  
  for (int i = 0; i < n - 2; ++i) {
    arma::rowvec p1 = points.row(i);
    arma::rowvec p2 = points.row(i + 1);
    arma::rowvec p3 = points.row(i + 2);
    
    arma::rowvec a = p1 - p2;
    arma::rowvec b = p3 - p2;
    dotprod[i] = arma::dot(a, b);
  }
  
  return dotprod;
}


// [[Rcpp::export]]
NumericVector distances_from_origin(NumericMatrix points) {
  int n = points.nrow();
  NumericVector distances(n);
  int k = 0;
  
  for (int i = 0; i < n; ++i) {
    double dist = sqrt(sum(pow(points(i, _), 2.0)));
    distances[k++] = dist;
  }
  return distances;
}

// [[Rcpp::export]]
NumericVector distances_from_point(NumericMatrix points,  NumericVector v = NumericVector::create(1.0, 0.0, 0.0, 0.0)) {
  int n = points.nrow();
  int m = points.ncol(); 
  NumericVector distances(n);
  int k = 0;
  
  if(m == 3){
    for (int i = 0; i < n; ++i) {
      double di = points(i, 0) - v[0];
      double dj = points(i, 1) - v[1];
      double dk = points(i, 2) - v[2];
      distances[k++] = sqrt(di * di + dj * dj + dk * dk);
    } 
  } else if(m == 4) {
    for (int i = 0; i < n; ++i) {
      double dr = points(i, 0) - v[0];
      double di = points(i, 1) - v[1];
      double dj = points(i, 2) - v[2];
      double dk = points(i, 3) - v[3];
      distances[k++] = sqrt(dr * dr + di * di + dj * dj + dk * dk);
    }
  } else{
    for (int i = 0; i < n; ++i) {
      double dist = 0;
      for (int p = 0; p < m; ++p) {
        double d = points(i, p) - v[p];
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
  
  for (int i = 0; i < n - 2; ++i) {
    arma::rowvec p1 = points.row(i);
    arma::rowvec p2 = points.row(i + 1);
    arma::rowvec p3 = points.row(i + 2);
    
    arma::rowvec a = p1 - p2;
    arma::rowvec b = p3 - p2;
    dotprod[i] = arma::dot(arma::cross(a, b), p2);
  }
  
  return dotprod;
}

// NumericVector by3rowscalartripleproduct(NumericMatrix points) {
//   int n = points.nrow();
//   int veclen = points.ncol();
//   NumericVector dotprod(n - 2);
//   NumericVector c(veclen, 0.0);
//   
//   c[0] = 1;
//   for (int i = 0; i < n - 2; ++i) {
//     NumericVector p1 = points.row(i);
//     NumericVector p2 = points.row(i + 1);
//     NumericVector p3 = points.row(i + 2);
//     
//     NumericVector a = p1 - p2;
//     NumericVector b = p3 - p2;
//     dotprod[i] = dot_product(cross_product(a, b), c);
//   }
//   
//   return dotprod;
// }

// [[Rcpp::export]]
NumericVector scalartripleproduct_pairs(arma::mat points, arma::rowvec v) {
  int n = points.n_rows;
  NumericVector dotprod(n - 1);
  
  for (int i = 0; i < n - 1 ; ++i) {
    arma::rowvec p1 = points.row(i);
    arma::rowvec  p2 = points.row(i + 1);
    
    dotprod[i] = arma::dot(arma::cross(p1, p2), v);
  }
  
  return dotprod;
}

// [[Rcpp::export]]
NumericVector scalartripleproduct0(arma::mat points, arma::rowvec v) {
  int n = points.n_rows;
  NumericVector dotprod(n - 1);
  
  for (int i = 0; i < n - 1 ; ++i) {
    arma::rowvec  p1 = points.row(i);
    arma::rowvec  p2 = points.row(i + 1);
    
    dotprod[i] = arma::dot(arma::cross(p1, p2), v);
  }
  
  return dotprod;
}

// [[Rcpp::export]]
NumericVector scalartripleproduct1(arma::mat points) {
  int n = points.n_rows;
  NumericVector dotprod(n - 2);

  for (int i = 0; i < n - 2; ++i) {
    arma::rowvec p1 = points.row(i);
    arma::rowvec p2 = points.row(i + 1);
    arma::rowvec p3 = points.row(i + 2);
    
    dotprod[i] = arma::dot(arma::cross(p1, p3), p2);
  }
  
  return dotprod;
}

// [[Rcpp::export]]
NumericVector scalartripleproduct2(arma::mat points, arma::rowvec c) {
  int n = points.n_rows;
  NumericVector dotprod(n - 1);
  
  for (int i = 0; i < n - 1; ++i) {
    arma::rowvec p1 = points.row(i);
    arma::rowvec p2 = points.row(i + 1);
    
    dotprod[i] = arma::dot(arma::cross(p1, p2), c);
  }
  
  return dotprod;
}

// [[Rcpp::export]]
arma::rowvec scalartripleproduct5(NumericMatrix points) {
  int n = points.nrow();
  int veclen = points.ncol();
  arma::vec dotprod(n - 1);
  arma::rowvec c = arma::zeros<arma::rowvec>(veclen);
  
  c[0] = -0.5 * 2 * std::sqrt(1.0/3.0);
  c[1] = -0.5 * 2 * std::sqrt(1.0/3.0);
  c[2] = -0.5 * 2 * std::sqrt(1.0/3.0);
  
  for (int i = 0; i < n - 1; ++i) {
    arma::rowvec p1 = points.row(i);
    arma::rowvec p2 = points.row(i + 1);
    
    dotprod[i] = arma::dot(arma::cross(p1, p2), c);
  }
  
  return dotprod;
}

// [[Rcpp::export]]
arma::vec scalartripleproduct6(NumericMatrix points) {
  int n = points.nrow();
  int veclen = points.ncol();
  arma::vec dotprod(n - 1);
  arma::rowvec c = arma::zeros<arma::rowvec>(veclen);
  
  c[0] = (1.0-0.5) * 2 * std::sqrt(1.0/3.0);
  c[1] = -0.5 * 2 * std::sqrt(1.0/3.0);
  c[2] = (1.0-0.5) * 2 * std::sqrt(1.0/3.0);
  
  for (int i = 0; i < n - 1; ++i) {
    arma::rowvec p1 = points.row(i);
    arma::rowvec p2 = points.row(i + 1);
    
    dotprod[i] = arma::dot(arma::cross(p1, p2), c);
  }
  
  return dotprod;
}

// [[Rcpp::export]]
arma::vec scalartripleproduct7(NumericMatrix points) {
  int n = points.nrow();
  int veclen = points.ncol();
  arma::vec dotprod(n - 1);
  arma::rowvec c = arma::zeros<arma::rowvec>(veclen);

  c[0] = -0.5 * 2 * std::sqrt(1.0/3.0);
  c[1] = (1.0-0.5) * 2 * std::sqrt(1.0/3.0);
  c[2] = (1.0-0.5) * 2 * std::sqrt(1.0/3.0);
  
  for (int i = 0; i < n - 1; ++i) {
    arma::rowvec p1 = points.row(i);
    arma::rowvec p2 = points.row(i + 1);
    
    dotprod[i] = arma::dot(arma::cross(p1, p2), c);
  }
  
  return dotprod;
}

// [[Rcpp::export]]
arma::vec scalartripleproduct8(NumericMatrix points) {
  int n = points.nrow();
  int veclen = points.ncol();
  arma::vec dotprod(n - 1);
  arma::rowvec c(veclen, 0.0);
  
  c[0] = (1.0-0.5) * 2 * std::sqrt(1.0/3.0);
  c[1] = (1.0-0.5) * 2 * std::sqrt(1.0/3.0);
  c[2] = -0.5 * 2 * std::sqrt(1.0/3.0);
  
  for (int i = 0; i < n - 1; ++i) {
    arma::rowvec p1 = points.row(i);
    arma::rowvec p2 = points.row(i + 1);
    
    dotprod[i] = arma::dot(arma::cross(p1, p2), c);
  }
  
  return dotprod;
}



// NumericVector vertices_dist(NumericMatrix points, const NumericMatrix &baseCoords) {
//   int n = points.nrow();
//   int m = baseCoords.ncol(); 
//   NumericMatrix distances(n, baseCoords.ncol());
//   int k = 0;
//   
//   if(m == 3){
//     for (int i = 0; i < n; ++i) {
//       for(int j = 0; j < m; ++j){
//         double di = points(i, 0) - baseCoords(j, 0);
//         double dj = points(i, 1) - baseCoords(j, 1);
//         double dk = points(i, 2) - baseCoords(j, 2);
//         distances[k++] = sqrt(di * di + dj * dj + dk * dk);
//       }
//     } 
//   } else if(m == 4) {
//     for (int i = 0; i < n; ++i) {
//       for(int j = 0; j < m; ++j){
//         double dr = points(i, 0) - baseCoords(j, 0);
//         double di = points(i, 1) - baseCoords(j, 1);
//         double dj = points(i, 2) - baseCoords(j, 2);
//         double dk = points(i, 3) - baseCoords(j, 3);
//         distances[k++] = sqrt(dr * dr + di * di + dj * dj + dk * dk);
//       }
//     }
//   } else{
//     for (int i = 0; i < n; ++i) {
//       j = i + 2;
//       double dist = 0;
//       for (int p = 0; p < m; ++p) {
//         double d = points(i, p) - points(j, p);
//         dist += d * d;
//       }
//       distances[k++] = sqrt(dist);
//     }
//   }
//   return distances;
// }








