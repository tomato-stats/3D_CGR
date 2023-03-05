#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix kmerdist(NumericMatrix input_hist) {
  // Input is a histogram with each row associated with an organism 
  // and and each column corresponds to a histogram bin
  int n_organisms = input_hist.nrow();
  NumericVector seq_lengths(n_organisms);
  NumericMatrix Fij(n_organisms, n_organisms);
  NumericMatrix res(n_organisms, n_organisms);
  double a = log(1.1);
  double b = log(0.1) - a;
  
  for (int i = 0; i < n_organisms; ++i) {
    seq_lengths[i] = sum(input_hist(i, _));
  }
  
  for (int i = 0; i < n_organisms; ++i) {
    for (int j = i ; j < n_organisms; ++j) {
      if(i == j){
        Fij(i, j) = 0.5;
        res(i, j) = 0.0;
      } else{
        double denom = fmin(seq_lengths[i], seq_lengths[j]);
        Fij(i, j) = sum(pmin(input_hist(i, _), input_hist(j, _)) / denom);
        res(i, j) = (log(0.1 + Fij(i, j)) - a) / b;
        if(res(i, j) < 0) res(i, j) = 0;
      }
    }
  }
  
  Fij += transpose(Fij);
  res += transpose(res);
  return res;
}

// [[Rcpp::export]]
NumericMatrix kmerdist2(NumericMatrix input_hist) {
  // Input is a histogram with each row associated with an organism 
  // and and each column corresponds to a histogram bin
  int n_organisms = input_hist.nrow();
  NumericVector seq_lengths(n_organisms);
  NumericMatrix Fij(n_organisms, n_organisms);
  NumericMatrix res(n_organisms, n_organisms);
  double a = log(1.1);
  double b = log(0.1) - a;
  
  for (int i = 0; i < n_organisms; ++i) {
    seq_lengths[i] = sum(input_hist(i, _));
  }
  
  for (int i = 0; i < n_organisms; ++i) {
    for (int j = i ; j < n_organisms; ++j) {
      if(i == j){
        Fij(i, j) = 0.5;
        res(i, j) = 0.0;
      } else{
        double denom = fmax(seq_lengths[i], seq_lengths[j]);
        Fij(i, j) = sum(pmin(input_hist(i, _), input_hist(j, _)) / denom);
        res(i, j) = (log(0.1 + Fij(i, j)) - a) / b;
        if(res(i, j) < 0) res(i, j) = 0;
      }
    }
  }
  
  Fij += transpose(Fij);
  res += transpose(res);
  return res;
}


double tancoef(NumericVector x) {
  // Input is a vector of length 2
  if((x(0) == 0) & (x(1) == 0)){
    return(NA_REAL);
  } else {
    return(min(x) / (sum(x) - min(x)));
  }
}

// [[Rcpp::export]]
NumericMatrix tandist(NumericMatrix input_hist) {
  // Input is a histogram with each row associated with an organism 
  // and and each column corresponds to a histogram bin
  int n_organisms = input_hist.nrow();
  int n_cols = input_hist.ncol(); 
  NumericVector seq_lengths(n_organisms);
  NumericVector each_tc(n_cols);
  NumericMatrix TC(n_organisms, n_organisms);
  NumericMatrix dist(n_organisms, n_organisms);
  double a = log(1.1);
  double b = log(0.1) - a;
  
  for (int i = 0; i < n_organisms; ++i) {
    seq_lengths[i] = sum(input_hist(i, _));
  }
  
  for (int i = 0; i < n_organisms; ++i) {
    for (int j = i ; j < n_organisms; ++j) {
      if(i == j){
        TC(i, j) = 0.5;
        dist(i, j) = 0.0;
      } else{
        for(int k = 0; k < n_cols; ++k){
          NumericVector x(2);
          x[0] = input_hist(i, k);
          x[1] = input_hist(j, k);
          each_tc[k] = tancoef(x);
          each_tc = each_tc[!is_na(each_tc)];
          TC(i, j) += mean(each_tc); 
        }
        dist(i, j) = (log(0.1 + TC(i, j)) - a) / b;
        if(dist(i, j) < 0) dist(i, j) = 0;
      }
    }
  }
  TC += transpose(TC);
  dist += transpose(dist);
  return dist;
}