#include <Rcpp.h>
using namespace Rcpp;

using Rcpp::log;

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


NumericMatrix tancoef(NumericVector x) {
  // Input is a vector of length 2
  if(x(1) == 0 & x(2) == 0){
    return(NA_REAL)
  } else{
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
        for(int k = 0; k < 
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
// [[Rcpp::export]]
tc <- function(input_hist){
  output <- matrix(0, ncol = nrow(input_hist), nrow = nrow(input_hist))
  for(i in 1:(nrow(input_hist) - 1)){
    for(j in (i+1) : nrow(input_hist)){
      output[i, j] <- apply(input_hist[c(i, j),], 2, function(x) min(x)/ (sum(x) - min(x))) |> mean(na.rm = T)
    }
  }
  output <- t(output) + output
    rownames(output) <- rownames(input_hist)
    colnames(output) <- rownames(input_hist)
    a <- log(1.1)
    b <- log(0.1) - a
    res <- (log(0.1 + output) - a) / b
    diag(res) <- 0
  diag(output) <- 1
  return(res)
}