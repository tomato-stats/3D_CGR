#include <Rcpp.h>
#include <algorithm>
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
  colnames(Fij) = rownames(input_hist);
  rownames(Fij) = rownames(input_hist);
  colnames(res) = rownames(input_hist);
  rownames(res) = rownames(input_hist);
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
  colnames(Fij) = rownames(input_hist);
  rownames(Fij) = rownames(input_hist);
  colnames(res) = rownames(input_hist);
  rownames(res) = rownames(input_hist);
  return res;
}

// [[Rcpp::export]]
double tancoef(NumericVector x) {
  // Input is a vector of length 2
  if((x[0] == 0) & (x[1] == 0)){
    return(NA_REAL);
  } else {
    if(x[0] * x[1] > 0){
      return(min(abs(x)) / (sum(abs(x)) - min(abs(x))));
    } else {
      return(0.0);
    }
  }
}


// [[Rcpp::export]]
int count_non_na(Rcpp::NumericVector x) {
  return std::count_if(x.begin(), x.end(), [](double val) { return !Rcpp::NumericVector::is_na(val); });
}

// [[Rcpp::export]]
NumericMatrix tandist(NumericMatrix input_hist) {
  // Input is a histogram with each row associated with an organism 
  // and and each column corresponds to a histogram bin
  int n_organisms = input_hist.nrow();
  int n_cols = input_hist.ncol(); 
  NumericVector seq_lengths(n_organisms);
  NumericVector each_tc(n_cols);
  NumericVector weights(n_cols);
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
      } else {
        each_tc = NumericVector(n_cols);
        each_tc.fill(NA_REAL);  // reset each_tc to all NAs
        double denominator = seq_lengths[i] + seq_lengths[j];
        if(denominator < 0.0001) denominator = 1.0;
        for(int k = 0; k < n_cols; ++k){
          NumericVector x(2);
          x[0] = input_hist(i, k);
          x[1] = input_hist(j, k);
          each_tc[k] = tancoef(x);
          weights[k] = (input_hist(i, k) + input_hist(j, k)) / denominator;
          if(denominator == 1.0) weights[k] = 1.0 / n_cols;
          // if(!NumericVector::is_na(each_tc[k])) Rcout << "K " << k << " " << each_tc[k] << std::endl;
          if(!NumericVector::is_na(each_tc[k])) TC(i, j) += weights[k] * each_tc[k];
        }
        // each_tc = na_omit(each_tc);
        // TC(i, j) = mean(each_tc); 
        dist(i, j) = (log(0.1 + TC(i, j)) - a) / (b);
        if(dist(i, j) < 0) dist(i, j) = 0;
      }
    }
  }
  TC += transpose(TC);
  dist += transpose(dist);
  colnames(TC) = rownames(input_hist);
  rownames(TC) = rownames(input_hist);
  colnames(dist) = rownames(input_hist);
  rownames(dist) = rownames(input_hist);
  return dist;
}

