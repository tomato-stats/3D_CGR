#include <Rcpp.h>
using namespace Rcpp;

// Function to calculate the next coordinates
NumericVector nextCoords(const char &base, const NumericVector &prevCoords,
                         const NumericMatrix &baseCoords) {
  int baseIndex;
  switch (base) {
  case 'A':
    baseIndex = 0;
    break;
  case 'C':
    baseIndex = 1;
    break;
  case 'G':
    baseIndex = 2;
    break;
  case 'T':
    baseIndex = 3;
    break;
  case 'R':
    baseIndex = 4;
    break;
  case 'Y':
    baseIndex = 5;
    break;
  case 'K':
    baseIndex = 6;
    break;
  case 'M':
    baseIndex = 7;
    break;
  case 'S':
    baseIndex = 8;
    break;
  case 'W':
    baseIndex = 9;
    break;
  case 'B':
    baseIndex = 10;
    break;
  case 'D':
    baseIndex = 11;
    break;
  case 'H':
    baseIndex = 12;
    break;
  case 'V':
    baseIndex = 13;
    break;
  case 'N':
    baseIndex = 14;
    break;
  default:
    baseIndex = -1;
  }
  
  NumericVector baseCoord = baseCoords(baseIndex, _);
  NumericVector newCoords(prevCoords.size());
  for (int i = 0; i < prevCoords.size(); i++) {
    newCoords[i] = (prevCoords[i] + baseCoord[i]) / 2;
  }
  return newCoords;
}

// [[Rcpp::export]]
DataFrame chaosGameDF(const CharacterVector &bases, const NumericMatrix &baseCoords) {
  int n = bases.size();
  NumericMatrix coords(n + 1, baseCoords.ncol());
  
  // Add a row of zeros at the top of the matrix
  for (int j = 0; j < coords.ncol(); j++) {
    coords(0, j) = 0.0;
  }
  
  NumericVector prevCoords(baseCoords.ncol());
  for (int i = 0; i < prevCoords.size(); i++) {
    prevCoords[i] = 0;
  }
  
  for (int i = 0; i < n; i++) {
    std::string base = Rcpp::as<std::string>(bases[i]);
    coords(i + 1, _) = nextCoords(base[0], prevCoords, baseCoords);
    prevCoords = coords(i + 1, _);
  }
  
  return DataFrame::create(_["i"] = coords(_, 0),
                           _["j"] = coords(_, 1),
                           _["k"] = coords(_, 2));

}

// [[Rcpp::export]]
NumericMatrix chaosGame(const CharacterVector &bases, const NumericMatrix &baseCoords) {
  int n = bases.size();
  NumericMatrix coords(n + 1, baseCoords.ncol());
  
  // Add a row of zeros at the top of the matrix
  for (int j = 0; j < coords.ncol(); j++) {
    coords(0, j) = 0.0;
  }
  
  NumericVector prevCoords(baseCoords.ncol());
  for (int i = 0; i < prevCoords.size(); i++) {
    prevCoords[i] = 0;
  }
  
  for (int i = 0; i < n; i++) {
    std::string base = Rcpp::as<std::string>(bases[i]);
    coords(i + 1, _) = nextCoords(base[0], prevCoords, baseCoords);
    prevCoords = coords(i + 1, _);
  }
  
  return coords;
}