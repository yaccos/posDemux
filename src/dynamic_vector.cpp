#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

typedef std::vector<int> ivec;

// [[Rcpp::export]]
Rcpp::XPtr<std::vector<int> > create_encoding_vector() {
  ivec *vec = new ivec();
  Rcpp::XPtr<ivec> ptr(vec, true);
  return ptr;
}

// [[Rcpp::export]]
void grow_encoding_vector(Rcpp::XPtr<std::vector<int>> vec, IntegerVector chunk) {
  vec->insert(vec->end(), chunk.begin(), chunk.end());
}

// [[Rcpp::export]]
IntegerVector get_encoding_vector(Rcpp::XPtr<std::vector<int>> vec) {
  return Rcpp::wrap(*vec);
}
