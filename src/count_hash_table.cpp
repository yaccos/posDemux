#include <Rcpp.h>
#include <unordered_map>
using namespace Rcpp;

typedef std::unordered_map<int, int> count_table_t;

// [[Rcpp::export]]
Rcpp::XPtr<std::unordered_map<int, int>> create_count_table() {
  count_table_t *table = new count_table_t();
  Rcpp::XPtr<count_table_t> ptr(table, true);
  return ptr;
}

// [[Rcpp::export]]
void add_table_entries(Rcpp::XPtr<std::unordered_map<int, int>> table, IntegerVector chunk) {
  for (int entry: chunk) {
    (*table)[entry]++;
  }
}

// [[Rcpp::export]]
List get_count_table(Rcpp::XPtr<std::unordered_map<int, int>> table) {
  int n_entries = table->size();
  IntegerVector encodings(n_entries);
  IntegerVector frequency(n_entries);
  int idx = 0;
  for (auto it = table->begin(); it != table->end(); it++) {
    encodings[idx] = it->first;
    frequency[idx] = it->second;
    idx++;
  }
  List res = List::create(Named("encoding") = encodings, Named("frequency") = frequency);
  return res;
}
