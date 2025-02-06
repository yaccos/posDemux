#include <Rcpp.h>
extern "C" {
  // These packages use C linkage, so we need to specify
  // this to avoid name mangling
#include <Biostrings_interface.h>
  // #include <IRanges_interface.h>
  // #include <XVector_interface.h>
  // #include <S4Vectors_interface.h>
}
using namespace Rcpp;

// [[Rcpp::export]]
List hamming_match(SEXP segment, SEXP segment_names, SEXP barcode,
                   CharacterVector barcode_names, int width) {
  XStringSet_holder segment_holder = hold_XStringSet(segment);
  XStringSet_holder barcode_holder = hold_XStringSet(barcode);
  int n_segments = get_length_from_XStringSet_holder(&segment_holder);
  int n_barcodes = get_length_from_XStringSet_holder(&barcode_holder);
  /*
   This array will hold the barcodes sequences in a linear order such that 
   letter j of barcode i is kept in barcode_holder_array[j*n_barcodes + i].
   This ensures optimal cache locality and allows for potential future SIMD
   vectorization
   */
  int n_barcode_letters = n_barcodes * width;
  char *barcode_holder_array = (char*) R_alloc(sizeof(char), n_barcode_letters);
  for (int i = 0; i < n_barcodes; i++) {
    Chars_holder this_barcode = get_elt_from_XStringSet_holder(&barcode_holder, i);
    for (int j = 0; j < width; j++) {
      barcode_holder_array[j*n_barcodes + i] = this_barcode.ptr[j];
    }
  }
  
  IntegerVector this_mismatches(n_barcodes);
  IntegerVector mismatches(n_segments);
  CharacterVector assigned_barcode(n_segments);
  for (int i = 0; i < n_segments; i++) {
    // Initialize the mismatches vector to zero. Otherwise, we would have the
    // iterations accumulate mismatches from previous segments
    this_mismatches.fill(0);
    Chars_holder this_segment = get_elt_from_XStringSet_holder(&segment_holder, i);
    for (int k = 0; k < width; k++) {
      for (int j = 0; j < n_barcodes; j++) {
        /* Increment the number of mismatches if the barcode letter does not match
         This was originally written as an if-statement, but with -O2 optimization 
         (the most common for package repositories),
         this creates a branch which likely causes considerable slowdowns as
         branch prediction will be pure guesswork.
         */
        this_mismatches[j] += this_segment.ptr[k] != barcode_holder_array[k*n_barcodes + j];
      }
      R_xlen_t assigned_barcode_idx = which_min(this_mismatches);
      assigned_barcode[i] = barcode_names[assigned_barcode_idx];
      mismatches[i] = this_mismatches[assigned_barcode_idx];
    }
  }
  
  /*
   * Note: segments_names might be NULL instead of a character vector. However,
   * even NULL values will be handled gracefully by the following assignments  
   */
  assigned_barcode.attr("names") = segment_names;
  mismatches.attr("names") = segment_names;
  return List::create(Named("assigned_barcode") = assigned_barcode,
                      Named("mismatches") = mismatches);
}
