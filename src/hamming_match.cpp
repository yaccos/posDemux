#include <Rcpp.h>
// [[Rcpp::depends(Biostrings)]]
// [[Rcpp::depends(IRanges)]]
// [[Rcpp::depends(XVector)]]
// [[Rcpp::depends(S4Vectors)]]
extern "C" {
  // These packages use C linkage, so we need to specify
  // this to avoid name mangling
#include <Biostrings_interface.h>
#include <IRanges_interface.h>
#include <XVector_interface.h>
#include <S4Vectors_defines.h>
}
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
List hamming_match(SEXP segment, CharacterVector segment_names, SEXP barcode, CharacterVector barcode_names, int width) {
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
    Chars_holder this_segment = get_elt_from_XStringSet_holder(&segment_holder, i);
    for (int k = 0; k < width; k++) {
      for (int j = 0; j < n_barcodes; j++) {
        if (this_segment.ptr[k] != barcode_holder_array[k*n_barcodes + j]) {
          this_mismatches[j]++;
        }
      }
    }
    R_xlen_t assigned_barcode_idx = which_min(this_mismatches);
    assigned_barcode[i] = barcode_names[assigned_barcode_idx];
    mismatches[i] = this_mismatches[assigned_barcode_idx];
  }
  
  assigned_barcode.attr("names") = segment_names;
  mismatches.attr("names") = segment_names;
  return List::create(Named("assigned_barcode") = assigned_barcode,
               Named("mismatches") = mismatches);
  // 
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
library(Biostrings)
barcodes <- DNAStringSet(c(a="ATCG", b="GGGA", c="CCTA"))
segments <- DNAStringSet(c(a="ATCG", b="GCGA", c="CCAT"))
hamming_match(segments, names(segments), barcodes, names(barcodes), 4L)
*/
