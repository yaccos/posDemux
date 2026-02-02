# Filter demultiplexed reads

Filters the demultiplexed reads from
[`combinatorial_demultiplex()`](https://yaccos.github.io/posDemux/reference/combinatorial_demultiplex.md)
such that any read exceeding the number of allowed mismatches for any of
the barcodes is removed. The function gives diagnostic information on
the number of reads removed per barcode and the total number of reads
removed.

## Usage

``` r
filter_demultiplex_res(demultiplex_res, allowed_mismatches)
```

## Arguments

- demultiplex_res:

  Unprocessed output from
  [`combinatorial_demultiplex()`](https://yaccos.github.io/posDemux/reference/combinatorial_demultiplex.md).

- allowed_mismatches:

  Integer vector of length one or the same length as the number of
  barcode segments; the threshold Hamming distance. All reads having a
  number of mismatches above this number in any of the barcodes will be
  filtered away.

## Value

A list with the following elements:

- `demultiplex_res`: The contents of the input argument
  `demultiplex_res` with the sequences filtered.

- `retained`: Logical vector with the same length as the number of reads
  in the input. `TRUE` if the corresponding read is retained. Useful for
  future filtering of paired-end reads.

- `summary_res`: Result of
  [`create_summary_res()`](https://yaccos.github.io/posDemux/reference/create_summary_res.md)
  called on the results of filtering.

## Details

The value of `n_removed` does not in general equal the sum of
`n_removed_per_barcode` since a read can have too many mismatches with
multiple barcodes.

## See also

[`create_summary_res()`](https://yaccos.github.io/posDemux/reference/create_summary_res.md)

## Examples

``` r
library(purrr)
library(Biostrings)
input_fastq <- system.file(
    "extdata", "PETRI-seq_forward_reads.fq.gz", package = "posDemux")
reads <- readDNAStringSet(input_fastq, format = "fastq")
barcode_files <- system.file(
    "extdata/PETRI-seq_barcodes",
    c(bc1 = "bc1.fa", bc2 = "bc2.fa", bc3 = "bc3.fa"), package = "posDemux"
    )
names(barcode_files) <- paste0("bc", 1L:3L)
barcode_index <- map(barcode_files, readDNAStringSet)
barcodes <- barcode_index[c("bc3", "bc2", "bc1")]
sequence_annotation <- c(UMI = "P", "B", "A", "B", "A", "B", "A")
segment_lengths <- c(7L, 7L, 15L, 7L, 14L, 7L, NA_integer_)
demultiplex_res <- posDemux::combinatorial_demultiplex(
    reads, barcodes = barcodes, segments = sequence_annotation,
    segment_lengths = segment_lengths
    )
filtered_res <- filter_demultiplex_res(demultiplex_res, allowed_mismatches = 1L)
freq_table <- create_freq_table(filtered_res$demultiplex_res$assigned_barcodes)
print(filtered_res$summary_res)
#> Total number of reads: 56895 
#> Number of reads failing to demultiplex: 4989 (8.77%) 
#> Observed number of unique barcode combinations: 978 
#> Number of possible barcode combinations: 884736 
#> Estimated number of features: 978.5 
#> Observed feature to barcode ratio: 0.001105 
#> Corrected feature to barcode ratio: 0.001106 
#> Estimated number of observed barcode combinations
#> corresponding to more than one feature: 0.5 (0.06%) 
#> Number of barcode sets: 3 
#> - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#> Barcode set: bc3 
#> Barcode width: 7 
#> Number of possible barcodes: 96 
#> Number of allowed mismatches: 1 
#> Number of reads with 0 mismatches: 51923 (91.26%) 
#> Number of reads with 1 mismatches: 671 (1.18%) 
#> Number of reads above mismatch threshold: 4301 (7.56%) 
#> - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#> Barcode set: bc2 
#> Barcode width: 7 
#> Number of possible barcodes: 96 
#> Number of allowed mismatches: 1 
#> Number of reads with 0 mismatches: 51918 (91.25%) 
#> Number of reads with 1 mismatches: 567 (1%) 
#> Number of reads above mismatch threshold: 4410 (7.75%) 
#> - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#> Barcode set: bc1 
#> Barcode width: 7 
#> Number of possible barcodes: 96 
#> Number of allowed mismatches: 1 
#> Number of reads with 0 mismatches: 51930 (91.27%) 
#> Number of reads with 1 mismatches: 642 (1.13%) 
#> Number of reads above mismatch threshold: 4323 (7.6%) 
#> - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# This also works, but is usually not necessary to call directly
alternative_summary_res <- create_summary_res(
    retained = filtered_res$retained, barcodes = barcodes,
    assigned_barcodes = filtered_res$demultiplex_res$assigned_barcodes,
    allowed_mismatches = 1L, mismatches = demultiplex_res$mismatches
    )
```
