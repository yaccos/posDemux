# Frequency table

Creates a sorted frequency table of each of the observed barcode
combinations. This function is indended to be used after running
[`filter_demultiplex_res()`](https://yaccos.github.io/posDemux/reference/filter_demultiplex_res.md)
and before creating frequency plots, knee plots, or selecting the number
of barcodes to include.

## Usage

``` r
create_freq_table(assigned_barcodes)
```

## Arguments

- assigned_barcodes:

  A character or integer matrix, corresponding to the field
  `assigned_barcodes` from
  [`combinatorial_demultiplex()`](https://yaccos.github.io/posDemux/reference/combinatorial_demultiplex.md)
  or the field `demultiplex_res$assigned_barcodes` from
  [`filter_demultiplex_res()`](https://yaccos.github.io/posDemux/reference/filter_demultiplex_res.md).

## Value

A data frame where each row corresponds to a unique observed barcode
combination. The rows are sorted in descending order of frequency. The
first columns specify the barcode assignment (e.g `bc3`, `bc2`, `bc1`)
and the last columns were the following:

- `frequency`: The number of reads with the barcode combination.

- `cumulative_frequency`: The cumulative frequency of the barcode
  combination counted from the top.

- `fraction`: The fraction of reads with the barcode combination.

- `cumulative_fraction`: The cumulative fraction of the barcode
  combination counted from the top.

## Examples

``` r
library(purrr)
library(Biostrings)
input_fastq <- system.file(
    "extdata", "PETRI-seq_forward_reads.fq.gz", package = "posDemux")
reads <- readDNAStringSet(input_fastq, format = "fastq")
barcode_files <- system.file(
    "extdata/PETRI-seq_barcodes",
    c(bc1 = "bc1.fa", bc2 = "bc2.fa", bc3 = "bc3.fa"),
    package = "posDemux"
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
```
