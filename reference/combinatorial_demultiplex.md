# Combinatorial demultiplexer

This function performs segmenting of sequences and combinatorial
demultiplexing and segmenting of sequences. The sequences have a
structure as defined by the argument `segments` with corresponding
lengths in `segment_lengths`. As segmentation and extraction can take
place from either end, a single middle segment can be variadic in
length. There are three types of segments: Adapter, Barcode and Payload.
The adapter is trimmed and ignored, the barcode is used for
demultiplexing, and the payload is kept after segmenting and
demultiplexing and returned from the function. If there are multiple
payload segments, then each segment constitutes its own segment in a
list. For type stability reasons, such a list is returned also when
there is zero or one payload segments. The barcodes can be positioned at
either end of the sequences, but no barcode can (for obvious reasons) be
variadic in length.

## Usage

``` r
combinatorial_demultiplex(sequences, barcodes, segments, segment_lengths)
```

## Arguments

- sequences:

  A
  [`XStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  object, the sequences to be demultiplexed.

- barcodes:

  A list of
  [`XStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  objects in the same order they appear in the `sequences`, the barcodes
  to be used for demultiplexing. All of the barcodes in each
  [`XStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  must have the same length as specified by the `segment_lengths`
  argument and be named. For computational reasons, the maximum possible
  length of an individual barcode is 127.

- segments:

  Character vector showing the segments of the sequences from 5' end to
  3' end. The code applied is as follows:

  - `'A'`: Adapter (often referred to as linker), is trimmed and ignored

  - `'B'`: Barcode, used for demultiplexing

  - `'P'`: Payload, sequence to be kept after trimming and
    demultiplexing (e.g. cDNA or UMI).

  If this vector is named, this will determine the names of the payload
  sets. Names of the barcode sets will be determined by the names of the
  argument `barcodes` (if any).

- segment_lengths:

  Integer vector with the same length as `segments`, lengths of the
  segments provided in the same order as in `segments`. Up to one of the
  non-barcode segments can have its length set to `NA` which means it is
  considered a variadic length segment.

## Value

A list with the following elements:

- `assigned_barcodes`: A `character` matrix with the names of the
  assigned barcodes as elements. The rows correspond to the sequences
  and the columns to the barcode segments.

- `mismatches`: An `integer` matrix with the number of mismatches
  between the assigned barcodes and the sequences. The rows correspond
  to the sequences and the columns to the barcode segments.

- `payload`: A list of
  [`XStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  objects, each containing the results for a payload segment.

- `barcodes`: The `barcodes` argument passed into the function. It is
  included in order to ease downstream processing.

## Details

If there are two barcodes both having the minimum number of mismatches
the first one will be selected. It is therefore important to choose the
error tolerance to be equal or less than the redundancy of the barcodes.
All sequences are assumed to be long enough for all segments to be
extracted. Otherwise, an error is raised.

## Examples

``` r
library(purrr)
library(Biostrings)
sequence_annotation <- c(UMI = "P", "B", "A", "B", "A", "B", "A")
segment_lengths <- c(7L, 7L, 15L, 7L, 14L, 7L, NA_integer_)
barcode_files <- system.file(
    "extdata/PETRI-seq_barcodes",
    c(bc1 = "bc1.fa", bc2 = "bc2.fa", bc3 = "bc3.fa"),
    package = "posDemux"
    )
names(barcode_files) <- paste0("bc", 1L:3L)
barcode_index <- map(barcode_files, readDNAStringSet)

barcodes <- barcode_index[c("bc3", "bc2", "bc1")]
input_fastq <- system.file(
    "extdata", "PETRI-seq_forward_reads.fq.gz", package = "posDemux")
reads <- readDNAStringSet(input_fastq, format = "fastq")
demultiplex_res <- combinatorial_demultiplex(
    reads, barcodes = barcodes, segments = sequence_annotation,
    segment_lengths = segment_lengths
    )
```
