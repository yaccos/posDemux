# Demultiplexing with streaming

This function provides an interface to
[`combinatorial_demultiplex()`](https://yaccos.github.io/posDemux/reference/combinatorial_demultiplex.md)
and
[`filter_demultiplex_res()`](https://yaccos.github.io/posDemux/reference/filter_demultiplex_res.md)
such that reads are streamed in chunks instead having to load everything
at once, hence reducing memory consumption. It accepts two functions
which are called once per chunk: A data loader function for producing
the sequences of the chunk and an archiver writing the results to file.

## Usage

``` r
streaming_demultiplex(
  state_init,
  loader,
  archiver,
  barcodes,
  allowed_mismatches,
  segments,
  segment_lengths
)
```

## Arguments

- state_init:

  The initial state to pass into `loader`.

- loader:

  Function loading the reads. It has the signature `f(state)`, where
  `state` is a user-defined object which is initialized to be
  `state_init` and for the subsequent iterations taken as the `state`
  field of the output of `archiver`. Its return value is a list with the
  following fields:

  - `state`: The state to be passed into `archiver`.

  - `sequences`: A
    [`XStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
    object, the sequences to be demultiplexed in the current chunk.

  - `should_terminate`: A scalar logical. If `TRUE`, the demultiplexing
    process terminates and the final results are returned. Notice that
    this termination happens before the sequences of the final call to
    `loader` are demultiplexed.

- archiver:

  Function taking care of archiving the demultiplexed results. Its
  arguments are:

  - `state`: The state of the process returned by `loader`.

  - `filtered_res`: The output from running
    [`combinatorial_demultiplex()`](https://yaccos.github.io/posDemux/reference/combinatorial_demultiplex.md)
    and
    [`filter_demultiplex_res()`](https://yaccos.github.io/posDemux/reference/filter_demultiplex_res.md)
    on the data expect that the field `summary_res` is missing.

  Its output is a state object fed into the next call to `loader`.

- barcodes:

  A list of
  [`XStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  objects in the same order they appear in the `sequences`, the barcodes
  to be used for demultiplexing. All of the barcodes in each
  [`XStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  must have the same length as specified by the `segment_lengths`
  argument and be named. For computational reasons, the maximum possible
  length of an individual barcode is 127.

- allowed_mismatches:

  Integer vector of length one or the same length as the number of
  barcode segments; the threshold Hamming distance. All reads having a
  number of mismatches above this number in any of the barcodes will be
  filtered away.

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

A list with three elements:

- `freq_table`: The frequency table for all reads, akin to the output of
  [`create_freq_table()`](https://yaccos.github.io/posDemux/reference/create_freq_table.md).

- `summary_res`: The summary result of match filtering of all reads per
  [`create_summary_res()`](https://yaccos.github.io/posDemux/reference/create_summary_res.md).

- `state_final`: The final state object returned from `loader`.

## Details

The data loader decides the size of each chunk. While this framework
does not provide any restriction on the `state` object, the loader and
archiver must be written such that the state objects they return are
compatible. Since the data loader alone decides when to terminate, bad
terminations crieria can cause a runaway loop. Usually, it will be
useful to have a progress tracker of how many reads are demultiplexed.
The framework itself does not implement this, so it is typically
implemented into the archiver or loader.

For technical reasons, it is not possible to do streaming when the
number of possible barcode combinations exceeds \\2^{32}-1\approx
2.1\cdot 10^{9}\\.

## See also

[`filter_demultiplex_res()`](https://yaccos.github.io/posDemux/reference/filter_demultiplex_res.md),
[`combinatorial_demultiplex()`](https://yaccos.github.io/posDemux/reference/combinatorial_demultiplex.md),
[`create_freq_table()`](https://yaccos.github.io/posDemux/reference/create_freq_table.md),
and
[`create_summary_res()`](https://yaccos.github.io/posDemux/reference/create_summary_res.md)
for the underlying processing.

## Examples

``` r
library(purrr)
library(Biostrings)
input_fastq <- system.file(
    "extdata", "PETRI-seq_forward_reads.fq.gz", package = "posDemux")
output_barcode_table <- tempfile(pattern = "barcode_table", fileext = ".txt")

callbacks <- streaming_callbacks(
    input_file = input_fastq, output_table_file = output_barcode_table,
    chunk_size = 10000, verbose = TRUE)
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
streaming_summary_res <- streaming_demultiplex(
    state_init = callbacks$state_init, loader = callbacks$loader,
    archiver = callbacks$archiver, barcodes = barcodes, allowed_mismatches = 1L,
    segments = sequence_annotation, segment_lengths = segment_lengths
    )
#> Tue Dec 16 10:16:19 2025 => Initializing FASTQ stream and output table
#> Tue Dec 16 10:16:19 2025 => Streaming FASTQ input file /home/runner/work/_temp/Library/posDemux/extdata/PETRI-seq_forward_reads.fq.gz
#> Tue Dec 16 10:16:19 2025 => Processed 10000 reads, successfully demultiplexed 9117 reads so far...
#> Tue Dec 16 10:16:19 2025 => Processed 20000 reads, successfully demultiplexed 18259 reads so far...
#> Tue Dec 16 10:16:19 2025 => Processed 30000 reads, successfully demultiplexed 27362 reads so far...
#> Tue Dec 16 10:16:19 2025 => Processed 40000 reads, successfully demultiplexed 36487 reads so far...
#> Tue Dec 16 10:16:19 2025 => Processed 50000 reads, successfully demultiplexed 45629 reads so far...
#> Tue Dec 16 10:16:19 2025 => Processed 56895 reads, successfully demultiplexed 51906 reads so far...
#> Tue Dec 16 10:16:19 2025 => Done demultiplexing
```
