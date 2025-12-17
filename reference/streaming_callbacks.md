# Suggested setup for FASTQ streaming

Even though the user can define the arguments `state_init`, `loader`,
and `archiver` for
[`streaming_demultiplex()`](https://yaccos.github.io/posDemux/reference/streaming_demultiplex.md),
this approach is only recommended for advanced users. This functions
defines a premade combinations of these three arguments which should be
suitable in most cases. The loader streams a FASTQ file (having multiple
files is also supported) and the archiver outputs a data frame to file
consisting of the read name (`read`), the sequences of all payloads
(e.g. `UMI`), and barcode assignments (`c('bc3','bc2','bc1')`).

## Usage

``` r
streaming_callbacks(
  input_file,
  output_table_file,
  chunk_size = 1e+06,
  verbose = TRUE,
  min_width = NULL
)
```

## Arguments

- input_file:

  A character vector containing the paths to the FASTQ files to be used
  for demultiplexing. Often this is only one file, but multiple files
  are supported such that demultiplexing data from multiple lanes does
  not require merging the lanes first.

- output_table_file:

  The path to which the output barcode table will be written.

- chunk_size:

  Integer, the number of reads to process in each chunk.

- verbose:

  Logical scalar: Should the progress be displayed?

- min_width:

  Optional integer scalar: Minimum width of the sequences to keep. For
  reads which are shorter than this, a warning is emitted and the reads
  are removed and ignored and thus not appear in any statistics. The
  data loader is **not** supposed to be used as a length filter, so this
  option is more like an escape hatch for being able to deal with
  sequences which have not been properly filtered beforehand.

## Value

A list with the following elements, all of which are intended to be used
as the corresponding arguments to
[`streaming_demultiplex()`](https://yaccos.github.io/posDemux/reference/streaming_demultiplex.md):

- `state_init`

- `loader`

- `archiver`

## Details

If the read names have any spaces in them, the loader will only keep the
portion of the read name preceding the first space. This is due to the
Illumina platform's behavior of encoding the sequencing direction
(forward or reverse) past the space. Keeping the read names with the
space is usually not desirable as it makes the resulting barcode table
more confusing and makes it more difficult to group the forward and
reverse reads together afterwards.

## See also

[`ShortRead::FastqStreamer()`](https://rdrr.io/pkg/ShortRead/man/Sampler-class.html)
which is used as a backend for the data loader.

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
#> Wed Dec 17 10:22:50 2025 => Initializing FASTQ stream and output table
#> Wed Dec 17 10:22:50 2025 => Streaming FASTQ input file /home/runner/work/_temp/Library/posDemux/extdata/PETRI-seq_forward_reads.fq.gz
#> Wed Dec 17 10:22:55 2025 => Processed 10000 reads, successfully demultiplexed 9117 reads so far...
#> Wed Dec 17 10:22:55 2025 => Processed 20000 reads, successfully demultiplexed 18259 reads so far...
#> Wed Dec 17 10:22:55 2025 => Processed 30000 reads, successfully demultiplexed 27362 reads so far...
#> Wed Dec 17 10:22:55 2025 => Processed 40000 reads, successfully demultiplexed 36487 reads so far...
#> Wed Dec 17 10:22:55 2025 => Processed 50000 reads, successfully demultiplexed 45629 reads so far...
#> Wed Dec 17 10:22:55 2025 => Processed 56895 reads, successfully demultiplexed 51906 reads so far...
#> Wed Dec 17 10:22:55 2025 => Done demultiplexing
```
