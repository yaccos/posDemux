# Convert between cutoff types

There are at least two ways to specify the cutoff to use when selecting
barcode combinations (cells) for further analysis. One way is to specify
the number of barcode combinations to keep, effectively keeping a given
number of barcode combinations with the highest frequencies. The other
way is to specify the frequency cutoff directly without regard to the
number of barcode combination to keep. In the former case,
`bc_to_freq_cutoff()` is used to find the corresponding frequency
cutoff, whereas in the latter case `freq_to_bc_cutoff()` is used to find
the corresponding barcode cutoff.

## Usage

``` r
bc_to_freq_cutoff(freq_table, cutoff)

freq_to_bc_cutoff(freq_table, cutoff)
```

## Arguments

- freq_table:

  The frequency table from
  [`create_freq_table()`](https://yaccos.github.io/posDemux/reference/create_freq_table.md).
  In case the table is derived from another source, it must be sorted in
  descending order of frequency.

- cutoff:

  Integer vector, the cutoff values to be converted.

## Value

Integer, the converted cutoff values.

## Details

In the edge case of the barcode threshold being zero, the frequency
cutoff is set to the maximum frequency in the table plus one. This
feature makes sure that the cutoff line is visible in the frequency
plot.

## Examples

``` r
library(purrr)
library(Biostrings)
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> Loading required package: stats4
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> 
#> Attaching package: ‘IRanges’
#> The following object is masked from ‘package:purrr’:
#> 
#>     reduce
#> Loading required package: XVector
#> 
#> Attaching package: ‘XVector’
#> The following object is masked from ‘package:purrr’:
#> 
#>     compact
#> Loading required package: Seqinfo
#> 
#> Attaching package: ‘Biostrings’
#> The following object is masked from ‘package:base’:
#> 
#>     strsplit
input_fastq <- system.file(
    "extdata", "PETRI-seq_forward_reads.fq.gz",
    package = "posDemux")
reads <- readDNAStringSet(input_fastq, format = "fastq")
barcode_files <- system.file(
    "extdata/PETRI-seq_barcodes", c(bc1 = "bc1.fa", bc2 = "bc2.fa",
    bc3 = "bc3.fa"), package = "posDemux")
names(barcode_files) <- paste0("bc", 1L:3L)
barcode_index <- map(barcode_files, readDNAStringSet)
barcodes <- barcode_index[c("bc3", "bc2", "bc1")]
sequence_annotation <- c(UMI = "P", "B", "A", "B", "A", "B", "A")
segment_lengths <- c(7L, 7L, 15L, 7L, 14L, 7L, NA_integer_)
demultiplex_res <- posDemux::combinatorial_demultiplex(
    reads, barcodes = barcodes, segments = sequence_annotation,
    segment_lengths = segment_lengths)
filtered_res <- filter_demultiplex_res(demultiplex_res, allowed_mismatches = 1L)
freq_table <- create_freq_table(filtered_res$demultiplex_res$assigned_barcodes)

bc_cutoff <- c(100L, 204L, 50L, 655L)
freq_cutoff <- bc_to_freq_cutoff(freq_table, bc_cutoff)
# Note: The reconstructed barcode cutoff is not equal to
# the original due to ties in the frequency table
reconstructed_bc_cutoff <- freq_to_bc_cutoff(freq_table, freq_cutoff)
# The frequency cutoff is still preserved through these conversions
reconstruced_freq_cutoff <- bc_to_freq_cutoff(
    freq_table, reconstructed_bc_cutoff)
```
