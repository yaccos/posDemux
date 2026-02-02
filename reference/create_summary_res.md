# Create a summary of match filtering

`create_summary_res()` is a helper function in order to create a summary
of the demultiplexing and following match filtering. It is not designed
to be invoked directly, but its results will be returned automatically
from
[`filter_demultiplex_res()`](https://yaccos.github.io/posDemux/reference/filter_demultiplex_res.md).
This returned object has it own method for printing the result in a
user-friendly manner.

## Usage

``` r
create_summary_res(
  retained,
  barcodes,
  assigned_barcodes,
  allowed_mismatches,
  mismatches
)

# S3 method for class 'demultiplex_filter_summary'
print(x, ...)
```

## Arguments

- retained:

  Logical vector with the same length as the number of reads in the
  input to the demultiplexer. `TRUE` if the corresponding read is
  retained. Corresponds to the field `retained` of the output of
  [`filter_demultiplex_res()`](https://yaccos.github.io/posDemux/reference/filter_demultiplex_res.md).

- barcodes:

  A list of
  [`XStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  objects, the barcodes which were used for demultiplexing.

- assigned_barcodes:

  Character matrix of the assigned barcodes only including the onces
  within the mismatch threshold. Corresponds to of the field
  `demultiplex_res$assigned_barcodes` of
  [`filter_demultiplex_res()`](https://yaccos.github.io/posDemux/reference/filter_demultiplex_res.md).

- allowed_mismatches:

  Integer vector of length one or the same length as the number of
  barcode segments; the threshold Hamming distance. All reads having a
  number of mismatches above this number in any of the barcodes will be
  filtered away.

- mismatches:

  Integer matrix of the number of mismatches of each assigned barcode.
  Corresponds to the field `mismatches` of
  [`combinatorial_demultiplex()`](https://yaccos.github.io/posDemux/reference/combinatorial_demultiplex.md).

- x:

  An object of class `demultiplex_filter_summary` from
  `create_summary_res()`.

- ...:

  Ignored

## Value

`create_summary_res()` returns a list of S3 class
`demultiplex_filter_summary` providing diagnostics for the filtering
process. It contains the the following fields:

- `n_reads`: The total number of reads in the dataset before filtering.

- `n_removed`: The number of reads removed because demultiplexing
  failed.

- `n_barcode_sets`: The number of barcode sets.

- `n_barcode_combinations`: The possible number of barcode combinations.

- `n_unique_barcodes`: The number of observed unique barcode
  combinations (i.e. features which may be cells) detected after
  filtering mismatches.

- `n_estimated_features`: The estimated number of features having a
  detected combination of barcodes. This number will always be greater
  or equal than `n_unique_barcodes` due to barcode collisions.

- `observed_collision_lambda`: The ratio of observed barcode
  combinations divided by the total number of possible barcode
  combinations.

- `corrected_collision_lambda`: The ratio of estimated number of
  features to the total number of possible barcode combinations.

- `expected_collisions`: The statistically expected number of barcode
  collisions or more precicely the expected number of observed barcodes
  which correspond to two or more features.

- `barcode_summary`: A list containing a summary for each barcode set.
  Each element contains the following:

  - `width`: The width (number of nucleotides) of the barcode set.

  - `n_barcodes`: Number of query barcodes.

  - `n_allowed_mismatches`: Number of allowed mismatches for the barcode
    set.

  - `n_removed`: Number of reads having too many mismatches for this
    barcode set.

  - `mismatch_frame`: A `data.frame` with the two columns,
    `n_mismatches` and `frequency` showing the number of reads for each
    of the allowed number of mismatches for the given barcode set.

The [`print()`](https://rdrr.io/r/base/print.html) method returns its
output invisibly.

## Details

Following a uniform distribution of barcodes, the expected number of
barcode collisions (observed barcodes combinations being composed of two
or more features) is given by \$\$N\left(1-e^{-\lambda}-\lambda
e^{-\lambda}\right),\$\$ where \\N\\ is the number of possible barcode
combinations and \\\lambda\\ is in this summary referred to as the
collision lambda: \$\$\lambda=\frac{n}{N},\$\$ where \\n\\ is the number
of features. However, \\n\\ is unknown as we cannot know how many
features there were originally due to potential collisions. Utilizing
the fact that the expected observed number of barcodes is given by
\$\$N\left(1-e^{-\lambda}\right),\$\$ we can correct the estimate for
\\\lambda\\ from the known value of the observed barcode combinations,
and thus estimate the number of features and barcode collisions.

While each unique feature can be conceptually thought of as single cell
with its transcripts, realistic datasets have many features with
relatively small numbers of reads which are artifacts and unlikely to
correspond to true cells.

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
