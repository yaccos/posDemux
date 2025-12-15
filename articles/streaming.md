# Streaming data into the demultiplexer

## Why streaming?

Given that sequence files can consist of gigabytes of data, storing
everything in memory at once can prove to be too much for the computer.
Fortunately, for demultiplexing of barcodes, we can demultiplex the
barcodes independently of each other. In theory, we could read the
entries from a FASTQ file one by one, and for each of them call the
demultiplexer and output the barcode assignment to a file. However, due
to overhead in file operations and the `R` interpreter, the most
efficient option is to process the files in chunks of a predetermined
size. If you have not read the main package vignette yet
([`vignette("demultiplexing")`](https://yaccos.github.io/posDemux/articles/demultiplexing.md)),
please do so before proceeding.

### Which chunk size to use?

The chunk size should be large enough for the constant overhead of
initiating a file read and interpreting the required `R` instructions to
be neglectable. Still, the chunk size should be small enough to only
make a modest memory footprint. In practice, with modern computers
having gigabytes of RAM, we find a chunk size of $2 \cdot 10^{5}$ reads
to be a reasonable choice. The example datasets used in this vignette
are far too small to require streaming (including a dataset of realistic
size would make the package exceed the maximum package size), so for
purely demonstrational purposes, we choose a chunk size of $10^{4}$
reads.

### Overview of the streaming API

The streaming API of `posDemux` is focused around two user-supplied
functions:

- The loader: This function is responsible for grabbing a chunk of the
  input file of the sequences and returning a `XStringSet` object of the
  sequences. Also, this function determines when to terminate the
  streaming, typically done when there are no sequences left.
- The archiver: This function is responsible for writing the results of
  the demultiplexer of a chunk to file, typically in the form of a
  barcode table.

Once the loader and archiver are specified and combined with the
demultiplexer, they can generate the barcode table chunk by chunk.
However, we still want to generate the summary statistics and frequency
table as if all the reads were demultiplexed simultaneously. If the
chunks are demultiplexed independently,
[`create_summary_res()`](https://yaccos.github.io/posDemux/reference/create_summary_res.md)
and
[`create_freq_table()`](https://yaccos.github.io/posDemux/reference/create_freq_table.md)
don’t work as intended. That is where the core function of the streaming
API;
[`streaming_demultiplex()`](https://yaccos.github.io/posDemux/reference/streaming_demultiplex.md)
comes into play. It accepts the loader and the archiver as callback
functions, communicates with the demultiplexer, and does the bookkeeping
to construct the summary statistics and frequency table. Hence, for each
chunk, the loader provides the sequences to demultiplex, the streaming
logic feeds these reads into the demultiplexer, updates is record of
observed barcodes and provides the results to the archiver which writes
to the barcode table. This process is illustrated in Figure
@ref(fig:structure).

For communicating between the loader and archiver, a user-defined state
object is passed through
[`streaming_demultiplex()`](https://yaccos.github.io/posDemux/reference/streaming_demultiplex.md).
The initial value of this object is provided to
[`streaming_demultiplex()`](https://yaccos.github.io/posDemux/reference/streaming_demultiplex.md)
as the argument `state_init`.

![Components of the streaming API](streaming_structure.svg)

Components of the streaming API

## Streaming using default callbacks

Although the option of specifying the loader and archiver gives great
flexibility, it does so at the cost of convenience as the loader and
archiver can be challenging to write. For this reason, `posDemux`
provides a default set of streaming callbacks available through
[`streaming_callbacks()`](https://yaccos.github.io/posDemux/reference/streaming_callbacks.md).
Actually, this is a function factory, which returns a list of the
following elements:

- `state_init`: The initial state to provide to
  [`streaming_demultiplex()`](https://yaccos.github.io/posDemux/reference/streaming_demultiplex.md)
- `loader`: Loader function
- `archiver`: Archiver function

### Assumptions

The callbacks are designed to be useful for a typical use case where:

- The reads are contained as a single FASTQ file on disk.
- The barcode table contains columns with the read names (`read`),
  sequences of all payloads (eg.`UMI`) and the barcode assignments (eg.
  `BC3`, `BC2`, `BC1`).
- Read names with spaces are truncated by the loader such that only the
  part of the read name before the space is kept. This is done because
  the information preceding the space is typically enough to ambiguously
  identify the read, while the information after the space usually
  contains irrelevant auxiliary information.
- Progress of how many reads have been processed can optionally be
  displayed by the archiver for each chunk.

If this does not fill your needs, please refer to the section [Streaming
using custom callbacks](#streaming-using-custom-callbacks).

### Obtaining the callbacks

As said,
[`streaming_callbacks()`](https://yaccos.github.io/posDemux/reference/streaming_callbacks.md)
is a function factory, we must specify the file path to the FASTQ input
file, and the barcode table output file for it to generate useful
callbacks. In addition, the chunk size is regulated through the
callbacks. We use the same dataset as in the main package vignette
([`vignette("demultiplexing")`](https://yaccos.github.io/posDemux/articles/demultiplexing.md))
and obtain the callbacks as follows:

``` r
library(posDemux)
input_fastq <- system.file("extdata",
    "PETRI-seq_forward_reads.fq.gz",
    package = "posDemux"
)
output_barcode_table <- tempfile(
    pattern = "barcode_table",
    fileext = ".txt"
)
default_callbacks <- streaming_callbacks(
    input_file = input_fastq,
    output_table_file = output_barcode_table,
    chunk_size = 1e+4,
    verbose = TRUE
)
```

### Running the demultiplexing

With the callbacks in place, we specify the barcodes, and the segments
and feed them into the streaming framework:

``` r
suppressPackageStartupMessages({
    library(purrr)
    library(Biostrings)
})
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
    state_init = default_callbacks$state_init,
    loader = default_callbacks$loader,
    archiver = default_callbacks$archiver,
    barcodes = barcodes,
    allowed_mismatches = 1L,
    segments = sequence_annotation,
    segment_lengths = segment_lengths
)
#> Mon Dec 15 17:15:59 2025 => Initializing FASTQ stream and output table
#> Mon Dec 15 17:15:59 2025 => Streaming FASTQ input file /home/runner/work/_temp/Library/posDemux/extdata/PETRI-seq_forward_reads.fq.gz
#> Mon Dec 15 17:16:05 2025 => Processed 10000 reads, successfully demultiplexed 9117 reads so far...
#> Mon Dec 15 17:16:05 2025 => Processed 20000 reads, successfully demultiplexed 18259 reads so far...
#> Mon Dec 15 17:16:05 2025 => Processed 30000 reads, successfully demultiplexed 27362 reads so far...
#> Mon Dec 15 17:16:05 2025 => Processed 40000 reads, successfully demultiplexed 36487 reads so far...
#> Mon Dec 15 17:16:05 2025 => Processed 50000 reads, successfully demultiplexed 45629 reads so far...
#> Mon Dec 15 17:16:05 2025 => Processed 56895 reads, successfully demultiplexed 51906 reads so far...
#> Mon Dec 15 17:16:05 2025 => Done demultiplexing
```

### Displaying the result from streaming

From the return value of
[`streaming_demultiplex()`](https://yaccos.github.io/posDemux/reference/streaming_demultiplex.md),
we obtain the frequency table

``` r
head(streaming_summary_res$freq_table)
#>      bc3    bc2    bc1 frequency cumulative_frequency    fraction
#> 1 bc3_37  bc2_4 bc1_30       270                  270 0.005201711
#> 2 bc3_45 bc2_36 bc1_65       246                  516 0.004739336
#> 3 bc3_25 bc2_34 bc1_25       235                  751 0.004527415
#> 4 bc3_40 bc2_62 bc1_42       221                  972 0.004257697
#> 5 bc3_69 bc2_95 bc1_37       219                 1191 0.004219165
#> 6 bc3_56 bc2_94 bc1_29       217                 1408 0.004180634
#>   cumulative_fraction
#> 1         0.005201711
#> 2         0.009941047
#> 3         0.014468462
#> 4         0.018726159
#> 5         0.022945324
#> 6         0.027125958
```

and summary printout:

``` r
streaming_summary_res$summary_res
#> Total number of reads: 56895 
#> Number of reads failing to demultiplex: 4989 (8.77%) 
#> Observed number of unique barcode combinations:978 
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
```

## Streaming the barcode table

### Considerations about size of tables

The frequency table from a realistic PETRI-seq dataset is relatively
small (just a few megabytes) and should therefore not pose any
substantial memory usage problems. On the other hand, the barcode table
will have one row for each read processed, making its memory footprint
closer to the size of the original FASTQ-file. To this end, the default
archiver writes to the barcode table on disk for each chunk, thus
avoiding keeping the table in memory. However, the barcode table will be
required by some downstream step of the pipeline. Typically, this will
be frequency filtering. The good news is that the cutoff determination
process only needs the light frequency table. The bad news is that the
full barcode table is required when determining which reads to keep
after the cutoff is selected.

If we want to avoid heavy load on RAM, we must stream the barcode table
in chunks, determine for each chunk which reads to keep, and appending
the results to another file. \## Subsetting the frequency table

We begin this example by subsetting the frequency table. In order not to
repeat ourselves, we refer to the main package vignette
([`vignette("demultiplexing")`](https://yaccos.github.io/posDemux/articles/demultiplexing.md))
for the considerations behind this process.

``` r
freq_table <- streaming_summary_res$freq_table
bc_cutoff <- 500L
selected_freq_table <- freq_table[seq_len(bc_cutoff), ]
```

### Writing to database

We typically want to create a copy of the barcode table including only
those reads which are selected by the barcode frequency and thus match a
row in `selected_freq_table`. We **could** use the most direct approach
and provide the selected barcode table in the same file format, just
with desired rows kept. However, this could become impractical for
downstream analysis as using this new table as a lookup would require us
to load it all into memory, causing the same problem one level down. For
memory-constrained environments, we therefore suggest another idea:
Dumping the data into a SQLite database which later can be queried. For
simplicity, we paste all the barcode designations together in one
column. We begin by setting up the database:

``` r
suppressPackageStartupMessages({
    library(DBI)
})
output_database_path <- tempfile(
    pattern = "selected_barcode_table",
    fileext = ".sqlite"
)
output_db <- dbConnect(RSQLite::SQLite(), output_database_path)
```

For doing the streaming of the barcode table and combining it with
outputting to the database, we use the `chunked` package:

``` r
chunk_size <- 1e4
suppressPackageStartupMessages({
    library(chunked)
    library(magrittr)
})
read_table_chunkwise(
    file = output_barcode_table,
    header = TRUE,
    sep = "\t",
    chunk_size = chunk_size
) %>%
    # For use within a pipeline, it is more convenient to use inner_join() 
    # than posDemux::row_match() and it achieves the save result
    inner_join(selected_freq_table) %>%
    mutate(
        celltag = do.call(
            paste,
            c(barcodes %>% names() %>% all_of() %>% across(), sep = "_")
        )
    ) %>%
    select(read, UMI, celltag) %>%
    write_chunkwise(dbplyr::src_dbi(output_db), "selected_barcodes")
#> Joining with `by = join_by(bc3, bc2, bc1)`
#> Joining with `by = join_by(bc3, bc2, bc1)`
#> Joining with `by = join_by(bc3, bc2, bc1)`
#> Joining with `by = join_by(bc3, bc2, bc1)`
#> Joining with `by = join_by(bc3, bc2, bc1)`
#> Joining with `by = join_by(bc3, bc2, bc1)`
#> # Source:   table<`selected_barcodes`> [?? x 3]
#> # Database: sqlite 3.51.1 [/tmp/Rtmp0p0lbR/selected_barcode_tablea893401f95a9.sqlite]
#>    read   UMI     celltag             
#>    <chr>  <chr>   <chr>               
#>  1 seq_1  GCCTAAC bc3_57_bc2_51_bc1_94
#>  2 seq_2  CCAAGCG bc3_72_bc2_94_bc1_95
#>  3 seq_4  CCTAACG bc3_85_bc2_81_bc1_17
#>  4 seq_5  GCTCGTC bc3_49_bc2_64_bc1_77
#>  5 seq_6  TGGAGAA bc3_3_bc2_45_bc1_33 
#>  6 seq_7  ACTTCGA bc3_17_bc2_57_bc1_71
#>  7 seq_8  GCCCGCA bc3_69_bc2_65_bc1_49
#>  8 seq_9  CCCATGT bc3_38_bc2_87_bc1_83
#>  9 seq_10 TGCCGGT bc3_76_bc2_36_bc1_85
#> 10 seq_11 CAACTCA bc3_15_bc2_84_bc1_16
#> # ℹ more rows
```

Finally, we make a index on the barcodes such that lookups are fast and
close the database connection:

``` r
dbExecute(output_db, "CREATE INDEX idx ON selected_barcodes(read)") %>%
    invisible()
dbDisconnect(output_db)
```

## Streaming using custom callbacks

### Desired properties

In cases where the default callbacks are not sufficient, it is possible
for the user to write custom callback functions. For our example, we
will use the same sample dataset as before, but pretend the reads for
some reason come in a FASTA file (which the default loader does not
support). In addition, we want the UMIs to go into its own FASTA files
instead of being a column in the barcode table. Unlike the default
callbacks, we don’t truncate read names.

### The state object

For progress printouts, we want the state object of the streaming to
contain the total number of reads and the number of reads where
demultiplexing was successful. Additionally, for our specific case, the
loader needs to know where to begin reading the FASTA file for each
chunk. For the final numbers, this is not really necessary because the
required information is provided by the streaming framework itself when
finishing. Thus, we want the state object to contain the following
fields which are both integers:

- `total_reads`
- `demultiplexed_reads`

`total_reads` will be updated by the loader, whereas
`demultiplexed_reads` will be updated by the archiver at the end of each
chunk. Furthermore, the archiver needs to write the column headers to
the barcode table when first writing to it, so it must know whether it
is working at its first chunk. This gives rise to the logical flag
`output_table_initialized`. We thus initialize the state as follows:

``` r
custom_state_init <- list(
    total_reads = 0L, demultiplexed_reads = 0L,
    output_table_initialized = FALSE
)
```

### The loader

The loader is a function which accepts the state as its sole parameter
and outputs a list with three fields:

- `state`: To be passed into the archiver.
- `sequences`: The sequences in the chunk as a
  [`Biostrings::DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  object.
- `should_terminate`: A logical flag to signal to the streaming
  framework that it should terminate. Note that this also means the
  loader is in charge of the termination criterion.

We decide to hard-code the loader to parse $10^{4}$ reads for each
chunk.

``` r
input_fasta <- system.file("extdata", "PETRI-seq_forward_reads.fa.gz",
    package = "posDemux"
)

loader <- function(state) {
    n_reads_per_chunk <- 1e4
    if (!state$output_table_initialized) {
        message("Starting to load sequences")
    }
    chunk <- Biostrings::readDNAStringSet(
        filepath = input_fasta,
        format = "fasta",
        nrec = n_reads_per_chunk,
        skip = state$total_reads
    )
    n_reads_in_chunk <- length(chunk)
    if (n_reads_in_chunk == 0L && state$output_table_initialized) {
        # The case when the initial chunk is empty is given
        # special treatment since
        # we want to create the table regardless
        # No more reads to process, make the outer framework return
        final_res <- list(
            state = state,
            sequences = NULL,
            should_terminate = TRUE
        )
        message("Done demultiplexing")
        return(final_res)
    }
    state$total_reads <- state$total_reads + n_reads_in_chunk
    list(
        state = state,
        sequences = chunk,
        should_terminate = FALSE
    )
}
```

### The archiver

The archiver has the state obtained from the loader as its first
argument and the results from
[`filter_demultiplex_res()`](https://yaccos.github.io/posDemux/reference/filter_demultiplex_res.md)
of the chunk as its second argument. Its return value is the state to
provide to the loader when it continues with its next chunk.

``` r
output_umi_file <- tempfile(pattern = "umi.fa", fileext = ".txt")
custom_output_table <- tempfile(
    pattern = "barcode_table_custom",
    fileext = ".txt"
)
archiver <- function(state, filtered_res) {
    barcode_matrix <- filtered_res$demultiplex_res$assigned_barcodes
    barcode_names <- colnames(barcode_matrix)
    read_names <- rownames(barcode_matrix)
    # If the table has no rows, we may risk getting a NULL value
    if (is.null(read_names)) {
        read_names <- character()
    }
    barcode_table <- as.data.frame(barcode_matrix)
    read_name_table <- data.frame(read = read_names)
    umis <- filtered_res$demultiplex_res$payload$UMI

    chunk_table <- cbind(read_name_table, barcode_table)
    if (!state$output_table_initialized) {
        append <- FALSE
        state$output_table_initialized <- TRUE
    } else {
        append <- TRUE
    }

    Biostrings::writeXStringSet(
        umis,
        filepath = output_umi_file, append = TRUE
    )
    readr::write_tsv(
        x = chunk_table,
        file = custom_output_table,
        append = append,
        col_names = !append,
        eol = "\n"
    )
    state <- within(state, {
        demultiplexed_reads <- demultiplexed_reads + nrow(barcode_matrix)
        paste0(
            "Processed {total_reads} reads,", " ",
            "successfully demultiplexed {demultiplexed_reads} reads so far..."
        ) %>%
            glue::glue() %>%
            message()
    })
    state
}
```

### Putting it all together

With the initial state and the callbacks being defined, we are now ready
to put them to use:

``` r
custom_summary_res <- streaming_demultiplex(
    state_init = custom_state_init,
    loader = loader,
    archiver = archiver,
    barcodes = barcodes,
    allowed_mismatches = 1L,
    segments = sequence_annotation,
    segment_lengths = segment_lengths
)
#> Starting to load sequences
#> Processed 10000 reads, successfully demultiplexed 9117 reads so far...
#> Processed 20000 reads, successfully demultiplexed 18259 reads so far...
#> Processed 30000 reads, successfully demultiplexed 27362 reads so far...
#> Processed 40000 reads, successfully demultiplexed 36487 reads so far...
#> Processed 50000 reads, successfully demultiplexed 45629 reads so far...
#> Processed 56895 reads, successfully demultiplexed 51906 reads so far...
#> Done demultiplexing
```

### Verifying results

We finally load the files we wrote to in order to verify that they
contain what we desire. We load the barcode table:

``` r
readr::read_tsv(custom_output_table, show_col_types = FALSE)
#> # A tibble: 51,906 × 4
#>    read   bc3    bc2    bc1   
#>    <chr>  <chr>  <chr>  <chr> 
#>  1 seq_1  bc3_57 bc2_51 bc1_94
#>  2 seq_2  bc3_72 bc2_94 bc1_95
#>  3 seq_4  bc3_85 bc2_81 bc1_17
#>  4 seq_5  bc3_49 bc2_64 bc1_77
#>  5 seq_6  bc3_3  bc2_45 bc1_33
#>  6 seq_7  bc3_17 bc2_57 bc1_71
#>  7 seq_8  bc3_69 bc2_65 bc1_49
#>  8 seq_9  bc3_38 bc2_87 bc1_83
#>  9 seq_10 bc3_76 bc2_36 bc1_85
#> 10 seq_11 bc3_15 bc2_84 bc1_16
#> # ℹ 51,896 more rows
```

and the UMI sequences:

``` r
Biostrings::readDNAStringSet(output_umi_file, format = "fasta")
#> DNAStringSet object of length 51906:
#>         width seq                                           names               
#>     [1]     7 GCCTAAC                                       seq_1
#>     [2]     7 CCAAGCG                                       seq_2
#>     [3]     7 CCTAACG                                       seq_4
#>     [4]     7 GCTCGTC                                       seq_5
#>     [5]     7 TGGAGAA                                       seq_6
#>     ...   ... ...
#> [51902]     7 ATGGATC                                       seq_56890
#> [51903]     7 ACGGCTT                                       seq_56891
#> [51904]     7 ATGAGCG                                       seq_56892
#> [51905]     7 ACCTGCG                                       seq_56893
#> [51906]     7 GCCTTGA                                       seq_56895
```

## Reproducibility

The *[posDemux](https://bioconductor.org/packages/3.23/posDemux)*
package was made possible thanks to:

- R
- *[BiocStyle](https://bioconductor.org/packages/3.23/BiocStyle)*
- *[knitr](https://CRAN.R-project.org/package=knitr)*
- *[RefManageR](https://CRAN.R-project.org/package=RefManageR)*
- *[rmarkdown](https://CRAN.R-project.org/package=rmarkdown)*
- *[sessioninfo](https://CRAN.R-project.org/package=sessioninfo)*
- *[testthat](https://CRAN.R-project.org/package=testthat)*

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.23/biocthis)*.

`R` session information:

    #> ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
    #>  setting  value
    #>  version  R Under development (unstable) (2025-12-12 r89163)
    #>  os       Ubuntu 24.04.3 LTS
    #>  system   x86_64, linux-gnu
    #>  ui       X11
    #>  language en
    #>  collate  C.UTF-8
    #>  ctype    C.UTF-8
    #>  tz       UTC
    #>  date     2025-12-15
    #>  pandoc   3.1.11 @ /opt/hostedtoolcache/pandoc/3.1.11/x64/ (via rmarkdown)
    #>  quarto   NA
    #> 
    #> ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
    #>  package              * version   date (UTC) lib source
    #>  abind                  1.4-8     2024-09-12 [1] RSPM
    #>  assertthat             0.2.1     2019-03-21 [1] RSPM
    #>  backports              1.5.0     2024-05-23 [1] RSPM
    #>  bibtex                 0.5.1     2023-01-26 [1] RSPM
    #>  Biobase                2.71.0    2025-10-30 [1] Bioconduc~
    #>  BiocGenerics         * 0.57.0    2025-10-30 [1] Bioconduc~
    #>  BiocManager            1.30.27   2025-11-14 [1] RSPM
    #>  BiocParallel           1.45.0    2025-10-30 [1] Bioconduc~
    #>  BiocStyle            * 2.39.0    2025-10-30 [1] Bioconduc~
    #>  Biostrings           * 2.79.2    2025-11-05 [1] Bioconduc~
    #>  bit                    4.6.0     2025-03-06 [1] RSPM
    #>  bit64                  4.6.0-1   2025-01-16 [1] RSPM
    #>  bitops                 1.0-9     2024-10-03 [1] RSPM
    #>  blob                   1.2.4     2023-03-17 [1] RSPM
    #>  bookdown               0.46      2025-12-05 [1] RSPM
    #>  bslib                  0.9.0     2025-01-30 [1] RSPM
    #>  cachem                 1.1.0     2024-05-16 [1] RSPM
    #>  chunked              * 0.6.1     2025-10-02 [1] RSPM
    #>  cigarillo              1.1.0     2025-10-31 [1] Bioconduc~
    #>  cli                    3.6.5     2025-04-23 [1] RSPM
    #>  codetools              0.2-20    2024-03-31 [3] CRAN (R 4.6.0)
    #>  crayon                 1.5.3     2024-06-20 [1] RSPM
    #>  DBI                  * 1.2.3     2024-06-02 [1] RSPM
    #>  dbplyr                 2.5.1     2025-09-10 [1] RSPM
    #>  DelayedArray           0.37.0    2025-10-31 [1] Bioconduc~
    #>  deldir                 2.0-4     2024-02-28 [1] RSPM
    #>  desc                   1.4.3     2023-12-10 [1] RSPM
    #>  digest                 0.6.39    2025-11-19 [1] RSPM
    #>  dplyr                * 1.1.4     2023-11-17 [1] RSPM
    #>  evaluate               1.0.5     2025-08-27 [1] RSPM
    #>  farver                 2.1.2     2024-05-13 [1] RSPM
    #>  fastmap                1.2.0     2024-05-15 [1] RSPM
    #>  fs                     1.6.6     2025-04-12 [1] RSPM
    #>  generics             * 0.1.4     2025-05-09 [1] RSPM
    #>  GenomicAlignments      1.47.0    2025-10-31 [1] Bioconduc~
    #>  GenomicRanges          1.63.1    2025-12-08 [1] Bioconduc~
    #>  ggplot2                4.0.1     2025-11-14 [1] RSPM
    #>  glue                   1.8.0     2024-09-30 [1] RSPM
    #>  gtable                 0.3.6     2024-10-25 [1] RSPM
    #>  hms                    1.1.4     2025-10-17 [1] RSPM
    #>  htmltools              0.5.9     2025-12-04 [1] RSPM
    #>  htmlwidgets            1.6.4     2023-12-06 [1] RSPM
    #>  httpuv                 1.6.16    2025-04-16 [1] RSPM
    #>  httr                   1.4.7     2023-08-15 [1] RSPM
    #>  hwriter                1.3.2.1   2022-04-08 [1] RSPM
    #>  interp                 1.1-6     2024-01-26 [1] RSPM
    #>  IRanges              * 2.45.0    2025-10-31 [1] Bioconduc~
    #>  jpeg                   0.1-11    2025-03-21 [1] RSPM
    #>  jquerylib              0.1.4     2021-04-26 [1] RSPM
    #>  jsonlite               2.0.0     2025-03-27 [1] RSPM
    #>  knitr                  1.50      2025-03-16 [1] RSPM
    #>  LaF                    0.8.6     2024-12-13 [1] RSPM
    #>  later                  1.4.4     2025-08-27 [1] RSPM
    #>  lattice                0.22-7    2025-04-02 [3] CRAN (R 4.6.0)
    #>  latticeExtra           0.6-31    2025-09-10 [1] RSPM
    #>  lifecycle              1.0.4     2023-11-07 [1] RSPM
    #>  lubridate              1.9.4     2024-12-08 [1] RSPM
    #>  magrittr             * 2.0.4     2025-09-12 [1] RSPM
    #>  Matrix                 1.7-4     2025-08-28 [3] CRAN (R 4.6.0)
    #>  MatrixGenerics         1.23.0    2025-10-30 [1] Bioconduc~
    #>  matrixStats            1.5.0     2025-01-07 [1] RSPM
    #>  memoise                2.0.1     2021-11-26 [1] RSPM
    #>  mime                   0.13      2025-03-17 [1] RSPM
    #>  otel                   0.2.0     2025-08-29 [1] RSPM
    #>  pillar                 1.11.1    2025-09-17 [1] RSPM
    #>  pkgconfig              2.0.3     2019-09-22 [1] RSPM
    #>  pkgdown                2.2.0     2025-11-06 [1] RSPM
    #>  plyr                   1.8.9     2023-10-02 [1] RSPM
    #>  png                    0.1-8     2022-11-29 [1] RSPM
    #>  posDemux             * 0.99.8    2025-12-15 [1] local
    #>  promises               1.5.0     2025-11-01 [1] RSPM
    #>  purrr                * 1.2.0     2025-11-04 [1] RSPM
    #>  pwalign                1.7.0     2025-10-31 [1] Bioconduc~
    #>  R6                     2.6.1     2025-02-15 [1] RSPM
    #>  ragg                   1.5.0     2025-09-02 [1] RSPM
    #>  RColorBrewer           1.1-3     2022-04-03 [1] RSPM
    #>  Rcpp                   1.1.0.8.1 2025-12-08 [1] RSPM
    #>  readr                  2.1.6     2025-11-14 [1] RSPM
    #>  RefManageR           * 1.4.0     2022-09-30 [1] RSPM
    #>  rlang                  1.1.6     2025-04-11 [1] RSPM
    #>  rmarkdown              2.30      2025-09-28 [1] RSPM
    #>  Rsamtools              2.27.0    2025-10-31 [1] Bioconduc~
    #>  RSQLite                2.4.5     2025-11-30 [1] RSPM
    #>  S4Arrays               1.11.1    2025-11-25 [1] Bioconduc~
    #>  S4Vectors            * 0.49.0    2025-10-30 [1] Bioconduc~
    #>  S7                     0.2.1     2025-11-14 [1] RSPM
    #>  sass                   0.4.10    2025-04-11 [1] RSPM
    #>  scales                 1.4.0     2025-04-24 [1] RSPM
    #>  Seqinfo              * 1.1.0     2025-10-31 [1] Bioconduc~
    #>  sessioninfo          * 1.2.3     2025-02-05 [1] RSPM
    #>  shiny                  1.12.1    2025-12-09 [1] RSPM
    #>  ShortRead              1.69.2    2025-11-05 [1] Bioconduc~
    #>  SparseArray            1.11.9    2025-12-10 [1] Bioconduc~
    #>  stringi                1.8.7     2025-03-27 [1] RSPM
    #>  stringr                1.6.0     2025-11-04 [1] RSPM
    #>  SummarizedExperiment   1.41.0    2025-10-31 [1] Bioconduc~
    #>  systemfonts            1.3.1     2025-10-01 [1] RSPM
    #>  textshaping            1.0.4     2025-10-10 [1] RSPM
    #>  tibble                 3.3.0     2025-06-08 [1] RSPM
    #>  tidyselect             1.2.1     2024-03-11 [1] RSPM
    #>  timechange             0.3.0     2024-01-18 [1] RSPM
    #>  tzdb                   0.5.0     2025-03-15 [1] RSPM
    #>  utf8                   1.2.6     2025-06-08 [1] RSPM
    #>  vctrs                  0.6.5     2023-12-01 [1] RSPM
    #>  vroom                  1.6.7     2025-11-28 [1] RSPM
    #>  withr                  3.0.2     2024-10-28 [1] RSPM
    #>  xfun                   0.54      2025-10-30 [1] RSPM
    #>  xml2                   1.5.1     2025-12-01 [1] RSPM
    #>  xtable                 1.8-4     2019-04-21 [1] RSPM
    #>  XVector              * 0.51.0    2025-10-31 [1] Bioconduc~
    #>  yaml                   2.3.12    2025-12-10 [1] RSPM
    #> 
    #>  [1] /home/runner/work/_temp/Library
    #>  [2] /opt/R/devel/lib/R/site-library
    #>  [3] /opt/R/devel/lib/R/library
    #>  * ── Packages attached to the search path.
    #> 
    #> ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

## Bibliography

This vignette was generated using
*[BiocStyle](https://bioconductor.org/packages/3.23/BiocStyle)* with
*[knitr](https://CRAN.R-project.org/package=knitr)* and
*[rmarkdown](https://CRAN.R-project.org/package=rmarkdown)* running
behind the scenes.
