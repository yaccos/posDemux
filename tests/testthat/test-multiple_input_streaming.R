library(Biostrings)
library(purrr)
input_fastq <- system.file("extdata", "PETRI-seq_forward_reads.fq.gz", package = "posDemux")
reads <- readQualityScaledDNAStringSet(input_fastq)
n_reads <- length(reads)
# Splits the reads into "odd" and "even" reads taking every second read
read_is_even <- seq_len(n_reads) %% 2 == 0
even_reads <- reads[read_is_even]
odd_reads <- reads[!read_is_even]
# Writes the split reads back to two temporary files
even_fastq <- tempfile()
odd_fastq <- tempfile()
writeQualityScaledXStringSet(even_reads, even_fastq)
writeQualityScaledXStringSet(odd_reads, odd_fastq)
barcode_files <- system.file("extdata/PETRI-seq_barcodes", c(
    bc1 = "bc1.fa", bc2 = "bc2.fa",
    bc3 = "bc3.fa"
), package = "posDemux")
names(barcode_files) <- paste0("bc", 1L:3L)
output_table_file <- tempfile(fileext = ".txt")
barcode_index <- map(barcode_files, readDNAStringSet)
barcodes <- barcode_index[c("bc3", "bc2", "bc1")]
sequence_annotation <- c(UMI = "P", "B", "A", "B", "A", "B", "A")
segment_lengths <- c(7L, 7L, 15L, 7L, 14L, 7L, NA_integer_)
callbacks <- streaming_callbacks(
    input_file = c(even_fastq, odd_fastq), output_table_file = output_table_file,
    chunk_size = 10000L, verbose = FALSE, min_width = sum(segment_lengths, na.rm = TRUE)
)
streaming_res <- rlang::exec(streaming_demultiplex, !!!callbacks,
                             barcodes = barcodes, allowed_mismatches = 1L,
                             segments = sequence_annotation, segment_lengths = segment_lengths
)
expected_filtered_res <- combinatorial_demultiplex(
    sequences = c(even_reads, odd_reads), barcodes = barcodes,
    segments = sequence_annotation, segment_lengths = segment_lengths
) %>%
    filter_demultiplex_res(allowed_mismatches = 1L)
expected_summary <- expected_filtered_res$summary_res
expected_assigned_barcodes <- expected_filtered_res$demultiplex_res$assigned_barcodes
expected_UMI <- expected_filtered_res$demultiplex_res$payload$UMI
expected_freq_table <- create_freq_table(expected_assigned_barcodes)

test_that("Demultiplex summary from streaming with multiple input files is correctly generated", {
    expect_equal(streaming_res$summary_res, expected_summary)
})

test_that("Frequency table is correctly generated from streaming with multiple input files", {
    test_freq_table(streaming_res$freq_table, expected_freq_table, names(barcodes))
})

test_that("Barcode table is correctly generated from streaming with multiple input files", {
    barcode_table <- read.table(output_table_file, header = TRUE)
    expected_barcode_table <- data.frame(read = rownames(expected_assigned_barcodes), UMI = expected_UMI %>%
                                             unname()) %>%
        cbind(as.data.frame(expected_assigned_barcodes, row.names = FALSE))
    expect_equal(barcode_table, expected_barcode_table)
})
