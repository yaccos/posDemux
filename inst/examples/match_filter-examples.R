library(purrr)
input_fastq <- system.file("extdata",
                           "PETRI-seq_forward_reads.fq.gz",
                           package = "posDemux")
reads <- readDNAStringSet(input_fastq, format = "fastq")
barcode_files <- system.file("extdata/PETRI-seq_barcodes",
                             c(bc1="bc1.fa",
                               bc2="bc2.fa",
                               bc3="bc3.fa"),
                             package = "posDemux")
names(barcode_files) <- paste0("bc", 1L:3L)
barcode_index <- map(barcode_files, readDNAStringSet)
barcodes <- barcode_index[c("bc3","bc2","bc1")]
sequence_annotation <- c(UMI = "P", "B", "A", "B", "A", "B", "A")
segment_lengths <- c(7L, 7L, 15L, 7L, 14L, 7L, NA_integer_)
demultiplex_res <- posDemux::combinatorial_demultiplex(reads, barcodes = barcodes,
                                                       segments = sequence_annotation,
                                                       segment_lengths = segment_lengths)
filtered_res <- filter_demultiplex_res(demultiplex_res, allowed_mismatches = 1L)
freq_table <- create_frequency_table(filtered_res$demultiplex_res$assigned_barcodes)
print(filtered_res$summary_res)

# This also works, but is usually not necessary to call directly
alternative_summary_res <- create_summary_res(
  retained_sequences = filtered_res$retained,
  barcodes = barcodes,
  assigned_barcodes = demultiplex_res$assigned_barcodes,
  allowed_mismatches = 1L,
  mismatches = demultiplex_res$mismatches
  )
