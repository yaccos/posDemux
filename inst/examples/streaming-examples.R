library(purrr)
library(Biostrings)
input_fastq <- system.file("extdata",
  "PETRI-seq_forward_reads.fq.gz",
  package = "posDemux"
)
output_barcode_table <- tempfile(
  pattern = "barcode_table",
  fileext = ".txt"
)

callbacks <- streaming_callbacks(
  input_file = input_fastq,
  output_table_file = output_barcode_table,
  chunk_size = 1e+4,
  verbose = TRUE
)
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
  state_init = callbacks$state_init,
  loader = callbacks$loader,
  archiver = callbacks$archiver,
  barcodes = barcodes,
  allowed_mismatches = 1L,
  segments = sequence_annotation,
  segment_lengths = segment_lengths
)
