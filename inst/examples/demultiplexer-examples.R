library(purrr)
sequence_annotation <- c(UMI = "P", "B", "A", "B", "A", "B", "A")
segment_lengths <- c(7L, 7L, 15L, 7L, 14L, 7L, NA_integer_)
barcode_files <- system.file("extdata/PETRI-seq_barcodes",
c(bc1="bc1.fa",
 bc2="bc2.fa",
 bc3="bc3.fa"),
package = "posDemux")
names(barcode_files) <- paste0("bc", 1L:3L)
barcode_index <- map(barcode_files, readDNAStringSet)

barcodes <- barcode_index[c("bc3","bc2","bc1")]
input_fastq <- system.file("extdata",
                           "PETRI-seq_forward_reads.fq.gz",
                           package = "posDemux")
reads <- readDNAStringSet(input_fastq, format = "fastq")
demultiplex_res <- combinatorial_demultiplex(
  reads,
  barcodes = barcodes,
  segments = sequence_annotation,
  segment_lengths = segment_lengths
)

