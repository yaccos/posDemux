library(Biostrings)
input_fasta <- system.file("extdata",
                           "PETRI-seq_forward_reads.fa.gz",
                           package = "posDemux")
reads <- readDNAStringSet(input_fasta)
barcode_files <- system.file("extdata/PETRI-seq_barcodes",
                             c(bc1="bc1.fa",
                               bc2="bc2.fa",
                               bc3="bc3.fa"),
                             package = "posDemux")
names(barcode_files) <- paste0("bc", 1L:3L)
output_table_file <- tempfile(fileext = ".txt")
barcode_index <- map(barcode_files, readDNAStringSet)
barcodes <- barcode_index[c("bc3","bc2","bc1")]
sequence_annotation <- c(UMI = "P", "B", "A", "B", "A", "B", "A")
segment_lengths <- c(7L, 7L, 15L, 7L, 14L, 7L, NA_integer_)
callbacks <- streaming_callbacks(input_file = input_fasta,
                                 output_table_file = output_table_file,
                                 chunk_size = 100L)
streaming_res <- rlang::exec(streaming_demultiplex, !!! callbacks, barcodes=barcodes, allowed_mismatches = 1L,
            segments = sequence_annotation, segment_lengths = segment_lengths)
