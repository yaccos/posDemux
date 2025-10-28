library(Biostrings)
# This script is used for the streaming vignette where we illustrate how we can
# use a custom callback for parsing FASTA files. Therefore, we must first convert
# the FASTQ file into FASTA

input_fastq <- system.file("extdata",
                           "PETRI-seq_forward_reads.fq.gz",
                           package = "posDemux")

output_fasta <- "inst/extdata/PETRI-seq_forward_reads.fa.gz"

stringset <- readDNAStringSet(filepath = input_fastq,format = "fastq")

writeXStringSet(stringset, filepath = output_fasta, compress = TRUE,
                format = "fasta")
