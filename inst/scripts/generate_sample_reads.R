# This script generates synthetic example PETRI-seq reads which are used for the
# vignette and the examples
set.seed(59252L)
library(Biostrings)
library(purrr)
library(tibble)
library(dplyr)
library(magrittr)
sample_DNA_sequences <- function(n, length) {
    vapply(seq_len(n), \(n) sample(Biostrings::DNA_BASES, length, replace = TRUE) |>
        paste0(collapse = ""), FUN.VALUE = character(1L))
}
n_cells <- 500L
n_barcodes_artifacts <- 500L
n_junk_sequences <- 5000L
mean_umis_per_cell_barcode <- 10L
mean_reads_per_artifact_barcode <- 3L
mean_reads_per_umi <- 10L
cell_rbinom_size <- 10L
artifact_rbinom_size <- 20L
bc3_to_bc2_linker <- "GGTCCTTGGCTTCGC"
bc2_to_bc1_linker <- "CCTCCTACGCCAGA"
segment_map <- c(UMI = "P", bc3 = "B", l_1 = "A", bc2 = "B", l_2 = "A", bc1 = "B")
segment_lengths <- segment_lengths <- c(7L, 7L, 15L, 7L, 14L, 7L)
seq_length <- sum(segment_lengths)


generate_junk_sequences <- function(n_junk_sequences, seq_length) {
    sample_DNA_sequences(n_junk_sequences, seq_length) |>
        DNAStringSet()
}

barcode_files <- system.file("extdata/PETRI-seq_barcodes", c(bc1 = "bc1.fa", bc2 = "bc2.fa",
    bc3 = "bc3.fa"), package = "posDemux")
names(barcode_files) <- paste0("bc", 1L:3L)
barcode_index <- map(barcode_files, readDNAStringSet)

concatenate_and_shuffle <- function(barcode_frame) {
    barcode_frame |>
        tidyr::unnest(cols = umi) %$%
        xscat(umi, barcode_index[["bc3"]][bc3], bc3_to_bc2_linker, barcode_index[["bc2"]][bc2],
            bc2_to_bc1_linker, barcode_index[["bc1"]][bc1])
}

generate_cell_sequences <- function(n_cells) {
    cell_barcode_frame <- map(barcode_index, \(barcode) sample(barcode |>
        names(), n_cells, replace = TRUE)) |>
        as_tibble()
    cell_barcode_frame$n_umis <- rnbinom(n_cells, mu = mean_umis_per_cell_barcode, size = cell_rbinom_size) |>
        as.integer()
    cell_barcode_frame$umi <- map(cell_barcode_frame$n_umis, \(n_umis) {
        umi <- sample_DNA_sequences(n_umis, segment_lengths[[1L]])
        reads_per_umi <- rnbinom(n_umis, mu = mean_reads_per_umi, size = cell_rbinom_size) |>
            as.integer()
        rep(umi, times = reads_per_umi)
    })
    cell_barcode_frame |>
        select(!n_umis) |>
        concatenate_and_shuffle()
}

generate_artifact_sequences <- function(n_artifacts) {
    artifact_barcode_frame <- map(barcode_index, \(barcode) sample(barcode |>
        names(), n_artifacts, replace = TRUE)) |>
        as_tibble()
    artifact_barcode_frame$n_reads <- rnbinom(n_artifacts, mu = mean_reads_per_artifact_barcode,
        size = artifact_rbinom_size) |>
        as.integer()
    artifact_barcode_frame$umi <- map(artifact_barcode_frame$n_reads, \(n_reads) {
            umi <- sample_DNA_sequences(n_reads, segment_lengths[[1L]])
            umi
        })
    artifact_barcode_frame |>
        select(!n_reads) |>
        concatenate_and_shuffle()
}

# Sequences without barcodes
junk_sequences <- generate_junk_sequences(n_junk_sequences, seq_length)
# The sequences we want to keep
cell_sequences <- generate_cell_sequences(n_cells)
# Sequences which have barcodes attaches to them, but do not have desireable properties
artifact_sequences <- generate_artifact_sequences(n_barcodes_artifacts)

combined_sequences <- c(junk_sequences, cell_sequences, artifact_sequences) |>
    sample()

names(combined_sequences) <- paste0("seq_", seq_along(combined_sequences))

# In order to make compliant, albeit not realistic PhreadQualities
Q_scores <- sample(15:40, size = seq_length * length(combined_sequences), replace = TRUE)
Phread_encoding <- character() %>%
    PhredQuality() %>%
    encoding()

seq_quality <- match(Q_scores, Phread_encoding) %>%
    {
        extract(names(Phread_encoding), .)
    } %>%
    matrix(ncol = length(combined_sequences)) %>%
    apply(2L, paste, collapse = "") %>%
    PhredQuality()

quality_combined_sequences <- QualityScaledDNAStringSet(combined_sequences, seq_quality)

writeQualityScaledXStringSet(quality_combined_sequences, "inst/extdata/PETRI-seq_forward_reads.fq.gz",
    compress = TRUE)
