mutation_p <- list(
  # Read as: p=0.1 for 1 mismatch, p=0.04 for 2 mismatches,
  # p=0.008 for 3 mismatches (above threshold)
  bc1=c(0.1,0.04,0.008),
  bc2=c(0.2,0.05,0.02,0.005),
  bc3=c(0.05,0.01)
)

# Three barcode sets
segment_map <- c(a_1="A",bc_1="B",a_2="A",p_1="P",
                 bc_2="B",p_2="P",a_4="A",bc_3="B",p_3="P")
barcode_n_allowed_mismatches <- c(bc_1=2L, bc_2=3L, bc_3=1L)
n_segments <- length(segment_map)
n_unique_barcodes <- 1000L
mean_reads_per_cell <- 100L
mean_reads_per_artifact <- 5L 
n_barcode_sets <- (segment_map == "B") %>% sum()
segment_length <- c(20L, 9L, 13L, 9L, 12L, NA_integer_, 18L, 5L, 14L)
rbinom_size <- 10L
