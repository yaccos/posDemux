library(posDemux)
library(magrittr)
library(Biostrings)
set.seed(5203)
width <- 10L
n_barcodes <- length(LETTERS)
barcodes_matrix <- sample(DNA_BASES, size = width * n_barcodes, replace = TRUE) %>%
  matrix(ncol = n_barcodes) %>%
  set_colnames(LETTERS)

barcodes_stringset <- barcodes_matrix %>% 
  apply(2L, FUN = paste0, collapse="") %>% 
  DNAStringSet()

# Up to three mismatches
mismatches <- seq_len(3L)
times <- c(10L,6L,3L)


mutated_barcodes <- posDemux:::mutate_barcodes(barcodes = barcodes_stringset,
                                               mismatches = mismatches,
                                               times = times)

