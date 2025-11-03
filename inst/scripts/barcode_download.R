# This scripts downloads the PETRI-seq barcodes from the GitHub repository of the
# PETRI-seq pipeline
library(Biostrings)
library(purrr)
BARCODE_WIDTH <- 7L
destdir <- "inst/extdata/PETRI-seq_barcodes/"
repository_prefix <- "https://raw.githubusercontent.com/tavalab/PETRI-seq-persistence"
commit <- "71032d8511a848f2c96ef3b19bf33fe744b3629f"
barcode_folder <- "PETRI_seq_scripts_v2/scripts/sc_barcodes_v2"

url_paths <- paste(repository_prefix, commit, barcode_folder,
  c(
    bc1 = "BC1_5p_anchor_v2.fa", bc2 = "BC2_anchored.fa",
    bc3 = "BC3_anchored.fa"
  ),
  sep = "/"
)
names(url_paths) <- paste0("bc", 1L:3L)
iwalk(url_paths, \(url, bc_name) readDNAStringSet(url) |>
  # These files also include the linkers between the barcodes,
  # these are being removed
  subseq(start = 1L, width = BARCODE_WIDTH) |>
  writeXStringSet(filepath = paste0(destdir, bc_name, ".fa")))
