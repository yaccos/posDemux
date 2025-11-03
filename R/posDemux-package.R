#' @name posDemux
#' 
#' @keywords internal
#' @aliases posDemux-package
#' @title Utilities for segmented sequences with combinatorial barcoding
#' @useDynLib posDemux, .registration = TRUE
#' 
#' @description
#' Provides tools for handling reads with combinatorial barcodes, and
#' multiple adapter regions. This includes utilities for demultiplexing,
#' sequence segmenting and filtering. The package is intended to work with
#' single-cell RNA-seq data with multiple barcoding
#' (e.g. SPLiT-seq and PETRI-seq).
#' 
#' 
#' @import assertthat
#' @importFrom magrittr %>% set_names %<>% extract equals subtract add
#' 
#' @author Jakob Peder Pettersen <jakobpeder.pettersen@gmail.com>
#' 
# Making sure R CMD CHECK does not complain over the magrittr placeholder
utils::globalVariables(c("."))
"_PACKAGE"
