#' @name posDemux
#' 
#' @title Utilities for segmented sequences with combinatorial barcoding
#' @useDynLib posDemux, .registration = TRUE
#' 
#' 
#' @docType package
#' 
#' @description
#' Provides tools for handling reads with combinatorial barcodes and
#' multiple adapter regions. This includes utilities for demultiplexing,
#' sequence trimming and filtering. The package is intended to work with
#' single-cell RNA-seq data with multiple barcoding
#' (e.g. SPLiT-seq and PETRI-seq)
#' 
#' 
#' @import assertthat
#' 
#' @author Jakob Peder Pettersen <jakobpeder.pettersen@gmail.com>
#' 
"_PACKAGE"