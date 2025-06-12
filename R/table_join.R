#' Row matching of tables
#'
#' @description
#' This function extends the functionality of \code{\link[base]{%in%}} for
#' finding which rows in the first argument is found in the second.
#' 
#' 
#' @param x A matrix or data frame which rows to be matched. 
#' Typically, this will be a matrix of assigned barcodes for each read.
#' @param table A matrix or data frame with the rows to matched against. Typically,
#' this will be the top portion of a frequency table.
#' character matrix of assigned barcodes with barcodes as
#' columns and reads as rows.
#' @seealso [create_frequency_table()] for how frequency tables are constructed,
#' [combinatorial_demultiplex()] for more information on the matrix of assigned
#'  barcodes
#'  and [dplyr::inner_join()] for a function with similar functionality.
#'  
#' @details
#' As this function is intended to be used for data frames containing more than
#' just the barcodes, the intersection of the column names is used for matching.
#' As opposed to [base::match()], this function is implemented more efficiently
#' by converting each row into a numeric encoding before matching.
#' 
#'
#' @returns Logical vector, for each row in \code{table}, is the same row found
#' row in \code{x}.
#' @export
#'
row_match <- function(x, table) {
  barcode_cols <- intersect(colnames(x), colnames(table))
  table_df <- as.data.frame(table[, barcode_cols])
  x_barcodes <- as.data.frame(x[,barcode_cols])
  n_barcodes <- ncol(table)
  # The order the barcodes appear in does not really matter, just that it stays
  # the same during execution of the function
  table_unique_barcodes <- map(table_df, unique)
  n_unique_barcodes <- map_int(table_unique_barcodes, length)
  cumulative_mapping <- c(1L, cumprod(n_unique_barcodes) %>% head(-1L))
  shifted_cumulative_mapping <- cumprod(n_unique_barcodes)
  encode <- . %>%
    # We need zero-based indexing for this to work
    map2(table_unique_barcodes, \(x, table) match(x, table) - 1L) %>%
    map2(cumulative_mapping, `*`) %>% 
    {Reduce(`+`, ., simplify = TRUE)}
  decode <- . %>% 
    {map2(shifted_cumulative_mapping, cumulative_mapping,
         \(x,y)  (. %% x) %/% y + 1L)} %>% 
    map2(table_unique_barcodes, `[`)
    
  table_encoded <- encode(table_df)
  x_encoded <- encode(x_barcodes)
  x_encoded %in% table
}

encode <- function(x, table_unique_values){
  
}

