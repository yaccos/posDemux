#' Row matching of tables
#'
#' @description
#' This function extends the functionality of \code{\link[base]{%in%}} for
#' finding which rows in the first argument is exists in the second.
#' 
#' 
#' @param x A matrix or data frame which rows to be matched. 
#' Typically, this will be a matrix of assigned barcodes for each read.
#' @param table A matrix or data frame with the rows to matched against. Typically,
#' this will be the top portion of a frequency table.
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
  mapping <- get_mapping(table_df)
  table_encoded <- encode(table_df, mapping)
  x_encoded <- encode(x_barcodes, mapping)
  x_encoded %in% table_encoded
}

get_mapping <- function(table){
  # The order the barcodes appear in does not really matter
  # as long as the same mapping object is used
  unique_values <- map(table, unique)
  n_unique_values <- map(unique_values, length)
  cumprod <- cumprod(n_unique_values)
  shifted_cumprod <- c(1L, cumprod(n_unique_values) %>% head(-1L))
  list(unique = unique_values, n_unique = n_unique_values, 
       cumprod = cumprod,
       shifted_cumprod =  shifted_cumprod)
}

encode <- function(x, mapping) {
   x %>% 
    map2(mapping$unique, \(x, table) match(x, table) - 1L) %>%
    map2(mapping$shifted_cumprod, `*`) %>% 
    {Reduce(`+`, ., simplify = TRUE)}
}

decode <- function(x_encoded, mapping) {
    {map2(mapping$cumprod, mapping$shifted_cumprod,
        \(modulus, dividend)  (x_encoded %% modulus) %/% dividend + 1L)} %>% 
    {map2(mapping$unique, ., `[`)}
}
