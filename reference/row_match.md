# Row matching of tables

This function extends the functionality of
[`%in%`](https://rdrr.io/r/base/match.html) for finding which rows in
the first argument exist in the second.

## Usage

``` r
row_match(x, table)
```

## Arguments

- x:

  A matrix or data frame which rows to be matched. Typically, this will
  be a matrix of assigned barcodes for each read.

- table:

  A matrix or data frame with the rows to be matched against. Typically,
  this will be the top portion of a frequency table.

## Value

Logical vector, for each row in `x`, is the same row found in `table`?

## Details

As this function is intended to be used for data frames containing more
than just the barcodes, the intersection of the column names is used for
matching. As opposed to
[`base::match()`](https://rdrr.io/r/base/match.html), this function is
implemented more efficiently by converting each row into a numeric
encoding before matching.

For technical reasons, it is not permitted for the product of the number
of the unique values of the columns in `table` to exceed
\\2^{32}-1\approx 2.1\cdot 10^{9}\\.

## See also

[`create_freq_table()`](https://yaccos.github.io/posDemux/reference/create_freq_table.md)
for how frequency tables are constructed,
[`combinatorial_demultiplex()`](https://yaccos.github.io/posDemux/reference/combinatorial_demultiplex.md)
for more information on the matrix of assigned barcodes, and
[`dplyr::inner_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html)
for a function with similar functionality.

## Examples

``` r
barcode_table <- data.frame(
    read = c("seq_1", "seq_2", "seq_3", "seq_4"),
    bc1 = c("A", "B", "C", "B"),
    bc2 = c("A", "C", "A", "A")
    )

freq_table <- data.frame(
    bc1 = c("B", "B", "C", "A"),
    bc2 = c("A", "C", "A", "A"),
    frequency = c(200L, 100L, 50L, 10L)
    )

freq_cutoff <- 100L

selected_freq_table <- freq_table[freq_table$frequency >= freq_cutoff, ]

selected_rows <- row_match(barcode_table, selected_freq_table)
selected_barcode_table <- barcode_table[selected_rows, ]
```
