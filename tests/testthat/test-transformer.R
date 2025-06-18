library(Biostrings)
barcode_names <- c("bc3","bc2","bc1")
input <- get_tables()
table_trimmed <- input$table[barcode_names]
x_df <- as.data.frame(input$x, row.names = FALSE)
mapping <- posDemux:::get_mapping(table_trimmed)
x_encoded <- posDemux:::encode(x_df, mapping)
table_encoded <- posDemux:::encode(table_trimmed, mapping)
x_decoded <- posDemux:::decode(x_encoded, mapping)
table_decoded <- posDemux:::decode(table_encoded, mapping)
test_that("Encode and decode work as inverses of each other",
          {
            expect_equal(x_df, x_decoded, ignore_attr = FALSE)
            expect_equal(table_trimmed, table_decoded, ignore_attr = FALSE)
          }
          )

test_that("Row matching work as expected",
          {
            x_labeled <- x_df %>% dplyr::mutate(read = rownames(input$x))
            joint_table <- dplyr::inner_join(x_labeled, table_trimmed,
                                               by = barcode_names)
            expected_mask <- x_labeled$read %in% joint_table$read
            realized_mask <- row_match(input$x, input$table)
            expect_equal(realized_mask, expected_mask)
          }
          
          )
