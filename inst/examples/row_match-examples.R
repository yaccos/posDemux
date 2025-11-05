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
