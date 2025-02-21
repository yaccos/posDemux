#' Interactive cutoff selection
#' 
#' @description Launches an interactive Shiny application for
#' determining cutoff from knee plot and barcode frequency plot. The user
#' will select the appropriate number of barcode combinations
#' (i.e. distinct cells) to keep for further analysis. Usually, this is done
#' by aiming for the "knee" of the knee plot in order to keep most reads,
#' while at the same time removing barcode combinations which are either
#' artifacts or broken cells.
#' @import shiny
#'
#' @param frequency_table The frequency table
#' from \code{\link{create_frequency_table}}
#' @return Integer scalar with the selected cutoff returned
#' invisibly
#' @export
interactive_bc_cutoff <- function(frequency_table) {
  n_barcode_combinations <- frequency_table %>% nrow()
  n_reads <- frequency_table %>% pull("frequency") %>% sum()
  ui <- pageWithSidebar(
    headerPanel("Barcode cutoff selector"),
    sidebarPanel = sidebarPanel(
      sliderInput("cutoff", "Number of barcodes combinations to keep", min = 0L,
                  max = n_barcode_combinations,
                  step = 1L,
                  value = as.integer(n_barcode_combinations / 2L)
                  ),
      glue("Total number of reads: {n_reads}") %>% h4(),
      textOutput("reads_kept") %>% h4(),
      actionButton("exit", "Confirm cutoff selection")
    ),
    mainPanel = mainPanel(
      shiny::h4("Knee plot"),
      plotOutput("knee_plot"),
      shiny::h4("Barcode frequency plot"),
      plotOutput("frequency_plot")
    )
    )
  
  server <- function(input, output, session) {
    output$knee_plot <- renderPlot(
      knee_plot(frequency_table, input$cutoff) %>% print()
    )
    
    filtered_table <- reactive(
      frequency_table[seq_len(input$cutoff),]
    )

    frequency_cutoff <- reactive(
      bc_to_frequency_cutoff(frequency_table,
                                             input$cutoff)
    )
    
    reads_kept <- reactive(
      filtered_table() %>%
        pull("frequency") %>%
        sum()
    )
    
    percentage_kept <- reactive(
      round(reads_kept() / n_reads*100, 2)
    )
    
    
    
    output$frequency_plot <- renderPlot(
      frequency_plot(frequency_table,
                     cutoff = frequency_cutoff()) %>%
        print()
    )
    output$frequency_cutoff <- renderPrint(
      frequency_cutoff()
    )
    output$reads_kept <- renderText(
        glue("Number of reads kept: {reads_kept()} \\
              ({percentage_kept()}%)")
    )
    observeEvent(input$exit,
                 stopApp(input$cutoff)
    )
  
  }
  cutoff_value <- runApp(list(ui = ui, server = server))
  return(cutoff_value)
}


#' Convert between cutoff types
#' 
#' @description There are at least two ways to specify the cutoff to use when
#' selecting barcode combinations (cells) for further analysis. One way is 
#' to specify the number of barcode combinations to keep, effectively
#' keepting a given barcode combinations with the highest frequency. The other
#' way is to specify the frequency cutoff directly without regard to the number
#' of barcode combination to keep. In the former case,
#' \code{bc_to_frequency_cutoff()} is used to find the frequency cutoff, whereas
#' in the latter case \code{frequency_to_bc_cutoff()} is used to find the barcode
#' cutoff.
#'
#' @param frequency_table The frequency table
#' from \code{\link{create_frequency_table}}. In case the table is derived
#' from another source, it must be sorted in descending
#' order of frequency.
#' @param cutoff Integer vector, the cutoff values to be converted.
#' 
#' @details
#' In the edge case of the barcode threshold being zero, the frequency cutoff
#' is set to the maximum frequency in the table plus one. This feature
#' makes sure that the cutoff line is visible in the frequency plot.
#' 
#' @returns Integer, the converted cutoff values.
#' @export
#'
bc_to_frequency_cutoff <- function(frequency_table, cutoff) {
  frequency_table %>%
    pull("frequency") %>%
    c(max(.) + 1L, .) %>%
    extract(cutoff + 1L)
}

#' @rdname bc_to_frequency_cutoff
#' @export
frequency_to_bc_cutoff <- function(frequency_table, cutoff) {
  frequency_table$frequency %>% outer(cutoff, `>=`) %>%
    colSums()
}
  