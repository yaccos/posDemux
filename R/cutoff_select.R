#' Interactive cutoff selection
#'
#' @description Returns an
#' interactive Shiny application for
#' determining cutoff from knee plot and barcode frequency plot. The user
#' will select the appropriate number of barcode combinations
#' (i.e. distinct cells) to keep for further analysis. Usually, this is done
#' by aiming for the 'knee' of the knee plot in order to keep most reads,
#' while at the same time removing barcode combinations which are either
#' artifacts or broken cells.
#' @import shiny
#'
#' @param freq_table The frequency table
#' from [create_freq_table()].
#' @return A \code{\link[shiny:shinyApp]{shiny.appobj}} which launches
#' when printed and returns the
#' last selected cutoff (invisibly) when it stops.
#'
#' @example inst/examples/interactive_bc_cutoff-examples.R
#'
#' @export
interactive_bc_cutoff <- function(freq_table) {
    n_reads <- freq_table %>%
        pull("frequency") %>%
        sum()
    ui <- cutoff_select_ui(freq_table, n_reads)
    server <- cutoff_select_server(freq_table, n_reads)
    shinyApp(ui, server)
}

cutoff_select_ui <- function(freq_table, n_reads) {
    n_barcode_combinations <- freq_table %>%
        nrow()
    n_reads <- freq_table %>%
        pull("frequency") %>%
        sum()
    pageWithSidebar(
        headerPanel("Barcode cutoff selector"),
        sidebarPanel = sidebarPanel(
            sliderInput(
                "cutoff",
                "Number of barcodes combinations to keep",
                min = 0L,
                max = n_barcode_combinations,
                step = 1L,
                value = as.integer(n_barcode_combinations / 2L)
            ),
            glue("Total number of reads: {n_reads}") %>%
                h4(),
            textOutput("reads_kept") %>%
                h4(),
            textOutput("frequency_cutoff") %>%
                h4(),
            radioButtons(
                "freq_plot_type",
                "Type of frequency plot",
                choices = c(Density = "density", Histogram = "histogram"),
                selected = "density"
            ),
            checkboxInput("log_scale_x", "Log scale on frequency plot x-axis"),
            checkboxInput("log_scale_y", "Log scale on frequency plot y-axis"),
            checkboxInput(
                "scale_by_reads",
                "Scale y-axis in frequency plot by number of reads"
            ),
            actionButton("exit", "Confirm cutoff selection")
        ),
        mainPanel = mainPanel(
            shiny::h4("Knee plot"),
            plotOutput("knee_plot"),
            shiny::h4("Barcode frequency plot"),
            plotOutput("frequency_plot")
        )
    )
}

cutoff_select_server <- function(freq_table, n_reads) {
    function(input, output, session) {
    output$knee_plot <- renderPlot(
        knee_plot(freq_table, input$cutoff) + 
            theme(text = element_text(size = 18))
        )
    
    filtered_table <- reactive(freq_table[seq_len(input$cutoff), ])
    
    frequency_cutoff <- reactive(
        bc_to_freq_cutoff(freq_table, input$cutoff)
    )
    reads_kept <- reactive(
        filtered_table() %>% pull("frequency") %>% sum()
    )
    
    percentage_kept <- reactive(round(reads_kept() / n_reads * 100, 2))
    
    
    output$frequency_plot <- renderPlot(
        freq_plot(
            freq_table,
            cutoff = frequency_cutoff(),
            type = input$freq_plot_type,
            log_scale_x = input$log_scale_x,
            log_scale_y = input$log_scale_y,
            scale_by_reads = input$scale_by_reads
        ) + theme(text = element_text(size = 18))
    )
    output$frequency_cutoff <- renderPrint(frequency_cutoff())
    output$reads_kept <- renderText(glue(
        "Number of reads kept: {reads_kept()} \\
            ({percentage_kept()}%)"
    ))
    output$frequency_cutoff <- renderText(glue(
        "Number of reads in last barcode included: {frequency_cutoff()}"
    ))
    observeEvent(input$exit, stopApp(input$cutoff))
    }
}




#' Convert between cutoff types
#'
#' @description There are at least two ways to specify the cutoff to use when
#' selecting barcode combinations (cells) for further analysis. One way is
#' to specify the number of barcode combinations to keep, effectively
#' keeping a given number of barcode combinations with the highest frequencies.
#' The other way is to specify the frequency cutoff directly
#' without regard to the number
#' of barcode combination to keep. In the former case,
#' \code{bc_to_freq_cutoff()} is used to find
#' the corresponding frequency cutoff,
#' whereas in the latter case \code{freq_to_bc_cutoff()}
#' is used to find the corresponding barcode cutoff.
#'
#' @param freq_table The frequency table
#' from [create_freq_table()]. In case the table is derived
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
#' @example inst/examples/cutoff_conversion-examples.R
#' @export
#'
bc_to_freq_cutoff <- function(freq_table, cutoff) {
    freq_table %>%
        pull("frequency") %>%
        c(max(.) + 1L, .) %>%
        extract(cutoff + 1L)
}

#' @rdname bc_to_freq_cutoff
#' @export
freq_to_bc_cutoff <- function(freq_table, cutoff) {
    freq_table$frequency %>%
        outer(cutoff, `>=`) %>%
        colSums()
}
