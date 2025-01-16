#' Launch a Shiny App for Exploring timer_transform Parameter Space
#'
#' This function launches a Shiny application that allows users to interactively explore
#' the parameter space of the `timer_transform` function. Users can adjust thresholds
#' and normalization methods to see how these changes affect the transformation of
#' flow cytometry data.
#'
#' @param prep A prep object containing file paths and variables, typically the output
#' from \code{\link{prep_tocky}}.
#' @param transformed_data A \code{TockyPrepData} object
#' @description The Shiny application provides a user interface for adjusting
#' the blue and red fluorescence thresholds and choosing between normalization methods.
#' It updates visualizations based on user inputs to aid in determining optimal parameters
#' for data analysis.
#'
#' @return Does not return a value; a Shiny app is launched in the default web browser.
#'
#' @examples
#' \dontrun{
#'   explore_timer_transform(prep_data, transformed_data)
#' }
#' @importFrom shiny fluidPage titlePanel sidebarLayout sidebarPanel sliderInput radioButtons actionButton mainPanel tabsetPanel tabPanel plotOutput verbatimTextOutput eventReactive renderPlot renderPrint shinyApp
#' @importFrom graphics screen split.screen close.screen
#' @export

explore_timer_transform <- function(prep, transformed_data) {
    blue_channel_name <- transformed_data@timer_fluorescence$original_blue_channel
    red_channel_name <- transformed_data@timer_fluorescence$original_red_channel
    
    blue_vals <- transformed_data@Data$Blue_log
    red_vals  <- transformed_data@Data$Red_log
    
    q_bl <- round(quantile(blue_vals, probs = c(0.1, 1.0), na.rm = TRUE), digits = 2)
    q_rd <- round(quantile(red_vals,  probs = c(0.1, 1.0), na.rm = TRUE), digits = 2)
    step_blue <- round((q_bl[2] - q_bl[1]) / 50, digits = 2)
    step_red  <- round((q_rd[2] - q_rd[1]) / 50, digits = 2)
    if(step_blue <= 0) step_blue <- abs(q_bl[2] - q_bl[1]) / 50 + 0.1
    if(step_red  <= 0) step_red  <- abs(q_rd[2] - q_rd[1]) / 50 + 0.1
    blue_mid <- transformed_data@normalization_parameters$blue_threshold
    red_mid  <- transformed_data@normalization_parameters$red_threshold
    
    ui <- fluidPage(
    
    titlePanel("Exploring Timer Data Transformation Parameters"),
    
    sidebarLayout(
    sidebarPanel(
    
    sliderInput("blueThresh", "Timer Blue Threshold:",
    min = q_bl[1], max = q_bl[2], value = blue_mid, step = step_blue),
    sliderInput("redThresh", "Timer Red Threshold:",
    min = q_rd[1], max = q_rd[2], value = red_mid, step = step_red),
    
    radioButtons("normalization", "Normalization:",
    choices = c("TRUE", "FALSE"),
    selected = "TRUE"),
    
    radioButtons("normMethod", "Normalization Method:",
    choices = c("MAD", "SD"),
    selected = "MAD"),
    
    actionButton("applyTransform", "Apply Transform")
    ),
    
    mainPanel(
    
    tabsetPanel(
    tabPanel("Plot", plotOutput("timerPlot", height = "800px")),
    tabPanel("Summary", verbatimTextOutput("summaryText"))
    )
    )
    )
    )
    
    server <- function(input, output, session) {
        
        timerReactive <- eventReactive(input$applyTransform, {
            
            result <- timer_transform(
            prep         = prep,
            select       = FALSE,
            blue_channel = blue_channel_name,
            red_channel  = red_channel_name,
            red_threshold  = input$redThresh,
            blue_threshold = input$blueThresh,
            normalization_method = input$normMethod,
            normalization = input$normalization,
            interactive_gating  = FALSE,
            verbose           = FALSE,
            use_negative_control = FALSE
            )
            return(result)
        })
        
        
        output$timerPlot <- renderPlot({
            transformedData <- timerReactive()
            
            split.screen(
            figs = matrix(
            c(
            0,    0.5,  0.5,  1,
            0.5,  1,    0.5,  1,
            0,    0.5,  0,    0.5,
            0.5,  1,    0,    0.5
            ),
            ncol = 4,
            byrow = TRUE
            )
            )
            
            screen(1)
            plot_tocky_shiny(
            transformedData,
            plot_mode = "Timer fluorescence"
            )
            # plot_timer_gating
            screen(2)
            plot_tocky_shiny(
            transformedData,
            plot_mode = "Normalized Timer fluorescence"
            )
            
            screen(3)
            plot_tocky_shiny(
            transformedData,
            plot_mode = "Timer Angle and Intensity"
            )
            screen(4)
            plot_timer_gating(prep, transformedData)
            
            close.screen(all.screens = TRUE)
        })
        
        output$summaryText <- renderPrint({
            transformedData <- timerReactive()
            
            cat("Current Parameter Selections:\n")
            cat(paste("Blue Threshold: ", input$blueThresh, "\n"))
            cat(paste("Red Threshold:  ", input$redThresh, "\n"))
            cat(paste("Normalization Method: ", input$normMethod, "\n\n"))
            
            cat("\nTimer Blue Fluorescence Channel:\n")
            print(blue_channel_name)
            
            cat("\nTimer Red Fluorescence Channel:\n")
            print(red_channel_name)
            
            cat("\nNormalization Parameters:\n")
            print(transformedData@normalization_parameters)
            
        })
    }
    
    shinyApp(ui = ui, server = server)
}

#' Plot Timer Data for Shiny Visualization
#'
#' This function provides visualizations for timer data within a Shiny application,
#' showing different modes such as "Timer fluorescence", "Normalized Timer fluorescence",
#' or "Timer Angle and Intensity". It is designed to be used internally by Shiny apps
#' to dynamically display changes based on user input parameters.
#'
#' @param x A `TockyPrepData` object from `timer_transform` containing processed timer data.
#' @param plot_mode A character string indicating the type of plot to display.
#'        Options are "Timer fluorescence", "Normalized Timer fluorescence", or "Timer Angle and Intensity".
#'        Defaults to "Timer fluorescence".
#' @param lower_quantile_cutoff A numeric value indicating the lower quantile cutoff
#'        for setting the plot ranges in fluorescence mode. This helps in adjusting the scale
#'        of the plot to focus on relevant data points. Defaults to 0.01.
#'
#' @details The function checks if the input data inherits from `TockyPrepData`.
#'          It then uses base R plotting functions to create plots according to the specified `plot_mode`.
#'          The function adjusts the x and y limits based on the provided `lower_quantile_cutoff`
#'          and uses density colors for better visualization of data point distributions.
#'
#' @return This function does not return a value but renders a plot directly to the Shiny application.
#'         It adjusts plot elements such as labels, axis limits, and coloring based on the input parameters.
#'
#' @examples
#' \dontrun{
#'   plot_tocky_shiny(x, plot_mode = "Timer fluorescence")
#' }
#'
#' @importFrom grDevices densCols colorRampPalette
#' @importFrom graphics plot abline
#' @importFrom stats quantile
#' @keywords internal

plot_tocky_shiny <- function(x, plot_mode = "Timer fluorescence", lower_quantile_cutoff = 0.01) {
    if(!inherits(x, "TockyPrepData")){
        stop("Use the output of timer_transform. \n")
    }
    data <- x@Data

    if (plot_mode == "Timer fluorescence") {
        plotting_params <- list(
        x_var = 'Red_log',
        y_var = 'Blue_log',
        x_label = "Timer Red",
        y_label = "Timer Blue",
        xlim = quantile(data[['Red_log']], c(lower_quantile_cutoff, 1)),
        ylim = quantile(data[['Blue_log']], c(lower_quantile_cutoff, 1))
        )
    } else if (plot_mode == "Normalized Timer fluorescence") {
        plotting_params <- list(
        x_var = 'Red_Normalized',
        y_var = 'Blue_Normalized',
        x_label = "Timer Red Normalized",
        y_label = "Timer Blue Normalized",
        xlim = quantile(data[['Red_Normalized']], c(lower_quantile_cutoff, 1)),
        ylim = quantile(data[['Blue_Normalized']], c(lower_quantile_cutoff, 1))
        )
    } else  {
        
        plotting_params <- list(
        x_var = 'Angle',
        y_var = 'Intensity',
        x_label = "Timer Angle",
        y_label = "Timer Intensity",
        xlim = c(0, 90),
        ylim = c(0, max(data[['Intensity']], na.rm = TRUE))
        )
    }
    tpdata <- data[, c(plotting_params$x_var, plotting_params$y_var)]
    colnames(tpdata) <- c('x', 'y')
    tpdata <- tpdata[is.finite(tpdata$x) & is.finite(tpdata$y), ]
    
    density <- suppressWarnings(densCols(tpdata$x, tpdata$y, colramp = colorRampPalette(rev(rainbow(10, end = 4/6)))))
    
    plot(tpdata, xlab = plotting_params$x_label, ylab = plotting_params$y_label, pch = 19, cex = 0.2,
    col = density, xlim = plotting_params$xlim, ylim = plotting_params$ylim,
    main = plot_mode)
    
    if (plot_mode == "Timer fluorescence") {
        abline(
        v = x@normalization_parameters$red_threshold,
        h = x@normalization_parameters$blue_threshold,
        col = 2
        )
    }
}
