# Copyright 2024 Masahiro Ono
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Generate basic QC plots for Tocky data
#'
#' This function visualizes either Timer fluorescence (Blue vs Red) or Timer dynamics
#' by the Tocky method (Angle vs Intensity) based on the specified mode.
#'
#' @param timer_transform_output A list object returned by `timer_transform`, containing `sample_definition`.
#' @param file File name.
#' @param pseudocolour A logical argument for whether to use pseudocolour in plots.
#' @param interactive Logical indicating whether to prompt the user to select plot_mode.
#'   Defaults to `TRUE`.
#' @param save A logical argument; if FALSE, plots are shown in an X window.
#' @param pdf A logical argument; if FALSE, a jpeg file is generated instead.
#' @param output The output directory name for output files.
#' @param n A number; n x n plots will be generated in the output Tocky plot file, max is 4 for 16 plots.
#' @param plot_mode Either "Timer fluorescence" for Blue vs Red plots, "Normalized Timer fluorescence" for normalized plots, or "Timer Angle and Intensity" for Angle vs Intensity plots.
#' @param group_order Optional character vector for specifying the order of the panels when using the group option.
#' @param lower_quantile_cutoff Lower quantile cutoff for setting the plot ranges in fluorescence mode.
#' @param select Logical indicating whether to manually select samples for plotting.
#' @param group Logical indicating whether to group plots based on the `group` field in `sample_definition`.
#' @param verbose Logical indicating whether to print progress messages.
#'        Default is `TRUE`.
#' @param samplefile Character vector specifying the sample files. Defaults to `NULL`.
#' @export
#' @examples
#' \dontrun{
#'   plot_tocky(x, plot_mode = "Timer fluorescence")
#'   plot_tocky(x, plot_mode = "Timer Angle and Intensity")
#' }
#' @importFrom grDevices jpeg densCols colorRampPalette rainbow dev.off
#' @importFrom stats quantile

plot_tocky <- function(timer_transform_output, file = 'PlotTocky', pseudocolour = TRUE, pdf = FALSE, output = 'QC', n = 4, plot_mode = "Timer fluorescence", lower_quantile_cutoff = 0.01, select = FALSE, group = TRUE, group_order = NULL, interactive = TRUE, save = FALSE,
    samplefile = NULL, verbose = TRUE) {
    n <- min(n, 4)
    
    if (!file.exists(output)) {
        dir.create(output)
    }
    
    if (interactive) {

        plot_mode_choices <-  c("Timer fluorescence", "Normalized Timer fluorescence", "Timer Angle and Intensity")
        plot_mode <- utils::select.list(choices = plot_mode_choices, title = "Select plot mode:", graphics = FALSE)
    }
    
    samples <- timer_transform_output$sample_definition$file

    if(is.null(samplefile)){
        if (select) {
            selected_samples <- utils::select.list(samples, graphics = TRUE, title = "Choose samples for plotting", multiple = TRUE)
            data <- timer_transform_output$transformed_data[timer_transform_output$transformed_data$file %in% selected_samples, ]
        } else {
            data <- timer_transform_output$transformed_data
        }
        
    }else{
        if(all(samplefile %in%  timer_transform_output$transformed_data$file)){
            data <- timer_transform_output$transformed_data[timer_transform_output$transformed_data$file %in% samplefile, ]
            
        }else{
            
            stop("Enter only sample files that exist in the data. \n")
        }
    }

    
    data <- merge(data, timer_transform_output$sample_definition, by = 'file')
    if(verbose){
        cat("Now plotting..\n")
        }
    
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

    
    data$file <- as.factor(data$file)
    data$group <- as.factor(data$group)
    
    if (group) {
        if (!is.null(group_order)) {
            data$group <- factor(data$group, levels = group_order)
        }
        plot_data <- split(data, data$group)
    } else {
        plot_data <- split(data, data$file)
    }

    if(save){
        if (pdf) {
            pdf(file = file.path(output, paste0(file, ".pdf")), width = 480 * 3, height = 480 * 3)
            par(mar = c(2, 2, 2, 2))
        } else {
            jpeg(filename = file.path(output, paste0(file, ".jpeg")), width = 480 * 3, height = 480 * 3)
            par(mar = c(5, 5, 5, 5))
        }
        
        par(mfrow = c(n, n))
        for (i in 1:min(16, length(plot_data))) {
            
            tpdata <- plot_data[[i]][, c(plotting_params$x_var, plotting_params$y_var)]
            colnames(tpdata) <- c('x', 'y')
            if (pseudocolour) {
                density <- densCols(tpdata$x, tpdata$y, colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
            } else {
                density <- 8  # A fixed color if pseudocolour is not used
            }
            plot(tpdata, xlab = plotting_params$x_label, ylab = plotting_params$y_label, pch = 19, cex = 0.2, col = density, xlim = plotting_params$xlim, ylim = plotting_params$ylim, cex.main = 2, cex.lab = 2, main = names(plot_data)[i])
            
            if(plot_mode == "Timer fluorescence"){
                abline(
                v = timer_transform_output$normalization_parameters$red_threshold,
                h = timer_transform_output$normalization_parameters$blue_threshold,
                col = 2
                )
            }
        }

    }else{

        par(mar = c(4, 4, 2.5, 2.5))
        
        par(mfrow = c(n, n))
        for (i in 1:min(16, length(plot_data))) {
            
            tpdata <- plot_data[[i]][, c(plotting_params$x_var, plotting_params$y_var)]
            colnames(tpdata) <- c('x', 'y')
            if (pseudocolour) {
                density <- densCols(tpdata$x, tpdata$y, colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
            } else {
                density <- 8  # A fixed color if pseudocolour is not used
            }
            plot(tpdata, xlab = plotting_params$x_label, ylab = plotting_params$y_label, pch = 19, cex = 0.2, col = density, xlim = plotting_params$xlim, ylim = plotting_params$ylim, cex.main = 0.8, cex.lab = 0.8, main = names(plot_data)[i])
            
            if(plot_mode == "Timer fluorescence"){
                abline(
                v = timer_transform_output$normalization_parameters$red_threshold,
                h = timer_transform_output$normalization_parameters$blue_threshold,
                col = 2
                )
            }
        }
    }



    if(save){
        dev.off()
    }
    if(verbose){
        cat("Plotting complete.\n")
    }
    
}
