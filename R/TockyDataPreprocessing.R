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

#' A class representing a TockyPrepData object for output of timer_transform
#'
#' This class is designed to encapsulate and structure the output of the
#' timer_transform function in the TockyPrep package.
#'
#' @slot Data A data.frame containing expression data.
#' @slot cell_counts A data.frame containing counts of cells per sample.
#' @slot sampledef A list including annotation data for sample grouping.
#' @slot timer_fluorescence A list containing channel names for fluorescence timer data.
#' @slot input A list of parameters used for creating TockyPrepData object.
#' @slot normalization_parameters A list of parameters used for data normalization.
#' @slot Tocky A list containing other Tocky-specific analysis data.
#'
#' @keywords classes
#' @export
#' @importFrom methods setClass
setClass(
  "TockyPrepData",
  slots = list(
    Data = "data.frame",
    cell_counts = "data.frame",
    sampledef = "list",
    timer_fluorescence = "list",
    input = "list",
    normalization_parameters = "list",
    Tocky = "list"
  ),
  validity = function(object) {
    problems <- NULL
    if (!is.data.frame(object@Data)) {
      problems <- c(problems, "Data must be a data.frame.")
    }
    if (!is.data.frame(object@cell_counts)) {
      problems <- c(problems, "cell_counts must be a data.frame.")
    }
    if (!is.list(object@sampledef)) {
      problems <- c(problems, "sampledef must be a list.")
    }
    if (!is.list(object@timer_fluorescence)) {
      problems <- c(problems, "timer_fluorescence must be a list.")
    }
    if (!is.list(object@input)) {
      problems <- c(problems, "input must be a list.")
    }
    if (!is.list(object@normalization_parameters)) {
      problems <- c(problems, "normalization_parameters must be a list.")
    }
    if (!is.list(object@Tocky)) {
      problems <- c(problems, "Tocky must be a list.")
    }
    if (length(problems) == 0) TRUE else problems
  }
)


#' Initialize a TockyPrepData object
#'
#' This method initializes a TockyPrepData with specific data for flow cytometry analysis.
#' It ensures all specific slots are set up.
#'
#' @param .Object A TockyPrepData object to be initialized.
#' @param Data A data.frame containing expression data.
#' @param cell_counts A data.frame containing counts of cells per sample.
#' @param sampledef A list including annotation data for sample grouping.
#' @param timer_fluorescence A list containing channel names for fluorescence timer data.
#' @param input A list of parameters used for creating TockyPrepData object.
#' @param normalization_parameters A list of parameters used for data normalization.
#' @param Tocky A list containing other Tocky-specific analysis data.
#'
#' @return A valid TockyPrepData that has been initialized with provided data.
#' @keywords internal
#' @importFrom methods initialize validObject
#' @export
setMethod("initialize", "TockyPrepData",
  function(.Object, Data, cell_counts, sampledef, timer_fluorescence, input, normalization_parameters, Tocky) {
    .Object@Data <- Data
    .Object@cell_counts <- cell_counts
    .Object@sampledef <- sampledef
    .Object@timer_fluorescence <- timer_fluorescence
    .Object@input <- input
    .Object@normalization_parameters <- normalization_parameters
    .Object@Tocky <- Tocky

    validObject(.Object)

    return(.Object)
  }
)


#' Show method for the TockyPrepData class
#'
#' Displays summary information for various slots of the TockyPrepData object.
#' Includes details such as total number of cells, variable names, sample numbers,
#' and group levels, providing a concise summary of the object.
#'
#' @param object An object of the TockyPrepData class
#' @export
#' @importFrom methods show
#' @importFrom utils head
setMethod("show", "TockyPrepData", function(object) {
    cat("TockyPrepData Object:\n")
    if (length(object@Data) > 0) {
        cat(paste("Total cell number:", nrow(object@Data), "\n"))
        cat("Variables: ", paste(head(colnames(object@Data)), collapse=", "), "\n")
    }
    
    if (length(object@sampledef) > 0) {
        cat(paste("Total sample number:", nrow(object@sampledef), "\n"))
        if ("group" %in% colnames(object@sampledef)) {
            cat("Groups: ", paste(levels(as.factor(object@sampledef[['group']])), collapse=", "), "\n")
        }
    }

    if (length(object@Tocky) > 0) {
         cat("Tocky Data: \n")
         print(names(object@Tocky))
     }
    cat("\n")
})




#' Prepare Data for Timer Transformation Using Flow Cytometric Data
#'
#' This function prepares the dataset for timer transformation analysis by identifying common variables
#' across sample files, configuring necessary control files, and setting up variables for downstream
#' analysis. The function supports both interactive and non-interactive file selection modes.
#'
#' @param path Character string specifying the directory where the data files are located.
#'   Defaults to the current directory `'.'`.
#' @param interactive Logical indicating whether to prompt the user to select sample files.
#'   Defaults to `TRUE`.
#' @param negfile Character string specifying the negative control file. If `NULL`, the user will
#'   be prompted to select a file. Defaults to `NULL`.
#' @param samplefile Character vector specifying the sample files. If `NULL` and `samplefilechoice`
#'   is `TRUE`, the user will be prompted to select files. Defaults to `NULL`.
#'
#' @return A list containing paths to the control file, selected sample files, and the standardized
#'   variables used in the analysis.
#' @export
#'
#' @examples
#' \dontrun{
#'   # Interactive file selection
#'   prep_data <- prep_tocky(path='data', output='output')
#'
#'   # Specifying files directly for non-interactive usage
#'   prep_data <- prep_tocky(
#'     path='data',
#'     output='output',
#'     negfile='neg_control.csv',
#'     samplefile=c('sample1.csv', 'sample2.csv')
#'   )
#' }

prep_tocky <- function(path = '.', interactive = TRUE,
negfile = NULL, samplefile = NULL) {
    files_in_path <- list.files(path, pattern = "\\.csv$", full.names = TRUE)
    

    
    if (interactive) {

            negfile <- utils::select.list(files_in_path, graphics = TRUE,
            title = "Choose a Timer negative control file", multiple = FALSE)

            samplefile <- utils::select.list(files_in_path, graphics = TRUE,
            title = "Choose sample CSV files to be analysed", multiple = TRUE)

    } else {
        if(is.null(negfile)|is.null(samplefile) ){
            stop("Include input files if you use interactive mode. \n")
        }else{
            negfile <- negfile
            samplefile <- samplefile
            
        }
    }
    
    neg <- read.table(file.path(path, negfile), nrows = 2, sep = ',', header = TRUE)
    lf <- samplefile
    common_variables <- colnames(neg)
    
    for (i in seq_along(samplefile)) {
        first_line <- read.table(file.path(path, lf[i]), nrows = 2, sep = ',', header = TRUE)
        common_variables <- intersect(common_variables, colnames(first_line))
    }
    
    vard <- colnames(neg)
    lg <- vard %in% common_variables
    vard[!lg] <- "Inconsistent"
    if (any(!lg)) {
        message("Warning - some variables are not consistent between samples!")
    }
    
    vard[grepl(pattern = 'Red', vard, ignore.case = TRUE)] <- "Timer_Red"
    vard[grepl(pattern = 'Blue', vard, ignore.case = TRUE)] <- "Timer_Blue"
    
    vardf <- data.frame(Channel.name = colnames(neg), Variable = vard, stringsAsFactors = FALSE)
    
    output_list <- list(
    neg = negfile,
    samplefile = samplefile,
    variables = vardf$Channel.name[vardf$Variable != "Inconsistent"],
    path = path
    )
    
    return(output_list)
}



#' Perform Timer Transformation on Flow Cytometry Data
#'
#' This function processes flow cytometry data by applying Timer thrsholding, normalization,
#' and trigonometric transformation to the Blue and Red fluorescence data.
#'
#' @param prep A list containing file paths and variables, typically the output from \code{\link{prep_tocky}}.
#' @param blue_channel Character string specifying the Blue fluorescence channel name.
#'        If `NULL`, the function attempts to determine it automatically.
#' @param red_channel Character string specifying the Red fluorescence channel name.
#'        If `NULL`, the function attempts to determine it automatically.
#' @param red_threshold Numeric specifying the Red channel gating threshold.
#'        If `NULL`, gating is performed automatically or interactively based on `interactive_gating`.
#' @param blue_threshold Numeric specifying the Blue channel gating threshold.
#'        If `NULL`, gating is performed automatically or interactively based on `interactive_gating`.
#' @param interactive_gating Logical indicating whether to perform interactive gating when thresholds are not provided.
#'        Default is `FALSE`.
#' @param verbose Logical indicating whether to print progress messages.
#'        Default is `TRUE`.
#' @param select Logical indicating whether to choose Timer fluorescence channels interactively.
#'        Default is `TRUE`.
#' @param normalization Logical indicating whether to apply Timer fluorescence normalization.
#'        Default is `TRUE`.
#' @param normalization_method Charcter string specifying the normalization method to be applied to Timer flfuorescence
#'        Default is `MAD`. The alternative is `SD`.
#' @param q Quantile value used for automatic Timer thresholds.
#'        Default is 0.998.
#'
#' @return The function returns a new TockyPrepData object containing:
#'   \itemize{
#'     \item \code{Data}: Data frame with angle, intensity, and other variables.
#'     \item \code{normalization_parameters}: List with normalization values and coefficients.
#'     \item \code{cell_counts}: Data frame with cell counts for each sample.
#'     \item \code{sampledef}: Data frame with sample file names and placeholder group column.
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#'   # Assuming `prep_data` is the output from `prep_tocky`
#'   result <- timer_transform(prep_data)
#' }
#'
#' @importFrom utils read.csv
#' @importFrom grDevices rgb
#' @importFrom methods new
#' @importFrom graphics abline locator par
#' @importFrom stats coef density lm sd quantile mad

timer_transform <- function(prep, select = TRUE, blue_channel = NULL, red_channel = NULL, normalization_method = 'MAD',
red_threshold = NULL, blue_threshold = NULL, interactive_gating = FALSE, verbose = TRUE, q = 0.998, normalization = TRUE) {
    
    input_list <- list(
        path = prep$path,
        select = select,
        blue_channel = blue_channel,
        red_channel = red_channel,
        xlim_red = red_threshold,
        ylim_blue = blue_threshold,
        interactive_gating = interactive_gating,
        verbose = verbose,
        quantile = q,
        normalization = normalization,
        normalization_method = normalization_method
    )
    
    LogSingleData <- function(x){
        x.log <- x
        lg <- x.log > 1
        x.log[lg] <- log10(x.log[lg])
        x.log[!lg] <- 0
        return(x.log)
    }
    
    if(!all(c("neg","path","samplefile", "variables") %in% names(prep))){
        stop("Use a correct prep input. \n")
    }
    
    negfile <- prep$neg
    path <- prep$path
    samplefile <- prep$samplefile
    variables <- prep$variables
    
    if(select){
        blue_channel <- utils::select.list(variables, graphics = TRUE, title = "Choose Timer Blue", multiple = FALSE)
        red_channel <- utils::select.list(variables, graphics = TRUE, title = "Choose Timer Red", multiple = FALSE)
    } else {
        if (is.null(blue_channel) || is.null(red_channel)) {
            choices <- variables
            if (is.null(blue_channel)) {
                blue_channel <- choices[grep("Timer_Blue", choices)]
                if (length(blue_channel) != 1) {
                    stop("Blue channel not specified and could not be determined automatically.")
                }
            }
            if (is.null(red_channel)) {
                red_channel <- choices[grep("Timer_Red", choices)]
                if (length(red_channel) != 1) {
                    stop("Red channel not specified and could not be determined automatically.")
                }
            }
        }
    }
    
    if (verbose) {
        message("Blue channel: ", blue_channel)
        message("Red channel: ", red_channel)
     }
    
    neg <- read.csv(file.path(path, negfile), header = TRUE)
    
    required_channels <- c(blue_channel, red_channel)
    if (!all(required_channels %in% colnames(neg))) {
        stop("Not all required channels are present in the negative control data.")
    }
    
    red_min_adjust <- -quantile(neg[[red_channel]], 0.005, na.rm = TRUE) + 100
    blue_min_adjust <- -quantile(neg[[blue_channel]], 0.005, na.rm = TRUE) + 100
    
    # New column names
    red_log_colname <- ifelse("Red_log" %in% colnames(neg), paste("Red_log", "transformed", sep = "_"), "Red_log")
    blue_log_colname <- ifelse("Blue_log" %in% colnames(neg), paste("Blue_log", "transformed", sep = "_"), "Blue_log")

    # Log transform and assign to new columns
    neg[[red_log_colname]] <- LogSingleData(neg[[red_channel]] + red_min_adjust)
    neg[[blue_log_colname]] <- LogSingleData(neg[[blue_channel]] + blue_min_adjust)

    if (verbose) {
        if(red_log_colname != "Red_log" || blue_log_colname != "Blue_log") {
            message("Column names 'Red_log' or 'Blue_log' already exist. Renamed to ", red_log_colname, " and ", blue_log_colname, ".")
        }
        message("Processed Timer Fluorescence: ", red_log_colname, " and ", blue_log_colname)
    }
        
    if (is.null(red_threshold) || is.null(blue_threshold)) {
        if (interactive_gating) {
            if (verbose) message("Interactive gating started. Please click on the plot to set thresholds.")
            plot(neg$Red_log, neg$Blue_log, xlab='Timer Red (log)', ylab='Timer Blue (log)',
            pch='.', col=rgb(0,0,0,alpha=0.2))
            cat("Set a threshold for Blue and Red by clicking on the plot: \n")
            scalegate <- locator(n = 1, type = 'p', col=2)
            if (!is.null(scalegate)) {
                red_threshold <- scalegate$x
                blue_threshold <- scalegate$y
                abline(v = red_threshold, h = blue_threshold, col=2)
                if (verbose) message(paste("Gating thresholds set at red_threshold =", red_threshold,
                "and blue_threshold =", blue_threshold))
            } else {
                stop("No point selected for gating thresholds.")
            }
        } else {
            red_threshold <- quantile(neg$Red_log, q, na.rm = TRUE)
            blue_threshold <- quantile(neg$Blue_log, q, na.rm = TRUE)
            if (verbose) message("Automatic gating applied.")
        }
    }
    if (verbose) {
        message("Gating thresholds:")
        message("red_threshold: ", red_threshold)
        message("blue_threshold: ", blue_threshold)
    }
    
    gate_filter <- (neg$Red_log < red_threshold) & (neg$Blue_log < blue_threshold)
    neg_gated <- neg[gate_filter, ]
    
    if (nrow(neg_gated) == 0) {
        stop("No data points remain after gating. Please adjust the gating thresholds.")
    }
    
    blue_channel_normalized <- "Blue_Normalized"
    red_channel_normalized <- "Red_Normalized"

    neg_normalized <- neg
    neg_normalized[[blue_channel_normalized]] <- neg$Blue_log
    neg_normalized[[red_channel_normalized]] <- neg$Red_log
    
    dataset_list <- list()
    cellcount_total <- integer(length(samplefile))
    names(cellcount_total) <- samplefile
    
    for (i in seq_along(samplefile)) {
        tmpdata <- read.csv(file.path(path, samplefile[i]),  header = TRUE)
        cellcount_total[i] <- nrow(tmpdata)
        
        if (!all(required_channels %in% colnames(tmpdata))) {
            stop(paste("Required channels not found in sample file:", samplefile[i]))
        }
        
        tmpdata$Red_log <- LogSingleData(tmpdata[[red_channel]] + red_min_adjust)
        tmpdata$Blue_log <- LogSingleData(tmpdata[[blue_channel]] + blue_min_adjust)
        
        tmpdata[[blue_channel_normalized]] <- tmpdata$Blue_log
        tmpdata[[red_channel_normalized]] <- tmpdata$Red_log
        
        dataset_list[[samplefile[i]]] <- tmpdata
        if (verbose) message(paste("Processed sample file:", samplefile[i]))
    }
    
    max_neg_blue <- max(neg_normalized[[blue_channel_normalized]], na.rm = TRUE)
    max_neg_red <- max(neg_normalized[[red_channel_normalized]], na.rm = TRUE)

    if (normalization_method == "MAD") {
        dispersion_blue <- mad(neg_normalized[[blue_channel_normalized]], na.rm = TRUE)
        dispersion_red <- mad(neg_normalized[[red_channel_normalized]], na.rm = TRUE)
        dispersion_label <- "MAD"
    } else {  # Default to SD
        dispersion_blue <- sd(neg_normalized[[blue_channel_normalized]], na.rm = TRUE)
        dispersion_red <- sd(neg_normalized[[red_channel_normalized]], na.rm = TRUE)
        dispersion_label <- "SD"
    }

    if (verbose) {
        message("Normalization parameters:")
        message("Max negative Blue: ", max_neg_blue)
        message(sprintf("%s negative Blue: ", dispersion_label), dispersion_blue)
        message("Max negative Red: ", max_neg_red)
        message(sprintf("%s negative Red: ", dispersion_label), dispersion_red)
    }

    process_fluorescence <- function(B, R, maxB, sdB, maxR, sdR, applyNormalization = TRUE) {
        B_normalized <- pmax(B, maxB)
        R_normalized <- pmax(R, maxR)

        if (applyNormalization) {
            B_normalized <- (B_normalized - maxB) / sdB
            R_normalized <- (R_normalized - maxR) / sdR
        } else {
             B_normalized <- B_normalized - maxB
            R_normalized <- R_normalized - maxR
        }
        
        return(data.frame(Blue = B_normalized, Red = R_normalized))
    }

    result_all <- data.frame()
    
    max_neg_blue <- max(neg_normalized[[blue_channel_normalized]], na.rm = TRUE)
    max_neg_red <- max(neg_normalized[[red_channel_normalized]], na.rm = TRUE)
    sd_blue <- sd(neg_normalized[[blue_channel_normalized]], na.rm = TRUE)
    sd_red <- sd(neg_normalized[[red_channel_normalized]], na.rm = TRUE)
    MAD_blue <- mad(neg_normalized[[blue_channel_normalized]], na.rm = TRUE)
    MAD_red <- mad(neg_normalized[[red_channel_normalized]], na.rm = TRUE)

    for (sample_name in names(dataset_list)) {
        data_norm <- dataset_list[[sample_name]]
        
        scale_parameters <- if (normalization_method == 'MAD') {
            list(sdB = MAD_blue, sdR = MAD_red)
        } else {
            list(sdB = sd_blue, sdR = sd_red)
        }
        
        norm_values <- process_fluorescence(
            B = data_norm[[blue_channel_normalized]],
            R = data_norm[[red_channel_normalized]],
            maxB = max_neg_blue,
            sdB = scale_parameters$sdB,
            maxR = max_neg_red,
            sdR = scale_parameters$sdR,
            applyNormalization = normalization
        )

        data_norm[[blue_channel_normalized]] <- norm_values$Blue
        data_norm[[red_channel_normalized]] <- norm_values$Red
        
        intensity <- sqrt(norm_values$Blue^2 + norm_values$Red^2)
        angle_rad <- acos(norm_values$Blue / intensity)
        angle_deg <- angle_rad * 180 / pi
        
        result_df <- data.frame(
            file = sample_name,
            Angle = angle_deg,
            Intensity = intensity,
            data_norm,
            norm_values
        )

        result_all <- rbind(result_all, result_df)
        if (verbose) message(paste("Computed angle and intensity for:", sample_name))
    }

    sample_files <- unique(result_all$file)
    
    sampledef_df <- data.frame(
    file = sample_files,
    group = rep("", length(sample_files)),
    stringsAsFactors = FALSE
    )
    sampledef <- list(sampledef = sampledef_df)

    if (normalization) {
        normalization_parameters <- list(
            blue_added = blue_min_adjust,
            red_added = red_min_adjust,
            red_threshold = red_threshold,
            blue_threshold = blue_threshold,
            max_neg_blue = max_neg_blue,
            max_neg_red = max_neg_red,
            SD_blue = sd_blue,
            SD_red = sd_red,
            MAD_blue = MAD_blue,
            MAD_red = MAD_red,
            blue_channel_normalized = blue_channel_normalized,
            red_channel_normalized = red_channel_normalized
        )
    } else {
        normalization_parameters <- list(
            blue_added = blue_min_adjust,
            red_added = red_min_adjust,
            red_threshold = red_threshold,
            blue_threshold = blue_threshold
        )
     }
    
    timer_fluorescence <- list(blue_channel = blue_log_colname, red_channel = red_log_colname)

    cell_counts <- data.frame(
            sample = names(cellcount_total),
            cell_count = cellcount_total,
            stringsAsFactors = FALSE
        )
    
    output <- new("TockyPrepData",
                  Data = result_all,
                  cell_counts = cell_counts,
                  sampledef = sampledef,
                  timer_fluorescence = timer_fluorescence,
                  input = input_list,
                  normalization_parameters = normalization_parameters,
                  Tocky = list())
    
    return(output)
}

#' Update sample definitions and group assignments
#'
#' This function takes the output from `timer_transform`, specifically the `sample_definition` data frame,
#' exports it to a CSV file for the user to edit group assignments, and then reads the updated file back into R.
#'
#' @param x A TockyPrepData object returned by `timer_transform`.
#' @param output_dir Character string specifying the directory to save the `sampledef.csv` file.
#'                   If `NULL`, the file is saved in the current working directory. Default is `NULL`.
#' @param filename Character string specifying the name of the sample definition file. Default is `"sampledef.csv"`.
#' @param sep Character string indicating the field separator used in the CSV file. Default is `","`.
#' @param verbose Logical indicating whether to display messages. Default is `TRUE`.
#' @param sample_definition (Optional) to use a data frame object as an annotation data for sample grouping. Defaul is `NULL`.
#' @param interactive Logical indicating whether to use an interactive session
#'  to export a file for sample grouping and enable user to edit it and import. Defaults to `TRUE`.

#'
#' @return An updated TockyPrepData with user-assigned groupings.
#' @export
#'
#' @examples
#' \dontrun{
#'   # Assuming `x` is the output from `timer_transform`
#'   x <- sampledef(x, output_dir = "output_directory")
#'   # The function will pause, allowing you to edit the 'group' column in the CSV file.
#'   # After editing and saving the file, press Enter in R to continue.
#'   # The updated sample definitions will be returned as a data frame.
#' }
#' @importFrom utils write.table read.table

sample_definition <- function(x, sample_definition = NULL, output_dir = NULL, filename = "sampledef.csv",
                              sep = ",", verbose = TRUE, interactive = FALSE) {
                                  
  if(!inherits(x, "TockyPrepData")){
      stop("Use the output of timer_transform. \n")
  }
  
  if(interactive){
      sn <- x@sampledef$sampledef
      
      if (!is.null(output_dir)) {
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        output_file <- file.path(output_dir, filename)
      } else {
        output_file <- filename
      }
      
      write.table(sn, file = output_file, sep = sep, row.names = FALSE, quote = FALSE)
      
      if (verbose) {
        message("Sample definition file written to: ", output_file)
        message("Please update the 'group' column in the file and save it.")
      }
      
      readline(prompt = "Press [Enter] when you have edited the group clumn updated the sample definition file and saved it.")
      
      sn_updated <- read.table(output_file, sep = sep, header = TRUE, stringsAsFactors = FALSE)
      x@sampledef$sampledef <- sn_updated
      
      if (verbose) {
        message("Updated sample definition has been read.")
      }
  }else{
      if(!is.null(sample_definition)){
          if(is.data.frame(sample_definition)){
              x@sampledef$sampledef <- sample_definition
          }else{
              stop("Use data frame for sampledef data. \n")
          }
      }else{
          stop("Include sampledef data frame or use interactive mode. \n")
      }
      
      
  }
  
  return(x)
}

