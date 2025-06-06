#' @title Compare Automatic and Manual Thresholds for VeDBA Data
#'
#' @description
#' This function analyzes log-transformed VeDBA values from `.parquet` files, computing automatic behavioral thresholds using a Gaussian Mixture Model (GMM),
#' and compares them against manually defined thresholds provided in a metadata table.
#' The function processes all files from the same individual/species, summarizes threshold differences, and saves visual diagnostics for quality control.
#'
#' @param parquet_path Character. Path to a `.parquet` file containing a `logvedba` column. This file should ideally be the baseline continuous VeDBA dataset
#' (i.e., the one without `burst_int` or `burst_len` in the filename).
#' @param k Integer. Number of components (clusters) for the Gaussian Mixture Model. Defaults to 2.
#' @param thres_value A numeric value indicating the minimum posterior probability cutoff to define a threshold. Default is 0.1.
#' @param metadata Data frame. Should contain columns `deployment_ID`, `individual_ID`, and the manually set threshold (`th1`).
#' @param fig_output_folder Character. Path to a directory where diagnostic plots will be saved.
#'
#' @return A data frame with automatic and manual thresholds per file, including:
#' \describe{
#'   \item{species}{Extracted species name.}
#'   \item{deployment_ID}{Deployment ID from file name.}
#'   \item{individual_ID}{Individual ID from file name.}
#'   \item{burst_int}{Burst interval extracted from file name, if present.}
#'   \item{burst_length}{Burst length extracted from file name, if present.}
#'   \item{manual_th}{Manual threshold from metadata.}
#'   \item{automatic_th}{Primary threshold from GMM.}
#'   \item{automatic_th2}{Secondary threshold from GMM, if `k > 2`.}
#'   \item{threshold_diff}{Difference between automatic and manual thresholds.}
#' }
#'
#' The function also produces:
#' \itemize{
#'   \item Histogram plots with posterior distributions and threshold lines.
#'   \item Correlation plots between burst interval and threshold difference.
#'   \item Linearity checks between burst regime and threshold.
#' }
#'
#' @note This function assumes that filenames follow a consistent naming pattern including species, deployment ID, individual ID, and optionally burst information.
#'
#' @examples
#' \dontrun{
#' metadata <- read.csv("cross_sleep_metadata.csv")
#' compare_auto_manual_thresholds(
#'   parquet_path = "VeDBA/hyena_2016_BORA_vedba_standard.parquet",
#'   k = 2,
#'   metadata = metadata,
#'   fig_output_folder = "Figures/thresholds/"
#' )
#' }
#'
#' @import arrow ggplot2 stringr dplyr cowplot
#' @author Cross-sleep team, MPIAB
#' @export
compare_auto_manual_thresholds <- function(parquet_path, thres_value = 0.1, k = 2, metadata, fig_output_folder) {
  # Load required libraries
  library(arrow)
  library(ggplot2)
  library(stringr)
  library(dplyr)
  library(cowplot)
  
  # Extract file name and folder path
  parquet_filename <- basename(parquet_path)         # Get just the filename from full path
  vedbadata_path <- dirname(parquet_path)            # Get directory path
  filename_core <- sub("\\.parquet$", "", parquet_filename)  # Remove extension
  parts <- strsplit(filename_core, "_")[[1]]         # Split name by underscores
  
  # Extract metadata from filename: species, deployment, individual
  curr_species <- parts[1]
  curr_deployment_ID <- parts[2]
  curr_individual_ID <- parts[3]
  
  # List all files in the same folder that match species and individual
  all_files <- list.files(vedbadata_path, full.names = TRUE)
  matching_files <- all_files[
    grepl(curr_individual_ID, all_files) & grepl(curr_species, all_files)
  ]
  
  # Extract burst intervals from file names and sort files accordingly
  burst_ints <- as.numeric(str_extract(matching_files, "(?<=burst_int)\\d+"))
  non_matching <- is.na(burst_ints)
  matching_files <- c(
    matching_files[non_matching],  # Put baseline file(s) first
    matching_files[!non_matching][order(burst_ints[!non_matching])]  # Sort burst files numerically
  )
  
  # Prepare an empty dataframe to store thresholds and file metadata
  n_rows <- length(matching_files)
  Auto_th <- data.frame(
    species       = rep(NA_character_, n_rows),
    deployment_ID = rep(NA_character_, n_rows),
    individual_ID = rep(NA_character_, n_rows),
    burst_int     = rep(NA_real_, n_rows),
    burst_length  = rep(NA_real_, n_rows),
    manual_th     = rep(NA_real_, n_rows),
    automatic_th  = rep(NA_real_, n_rows),
    automatic_th2 = rep(NA_real_, n_rows),
    stringsAsFactors = FALSE
  )
  
  # Read baseline file to determine plotting limits
  d1 <- read_parquet(parquet_path)
  x_limits <- c(floor(min(d1$logvedba, na.rm = TRUE) - 1),
                ceiling(max(d1$logvedba, na.rm = TRUE) + 1))  # Extend range for cleaner plots
  
  plot <- list()  # Container for distribution plots
  
  for (i in seq_along(matching_files)) {
    # Read each file and extract relevant info
    print(paste0("Working on ", matching_files[i]))
    d1 <- read_parquet(matching_files[i])
    parquet_filename <- basename(matching_files[i])
    filename_core <- sub("\\.parquet$", "", parquet_filename)
    parts <- strsplit(filename_core, "_")[[1]]
    curr_species <- parts[1]
    curr_deployment_ID <- parts[2]
    curr_individual_ID <- parts[3]
    curr_burst_int <- as.numeric(str_extract(parquet_filename, "(?<=burst_int)\\d+"))
    curr_burst_len <- as.numeric(str_extract(parquet_filename, "(?<=_len)\\d+"))
    
    # Retrieve manual threshold from metadata
    manual_threshold <- metadata %>%
      filter(deployment_ID == curr_deployment_ID & individual_ID == curr_individual_ID) %>%
      pull(th1)
    
    # Run threshold estimation, retry if infinite value returned
    attempt <- 1
    print(paste0("#", attempt, " attempt"))
    results <- calculate_threshold(d1, thres_value, k)
    while (is.infinite(results$threshold[1])) {
      attempt <- attempt + 1
      print(paste0("#", attempt, " attempt"))
      results <- calculate_threshold(d1, thres_value, k)
    }
    
    # Save metadata and thresholds
    Auto_th$species[i] <- curr_species
    Auto_th$deployment_ID[i] <- curr_deployment_ID
    Auto_th$individual_ID[i] <- curr_individual_ID
    Auto_th$burst_int[i] <- curr_burst_int
    Auto_th$burst_length[i] <- curr_burst_len
    Auto_th$manual_th[i] <- manual_threshold
    Auto_th$automatic_th[i] <- results$threshold[1]
    Auto_th$automatic_th2[i] <- if (length(results$threshold) > 1) results$threshold[2] else NA
    
    # Create a plot of posterior distribution and thresholds
    plot[[i]] <- ggplot(results$vedba_table_posterior, aes(x = logvedba)) +
      geom_histogram(aes(y = ..density..), bins = 100, fill = "blue", alpha = 0.5) +
      geom_point(aes(y = posterior_comp1, color = "Component 1"), size = 1, alpha = 0.7) +
      labs(title = filename_core, x = "log(VeDBA)", y = "Density / Probability") +
      xlim(x_limits) +
      geom_vline(xintercept = results$threshold, linetype = "dotted", color = "red", size = 2) +
      geom_vline(xintercept = manual_threshold, linetype = "dotted", color = "blue", size = 2) +
      theme_minimal() +
      theme(legend.position = "none")
    num_components <- sum(grepl("^posterior_comp", names(results$vedba_table_posterior)))
    if (num_components > 2) {
      for (ii in 2:num_components-1) {
        plot[[i]] <- plot[[i]] + geom_point(aes_string(y = paste0("posterior_comp", ii), color = shQuote(paste0("Component ", ii))), size = 1, alpha = 0.7)
      }
    }
  }
  
  # Combine individual distribution plots into a single grid
  grid_plot <- do.call(plot_grid, c(plot[1:length(plot)], ncol = 3, labels = "AUTO", align = "hv"))
  ggsave(paste0(fig_output_folder, curr_species, "_", curr_deployment_ID, "_", curr_individual_ID, "_vedba_disttrib.png"), 
         plot = grid_plot, width = 24, height = 16, dpi = 100, bg = "white")
  
  # Compute and store threshold difference
  Auto_th$threshold_diff <- Auto_th$automatic_th - Auto_th$manual_th
  
  # Coefficient of variation (CV) of automatic thresholds
  cv_percent <- (sd(Auto_th$automatic_th, na.rm = TRUE) / abs(mean(Auto_th$automatic_th, na.rm = TRUE))) * 100
  cv_label <- paste0("CV automatic ths = ", round(cv_percent, 2), "%")
  
  # Convert log-thresholds back to linear scale and compute % difference
  manual_linear <- exp(Auto_th$manual_th)
  auto_linear <- exp(Auto_th$automatic_th)
  percent_diff <- (auto_linear - manual_linear) / manual_linear * 100
  mean_percent_diff <- mean(percent_diff, na.rm = TRUE)
  
  # Plot difference between auto/manual threshold as a function of burst interval
  diff_plot <- ggplot(Auto_th, aes(x = burst_int, y = threshold_diff)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    labs(title = "Correlation between Burst Interval and Threshold Difference",
         x = "Burst Interval", y = "Automatic - Manual Threshold") +
    annotate("text", x = -Inf, y = Inf, label = cv_label, hjust = -0, vjust = 1.5, size = 4) +
    theme_minimal() +
    geom_hline(yintercept = Auto_th$automatic_th[is.na(Auto_th$burst_int)] - Auto_th$manual_th[1],
               linetype = "dotted", color = "red", size = 2)
  
  ggsave(paste0(fig_output_folder, Auto_th$species[1], "_", 
                Auto_th$deployment_ID[1], "_", Auto_th$individual_ID[1], "_th-int.png"), 
         plot = diff_plot, width = 8, height = 5, dpi = 100, bg = "white")
  
  # Plot automatic thresholds over burst intervals
  threshold_plot <- ggplot(Auto_th, aes(x = burst_int, y = automatic_th)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    labs(x = "Burst Interval", y = "Thresholds") +
    geom_hline(yintercept = Auto_th$manual_th[1], linetype = "dotted", color = "blue", size = 2) +
    geom_hline(yintercept = Auto_th$automatic_th[is.na(Auto_th$burst_int)], linetype = "dotted", color = "red", size = 2) +
    theme_minimal() + ylim(x_limits)
  
  ggsave(paste0(fig_output_folder, Auto_th$species[1], "_", 
                Auto_th$deployment_ID[1], "_", Auto_th$individual_ID[1], "_th-int2.png"), 
         plot = threshold_plot, width = 8, height = 5, dpi = 100, bg = "white")
  
  return(Auto_th)
}

