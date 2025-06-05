library(arrow)
library(dplyr)
library(stringr)
library(zoo)
library(readr)      # if needed for reading data
library(mixtools)   # for normalmixEM
library(hms)
library(lubridate)
library(purrr)

## functions

source('functions/calculate_threshold.R')
source('functions/classify_sleep_windows.R')
source('functions/calculate_sleep_metrics.R')
source('functions/plot_interactive_spt.R')

# Define parameter grids; adjust as needed
mov_window_grid <- c(9, 25) # e.g., 541 and 1501 for spider monkey and hyena (in minutes*2+1 for odd window)
block_size_grid <- c(15, 30, 45)    # e.g., 1800, 3600 (in seconds)
gap_size_grid <- c(15, 30, 45)                # e.g., 15 (in minutes)
waso_block_grid <- c(3*60)            # e.g., 180 (in seconds)
frag_block_grid <- c(2*60)            # e.g., 120 (in seconds)

# Paths
input_dir <- "/mnt/EAS_shared/cross_sleep/working/Data/VeDBA_cont/"
output_classified_dir <- "/mnt/EAS_shared/cross_sleep/working/Data/ClassifiedSleep/"
output_metrics_dir <- "/mnt/EAS_shared/cross_sleep/working/Data/SleepMetrics"
output_plot_dir <- "/mnt/EAS_shared/cross_sleep/working/Figures/SPT_sample_check"


# Get list of .parquet files
parquet_files <- list.files(input_dir, pattern = "\\.parquet$", full.names = TRUE)

#sleep_metrics_all <- list()

# Grid over all parameter combinations
param_grid <- expand.grid(
  mov_window = mov_window_grid,
  block_size = block_size_grid,
  gap_size = gap_size_grid,
  waso_block = waso_block_grid,
  frag_block = frag_block_grid,
  stringsAsFactors = FALSE
)

sample_files = c(3, 5, 8, 14, 17, 21, 29, 30)  # parquet_files
for (parquet_path in parquet_files[sample_files]) {
  try({
    
    sleep_metrics_tag <- list()
    
    parquet_filename <- basename(parquet_path)
    filename_core <- sub("\\.parquet$", "", parquet_filename)
    parts <- strsplit(filename_core, "_")[[1]]
    species <- parts[1]
    deployment_ID <- parts[2]
    individual_ID <- parts[3]
    tag_name <- paste(species, deployment_ID, individual_ID, sep = "_")
    
    # Read data
    d1 <- read_parquet(parquet_path) |> as.data.frame()
    
    # Time prep
    time_factor <- 12*60*60
    time_shift <- d1$Timestamp - time_factor
    start_date <- as.Date(min(d1$Timestamp) - time_factor)
    
    d1$night <- as.numeric(as.Date(time_shift) - start_date + 1)
    d1$night_date <- as.Date(d1$Timestamp - time_factor)
    d1$Time <- format(d1$Timestamp, format = "%H:%M:%S")
    
    # Threshold + classification
    min_vedba <- min(d1$log_vedba[is.finite(d1$log_vedba)], na.rm = TRUE)
    d1$log_vedba[!is.finite(d1$log_vedba)] <- min_vedba
    thresh <- calculate_threshold(d1$log_vedba)
    
    # Loop over parameter grid
    for (i in seq_len(nrow(param_grid))) {
      mov_window <- param_grid$mov_window[i]
      block_size <- param_grid$block_size[i]
      gap_size <- param_grid$gap_size[i]
      waso_block <- param_grid$waso_block[i]
      frag_block <- param_grid$frag_block[i]
      
      # You may need to pass these to your functions explicitly if they are not global
      classified_sleep <- classify_sleep_windows(
        d1, tag_name, thresh,
        mov_window = mov_window,
        block_size = block_size,
        waso_block = waso_block,
        gap_size = gap_size
      )
      
      # Save classified sleep
      out_file_classified <- file.path(
        output_classified_dir,
        sprintf("%s_classified_sleep_MW%d_BS%d_GS%d_WB%d_FB%d.parquet",
                tag_name, mov_window, block_size, gap_size, waso_block, frag_block
        )
      )
      write_parquet(classified_sleep, out_file_classified)
      
      # Sleep metrics
      sleep_metrics <- calculate_sleep_metrics(classified_sleep, tag_name, frag_block = frag_block)
      # Add metadata for traceability
      sleep_metrics$mov_window <- mov_window
      sleep_metrics$block_size <- block_size
      sleep_metrics$gap_size <- gap_size
      sleep_metrics$waso_block <- waso_block
      sleep_metrics$frag_block <- frag_block
      sleep_metrics$input_file <- parquet_filename
      
      sleep_metrics_tag[[paste(tag_name, mov_window, block_size, gap_size, waso_block, frag_block, sep = "_")]] <- sleep_metrics
    }
  }, silent = TRUE)

  # Combine all metrics into one long data frame
  sleep_metrics_tag <- bind_rows(sleep_metrics_tag, .id = "tag_param")

  out_file_plot <- file.path(
    output_classified_dir,
    sprintf("%s_classified_sleep_MW%d_BS%d_GS%d_WB%d_FB%d.html",
            tag_name, mov_window, block_size, gap_size, waso_block, frag_block
    ))
  interactive_spt_plot(classified_sleep, sleep_metrics_tag, night_id = 2, output_plot_dir, out_file_plot)
  
}
