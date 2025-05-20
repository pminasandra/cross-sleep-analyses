#' Compute 1s VeDBA and select sleep detection threshold
#'
#' This function reads a .parquet accelerometer file, computes VeDBA at 1Hz,
#' fills in any missing seconds, and allows interactive threshold selection
#' on log-transformed VeDBA values for sleep detection purposes.
#'
#' @param parquet_path Path to the input .parquet file.
#' @return A list with:
#'   - df_1s: the downsampled 1Hz data table with complete timestamps
#'   - thresholds: a data.frame of selected thresholds with metadata
#' @author Your Name
compute_vedba_1s_and_select_sleep_threshold <- function(parquet_path) {
  # Load required libraries
  library(arrow)
  library(data.table)
  library(lubridate)
  library(dplyr)
  
  # Read the file
  df <- read_parquet(parquet_path)
  setDT(df)
  
  # Ensure Timestamp is POSIXct
  if (!inherits(df$Timestamp, "POSIXt")) {
    df$Timestamp <- as.POSIXct(df$Timestamp, tz = "UTC")
  }
  
  # Extract metadata from filename
  parquet_filename <- basename(parquet_path)
  filename_core <- sub("\\.parquet$", "", parquet_filename)
  parts <- strsplit(filename_core, "_")[[1]]
  species <- parts[1]
  deployment_ID <- parts[2]
  individual_ID <- parts[3]
  
  # Step 1: Estimate sampling rate (Hz)
  time_diffs <- diff(as.numeric(df$Timestamp))
  sampling_rate <- round(1 / median(time_diffs, na.rm = TRUE))
  cat("Estimated sampling rate (Hz):", sampling_rate, "\n")
  
  #Adjust the timestamp so that each sample always fall on the nearest 1/sampling rate timepoint
  df$Timestamp_adj <- round_to_nearest_sample(df$Timestamp, sampling_rate)
  
  
  # Generate the full sequence of timestamps from start to end
  # This will be used to make sure that the table is regularly sampled (NA if no data)
  full_timestamps <- seq(
    from = as.numeric(min(df$Timestamp_adj)),
    to   = as.numeric(max(df$Timestamp_adj)),
    by   = 1/sampling_rate  # 20 Hz, i.e., every 50 ms
  )
  full_timestamps <- as.POSIXct(full_timestamps, origin = "1970-01-01", tz = "UTC")
  full_timestamps <- round_to_nearest_sample(full_timestamps, sampling_rate)
  format(full_timestamps[1:10], "%Y-%m-%d %H:%M:%OS3")
  
  
  # Identify missing timestamps
  missing_seconds <- setdiff(full_timestamps, df$Timestamp_adj)
  missing_seconds <- as.POSIXct(missing_seconds, origin = "1970-01-01", tz = "UTC")
  
  # If there are missing timestamps, fill them in with NA
  if (length(missing_seconds) > 0) {
    # Create a data frame with missing timestamps and NA in other columns
    missing_df <- data.frame(Timestamp_adj = missing_seconds)
    other_cols <- setdiff(names(df), "Timestamp_adj")
    for (col in other_cols) {
      missing_df[[col]] <- NA
    }
    
    # Combine original and missing rows, then sort
    missing_df <- missing_df[, names(df)]  # Reorder columns to match df
    missing_df$Timestamp <- as.POSIXct(missing_df$Timestamp, tz = "UTC")
    
    df_filled <- rbind(df, missing_df)
    df_filled <- df_filled[order(df_filled$Timestamp_adj), ]
  } else {
    df_filled <- df
  }
  
  #just check original timestamps vs. adjusted
  ii = 3
  format(df_filled$Timestamp[ii], "%Y-%m-%d %H:%M:%OS3");format(df_filled$Timestamp_adj[ii], "%Y-%m-%d %H:%M:%OS3");
  
  #put back in df if good
  df <- df_filled
  
  # Step 3: Compute VeDBA
  present_axes <- c("X", "Y", "Z")
  rolling_mean_width <- sampling_rate
  
  df[, (paste0(present_axes, "_static")) := lapply(present_axes, function(axis) {
    frollmean(get(axis), n = rolling_mean_width, align = "center", na.rm = FALSE)
  })]
  
  df[, (paste0(present_axes, "_dynamic")) := Map(
    function(raw, stc) raw - stc,
    lapply(present_axes, function(axis) get(axis)),
    .SD
  ), .SDcols = paste0(present_axes, "_static")]
  
  df[, VeDBA := sqrt(Reduce(`+`, lapply(.SD, function(d) d^2))),
     .SDcols = paste0(present_axes, "_dynamic")]
  
  df <- df[, !c(paste0(present_axes, "_static"), paste0(present_axes, "_dynamic")), with = FALSE]
  
  # Step 4: Downsample to 1Hz
  df_1s <- df[df$Timestamp_adj == floor_date(df$Timestamp_adj, unit = "second")]
  cat("Downsampling ratio (raw/1s):", nrow(df) / nrow(df_1s), "\n")
  
  # Log-transform VeDBA
  df_1s[, log_VeDBA := log(VeDBA)]
  
  #View(df_1s)
  ii = 3
  format(df_1s$Timestamp[ii], "%Y-%m-%d %H:%M:%OS3");format(df_1s$Timestamp_adj[ii], "%Y-%m-%d %H:%M:%OS3");
  
  
  # Step 6: Select thresholds interactively
  thresholds <- select_thresholds_on_hist(df_1s$log_VeDBA)
  #print(thresholds)
  
  thresholds <- as.data.frame(as.list(thresholds)) %>%
    mutate(
      species = species,
      deployment_ID = deployment_ID,
      individual_ID = individual_ID
    )
  
  return(list(
    df_1s = df_1s,
    thresholds = thresholds
  ))
}

