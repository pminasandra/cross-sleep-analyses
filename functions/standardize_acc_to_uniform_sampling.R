#' Convert ACC data that may have gaps to uniformly sampled ACC data with NAs
#' where gaps exist
#'
#' @param df data frame containing columns Timestamp, X, Y, Z
#' 
#' @returns df_standard: data frame in the same format but with additional timestamps where needed with NAs for X,Y,Z

standardize_acc_to_uniform_sampling <- function(df){
  
  library(arrow)
  library(data.table)
  library(lubridate)
  library(dplyr)
  
  setDT(df)
  
  # Ensure Timestamp is POSIXct
  if (!inherits(df$Timestamp, "POSIXt")) {
    df$Timestamp <- as.POSIXct(df$Timestamp, tz = "UTC")
  }
  
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
  
  #put back in df
  df <- df_filled
  
  out <- list()
  out$df <- df
  out$sampling_rate <- sampling_rate
  
  return(out)
}
