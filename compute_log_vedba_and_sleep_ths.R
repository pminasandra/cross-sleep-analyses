# This script processes pre-prepared triaxial acceleration data with columns:
# Timestamp, X, Y, and Z. It computes Vectorial Dynamic Body Acceleration (VeDBA),
# log-transforms it, downsamples it to 1 Hz, and provides an interactive function
# to visually select thresholds from a histogram using mouse clicks.

library(data.table)  # For fast, efficient, in-place data manipulation (:=, .SD, .SDcols)
library(lubridate)   # For handling date-time values
library(stringr)     # For string manipulation (not currently used here)
library(arrow)       # For reading .parquet files
library(dplyr)

# ----------------------------------------------------------
# Custom function by Rich Gunner
# ----------------------------------------------------------
# Downsample a data.table to a lower frequency by selecting one row every 
# (original_rate / target_rate) rows. Ensures the resulting data keeps 
# consistent spacing. Only works if target_rate <= original_rate.
downsample_data_table <- function(data, original_rate, target_rate) {
  if (length(original_rate) != 1) {
    stop("Multiple sampling frequencies found. Please check data integrity.")
  }
  if (target_rate > original_rate) {
    stop("target_rate must be <= the original sampling frequency.")
  }
  
  factor <- original_rate / target_rate
  indices <- seq(1, nrow(data), by = factor)
  rounded_indices <- round(indices)
  valid_indices <- rounded_indices[rounded_indices > 0 & rounded_indices <= nrow(data)]
  
  # Return the downsampled rows
  data[valid_indices, ]
}

# ----------------------------------------------------------
# Interactive function to define thresholds via mouse clicks
# ----------------------------------------------------------
# This function opens a histogram of a numeric vector and lets the user click 
# once or twice to set thresholds (e.g., for activity classification).
# - A red vertical line is drawn for the first click (th1)
# - A blue vertical line is drawn for the second click (th2)
# - Clicking outside the x-axis range or pressing Esc results in NA
select_thresholds_on_hist <- function(x, breaks = 100) {
  if (!is.numeric(x)) stop("x must be a numeric vector")
  
  dev.new()  # Open a new window for better interactivity
  
  # Plot histogram and determine x-axis range for valid clicks
  hist_data <- hist(x, breaks = breaks, col = "lightgray", border = "white",
                    main = "Click up to 2 thresholds\n(Esc or outside = skip)",
                    xlab = "Value")
  x_range <- range(hist_data$breaks)
  thresholds <- c(th1 = NA, th2 = NA)
  
  for (i in 1:2) {
    cat(sprintf("Click to set threshold th%d (Esc or click outside to skip)\n", i))
    click <- try(locator(1), silent = TRUE)
    
    # If user presses Esc or doesn't click
    if (inherits(click, "try-error") || length(click$x) == 0) {
      cat(sprintf("Threshold th%d skipped (Esc pressed)\n", i))
      next
    }
    
    xval <- click$x
    if (xval < x_range[1] || xval > x_range[2]) {
      cat(sprintf("Threshold th%d ignored (%.3f out of range)\n", i, xval))
      next
    }
    
    thresholds[i] <- xval
    abline(v = xval, col = ifelse(i == 1, "red", "blue"), lwd = 2, lty = 2)
    cat(sprintf("Threshold th%d = %.3f\n", i, xval))
  }
  
  dev.off()
  cat("Final thresholds:\n")
  print(thresholds)
  return(thresholds)
}

# --------------------------------------------------------------------
#' Round POSIXct timestamps to the nearest sampling interval
#'
#' @param timestamps A vector of POSIXct timestamps
#' @param sampling_rate Numeric, in Hz (e.g., 20 means 20 samples per second)
#' @return A vector of rounded POSIXct timestamps
#' @author Your Name
round_to_nearest_sample <- function(timestamps, sampling_rate) {
  # Compute the interval length in seconds (e.g., 0.05 for 20 Hz)
  interval <- 1 / sampling_rate
  
  # Convert to numeric (seconds), round, then convert back to POSIXct
  rounded <- as.POSIXct(
    round(as.numeric(timestamps) / interval) * interval,
    origin = "1970-01-01", tz = "UTC"
  )
  
  return(rounded)
}


# ----------------------------------------------------------
# Read and prepare acceleration data (.parquet file)
# Each file corresponds to one individual and one deployment
# ----------------------------------------------------------

# Define the full path to the main Parquet folder
parquet_directory <- "//10.126.19.90/EAS_shared/cross_sleep/working/Data/Acc/"
metadata_directory <- "//10.126.19.90/EAS_shared/cross_sleep/working/Data/cross_sleep_metadata.csv"


file_list <- list.files(
  path = parquet_directory,
  pattern = "\\.parquet$",   # match files ending in .parquet
  recursive = TRUE,
  full.names = TRUE
)


metadata <- read.csv(metadata_directory,colClasses = "character")

# Filenames (e.g., "spidermonkey_1_Albus.parquet")

for (i in 1:length(file_list)){
  print(i)
  print(paste0("Processing ", file_list[i]))
  result <- compute_vedba_1s_and_select_sleep_threshold(file_list[i])
  
  # Check if the row exists in metadata and then add the selected th
  exists <- any(
    metadata$species == result$thresholds$species & 
      metadata$individual_ID == result$thresholds$individual_ID &
      metadata$deployment_ID == result$thresholds$deployment_ID
  )
  
  if (exists) {
    metadata <- metadata %>%
      mutate(across(everything(), ~ .)) %>%  # noop, keeps metadata a tibble
      rows_update(result$thresholds, by = c("species", "deployment_ID", "individual_ID"))
  } else {
    metadata <- bind_rows(metadata, result$thresholds)
  }
  

}
print(metadata)
#write the updated metadata

# Today's date
{date_suffix <- format(Sys.Date(), "%Y-%m-%d")

# New path with date inserted
new_metadata_path <- sub(
  pattern = "\\.csv$",
  replacement = paste0("_", date_suffix, ".csv"),
  x = metadata_directory
)

# Prompt user for confirmation
cat("Are you sure you want to write metadata to:\n", new_metadata_path, "\n"); response <- readline(prompt = "Type 'yes' to confirm: ")

if (tolower(response) == "yes") {
  write.csv(metadata, new_metadata_path, row.names = FALSE)
  cat("✅ Metadata successfully written to:\n", new_metadata_path, "\n")
} else {
  cat("❌ Operation cancelled. Metadata not written.\n")
}}
