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
  
  cat("Final thresholds:\n")
  print(thresholds)
  return(thresholds)
}

# ----------------------------------------------------------
# Read and prepare acceleration data (.parquet file)
# Each file corresponds to one individual and one deployment
# ----------------------------------------------------------

# Define the full path to the main Parquet folder
parquet_path <- "//10.126.19.90/EAS_shared/cross_sleep/working/Data/Acc/spidermonkey/"

# Filename (e.g., "spidermonkey_1_Albus.parquet")
parquet_filename <- "spidermonkey_1_Albus.parquet"

# Read the file
df <- read_parquet(parquet_path)

# Convert to data.table for fast, memory-efficient operations
setDT(df)

# Ensure the Timestamp column is POSIXct (datetime)
if (!inherits(df$Timestamp, "POSIXt")) {
  df$Timestamp <- as.POSIXct(df$Timestamp, tz = "UTC")
}

# Retrieve ID
# Remove ".parquet" manually
filename_core <- sub("\\.parquet$", "", parquet_filename)
# Split by underscores
parts <- strsplit(filename_core, "_")[[1]]
# Assign to variables
species <- parts[1]
deployment_ID <- parts[2]
individual_ID <- parts[3]

# ----------------------------------------------------------
# Step 1: Determine sampling rate (Hz)
# ----------------------------------------------------------
time_diffs <- diff(as.numeric(df$Timestamp))
sampling_rate <- round(1 / median(time_diffs, na.rm = TRUE))
print(sampling_rate)

# ----------------------------------------------------------
# Step 2: Check how many rows fall exactly on the second
# ----------------------------------------------------------
is_exact_second <- df$Timestamp == floor_date(df$Timestamp, unit = "second")
cat("Samples on exact seconds: ", sum(is_exact_second), "\n")
cat("Unique seconds: ", length(unique(floor_date(df$Timestamp, unit = "second"))), "\n")

# ----------------------------------------------------------
# Step 3: Compute 1-second VeDBA (no post-smoothing)
# ----------------------------------------------------------

# Axes to process
present_axes <- c("X", "Y", "Z")
rolling_mean_width <- sampling_rate  # For example, 20 for 20 Hz

# (a) Compute static acceleration columns:
#     For each axis (X, Y, etc.), we do a rolling mean with width = rolling_mean_width = 1s.
df[,(paste0(present_axes, "_static")) := lapply(present_axes, function(axis) {
    frollmean(
      get(axis),               # <-- the raw acceleration column "X", "Y", etc.
      n = rolling_mean_width,  
      align = "center",
      na.rm = FALSE)
  })]

# (b) Compute dynamic acceleration columns:
#     Subtract the static from the raw signal, axis by axis.
df[,(paste0(present_axes, "_dynamic")) := Map(
    function(raw, stc) raw - stc,
    lapply(present_axes, function(axis) get(axis)),
    .SD), # .SD refers to the newly created static columns in the same row subset
  .SDcols = paste0(present_axes, "_static")]

# (c) Compute VeDBA (based on dynamic columns just created):
#     VeDBA = sqrt(X_dyn^2 + Y_dyn^2 + Z_dyn^2).
df[,VeDBA := sqrt(Reduce(`+`, lapply(.SD, function(d) d^2)) # sums them, producing a single vector of “X^2 + Y^2 + Z^2”
  ),.SDcols = paste0(present_axes, "_dynamic")]

# (d) Remove temporary static and dynamic columns
df <- df[, !c(paste0(present_axes, "_static"), paste0(present_axes, "_dynamic")), with = FALSE]

# ----------------------------------------------------------
# Step 4: Downsample to 1 Hz (based on exact-second timestamps)
# ----------------------------------------------------------
df_1s <- df[df$Timestamp == floor_date(df$Timestamp, unit = "second")]
cat("Downsampling ratio (raw/1s): ", nrow(df) / nrow(df_1s), "\n")

# Can also simply use
# df_1s <- downsample_data_table(df, original_rate = sampling_rate, target_rate = 1)

# Log-transform the VeDBA
df_1s[, log_VeDBA := log(VeDBA)]

# ----------------------------------------------------------
# Step 5: Interactive threshold selection on histogram
# ----------------------------------------------------------
thresholds <- select_thresholds_on_hist(df_1s$log_VeDBA)
print(thresholds)

thresholds <- as.data.frame(as.list(thresholds))
thresholds <- thresholds %>%
  mutate(
    species = species,
    deployment_ID = deployment_ID,
    individual_ID = individual_ID
  )
print(thresholds)
