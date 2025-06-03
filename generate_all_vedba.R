############ Script for generating burst sampled and continuous VeDBA #########

############# Functions #############

# standardize_acc_to_uniform_sampling = Estimating a single sampling rate from the median timestamp difference.
# = Rounding each row’s timestamp (Timestamp_adj) to that sampling grid.
# = Creating a full grid from the earliest to the latest adjusted timestamp.
# = Filling in missing timestamps with rows of NA in the acceleration columns.

# burst_sample = Simulates burst sampling of acc data at a given interval burst_interval and for a given burst_length (in seconds)
# round_to_nearest_sample = Helper function to round timestamps to the nearest sample
# acc_to_vedba = Compute rolling VeDBA (optiionally via grouping)

#Params
rolling_width = 20
fracNA_thresh = .5
outdir <- '~/EAS_shared/cross_sleep/working/Data/'
codedir <- '~/EAS_ind/astrandburg/code/'


# Required libraries
library(arrow)
library(data.table)
library(lubridate)
library(dplyr)

# Source required functions
source(paste0(codedir, 'vedba-general/functions/acc_to_vedba.R'))
source(paste0(codedir, 'vedba-general/functions/round_to_nearest_sample.R'))
source(paste0(codedir, 'vedba-general/functions/standardize_acc_to_uniform_sampling.R'))
source(paste0(codedir, 'burst-sampling-regimes/burst_sample.R'))

#specify default time zone as UTC
Sys.setenv(TZ='UTC')

#set working directory
setwd(outdir)

#Ensure decimal secs are displayed
options(digits.secs = 3) 

# Spider monkey - Albus
df <- read_parquet("Acc/spidermonkey/spidermonkey_1_Albus.parquet")

df <- as.data.table(df) # Make data.table
df[, Timestamp := as.POSIXct(Timestamp, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")] # Ensure UTC for consistency 

# 1)  Standardize acc to uniform_sampling - standardizing the time series to a uniform sampling grid with NA for missing records
df_stand = standardize_acc_to_uniform_sampling(df)

# 2) Clean data frame - make adjusted (rounded) time stamp the main one we use
df_stand$Timestamp = df_stand$Timestamp_adj ; df_stand$Timestamp_adj = NULL

# Write standarsised data as parquet file
write_parquet(df_stand, "Acc_standard/spidermonkey_1_Albus_standard.parquet")

# 3) Raw VeDBA all standardised data 
df_Ved = acc_to_vedba(df_stand, rolling_mean_width = rolling_width, group_col = NULL)

# 3b) Optional - remove VeDBA values if fracNA >= 50 % #Arbitrary!
# Identify non-POSIXct columns (i.e., not Timestamp)
non_time_cols <- names(df_Ved)[!sapply(df_Ved, inherits, what = "POSIXct")]
non_time_cols = subset(non_time_cols, non_time_cols != "fracNA") # Keep value in this column to aid diagnostics
# Set values to NA where fracNA >= 0.5
df_Ved[fracNA >= fracNA_thresh, (non_time_cols) := lapply(.SD, function(x) NA), .SDcols = non_time_cols]

# 4) Log-transform VeDBA
df_Ved[, log_VeDBA := log(VeDBA)]

# 5) Median VeDBA value per second
# Create a column that rounds down to the nearest second
df_Ved[, Timestamp_sec := as.POSIXct(floor(as.numeric(Timestamp)), origin = "1970-01-01", tz = attr(Timestamp, "tzone"))]
# Now aggregate by second
df_Ved_1Hz <- df_Ved[
  ,
  .(
    VeDBA = median(VeDBA, na.rm = TRUE),
    logvedba = median(log_VeDBA, na.rm = TRUE)
  ),
  by = Timestamp_sec
]

# 6) Export median VeDBA per second - continous file
colnames(df_Ved_1Hz) <- c("Timestamp", "vedba", "log_vedba")
write_parquet(df_Ved_1Hz, "VeDBA_cont/spidermonkey_1_Albus_VeDBA_standard.parquet")

#####################################################################################################################################

# 7) Pseudo undersample (burst lengths and inter-intervals)
# a) 2 s burst length, 1 min burst interval

#Params
burst_interval <- 60 #burst interval in seconds 
burst_length <- 2 #burst length in seconds
acc_sample_rate <- 20 #samples per second in the raw ACC data
start_timestamp <- as.POSIXct('2024-05-04 11:00:00.00', tz = 'UTC')
filename <- 'Acc_standard/spidermonkey/spidermonkey_1_Albus_standard.parquet'
savename <- paste0('VeDBA_burst/spidermonkey_1_Albus_VeDBA_burst_int',burst_interval,'_len',burst_length,'.parquet')

run_and_save_burst_data <- function(filename, savename, burst_interval, burst_length, acc_sample_rate, start_timestamp){

  df_stand = read_parquet(filename)
  
  # Burst IDs
  df_Ved_2s_60s <- burst_sample(accdat = as.data.frame(df_stand), burst_interval = burst_interval, burst_length = burst_length, acc_sample_rate = acc_sample_rate, start_timestamp = start_timestamp)
  # Remove inter-burst interval data
  df_Ved_2s_60s = subset(df_Ved_2s_60s, is.na(df_Ved_2s_60s$Burst_ID) != TRUE)
  # Compute VeDBA per burst ID 
  df_Ved_2s_60s <- as.data.table(df_Ved_2s_60s) # Make data.table
  df_Ved_2s_60s = acc_to_vedba(df_Ved_2s_60s, rolling_mean_width = rolling_width, group_col = "Burst_ID")
  
  # Optional - remove VeDBA values if fracNA >= 50 % #Arbitrary!
  # Identify non-POSIXct columns (i.e., not Timestamp)
  non_time_cols <- names(df_Ved_2s_60s)[!sapply(df_Ved_2s_60s, inherits, what = "POSIXct")]
  non_time_cols = subset(non_time_cols, non_time_cols != "fracNA") # Keep value in this column to aid diagnostics
  # Set values to NA where fracNA >= 0.5
  df_Ved_2s_60s[fracNA >= fracNA_thresh, (non_time_cols) := lapply(.SD, function(x) NA), .SDcols = non_time_cols]
  
  # Log-transform VeDBA
  df_Ved_2s_60s[, log_VeDBA := log(VeDBA)]
  
  # Create a column that rounds down to the nearest second
  df_Ved_2s_60s[, Timestamp_sec := as.POSIXct(floor(as.numeric(Timestamp)), origin = "1970-01-01", tz = attr(Timestamp, "tzone"))]
  
  # Mean per ID
  df_Ved_2s_60s_seq <- df_Ved_2s_60s[
    ,
    .(
      Timestamp = first(Timestamp_sec),  # retain first timestamp
      VeDBA = mean(VeDBA, na.rm = TRUE),
      logvedba = mean(log_VeDBA, na.rm = TRUE)
    ),
    by = Burst_ID
  ]

  # Diagnostics
  #plot(df_Ved_2s_60s_seq$Timestamp[-1], cumsum(as.numeric(diff(df_Ved_2s_60s_seq$Timestamp))), type = "l", xlab = "Time (s)", ylab = "Cumulated time difference (s)")
  #hist(as.numeric(diff(df_Ved_2s_60s_seq$Timestamp)))
  print(table(diff(df_Ved_2s_60s_seq$Timestamp)))

  # Export mean VeDBA per burst ID- sequence file
  df_Ved_2s_60s_seq = df_Ved_2s_60s_seq[, c("Timestamp", "VeDBA", "logvedba")] # Remove Burst ID column
  colnames(df_Ved_2s_60s_seq) <- c("Timestamp", "vedba", "log_vedba")
  write_parquet(df_Ved_2s_60s_seq, savename)
}
