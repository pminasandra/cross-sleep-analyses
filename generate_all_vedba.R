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
rolling_width = 1 #width of the vedba static computation window - default = 1 sec, rolling_width*sampling_rate must be an integer
fracNA_thresh = .5
outdir <- '~/EAS_shared/cross_sleep/working/Data/'
codedir <- '~/EAS_ind/astrandburg/code/'

#Burst params
burst_intervals <- c(4, 10, 30, 60, 2*60, 5*60, 10*60, 15*60)
burst_length <- 2 #burst length in seconds

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

#Wrapper function
run_and_save_burst_data <- function(filename, savename, burst_interval, burst_length, acc_sample_rate, start_timestamp = NULL){
  
  df_stand = read_parquet(filename)
  
  #get start timestamp if not specified - ceiling to the nearest hour
  if(is.null(start_timestamp)){
    first_non_na_timestamp <- df_stand$Timestamp[which(!is.na(df_stand$Timestamp))[1]]
    start_timestamp <- ceiling_date(first_non_na_timestamp, unit = 'hours')
  }
  
  # Burst IDs
  df_burst <- burst_sample(accdat = as.data.frame(df_stand), burst_interval = burst_interval, burst_length = burst_length, acc_sample_rate = acc_sample_rate, start_timestamp = start_timestamp)
  # Remove inter-burst interval data
  df_burst = subset(df_burst, is.na(df_burst$Burst_ID) != TRUE)
  # Compute VeDBA per burst ID 
  df_burst <- as.data.table(df_burst) # Make data.table
  df_burst = acc_to_vedba(df_burst, rolling_mean_width = rolling_width, group_col = "Burst_ID")
  
  # Optional - remove VeDBA values if fracNA >= 50 % #Arbitrary!
  # Identify non-POSIXct columns (i.e., not Timestamp)
  non_time_cols <- names(df_burst)[!sapply(df_burst, inherits, what = "POSIXct")]
  non_time_cols = subset(non_time_cols, non_time_cols != "fracNA") # Keep value in this column to aid diagnostics
  # Set values to NA where fracNA >= 0.5
  df_burst[fracNA >= fracNA_thresh, (non_time_cols) := lapply(.SD, function(x) NA), .SDcols = non_time_cols]
  
  # Log-transform VeDBA
  df_burst[, log_VeDBA := log(VeDBA)]
  
  # Create a column that rounds down to the nearest second
  df_burst[, Timestamp_sec := as.POSIXct(floor(as.numeric(Timestamp)), origin = "1970-01-01", tz = attr(Timestamp, "tzone"))]
  
  # Median per ID
  df_burst_median <- df_burst[
    ,
    .(
      timestamp = first(Timestamp_sec),  # retain first timestamp
      vedba = median(VeDBA, na.rm = TRUE),
      logvedba = median(log_VeDBA, na.rm = TRUE)
    ),
    by = Burst_ID
  ]
  
  # Export mean VeDBA per burst ID- sequence file
  df_burst_median = df_burst_median[, c("timestamp", "vedba", "logvedba")] # Remove Burst ID column
  write_parquet(df_burst_median, savename)
}

#specify default time zone as UTC
Sys.setenv(TZ='UTC')

#set working directory
setwd(outdir)

#Ensure decimal secs are displayed
options(digits.secs = 3) 

#get directory where raw data is stored
rawdir <- paste0(outdir, 'Acc/')

#list all files that are parquets in the raw acc directory
files <- list.files(rawdir, recursive = T, pattern = '.parquet$', full.names = T)

#get basenames for the files
basenames <- gsub('[.]parquet','',basename(files))

#loop over files and generate + save all data
for(i in 2:length(files)){
  print(paste('Running all processing on file',i))
  print(paste('File path =', files[i]))
  df <- read_parquet(files[i])
  
  df <- as.data.table(df) # Make data.table
  df[, Timestamp := as.POSIXct(Timestamp, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")] # Ensure UTC for consistency 
  
  # 1)  Standardize acc to uniform_sampling - standardizing the time series to a uniform sampling grid with NA for missing records
  out = standardize_acc_to_uniform_sampling(df)
  df_stand <- out$df
  sampling_rate <- out$sampling_rate #also get the sampling rate from this function
  
  # 2) Clean data frame - make adjusted (rounded) time stamp the main one we use
  print("Getting standardized acc data")
  df_stand$Timestamp = df_stand$Timestamp_adj ; df_stand$Timestamp_adj = NULL
  
  # Write standarsised data as parquet file
  savename_standard <- paste0(outdir, 'Acc_standard/',basenames[i],'_standard.parquet')
  write_parquet(df_stand, savename_standard)
  
  # 3) Raw VeDBA all standardised data 
  print('Computing vedba from continuous data')
  if(!(round(rolling_width*sampling_rate) == rolling_width*sampling_rate)){
    stop('rolling_width*sampling_rate must be a whole number value')
  }
  df_Ved = acc_to_vedba(df_stand, rolling_mean_width = rolling_width*sampling_rate, group_col = NULL)
  
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
  savename_vedba_standard <- paste0(outdir,'VeDBA_cont/',basenames[i],'_VeDBA_standard.parquet')
  write_parquet(df_Ved_1Hz, savename_vedba_standard)
  
  #####################################################################################################################################
  
  # 7) Pseudo undersample (burst lengths and inter-intervals) and regenerate vedba
  print('Generating vedba for burst sampled data')
  for(j in 1:length(burst_intervals)){
    print(paste('Burst interval =', burst_intervals[j]))
    savename_burst <- paste0(outdir,'VeDBA_burst/',basenames[i],'_vedba_burst_int',burst_intervals[j],'_len',burst_length,'.parquet')
  run_and_save_burst_data(filename = savename_standard, savename = savename_burst, burst_intervals[j], burst_length, acc_sample_rate = sampling_rate)
  }
}