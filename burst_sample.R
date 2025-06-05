#Test code for burst sampling

#-----Burst sampling function-----
#Simulates burst sampling of acc data at a given interval burst_interval and for a given burst_length (in seconds)
#TODO: deal with when start time has NA data - not yet implemented

#INPUTS:
# accdat: data frame containing columns:
#   $Timestamp: timestamp in UTC (POSIXct)
#   $X, $Y, $Z: X, Y, and Z values of ACC data (numeric)
# burst_interval: how often a burst should be sampled, in seconds (must be an integer)
# burst_length: how long a burst should be, in seconds (must be an integer)
# acc_sample_rate: sample rate of the ACC data in Hz (must be an integer)
# start_timestamp: what timestamp to start at (POSIXct)
#
#OUTPUTS:
# accdat: same table as before, but with a column Burst_ID (numeric) specifying which burst each row is part of (if not part of a burst, Burst_ID is NA)
#REQUIRES:
# lubridate
burst_sample <- function(accdat, burst_interval, burst_length, acc_sample_rate, start_timestamp = NULL){
  
  #---Input checking---
  #check that data are the right format - must contain colnames 'Timestamp' 'X' 'Y' 'Z'
  expected_cols <- c('Timestamp','X','Y','Z')
  if(!setequal(intersect(colnames(accdat),expected_cols),expected_cols)){
    stop('accdat must include columns named Timestamp, X, Y, and Z')
  }
  
  #Check that timestamps are POSIXct
  if(!('POSIXct' %in% class(accdat$Timestamp))){
    stop('Timestamp column must be in POSIXct format')
  }
  #Check that X, Y, and Z are numeric (at least first element)
  for(col in c('X','Y','Z')){
    if(!('numeric' %in% class(accdat[1,col][[1]]))){
      stop('Timestamp column must be in POSIXct format')
    }
  }
  
  #Check that start_timestamp is either NULL or POSIXct
  if(!is.null(start_timestamp) & !is.POSIXct(start_timestamp)){
    stop('start_timestamp must be in POSIXct format')
  }
  
  #Check that other parameters are integers
  if(!(round(burst_interval)==burst_interval)){
    stop('burst_interval must be an integer')
  }
  if(!(round(burst_length)==burst_length)){
    stop('burst_length must be an integer')
  }
  if(!(round(acc_sample_rate)==acc_sample_rate)){
    stop('acc_sample_rate must be an integer')
  }
  
  #----Main----
  
  #Get start_timestamp and end_timestamp, if not specified
  #Take first full-minute timestamp for start timestamp
  if(is.null(start_timestamp)){
    first_timestamp <- accdat$Timestamp[1]
    start_timestamp <- lubridate::ceiling_date(first_timestamp, unit = 'min')
  }
  #Take last timestamp in the table as the end_timestamp (not rounded)
  end_timestamp <- accdat$Timestamp[nrow(accdat)]
  
  #Get index of first sample and last index
  idx0 <- which(accdat$Timestamp == start_timestamp)
  idxf <- nrow(accdat)
  
  #Get indexes to burst starts and ends
  idxs_starts <- seq(from = idx0, to = idxf, by = burst_interval*acc_sample_rate)
  idxs_ends <- idxs_starts + burst_length*acc_sample_rate
  
  #If last burst goes over the end of the table, remove it
  if(idxs_ends[length(idxs_ends)] > nrow(accdat)){
    idxs_starts <- idxs_starts[1:(length(idxs_starts)-1)]
    idxs_ends <- idxs_ends[1:(length(idxs_ends)-1)]
  }
  
  #Label bursts in the table
  accdat$Burst_ID <- NA
  n_bursts <- length(idxs_starts)
  burst_IDs <- rep(1:n_bursts, each = burst_length*acc_sample_rate)
  row_IDs <- rep(idxs_starts, each = burst_length*acc_sample_rate) + seq_len(burst_length*acc_sample_rate) -1
  accdat$Burst_ID[row_IDs] <- burst_IDs
  
  #TODO: deal with NA start times - not yet dealt with (but might not be needed? could just consider them missing data)
  na_starts <- which(is.na(accdat$X[idxs_starts]))
  
  return(accdat)
}

#----Test script----

run_test <- F

if(run_test){

  #Set computer time zone to UTC - simplifies everything
  Sys.setenv(TZ='UTC')
  
  library(arrow) #needed for read_parquet
  library(lubridate) #for round_date
  
  #read in an example file
  accdat <- arrow::read_parquet('~/EAS_shared/cross_sleep/working/Data/Acc/hyrax/hyrax_bachelor_CF.parquet')
  
  #Params
  burst_interval <- 300 #burst interval in seconds
  burst_length <- 6 #burst length in seconds
  acc_sample_rate <- 50 #samples per second in the raw ACC data
  start_timestamp <- as.POSIXct('2023-07-02 8:20:00', tz = 'UTC')
  
  #test the function
  test_out <- burst_sample(accdat = accdat, burst_interval = burst_interval, burst_length = burst_length, acc_sample_rate = acc_sample_rate, start_timestamp = start_timestamp)
  
}
    
    
