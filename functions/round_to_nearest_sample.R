#Helper function to round timestamps to the nearest sample
round_to_nearest_sample <- function(timestamps, sampling_rate){
  
  interval <- 1 / sampling_rate
  
  rounded <- as.POSIXct(
    round(as.numeric(timestamps) / interval) * interval,
    origin = "1970-01-01", tz = "UTC"
  )
  
  return(rounded)
}
