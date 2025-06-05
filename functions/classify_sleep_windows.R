classify_sleep_windows <- function(df_vedba, tag, thresh, mov_window = 9, block_size = 3, waso_block = 3, gap_size = 60) {
  df_vedba <- df_vedba %>%
    dplyr::mutate(
      roll_log_vedba = zoo::rollmedian(log_vedba, mov_window, fill = NA, align = "center")
    )
  
  df_vedba <- df_vedba %>%
    dplyr::mutate(
      sleep_per = 0,
      pot_sleep = as.numeric(log_vedba < thresh),
      sleep_bouts = 0,
      n_bursts = n(),
      max_time_diff = max(diff(as.numeric(Timestamp)), na.rm = TRUE)
    )
  
  for (night in unique(df_vedba$night_date)) {
    night_dat <- df_vedba[df_vedba$night_date == night, ]
    if (nrow(night_dat) < 10) next
    
    temp <- rle(as.numeric(night_dat$roll_log_vedba < thresh))
    sleep_per_runs <- rep(temp$lengths > block_size, times = temp$lengths)
    sleep_per_sleep_bouts <- as.numeric(night_dat$roll_log_vedba < thresh & sleep_per_runs)
    
    diffs <- diff(c(0, sleep_per_sleep_bouts))
    starts <- which(diffs == 1)[-1]
    ends <- which(diffs == -1)
    
    if (length(starts) > 0 && length(ends) > 0) {
      ends <- ends[1:length(starts)]
      gaps <- as.numeric(night_dat$Timestamp[starts] - night_dat$Timestamp[ends], units = "mins")
      inds_to_remove <- which(gaps < gap_size)
      
      onset <- if (length(inds_to_remove) == 0) which(diffs == 1)
      else which(diffs == 1)[-(inds_to_remove + 1)]
      wake <- if (length(inds_to_remove) == 0) ends
      else ends[-inds_to_remove]
      
      if (length(onset) > 0 && length(wake) > 0) {
        longest <- which.max(as.numeric(night_dat$Timestamp[wake] - night_dat$Timestamp[onset], units = "secs"))
        start_time <- night_dat$Timestamp[onset[longest]]
        end_time <- night_dat$Timestamp[wake[longest]]
        
        df_vedba[df_vedba$night_date == night & df_vedba$Timestamp >= start_time & df_vedba$Timestamp <= end_time, "sleep_per"] <- 1
      }
    }
    
    temp2 <- rle(as.numeric(night_dat$log_vedba < thresh))
    runs <- rep(temp2$lengths >= waso_block, times = temp2$lengths)
    df_vedba[df_vedba$night_date == night, "sleep_bouts"] <- as.numeric(night_dat$log_vedba < thresh & runs)
  }
  
  return(df_vedba)
  
  
}
