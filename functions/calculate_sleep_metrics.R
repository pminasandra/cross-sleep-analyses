calculate_sleep_metrics <- function(df_classified, tag, frag_block = 3) {
  # Create study night range
  nights <- unique(df_classified$night_date)
  
  # Initialize output dataframe
  sleep_per <- data.frame(
    tag = rep(tag, each = length(nights)),
    night = rep(seq_along(nights), times = length(tag)),  # running number
    night_date = rep(nights, times = length(tag)),
    total_pot_sleep = NA,
    total_sleep_bouts = NA,
    onset = NA,
    waking = NA,
    SPT = NA,
    WASO = NA,
    TST = NA,
    sleep_eff = NA,
    wake_bouts = NA,
    frag_wake_bouts = NA,
    summed_VeDBA = NA,
    night_VeDBA_corr = NA,
    ave_vedba = NA,
    dark_pot_sleep = NA,
    dark_ave_vedba = NA,
    max_time_diff = NA,
    n_bursts = NA
  )
  
  sleep_durs <- c()
  wake_durs <- c()
  
  for (night in nights) {
    message(paste(tag, night))
    night_dat <- df_classified[df_classified$night_date == night, ]
    night_dat <- night_dat[order(night_dat$Timestamp), ]
    
    idx <- sleep_per$night_date == night
    sleep_per[idx, "n_bursts"] <- unique(night_dat$n_bursts)
    sleep_per[idx, "max_time_diff"] <- unique(night_dat$max_time_diff)
    sleep_per[idx, "total_pot_sleep"] <- sum(night_dat$pot_sleep, na.rm = TRUE)
    sleep_per[idx, "total_sleep_bouts"] <- sum(night_dat$sleep_bouts, na.rm = TRUE)
    sleep_per[idx, "ave_vedba"] <- mean(night_dat$log_vedba, na.rm = TRUE)
    
    SPT_dat <- night_dat[night_dat$sleep_per == 1, ]
    if (nrow(SPT_dat) > 0) {
      onset <- min(SPT_dat$Timestamp)
      waking <- max(SPT_dat$Timestamp)
      SPT <- as.numeric(difftime(waking, onset, units = "mins")) + 1
      WASO <- sum(SPT_dat$sleep_bouts == 0, na.rm = TRUE)
      TST <- sum(SPT_dat$sleep_bouts == 1, na.rm = TRUE)
      sleep_eff <- TST / nrow(SPT_dat)
      summed_VeDBA <- sum(SPT_dat$log_VeDBA, na.rm = TRUE)
      night_VeDBA_corr <- summed_VeDBA / SPT
      
      sleep_per[idx, "onset"] <- onset
      sleep_per[idx, "waking"] <- waking
      sleep_per[idx, "SPT"] <- SPT
      sleep_per[idx, "WASO"] <- WASO
      sleep_per[idx, "TST"] <- TST
      sleep_per[idx, "sleep_eff"] <- sleep_eff
      sleep_per[idx, "summed_VeDBA"] <- summed_VeDBA
      sleep_per[idx, "night_VeDBA_corr"] <- night_VeDBA_corr
      
      # Fragmented wake bouts
      temp_rle <- rle(SPT_dat$sleep_bouts)
      frag_runs <- as.numeric(rep(temp_rle$lengths >= frag_block, times = temp_rle$lengths))
      frag_wake_bouts <- as.numeric(SPT_dat$sleep_bouts == 0 & frag_runs == 1)
      frag_diffs <- diff(c(1, frag_wake_bouts))
      sleep_per[idx, "frag_wake_bouts"] <- sum(frag_diffs == 1, na.rm = TRUE)
      
      # Wake bouts
      wake_diffs <- diff(c(0, SPT_dat$sleep_bouts))
      sleep_per[idx, "wake_bouts"] <- sum(wake_diffs == -1, na.rm = TRUE)
      
      # Save durations
      rle_bouts <- rle(SPT_dat$sleep_bouts)
      sleep_durs <- c(sleep_durs, rle_bouts$lengths[rle_bouts$values == 1])
      wake_durs <- c(wake_durs, rle_bouts$lengths[rle_bouts$values == 0])
    }
  }
  
  # Reformat timestamps and add time-only columns
  sleep_per$onset <- as.POSIXct(sleep_per$onset, origin = "1970-01-01", tz = tz)
  sleep_per$waking <- as.POSIXct(sleep_per$waking, origin = "1970-01-01", tz = tz)
  sleep_per$onset_time <- as_hms(sleep_per$onset)
  sleep_per$waking_time <- as_hms(sleep_per$waking)
  
  return(sleep_per)
  
}
