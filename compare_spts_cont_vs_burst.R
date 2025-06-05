#Compare computed SPTs across burst lengths

#directory where inactivity periods / SPT data are stored
dir <- '~/EAS_shared/cross_sleep/working/Data/inactivity_periods/'

#read in data across all files and store in data frame
all_files <- list.files(path = dir)
spt_dat_all <- data.frame()
for(f in 1:length(all_files)){
  
  #load file
  file <- all_files[f]
  load(paste0(dir,'/',file))
  
  #get base id
  base_id <- gsub('spt_','',strsplit(file,'_vedba')[[1]][1])
  
  #if continuous, type is continuous
  if(grepl('vedba_standard',file)){
    inactive_periods$type <- 'cont'
    inactive_periods$burst_interval <- NA
    inactive_periods$burst_length <- NA
  } else{ #otherwise it's a burst file, and get also burst parameters
    inactive_periods$type <- 'burst'
    splitfile <- strsplit(gsub('[.]RData','',file),'_')[[1]]
    inactive_periods$burst_interval <- gsub('int','',splitfile[length(splitfile)-1])
    inactive_periods$burst_length <- gsub('len','',splitfile[length(splitfile)])
  }
  inactive_periods$base_id <- base_id
  
  spt_dat_all <- rbind(spt_dat_all, inactive_periods)
}

spt_dat_all <- spt_dat_all[which(spt_dat_all$spt==T),]

#for each file id and each date, get the difference in onset / offset of SPT from continuous
base_ids <- unique(spt_dat_all$base_id)
spt_dat_all$spt_start_offset <- spt_dat_all$spt_end_offset <- NA
for(i in 1:length(base_ids)){
  
  #get current data for that individual
  curr_idxs <- which(spt_dat_all$base_id == base_ids[i])
  curr <- spt_dat_all[curr_idxs,]
  
  #get all dates
  dates <- unique(curr$date)
  
  for(j in 1:length(dates)){
    cont_idx <- which(curr$base_id == base_ids[i] & curr$type == 'cont' & curr$date==dates[j] & curr$spt==T)
    if(length(cont_idx)>0){
      spt_start_cont <- curr$start_time_UTC[cont_idx]
      spt_end_cont <- curr$end_time_UTC[cont_idx]
      burst_idxs <- which(spt_dat_all$base_id == base_ids[i] & spt_dat_all$type == 'burst' & spt_dat_all$date == dates[j])
      spt_dat_all$spt_start_offset[burst_idxs] <- difftime(spt_dat_all$start_time_UTC[burst_idxs], spt_start_cont, unit = 'secs')
      spt_dat_all$spt_end_offset[burst_idxs] <- difftime(spt_dat_all$end_time_UTC[burst_idxs], spt_end_cont, unit = 'secs')
    }
  }
}

#plot as a function of burst interval
#SPT onset error
burst_intervals <- as.numeric(unique(spt_dat_all$burst_interval))
plot(NULL, xlim = c(.1,max(burst_intervals,na.rm=T)/60), ylim = c(-10000/60,10000/60),log='x', xlab = 'Burst interval (min)',ylab = 'SPT onset error (min)')
for(i in 1:length(burst_intervals)){
  idxs <- which(spt_dat_all$burst_interval == burst_intervals[i])
  jitter <- rnorm(length(idxs))*burst_intervals[i]/1000
  points(rep(burst_intervals[i],length(idxs))/60 + jitter, spt_dat_all$spt_start_offset[idxs]/60, pch = 19, cex = 0.5, col = '#0000FF44')
  
}
abline(h=0, lty = 2)
abline(h=60,lty = 2, col='gray')
abline(h=-60,lty = 2, col='gray')

#SPT offset error
plot(NULL, xlim = c(.1,max(burst_intervals,na.rm=T)/60), ylim = c(-10000/60,10000/60),log='x', xlab = 'Burst interval (min)',ylab = 'SPT offset error (min)')
for(i in 1:length(burst_intervals)){
  idxs <- which(spt_dat_all$burst_interval == burst_intervals[i])
  jitter <- rnorm(length(idxs))*burst_intervals[i]/1000
  points(rep(burst_intervals[i],length(idxs))/60 + jitter, spt_dat_all$spt_end_offset[idxs]/60, pch = 19, cex = 0.5, col = '#FF000044')
  
}
abline(h=0, lty = 2)
abline(h=60,lty = 2, col='gray')
abline(h=-60,lty = 2, col='gray')
