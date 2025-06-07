#Compare inactive intervals across burst lengths

#directory where inactivity periods / SPT data are stored
dir <- '~/EAS_shared/cross_sleep/working/Data/inactive_periods'

#directory to store plots
plotdir <-'~/EAS_shared/cross_sleep/working/Figures/inactive_period_comparisons/'

#directory to save comparison tables
savedir <- '~/EAS_shared/cross_sleep/working/Data/inactive_period_comparisons/'

#what species to use (or NULL for all)
species_to_use <- 'hyena'
make_plots <- T
utc_offset <- -5 #this is just a hacky solution to avoid SPTs that start around midnight UTC which leads to errouneously bad results - only relevant for spider monkeys

#read in data across all files and store in data frame
all_files <- list.files(path = dir)

if(!is.null(species_to_use)){
  all_files <- all_files[grepl(species_to_use,all_files)]
}
inactive_periods_all <- data.frame()
for(f in 1:length(all_files)){
  
  #load file
  file <- all_files[f]
  load(paste0(dir,'/',file))
  
  for(i in 1:nrow(inactive_periods)){
    #if continuous, type is continuous
    if(grepl('vedba_standard',inactive_periods$filename[i])){
      inactive_periods$type[i] <- 'cont'
      inactive_periods$burst_interval[i] <- NA
      inactive_periods$burst_length[i] <- NA
    } else{ #otherwise it's a burst file, and get also burst parameters
      inactive_periods$type[i] <- 'burst'
      splitfile <- strsplit(gsub('[.]parquet','',inactive_periods$filename[i]),'_')[[1]]
      inactive_periods$burst_interval[i] <- gsub('int','',splitfile[length(splitfile)-1])
      inactive_periods$burst_length[i] <- gsub('len','',splitfile[length(splitfile)])
    }
  }
  inactive_periods$base_id <- paste(inactive_periods$species, inactive_periods$deploy_id, inactive_periods$ind_id,sep='_')
    
  inactive_periods_all <- rbind(inactive_periods_all, inactive_periods)
}

#hacky solution to SPTs starting around midnight being erroneously linked with the wrong day
if(species_to_use=='spidermonkey'){
  inactive_periods_all$start_time_UTC <- inactive_periods_all$start_time_UTC + utc_offset*60*60
  inactive_periods_all$end_time_UTC <- inactive_periods_all$end_time_UTC + utc_offset*60*60
  #also remove on bad file
  inactive_periods_all <- inactive_periods_all[which(inactive_periods_all$filename!='spidermonkey_1_Albus_vedba_burst_int180_len2.parquet'),]
}

#for each burst file compute:
# total time in inactive periods (burst - continuous)
# IoU of inactive periods between burst data and continuous data
# number of inactive periods (burst - continuous)
cont_dat <- inactive_periods_all[which(inactive_periods_all$type == 'cont'),]
burst_dat <- inactive_periods_all[which(inactive_periods_all$type == 'burst'),]

filenames <- unique(burst_dat$filename) #comp_dat will have one row per burst filename
comp_dat <- data.frame(filename = filenames)
comp_dat$tot_time_cont <- NA
comp_dat$tot_time_burst <- NA
comp_dat$iou <- NA
comp_dat$n_periods_cont <- NA
comp_dat$n_periods_burst <- NA
for(i in 1:nrow(comp_dat)){
  
  print(paste(i,'/',nrow(comp_dat)))
  filename <- comp_dat$filename[i]
  print(filename)
  
  #get continuous and burst sampled inactive periods for current file
  burst_curr <- burst_dat[which(burst_dat$filename == filename),]
  cont_curr <- cont_dat[which(cont_dat$base_id == burst_curr$base_id[1]),]
  
  comp_dat$base_id[i] <- cont_curr$base_id[1]
  comp_dat$species[i] <- cont_curr$species[1]
  comp_dat$type[i] <- 'burst'
  comp_dat$burst_interval[i] <- as.numeric(burst_curr$burst_interval[1])
  comp_dat$burst_length[i] <- as.numeric(burst_curr$burst_length[1])
  
  #---total time in an inactive period--
  comp_dat$tot_time_cont[i] <- sum(cont_curr$duration_hr, na.rm=T)
  comp_dat$tot_time_burst[i] <- sum(burst_curr$duration_hr, na.rm=T)
  
  #---IoU---
  min_time <- floor_date(min(min(burst_curr$start_time_UTC, na.rm=T),min(cont_curr$start_time_UTC, na.rm=T)), unit = 'secs')
  max_time <- ceiling_date(max(max(burst_curr$end_time_UTC, na.rm=T),max(cont_curr$end_time_UTC, na.rm=T)), unit = 'secs')
  
  #get time sequences of 1s and 0s (1 if in inactive period)
  timestamps_curr <- seq.POSIXt(min_time, max_time, by = 1)
  burst_seq <- cont_seq <- rep(0, length(timestamps_curr))
  for(j in 1:nrow(burst_curr)){
    idx0 <- match(burst_curr$start_time_UTC[j], timestamps_curr)
    idxf <- match(burst_curr$end_time_UTC[j], timestamps_curr)
    burst_seq[idx0:idxf] <- 1
  }
  for(j in 1:nrow(cont_curr)){
    idx0 <- match(cont_curr$start_time_UTC[j], timestamps_curr)
    idxf <- match(cont_curr$end_time_UTC[j], timestamps_curr)
    cont_seq[idx0:idxf] <- 1
  }
  comp_dat$iou[i] <- sum(cont_seq & burst_seq) / sum(cont_seq | burst_seq)
  
  #---number of time periods---
  comp_dat$n_periods_burst[i] <- nrow(burst_curr)
  comp_dat$n_periods_cont[i] <- nrow(cont_curr)
  
}

#---for each file and each date, define the "SPT" as the longest interval starting on that date
#compare SPT onset and offset between continuous and burst data

#create a data frame to hold each date within each file
comp_dat_date <- data.frame()
for(i in 1:nrow(comp_dat)){
  base_id <- comp_dat$base_id[i]
  species <- comp_dat$species[i]
  burst_interval <- comp_dat$burst_interval[i]
  burst_length <- comp_dat$burst_length[i]
  cont_curr <- cont_dat[which(cont_dat$base_id == base_id),]
  burst_curr <- burst_dat[which(burst_dat$base_id == base_id),]
  dates <- unique(date(c(cont_curr$start_time_UTC, burst_curr$start_time_UTC)))
  filename <- comp_dat$filename[i]
  tmp <- data.frame(filename = rep(filename, length(dates)),
                                  date = dates)
  tmp$base_id <- base_id
  tmp$species <- species
  tmp$type <- 'burst'
  tmp$burst_interval <- burst_interval
  tmp$burst_length <- burst_length
  comp_dat_date <- rbind(comp_dat_date, tmp)
}

#add in continuous data rows
base_ids <- unique(comp_dat_date$base_id)
for(i in 1:length(base_ids)){
  dates <- unique(comp_dat_date$date[which(comp_dat_date$base_id == base_ids[i])])
  cont_filename <- inactive_periods_all$filename[which(inactive_periods_all$type == 'cont' & inactive_periods_all$base_id == base_ids[i])][1]
  tmp <- data.frame(filename = rep(cont_filename, length(dates)),
                    date = dates)
  tmp$base_id <- base_ids[i]
  tmp$species <- strsplit(base_ids[i],'_')[[1]][1]
  tmp$type <- 'cont'
  tmp$burst_interval <- NA
  tmp$burst_length <- NA
  comp_dat_date <- rbind(comp_dat_date, tmp)
}

#get the spt for each file/date combo
comp_dat_date$spt_onset_UTC <- comp_dat_date$spt_offset_UTC <- as.POSIXct(NA, tz = 'UTC')
for(i in 1:nrow(comp_dat_date)){
  curr <- inactive_periods_all[which(inactive_periods_all$filename == comp_dat_date$filename[i] & date(inactive_periods_all$start_time_UTC) == comp_dat_date$date[i]),]
  
  if(nrow(curr)>0){
    longest <- which(curr$duration_hr == max(curr$duration_hr,na.rm=T))
    comp_dat_date$spt_onset_UTC[i] <- curr$start_time_UTC[longest]
    comp_dat_date$spt_offset_UTC[i] <- curr$end_time_UTC[longest]
  }
  
}

#compare to continuous data SPT
comp_dat_date$spt_onset_error <- comp_dat_date$spt_offset_error <- NA
for(i in 1:nrow(comp_dat_date)){
  cont_idx <- which(comp_dat_date$base_id == comp_dat_date$base_id[i] & 
                      comp_dat_date$type == 'cont' &
                      comp_dat_date$date == comp_dat_date$date[i])
  comp_dat_date$spt_onset_error[i] <- difftime(comp_dat_date$spt_onset_UTC[cont_idx],comp_dat_date$spt_onset_UTC[i], units = 'secs')
  comp_dat_date$spt_offset_error[i] <- difftime(comp_dat_date$spt_offset_UTC[cont_idx],comp_dat_date$spt_offset_UTC[i], units = 'secs')
}


#save data
if(is.null(species_to_use)){
  paste0(savedir,'inactive_period_comparisons_all.RData')
} else{
  savename <- paste0(savedir,'inactive_period_comparisons_',species_to_use,'.RData')
}
save(list = c('comp_dat','comp_dat_date'), file = savename)


#-----PLOTTING---
burst_dat_date <- comp_dat_date[which(comp_dat_date$type=='burst'),]

if(make_plots){
  #SPT onset error
  if(is.null(species_to_use)){
    plotname <- paste0(plotdir, 'spt_onset_error_all.png')
  } else{
    plotname <- paste0(plotdir, 'spt_onset_error_',species_to_use,'.png')
  }
  png(filename = plotname, width = 10, height = 8, units = 'in', res = 300)
  burst_intervals <- as.numeric(unique(burst_dat_date$burst_interval))
  plot(NULL, xlim = c(.1,max(burst_intervals,na.rm=T)/60), ylim = c(-120,120),log='x', xlab = 'Burst interval (min)',ylab = 'SPT onset error (min)', main = species_to_use)
  abline(h=0, lty = 2)
  abline(h=seq(30,120,30),lty = 2, col='gray')
  abline(h=-seq(30,120,30),lty = 2, col='gray')
  for(i in 1:length(burst_intervals)){
    idxs <- which(burst_dat_date$burst_interval == burst_intervals[i])
    jitter <- rnorm(length(idxs))*burst_intervals[i]/1000
    points(rep(burst_intervals[i],length(idxs))/60 + jitter, burst_dat_date$spt_onset_error[idxs]/60, pch = 19, cex = 0.5, col = '#0000FF22')
    upper <- quantile(burst_dat_date$spt_onset_error[idxs]/60, 0.975,na.rm=T)
    lower <- quantile(burst_dat_date$spt_onset_error[idxs]/60, 0.025,na.rm=T)
    arrows(burst_intervals[i]/60,upper,burst_intervals[i]/60,lower, col = 'blue', angle = 90, len = 0.05, code = 3)
  }
  dev.off()
  
  #SPT offset error
  if(is.null(species_to_use)){
    plotname <- paste0(plotdir, 'spt_offset_error_all.png')
  } else{
    plotname <- paste0(plotdir, 'spt_offset_error_',species_to_use,'.png')
  }
  png(filename = plotname, width = 10, height = 8, units = 'in', res = 300)
  burst_intervals <- as.numeric(unique(burst_dat_date$burst_interval))
  plot(NULL, xlim = c(.1,max(burst_intervals,na.rm=T)/60), ylim = c(-120,120),log='x', xlab = 'Burst interval (min)',ylab = 'SPT offset error (min)', main = species_to_use)
  abline(h=0, lty = 2)
  abline(h=seq(30,120,30),lty = 2, col='gray')
  abline(h=-seq(30,120,30),lty = 2, col='gray')
  for(i in 1:length(burst_intervals)){
    idxs <- which(burst_dat_date$burst_interval == burst_intervals[i])
    jitter <- rnorm(length(idxs))*burst_intervals[i]/1000
    points(rep(burst_intervals[i],length(idxs))/60 + jitter, burst_dat_date$spt_offset_error[idxs]/60, pch = 19, cex = 0.5, col = '#FF000022')
    upper <- quantile(burst_dat_date$spt_offset_error[idxs]/60, 0.975,na.rm=T)
    lower <- quantile(burst_dat_date$spt_offset_error[idxs]/60, 0.025,na.rm=T)
    arrows(burst_intervals[i]/60,upper,burst_intervals[i]/60,lower, col = 'red', angle = 90, len = 0.05, code = 3)
  }
  dev.off()
  
  #Inactive intervals IoU
  if(is.null(species_to_use)){
    plotname <- paste0(plotdir, 'iou_all.png')
  } else{
    plotname <- paste0(plotdir, 'iou_',species_to_use,'.png')
  }
  png(filename = plotname, width = 10, height = 8, units = 'in', res = 300)
  burst_intervals <- as.numeric(unique(comp_dat$burst_interval))
  plot(NULL, xlim = c(.1,max(burst_intervals,na.rm=T)/60), ylim = c(0.5,1),log='x', xlab = 'Burst interval (min)',ylab = 'Inactive intervals IoU (burst vs cont)', main = species_to_use)
  abline(h=1, lty = 2)
  abline(h=seq(0,.9,.1), lty = 2, col = 'gray')
  for(i in 1:length(burst_intervals)){
    idxs <- which(comp_dat$burst_interval == burst_intervals[i])
    jitter <- rnorm(length(idxs))*burst_intervals[i]/1000
    points(rep(burst_intervals[i],length(idxs))/60 + jitter, comp_dat$iou[idxs], pch = 19, cex = 2, col = '#00000044')
  }
  dev.off()
}
