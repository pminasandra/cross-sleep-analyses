#This script gets contiguous periods of inactivity based on log vedba and idenfies the "spt", i.e. the longest inactive period
#It works as follows:
# 1. For each time step, identify the individual as active or inactive using a threshold on log vedba (manual_set_thresh or pulled from metadata table)
# 2. Run a rolling window of size time_win seconds across the data from start_hour_UTC to start_hour_UTC on the next day
# For each time interval, compute the majority state (active or inactive) and assign it to that time point.
# The rolling window is centered on the associated time point. In addition, NAs within this window are ignored.
# 3. Get contiguous sequences of inactivity. These get saved in the data frame inactive_periods
# 4. Find the longest duration sequence of inactivity - this is defined as the SPT (sleep period time)
#The code currently does not save the inactivity_periods data frame, but it does save plots in plotdir
#In the plots, the SPT is colored red, and other inactivity periods are colored blue.

#----PARAMS----
time_win <- 60*60 #time window in seconds 
manual_set_thresh <- -2.5 #logvedba threshold to use (if NULL, this will be pulled from metadata table)
start_hour_UTC <- 11 #hour (in UTC) to use as the start of the window (should be somewhere in the middle of the 'active period')
make_plot <- T #T or F, whether to make plots
species_to_use <- 'meerkat' #can be hyena, meerkat, spidermonkey, or hyrax

#directories and files
datadir <- '~/EAS_shared/cross_sleep/working/Data/'
plotdir <- '~/EAS_shared/cross_sleep/working/Figures/spt_majority_windowing/'
metadata_file <- 'cross_sleep_metadata_2025-05-20.csv'

#----MAIN----
#libraries
library(arrow)
library(lubridate)

#get all files for the target species
files <- list.files(paste0(datadir,'VeDBA/'))
files <- files[grep(species_to_use,files)]

#loop over files
for(f in 1:length(files)){
  
  filename <- files[f]
  print(filename)

  #get info from the filename
  splitname <- strsplit(basename(filename),split='_')[[1]]
  species <- splitname[1]
  deploy_id <- splitname[2]
  ind_id <- splitname[3]
  
  #main
  setwd(datadir)
  
  #Load in file
  dat <- arrow::read_parquet(paste0(datadir,'VeDBA/',filename))
  
  #replace infinity values with NAs (fix in original files later)
  dat$logvedba[which(is.infinite(dat$logvedba))] <- NA
  
  #get threshold from metadata file if needed (or set manually)
  if(is.null(manual_set_thresh)){
    metadata <- read.csv(metadata_file)
    row <- which(metadata$individual_ID==ind_id & metadata$deployment_ID==deploy_id & metadata$species==species)
    metadata_ind <- metadata[row,]
    thresh <- metadata_ind$th1
  } else{
    thresh <- manual_set_thresh
  }
  
  #add column 'active'
  dat$active <- dat$logvedba > thresh
  
  #get timestep from the data
  dt <- as.numeric(difftime(dat$timestamp[2], dat$timestamp[1], units = 'secs'))
  
  #get window to use in units of samples (rather than times)
  win <- ceiling(time_win / dt)
  
  #rolling mean
  smooth_mean <- zoo::rollmean(dat$active, k = win, align = 'center', fill = NA, na.rm=T)
  
  #take majority state
  dat$active_smoothed <- smooth_mean > 0.5
  
  #get indexes of day starts and ends
  first_start <- min(which(hour(dat$timestamp)==start_hour_UTC))
  day_starts <- seq(first_start, nrow(dat), by = floor(24*60*60/dt))
  
  #save parameters used in an object
  params <- list()
  params$time_win <- time_win
  params$vedba_thresh <- thresh
  params$start_hour_UTC <- start_hour_UTC
  params$species <- species
  params$deploy_id <- deploy_id
  params$ind_id <- ind_id
  params$filename <- filename
  
  #data for each day - plot and process
  inactive_periods_all <- data.frame()
  for(day in 1:length(day_starts)){
    row_start <- day_starts[day]
    day_duration <- 24*60*60 / dt 
    row_end <- row_start + day_duration - 1
    idxs <- row_start:row_end
    
    #get contiguous segments of inactivity
    runs <- rle(dat$active_smoothed[idxs])
    runs$starts <- cumsum(runs$lengths)
    starts <- c(1,runs$starts[1:length(runs$starts)-1])
    ends <- runs$starts-1
    segments <- data.frame(start_row = starts + row_start - 1, end_row = ends + row_start - 1, state = runs$values)
    inactive_periods <- segments[which(segments$state==F),]
    
    #remove half a window length on either side
    inactive_periods$start_row <- inactive_periods$start_row + floor(win/2)
    inactive_periods$end_row <- inactive_periods$end_row - floor(win/2)
    
    inactive_periods$start_time_UTC <- dat$timestamp[inactive_periods$start_row]
    inactive_periods$end_time_UTC <- dat$timestamp[inactive_periods$end_row]
    inactive_periods$duration_hr <- as.numeric(difftime(inactive_periods$end_time_UTC, inactive_periods$start_time_UTC, units = 'hours'))
    inactive_periods$spt <- F
    inactive_periods$spt[which(inactive_periods$duration_hr == max(inactive_periods$duration_hr, na.rm=T))] <- T
    inactive_periods$date <- as.Date(dat$timestamp[idxs[1]])
    
    #save in inactive_periods_all
    inactive_periods_all <- rbind(inactive_periods_all, inactive_periods)
    
    #make plot
    if(make_plot){
      
      #make directory if needed
      plotdir_specific <- paste0(plotdir,species)
      if(!dir.exists(plotdir_specific)){
        dir.create(plotdir_specific)
      }
      
      plotname <- paste0(plotdir_specific,'/spt_',time_win/60,'min_',gsub('[.]parquet','',filename),'_',as.Date(dat$timestamp[idxs[1]]),'.png')
      png(filename = plotname, width = 800, height = 400, units='px')
      ymin <- quantile(dat$logvedba,.05,na.rm=T)
      ymax <- quantile(dat$logvedba,.95,na.rm=T)
      plot(idxs, dat$logvedba[idxs],type='l',col='#00000033',ylim=c(ymin,ymax),xaxt='n',ylab = 'Log VeDBA',xlab='Time UTC (hr)')
      abline(h=thresh,lty=2)
      axis(1, at = seq(idxs[1],idxs[length(idxs)],length.out=2),labels=rep(start_hour_UTC,2))
      for(i in 1:nrow(inactive_periods)){
        if(inactive_periods$spt[i]){
          col <- '#FF000044'
        } else{
          col <- '#0000FF22'
        }
        polygon(c(inactive_periods$start_row[i], inactive_periods$start_row[i],inactive_periods$end_row[i],inactive_periods$end_row[i]),
                c(ymin-1,ymax+1,ymax+1,ymin-1), col = col, border=NA)
      }
      dev.off()
    }
  }
  
  inactive_periods <- inactive_periods_all
  
  savename <- paste0(datadir,'inactivity_periods/spt_',gsub('[.]parquet','',filename),'.RData')
  save(list = c('inactive_periods','params'), file = savename)
}
