library(lubridate)

#This function gets contiguous periods of inactivity based on log vedba using a windowing approach.
#It works as follows:
# 1. For each time step, identify the individual as active or inactive using a threshold on log vedba (thresh)
# 2. Run a rolling window of size time_win seconds across the entire datastream, which is assumed to be continuous and uniformly sampled in time (but can have NAs)
#    For each time interval, compute the majority state (active or inactive) and assign it to that time point.
#    The rolling window is centered on the associated time point. In addition, NAs within this window are ignored.
# 3. Get contiguous sequences of inactivity. These get saved in the data frame inactive_periods.
#    The beginning of a contiguous sequence of activity is adjusted by adding time_win / 2 seconds and the end is
#    adjusted by subtracting time_win / 2 seconds. This is to make the estimates slightly more conservative and to account
#    for the time window's width.
# 4. Outputs a data frame called inactive_periods. If specified, also outputs plots
# INPUTS:
#   dat: data frame containing columns 'timestamp' (UTC time), 'vedba', and 'logvedba'
#   time_win: time window over which to do the smoothing (in seconds)
#   logvedba_thresh: threshold on logvedba below which the individual is considered inactive
#   make_plot: T or F, whether to make a plot
#   plot_savepath: full path to the file in which to save the plot (with .png extension)
# OUTPUTS:
#   inactive_periods: data frame containing the following columns:
#      $start_row: index of the starting row of the inactive period (in the original dat dataframe)
#      $end_row: index of the ending row of the inactive period (in the original dat dataframe)
#      $start_time_UTC: starting timestamp of the inactive period
#      $end_time_UTC: ending timestamp of the inactive period
#      $duration_hr: duration of the inactive period in hours
# If specified a plot is also produced and saved at the path plot_savepath
get_inactive_periods <- function(dat, time_win, logvedba_thresh, make_plot = F, plot_savepath = NULL) { 

  #add column 'active'
  dat$active <- dat$logvedba > logvedba_thresh
  
  #get timestep from the data
  dt <- as.numeric(difftime(dat$timestamp[2], dat$timestamp[1], units = 'secs'))
  
  #get window to use in units of samples (rather than times)
  win <- ceiling(time_win / dt)
  
  #rolling mean
  smooth_mean <- zoo::rollmean(dat$active, k = win, align = 'center', fill = NA, na.rm=T)
  
  #take majority state
  dat$active_smoothed <- smooth_mean > 0.5
    
  #get contiguous segments of inactivity
  runs <- rle(dat$active_smoothed)
  runs$starts <- cumsum(runs$lengths)
  starts <- c(1,runs$starts[1:length(runs$starts)-1]+1)
  ends <- c(starts[2:length(starts)]-1, runs$starts[length(runs$starts)])
  segments <- data.frame(start_row = starts, end_row = ends, state = runs$values)
  inactive_periods <- segments[which(segments$state==F),]
    
  #remove half a window length on either side
  inactive_periods$start_row <- inactive_periods$start_row + floor(win/2)
  inactive_periods$end_row <- inactive_periods$end_row - floor(win/2)
  
  #remove any negative or zero length durations (created due to the window adjustment at the end)
  inactive_periods <- inactive_periods[which(inactive_periods$end_row > inactive_periods$start_row),]
  
  inactive_periods$start_time_UTC <- dat$timestamp[inactive_periods$start_row]
  inactive_periods$end_time_UTC <- dat$timestamp[inactive_periods$end_row]
  inactive_periods$duration_hr <- as.numeric(difftime(inactive_periods$end_time_UTC, inactive_periods$start_time_UTC, units = 'hours'))
    
  #make plot
  if(make_plot){

    png(filename = plot_savepath, width = 1600, height = 1200, units='px')
    n_days <- as.numeric(difftime(max(dat$timestamp),min(dat$timestamp), units = 'days'))
    par(mfrow=c(ceiling(n_days/2),1),mar=c(0,0,0,0))
    plot_start_intervals <- seq(1, nrow(dat), by = floor(24*60*60*2/dt)) #each row is 2 days
    for(p in 1:length(plot_start_intervals)){
      idxs <- plot_start_intervals[p]:(plot_start_intervals[p] + floor(24*60*60*2/dt))
      plot(NULL, xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(min(idxs),max(idxs)),ylim=c(0,1))
      abline(v = idxs[which(!dat$active[idxs])], col = 'gray',lwd=.1)
      abline(v = idxs[which(dat$active[idxs])], col = '#00AA00',lwd=.1)
      abline(v = which(diff(as.Date(dat$timestamp))==1)+1, col = 'yellow',lwd=2,lty=2) #midnight UTC
      if(nrow(inactive_periods)>0){
        arrows(inactive_periods$start_row, rep(.5, nrow(inactive_periods)), inactive_periods$end_row, rep(.5, nrow(inactive_periods)), length = .4, col = 'red',lwd=3, code = 3, angle = 90)
      }
    }
    dev.off()

  }
  
  return(inactive_periods)
  
}


