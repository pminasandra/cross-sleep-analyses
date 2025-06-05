
{
  library(httr)
  library( hms )
  library( data.table )
  library( stringr )
  library( lubridate )
  library( suncalc )
  library( dplyr )
  library( purrr )
  library( HDInterval )
  library(tidyr) 
  library( sp )
  library( stats )
  library( reshape2 )
  library( plyr )
  library(mixtools)
  library(arrow)
  library(zoo)
}

################# Functions #########################

## function for normalizing a vector
normalize_func <- function( x ) return( (x - mean( x, na.rm = T ) )/ sd( x, na.rm = T ) )

## function for plotting times from noon to noon. It will make correct 12:00 - 24:00 to be 00:00 - 12:00 and 00:00 - 12:00 to be 12:00 to 24:00
calculate_threshold <- function(log_vedba_vector, thres_value = 0.1) {
  # Remove NA values
  log_vedba_vector <- log_vedba_vector[!is.na(log_vedba_vector)]
  
  # Standardize the data
  #log_vedba_vector <- scale(log_vedba_vector)
  
  # Check if the data is valid
  if (length(log_vedba_vector) < 2) {
    stop("Insufficient data to calculate threshold.")
  }
  
  # Fit a Gaussian mixture model with 2 components
  fit <- tryCatch({
    normalmixEM(log_vedba_vector, k = 2, arbmean = TRUE, arbvar = TRUE)
  }, error = function(e) {
    stop("Error fitting model: ", e$message)
  })
  
  # Identify the leftmost distribution
  sorted_indices <- order(fit$mu)
  leftmost_index <- sorted_indices[1]
  
  # Calculate the posterior probabilities for the leftmost component
  posterior_probs_leftmost <- fit$posterior[, leftmost_index]
  
  # Calculate the leftmost point where the probability of being from the left distribution is under 25%
  threshold <- min(log_vedba_vector[
    log_vedba_vector > max(log_vedba_vector[posterior_probs_leftmost > thres_value]) &
      posterior_probs_leftmost < thres_value])
  
  return(threshold)
}

################# Determining sleep periods with modification of Van Hees et al. 2018 method ###################

################## Read in the d1 (accelerometer burst) data ###################

## d1 is a dataframe with a row for each minute for each baboon. Each row contains the raw (or interpolated) GPS acc burst, and several different measures calculated from bursts (like VeDBA)
d1 <- read_parquet("/mnt/EAS_shared/cross_sleep/working/Data/VeDBA_cont/spidermonkey_1_Albus_vedba_standard.parquet")
## turn the data table into a dataframe

d1$Timestamp <- with_tz(d1$Timestamp, "America/Panama")

d1 <- d1[d1$Timestamp == floor_date(d1$Timestamp, unit = "second")]
d1[, log_VeDBA := log(VeDBA)]
d1 <- as.data.frame( d1 )
d1$tag <- "spidermonkey_1_Albus"

# Define the coordinates for Mpala Research Centre
site_lat <- 47.6667
site_lon <- 9.1667

## make a column for local time
#d1 <- d1[ !duplicated( d1[ , c( 'tag', 'Timestamp' ) ] ), ]

## assign each minute of data to a given night. A night lasts from noon to noon. First, apply a time shift so that each night is a unit, and not each day
time_shift <- d1$Timestamp - 17*60*60

## save the date of the first night of the study (the date of the night is always the date of the evening at the beginning of that night; so the first night of the study is 2012-07-31, although the data starts on 2012-08-01, because the data on that first morning is still technically part of the data for the previous night, as a night is noon to noon)
start_date <- as.Date(min(d1$Timestamp)- 17*60*60)

## assign night as number of nights from the start of the study, with all data before the first noon representing night 1
d1$night <- as.numeric( as.Date(time_shift) - start_date + 1 )
d1$night_date <- as.Date( d1$Timestamp - 17*60*60 )
d1$Time <- format(d1$Timestamp, format = "%H:%M:%S")

## show how many nights there are
#nrow( unique( d1[ , c( 'tag', 'night' ) ] ) )

## save a variable denoting the total number of minutes in the day
thresh <- -4.511

mins_in_day <- 60*24 # there are 14 hours between 17:00:00 and 07:00:00 
secs_in_day <- 60*24*60 # there are 14 hours between 17:00:00 and 07:00:00 

missing_mins <- 45*60 ## this is the maximum total number of minutes of data that can be missing from a day and still have that day included in the analysis (for sleep period time and sleep based analyses; i.e. not ave_vedba)

time_gap <- 20*60 ## this is the maximum allowable time gap between two accelerometer bursts (in seconds) that can exist in a noon-to-noon period without removing this noon-to-noon period from the data

mov_window <- 9*60 ## this is the size of the moving window (in minutes) used in calculating the rolling median of the average VeDBA

block_size <- 30*60 ## duration in minutes of the blocks of continuous inactivity that will be considered sleep

gap_size <- 45*60 ## maximum duration between sleep blocks that will be merged

thres <- 0.7 # this is the percentile threshold of the log VeDBA within the 18:00 to 06:30 period used to classify activity vs. inactivity (without multiplying by a multiplier)

waso_block <- 3*10 ## this is the number of consecutive minutes of inactivity needed to classify a period as sleep. A waso_block of 1 means that anytime the value is below the threshold, the baboon in considered sleeping and anytime the value is above the threshold the baboon is considered awake

frag_block <- 2*10 ## this is the number of minutes of waking that need to be consecutive to be considered a wake bout during the night (other epochs of wake that do not meet this criterion will still be considered wake for WASO and wake_bouts, but not frag_wake_bouts)


## create a vector containing the names of each baboon
tag_names <- unique( d1$tag )

## make a copy of d1. We will fill in this new dataframe with information about if the baboon was asleep in each epoch
full_dat <- d1[ d1$Time > "12:00:00" | d1$Time < "12:00:00", ]

full_dat$sleep_per <- NA ## binary indicating whether a row belongs to the sleep period window
full_dat$pot_sleep <- NA ## binary indicating whether the row is below the VeDBA threshold, making it a potential sleep bout. Three or more of these in a row get labeled as sleep bouts
full_dat$sleep_bouts <- NA ## binary indicating whether the row is considered sleep, based on the waso or nap requirements
full_dat$n_bursts <- NA ## the number of bursts collected in a given noon-to-noon period (to be compared to the total number of minutes in the day). This column will indicate whether the data for a given night is insufficient to calculate the sleep period (and thus: onset, waking, SPT, sleep_eff, TST, sleep_bouts -- because this depends on a percentile of bursts' log vedba, WASO, wake_bouts, summed_VeDBA, night_VeDBA_corr, dark_TST, prev_naps, prev_day_sleep)
full_dat$max_time_diff <- NA ## the maximum difference between consecutive fixes in a given noon-to-noon period. With the previous column, this column will indicate whether the data is insufficient to calculate the sleep period (and thus: onset, waking, SPT, sleep_eff, TST, WASO, wake_bouts, summed_VeDBA, night_VeDBA_corr, prev_naps )
full_dat$roll_log_vedba <- NA ## the maximum difference between consecutive fixes in a given noon-to-noon period. With the previous column, this column will indicate whether the data is insufficient to calculate the sleep period (and thus: onset, waking, SPT, sleep_eff, TST, WASO, wake_bouts, summed_VeDBA, night_VeDBA_corr, prev_naps )


## create a vector containing the names of each baboon
#tag_names <- ( unique( full_dat$tag ) )
## for each individual...
for( tag in tag_names ){
  
  ## subset the data to just this individual's data
  #id_dat <- full_dat[ full_dat$tag == tag, ]
  id_dat <- full_dat
  ## create a vector the nights for which this individual has data
  nights <- unique( full_dat$night_date )
  
  #if( length( nights ) > 15 ){ # unhash this if you only want to run sleep classification for individuals that have many nights of data, and thus reliable thresholds (as their thresholds are based on a percentile of their entire 18:00 - 06:30 log VeDBA data, if class_meth is set to 'percentile_thresh' )
  
  # we are going to classify sleep based on a percentile of the individual's rolling log VeDBA for the study period, let's first find the individual's rolling log VeDBA for the full study period, and save the relevant percentile as the threshold
  # create an empty vector to fill with the rolling log VeDBAs from each night
  # full_roll <- c()
  # 
  # # for each night on which this individual has data
  # for( night in nights ){
  #   # subset the individual's data to this night
  #   night_dat <- id_dat[ id_dat$night_date == night, ]
  #   
  #   # take the rolling median of the log VeDBA
  #   roll_log_vedba <- rollmedian( night_dat$log_vedba, mov_window, fill = NA, align = 'center' )
  #   
  #   # add the rolling medians to the vector of the individuals rolling medians for the whole study period
  #   full_roll <- c( full_roll, roll_log_vedba )
  # }
  
  ## determine the threshold activity vs. inactivity threshold based on the percentile, multiplier, and the rolling median just produced
  #thresh <- calculate_threshold(full_roll, thres)
  
  for(night in nights ){
    print(paste(tag, night))
    
    ## subset this individual's data to just that night
    night_dat <- id_dat[ id_dat$night_date == night, ]
    
    ## create empty columns for the sleep period, potential sleep bouts, and sleep bout binary variables
    night_dat$sleep_per <- NA
    night_dat$pot_sleep <- NA
    night_dat$sleep_bouts <- NA
    
    ## save a column of the total number of bursts for that day. This will also make it easier to remove these days from the dataframe later
    night_dat$n_bursts <- nrow( night_dat )
    
    ## sort the timestamps, and book end them with the beginning and end of the night
    #sorted_times <- c( as.POSIXct( paste( as.Date( night, origin = "1970-01-01", tz = 'UTC' ), '18:00:00' ), tz = 'UTC' ), sort( night_dat$Timestamp ), as.POSIXct( paste( as.Date( ( night + 1 ), origin = "1970-01-01", tz = 'UTC' ), '06:30:00' ), tz = 'UTC' ) )
    sorted_times <- night_dat$Timestamp
    
    ## find the time difference in seconds between each burst
    time_diffs <- as.numeric( diff( sorted_times, units = 'secs' ) )
    
    ### find blocks of continuous inactivity
    
    ## take the rolling median of the log VeDBA and save it as a column
    night_dat$roll_log_vedba <- rollmedian( night_dat$log_VeDBA, mov_window, fill = NA, align = 'center' )
    #roll_log_vedba <- night_dat$log_VeDBA
    
    ## find the run length encoding of periods above and below the threshold
    temp <- rle( as.numeric( night_dat$roll_log_vedba < thresh ) ) 
    
    ## mark the rows that are part of runs (i.e. part of chunks that are greater than the block_size of either continuous activity or continuous inactivity )
    sleep_per_runs <- as.numeric( rep( temp$lengths > block_size, times = temp$lengths ) )
    
    ## mark the rows corresponding to sleep bouts. These sleep bouts are runs of inactivity
    sleep_per_sleep_bouts <- as.numeric( night_dat$roll_log_vedba < thresh & sleep_per_runs == 1 )
    
    ## find when sleep bouts start and end
    diffs <- diff( c(0, sleep_per_sleep_bouts ) )
    starts <- which( diffs == 1 ) [ -1 ]
    ends <- which( diffs == -1 )
    
    ## if there are any sleep bouts...
    if( length( which( diffs == 1 ) ) != 0 ){
      
      ## find the duration of the gaps between each sleep bout (the end of one sleep bout and the start of the next)
      gaps <- as.numeric( night_dat$Timestamp [ starts ] - night_dat$Timestamp [ ends[ 1: length( starts ) ] ], units = 'secs' )
      
      ## sleep bouts separated by gaps that are shorter than that specified by gap_size will be merged. Note which of these gaps are shorter than the gap_size
      inds_to_remove <- which( gaps < gap_size ) 
      
      ## if there are NO gaps between sleep bouts that are to be removed...
      if( length( inds_to_remove ) == 0 ){
        
        ## set sleep onset index to be the start of sleep bouts
        onset <- which( diffs == 1 ) 
        
        ## set waking index to be the end of sleep bouts
        wake <- ends
        
      }else{ ## if there ARE gaps between sleep bouts that are to be removed...
        
        ## set sleep onset index to be the start of sleep bouts that do not correspond to the gaps to be removed (because these will be within sleep periods, not a start of a new bout)
        onset <- which( diffs == 1 ) [ - (inds_to_remove + 1) ]
        
        ## set waking index to be the end of sleep bouts that do not correspond to the gaps to be removed
        wake <- ends [ - inds_to_remove ]
        
      }
      
      ## determine which sleep period is the longest
      per_ind <- which.max( as.numeric( night_dat$Timestamp[ wake ] - night_dat$Timestamp[ onset ], units = 'secs' ) )
      
      ## fill in the sleep period data frame with the sleep onset and waking time associated with the longest sleep period in the day (noon to noon)
      
      night_dat$sleep_per <- as.numeric( night_dat$Timestamp >= night_dat$Timestamp[ onset[ per_ind ] ] & night_dat$Timestamp <= night_dat$Timestamp[ wake[ per_ind ] ] )
      
    }else{ ## if there aren't any sleep bouts, record all rows as a 0 in the sleep_per column
      night_dat$sleep_per <- 0
    }
    
    # I was including the possibility of taking another rolling median here of the log VeDBA to determine the epoch-by-epoch sleep-wake classification. But instead I am going to do this with the raw log VeDBA (I was doing that before anyway, as waso_window as set to 1, so I wasn't actually taking a rolling median)
    # ## take the rolling median of the log VeDBA
    
    night_dat$pot_sleep <- as.numeric( night_dat$log_VeDBA < thresh )
    
    ## find the run length encoding of periods above and below the threshold
    temp <- rle( as.numeric( night_dat$log_VeDBA < thresh ) ) 
    
    ## mark the rows that are part of runs (i.e. part of chunks that are greater than the block_size of either continuous activity or continuous inactivity )
    runs <- as.numeric( rep( temp$lengths >= waso_block, times = temp$lengths ) )
    
    ## mark the rows corresponding to sleep bouts. These sleep bouts are runs of inactivity
    sleep_bouts <- as.numeric( night_dat$log_VeDBA < thresh & runs == 1 )
    
    
    ## make which rows are part of runs of inactivity. These are the periods of sleep within and outside of the sleep period
    night_dat$sleep_bouts <- sleep_bouts
    
    ### put the night data back into full_dat
    full_dat[ full_dat$tag == tag & full_dat$night_date == night, ] <- night_dat
  }
  
}

full_dat <- full_dat %>%
  dplyr::group_by(night_date) %>%
  dplyr::mutate(running_number = row_number())

ggplot(full_dat[full_dat$night == 3,], aes(x = running_number, y = roll_log_vedba, group = night)) +
geom_line()
  
filtered_data <- full_dat %>%
    filter(tag == tag) %>%
    dplyr::mutate(timestamp = ifelse(nchar(as.character(Timestamp)) == 10,
                                     paste0(Timestamp, " 00:00:00"),
                                     as.character(Timestamp)),
                  timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d %H:%M:%S"),
                  activity = case_when(
                    sleep_bouts == 1 ~ "Inactive",
                    sleep_bouts == 0 ~ "Active"
                  )) %>%
    dplyr::group_by(night_date) %>%
    dplyr::mutate(time_difference = as.numeric(difftime(timestamp, min(timestamp, na.rm = TRUE), units = "secs")) - 43200) %>%
    ungroup()
  
    ggplot(filtered_data, aes(x = time_difference, y = night_date, fill = activity)) +
    geom_tile() +
    labs(x = tag, y = "") +
    scale_fill_manual(values = c("Inactive" = "white", "Active" = "gray80")) +
    scale_x_continuous(breaks = seq(-43200, 43200, by = 14400), labels = time_labels) +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  


study_nights <- min( d1$night ):max( d1$night )
study_night_dates <- as.Date( min( d1$night_date ):max( d1$night_date ), origin = '1970-01-01' )
sleep_per <- data.frame( tag = rep( unique( d1$tag ), each = length( study_nights ) ), night = rep( study_nights, times = length( tag ) ), night_date = rep( study_night_dates, times = length( tag ) ), total_pot_sleep = NA, total_sleep_bouts = NA, onset = NA, waking = NA, SPT = NA, WASO = NA, TST = NA, sleep_eff = NA, wake_bouts = NA, frag_wake_bouts = NA, summed_VeDBA = NA, night_VeDBA_corr = NA, ave_vedba = NA, dark_pot_sleep = NA, dark_ave_vedba = NA, max_time_diff = NA, n_bursts= NA )

## create empty vectors for the durations of sleep and wake bouts. We will fill these in to see if the distributions of the durations of these bouts later
sleep_durs <- c()
wake_durs <- c() 

## for each individual...
for( tag in tag_names ){
  
  ## subset the data to just this individual's data
  id_dat <- full_dat[ full_dat$tag == tag, ]

  ## create a vector the nights for which this individual has data
  nights <- unique( id_dat$night_date )
  
  ## for each night on which this individual has data
  for( n in 1:length(nights) ){
    print(paste(tag, night))
    night <- nights[n]
    
    ## subset this individual's data to just that night
    night_dat <- id_dat[ id_dat$night_date == night, ]
    
    ## should already be in order, but just in case
    night_dat <- night_dat[ order( night_dat$Timestamp), ]
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night_date == night, ]$n_bursts <- unique( night_dat$n_bursts )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night_date == night, ]$max_time_diff <- unique( night_dat$max_time_diff )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night_date == night, ]$total_pot_sleep <- sum( night_dat$pot_sleep )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night_date == night, ]$total_sleep_bouts <- sum( night_dat$sleep_bouts )
    
    sleep_per[ sleep_per$tag == tag & sleep_per$night_date == night, ]$ave_vedba <- mean( night_dat$log_VeDBA, na.rm = T )
    
    # Get the sun times for the specified location and date
    #sun_times <- getSunlightTimes(date = night, lat = site_lat, lon = site_lon, keep = c("nightEnd", "night"))
    # Extract dark_start and dark_end
    #dark_start <- str_split_fixed(sun_times$night, " ", 2)[,2] 
    #dark_end <- str_split_fixed(sun_times$nightEnd, " ", 2)[,2]
    
    #sleep_per[ sleep_per$tag == tag & sleep_per$night_date == night, ]$dark_pot_sleep <- sum( night_dat$pot_sleep[ night_dat$local_time > dark_start | night_dat$local_time < dark_end ] )
    #sleep_per[ sleep_per$tag == tag & sleep_per$night_date == night, ]$dark_ave_vedba <- mean( night_dat$log_vedba[ night_dat$local_time > dark_start | night_dat$local_time < dark_end ] )
    
    SPT_dat <- night_dat[ night_dat$sleep_per == 1, ]
    
    if( nrow( SPT_dat ) > 0 ){
      
      onset <- min( SPT_dat$Timestamp )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night_date == night, ]$onset <- as.POSIXct(onset, origin = "1970-01-01")
      
      waking <- max( SPT_dat$Timestamp ) 
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night_date == night, ]$waking <- as.POSIXct(waking, origin = "1970-01-01")
      
      SPT <- as.numeric( waking - onset, units = 'mins' ) + 1
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night_date == night, ]$SPT <- SPT
      
      WASO <- sum( SPT_dat$sleep_bouts == 0 )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night_date == night, ]$WASO <- WASO
      
      TST <- sum( SPT_dat$sleep_bouts == 1 )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night_date == night, ]$TST <- TST
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night_date == night, ]$sleep_eff <- TST/ nrow( SPT_dat )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night_date == night, ]$summed_VeDBA <- sum( SPT_dat$log_VeDBA )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night_date == night, ]$night_VeDBA_corr <- sum( SPT_dat$log_VeDBA ) / SPT
      
      temp <- rle( SPT_dat$sleep_bouts )
      
      runs <- as.numeric( rep( temp$lengths >= frag_block, times = temp$lengths ) )
      
      frag_wake_bouts <- as.numeric( SPT_dat$sleep_bouts == 0 & runs == 1 )
      
      diffs <- diff( c( 1, frag_wake_bouts ) )
      
      sleep_per[ sleep_per$tag == tag & sleep_per$night_date == night, ]$frag_wake_bouts <- sum( diffs == 1 )
      
      ## find the distinct sleep bouts (i.e. epochs of sleep separated by waking)
      diffs <- diff( c( 0, SPT_dat$sleep_bouts ) )
      
      ## save the number of distinct wake bouts
      sleep_per[ sleep_per$tag == tag & sleep_per$night_date == night, ]$wake_bouts <- sum( diffs == -1 )
      
      ## find durations of sleep and wake bouts
      temp <- rle( SPT_dat$sleep_bouts )
      
      ## add the duration of sleep bouts to the sleep bout duration vector
      sleep_durs <- c( sleep_durs, temp$lengths[ temp$values == 1 ] )
      ## add the duration of wake bouts to the wake bout duration vector
      wake_durs <- c( wake_durs, temp$lengths[ temp$values == 0 ] )
      
      
    }
  }
}

## reformat sleep timestamp
sleep_per$onset <- as.POSIXct( sleep_per$onset, origin = "1970-01-01 00:00:00", tz = "America/Panama" )

## reformat waking timestamp
sleep_per$waking <- as.POSIXct( sleep_per$waking, origin = "1970-01-01 00:00:00", tz = "America/Panama" )

## make columns for just the time part of the sleep onset and waking timestamps
sleep_per$onset_time <- as_hms( sleep_per$onset )
sleep_per$waking_time <- as_hms( sleep_per$waking )

write.csv( full_dat, paste0( 'data/full_dat.csv' ), row.names = F )
write.csv( sleep_per, paste0( 'data/sleep_per.csv' ), row.names = F )

