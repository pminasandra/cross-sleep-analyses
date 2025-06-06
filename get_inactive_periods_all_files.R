#Get inactivity periods across all files

#libraries
library(arrow)

#----PARAMS----
time_win <- 60*60 #time window in seconds 
manual_set_thresh <- -2.5 #logvedba threshold to use (if NULL, this will be pulled from metadata table)
start_hour_UTC <- 11 #hour (in UTC) to use as the start of the window (should be somewhere in the middle of the 'active period')
make_plots <- T #T or F, whether to make plots
species_to_use <- 'meerkat' #can be hyena, meerkat, spidermonkey, or hyrax -or can be NULL and all files in the folder are used

#directories and files
datadir <- '~/EAS_shared/cross_sleep/working/Data/'
plotdir <- '~/EAS_shared/cross_sleep/working/Figures/inactive_periods/'
metadata_file <- 'cross_sleep_metadata_2025-05-20.csv'

#source needed function
source('functions/get_inactive_periods.R')

#----MAIN----

#get all files for the target species
files <- list.files(paste0(datadir,'VeDBA/'))

if(!is.null(species_to_use)){
  files <- files[grep(species_to_use,files)]
}

#loop over files
print('getting inactivity periods over all files...')
inactive_periods <- data.frame()
for(f in 1:length(files)){
  
  filename <- files[f]
  print(filename)
  
  #get info from the filename
  splitname <- strsplit(basename(filename),split='_')[[1]]
  species <- splitname[1]
  deploy_id <- splitname[2]
  ind_id <- splitname[3]
  
  #get threshold from metadata file if needed (or set manually)
  if(is.null(manual_set_thresh)){
    metadata <- read.csv(metadata_file)
    row <- which(metadata$individual_ID==ind_id & metadata$deployment_ID==deploy_id & metadata$species==species)
    metadata_ind <- metadata[row,]
    thresh <- metadata_ind$th1
  } else{
    thresh <- manual_set_thresh
  }
  
  #Load in file
  dat <- arrow::read_parquet(paste0(datadir,'VeDBA/',filename))
  
  #replace infinity values with values below threshold (fix in original files later)
  dat$logvedba[which(is.infinite(dat$logvedba))] <- thresh - 1
  
  plot_savepath <- paste0(plotdir, gsub('[.]parquet','', filename),'_inactive_periods.png')
  
  inactive_periods_file <- get_inactive_periods(dat = dat, 
                                                time_win = time_win,
                                                logvedba_thresh = thresh,
                                                make_plot = make_plots,
                                                plot_savepath = plot_savepath)
  
  inactive_periods_file$filename <- filename
  inactive_periods_file$species <- species
  inactive_periods_file$deploy_id <- deploy_id
  inactive_periods_file$ind_id <- ind_id
  inactive_periods_file$time_win <- time_win
  inactive_periods_file$logvedba_thresh <- thresh
  
  inactive_periods <- rbind(inactive_periods, inactive_periods_file)
  
}

if(is.null(species)){
  savepath <- paste0(datadir,'inactive_periods/inactive_periods_allspecies.RData')
} else{
  savepath <- paste0(datadir,'inactive_periods/inactive_periods_',species,'.RData')
}

save(list = 'inactive_periods', file = savepath)


