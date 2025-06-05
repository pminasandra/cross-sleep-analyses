library(dplyr)
library(ggplot2)
library(zoo)
library(arrow)
library(lubridate)
library(tools)
library(viridis)
library(hms)

# Set paths
parquet_path <- "/mnt/EAS_shared/cross_sleep/working/Data/VeDBA_burst/"
output_dir_data <- "/mnt/EAS_shared/cross_sleep/working/Figures/data_collection/"
output_dir_night <- "/mnt/EAS_shared/cross_sleep/working/Figures/rolling_median/"

# List all .parquet files in the folder
parquet_files <- list.files(parquet_path, pattern = "\\.parquet$", full.names = TRUE)

# Define moving window values
mov_windows <- c(1, 300, 600, 1200, 3600) # in seconds - cont ?
mov_windows <- c(1, 5, 10, 20, 60) # in minutes - burst ?

# Loop through each file
for (file in parquet_files[10:length(parquet_files)]) {
  
  # Extract tag name from file name
  tag_name <- file_path_sans_ext(basename(file))
  
  # Read parquet file
  d1 <- read_parquet(file)
  d1 <- as.data.frame(d1)
  
  if ("timestamp" %in% names(d1)) {
    names(d1)[names(d1) == "timestamp"] <- "Timestamp"
  }
  
  if ("logvedba" %in% names(d1)) {
    names(d1)[names(d1) == "logvedba"] <- "log_vedba"
  }
  
  ## prep time
  time_factor <- 12*60*60
  time_shift <- d1$Timestamp - time_factor
  
  ## save the date of the first night of the study (the date of the night is always the date of the evening at the beginning of that night; for example the first night of the study is 2012-07-31, although the data starts on 2012-08-01, because the data on that first morning is still technically part of the data for the previous night, as a night is noon to noon)
  start_date <- as.Date(min(d1$Timestamp)- time_factor)
  
  ## assign night as number of nights from the start of the study, with all data before the first noon representing night 1
  d1$night <- as.numeric( as.Date(time_shift) - start_date + 1 )
  d1$night_date <- as.Date( d1$Timestamp - time_factor)
  d1$Time <- format(d1$Timestamp, format = "%H:%M:%S")
  d1$Time <- as_hms(d1$Timestamp)
  unique_nights <- unique(d1$night)
  
  print(paste(tag_name,length(unique_nights), sep = "_"))
  
  # p <- ggplot(d1, aes(x = Time, y = night_date)) + 
  #   geom_point(shape = ".", alpha = 0.3) +
  #   theme_minimal() +
  #   theme(
  #     panel.background = element_rect(fill = "white", color = NA),
  #     plot.background = element_rect(fill = "white", color = NA),
  #     panel.grid = element_blank(),
  #     axis.line = element_blank(),
  #     axis.ticks = element_blank(),
  #     axis.text = element_text(size = 8),
  #     plot.title = element_text(size = 10)
  #   ) 
  
  #file_name <- paste0(tag_name, ".jpg")
  #ggsave(filename = file.path(output_dir_data, file_name), plot = p, width = 10, height = 4)  
    
  # Loop over each night
  night_sample <- sample(unique_nights[2:length(unique_nights)-1],min(5, length(unique_nights)-2),replace = TRUE)# unique_nights[2:length(unique_nights)-1]
  for (night_num in night_sample) {
    
    night_data <- d1 %>%
      dplyr::filter(night == night_num) %>%
      dplyr::arrange(Timestamp)
    
    # Skip if too few rows
    if (nrow(night_data) < max(mov_windows)) next
    
    all_rolls <- data.frame()
    
    for (mw in mov_windows) {
      temp <- night_data %>%
        dplyr::mutate(
          mov_window = mw,
          roll_log_vedba = zoo::rollapply(log_vedba, width = mw, FUN = median, na.rm = TRUE, fill = NA, align = "center")
          #roll_log_vedba = zoo::rollmedian(log_vedba, mw, fill = NA, align = "center")
        )
      
      all_rolls <- bind_rows(all_rolls, temp)
    }
    
    alpha_vals <- seq(0.3, 1, length.out = length(mov_windows))
    names(alpha_vals) <- as.character(mov_windows)
    
    # Plot - mov_window overlays per night
    p <- ggplot(all_rolls, aes(x = Timestamp, y = roll_log_vedba, group = factor(mov_window))) +
      geom_line(aes(color = factor(mov_window), alpha = factor(mov_window)), linewidth = 0.6) +
      scale_color_viridis_d(option = "D", name = "mov_window") +  # Discrete viridis
      scale_alpha_manual(values = alpha_vals, guide = "none") +   # Discrete alpha
      # scale_color_manual(
      #   name = "mov_window",
      #   values = gray.colors(length(mov_windows), start = 0.9, end = 0.1),
      #   labels = mov_windows
      # ) +
      geom_line(data = night_data, aes(x = Timestamp, y = log_vedba), inherit.aes = FALSE, color = "lightgrey", alpha = .1,size = 0.1) +
      theme_minimal() +
      labs(title = paste(tag_name, "- night:", night_num, "- rolling median"),
           x = "", y = "log_vedba") +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 10)
      )    
    # File name
    file_name <- paste0(tag_name, "_night-", night_num, "_all_mov_windows.jpg")
    ggsave(filename = file.path(output_dir_night, file_name), plot = p, width = 10, height = 4)
  }
  
}


