library(arrow)
library(dplyr)
library(ggplot2)
library(stringr)

# Define paths
input_dir <- "/mnt/EAS_shared/cross_sleep/working/Data/ClassifiedSleep/"
output_dir <- "/mnt/EAS_shared/cross_sleep/working/Figures/activity_noon_noon/"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir)

# Get list of all classified sleep parquet files
files <- list.files(input_dir, pattern = ".parquet$", full.names = TRUE)

# Time labels for plotting
time_labels <- c("1200", "1600", "2000", "2400", "0400", "0800", "1200")

for (file in files) {
  # Read parquet file
  parquet_filename <- basename(file)
  filename_core <- sub("\\.parquet$", "", parquet_filename)
  
  classified_sleep <- read_parquet(file) |> as.data.frame()
  
  nights <- unique(classified_sleep$night_date)
  
    filtered_data <- classified_sleep %>%
      dplyr::mutate(
        Timestamp = ifelse(nchar(as.character(Timestamp)) == 10,
                           paste0(Timestamp, " 00:00:00"),
                           as.character(Timestamp)),
        Timestamp = as.POSIXct(Timestamp, format = "%Y-%m-%d %H:%M:%S"),
        activity = dplyr::case_when(
          pot_sleep == 1 ~ "Inactive",
          pot_sleep == 0 ~ "Active"
        )
      ) %>%
      dplyr::group_by(night_date) %>%
      dplyr::mutate(
        time_difference = as.numeric(difftime(Timestamp, min(Timestamp, na.rm = TRUE), units = "secs")) - 43200
      ) %>%
      ungroup()
    
    p <- ggplot(filtered_data, aes(x = time_difference, y = night_date, fill = activity)) +
      geom_tile() +
      scale_fill_manual(values = c("Inactive" = "gray90", "Active" = "gray60")) +
      scale_x_continuous(breaks = seq(-43200, 43200, by = 14400), labels = time_labels) +
      theme_minimal() +
      theme(
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 8),
        plot.title = element_text(size = 10)
      )    
    
    # Save plot
    ggsave(
      filename = file.path(output_dir, paste0(filename_core,".jpg")),
      plot = p, width = 8, height = 6, dpi = 150
    )
  }
