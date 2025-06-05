## function for plotting times from noon to noon. It will make correct 12:00 - 24:00 to be 00:00 - 12:00 and 00:00 - 12:00 to be 12:00 to 24:00
calculate_threshold <- function(vedba_table, thres_value = 0.1) {
  # Remove NA values
  vedba_table <- vedba_table[!is.na(vedba_table$logvedba), ]
  
  # Standardize the data
  #log_vedba_vector <- scale(log_vedba_vector)
  
  # Check if the data is valid
  if (nrow(vedba_table) < 2) {
    stop("Insufficient data to calculate threshold.")
  }
  
  # Fit a Gaussian mixture model with 2 components
  fit <- tryCatch({
    normalmixEM(vedba_table$logvedba, k = 2, arbmean = TRUE, arbvar = TRUE)
  }, error = function(e) {
    stop("Error fitting model: ", e$message)
  })
  
  # Identify the leftmost distribution
  sorted_indices <- order(fit$mu)
  leftmost_index <- sorted_indices[1]
  
  # Calculate the posterior probabilities for the leftmost component
  posterior_probs_leftmost <- fit$posterior[, leftmost_index]
  
  # Calculate the leftmost point where the probability of being from the left distribution is under 25%
  threshold <- min(vedba_table$logvedba[
    vedba_table$logvedba > max(vedba_table$logvedba[posterior_probs_leftmost > thres_value]) &
      posterior_probs_leftmost < thres_value])
  
  vedba_table$posterior_prob_left <- posterior_probs_leftmost

  return(list(threshold = threshold, vedba_table_posterior = vedba_table))

}


########## Test

library(arrow)
library(mixtools)
library(ggplot2)
library(stringr)
library(dplyr)
library(cowplot)



metadata <- read.csv("/mnt/EAS_shared/cross_sleep/working/Data/cross_sleep_metadata_2025-05-20.csv")
vedbadata_path <- "/mnt/EAS_shared/cross_sleep/working/Data/VeDBA/"
fig_output_folder <- "/mnt/EAS_shared/cross_sleep/working/Figures/thresholds/"

# Extract metadata from filename
#first, just pick one file

parquet_path<- "/mnt/EAS_shared/cross_sleep/working/Data/VeDBA/hyena_2016_BORA_vedba_standard.parquet"
parquet_filename <- basename(parquet_path)
filename_core <- sub("\\.parquet$", "", parquet_filename)
parts <- strsplit(filename_core, "_")[[1]]
curr_species <- parts[1]
curr_deployment_ID <- parts[2]
curr_individual_ID <- parts[3]
curr_burst_int <- as.numeric(str_extract(parquet_filename, "(?<=burst_int)\\d+"))
curr_len_val   <- as.numeric(str_extract(parquet_filename, "(?<=_len)\\d+"))


# List all files in the folder
all_files <- list.files(vedbadata_path, full.names = TRUE)

# Filter files that are corresponding to that individual
matching_files <- all_files[
  grepl(curr_individual_ID, all_files) & grepl(curr_species, all_files)
]

# sort the file by increasing burst_int value

# Extract burst_int (NA if not found)
burst_ints <- as.numeric(str_extract(matching_files, "(?<=burst_int)\\d+"))
# Create logical index: TRUE if burst_int is NA (non-matching files)
non_matching <- is.na(burst_ints)
# Sort: non-matching first, then by burst_int ascending
matching_files <- c(
  matching_files[non_matching],
  matching_files[!non_matching][order(burst_ints[!non_matching])]
)


# Create the empty data frame with predefined column types
n_rows <- length(matching_files)
Auto_th <- data.frame(
  species       = rep(NA_character_, n_rows),
  deployment_ID = rep(NA_character_, n_rows),
  individual_ID = rep(NA_character_, n_rows),
  burst_int     = rep(NA_real_, n_rows),
  burst_length  = rep(NA_real_, n_rows),
  manual_th     = rep(NA_real_, n_rows),
  automatic_th  = rep(NA_real_, n_rows),
  stringsAsFactors = FALSE
)

#just read this one first to define axis limits
d1 <- read_parquet(parquet_path)
x_limits = c(floor(min(d1$logvedba, na.rm = T)-1), ceiling(max(d1$logvedba, na.rm = T)+1))

for (i in 1:length(matching_files)) {

  d1 <- read_parquet(matching_files[i])
  # Extract metadata from filename
  parquet_filename <- basename(matching_files[i])
  filename_core <- sub("\\.parquet$", "", parquet_filename)
  parts <- strsplit(filename_core, "_")[[1]]
  curr_species <- parts[1]
  curr_deployment_ID <- parts[2]
  curr_individual_ID <- parts[3]
  curr_burst_int <- as.numeric(str_extract(parquet_filename, "(?<=burst_int)\\d+"))
  curr_burst_len   <- as.numeric(str_extract(parquet_filename, "(?<=_len)\\d+"))
  

  manual_threshold <- metadata %>%
    filter(deployment_ID == curr_deployment_ID & individual_ID == curr_individual_ID) %>%
    pull(th1)
  
  results <- calculate_threshold(d1)
  attempt = 1
  while (is.infinite(results$threshold)) {
    print("Returning infinite threshold. Trying again")
    attempt = attempt+1
    print(paste0("#",attempt," attempt"))
    results <- calculate_threshold(d1)
  } 
  
  Auto_th$species[i] = curr_species
  Auto_th$deployment_ID[i] = curr_deployment_ID
  Auto_th$individual_ID[i] = curr_individual_ID
  Auto_th$burst_int[i] = curr_burst_int
  Auto_th$burst_length[i] = curr_burst_len
  Auto_th$manual_th[i] = manual_threshold
  Auto_th$automatic_th[i] = results$threshold
  
  # Create the plot
  plot[[i]] <- ggplot(results$vedba_table_posterior, aes(x = logvedba)) +
    geom_histogram(aes(y = ..density..), bins = 100, fill = "blue", alpha = 0.5) +
    geom_point(aes(y = posterior_prob_left, color = "Leftmost"), size = 1, alpha = 0.7) +
    labs(title = filename_core,
         x = "log(VeDBA)",
         y = "Density / Probability") + xlim(x_limits)+
    theme_minimal() +
    geom_vline(xintercept = results$threshold, linetype = "dotted", color = "red", size = 2) +
    geom_vline(xintercept = manual_threshold, linetype = "dotted", color = "blue", size = 2)+
    theme(legend.position = "none")

}

grid_plot <- do.call(plot_grid, c(plot[1:9], ncol = 3, labels = "AUTO", align = "hv"))

ggsave(paste0(fig_output_folder,
              curr_species ,"_",curr_deployment_ID ,"_",curr_individual_ID, "_vedba_disttrib.png"), 
       plot = grid_plot, width = 24, height = 16, dpi = 100, bg = "white")

# Create a column with the difference
Auto_th$threshold_diff <- Auto_th$automatic_th - Auto_th$manual_th

# calculate Coefficient of variation
cv_percent <- (sd(Auto_th$automatic_th) / abs(mean(Auto_th$automatic_th))) * 100
cv_label <- paste0("CV automatic ths = ", round(cv_percent, 2), "%")
cv_label

manual_linear <- exp(Auto_th$manual_th)
auto_linear   <- exp(Auto_th$automatic_th)
percent_diff <- (auto_linear - manual_linear) / manual_linear * 100
mean_percent_diff <- mean(percent_diff, na.rm = TRUE)

mean_percent_diff
mean(Auto_th$threshold_diff)


# Plot the correlation
plot <- ggplot(Auto_th, aes(x = burst_int, y = threshold_diff)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(
    title = "Correlation between Burst Interval and Threshold Difference",
    x = "Burst Interval",
    y = "Automatic - Manual Threshold"
  ) +
  annotate("text", x = -Inf, y = Inf, label = cv_label, hjust = -0, vjust = 1.5, size = 4) +
  theme_minimal()+
  geom_hline(yintercept = Auto_th$automatic_th[is.na(Auto_th$burst_int)]-Auto_th$manual_th[1], linetype = "dotted", color = "red", size = 2) +
  theme_minimal()

ggsave(paste0(fig_output_folder, Auto_th$species[1],"_", 
              Auto_th$deployment_ID[1],"_", Auto_th$individual_ID[1],"_th-int.png"), 
       plot = plot, width = 8, height = 5, dpi = 100, bg = "white")


plot <- ggplot(Auto_th, aes(x = burst_int, y = automatic_th)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(
    title = "",
    x = "Burst Interval",
    y = "Thresholds"
  ) +
  geom_hline(yintercept = Auto_th$manual_th[1], linetype = "dotted", color = "blue", size = 2) +
  geom_hline(yintercept = Auto_th$automatic_th[is.na(Auto_th$burst_int)], linetype = "dotted", color = "red", size = 2) +
  theme_minimal()+ylim(x_limits)

ggsave(paste0(fig_output_folder, Auto_th$species[1],"_", 
              Auto_th$deployment_ID[1],"_", Auto_th$individual_ID[1],"_th-int2.png"), 
       plot = plot, width = 8, height = 5, dpi = 100, bg = "white")

Auto_th

