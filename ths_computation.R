# =============================================================================
# Title       : Compute and Compare Automatic Behavioral Thresholds with Manual Metadata
# Description : This script sources all R functions from a specified project directory,
#               loads the most recent metadata file, runs a behavioral thresholding function
#               (`compare_auto_manual_thresholds`) on a specified .parquet file containing
#               VeDBA (Vectorial Dynamic Body Acceleration) data, and compares the automatic
#               thresholds with those manually defined in the metadata. It generates visual
#               outputs, updates the metadata table, saves results, and ensures reproducibility
#               by timestamping output files.
#
# Author      : 
# Date        : 
# License     : 
# =============================================================================

# ==== Load and source external R function files =================================

# Path to folder containing necessary custom R functions
func_folder <- "/mnt/EAS_ind/mratsimbazafindranah/code/cross-sleep-analyses/functions"
# List all .R script files in the folder
func_files <- list.files(func_folder, pattern = "\\.R$", full.names = TRUE)
# Source each function file to load them into the current R session
sapply(func_files, source)

# ==== Load metadata (selecting most recent version based on date in filename) ===

# Define folder where metadata CSV files are stored
metadata_folder <- "/mnt/EAS_shared/cross_sleep/working/Data/"

# List all metadata files matching naming convention with a date suffix
meta_files <- list.files(metadata_folder, 
                         pattern = "^cross_sleep_metadata_\\d{4}-\\d{2}-\\d{2}\\.csv$", 
                         full.names = TRUE)
# Extract the date from each filename using regex and convert to Date format
dates <- as.Date(stringr::str_extract(meta_files, "\\d{4}-\\d{2}-\\d{2}"))
# Identify the most recent metadata file based on the date extracted
metadata_path <- meta_files[which.max(dates)]
metadata_path

# ==== Define output folders for figures and threshold tables ===================

fig_output_folder <- "/mnt/EAS_shared/cross_sleep/working/Figures/thresholds/"
ths_table_output_folder <- "/mnt/EAS_shared/cross_sleep/working/Data/thresholds/"

# Read the selected metadata CSV into a data frame
metadata <- read.csv(metadata_path)

# ==== Perform automatic threshold computation and compare to manual ============

# Path to VeDBA data (parquet format); continuous or lowest burst interval preffered
parquet_path <- "/mnt/EAS_shared/cross_sleep/working/Data/VeDBA/baboon_4707657037_14556_vedba_standard.parquet"

# Run custom function to compute thresholds using GMM (Gaussian Mixture Model)
# `thres_value` is the lower quantile used as a baseline; `k` is number of GMM components
result <- compare_auto_manual_thresholds(parquet_path, thres_value = 0.1, k = 2,
                                         metadata, fig_output_folder)

# ==== Extract and integrate automatic threshold from continuous (or lowest burst) ===

# Keep only the first row from result table (typically the continuous data)
first_result <- result %>%
  slice(1) %>%
  select(species, deployment_ID, individual_ID, automatic_th, automatic_th2)

# Add columns if they don’t exist yet
{
  if (!"automatic_th" %in% names(metadata)) {
    metadata$automatic_th <- NA_real_
  }
  if (!"automatic_th2" %in% names(metadata)) {
    metadata$automatic_th2 <- NA_real_
  }
  
  # Option 1: Manually update only the matching row in metadata
  metadata <- metadata %>%
    mutate(
      automatic_th = if_else(
        species == first_result$species &
          deployment_ID == first_result$deployment_ID &
          individual_ID == first_result$individual_ID,
        first_result$automatic_th,
        automatic_th  # keep original
      ),
      automatic_th2 = if_else(
        species == first_result$species &
          deployment_ID == first_result$deployment_ID &
          individual_ID == first_result$individual_ID,
        first_result$automatic_th2,
        automatic_th2
      )
    )
  # Print updated metadata to console for verification
  print(metadata) }

# ==== Write full result table to file ==========================================

# Construct output filename using subject identifiers from result
write.csv(result, file = file.path(ths_table_output_folder, 
                                   paste0(result$species[1], "_", 
                                          result$deployment_ID[1], "_",
                                          result$individual_ID[1], "_ths.csv")), 
          row.names = FALSE)



# ==== Write updated metadata file with today's date =============================

# Get today's date in yyyy-mm-dd format for filename versioning
{date_suffix <- format(Sys.Date(), "%Y-%m-%d")

# Construct new file path using today's date
new_metadata_path <- paste0(dirname(metadata_path), "/cross_sleep_metadata_", date_suffix, ".csv")

# Save the updated metadata as a new CSV
write.csv(metadata, new_metadata_path, row.names = FALSE)

# Print confirmation message with the file path
cat("✅ Metadata successfully written to:\n", new_metadata_path, "\n")
}
