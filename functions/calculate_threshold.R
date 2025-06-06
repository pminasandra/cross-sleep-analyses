#' @title Calculate Behavioral Thresholds Using Gaussian Mixture Model
#'
#' @description
#' This function fits a Gaussian Mixture Model (GMM) with `k` components to a vector of log-transformed VeDBA (Vectorial Dynamic Body Acceleration) values
#' to identify behavioral thresholds. These thresholds help distinguish between different activity states (e.g., resting vs. active).
#' The function returns threshold values based on posterior probabilities from the fitted model, specifically where posterior drops below a defined threshold.
#'
#' The function also appends the posterior probabilities for each mixture component to the input data frame, useful for visual inspection or further classification.
#'
#' @param vedba_table A data frame containing a column named `logvedba` with log-transformed VeDBA values.
#' @param thres_value A numeric value indicating the minimum posterior probability cutoff to define a threshold. Default is 0.1.
#' @param k Integer specifying the number of components (clusters) to use in the Gaussian Mixture Model. Default is 2.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{threshold}{A numeric vector of threshold values (length = k - 1) indicating transitions between behavioral states.}
#'   \item{vedba_table_posterior}{The input data frame `vedba_table` augmented with additional columns for each posterior probability.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- calculate_threshold(mydata, thres_value = 0.1, k = 2)
#'
#' @author Cross-sleep team, MPIAB
#' @export

library(mixtools)

calculate_threshold <- function(vedba_table, thres_value = 0.1, k = 2) {
  # Remove rows with missing (NA) values in the logvedba column
  vedba_table <- vedba_table[!is.na(vedba_table$logvedba), ]
  
  # If fewer than 2 data points remain, terminate with an error
  if (nrow(vedba_table) < 2) {
    stop("Insufficient data to calculate threshold.")
  }
  
  # Fit a Gaussian Mixture Model with k components to the logvedba values
  fit <- tryCatch({
    normalmixEM(vedba_table$logvedba, k, arbmean = TRUE, arbvar = TRUE)
  }, error = function(e) {
    # Catch fitting errors and stop the function with a descriptive message
    stop("Error fitting model: ", e$message)
  })
  
  # Sort the mixture components by increasing mean value (for logical ordering)
  sorted_indices <- order(fit$mu)
  
  # Reorder posterior probabilities accordingly
  posterior_sorted <- fit$posterior[, sorted_indices]
  
  # Determine how many components were actually returned (could be < k)
  num_components <- ncol(posterior_sorted)
  if (num_components < k) {
    message(sprintf(" ⚠️ Model returned only %d components instead of requested k = %d.", num_components, k))
  }
  
  # Add each posterior probability column to the data frame
  # Useful for visualization or additional downstream analysis
  for (i in 1:num_components) {
    vedba_table[[paste0("posterior_comp", i)]] <- posterior_sorted[, i]
  }
  
  # --- Helper function to compute thresholds ---
  
  # This function identifies the logvedba value after which a component's posterior
  # probability drops below the specified threshold value (thres_value)
  get_tail_threshold <- function(posterior, logvedba) {
    # Identify maximum logvedba with posterior above threshold
    upper <- max(logvedba[posterior > thres_value])
    # Find the smallest logvedba greater than 'upper' where posterior drops below threshold
    tail_thresh <- min(logvedba[logvedba > upper & posterior < thres_value])
    return(tail_thresh)
  }
  
  # This function identifies the logvedba value where two components
  # have approximately equal posterior probabilities (i.e., intersection)
  get_midpoint_threshold <- function(p1, p2, logvedba) {
    # Compute the absolute difference in posterior probabilities
    delta <- abs(p1 - p2)
    # Return the logvedba value with the smallest difference
    return(logvedba[which.min(delta)])
  }
  
  # --- Threshold calculation across components ---
  
  # Initialize vectors to store tail and midpoint thresholds between adjacent components
  tail_thresholds <- numeric(num_components - 1)
  midpoint_thresholds <- numeric(num_components - 1)
  
  # Loop over all adjacent pairs of components (1 vs 2, 2 vs 3, ...)
  for (i in 1:(num_components - 1)) {
    # Compute tail threshold: where posterior of component i drops below thres_value
    tail_thresholds[i] <- get_tail_threshold(
      posterior_sorted[, i],
      vedba_table$logvedba
    )
    
    # Compute midpoint threshold: where components i and i+1 are equally likely
    midpoint_thresholds[i] <- get_midpoint_threshold(
      posterior_sorted[, i],
      posterior_sorted[, i + 1],
      vedba_table$logvedba
    )
  }
  
  # Return both the vector of tail thresholds and the enriched table with posterior probabilities
  return(list(threshold = tail_thresholds, vedba_table_posterior = vedba_table))
}
