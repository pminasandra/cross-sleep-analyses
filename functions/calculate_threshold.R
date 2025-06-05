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
