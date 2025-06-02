#Calculate VeDBA from standardized ACC data
acc_to_vedba <- function(df, rolling_mean_width){
  
  present_axes <- c("X", "Y", "Z")
  
  #static component
  df[, (paste0(present_axes, "_static")) := lapply(present_axes, function(axis) {
    frollmean(get(axis), n = rolling_mean_width, align = "center", na.rm = FALSE)
  })]
  
  #dynamic component
  df[, (paste0(present_axes, "_dynamic")) := Map(
    function(raw, stc) raw - stc,
    lapply(present_axes, function(axis) get(axis)),
    .SD
  ), .SDcols = paste0(present_axes, "_static")]
  
  #vedba
  df[, VeDBA := sqrt(Reduce(`+`, lapply(.SD, function(d) d^2))),
     .SDcols = paste0(present_axes, "_dynamic")]
  
  #number of NAs
  df$fracNA <- frollmean(is.na(df[,'X'])[,1], n = rolling_mean_width, align = 'center', na.rm=T)
  
  df <- df[, !c(paste0(present_axes, "_static"), paste0(present_axes, "_dynamic")), with = FALSE]
  
  # Step 4: Downsample to 1Hz
  #df_1s <- df[df$Timestamp_adj == floor_date(df$Timestamp_adj, unit = "second")]
  #cat("Downsampling ratio (raw/1s):", nrow(df) / nrow(df_1s), "\n")

  #get only relevant columns and standardize naming
  df <- df[,c('Timestamp_adj','VeDBA','fracNA')]
  colnames(df) <- c('Timestamp','VeDBA','fracNA')
  
  return(df)
}
