#Calculate VeDBA from standardized ACC data
acc_to_vedba <- function(df, rolling_mean_width, group_col = NULL) {
  df <- copy(df)  # <-- prevents modifying the original data.table
  
  present_axes <- c("X", "Y", "Z")
  by_expr <- if (!is.null(group_col)) group_col else NULL
  
  # STATIC
  df[, (paste0(present_axes, "_static")) := lapply(present_axes, function(axis) {
    frollmean(get(axis), n = rolling_mean_width, align = "center", na.rm = TRUE)
  }), by = by_expr]
  
  # DYNAMIC
  df[, (paste0(present_axes, "_dynamic")) := Map(
    function(raw, stc) raw - stc,
    lapply(present_axes, function(axis) get(axis)),
    .SD
  ), by = by_expr, .SDcols = paste0(present_axes, "_static")]
  
  # VeDBA
  df[, VeDBA := sqrt(Reduce(`+`, lapply(.SD, function(d) d^2))),
     by = by_expr,
     .SDcols = paste0(present_axes, "_dynamic")]
  
  # fracNA
  df[, fracNA := frollmean(is.na(get("X")), n = rolling_mean_width, align = 'center', na.rm = TRUE),
     by = by_expr]
  
  # Cleanup
  df[, c(paste0(present_axes, "_static"), paste0(present_axes, "_dynamic")) := NULL]
  
  return(df)
}
