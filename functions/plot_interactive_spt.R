#' Create an interactive rolling log VeDBA plot with paired onset/waking lines and legend
#'
#' @param classified_sleep Data frame with at least columns: night, night_date, Timestamp, roll_log_vedba
#' @param sleep_metrics_long Data frame with sleep metrics
#' @param night_id The night number to plot
#' @param folder Folder to save the HTML plot to
#' @param file_prefix 
#' @return Invisibly returns the plotly object
#' @import dplyr
#' @import plotly
#' @importFrom htmlwidgets saveWidget
interactive_spt_plot <- function(
    classified_sleep,
    sleep_metrics_long,
    night_id = 2,
    folder,
    file_prefix = "rolling_log_VeDBA_night"
) {
  # Ensure required packages are loaded
  requireNamespace("dplyr")
  requireNamespace("plotly")
  requireNamespace("htmlwidgets")
  
  # Prepare your data
  df_night <- classified_sleep[classified_sleep$night == night_id,]
  night_dates <- unique(df_night$night_date)
  metrics_night <- sleep_metrics_long %>%
    dplyr::filter(night_date %in% night_dates)
  
  # Pair onset and waking by row
  pairs_df <- metrics_night %>%
    dplyr::filter(!is.na(onset) & !is.na(waking)) %>%
    dplyr::mutate(
      onset_time = as.POSIXct(onset, tz = "UTC"),
      waking_time = as.POSIXct(waking, tz = "UTC")
    )
  
  # Tooltips for onset and waking
  onset_text <- apply(pairs_df, 1, function(row) paste0(
    "<b>Onset:</b> ", row["onset"], "<br>",
    "<b>Night Date:</b> ", row["night_date"], "<br>",
    "<b>SPT:</b> ", row["SPT"], "<br>",
    "<b>WASO:</b> ", row["WASO"], "<br>",
    "<b>TST:</b> ", row["TST"], "<br>",
    "<b>sleep_eff:</b> ", row["sleep_eff"], "<br>",
    "<b>wake_bouts:</b> ", row["wake_bouts"], "<br>",
    "<b>frag_wake_bouts:</b> ", row["frag_wake_bouts"], "<br>",
    "<b>summed_VeDBA:</b> ", row["summed_VeDBA"], "<br>",
    "<b>block_size:</b> ", row["block_size"], "<br>",
    "<b>gap_size:</b> ", row["gap_size"], "<br>",
    "<b>mov_window:</b> ", row["mov_window"]
  ))
  waking_text <- apply(pairs_df, 1, function(row) paste0(
    "<b>Waking:</b> ", row["waking"], "<br>",
    "<b>Night Date:</b> ", row["night_date"], "<br>",
    "<b>SPT:</b> ", row["SPT"], "<br>",
    "<b>WASO:</b> ", row["WASO"], "<br>",
    "<b>TST:</b> ", row["TST"], "<br>",
    "<b>sleep_eff:</b> ", row["sleep_eff"], "<br>",
    "<b>wake_bouts:</b> ", row["wake_bouts"], "<br>",
    "<b>frag_wake_bouts:</b> ", row["frag_wake_bouts"], "<br>",
    "<b>summed_VeDBA:</b> ", row["summed_VeDBA"], "<br>",
    "<b>block_size:</b> ", row["block_size"], "<br>",
    "<b>gap_size:</b> ", row["gap_size"], "<br>",
    "<b>mov_window:</b> ", row["mov_window"]
  ))
  
  # Base plot
  p <- plotly::plot_ly() %>%
    plotly::add_lines(
      data = df_night,
      x = ~Timestamp,
      y = ~roll_log_vedba,
      name = "roll_log_vedba",
      line = list(color = "black", width = 0.5),
      opacity = 0.7,
      hoverinfo = "x+y",
      showlegend = FALSE
    )
  
  # Add each pair as a group with interactive legend
  for (i in seq_len(nrow(pairs_df))) {
    legend_label <- paste0("gap ", pairs_df$gap_size[i], ", window ", pairs_df$mov_window[i])
    legend_grp <- paste0("pair_", i)
    
    # Onset (red, legend shown)
    p <- p %>% plotly::add_segments(
      x = pairs_df$onset_time[i],
      xend = pairs_df$onset_time[i],
      y = min(df_night$roll_log_vedba, na.rm = TRUE),
      yend = max(df_night$roll_log_vedba, na.rm = TRUE),
      name = legend_label,
      legendgroup = legend_grp,
      showlegend = TRUE,
      line = list(color = "red", width = 2, dash = "dash"),
      hoverinfo = "text",
      text = onset_text[i]
    )
    # Waking (blue, legend not shown but in same group)
    p <- p %>% plotly::add_segments(
      x = pairs_df$waking_time[i],
      xend = pairs_df$waking_time[i],
      y = min(df_night$roll_log_vedba, na.rm = TRUE),
      yend = max(df_night$roll_log_vedba, na.rm = TRUE),
      name = legend_label,
      legendgroup = legend_grp,
      showlegend = FALSE,
      line = list(color = "blue", width = 2, dash = "dash"),
      hoverinfo = "text",
      text = waking_text[i]
    )
  }
  
  p <- p %>% plotly::layout(
    xaxis = list(title = "Timestamp"),
    yaxis = list(title = "roll_log_vedba"),
    hovermode = "closest"
  )
  
  file_path <- file.path(folder, paste0(file_prefix, night_id, ".html"))
  htmlwidgets::saveWidget(p, file = file_path)
  
  message("Plot saved to: ", file_path)
  invisible(p)
}