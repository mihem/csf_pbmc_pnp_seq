plotPropeller <- function(data, color, filename, width = 5, height = 5, FDR = NULL, use_permFDP = FALSE, dir_output = ".") {
  # Determine the threshold to use and y-axis variable
  if (use_permFDP) {
    if (!"permFDP_threshold" %in% names(data)) {
      stop("use_permFDP = TRUE but permFDP_threshold column not found in data")
    }
    threshold <- unique(data$permFDP_threshold)[1]
    y_var <- -log10(data$P.Value)
    message(paste0("Using permFDP threshold: ", round(threshold, 4)))
  } else {
    if (is.null(FDR)) {
      stop("FDR must be specified when use_permFDP = FALSE")
    }
    threshold <- FDR
    y_var <- data$FDR_log
  }

  # label only significant points (above p threshold and |log2FC| >= 0.5)
  label_var <- dplyr::if_else(
    y_var >= -log10(threshold) & abs(data$log2ratio) >= 0.5,
    data$cluster,
    NA_character_
  )

  plot <- ggplot(data, aes(x = log2ratio, y = y_var, color = cluster, size = 3, label = label_var)) +
    geom_point(size = 3) +
    scale_color_manual(values = color) +
    theme_classic() +
    ggrepel::geom_text_repel(nudge_y = 0.07, max.overlaps = 20, na.rm = TRUE) +
    geom_hline(yintercept = -log10(threshold), color = "blue", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "red", linetype = "solid") +
    geom_vline(xintercept = -0.5, color = "red", linetype = "dashed") +
    geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") +
    xlab(bquote(~ Log[2] ~ "fold change")) +
    ylab(bquote(~ -Log[10] ~ "adjusted p value")) +
    theme(legend.position = "none")
  ggsave(file.path(dir_output, glue::glue("propeller_{filename}.pdf")), plot, width = width, height = height)
  return(plot)
}
