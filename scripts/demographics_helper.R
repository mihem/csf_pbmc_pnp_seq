########################################################
# Helper functions for demographics.R
########################################################

# Create boxplot with statistics ----
create_boxplot <- function(
  data,
  x_var,
  y_var,
  group_var,
  title,
  color_palette,
  geom_type = "point",
  filter_na = TRUE
) {
  # Filter NA values if requested
  plot_data <- if (filter_na) {
    data |> dplyr::filter(!is.na(.data[[y_var]]))
  } else {
    data
  }

  # Calculate statistics
  stats <- scMisc:::compStat(
    x_var = y_var,
    group = group_var,
    data = plot_data,
    paired = FALSE
  )

  # Create base plot
  p <- plot_data |>
    ggplot(aes(
      x = .data[[x_var]],
      y = .data[[y_var]],
      fill = .data[[x_var]]
    )) +
    geom_boxplot()

  # Add points based on geom_type
  if (geom_type == "point") {
    p <- p + geom_point()
  } else if (geom_type == "jitter") {
    p <- p + geom_jitter(height = 0, width = 0.3)
  }

  # Add significance annotations
  p <- p +
    ggsignif::geom_signif(
      comparisons = stats$comparisons,
      annotation = stats$annotation,
      textsize = 5,
      step_increase = 0.05,
      vjust = 0.7
    ) +
    theme_bw() +
    scale_fill_manual(values = color_palette) +
    xlab("") +
    ylab("") +
    ggtitle(title) +
    theme(legend.position = "none")

  return(p)
}

# Create barplot ----
create_barplot <- function(
  data,
  x_var,
  fill_var,
  title,
  color_palette
) {
  p <- data |>
    ggplot(aes(x = .data[[x_var]], fill = .data[[fill_var]])) +
    geom_bar() +
    theme_bw() +
    scale_fill_manual(values = color_palette) +
    xlab("") +
    ylab("") +
    ggtitle(title)

  return(p)
}

# Batch create and save boxplots ----
create_and_save_boxplots <- function(
  data,
  configs,
  x_var,
  group_var,
  color_palette,
  output_dir
) {
  for (name in names(configs)) {
    config <- configs[[name]]
    plot <- create_boxplot(
      data = data,
      x_var = x_var,
      y_var = config$y_var,
      group_var = group_var,
      title = config$title,
      color_palette = color_palette,
      geom_type = config$geom_type
    )
    
    ggsave(
      file.path(output_dir, paste0("boxplot_", name, ".pdf")),
      plot = plot,
      width = config$width,
      height = config$height
    )
  }
}

# Batch create and save barplots ----
create_and_save_barplots <- function(
  data,
  configs,
  x_var,
  output_dir
) {
  for (name in names(configs)) {
    config <- configs[[name]]
    plot <- create_barplot(
      data = data,
      x_var = x_var,
      fill_var = config$fill_var,
      title = config$title,
      color_palette = config$color_palette
    )
    
    ggsave(
      file.path(output_dir, paste0("barplot_", name, ".pdf")),
      plot = plot,
      width = config$width,
      height = config$height
    )
  }
}

# Create overview table from patient data ----
create_overview_table <- function(patient_data) {
  overview_table <- patient_data |>
    dplyr::mutate(sex_cat = if_else(sex == "male", 1, 0)) |>
    dplyr::group_by(diagnosis) |>
    dplyr::summarize(
      samples = n(),
      patients = n_distinct(orbis_id),
      age = mean(age, na.rm = TRUE),
      female = (1 - mean(sex_cat, na.rm = TRUE)) * 100
    )
  
  return(overview_table)
}

