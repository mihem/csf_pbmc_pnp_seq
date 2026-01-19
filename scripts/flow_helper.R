ggboxplotFun <- function(var, data, group, stats, cols) {
  # plot boxplot
  ggplot(data, aes(x = .data[[group]], y = .data[[var]])) +
    geom_boxplot(aes(fill = .data[[group]])) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 0.5) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        size = 12,
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 20)
    ) +
    scale_fill_manual(values = cols) +
    ggtitle(var)
}

# calculate statistics for volcano plot
# using linear regression
statVolcano <- function(
  vars,
  reference,
  data,
  fdr_threshold = 0.1
) {
  result <- vector("list")
  for (var in vars) {
    # Create formula with sex and age as covariates
    f_str <- paste0(var, "~", reference, " + sex + age")

    # Use lm for linear regression with covariates
    model <- lm(as.formula(f_str), data = data)

    # Extract p-value for the reference variable (group/diagnosis)
    p_value <-
      broom::tidy(model) |>
      dplyr::slice(2) |>
      dplyr::pull(p.value)

    result[[var]] <- tibble(
      var = var,
      p.value = p_value
    )
  }

  result <- do.call(rbind, result)

  # Calculate adjusted p-values using Benjamini-Hochberg
  result <- result |>
    dplyr::mutate(
      p.adj = p.adjust(p.value, method = "BH"),
      significant = p.adj <= fdr_threshold,
      neg_log10_p = -log10(p.adj)
    )
  # browser()
  return(result)
}

# volcano plot function for cell abundancies
VolPlot <- function(data, cols, n, fdr_threshold = 0.1) {
  # Horizontal line at FDR threshold
  threshold_line <- -log10(fdr_threshold)

  data |>
    ggplot(aes(
      x = log2_ratio,
      y = neg_log10_p,
      color = var,
      label = var
    )) +
    geom_point(size = 3) +
    geom_hline(
      yintercept = threshold_line,
      color = "blue",
      linetype = "dashed"
    ) +
    geom_vline(xintercept = 0, color = "red", linetype = "solid") +
    geom_vline(xintercept = -0.5, color = "red", linetype = "dashed") +
    geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") +
    ggrepel::geom_text_repel() +
    theme_classic() +
    theme(legend.position = "none") +
    xlab(bquote(~ Log[2] ~ "fold change")) +
    ylab(bquote(~ -Log[10] ~ "adjusted p value")) +
    scale_color_manual(values = cols)
}

logfcVolcano <- function(data, group, group1, group2) {
  data <- select(data, all_of(group), all_of(flow_vars))
  data <- group_by(data, .data[[group]])
  data <- summarize(
    data,
    across(all_of(flow_vars), function(x) mean(x, na.rm = TRUE))
  )
  data <- pivot_longer(data, all_of(flow_vars), names_to = "var")
  data <- pivot_wider(data, names_from = all_of(group), values_from = value)
  data <- mutate(data, log2_ratio = log2(.data[[group1]] / .data[[group2]]))
  return(data)
}

# Create boxplots for all flow variables
createBoxplots <- function(data, flow_vars, group, cols) {
  plots <- lapply(
    flow_vars,
    FUN = function(x) {
      ggboxplotFun(
        var = x,
        data = data,
        group = group,
        cols = cols
      )
    }
  )
  return(plots)
}

# Create and save patchwork of boxplots
createAndSaveBoxplots <- function(
  data,
  flow_vars,
  group,
  cols,
  output_file,
  ncol = 4,
  width = 5,
  height = 15
) {
  plots <- createBoxplots(data, flow_vars, group, cols)
  patch <- patchwork::wrap_plots(plots, ncol = ncol)
  ggsave(
    file.path("results", "flow", output_file),
    width = width,
    height = height,
    plot = patch
  )
  return(patch)
}

# Function to create volcano plots for different comparisons
createVolcanoPlot <- function(
  data,
  group_column,
  group1,
  group2,
  tissue,
  output_dir,
  width = 5,
  height = 5,
  top_n = NULL
) {
  # Generate output file name based on parameters
  output_file <- paste0(
    "volcano_",
    group1,
    "_",
    group2,
    "_",
    tissue,
    ".pdf"
  )

  data <- data[data[[group_column]] %in% c(group1, group2), ]

  # Get fold changes first
  fc_data <- logfcVolcano(
    data = data,
    group = group_column,
    group1 = group1,
    group2 = group2
  )

  # Filter to top N variables by absolute fold change if specified
  if (!is.null(top_n)) {
    fc_data <- fc_data |>
      dplyr::arrange(desc(abs(log2_ratio))) |>
      dplyr::slice_head(n = top_n)
  }

  # Get variables to test (after filtering)
  vars_to_test <- fc_data$var

  # Get p-values only for filtered variables
  pval_data <- statVolcano(
    vars_to_test,
    reference = group_column,
    data = data
  )

  # Join data
  vol_data <- left_join(fc_data, pval_data, join_by(var))

  # Create plot
  vol_plot <- VolPlot(
    data = vol_data,
    cols = flow_vars_cols,
    n = length(vars_to_test)
  )

  # Save plot
  ggsave(
    file.path("results", output_dir, output_file),
    plot = vol_plot,
    width = width,
    height = height
  )

  # Return the data and plot for potential further use
  return(list(data = vol_data, plot = vol_plot))
}

# Helper function to create pairwise combinations
createCombinations <- function(
  conditions,
  group_column_name,
  tissues = c("CSF", "blood")
) {
  expand.grid(
    tissue = tissues,
    condition1 = conditions,
    condition2 = conditions
  ) |>
    dplyr::filter(
      condition1 != condition2,
      as.numeric(condition1) < as.numeric(condition2)
    ) |>
    dplyr::mutate(across(everything(), as.character)) |>
    dplyr::mutate(group_column = group_column_name)
}
