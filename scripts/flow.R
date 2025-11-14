#####################################
# flow analysis
#####################################

# load libraries ------
library(tidyverse)
library(readxl)
library(patchwork)
library(broom)
library(pals)
library(janitor)

# define order and colors
diagnosis_order <- c(
  "CTRL",
  "CIAP",
  "CIDP",
  "GBS",
  "MAG",
  "MFS",
  "PNC",
  "CAN",
  "PPN"
)

diagnosis_color <- setNames(
  pals::cols25(length(diagnosis_order)),
  diagnosis_order
)

group_order <- c("CTRL", "PNP")
group_color <- setNames(pals::cols25(length(group_order)), group_order)

# load flow and lookup file ----
# load flow data, keep only first measurement if multiple are available
flow_pre <-
  read_excel(file.path("raw", "flow", "flowbasic_v3.xlsx")) |>
  mutate(date = as_date(date)) |>
  group_by(last_name, first_name, tissue) |>
  filter(date == min(date)) |>
  ungroup()

lookup <-
  read_excel(file.path("lookup", "SEED_lookup_v8.xlsx")) |>
  janitor::clean_names() |>
  mutate(age = lubridate::time_length(difftime(date, birth_date), "years")) |>
  mutate(diagnosis = factor(diagnosis, levels = diagnosis_order)) |>
  mutate(group = factor(group, levels = group_order))

# sanity checks
lookup |>
  anti_join(flow_pre, join_by(last_name, first_name)) |>
  select(last_name, first_name, birth_date, date)

flow_pre |>
  anti_join(lookup, join_by(last_name, first_name, date)) |>
  print(n = Inf)

# join flow and lookup
flow <-
  flow_pre |>
  inner_join(lookup, join_by(last_name, first_name)) |>
  (function(df) split(df, df$tissue))()

# #function to show significant comparisons using dunn test
# compStat <- function(x_var, group, data) {
#   # initalize stats
#   stats <- vector("list")

#   # for character run pairwise fisher test for all parameters, only keep important columns so they match
#   for (par in x_var) {
#     if (is.character(data[[par]])) {
#       stats[[par]] <- rstatix::pairwise_fisher_test(
#         table(data[[group]], data[[par]]),
#         p.adjust.method = "BH"
#       ) |>
#         mutate(.y. = par) |>
#         select(.y., group1, group2, p, p.adj, p.adj.signif)
#     }
#     # for numeric run dunn if more than two or wilcox otherwise for all parameters
#     if (is.numeric(data[[par]])) {
#       f_str <- paste0(par, "~", group)
#       if (length(unique(data[[group]])) > 2) {
#         stats[[par]] <- rstatix::dunn_test(
#           as.formula(f_str),
#           data = data,
#           p.adjust.method = "BH"
#         ) |>
#           select(.y., group1, group2, p, p.adj, p.adj.signif)
#       } else {
#         stats[[par]] <- rstatix::wilcox_test(as.formula(f_str), data = data) |>
#           select(.y., group1, group2, p) |>
#           mutate(p.adj = p)
#       }
#     }
#   }
#   # combine these lists into a dataframe and extract significant comparions
#   stats_df <- do.call(rbind, stats) |>
#     # dplyr::mutate(p.adj = p.adjust(p.adj, method = "BH" )) |>
#     dplyr::filter(p.adj < 0.05) |>
#     mutate(
#       p.adj.signif = as.character(symnum(
#         p.adj,
#         corr = FALSE,
#         na = FALSE,
#         cutpoints = c(0, 0.001, 0.01, 0.05, 1),
#         symbols = c("***", "**", "*", " ")
#       ))
#     )
#   return(stats_df)
# }

ggboxplotFun <- function(var, data, group, stats, cols) {
  # # extract stats
  # stats_df <- dplyr::filter(stats, .y. == var)
  # stats_list <- vector("list")
  # if (nrow(stats_df) != 0) {
  #   stats_list$annotation <- stats_df$p.adj.signif
  #   for (i in 1:nrow(stats_df)) {
  #     stats_list$comparisons[[i]] <- c(stats_df$group1[i], stats_df$group2[i])
  #   }
  # }

  # plot boxplot
  ggplot(data, aes(x = .data[[group]], y = .data[[var]])) +
    geom_boxplot(aes(fill = .data[[group]])) +
    # ggsignif::geom_signif(
    #   comparisons = stats_list$comparisons,
    #   annotation = stats_list$annotation,
    #   textsize = 5,
    #   step_increase = 0.05,
    #   vjust = 0.7
    # ) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 0.5) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 20)
    ) +
    scale_fill_manual(values = cols) +
    ggtitle(var)
}

# BarplotPercent <- function(var, data, group, stats) {
#   # extract stats
#   stats_df <- dplyr::filter(stats, .y. == var)
#   stats_list <- vector("list")
#   if (nrow(stats_df) != 0) {
#     stats_list$annotation <- stats_df$p.adj.signif
#     for (i in 1:nrow(stats_df)) {
#       stats_list$comparisons[[i]] <- c(stats_df$group1[i], stats_df$group2[i])
#     }
#   }
#   # plot
#   data |>
#     drop_na(.data[[var]]) |>
#     group_by(.data[[group]]) |>
#     count(.data[[var]]) |>
#     mutate(freq = n / sum(n)) |>
#     ungroup() |>
#     complete(.data[[group]], .data[[var]], fill = list(n = 0, freq = 0)) |>
#     dplyr::filter(.data[[var]] == "yes") |>
#     ggplot(aes(x = .data[[group]], y = freq, fill = .data[[group]])) +
#     geom_col() +
#     scale_y_continuous(labels = scales::percent_format()) +
#     ggsignif::geom_signif(
#       comparisons = stats_list$comparisons,
#       annotation = stats_list$annotation,
#       textsize = 5,
#       step_increase = 0.05,
#       vjust = 0.7
#     ) +
#     theme_bw() +
#     theme(
#       legend.position = "none",
#       axis.title.x = element_blank(),
#       axis.text.x = element_text(size = 12),
#       axis.title.y = element_blank(),
#       plot.title = element_text(size = 20)
#     ) +
#     # scale_fill_brewer(palette = "Set2")+
#     scale_fill_manual(values = colorset_dutch) +
#     ggtitle(var)
# }

# CSF flow boxplots -----
flow_vars <-
  flow$CSF |>
  select(Gran:intMono) |>
  names()

# PNP vs CTRL
# stats_csf_con_group <- compStat(
#   x_var = flow_vars,
#   group = "group",
#   data = flow$CSF
# )

csf_con_plots_group <- lapply(
  flow_vars,
  FUN = function(x) {
    ggboxplotFun(
      var = x,
      data = flow$CSF,
      group = "group",
      # stats = stats_csf_con_group,
      cols = group_color
    )
  }
)

csf_con_plots_patch_group <- patchwork::wrap_plots(
  csf_con_plots_group,
  ncol = 4
)

ggsave(
  file.path("results", "flow", "csf_con_plots_group.pdf"),
  width = 5,
  height = 15,
  plot = csf_con_plots_patch_group
)

# subgroups
# stats_csf_con_dx <- compStat(
#   x_var = flow_vars,
#   group = "diagnosis",
#   data = flow$CSF
# )

csf_con_plots_dx <- lapply(
  flow_vars,
  FUN = function(x) {
    ggboxplotFun(
      var = x,
      data = flow$CSF,
      group = "diagnosis",
      # stats = stats_csf_con_dx,
      cols = diagnosis_color
    )
  }
)

csf_con_plots_patch_dx <- patchwork::wrap_plots(csf_con_plots_dx, ncol = 4)

ggsave(
  file.path("results", "flow", "csf_con_plots_dx.pdf"),
  width = 10,
  height = 15,
  plot = csf_con_plots_patch_dx
)


# blood flow boxplots -----

# # PNP vs CTRL
# stats_blood_con_group <- compStat(
#   x_var = flow_vars,
#   group = "group",
#   data = flow$blood
# )

blood_con_plots_group <- lapply(
  flow_vars,
  FUN = function(x) {
    ggboxplotFun(
      var = x,
      data = flow$blood,
      group = "group",
      # stats = stats_blood_con_group,
      cols = group_color
    )
  }
)

blood_con_plots_patch_group <- patchwork::wrap_plots(
  blood_con_plots_group,
  ncol = 4
)

ggsave(
  file.path("results", "flow", "blood_con_plots_group.pdf"),
  width = 5,
  height = 15,
  plot = blood_con_plots_patch_group
)

# subgroups
# stats_blood_con_dx <- compStat(
#   x_var = flow_vars,
#   group = "diagnosis",
#   data = flow$blood
# )

blood_con_plots_dx <- lapply(
  flow_vars,
  FUN = function(x) {
    ggboxplotFun(
      var = x,
      data = flow$blood,
      group = "diagnosis",
      # stats = stats_blood_con_dx,
      cols = diagnosis_color
    )
  }
)

blood_con_plots_patch_dx <- patchwork::wrap_plots(blood_con_plots_dx, ncol = 4)

ggsave(
  file.path("results", "flow", "blood_con_plots_dx.pdf"),
  width = 10,
  height = 15,
  plot = blood_con_plots_patch_dx
)

# continous routine variables -----
# routine_con <-
#     flow$CSF |>
#     select(cell_count:csf_protein, glucose_csf:lactate_csf) |>
#     names()

# stats_routine_con <- compStat(x_var = routine_con, group = "group", data = flow$CSF)

# flow_cat_plots <- lapply(
#     flow_cat,
#     FUN = function(x) {
#         BarplotPercent(var = x, data = flow, group = "Group", stats = stats_flow_cat)
#     }
# )

# BarplotPercent(var = "BCBD", data = flow, group = "Group", stats = stats_flow_cat)

# flow_cat_plots_patch <- patchwork::wrap_plots(flow_cat_plots, ncol = 4)

# ggsave(file.path("results", "flow", "flow_cat_plots.pdf"), width = 10, height = 5, plot = flow_cat_plots_patch)

# abundance volcano plot ----

# calculate statistics for volcano plot
# Original version using Wilcoxon test
statVolcano <- function(vars, reference, data) {
  result <- vector("list")
  for (var in vars) {
    f_str <- paste0(var, "~", reference)
    result[[var]] <- broom::tidy(wilcox.test(as.formula(f_str), data = data)) |>
      dplyr::mutate(var = var) |>
      dplyr::select(var, p.value)
  }
  result <-
    do.call(rbind, result) |>
    dplyr::mutate(p.adj = p.adjust(p.value, method = "BH")) |>
    mutate(neg_log10_p = -log10(p.value)) |>
    mutate(neg_log10_p_adj = -log10(p.adj))
  return(result)
}

# calculate statistics for volcano plot
# using linear regression
statVolcano <- function(vars, reference, data) {
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

  result <-
    do.call(rbind, result) |>
    dplyr::mutate(p.adj = p.adjust(p.value, method = "BH")) |>
    mutate(neg_log10_p = -log10(p.value)) |>
    mutate(neg_log10_p_adj = -log10(p.adj))

  return(result)
}


old <- statVolcano(
  flow_vars,
  reference = "group",
  data = flow$CSF
)

new <- statVolcano(
  flow_vars,
  reference = "group",
  data = flow$CSF
)

flow$CSF$age


# volcano plot function for cell abundancies
VolPlot <- function(data, cols, n) {
  data |>
    ggplot(aes(x = log2_ratio, y = neg_log10_p_adj, color = var, label = var)) +
    geom_point(size = 3) +
    geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed") +
    # geom_hline(
    #   yintercept = -log10(0.05 / n),
    #   color = "blue",
    #   linetype = "solid"
    # ) +
    geom_vline(xintercept = 0, color = "red", linetyp = "solid") +
    geom_vline(xintercept = -1, color = "red", linetype = "dashed") +
    geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
    ggrepel::geom_text_repel() +
    # ggrepel::geom_label_repel() +
    # geom_text(nudge_y = nudge_y) +
    theme_classic() +
    theme(legend.position = "none") +
    xlab(bquote(~ Log[2] ~ "fold change")) +
    ylab(bquote(~ -Log[10] ~ "p value")) +
    scale_color_manual(values = cols)
}

logfcVolcano <- function(data, group, group1, group2) {
  data <- select(data, .data[[group]], all_of(flow_vars))
  data <- group_by(data, .data[[group]])
  data <- summarize(data, across(flow_vars, function(x) mean(x, na.rm = TRUE)))
  data <- pivot_longer(data, flow_vars, names_to = "var")
  data <- pivot_wider(data, names_from = group, values_from = value)
  data <- mutate(data, log2_ratio = log2(.data[[group1]] / .data[[group2]]))
  return(data)
}

# set colors for flow variables
flow_vars_cols <- setNames(pals::cols25(length(flow_vars)), flow_vars)

# Function to create volcano plots for different comparisons
createVolcanoPlot <- function(
  data,
  group_column,
  group1,
  group2,
  tissue
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

  # Filter data if needed (when comparing specific diagnoses)
  if (group_column == "diagnosis" && group1 != "PNP" && group2 != "PNP") {
    data <- filter(data, .data[[group_column]] %in% c(group1, group2))
  }

  # Get p-values
  pval_data <- statVolcano(
    flow_vars,
    reference = group_column,
    data = data
  )

  # Get fold changes
  fc_data <- logfcVolcano(
    data = data,
    group = group_column,
    group1 = group1,
    group2 = group2
  )

  # Join data
  vol_data <- left_join(fc_data, pval_data, join_by(var))

  # Create plot
  vol_plot <- VolPlot(
    data = vol_data,
    cols = flow_vars_cols,
    n = length(flow_vars)
  )

  # Save plot
  ggsave(
    file.path("results", "flow", output_file),
    plot = vol_plot,
    width = 5,
    height = 5
  )

  # Return the data and plot for potential further use
  return(list(data = vol_data, plot = vol_plot))
}

diagnosis <- factor(
  c("CIDP", "GBS", "CIAP", "CTRL"),
  levels = c("CIDP", "GBS", "CIAP", "CTRL")
)

combinations <- crossing(
  tissue = c("CSF", "blood"),
  condition1 = diagnosis,
  condition2 = diagnosis
) |>
  dplyr::filter(
    condition1 != condition2,
    as.numeric(condition1) < as.numeric(condition2)
  ) |>
  dplyr::mutate(across(everything(), as.character)) |>
  dplyr::mutate(group_column = "diagnosis")


# Create a configuration table for all volcano plot comparisons
volcano_configs <- tibble::tibble(
  group_column = c(rep("group", 2), combinations$group_column),
  group1 = c(rep("PNP", 2), combinations$condition1),
  group2 = c(rep("CTRL", 2), combinations$condition2),
  tissue = c("CSF", "blood", combinations$tissue)
)

# Generate all volcano plots
volcano_results_list <- pmap(
  volcano_configs,
  function(group_column, group1, group2, tissue) {
    createVolcanoPlot(
      data = flow[[tissue]],
      group_column = group_column,
      group1 = group1,
      group2 = group2,
      tissue = tissue
    )
  }
)
