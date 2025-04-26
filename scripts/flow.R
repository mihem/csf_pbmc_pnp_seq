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
flow_pre <-
  read_excel(file.path("raw", "flow", "SEED_flow_v6.xlsx")) |>
  mutate(date = as_date(date))

lookup <-
  read_excel(file.path("lookup", "SEED_lookup_v6.xlsx")) |>
  janitor::clean_names() |>
  mutate(diagnosis = factor(diagnosis, levels = diagnosis_order)) |>
  mutate(group = factor(group, levels = group_order))

# sanity check
lookup |>
  anti_join(flow_pre, join_by(last_name, first_name, date)) |>
  select(last_name, first_name, birth_date, date)
# 5 missing because either no flow cytometry at all or only at another date

# join flow and lookup
flow <-
  flow_pre |>
  inner_join(lookup, join_by(last_name, first_name, date)) |>
  (function(df) split(df, df$tissue))()

#function to show significant comparisons using dunn test
compStat <- function(x_var, group, data) {
  # initalize stats
  stats <- vector("list")

  # for character run pairwise fisher test for all parameters, only keep important columns so they match
  for (par in x_var) {
    if (is.character(data[[par]])) {
      stats[[par]] <- rstatix::pairwise_fisher_test(
        table(data[[group]], data[[par]]),
        p.adjust.method = "BH"
      ) |>
        mutate(.y. = par) |>
        select(.y., group1, group2, p, p.adj, p.adj.signif)
    }
    # for numeric run dunn if more than two or wilcox otherwise for all parameters
    if (is.numeric(data[[par]])) {
      f_str <- paste0(par, "~", group)
      if (length(unique(data[[group]])) > 2) {
        stats[[par]] <- rstatix::dunn_test(
          as.formula(f_str),
          data = data,
          p.adjust.method = "BH"
        ) |>
          select(.y., group1, group2, p, p.adj, p.adj.signif)
      } else {
        stats[[par]] <- rstatix::wilcox_test(as.formula(f_str), data = data) |>
          select(.y., group1, group2, p) |>
          mutate(p.adj = p)
      }
    }
  }
  # combine these lists into a dataframe and extract significant comparions
  stats_df <- do.call(rbind, stats) |>
    # dplyr::mutate(p.adj = p.adjust(p.adj, method = "BH" )) |>
    dplyr::filter(p.adj < 0.05) |>
    mutate(
      p.adj.signif = as.character(symnum(
        p.adj,
        corr = FALSE,
        na = FALSE,
        cutpoints = c(0, 0.001, 0.01, 0.05, 1),
        symbols = c("***", "**", "*", " ")
      ))
    )
  return(stats_df)
}


ggboxplotFun <- function(var, data, group, stats, cols) {
  # extract stats
  stats_df <- dplyr::filter(stats, .y. == var)
  stats_list <- vector("list")
  if (nrow(stats_df) != 0) {
    stats_list$annotation <- stats_df$p.adj.signif
    for (i in 1:nrow(stats_df)) {
      stats_list$comparisons[[i]] <- c(stats_df$group1[i], stats_df$group2[i])
    }
  }

  # plot boxplot
  ggplot(data, aes(x = .data[[group]], y = .data[[var]])) +
    geom_boxplot(aes(fill = .data[[group]])) +
    ggsignif::geom_signif(
      comparisons = stats_list$comparisons,
      annotation = stats_list$annotation,
      textsize = 5,
      step_increase = 0.05,
      vjust = 0.7
    ) +
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

BarplotPercent <- function(var, data, group, stats) {
  # extract stats
  stats_df <- dplyr::filter(stats, .y. == var)
  stats_list <- vector("list")
  if (nrow(stats_df) != 0) {
    stats_list$annotation <- stats_df$p.adj.signif
    for (i in 1:nrow(stats_df)) {
      stats_list$comparisons[[i]] <- c(stats_df$group1[i], stats_df$group2[i])
    }
  }
  # plot
  data |>
    drop_na(.data[[var]]) |>
    group_by(.data[[group]]) |>
    count(.data[[var]]) |>
    mutate(freq = n / sum(n)) |>
    ungroup() |>
    complete(.data[[group]], .data[[var]], fill = list(n = 0, freq = 0)) |>
    dplyr::filter(.data[[var]] == "yes") |>
    ggplot(aes(x = .data[[group]], y = freq, fill = .data[[group]])) +
    geom_col() +
    scale_y_continuous(labels = scales::percent_format()) +
    ggsignif::geom_signif(
      comparisons = stats_list$comparisons,
      annotation = stats_list$annotation,
      textsize = 5,
      step_increase = 0.05,
      vjust = 0.7
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 20)
    ) +
    # scale_fill_brewer(palette = "Set2")+
    scale_fill_manual(values = colorset_dutch) +
    ggtitle(var)
}


# CSF flow boxplots -----
flow_vars <-
  flow$CSF |>
  select(Gran:intMono) |>
  names()

# PNP vs CTRL
stats_csf_con_group <- compStat(
  x_var = flow_vars,
  group = "group",
  data = flow$CSF
)

csf_con_plots_group <- lapply(
  flow_vars,
  FUN = function(x) {
    ggboxplotFun(
      var = x,
      data = flow$CSF,
      group = "group",
      stats = stats_csf_con_group,
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
stats_csf_con_dx <- compStat(
  x_var = flow_vars,
  group = "diagnosis",
  data = flow$CSF
)

csf_con_plots_dx <- lapply(
  flow_vars,
  FUN = function(x) {
    ggboxplotFun(
      var = x,
      data = flow$CSF,
      group = "diagnosis",
      stats = stats_csf_con_dx,
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

# PNP vs CTRL
stats_blood_con_group <- compStat(
  x_var = flow_vars,
  group = "group",
  data = flow$blood
)

blood_con_plots_group <- lapply(
  flow_vars,
  FUN = function(x) {
    ggboxplotFun(
      var = x,
      data = flow$blood,
      group = "group",
      stats = stats_blood_con_group,
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
stats_blood_con_dx <- compStat(
  x_var = flow_vars,
  group = "diagnosis",
  data = flow$blood
)

blood_con_plots_dx <- lapply(
  flow_vars,
  FUN = function(x) {
    ggboxplotFun(
      var = x,
      data = flow$blood,
      group = "diagnosis",
      stats = stats_blood_con_dx,
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

# calculcate statistics for volcano plot
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

# volcano plot function for cell abundancies
VolPlot <- function(data, cols, n) {
  data |>
    ggplot(aes(x = log2_ratio, y = neg_log10_p, color = var, label = var)) +
    geom_point(size = 3) +
    # geom_hline(yintercept = -log10(0.05 / n), color = "blue", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    # geom_vline(xintercept = -0.5, color = "red", linetype = "dashed") +
    # geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") +
    # ggrepel::geom_text_repel() +
    ggrepel::geom_label_repel() +
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

# volcano plots PNP vs CTRL csf ----
flow_pnp_ctrl_csf_pval <- statVolcano(
  flow_vars,
  reference = "group",
  data = flow$CSF
)

flow_pnp_ctrl_csf_fc <- logfcVolcano(
  data = flow$CSF,
  group = "group",
  group1 = "PNP",
  group2 = "CTRL"
)

flow_pnp_ctrl_csf_vol_data <- left_join(
  flow_pnp_ctrl_csf_fc,
  flow_pnp_ctrl_csf_pval,
  join_by(var)
)

flow_pnp_ctrl_csf_vol_plot <- VolPlot(
  data = flow_pnp_ctrl_csf_vol_data,
  cols = flow_vars_cols,
  n = length(flow_vars)
)

ggsave(
  file.path("results", "flow", "volcano_pnp_ctrl_csf.pdf"),
  plot = flow_pnp_ctrl_csf_vol_plot,
  width = 5,
  height = 5
)

# volcano plots PNP vs CTRL blood ----
flow_pnp_ctrl_blood_pval <- statVolcano(
  flow_vars,
  reference = "group",
  data = flow$blood
)

flow_pnp_ctrl_blood_fc <- logfcVolcano(
  data = flow$blood,
  group = "group",
  group1 = "PNP",
  group2 = "CTRL"
)

flow_pnp_ctrl_blood_vol_data <- left_join(
  flow_pnp_ctrl_blood_fc,
  flow_pnp_ctrl_blood_pval,
  join_by(var)
)

flow_pnp_ctrl_blood_vol_plot <- VolPlot(
  data = flow_pnp_ctrl_blood_vol_data,
  cols = flow_vars_cols,
  n = length(flow_vars)
)

ggsave(
  file.path("results", "flow", "volcano_pnp_ctrl_blood.pdf"),
  plot = flow_pnp_ctrl_blood_vol_plot,
  width = 5,
  height = 5
)

# CSF volcano plots subgroups ----
# ciap vs ctrl csf
flow_ciap_ctrl <- filter(flow$CSF, diagnosis %in% c("CIAP", "CTRL"))
flow_ciap_ctrl_csf_pval <- statVolcano(
  flow_vars,
  reference = "diagnosis",
  data = flow_ciap_ctrl
)

flow_ciap_ctrl_csf_fc <- logfcVolcano(
  data = flow_ciap_ctrl,
  group = "diagnosis",
  group1 = "CIAP",
  group2 = "CTRL"
)

flow_ciap_ctrl_csf_vol_data <- left_join(
  flow_ciap_ctrl_csf_fc,
  flow_ciap_ctrl_csf_pval,
  join_by(var)
)

flow_ciap_ctrl_csf_vol_plot <- VolPlot(
  data = flow_ciap_ctrl_csf_vol_data,
  cols = flow_vars_cols,
  n = length(flow_vars)
)

ggsave(
  file.path("results", "flow", "volcano_ciap_ctrl_csf.pdf"),
  plot = flow_ciap_ctrl_csf_vol_plot,
  width = 5,
  height = 5
)

# cidp vs ctrl csf
flow_cidp_ctrl <- filter(flow$CSF, diagnosis %in% c("CIDP", "CTRL"))

flow_cidp_ctrl_csf_pval <- statVolcano(
  flow_vars,
  reference = "diagnosis",
  data = flow_cidp_ctrl
)

flow_cidp_ctrl_csf_fc <- logfcVolcano(
  data = flow_cidp_ctrl,
  group = "diagnosis",
  group1 = "CIDP",
  group2 = "CTRL"
)

flow_cidp_ctrl_csf_vol_data <- left_join(
  flow_cidp_ctrl_csf_fc,
  flow_cidp_ctrl_csf_pval,
  join_by(var)
)

flow_cidp_ctrl_csf_vol_plot <- VolPlot(
  data = flow_cidp_ctrl_csf_vol_data,
  cols = flow_vars_cols,
  n = length(flow_vars)
)

ggsave(
  file.path("results", "flow", "volcano_cidp_ctrl_csf.pdf"),
  plot = flow_cidp_ctrl_csf_vol_plot,
  width = 5,
  height = 5
)

# gbs vs ctrl csf
flow_gbs_ctrl <- filter(flow$CSF, diagnosis %in% c("GBS", "CTRL"))

flow_gbs_ctrl_csf_pval <- statVolcano(
  flow_vars,
  reference = "diagnosis",
  data = flow_gbs_ctrl
)

flow_gbs_ctrl_csf_fc <- logfcVolcano(
  data = flow_gbs_ctrl,
  group = "diagnosis",
  group1 = "GBS",
  group2 = "CTRL"
)

flow_gbs_ctrl_csf_vol_data <- left_join(
  flow_gbs_ctrl_csf_fc,
  flow_gbs_ctrl_csf_pval,
  join_by(var)
)

flow_gbs_ctrl_csf_vol_plot <- VolPlot(
  data = flow_gbs_ctrl_csf_vol_data,
  cols = flow_vars_cols,
  n = length(flow_vars)
)

ggsave(
  file.path("results", "flow", "volcano_gbs_ctrl_csf.pdf"),
  plot = flow_gbs_ctrl_csf_vol_plot,
  width = 5,
  height = 5
)


# cidp vs ciap csf
flow_cidp_ciap <- filter(flow$CSF, diagnosis %in% c("CIDP", "CIAP"))
flow_cidp_ciap_csf_pval <- statVolcano(
  flow_vars,
  reference = "diagnosis",
  data = flow_cidp_ciap
)

flow_cidp_ciap_csf_fc <- logfcVolcano(
  data = flow_cidp_ciap,
  group = "diagnosis",
  group1 = "CIDP",
  group2 = "CIAP"
)

flow_cidp_ciap_csf_vol_data <- left_join(
  flow_cidp_ciap_csf_fc,
  flow_cidp_ciap_csf_pval,
  join_by(var)
)

flow_cidp_ciap_csf_vol_plot <- VolPlot(
  data = flow_cidp_ciap_csf_vol_data,
  cols = flow_vars_cols,
  n = length(flow_vars)
)

ggsave(
  file.path("results", "flow", "volcano_cidp_ciap_csf.pdf"),
  plot = flow_cidp_ciap_csf_vol_plot,
  width = 5,
  height = 5
)

# gbs vs ciap csf
flow_gbs_ciap <- filter(flow$CSF, diagnosis %in% c("GBS", "CIAP"))
flow_gbs_ciap_csf_pval <- statVolcano(
  flow_vars,
  reference = "diagnosis",
  data = flow_gbs_ciap
)

flow_gbs_ciap_csf_fc <- logfcVolcano(
  data = flow_gbs_ciap,
  group = "diagnosis",
  group1 = "GBS",
  group2 = "CIAP"
)

flow_gbs_ciap_csf_vol_data <- left_join(
  flow_gbs_ciap_csf_fc,
  flow_gbs_ciap_csf_pval,
  join_by(var)
)

flow_gbs_ciap_csf_vol_plot <- VolPlot(
  data = flow_gbs_ciap_csf_vol_data,
  cols = flow_vars_cols,
  n = length(flow_vars)
)

ggsave(
  file.path("results", "flow", "volcano_gbs_ciap_csf.pdf"),
  plot = flow_gbs_ciap_csf_vol_plot,
  width = 5,
  height = 5
)

# blood volcano plots subgroups ----
# ciap vs ctrl blood
flow_ciap_ctrl <- filter(flow$blood, diagnosis %in% c("CIAP", "CTRL"))

flow_ciap_ctrl_blood_pval <- statVolcano(
  flow_vars,
  reference = "diagnosis",
  data = flow_ciap_ctrl
)

flow_ciap_ctrl_blood_fc <- logfcVolcano(
  data = flow_ciap_ctrl,
  group = "diagnosis",
  group1 = "CIAP",
  group2 = "CTRL"
)

flow_ciap_ctrl_blood_vol_data <- left_join(
  flow_ciap_ctrl_blood_fc,
  flow_ciap_ctrl_blood_pval,
  join_by(var)
)

flow_ciap_ctrl_blood_vol_plot <- VolPlot(
  data = flow_ciap_ctrl_blood_vol_data,
  cols = flow_vars_cols,
  n = length(flow_vars)
)

ggsave(
  file.path("results", "flow", "volcano_ciap_ctrl_blood.pdf"),
  plot = flow_ciap_ctrl_blood_vol_plot,
  width = 5,
  height = 5
)

# cidp vs ctrl blood
flow_cidp_ctrl <- filter(flow$blood, diagnosis %in% c("CIDP", "CTRL"))

flow_cidp_ctrl_blood_pval <- statVolcano(
  flow_vars,
  reference = "diagnosis",
  data = flow_cidp_ctrl
)

flow_cidp_ctrl_blood_fc <- logfcVolcano(
  data = flow_cidp_ctrl,
  group = "diagnosis",
  group1 = "CIDP",
  group2 = "CTRL"
)

flow_cidp_ctrl_blood_vol_data <- left_join(
  flow_cidp_ctrl_blood_fc,
  flow_cidp_ctrl_blood_pval,
  join_by(var)
)

flow_cidp_ctrl_blood_vol_plot <- VolPlot(
  data = flow_cidp_ctrl_blood_vol_data,
  cols = flow_vars_cols,
  n = length(flow_vars)
)

ggsave(
  file.path("results", "flow", "volcano_cidp_ctrl_blood.pdf"),
  plot = flow_cidp_ctrl_blood_vol_plot,
  width = 5,
  height = 5
)

# gbs vs ctrl blood
flow_gbs_ctrl <- filter(flow$blood, diagnosis %in% c("GBS", "CTRL"))

flow_gbs_ctrl_blood_pval <- statVolcano(
  flow_vars,
  reference = "diagnosis",
  data = flow_gbs_ctrl
)

flow_gbs_ctrl_blood_fc <- logfcVolcano(
  data = flow_gbs_ctrl,
  group = "diagnosis",
  group1 = "GBS",
  group2 = "CTRL"
)

flow_gbs_ctrl_blood_vol_data <- left_join(
  flow_gbs_ctrl_blood_fc,
  flow_gbs_ctrl_blood_pval,
  join_by(var)
)

flow_gbs_ctrl_blood_vol_plot <- VolPlot(
  data = flow_gbs_ctrl_blood_vol_data,
  cols = flow_vars_cols,
  n = length(flow_vars)
)

ggsave(
  file.path("results", "flow", "volcano_gbs_ctrl_blood.pdf"),
  plot = flow_gbs_ctrl_blood_vol_plot,
  width = 5,
  height = 5
)
