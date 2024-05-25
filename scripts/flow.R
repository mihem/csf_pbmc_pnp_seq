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

# load xlsx file ----
flow <-
 read_excel(file.path("raw", "flow", "SEED_flow.xlsx")) |>
 janitor::clean_names() |>
 mutate(birth_date = as_date(birth_date)) |>
 mutate(date = as_date(date))

flow_first <- 
  flow  |>
  group_by(last_name, first_name) |>
  dplyr::filter(date == min(date))

lookup <-
  read_excel(file.path("lookup", "SEED_lookup_v3.xlsx")) |>
  janitor::clean_names() 

flow_first |>
  inner_join(lookup)

lookup |>
  anti_join(flow_first, join_by(last_name, first_name, birth_date, date))


compStat <- function(x_var, group, data, paired) {
# initalize stats
  stats <- vector("list")

  # for character run pairwise fisher test for all parameters, only keep important columns so they match
  for (par in x_var) {
    if (is.character(data[[par]])) {
      stats[[par]] <- rstatix::pairwise_fisher_test(table(data[[group]], data[[par]]), p.adjust.method = "none") |>
        mutate(.y. = par) |>
        select(.y., group1, group2, p, p.adj, p.adj.signif)
    }
    # for numeric run dunn if more than two or wilcox otherwise for all varameters
    if (is.numeric(data[[par]])) {
    f_str <- paste0(par, "~", group)
    if (length(unique(data[[par]])) > 2) {
      if(paired == FALSE) {
        stats[[par]] <- rstatix::dunn_test(as.formula(f_str), data = data, p.adjust.method = "none") |>
          dplyr::select(.y., group1, group2, p, p.adj, p.adj.signif)
      }
      if (paired == TRUE) {
        stats[[par]] <- rstatix::wilcox_test(as.formula(f_str), data = data, p.adjust.method = "none", paired = paired) |>
          dplyr::select(.y., group1, group2, p) |>
          dplyr::mutate(p.adj = NA,
                        p.adj.signif = NA)
      }
    } else {
      stats[[par]] <- rstatix::wilcox_test(as.formula(f_str), data = data, p.adjust.method = "none", paired = paired) |>
        dplyr::select(.y., group1, group2, p) |>
        dplyr::mutate(p.adj = NA,
                      p.adj.signif = NA)
    }
  }
  }
  # combine these lists into a dataframe, do p value adjustment with BH, and then extract only significant values
  stats_df <-
    do.call(rbind, stats) |>
    dplyr::mutate(p.adj = p.adjust(p, method = "BH")) |>
    dplyr::filter(p.adj < 0.05) |>
    dplyr::mutate(p.adj.signif = as.character(symnum(p.adj, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " "))))
  return(stats_df)
}


ggboxplotFun <- function(var, data, group, stats) {
  # extract stats
  stats_df <- dplyr::filter(stats, .y. == var)
  stats_list <- vector("list")
  if(nrow(stats_df) != 0) {
    stats_list$annotation <- stats_df$p.adj.signif
    for (i in 1:nrow(stats_df)) {
      stats_list$comparisons[[i]] <- c(stats_df$group1[i], stats_df$group2[i])
    }
  }

  # plot boxplot
  ggplot(data, aes(x = .data[[group]], y = .data[[var]])) +
    geom_boxplot(aes(fill = .data[[group]])) +
    ggsignif::geom_signif(comparisons = stats_list$comparisons, annotation = stats_list$annotation, textsize = 5, step_increase = 0.05, vjust = 0.7)+
    geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 0.5) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=12),
          #          axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 20)
          )+
    scale_fill_brewer(palette = "Set2")+
    ggtitle(var)
}

BarplotPercent <- function(var, data, group, stats) {
  # extract stats
  stats_df <- dplyr::filter(stats, .y. == var)
  stats_list <- vector("list")
  if(nrow(stats_df) != 0) {
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
    mutate(freq = n/sum(n)) |>
    ungroup() |>
    complete(.data[[group]], .data[[var]], fill = list(n = 0, freq = 0)) |>
    dplyr::filter(.data[[var]] == "yes") |>
    ggplot(aes(x = .data[[group]], y = freq, fill = .data[[group]])) +
    geom_col() + 
     scale_y_continuous(labels = scales::percent_format()) + 
    ggsignif::geom_signif(comparisons = stats_list$comparisons, annotation = stats_list$annotation, textsize = 5, step_increase = 0.05, vjust = 0.7) +
    theme_bw()+
    theme(legend.position = "none",
           axis.title.x = element_blank(),
           axis.text.x = element_text(size=12),
           axis.title.y = element_blank(),
           plot.title = element_text(size = 20)
           )+
    scale_fill_brewer(palette = "Set2")+
    ggtitle(var)
}


# continous variables -----
flow_con <-
    flow |>
    select(cell_count, gluc_ratio, lactate, Bcell:CD56brightNK) |>
    names()

stats_flow_con <- compStat(x_var = flow_con, group = "Group", data = flow, paired = FALSE)

flow_con_plots <- lapply(
    flow_con,
    FUN = function(x) {
        ggboxplotFun(var = x, data = flow, group = "Group", stats = stats_flow_con)
    }
)

flow_con_plots_patch <- patchwork::wrap_plots(flow_con_plots, ncol = 4)

ggsave(file.path("results", "flow", "flow_con_plots.pdf"), width = 10, height = 15, plot = flow_con_plots_patch)

# categorical variables -----

flow_cat <-
    flow |>
    select(BCBD, IgG:IgM) |>
    names()

stats_flow_cat <- compStat(x_var = flow_cat, group = "Group", data = flow, paired = FALSE)

flow_cat_plots <- lapply(
    flow_cat,
    FUN = function(x) {
        BarplotPercent(var = x, data = flow, group = "Group", stats = stats_flow_cat)
    }
)

BarplotPercent(var = "BCBD", data = flow, group = "Group", stats = stats_flow_cat)

flow_cat_plots_patch <- patchwork::wrap_plots(flow_cat_plots, ncol = 4)

ggsave(file.path("results", "flow", "flow_cat_plots.pdf"), width = 10, height = 5, plot = flow_cat_plots_patch)

# abundance volcano plot ----

# calculcate statistics for volcano plot
statVolcano <- function(vars, reference, data) {
  result <- vector("list")
  for (var in vars) {
    f_str <- paste0(var, "~", reference)
    result[[var]] <- broom::tidy(wilcox.test(as.formula(f_str), data = data)) |>
      mutate(var = var) |>
      select(var, p.value)
  }
  result <-
    do.call(rbind, result) |>
    mutate(p.adj = p.adjust(p.value, method = "BH"))
  return(result)
}


# volcano plot function for cell abundancies
VolPlot <- function(data) {
  data |>
    ggplot(aes(x = log2_ratio, y = neg_log10_p, color = var, label = var)) +
    geom_point(size = 3) +
    geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05 / length(vars)), color = "blue", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    # geom_vline(xintercept = -0.5, color = "red", linetype = "dashed") +
    # geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") +
    ggrepel::geom_text_repel() +
    # ggrepel::geom_label_repel() +
    # geom_text(nudge_y = nudge_y) +
    theme_classic() +
    theme(legend.position = "none") +
    xlab(bquote(~ Log[2] ~ "fold change")) +
    ylab(bquote(~ -Log[10] ~ "adjusted p value")) +
    scale_color_manual(values = flow_cols)
}

my_cols_25 <- pals::cols25()

flow_cols <- my_cols_25
names(flow_cols) <- flow_con

# pms vs rrms
flow_pms_rrms <- 
  flow |>
  dplyr::filter(Group %in% c("RRMS", "PMS"))

flow_pms_rrms_pval <- statVolcano(flow_con, reference = "Group", data = flow_pms_rrms)

flow_pms_rrms_fc <-
  flow_pms_rrms |>
  select(Group, all_of(flow_con)) |>
  group_by(Group) |>
  summarize(across(flow_con, function(x) mean(x, na.rm = TRUE))) |>
  pivot_longer(flow_con, names_to = "var") |>
  pivot_wider(names_from = Group, values_from = value) |>
  mutate(log2_ratio = log2(PMS / RRMS))

# join p values and fc
flow_pms_rrms_volcano <-
  flow_pms_rrms_fc |>
  left_join(flow_pms_rrms_pval) |>
  mutate(neg_log10_p = -log10(p.adj)) |>
  mutate(var = if_else(neg_log10_p < -log10(0.05), "", var))
  # mutate(var = if_else(neg_log10_p < -log10(0.05) | abs(log2_ratio) < 0.5, "", var))

VolPlot(data = flow_pms_rrms_volcano)
ggsave(file.path("results", "flow", "volcano_pms_rrms.pdf"), width = 5, height = 5)

