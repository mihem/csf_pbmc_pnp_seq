##################################################
# Helper function to filter red cells
##################################################

filterRedCells <- function(data, excel_file, value_col) {
    # Retrieve formatting information from the Excel file
    cells <- tidyxl::xlsx_cells(excel_file)
    formats <- tidyxl::xlsx_formats(excel_file)

    # Find the column number for the value column
    val_col_num <- cells |>
        filter(row == 1, character == value_col) |>
        pull(col)

    # Filter cells to only include the value column (excluding header)
    value_cells <- cells |>
        filter(col == val_col_num, row > 1)

    # Identify red-marked cells
    red_cells <- value_cells |>
        filter(!is.na(local_format_id)) |>
        mutate(
            is_red = map_lgl(
                local_format_id,
                function(x) {
                    if (
                        x > 0 &&
                            x <=
                                length(
                                    formats$local$fill$patternFill$fgColor$rgb
                                )
                    ) {
                        fg_rgb <- formats$local$fill$patternFill$fgColor$rgb[x]
                        if (!is.na(fg_rgb)) {
                            return(
                                fg_rgb %in%
                                    c("FFF5AA9C", "FFFF0000", "FFFF6B6B")
                            )
                        }
                    }
                    return(FALSE)
                }
            )
        ) |>
        select(row, col, is_red)

    # Join red cell information with the data
    value_data <- value_cells |>
        filter(!is_blank) |>
        left_join(red_cells, by = c("row", "col")) |>
        select(row, content, is_red) |>
        rename(uncertain = is_red) |>
        mutate(
            row_index = row - 1,
            uncertain = ifelse(is.na(uncertain), FALSE, uncertain)
        ) |>
        select(row_index, uncertain)

    # Filter data to remove uncertain values
    data_filtered <- data |>
        mutate(row_index = row_number()) |>
        left_join(value_data, by = "row_index") |>
        select(-row_index) |>
        dplyr::filter(uncertain == FALSE | is.na(uncertain))

    return(data_filtered)
}


##################################################
# Helper functions for analysis
##################################################

calcStats <- function(data_wide, vars, group_var) {
    n_group <- length(unique(data_wide[[group_var]]))
    contr <- list()
    for (var in vars) {
        # Check if there are at least four non-NA values for all groups
        var_counts <- data_wide |>
            dplyr::filter(!is.na(.data[[var]])) |>
            dplyr::count(.data[[group_var]]) |>
            dplyr::filter(n >= 4)
        if (nrow(var_counts) >= n_group) {
            # Use mixed-effects model with orbis_id as random effect to account for repeated measurements
            formula <- paste0(var, " ~ ", group_var, " + sex + age + (1|orbis_id)")
            fit <- lme4::lmer(as.formula(formula), data = data_wide)
            contr[[var]] <- emmeans::emmeans(fit, group_var, adjust = "none")
            contr[[var]] <- pairs(contr[[var]], adjust = "none")
            contr[[var]] <- broom::tidy(contr[[var]])
            contr[[var]] <- tidyr::separate_wider_regex(
                contr[[var]],
                contrast,
                pattern = c(group1 = "\\w+", " - ", group2 = "\\w+")
            )
        } else {
            message("Skipping variable ", var, " due to insufficient data.")
        }
    }

    stats_df <- bind_rows(contr, .id = "var") |>
        dplyr::mutate(p.adj = p.adjust(p.value, method = "BH")) |>
        dplyr::mutate(
            p.adj.signif = as.character(symnum(
                p.adj,
                corr = FALSE,
                na = FALSE,
                cutpoints = c(0, 0.001, 0.01, 0.1, 1),
                symbols = c("***", "**", "*", " ")
            ))
        )
    return(stats_df)
}

createBoxplot <- function(var, data_wide, stats, unit_label, group_var, colors) {
    stats_var <- dplyr::filter(stats, var == !!var, p.adj < 0.1)
    if (nrow(stats_var) != 0) {
        stats_list <- list()
        stats_list$annotation <- stats_var$p.adj.signif
        for (i in 1:nrow(stats_var)) {
            stats_list$comparisons[[i]] <- c(
                stats_var$group1[i],
                stats_var$group2[i]
            )
        }
    }
    boxplot <-
        data_wide |>
        ggplot(aes(x = .data[[group_var]], y = .data[[var]], fill = .data[[group_var]])) +
        geom_boxplot() +
        geom_jitter(width = 0.2) +
        theme_bw() +
        theme(legend.position = "none") +
        scale_fill_manual(values = colors) +
        ggtitle(paste0(var, " (", unit_label, ")")) +
        xlab(NULL) +
        ylab(NULL)

    if (nrow(stats_var) != 0) {
        boxplot <- boxplot +
            ggsignif::geom_signif(
                comparisons = stats_list$comparisons,
                annotation = stats_list$annotation,
                step_increase = 0.05,
                vjust = 0.7
            )
    }
    return(boxplot)
}

processOlinkData <- function(data, value_col, unit_label, group_var, group_levels, colors) {
    # Process data with metadata
    data_metadata <- data |>
        left_join(olink_metadata, by = "SampleID") |>
        relocate(Assay, all_of(value_col), orbis_id, all_of(group_var), age, sex)

    # Convert to wide format
    id_cols <- c("SampleID", "orbis_id", group_var, "age", "sex")
    data_wide <- data_metadata |>
        pivot_wider(
            id_cols = all_of(id_cols),
            names_from = Assay,
            values_from = all_of(value_col)
        ) |>
        mutate(!!group_var := factor(.data[[group_var]], levels = group_levels))

    # Get assay names and sort alphabetically
    assays <- sort(unique(data_metadata$Assay))

    # Calculate statistics
    stats <- calcStats(data_wide, assays, group_var)

    # Create boxplots in alphabetical order
    boxplots_list <- lapply(
        assays,
        createBoxplot,
        data_wide = data_wide,
        stats = stats,
        unit_label = unit_label,
        group_var = group_var,
        colors = colors
    )

    boxplots <- patchwork::wrap_plots(boxplots_list)

    return(list(
        data_wide = data_wide,
        stats = stats,
        boxplots = boxplots,
        assays = assays
    ))
}


##################################################
# Volcano plot functions for Olink data
##################################################

# Calculate statistics for volcano plot using linear mixed-effects model
statVolcanoOlink <- function(
    vars,
    reference,
    data,
    fdr_threshold = 0.1
) {
    result <- vector("list")
    n_group <- length(unique(data[[reference]]))
    
    for (var in vars) {
        # Check if there are at least four non-NA values for all groups
        var_counts <- data |>
            dplyr::filter(!is.na(.data[[var]])) |>
            dplyr::count(.data[[reference]]) |>
            dplyr::filter(n >= 4)
        
        if (nrow(var_counts) >= n_group) {
            # Create formula with sex and age as covariates and orbis_id as random effect
            f_str <- paste0(var, " ~ ", reference, " + sex + age + (1|orbis_id)")

            # Use mixed-effects model with orbis_id as random effect to account for repeated measurements
            model <- lme4::lmer(as.formula(f_str), data = data)

            # Extract p-value for the reference variable (group/diagnosis)
            emm <- emmeans::emmeans(model, reference, adjust = "none")
            contr <- pairs(emm, adjust = "none")
            contr_tidy <- broom::tidy(contr)
            
            # Get p-value from the first contrast (assuming two groups)
            p_value <- contr_tidy$p.value[1]

            result[[var]] <- tibble(
                var = var,
                p.value = p_value
            )
        } else {
            message("Skipping variable ", var, " due to insufficient data.")
        }
    }

    # Check if we have any results
    if (length(result) == 0) {
        warning("No variables with sufficient data for analysis")
        return(NULL)
    }

    result <- do.call(rbind, result)

    # Apply BH correction (same as in calcStats for boxplots)
    result <- result |>
        dplyr::mutate(
            p.adj = p.adjust(p.value, method = "BH"),
            significant = p.adj < fdr_threshold,
            p.adj.threshold = fdr_threshold
        ) |>
        mutate(neg_log10_p = -log10(p.value))

    return(result)
}

# Calculate log2 fold change for volcano plot
logfcVolcanoOlink <- function(data, group, group1, group2, vars) {
    data <- select(data, .data[[group]], all_of(vars))
    data <- group_by(data, .data[[group]])
    data <- summarize(
        data,
        across(all_of(vars), function(x) mean(x, na.rm = TRUE))
    )
    data <- pivot_longer(data, all_of(vars), names_to = "var")
    data <- pivot_wider(data, names_from = group, values_from = value)
    data <- mutate(data, log2_ratio = log2(.data[[group1]] / .data[[group2]]))
    return(data)
}

# Volcano plot visualization function
VolPlotOlink <- function(data, colors, n) {
    # Get the corrected threshold for the horizontal line
    threshold_line <- -log10(data$p.adj.threshold[1])

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
        geom_vline(xintercept = -1, color = "red", linetype = "dashed") +
        geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
        ggrepel::geom_text_repel() +
        theme_classic() +
        theme(legend.position = "none") +
        xlab(bquote(~ Log[2] ~ "fold change")) +
        ylab(bquote(~ -Log[10] ~ "p value")) +
        scale_color_manual(values = colors)
}

# Function to create volcano plots for different comparisons
createVolcanoPlotOlink <- function(
    data_wide,
    assays,
    group_column,
    group1,
    group2,
    output_dir,
    suffix = "",
    width = 5,
    height = 5,
    top_n = NULL,
    colors = NULL
) {
    # Generate output file name based on parameters
    output_file <- paste0(
        "volcano_",
        group1,
        "_vs_",
        group2,
        suffix,
        ".pdf"
    )

    # Filter data to only include the two groups being compared
    data <- data_wide[data_wide[[group_column]] %in% c(group1, group2), ]

    # Get fold changes first
    fc_data <- logfcVolcanoOlink(
        data = data,
        group = group_column,
        group1 = group1,
        group2 = group2,
        vars = assays
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
    pval_data <- statVolcanoOlink(
        vars_to_test,
        reference = group_column,
        data = data
    )

    # Join data
    vol_data <- left_join(fc_data, pval_data, join_by(var))

    # Set up colors for proteins if not provided
    if (is.null(colors)) {
        colors <- setNames(
            rep("black", length(vars_to_test)),
            vars_to_test
        )
    }

    # Create plot
    vol_plot <- VolPlotOlink(
        data = vol_data,
        colors = colors,
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
