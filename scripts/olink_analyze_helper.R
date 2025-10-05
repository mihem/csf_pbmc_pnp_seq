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

calcStats <- function(data_wide, vars) {
    contr <- list()
    for (var in vars) {
        # Check if there are at least four non-NA values for three groups
        var_counts <- data_wide |>
            dplyr::filter(!is.na(.data[[var]])) |>
            dplyr::count(group) |>
            dplyr::filter(n >= 4)
        if (nrow(var_counts) >= 3) {
            formula <- paste0(var, "~ group + sex + age")
            # formula <- paste0(var, "~ group")
            fit <- lm(as.formula(formula), data = data_wide)
            contr[[var]] <- emmeans::emmeans(fit, "group", adjust = "none")
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

createBoxplot <- function(var, data_wide, stats, unit_label) {
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
        ggplot(aes(x = group, y = .data[[var]], fill = group)) +
        geom_boxplot() +
        geom_jitter(width = 0.2) +
        theme_bw() +
        theme(legend.position = "none") +
        scale_fill_manual(values = color_olink_diagnosis) +
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

processOlinkData <- function(data, value_col, unit_label) {
    # Process data with metadata
    data_metadata <- data |>
        left_join(olink_metadata, by = "SampleID") |>
        relocate(Assay, all_of(value_col), group, age, sex)

    # Convert to wide format
    data_wide <- data_metadata |>
        pivot_wider(
            id_cols = c(SampleID, group, age, sex),
            names_from = Assay,
            values_from = all_of(value_col)
        ) |>
        mutate(group = factor(group, levels = c("CTRL", "GBS", "CIDP")))

    # Get assay names and sort alphabetically
    # assays <- sort(unique(data_metadata$Assay))

    # Define assays
    assays <- c("CCL2", "CSF1", "CXCL10", "FASLG", "GZMA", "KLRD1", "TREM2")

    # Calculate statistics
    stats <- calcStats(data_wide, assays)

    # Create boxplots in alphabetical order
    boxplots_list <- lapply(
        assays,
        createBoxplot,
        data_wide = data_wide,
        stats = stats,
        unit_label = unit_label
    )

    boxplots <- patchwork::wrap_plots(boxplots_list)

    return(list(
        data_wide = data_wide,
        stats = stats,
        boxplots = boxplots,
        assays = assays
    ))
}
