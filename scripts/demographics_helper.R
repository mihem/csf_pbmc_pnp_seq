
########################################################
# Helper function for demographics.R
########################################################

# Helper function for creating boxplots with statistics ----
create_boxplot <- function(data, x_var, y_var, group_var, title, 
                          color_palette, geom_type = "point", 
                          ylim_range = NULL, filter_na = TRUE) {
    
    # Filter NA values if requested
    plot_data <- if (filter_na && y_var != "sex") {
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
        ggplot(aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]])) +
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
    
    # Add y-axis limits if specified
    if (!is.null(ylim_range)) {
        p <- p + ylim(ylim_range)
    }
    
    return(p)
}
