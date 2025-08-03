##################################################
# analyze correlations
##################################################
library(tidyverse)
library(qs)
library(Seurat)
library(rcna)
library(glue)
library(scMisc)
library(viridis)

# load data ----
lookup <- qs::qread(file.path("objects", "lookup.qs"))
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"))

# subset those that have the severity variable

sc_merge$incat_at_lumbar_puncture
sc_merge[["incat_at_lumbar_puncture"]]

rcna_plot <- function(seu_obj, param) {
    # check if file exists
    if (
        file.exists(file.path(
            "results",
            "correlation",
            paste0("rcna_", param, ".pdf")
        ))
    ) {
        return(paste0("File ", param, " already exists, skipping."))
    }
    seu_obj_filter <- seu_obj[, !is.na(seu_obj[[param]])]
    seu_obj_filter$sex_numeric <- as.numeric(seu_obj_filter$sex == "male")

    seu_obj_filter <- association.Seurat(
        seurat_object = seu_obj_filter,
        test_var = param,
        samplem_key = "sample",
        graph_use = "RNA_nn",
        verbose = TRUE,
        batches = NULL,
        # no batch variables to include, only works with matched design https://github.com/immunogenomics/cna/issues/11
        covs = c("sex_numeric", "age")
    )

    title <- gsub(pattern = "_", replacement = " ", x = param)
    rcna_fplot <-
        FeaturePlot(
            seu_obj_filter,
            features = c("cna_ncorrs"),
            reduction = "umap.stacas.ss.all",
            pt.size = 0.1,
            order = FALSE,
            coord.fixed = TRUE,
            raster = FALSE,
            alpha = 0.2
        ) +
        viridis::scale_color_viridis(option = "inferno") +
        theme_rect() +
        labs(
            title = title,
            color = "Correlation",
            x = "UMAP1",
            y = "UMAP2"
        )

    ggsave(
        plot = rcna_fplot,
        filename = file.path(
            "results",
            "correlation",
            paste0("rcna_", param, ".pdf")
        ),
        width = 8,
        height = 6
    )
}

rcna_params <- c(
    "incat_at_lumbar_puncture",
    "incat_follow_up",
    "incat_progress",
    "onls_at_lumbar_puncture",
    "onls_follow_up",
    "onls_progress",
    "mrc_sum_score_60_at_lumbar_puncture",
    "mrc_sum_score_60_follow_up",
    "mrc_sum_score_60_progress",
    "csf_protein",
    "dml_ulnar_motoric",
    "cmap_ulnar_motoric",
    "ncv_ulnar_motoric",
    "min_f_latency_ulnar_motoric",
    "cmap_tibial_motoric",
    "ncv_tibial_motoric"
)

lapply(
    rcna_params,
    function(x) {
        rcna_plot(
            seu_obj = sc_merge,
            param = x
        )
    }
)

# analyze subgroup with high NKCD56bright_1 cells ----
nkbright_high_group <- c("P16", "P17", "P18", "P21", "P23")

lookup_subgroup <- lookup |>
    dplyr::mutate(
        subgroup = dplyr::case_when(
            patient %in% nkbright_high_group ~ "NKhi",
            TRUE ~ as.character(diagnosis)
        )
    )

plot_severity_subgroup <-
    lookup_subgroup |>
    ggplot(
        aes(x = subgroup, y = incat_at_lumbar_puncture, fill = subgroup)
    ) +
    geom_boxplot() +
    geom_jitter(
        width = 0.2,
        height = 0,
        alpha = 0.5
    ) +
    theme_classic() +
    theme(legend.position = "none")

ggsave(
    plot = plot_severity_subgroup,
    filename = file.path("results", "severity", "severity_subgroup.pdf"),
    width = 6,
    height = 4
)
