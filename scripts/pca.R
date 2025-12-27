##################################################
# PCA plots of cluster abundance across samples
# requires running annotate.R first
##################################################

# libraries ---
library(Seurat)
library(tidyverse)
library(scMisc)
library(qs)
library(FactoMineR)
library(factoextra)

# load preprocessed data ----
sc_merge <- qs::qread(file.path("objects", "sc_merge.qs"), nthread = 4)

sc_merge_csf <- subset(sc_merge, tissue %in% "CSF")
sc_merge_pbmc <- subset(sc_merge, tissue %in% "PBMC")

sc_merge_csf_main <- subset(
  sc_merge_csf,
  subset = diagnosis %in% c("CIDP", "GBS", "CIAP", "CTRL")
)

sc_merge_pbmc_main <- subset(
  sc_merge_pbmc,
  subset = diagnosis %in% c("CIDP", "GBS", "CIAP", "CTRL")
)

pcaSeurat1 <- function(
  object_list,
  cluster,
  sample,
  patient,
  condition,
  tissue = NULL,
  object_name = NULL,
  width = 20,
  height = 5,
  dir_output = "."
) {
  # If single object, convert to list
  if (inherits(object_list, "Seurat")) {
    object_list <- list(object_list)
  }

  # Validate all objects are Seurat
  if (!all(sapply(object_list, inherits, "Seurat"))) {
    stop("All objects must be Seurat objects")
  }

  # Auto-capture object name if not provided
  if (is.null(object_name)) {
    object_name <- deparse1(substitute(object_list))
  }

  # Calculate proportions separately for each object
  if (length(object_list) == 1) {
    # Single tissue: patients as rows, clusters as columns
    obj <- object_list[[1]]
    props <- speckle::getTransformedProps(
      clusters = obj@meta.data[[cluster]],
      sample = obj@meta.data[[sample]],
      transform = "logit"
    )
    cl_size <- t(as.data.frame.matrix(props$TransformedProps))

    # Create lookup
    lookup_pre <- data.frame(
      sample = obj@meta.data[[sample]],
      patient = obj@meta.data[[patient]],
      condition = obj@meta.data[[condition]]
    ) |>
      distinct()

    lookup_final <- data.frame(sample = rownames(cl_size)) |>
      left_join(lookup_pre, by = "sample")
  } else {
    # Multiple tissues: calculate proportions per tissue, add tissue suffix
    prop_list <- lapply(object_list, function(obj) {
      tissue_name <- unique(obj@meta.data[[tissue]])[1]

      # Calculate proportions per sample (which is per patient)
      props <- speckle::getTransformedProps(
        clusters = obj@meta.data[[cluster]],
        sample = obj@meta.data[[sample]],
        transform = "logit"
      )

      # Transpose: samples as rows, clusters as columns
      cl_size_sample <- t(as.data.frame.matrix(props$TransformedProps))
      
      # Extract patient ID from sample names
      rownames(cl_size_sample) <- gsub(
        pattern = ".*_(P\\d+)$",
        replacement = "\\1",
        rownames(cl_size_sample)
      )

      # Add tissue suffix to cluster names
      props_df <- as.data.frame(cl_size_sample)
      names(props_df) <- paste0(names(props_df), "_", tissue_name)
      props_df$patient <- rownames(cl_size_sample)

      props_df
    })

    # Merge all tissues by patient
    cl_size_df <- Reduce(
      function(x, y) full_join(x, y, by = "patient"),
      prop_list
    ) |>
      as.data.frame()
    
    rownames(cl_size_df) <- cl_size_df$patient
    cl_size <- cl_size_df |> select(-patient) |> as.matrix()

    # Create lookup
    lookup_df <- do.call(
      rbind,
      lapply(object_list, function(obj) {
        obj@meta.data |> select(all_of(c(patient, condition))) |> distinct()
      })
    ) |>
      distinct() |>
      as.data.frame()
    names(lookup_df) <- c("patient", "condition")

    lookup_final <- data.frame(patient = rownames(cl_size)) |>
      left_join(lookup_df, by = "patient")
  }

  pca_result <- FactoMineR::PCA(
    cl_size,
    scale.unit = FALSE,
    ncp = 30,
    graph = FALSE
  )
  pca_eigen <- factoextra::fviz_eig(
    pca_result,
    addlabels = TRUE,
    ylim = c(0, 50),
    ncp = 7
  )

  pca_var_plot <-
    factoextra::fviz_pca_var(
      pca_result,
      col.var = "contrib",
      gradient.cols = viridis::viridis(100),
      repel = TRUE,
      select.var = list(contrib = 20)
    ) +
    labs(title = "") +
    theme_classic()

  factoextra::fviz_pca_ind(pca_result)

  pca_plot_ind <-
    factoextra::fviz_pca_ind(
      pca_result,
      pointsize = 5,
      pointshape = 21,
      fill = "#E7B800",
      col.ind = "black",
      palette = "Set2",
      axes.linetype = "solid"
    )

  pca_ggplot_ind <-
    ggpubr::ggpar(
      pca_plot_ind,
      title = "",
      xlab = "PC1",
      ylab = "PC2",
      ggtheme = theme_bw() +
        theme(
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(size = 25)
        )
    ) +
    theme_classic()

  pca_plot_group <-
    factoextra::fviz_pca_ind(
      pca_result,
      pointsize = 5,
      pointshape = 21,
      geom.ind = "point",
      fill.ind = lookup_final$condition,
      col.ind = "black",
      palette = "Set2",
      addEllipses = TRUE,
      ellipse.type = "confidence",
      legend.title = "group",
      axes.linetype = "solid"
    )

  pca_ggplot_group <-
    ggpubr::ggpar(
      pca_plot_group,
      title = "",
      xlab = "PC1",
      ylab = "PC2",
      ggtheme = theme_bw() +
        theme(
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(size = 25)
        )
    ) +
    theme_classic()

  pca_plots <- patchwork::wrap_plots(
    pca_eigen,
    pca_var_plot,
    pca_ggplot_ind,
    pca_ggplot_group,
    ncol = 4
  )
  ggsave(
    file.path(
      dir_output,
      paste0(object_name, "_", condition, "_", cluster, ".pdf")
    ),
    width = width,
    height = height,
    plot = pca_plots
  )
  return(pca_plots)
}

# Run PCA on combined CSF and PBMC data

# PCA for CSF only
pcaSeurat1(
  object_list = sc_merge_csf_main,
  cluster = "cluster",
  sample = "sample",
  patient = "patient",
  condition = "diagnosis",
  object_name = "csf",
  width = 20,
  height = 5,
  dir_output = file.path("results", "abundance")
)

# PCA for PBMC only
pcaSeurat1(
  object_list = sc_merge_pbmc_main,
  cluster = "cluster",
  sample = "sample",
  patient = "patient",
  condition = "diagnosis",
  object_name = "pbmc",
  width = 20,
  height = 5,
  dir_output = file.path("results", "abundance")
)

# PCA for combined CSF and PBMC
pcaSeurat1(
  object_list = list(sc_merge_csf_main, sc_merge_pbmc_main),
  cluster = "cluster",
  sample = "sample",
  patient = "patient",
  condition = "diagnosis",
  tissue = "tissue",
  object_name = "combined_csf_pbmc",
  width = 20,
  height = 5,
  dir_output = file.path("results", "abundance")
)
