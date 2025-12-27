##################################################
# PCA plots of cluster abundance across samples
# requires running annotate.R first
# not part of final analysis
##################################################
sc_merge_csf_main <- subset(
  sc_merge_csf,
  subset = diagnosis %in% c("CIDP", "GBS", "CIAP", "CTRL")
)

scMisc::pcaSeurat(
  object = sc_merge_csf_main,
  cluster = "cluster",
  sample = "sample",
  condition = "diagnosis",
  width = 20,
  height = 5,
  dir_output = file.path("results", "abundance")
)


props1 <- speckle::getTransformedProps(
  clusters = sc_merge_csf_main$cluster,
  sample = sc_merge_csf_main$sample,
  transform = "logit"
)

pcaSeurat1 <- function(
  object,
  cluster,
  sample,
  condition,
  width = 20,
  height = 5,
  dir_output = "."
) {
  if (!inherits(object, "Seurat")) {
    stop("Object must be a Seurat object")
  }
  object_parse <- deparse(substitute(object))

  props <- speckle::getTransformedProps(
    clusters = object@meta.data[[cluster]],
    sample = object@meta.data[[sample]],
    transform = "logit"
  )

  cl_size <- t(as.data.frame.matrix(props$TransformedProps))

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
      repel = TRUE
    ) +
    labs(title = "") +
    theme_classic()

  lookup_pre <-
    data.frame(
      cluster = object@meta.data[[sample]],
      condition = object@meta.data[[condition]]
    ) |>
    distinct()

  lookup <-
    data.frame(cluster = rownames(cl_size)) |>
    left_join(lookup_pre)

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
      fill.ind = lookup$condition,
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
      paste0(object_parse, "_", condition, "_", cluster, ".pdf")
    ),
    width = width,
    height = height,
    plot = pca_plots
  )
  return(pca_plots)
}

pcaSeurat1(
  object = sc_merge_csf_main,
  cluster = "cluster",
  sample = "sample",
  condition = "diagnosis",
  width = 20,
  height = 5,
  dir_output = file.path("results", "abundance")
)

lookup |>
  select(patient, diagnosis) |>
  print(n = Inf)

