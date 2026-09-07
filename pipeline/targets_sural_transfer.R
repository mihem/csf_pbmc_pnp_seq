targets_sural_transfer <- list(
  tar_target(
    sural_study_transfer_config,
    list(
      max_cells_per_patient_group = 50L,
      max_cells_per_cluster_tissue = 500L,
      nfeatures = 3000L,
      dimensions = 30L,
      confidence_threshold = 0.5,
      margin_threshold = 0.1,
      sural_tnk_clusters = c(
        "CD4", "Treg", "MAIT", "CD4_CD8", "CD8", "NK_CD8", "NK"
      ),
      sural_cd8_clusters = c("CD8", "NK_CD8"),
      study_tnk_clusters = c(
        "CD4naive_1", "CD4TCM_1", "CD4TCM_2", "CD4TEM", "CD4CTL",
        "Treg", "CD8naive", "CD8TCM", "CD8TEM_1", "CD8TEM_2",
        "CD8TEM_3", "MAIT", "gdT", "NKCD56bright_1",
        "NKCD56bright_2", "NKCD56dim"
      ),
      seed = 20260906L
    )
  ),
  tar_target(
    sural_study_transfer_input_file,
    "raw/sural/ic_study_transfer.qs",
    format = "file"
  ),
  tar_target(
    sural_study_transfer_query,
    read_sural_study_transfer_query(sural_study_transfer_input_file)
  ),
  tar_target(
    sural_tnk_study_transfer,
    run_sural_study_label_transfer(
      sc_annotated[
        ,
        as.character(sc_annotated$cluster) %in%
          sural_study_transfer_config$study_tnk_clusters
      ],
      sural_study_transfer_query[
        ,
        as.character(sural_study_transfer_query$ic_cluster) %in%
          sural_study_transfer_config$sural_tnk_clusters
      ],
      sural_study_transfer_config$max_cells_per_patient_group,
      sural_study_transfer_config$max_cells_per_cluster_tissue,
      sural_study_transfer_config$nfeatures,
      sural_study_transfer_config$dimensions,
      sural_study_transfer_config$confidence_threshold,
      sural_study_transfer_config$margin_threshold,
      sural_study_transfer_config$seed
    )
  ),
  tar_target(
    sural_cd8_study_transfer,
    run_sural_study_label_transfer(
      sc_annotated[
        ,
        as.character(sc_annotated$cluster) %in%
          sural_study_transfer_config$study_tnk_clusters
      ],
      sural_study_transfer_query[
        ,
        as.character(sural_study_transfer_query$ic_cluster) %in%
          sural_study_transfer_config$sural_cd8_clusters
      ],
      sural_study_transfer_config$max_cells_per_patient_group,
      sural_study_transfer_config$max_cells_per_cluster_tissue,
      sural_study_transfer_config$nfeatures,
      sural_study_transfer_config$dimensions,
      sural_study_transfer_config$confidence_threshold,
      sural_study_transfer_config$margin_threshold,
      sural_study_transfer_config$seed
    )
  ),
  tar_target(
    sural_study_transfer_workbook_file,
    write_sural_study_transfer_workbook(
      sural_tnk_study_transfer, sural_cd8_study_transfer
    ),
    format = "file"
  ),
  tar_target(
    sural_study_transfer_plot_files,
    write_sural_study_transfer_plots(
      sural_tnk_study_transfer,
      sural_cd8_study_transfer,
      sural_study_transfer_query
    ),
    format = "file"
  )
)
