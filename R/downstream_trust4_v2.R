prepare_sural_trust4_v2_analysis <- function(
  tcr_contigs, combined_tcr, mapped, patient_map, sc_tcr,
  enrichment_permutations = 10000L, enrichment_seed = 20260906L,
  enrichment_min_shared_clonotypes = 3L
) {
  beta <- prepare_sural_tcr_comparison(
    tcr_contigs, mapped$records, patient_map
  )
  chains <- prepare_sural_tcr_chain_comparisons(
    tcr_contigs, combined_tcr, mapped$records, patient_map
  )
  list(
    beta = beta,
    chains = chains,
    normalized_intersections = prepare_sural_tcr_normalized_intersections(
      beta, chains
    ),
    expansion = prepare_sural_tcr_expansion(beta, chains),
    enrichment = prepare_sural_tcr_cluster_enrichment(
      beta,
      chains,
      sc_tcr,
      enrichment_permutations,
      enrichment_seed,
      enrichment_min_shared_clonotypes
    ),
    cluster_compositions = list(
      trb = prepare_tcr_intersection_cluster_composition(
        beta, "TRB_CDR3aa", sc_tcr, mapped$records
      ),
      tra = prepare_tcr_intersection_cluster_composition(
        chains$alpha, "clonotype", sc_tcr, mapped$records
      )
    )
  )
}

prepare_trust4_version_comparison <- function(
  v1_mapped, v1_ic_mapped, v1_beta, v1_chains,
  v2_mapped, v2_ic_mapped, v2_analysis
) {
  mapping_summary <- dplyr::bind_rows(
    dplyr::mutate(v1_mapped$match_summary, version = "V1", scope = "all sural"),
    dplyr::mutate(v1_ic_mapped$match_summary, version = "V1", scope = "immune"),
    dplyr::mutate(v2_mapped$match_summary, version = "V2 R2-only", scope = "all sural"),
    dplyr::mutate(v2_ic_mapped$match_summary, version = "V2 R2-only", scope = "immune")
  ) |>
    dplyr::mutate(
      matched = dplyr::coalesce(.data$matched, 0L),
      unmatched = dplyr::coalesce(.data$unmatched, 0L),
      total = .data$matched + .data$unmatched,
      match_percent = 100 * .data$matched / .data$total,
      .before = "library_id"
    ) |>
    dplyr::arrange(.data$scope, .data$library_id, .data$receptor, .data$version)

  bind_intersections <- function(beta, chains, version) {
    dplyr::bind_rows(
      dplyr::mutate(beta$intersection_summary, chain = "TRB"),
      dplyr::mutate(chains$alpha$intersections, chain = "TRA"),
      dplyr::mutate(chains$paired$intersections, chain = "paired TRA + TRB")
    ) |>
      dplyr::mutate(version = version, .before = 1L)
  }
  intersections <- dplyr::bind_rows(
    bind_intersections(v1_beta, v1_chains, "V1"),
    bind_intersections(
      v2_analysis$beta, v2_analysis$chains, "V2 R2-only"
    )
  ) |>
    dplyr::arrange(
      .data$patient, .data$chain, .data$presence, .data$version
    )

  bind_sural_shared <- function(beta, chains, version) {
    dplyr::bind_rows(
      dplyr::transmute(
        beta$sural_shared,
        patient = as.character(.data$patient),
        chain = "TRB",
        clonotype = .data$TRB_CDR3aa
      ),
      dplyr::transmute(
        chains$alpha$sural_shared,
        patient = as.character(.data$patient),
        chain = "TRA",
        clonotype = .data$clonotype
      ),
      dplyr::transmute(
        chains$paired$sural_shared,
        patient = as.character(.data$patient),
        chain = "paired TRA + TRB",
        clonotype = .data$clonotype
      )
    ) |>
      dplyr::mutate(version = version, .before = 1L)
  }
  sural_shared <- dplyr::bind_rows(
    bind_sural_shared(v1_beta, v1_chains, "V1"),
    bind_sural_shared(
      v2_analysis$beta, v2_analysis$chains, "V2 R2-only"
    )
  ) |>
    dplyr::arrange(.data$patient, .data$chain, .data$clonotype, .data$version)
  patient_version_availability <- dplyr::bind_rows(
    tibble::tibble(
      version = "V1",
      patient = as.character(v1_beta$patient_mapping$patient)
    ),
    tibble::tibble(
      version = "V2 R2-only",
      patient = as.character(v2_analysis$beta$patient_mapping$patient)
    )
  ) |>
    dplyr::distinct()
  sural_shared_summary <- patient_version_availability |>
    tidyr::crossing(chain = c("TRA", "TRB", "paired TRA + TRB")) |>
    dplyr::left_join(
      sural_shared |>
        dplyr::count(
          .data$version, .data$patient, .data$chain,
          name = "clonotypes"
        ),
      by = c("version", "patient", "chain")
    ) |>
    dplyr::mutate(clonotypes = dplyr::coalesce(.data$clonotypes, 0L)) |>
    dplyr::arrange(.data$patient, .data$chain, .data$version)

  list(
    mapping_summary = mapping_summary,
    patient_version_availability = patient_version_availability,
    intersections = intersections,
    sural_shared_summary = sural_shared_summary,
    sural_shared_clonotypes = sural_shared
  )
}

write_trust4_version_comparison <- function(comparison) {
  path <- file.path(
    "results", "targets", "trust4_v2_r2", "trust4_v1_v2_comparison.xlsx"
  )
  ensure_parent_dir(path)
  writexl::write_xlsx(comparison, path)
  path
}

write_sural_trust4_v2_outputs <- function(
  mapped, ic_mapped, cells, ic_cells, analysis
) {
  temporary_root <- tempfile("trust4_v2_")
  dir.create(temporary_root, recursive = TRUE)
  on.exit(unlink(temporary_root, recursive = TRUE), add = TRUE)

  generated <- withr::with_dir(
    temporary_root,
    c(
      write_sural_trust4_table(mapped),
      write_sural_trust4_umap(cells, mapped),
      write_sural_trust4_cluster_summary(mapped),
      write_sural_trust4_table(
        ic_mapped, "trust4_immune_cells_with_subclusters.xlsx"
      ),
      write_sural_trust4_umap(
        ic_cells,
        ic_mapped,
        "trust4_immune_cells_umap.png",
        "Immune cells with TRUST4 receptor calls"
      ),
      write_sural_trust4_cluster_summary(
        ic_mapped,
        "trust4_receptor_cells_by_immune_subcluster.pdf",
        "TRUST4 receptor-positive cells by immune subcluster"
      ),
      write_sural_tcr_comparison_workbook(analysis$beta),
      write_sural_tcr_intersection_plot(analysis$beta),
      write_sural_tcr_tracking_plot(analysis$beta),
      write_sural_chain_comparison_workbook(analysis$chains),
      write_sural_chain_intersection_plot(
        analysis$chains$alpha,
        "TRA",
        "tra_shared_tissue_intersections.pdf"
      ),
      write_sural_chain_tracking_plot(
        analysis$chains$alpha,
        "TRA",
        "sural_shared_tra_clonotype_tracking.pdf"
      ),
      write_sural_chain_intersection_plot(
        analysis$chains$paired,
        "paired TRA + TRB",
        "paired_tra_trb_shared_tissue_intersections.pdf"
      ),
      write_sural_chain_tracking_plot(
        analysis$chains$paired,
        "paired TRA + TRB",
        "sural_shared_paired_tra_trb_tracking.pdf"
      ),
      write_sural_tcr_normalized_intersection_outputs(
        analysis$normalized_intersections
      ),
      write_sural_tcr_expansion_workbook(analysis$expansion),
      write_sural_tcr_expansion_plot(analysis$expansion),
      write_sural_tcr_cluster_enrichment_workbook(analysis$enrichment),
      write_sural_tcr_cluster_enrichment_plot(analysis$enrichment),
      write_tcr_intersection_cluster_plot(
        analysis$cluster_compositions$trb,
        "TRB clonotype",
        "trb_shared_tissue_intersections_by_cluster.pdf"
      ),
      write_tcr_intersection_cluster_plot(
        analysis$cluster_compositions$tra,
        "TRA clonotype",
        "tra_shared_tissue_intersections_by_cluster.pdf"
      )
    )
  )

  source_root <- file.path(temporary_root, "results", "targets", "trust4")
  source_paths <- file.path(temporary_root, generated)
  relative_paths <- substring(
    normalizePath(source_paths),
    nchar(normalizePath(source_root)) + 2L
  )
  destination_paths <- file.path(
    "results", "targets", "trust4_v2_r2", relative_paths
  )
  purrr::walk(destination_paths, ensure_parent_dir)
  copied <- file.copy(source_paths, destination_paths, overwrite = TRUE)
  stopifnot(all(copied), all(file.exists(destination_paths)))
  destination_paths
}
