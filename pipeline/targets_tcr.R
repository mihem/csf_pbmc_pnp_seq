targets_tcr <- list(
  tar_target(tcr_manifest, discover_tcr_inputs(file.path("raw", "tcr"))),
  tar_target(tcr_files, tcr_input_files(tcr_manifest), format = "file"),
  tar_target(
    tcr_contigs,
    import_tcr_contigs(tcr_manifest, tcr_files, donor_assignments, lookup)
  ),
  tar_target(combined_tcr, combine_tcr_contigs(tcr_contigs)),
  tar_target(sc_tcr, annotate_tcr_cells(sc_annotated, combined_tcr)),
  tar_target(tcr_report_tables, make_tcr_report_tables(sc_tcr)),
  tar_target(
    tcr_report_table_files,
    write_tcr_report_tables(tcr_report_tables),
    format = "file"
  ),
  tar_target(
    tcr_basic_plot_files,
    write_tcr_basic_plots(combined_tcr, sc_tcr, seed = 42L),
    format = "file"
  ),
  tar_target(
    tcr_shared_clone_plot_files,
    write_tcr_shared_clone_plots(sc_tcr, tcr_report_tables),
    format = "file"
  ),
  tar_target(
    tcr_alluvial_plot_files,
    write_tcr_alluvial_plots(sc_tcr, seed = 42L),
    format = "file"
  ),
  tar_target(
    tcr_abundance_table_files,
    write_tcr_abundance_tables(sc_tcr),
    format = "file"
  ),
  tar_target(
    tcr_abundance_plot_files,
    write_tcr_abundance_plots(sc_tcr),
    format = "file"
  ),
  tar_target(
    tcr_invariant_plot_files,
    write_tcr_invariant_plots(sc_tcr, seed = 42L),
    format = "file"
  ),
  tar_target(tcr_comparison_manifest, tcr_comparison_input_manifest()),
  tar_target(tcr_comparison_input_files, tcr_comparison_manifest, format = "file"),
  tar_target(
    tcr_comparison,
    prepare_tcr_comparison(sc_tcr, tcr_comparison_input_files, seed = 42L)
  ),
  tar_target(
    tcr_comparison_plot_files,
    write_tcr_comparison_plots(tcr_comparison, seed = 42L),
    format = "file"
  ),

  # Reviewer 4 comment 4: select clones on disease enrichment
  # rather than raw size, then reverse phenotype them.
  #
  # Two selections are computed. The specificity group is the primary one. The
  # per-clone Fisher test is retained because it is the evidence for why the
  # unit was moved: on this cohort it selects only single-patient clones and its
  # ranking is clone size, which is the selection the reviewer objected to.
  tar_target(
    tcr_disease_enriched_clones,
    select_gliph_enriched_clones(
      sc_tcr,
      cluster_membership = tcr_comparison$cluster_membership,
      diagnosis_enrichment = tcr_comparison$diagnosis_enrichment,
      targets = c("GBS", "CIDP"),
      fdr = 0.1,
      min_patients = 2L
    )
  ),
  tar_target(
    tcr_per_clone_enriched_clones,
    select_disease_enriched_clones(
      sc_tcr,
      targets = c("GBS", "CIDP"),
      background = c("CTRL", "CIAP"),
      min_cells = 3L,
      fdr = 0.1
    )
  ),
  tar_target(
    tcr_all_clones,
    reverse_phenotype_clone_universe(sc_tcr, tcr_disease_enriched_clones)
  ),
  tar_target(
    tcr_size_matched_clones,
    match_clones_by_size(tcr_all_clones, seed = 42L)
  ),
  tar_target(
    tcr_clone_cluster_composition,
    reverse_phenotype_cluster_composition(sc_tcr, tcr_size_matched_clones)
  ),
  # Estimated over repeated size-matched draws. A single draw is one
  # arbitrary sample: CD8TEM_3 lands at zero matched cells in 7 of 20 draws,
  # which sends its odds ratio to infinity and its p-value from 0.003 to 0.31
  # depending only on the seed.
  tar_target(
    tcr_clone_cluster_enrichment,
    reverse_phenotype_cluster_enrichment_resampled(
      sc_tcr,
      tcr_all_clones,
      n_resamples = 200L,
      seed = 42L
    )
  ),
  # The contrast runs over the CD8 effector clusters rather than CD8TEM_3
  # alone. The selected clones sit in CD8TEM_1, and inside CD8TEM_3 no patient
  # has both arms populated, so the narrower contrast is not estimable.
  tar_target(
    tcr_reverse_phenotype,
    reverse_phenotype_pseudobulk(
      sc_tcr,
      tcr_disease_enriched_clones,
      cluster_name = c("CD8TEM_1", "CD8TEM_2", "CD8TEM_3"),
      min_cells_per_arm = 5L
    )
  ),
  tar_target(
    tcr_reverse_phenotype_panel,
    reverse_phenotype_panel_scores(
      sc_tcr,
      tcr_disease_enriched_clones,
      cluster_name = c("CD8TEM_1", "CD8TEM_2", "CD8TEM_3"),
      min_cells_per_arm = 5L
    )
  ),
  # Kept because the reply needs it as the power statement: the CD8TEM_3-only
  # contrast has no eligible patients.
  tar_target(
    tcr_reverse_phenotype_cd8tem3_cohort,
    reverse_phenotype_cohort(
      sc_tcr,
      tcr_disease_enriched_clones,
      cluster_name = "CD8TEM_3",
      min_cells_per_arm = 10L
    )
  ),
  tar_target(
    tcr_reverse_phenotype_file,
    write_reverse_phenotype_workbook(
      list(
        specificity_group = tcr_disease_enriched_clones,
        per_clone_fisher = tcr_per_clone_enriched_clones,
        cluster_composition = tcr_clone_cluster_composition,
        cluster_enrichment = tcr_clone_cluster_enrichment,
        pseudobulk = tcr_reverse_phenotype$results,
        cohort = tcr_reverse_phenotype$cohort,
        cd8tem3_cohort = tcr_reverse_phenotype_cd8tem3_cohort,
        panel_tests = tcr_reverse_phenotype_panel$tests,
        panel_paired = tcr_reverse_phenotype_panel$paired
      ),
      file.path(reverse_phenotype_result_dir(), "reverse_phenotype.xlsx")
    ),
    format = "file"
  ),
  tar_target(
    tcr_specificity_group_plot_file,
    write_specificity_group_plot(
      tcr_disease_enriched_clones,
      file.path(
        reverse_phenotype_result_dir(),
        "fig_specificity_group_patient_support.pdf"
      )
    ),
    format = "file"
  ),
  # Plots the per-clone test, not the specificity-group selection: the point of
  # the panel is that significance tracks clone size.
  tar_target(
    tcr_per_clone_selection_plot_file,
    write_disease_enriched_clone_plot(
      tcr_per_clone_enriched_clones,
      sc_tcr@misc$diagnosis_col,
      file.path(
        reverse_phenotype_result_dir(),
        "fig_per_clone_selection_tracks_size.pdf"
      )
    ),
    format = "file"
  ),
  tar_target(
    tcr_clone_cluster_composition_plot_file,
    write_clone_cluster_composition_plot(
      tcr_clone_cluster_composition,
      sc_tcr@misc$cluster_col,
      file.path(
        reverse_phenotype_result_dir(), "fig_clone_cluster_composition.pdf"
      )
    ),
    format = "file"
  ),
  tar_target(
    tcr_clone_cluster_enrichment_plot_file,
    write_cluster_enrichment_plot(
      tcr_clone_cluster_enrichment,
      file.path(
        reverse_phenotype_result_dir(), "fig_clone_cluster_enrichment.pdf"
      ),
      highlight = "CD8TEM_3"
    ),
    format = "file"
  ),
  tar_target(
    tcr_reverse_phenotype_plot_file,
    write_reverse_phenotype_plot(
      tcr_reverse_phenotype,
      file.path(
        reverse_phenotype_result_dir(), "fig_reverse_phenotype_volcano.pdf"
      ),
      seed = 42L
    ),
    format = "file"
  ),
  tar_target(
    tcr_reverse_phenotype_panel_plot_file,
    write_reverse_phenotype_panel_plot(
      tcr_reverse_phenotype_panel,
      file.path(
        reverse_phenotype_result_dir(), "fig_reverse_phenotype_panel.pdf"
      )
    ),
    format = "file"
  ),
  # Brain revision, Reviewer 4 comment 5: whether clones from a systemic trigger
  # are detectable in blood at all. Figure 5C cannot answer this because it is
  # computed only over clones already found in both compartments.
  tar_target(
    tcr_blood_expansion,
    blood_clonal_expansion(sc_tcr, top_n = 10L, min_cells = 50L)
  ),
  tar_target(
    tcr_blood_expansion_summary,
    blood_expansion_summary(tcr_blood_expansion)
  ),
  tar_target(
    tcr_blood_expansion_file,
    write_blood_expansion_workbook(
      tcr_blood_expansion,
      tcr_blood_expansion_summary,
      file.path(tcr_blood_result_dir(), "blood_clonal_expansion.xlsx")
    ),
    format = "file"
  ),
  tar_target(
    tcr_blood_expansion_plot_file,
    write_blood_expansion_plot(
      tcr_blood_expansion,
      sc_tcr@misc$diagnosis_col,
      file.path(tcr_blood_result_dir(), "fig_blood_clonal_expansion.pdf"),
      seed = 42L
    ),
    format = "file"
  ),
  tar_target(
    tcr_blood_expansion_gbs_plot_file,
    write_blood_expansion_plot(
      tcr_blood_expansion,
      sc_tcr@misc$diagnosis_col,
      file.path(tcr_blood_result_dir(), "fig_blood_clonal_expansion_gbs.pdf"),
      diagnoses = c("CTRL", "CIAP", "GBS", "CIDP"),
      seed = 42L
    ),
    format = "file"
  )
)
