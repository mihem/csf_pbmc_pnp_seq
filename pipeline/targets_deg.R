targets_deg <- tarchetypes::tar_map(
  values = deg_comparisons(),
  names = comparison,
  unlist = TRUE,
  targets::tar_target(
    deg_subset,
    subset_deg_object(sc_annotated, condition_1, condition_2, tissue)
  ),
  targets::tar_target(
    deg_cluster_pseudobulk,
    make_deg_cluster_pseudobulk(deg_subset)
  ),
  targets::tar_target(
    deg_cluster_name,
    deg_cluster_names(deg_cluster_pseudobulk, deg_subset@misc$cluster_order),
    iteration = "vector"
  ),
  targets::tar_target(
    deg_cluster_result,
    run_deg_limma(
      deg_cluster_pseudobulk[[deg_cluster_name]],
      lookup,
      condition_1,
      condition_2
    ),
    pattern = map(deg_cluster_name),
    iteration = "list"
  ),
  targets::tar_target(
    deg_cluster_results,
    collect_deg_cluster_results(deg_cluster_name, deg_cluster_result)
  ),
  targets::tar_target(
    deg_combined_result,
    run_deg_combined(deg_subset, lookup, condition_1, condition_2)
  ),
  targets::tar_target(
    deg_eligibility,
    deg_cluster_eligibility(deg_subset, deg_cluster_pseudobulk)
  ),
  targets::tar_target(deg_paths, deg_output_paths(comparison)),
  targets::tar_target(
    deg_cluster_workbook,
    write_deg_cluster_workbook(
      deg_cluster_results,
      deg_paths[["cluster_workbook"]]
    ),
    format = "file"
  ),
  targets::tar_target(
    deg_combined_workbook,
    write_deg_combined_workbook(
      deg_combined_result,
      deg_paths[["combined_workbook"]]
    ),
    format = "file"
  ),
  targets::tar_target(
    deg_count_plot,
    write_deg_count_plot(
      deg_cluster_results,
      deg_subset@misc$cluster_col,
      paste(condition_1, "vs", condition_2, tissue),
      deg_paths[["plot"]]
    ),
    format = "file"
  ),
  targets::tar_target(
    deg_eligibility_file,
    write_deg_eligibility(deg_eligibility, deg_paths[["eligibility"]]),
    format = "file"
  )
)
