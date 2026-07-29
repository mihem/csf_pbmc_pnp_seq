targets_liana <- list(
  targets::tar_target(liana_markers, liana_olink_markers()),
  targets::tar_target(
    liana_resource,
    make_liana_resource(liana_markers),
    packages = c("dplyr", "liana")
  ),
  targets::tar_target(
    liana_object,
    prepare_liana_object(sc_annotated, seed = 123L, max_cells = 1000L),
    packages = c("Seurat", "SeuratObject", "withr")
  ),
  targets::tar_target(
    liana_results,
    run_liana_consensus(liana_object, liana_resource, seed = 123L),
    packages = c("liana", "withr")
  ),
  targets::tar_target(
    liana_aggregate,
    aggregate_liana_results(liana_results),
    packages = c("dplyr", "liana")
  ),
  targets::tar_target(
    liana_selected,
    selected_liana_results(liana_aggregate),
    packages = "dplyr"
  ),
  targets::tar_target(
    liana_results_file,
    write_liana_qs(
      liana_results,
      file.path(liana_result_dir(), "liana_results.qs")
    ),
    format = "file",
    packages = "qs"
  ),
  targets::tar_target(
    liana_aggregate_file,
    write_liana_aggregate(
      liana_aggregate,
      file.path(liana_result_dir(), "liana_aggregate.csv")
    ),
    format = "file",
    packages = "readr"
  ),
  targets::tar_target(
    liana_all_dotplot_file,
    write_liana_all_dotplot(
      liana_aggregate,
      file.path(liana_result_dir(), "liana_all_dotplot.pdf")
    ),
    format = "file",
    packages = c("ggplot2", "liana")
  ),
  targets::tar_target(
    liana_selected_dotplot_file,
    write_liana_selected_dotplot(
      liana_selected,
      file.path(liana_result_dir(), "liana_selected_dotplot.pdf")
    ),
    format = "file",
    packages = c("ggplot2", "liana")
  ),
  targets::tar_target(
    liana_selected_chord_file,
    write_liana_selected_chord(
      liana_aggregate,
      file.path(liana_result_dir(), "liana_chord_freq_selected.pdf")
    ),
    format = "file",
    packages = c("dplyr", "liana")
  ),
  targets::tar_target(
    liana_legacy_results_file,
    file.path("objects", "liana_results.qs"),
    format = "file"
  ),
  targets::tar_target(
    liana_regression_validation,
    validate_liana_results(liana_results, liana_legacy_results_file),
    packages = c("dplyr", "qs", "tibble")
  ),
  targets::tar_target(
    liana_regression_validation_file,
    write_liana_validation(
      liana_regression_validation,
      file.path(liana_result_dir(), "legacy_regression.csv")
    ),
    format = "file",
    packages = "readr"
  )
)
