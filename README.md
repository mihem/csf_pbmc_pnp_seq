# CSF/PBMC Polyneuropathy Sequencing Study

The reproducible workflow is a single hierarchical `targets` graph. Pipeline
functions live in `R/` and stage definitions live in `pipeline/`.

## Overview

The analysis pipeline covers:

| Script | Description |
|---|---|
| `preprocess.R` | scRNA-seq preprocessing and QC |
| `integrate.R` | Batch correction and data integration |
| `annotate.R` | Cell type annotation and clustering |
| `azimuth.R` | Reference mapping with Azimuth |
| `abundance.R` | Differential cell abundance |
| `deg.R` | Differential gene expression (limma) |
| `tcr.R` / `tcr_comparison.R` | TCR repertoire analysis and GLIPH2 motif clustering |
| `liana.R` | Cell-cell communication (LIANA) |
| `olink_prepare.R` / `olink_analyze.R` | Olink proteomics preparation and differential analysis |
| `flow.R` | Flow cytometry preprocessing and analysis |
| `correlation.R` | Correlation of immune cell abundance with clinical scores |
| `enrichment.R` | Pathway enrichment (MSigDB, Enrichr) |
| `demographics.R` | Patient demographics analysis |
| `pca.R` | PCA of proteomics/flow data |
| `projectil.R` | ProjecTILs reference mapping |
| `sukenikova.R` | Comparison with Sukenikova et al. dataset |

## Pipeline

The default graph covers preprocessing, Azimuth mapping, STACAS integration,
annotation, abundance, DEG, enrichment, LIANA, TCR, flow cytometry, Olink,
clinical correlations, demographics, Sukenikova comparison, and ProjectIL
readiness. Outputs are written below `results/targets/`; large intermediate
objects are managed by the `_targets/` store.

Inspect or run the graph from the project root:

```bash
Rscript -e 'targets::tar_manifest()'
Rscript -e 'targets::tar_make()'
Rscript -e 'targets::tar_visnetwork()'
```

Export the complete target DAG as a self-contained interactive HTML file:

```bash
LOCAL_UID=$(id -u) LOCAL_GID=$(id -g) docker compose run --rm dag
```

The graph is written to `results/pipeline_dag.html`. Nodes are colored by
their current `targets` status, and generating the graph does not run the
pipeline. To include the optional Cerebro branch, set
`TAR_INCLUDE_HEAVY=true` for the Compose command.

Run a bounded stage instead of the full graph:

```bash
Rscript -e 'targets::tar_make(names = tidyselect::starts_with("deg_"))'
Rscript -e 'targets::tar_make(names = tidyselect::starts_with("flow_"))'
```

Historical regression checks are not part of the active graph. The graph has
two explicit reference dependencies:

- `references/checkpoints/sc_merge_batch.qs` supplies the validated historical
  STACAS and UMAP reductions required for exact historical annotation.
- `references/projectil/conventional_cd8_projectil_reference.qs` is the
  external Terekhova CD8 atlas used by ProjecTILs.

Their provenance and checksums are documented in `references/README.md`.
The unused `legacy` block in `config/analysis.yml` is retained only to avoid
invalidating the cached memory-intensive STACAS target; no active command reads
those paths. It can be removed when `_targets/` is intentionally rebuilt from
scratch.

Memory-intensive Cerebro export targets are excluded by default. Enable them
explicitly with `TAR_INCLUDE_HEAVY=true`. The overlay-heavy TCR plot bundle is
not in the active graph because `clonalOverlay()` exceeded workstation memory;
the bounded TCR tables, comparisons, alluvial plots, basic plots, and invariant
plots remain active.

## Reproducibility

### renv
To ensure reproducibility, we used the *renv* package. To restore the environment from the *renv.lock* file, use:

```r
renv::restore()
```

## Docker
The Docker image includes the locked R environment and Azimuth PBMC reference.

```bash
docker pull mihem/csf_pbmc_pnp_seq:v0.6
LOCAL_UID=$(id -u) LOCAL_GID=$(id -g) docker compose run --rm csf_pbmc_pnp_seq
```

Compose limits the container to 48 GiB by default so a single analysis cannot
exhaust the workstation. Override this only when necessary:

```bash
PIPELINE_MEMORY_LIMIT=56g LOCAL_UID=$(id -u) LOCAL_GID=$(id -g) \
  docker compose run --rm csf_pbmc_pnp_seq
```

Run a specific target set in a 16 GiB container:

```bash
docker run --rm --memory=16g --memory-swap=16g \
  --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -e RENV_CONFIG_SANDBOX_ENABLED=FALSE -v "$PWD:/csf_pbmc_pnp_seq" \
  -w /csf_pbmc_pnp_seq mihem/csf_pbmc_pnp_seq:v0.6 \
  Rscript -e 'targets::tar_make(names = "tcr_alluvial_plot_files")'
```

## Questions?
If you have any questions, please contact us at [mheming.com](https://.mheming.com).
