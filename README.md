# CSF/PBMC Polyneuropathy Study

This is the code of the corresponding manuscript (preprint/in preparation).

All scripts can be found in the `scripts/` folder.

## Overview

The analysis pipeline covers:

| Script | Description |
|---|---|
| `preprocess.R` | scRNA-seq preprocessing and QC |
| `integrate.R` | Batch correction and data integration |
| `annotate.R` | Cell type annotation and clustering |
| `azimuth.R` | Reference mapping with Azimuth |
| `abundance.R` | Differential cell abundance (MiloDE) |
| `deg.R` | Differential gene expression (limma) |
| `tcr.R` / `tcr_comparison.R` | TCR repertoire analysis and GLIPH2 motif clustering |
| `liana.R` | Cell-cell communication (LIANA) |
| `olink_prepare.R` / `olink_analyze.R` | Olink proteomics preparation and differential analysis |
| `flow_pre.R` / `flow.R` | Flow cytometry preprocessing and analysis |
| `correlation.R` | Correlation of immune cell abundance with clinical scores |
| `enrichment.R` | Pathway enrichment (MSigDB, Enrichr) |
| `demographics.R` | Patient demographics analysis |
| `pca.R` | PCA of proteomics/flow data |
| `projectil.R` | ProjecTILs reference mapping |
| `sukenikova.R` | Comparison with Sukenikova et al. dataset |

## Reproducibility

### renv
To ensure reproducibility, we used the *renv* package. To restore the environment from the *renv.lock* file, use:

```R
renv::restore()
```

## Docker
If you have trouble restoring the environment via *renv*, you can also use the *Docker* image.

```bash
docker pull mihem/csf_pbmc_pnp_seq:v0.5
```

## Questions?
If you have any questions, please contact us at [mheming.com](https://.mheming.com).