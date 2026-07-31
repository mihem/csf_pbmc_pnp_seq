# CSF/PBMC Polyneuropathy Sequencing Study

Reproducible analysis of paired CSF and PBMC single-cell RNA/TCR sequencing,
flow cytometry, and Olink proteomics using an R `targets` pipeline.

## Run

```bash
docker pull mihem/csf_pbmc_pnp_seq:v0.6
docker compose run --rm csf_pbmc_pnp_seq
```

The container uses a 48 GiB memory limit by default. Pipeline outputs are
written to `results/targets/`, while intermediate objects are stored in
`_targets/`.

## Project Structure

- `R/`: analysis and plotting functions
- `pipeline/`: target definitions grouped by analysis stage
- `config/`: analysis parameters
- `_targets.R`: pipeline entry point
- `docker-compose.yml`: containerized pipeline services

## Pipeline DAG

```bash
docker compose run --rm dag
```

This writes the self-contained graph to `results/pipeline_dag.html`.
