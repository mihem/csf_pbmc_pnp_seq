# Legacy Archive

This directory preserves the implementation and artifacts that predate the
hierarchical `targets` pipeline. Nothing below `legacy/` is read by the active
graph.

## Contents

- `scripts/`: original script-based analysis and exploratory code.
- `objects/`: historical serialized checkpoints not used by active targets.
- `results/`: historical outputs, former regression baselines, and files removed
  from the active `results/targets/` manifest under `targets-stale/`.
- `docker/`: obsolete Docker launchers.
- `runtime/`: generated logs and stale local tooling configuration.

The active implementation is in `R/`, `pipeline/`, and `_targets.R`. Current
outputs are written exclusively below `results/targets/`.

Two files that remain active dependencies were deliberately not archived:

- The exact annotation checkpoint is in `references/checkpoints/`.
- The external ProjecTILs atlas is in `references/projectil/`.

This archive was organized on 2026-07-29 after manuscript panel reproduction
was completed. Preserve it for provenance; do not source the scripts or mix the
archived outputs with current target outputs.
