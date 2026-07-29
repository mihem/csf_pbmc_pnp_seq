# Analysis References

The active graph requires two large binary references that are intentionally
ignored by Git. They must be restored at the paths below before running the
pipeline in a fresh checkout.

## Annotation Checkpoint

`checkpoints/sc_merge_batch.qs` is a historical Seurat checkpoint created after
STACAS integration. The pipeline still recomputes preprocessing, Azimuth, and
STACAS from study inputs, then copies the checkpoint's `stacas.ss.all` and
`umap.stacas.ss.all` reductions plus `group` and `age` metadata onto the
recomputed object before clustering.

This boundary is required for exact manuscript annotation. A clean STACAS
recomputation was biologically equivalent but produced 33 resolution-1 clusters
instead of the historical 32.

- Size: approximately 1.9 GiB
- SHA-256: `d583a81efcc3220375a401e0e73374507244ba7df7f18aeb92076298c63c30d8`

## ProjecTILs Reference

`projectil/conventional_cd8_projectil_reference.qs` is the external Terekhova
conventional CD8 reference atlas used to project the study's `CD8TEM_3` cells.
It is not derivable from this study's raw data. The historical repository did
not retain a download or construction recipe, so this exact serialized atlas is
currently the reproducible input.

- Size: approximately 246 MiB
- SHA-256: `e5d800ca987ddf58c3494320384f69ec7d7e89cdd7dcc8fd2db9b4bfe375f93e`

The Azimuth PBMC reference is a separate external dependency bundled directly
in Docker image `mihem/csf_pbmc_pnp_seq:v0.6`.
