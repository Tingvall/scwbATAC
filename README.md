# scwbATAC

`scwbATAC` owns the `scwbatac::chromvar` and `scwbatac::hvp_select` methods for
scworkbench. Install the
pack and its Bioconductor dependencies, then activate its installed manifest
with `scworkbench::wb_operation_pack_activate()` before execution. The engine
retains the matrix, genome-provider, motif-library and optional cached-match
registries; the pack receives only bounded authenticated readers and writers.

`hvp_select` selects highly variable peaks from a count-like peak matrix. It
binarizes accessibility for scoring, optionally removes the linear association
with `log1p(detected peaks per cell)`, and supports inclusive minimum and
maximum detection-percentage filters. Depth adjustment changes only the
ranking statistic; the operation never residualizes or replaces the source
matrix. It emits the complete scoring ledger as `selected_features` and an
auditable `fit_summary`.

```r
result <- scworkbench::wb_session_run_op(
  session, bundle, "scwbatac::hvp_select",
  inputs = list(
    matrix = "ATAC_counts",
    n_features = 25000L,
    depth_adjustment = TRUE,
    min_detection_pct = 0.5,
    max_detection_pct = 95
  )
)
```

`chromvar` emits `deviation_z` and `deviation_raw` as complete cell-major
`cell_vector_family` artifacts and `motif_summary` as the third terminal in
one atomic closure. Finite values and `NA_real_` retain exact cell/key
positions across promotion and recipient reopen. The pack rejects provider,
library, match-set and genome-build disagreements before terminal publication.
The authenticated runtime-family grant is qualified through 12,000,000 rows
per family, covering the 3,730-cell by 2,957-motif (11,029,610-row) shape.
Finite subnormal chromVAR p-values remain numeric doubles and use the core
typed-table content identity when the narrower canonical-JSON identity cannot
represent them; they are not rounded, floored or converted to strings.

This is not a streaming method. It reads bounded typed matrix chunks and then
materializes the selected sparse peak-by-cell matrix for chromVAR. The
summary retains its materialized dimensions, nonzero count and R object size
as a bounded-resource observation, not a streaming claim. The
effective dependency versions and provider are recorded in transition
evidence; a serial BiocParallel setting does not guarantee deterministic
reruns across different dependency versions.

Historical unqualified core `chromvar` runs retain their recorded operation
identity. New executions use `scwbatac::chromvar`, and recipient artifact and
provenance inspection remains possible without activating this pack.

Run the source-worktree contract suite with `SCWBATAC_RUN_SUITE=1` and
`SCWBATAC_CORE_SOURCE` set to the exact compatible scworkbench worktree:

```sh
SCWBATAC_RUN_SUITE=1 SCWBATAC_CORE_SOURCE=/path/to/scworkbench \
  Rscript tests/run_contract_tests.R
```
