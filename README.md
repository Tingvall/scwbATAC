# scwbATAC

`scwbATAC` owns the `scwbatac::chromvar` method for scworkbench. Install the
pack and its Bioconductor dependencies, then activate its installed manifest
with `scworkbench::wb_operation_pack_activate()` before execution. The engine
retains the matrix, genome-provider, motif-library and optional cached-match
registries; the pack receives only bounded authenticated readers and writers.

The method emits `deviation_z` and `deviation_raw` as complete cell-major
`cell_vector_family` artifacts and `motif_summary` as the third terminal in
one atomic closure. Finite values and `NA_real_` retain exact cell/key
positions across promotion and recipient reopen. The pack rejects provider,
library, match-set and genome-build disagreements before terminal publication.

This is not a streaming method. It reads bounded typed matrix chunks and then
materializes the selected sparse peak-by-cell matrix for chromVAR. The
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
