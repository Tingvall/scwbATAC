# Prospective chromVAR comparison policy

Frozen before implementing the scwbATAC runner. Oracle: old-core scworkbench
`26f14667128339567275b9d11e005e815aae6c76`, chromVAR `mini_counts`,
JASPAR2020 human CORE motifs MA0139.1/MA1102.2/MA0140.2, seed 1701,
35 selected cells and three background iterations.

- Require exactly the three terminals `deviation_z`, `deviation_raw`, and
  `motif_summary`, in one promotion closure; no missing or extra output.
- Require exact ordered cell IDs and motif keys. Both families remain
  **cell-major**: for each cell, the full ordered motif-key sequence appears.
  Require the same continuous value type and row counts. The frozen old-core
  oracle contains five `NA` values in each family; require their exact
  positions and reject `NaN`/infinity. This sentence corrects the initial
  no-nonfinite draft after inspecting the frozen oracle, before accepting any
  pack result. No value is imputed merely to satisfy a new ABI.
- Compare each paired deviation value within absolute tolerance `1e-6` and
  relative tolerance `1e-6`. This is numerical equivalence, not a promise of
  byte-identical payloads or deterministic reruns on other dependency versions.
- Require exact summary motif identifiers, library/source labels, matched-peak
  counts, matrix/genome labels and row order. Compare floating summary columns
  (`fraction_matches`, background overlap, variability, adjusted/unadjusted
  P-values, mean absolute z) with the same `1e-6` absolute/relative tolerance;
  compare missingness exactly.
- Require effective parameter meaning and values, selected matrix/cell/motif
  bindings, output closure and source ancestry. Operation identity is expected
  to change from the legacy core name to `scwbatac::chromvar`; old content
  hashes are not an acceptance criterion for newly executed pack outputs.
- Compare a live motifmatchr path and an authenticated cached match-set path
  separately. The old-core oracle freezes the live positive; cached-path
  equivalence uses the same prepared source fixture and policy.

Resource observation: record dimensions, nonzeros and `object.size()` for the
selected sparse matrix after materialization. It is **not** a streaming test.
