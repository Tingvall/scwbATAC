#!/usr/bin/env Rscript

core <- Sys.getenv("SCWBATAC_OLD_CORE_SOURCE", unset = "")
output <- Sys.getenv("SCWBATAC_ORACLE_PATH", unset = "")
stopifnot(nzchar(core), nzchar(output), dir.exists(core))
revision <- trimws(system2("git", c("-C", shQuote(core), "rev-parse", "HEAD"),
  stdout = TRUE))
stopifnot(identical(revision,
  "26f14667128339567275b9d11e005e815aae6c76"))
devtools::load_all(core, quiet = TRUE)

root <- tempfile("scwbatac-old-core-oracle-")
bundle <- wb_init_bundle(root)
on.exit({try(wb_close_bundle(bundle), silent = TRUE)
  unlink(root, recursive = TRUE, force = TRUE)}, add = TRUE)
data("mini_counts", package = "chromVAR", envir = environment())
rr <- SummarizedExperiment::rowRanges(mini_counts)
cells <- colnames(mini_counts)
peak_ids <- paste0(as.character(GenomicRanges::seqnames(rr)), ":",
  as.integer(GenomicRanges::start(rr)), "-",
  as.integer(GenomicRanges::end(rr)))
counts <- SummarizedExperiment::assay(mini_counts, "counts")
dimnames(counts) <- list(peak_ids, cells)
wb_define_cells(bundle, cells)
wb_attach_matrix(bundle, "ATAC_counts", Matrix::Matrix(counts,
  sparse = TRUE), feature_type = "peak")
wb_matrix_set_genome(bundle, "ATAC_counts", organism = "human",
  genome_build = "hg19", genome_provider = "BSgenome.Hsapiens.UCSC.hg19")
wb_motif_library_import_jaspar(bundle,
  motif_library = "jaspar2020_core_human", jaspar_pkg = "JASPAR2020",
  species = "9606", collection = "CORE", organism = "human",
  genome_build = "hg19", motif_ids = c("MA0139.1", "MA1102.2", "MA0140.2"))
inputs <- list(matrix = "ATAC_counts", cell_set = "all", cell_ids = cells,
  motif_library = "jaspar2020_core_human", min_motif_peaks = 3L,
  motif_p_cutoff = 5e-05, motif_w = 7L,
  background_iterations = 3L, background_bs = 10L, background_w = 0.1)
set.seed(1701L)
result <- wb_session_run_op(wb_session_new(), bundle, "chromvar",
  inputs = inputs)
outputs <- result$run$outputs
stopifnot(identical(sort(names(outputs)),
  c("deviation_raw", "deviation_z", "motif_summary")))
selected <- wb_slice(bundle, "ATAC_counts", cells = cells)
oracle <- list(schema = "scwbatac.chromvar.old_core_oracle.v1",
  source_commit = revision, seed = 1701L, input = inputs,
  selected_matrix = list(dimensions = as.integer(dim(selected)),
    nonzero = as.integer(Matrix::nnzero(selected)),
    materialized_object_bytes = as.numeric(object.size(selected))),
  outputs = list(
    deviation_z = list(key = as.character(outputs$deviation_z$key),
      cell_id = as.character(outputs$deviation_z$cell_id),
      value = as.numeric(outputs$deviation_z$value),
      value_type = as.character(outputs$deviation_z$value_type)),
    deviation_raw = list(key = as.character(outputs$deviation_raw$key),
      cell_id = as.character(outputs$deviation_raw$cell_id),
      value = as.numeric(outputs$deviation_raw$value),
      value_type = as.character(outputs$deviation_raw$value_type)),
    motif_summary = outputs$motif_summary$table,
    effective_params = result$run$params,
    ancestry = result$run$inputs$parents))
saveRDS(oracle, output, version = 3)
cat("OLD_CORE_CHROMVAR_ORACLE cells=", length(cells),
  " rows=", length(oracle$outputs$deviation_z$key),
  " selected_bytes=", oracle$selected_matrix$materialized_object_bytes,
  "\n", sep = "")
