.atac_test_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
.atac_core_root <- Sys.getenv("SCWBATAC_CORE_SOURCE", unset = "")
if (nzchar(.atac_core_root)) {
  devtools::load_all(.atac_core_root, quiet = TRUE)
  devtools::load_all(.atac_test_root, quiet = TRUE)
} else {
  library(scworkbench)
  library(scwbATAC)
}

.atac_test_fixture <- function() {
  root <- tempfile("scwbatac-pack-")
  bundle <- scworkbench::wb_init_bundle(root)
  data("mini_counts", package = "chromVAR", envir = environment())
  rr <- SummarizedExperiment::rowRanges(mini_counts)
  cells <- colnames(mini_counts)
  ids <- paste0(as.character(GenomicRanges::seqnames(rr)), ":",
    as.integer(GenomicRanges::start(rr)), "-",
    as.integer(GenomicRanges::end(rr)))
  counts <- SummarizedExperiment::assay(mini_counts, "counts")
  dimnames(counts) <- list(ids, cells)
  scworkbench::wb_define_cells(bundle, cells)
  scworkbench::wb_attach_matrix(bundle, "ATAC_counts",
    Matrix::Matrix(counts, sparse = TRUE), feature_type = "peak")
  scworkbench::wb_matrix_set_genome(bundle, "ATAC_counts",
    organism = "human", genome_build = "hg19",
    genome_provider = "BSgenome.Hsapiens.UCSC.hg19")
  scworkbench::wb_motif_library_import_jaspar(bundle,
    motif_library = "jaspar2020_core_human", jaspar_pkg = "JASPAR2020",
    species = "9606", collection = "CORE", organism = "human",
    genome_build = "hg19", motif_ids = c("MA0139.1", "MA1102.2",
      "MA0140.2"))
  list(root = root, bundle = bundle, cells = cells, features = ids)
}

.atac_test_close <- function(fixture) {
  try(scworkbench::wb_close_bundle(fixture$bundle), silent = TRUE)
  unlink(fixture$root, recursive = TRUE, force = TRUE)
}

.atac_compare_oracle <- function(result, oracle) {
  stopifnot(identical(sort(names(result$run$outputs)),
    c("deviation_raw", "deviation_z", "motif_summary")))
  for (name in c("deviation_z", "deviation_raw")) {
    observed <- result$run$outputs[[name]]
    expected <- oracle$outputs[[name]]
    stopifnot(identical(as.character(observed$key), expected$key),
      identical(as.character(observed$cell_id), expected$cell_id),
      identical(which(is.na(observed$value)),
        which(is.na(expected$value))),
      !any(is.nan(observed$value) | is.infinite(observed$value)),
      isTRUE(all.equal(as.numeric(observed$value), expected$value,
        tolerance = 1e-6)))
  }
  observed <- result$run$outputs$motif_summary$table
  expected <- oracle$outputs$motif_summary
  labels <- c("motif_key", "motif_id", "motif_name", "motif_library",
    "collection", "source_id", "source_version", "genome_pkg",
    "genome_build", "organism", "matched_peaks")
  numerics <- c("fraction_matches", "fraction_background_overlap",
    "variability", "variability_p_value", "variability_p_value_adj",
    "mean_abs_z")
  for (field in labels) if (!identical(as.character(observed[[field]]),
      as.character(expected[[field]]))) stop("summary label mismatch: ",
      field, " observed=", paste(observed[[field]], collapse = ","),
      " expected=", paste(expected[[field]], collapse = ","))
  for (field in numerics) stopifnot(identical(
    which(is.na(observed[[field]])), which(is.na(expected[[field]]))),
    isTRUE(all.equal(as.numeric(observed[[field]]),
      as.numeric(expected[[field]]), tolerance = 1e-6)))
  TRUE
}

.atac_cached_match_set <- function(fixture) {
  data("mini_counts", package = "chromVAR", envir = environment())
  ranges <- SummarizedExperiment::rowRanges(mini_counts)
  names(ranges) <- fixture$features
  motifs <- scworkbench:::.wb_motifs_load_tfbstools(fixture$bundle,
    motif_library = "jaspar2020_core_human")
  provider <- getExportedValue("BSgenome.Hsapiens.UCSC.hg19",
    "BSgenome.Hsapiens.UCSC.hg19")
  matches <- motifmatchr::matchMotifs(pwms = motifs, subject = ranges,
    genome = provider, p.cutoff = 5e-5, w = 7L)
  assay <- SummarizedExperiment::assay(matches)
  positions <- which(as.matrix(assay), arr.ind = TRUE)
  stopifnot(nrow(positions) > 0L)
  scworkbench::wb_motif_match_replace(fixture$bundle,
    motif_match_set = "cached_jaspar", matrix_id = "ATAC_counts",
    motif_library = "jaspar2020_core_human",
    matches = data.frame(
      motif_key = colnames(assay)[positions[, 2L]],
      feature_id = rownames(assay)[positions[, 1L]],
      stringsAsFactors = FALSE),
    organism = "human", genome_build = "hg19",
    genome_provider = "BSgenome.Hsapiens.UCSC.hg19",
    match_params = list(p_cutoff = 5e-5, w = 7L))
  invisible(TRUE)
}

test_chromvar_live_oracle <- function() {
  fixture <- .atac_test_fixture()
  on.exit(.atac_test_close(fixture), add = TRUE)
  manifest <- system.file("scworkbench", "operation-pack.json",
    package = "scwbATAC")
  try(scworkbench::wb_operation_pack_unload("scwbatac"), silent = TRUE)
  missing_pack <- try(scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), fixture$bundle,
    "scwbatac::chromvar", inputs = list(matrix = "ATAC_counts",
      motif_library = "jaspar2020_core_human", cell_set = "all")),
    silent = TRUE)
  stopifnot(inherits(missing_pack, "try-error"),
    nrow(scworkbench::wb_art_list(fixture$bundle)) == 0L)
  scworkbench::wb_operation_pack_activate(manifest)
  stopifnot(sum(scworkbench::wb_list_ops() ==
    "scwbatac::chromvar") == 1L)
  set.seed(1701L)
  result <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    fixture$bundle, "scwbatac::chromvar", inputs = list(
      matrix = "ATAC_counts", motif_library = "jaspar2020_core_human",
      min_motif_peaks = 3L, background_iterations = 3L,
      background_bs = 10L, chunk_size = 35L, cell_set = "all"))
  oracle <- readRDS(file.path("tests", "fixtures",
    "old_core_chromvar_oracle.rds"))
  stopifnot(.atac_compare_oracle(result, oracle))
  stopifnot(is.finite(result$run$outputs$motif_summary$meta$
    selected_materialized_matrix_bytes),
    result$run$outputs$motif_summary$meta$
      selected_materialized_matrix_bytes > 0)
  source_counts <- c(artifacts = nrow(scworkbench::wb_art_list(fixture$bundle)),
    runs = nrow(DBI::dbGetQuery(fixture$bundle$con,
      "SELECT * FROM prov_runs;")),
    commits = nrow(DBI::dbGetQuery(fixture$bundle$con,
      "SELECT * FROM prov_commits;")))
  source_files <- sort(list.files(fixture$root, all.files = TRUE,
    recursive = TRUE, no.. = TRUE))
  incomplete <- try(scworkbench::wb_session_promote_current(fixture$bundle,
    result$session, outputs = "deviation_z"), silent = TRUE)
  stopifnot(inherits(incomplete, "try-error"),
    grepl("complete terminal closure", as.character(incomplete),
      fixed = TRUE))
  old_hook <- getOption("scworkbench.promotion_transaction_test_hook")
  options(scworkbench.promotion_transaction_test_hook = function(event, ...) {
    if (identical(event, "after_first_private_terminal_write"))
      stop("injected chromVAR first-terminal fault")
  })
  failed <- try(scworkbench::wb_session_promote_current(fixture$bundle,
    result$session, outputs = c("deviation_z", "deviation_raw",
      "motif_summary")), silent = TRUE)
  options(scworkbench.promotion_transaction_test_hook = old_hook)
  retained_counts <- c(artifacts = nrow(scworkbench::wb_art_list(fixture$bundle)),
    runs = nrow(DBI::dbGetQuery(fixture$bundle$con,
      "SELECT * FROM prov_runs;")),
    commits = nrow(DBI::dbGetQuery(fixture$bundle$con,
      "SELECT * FROM prov_commits;")))
  stopifnot(inherits(failed, "try-error"),
    grepl("injected chromVAR first-terminal fault", as.character(failed),
      fixed = TRUE), identical(retained_counts, source_counts),
    identical(sort(list.files(fixture$root, all.files = TRUE,
      recursive = TRUE, no.. = TRUE)), source_files))
  promoted <- scworkbench::wb_session_promote_current(fixture$bundle,
    result$session, outputs = c("deviation_z", "deviation_raw",
      "motif_summary"),
    created_by = "phase7-chromvar-pack")
  stopifnot(all(promoted$promoted$status == "ok"),
    setequal(as.character(promoted$promoted$output),
      c("deviation_z", "deviation_raw", "motif_summary")))
  artifact_ids <- stats::setNames(as.character(promoted$promoted$artifact_id),
    as.character(promoted$promoted$output))
  for (name in c("deviation_z", "deviation_raw")) {
    retained <- scworkbench::wb_art_load(fixture$bundle,
      artifact_ids[[name]], expected_kind = "cell_vector_family")
    stopifnot(identical(as.character(retained$key),
      oracle$outputs[[name]]$key),
      identical(as.character(retained$cell_id),
        oracle$outputs[[name]]$cell_id),
      identical(which(is.na(retained$value)),
        which(is.na(oracle$outputs[[name]]$value))),
      identical(retained$value_type, result$run$outputs[[name]]$value_type),
      identical(retained$meta, result$run$outputs[[name]]$meta),
      identical(scworkbench:::.wb_prov_payload_identity(
        fixture$bundle, retained)$content_sha256,
        scworkbench:::.wb_prov_payload_identity(fixture$bundle,
          result$run$outputs[[name]])$content_sha256))
  }
  stopifnot(isTRUE(scworkbench::wb_finalize_recipient_bundle(
    fixture$bundle)$ready))
  scworkbench::wb_close_bundle(fixture$bundle)
  relocated <- tempfile("scwbatac-relocated-")
  on.exit(unlink(relocated, recursive = TRUE, force = TRUE), add = TRUE)
  fs::dir_copy(fixture$root, relocated)
  fixture$bundle <- scworkbench::wb_open_bundle(relocated,
    read_only = TRUE)
  recipient <- scworkbench::wb_bundle_validation_report(relocated,
    level = "recipient")
  stopifnot(isTRUE(recipient$ok))
  for (name in c("deviation_z", "deviation_raw", "motif_summary")) {
    retained <- scworkbench::wb_art_load(fixture$bundle,
      artifact_ids[[name]])
    stopifnot(!is.null(retained))
    if (name != "motif_summary") stopifnot(
      identical(which(is.na(retained$value)),
        which(is.na(oracle$outputs[[name]]$value))),
      identical(scworkbench:::.wb_prov_payload_identity(
        fixture$bundle, retained)$content_sha256,
        scworkbench:::.wb_prov_payload_identity(fixture$bundle,
          result$run$outputs[[name]])$content_sha256))
  }
  pack_free_code <- paste(
    "args <- commandArgs(TRUE)",
    "devtools::load_all(args[[1L]], quiet = TRUE)",
    "stopifnot(!'scwbatac::chromvar' %in% wb_list_ops())",
    "bundle <- wb_open_bundle(args[[2L]], read_only = TRUE)",
    "stopifnot(nrow(wb_art_list(bundle)) == 3L)",
    "for (id in args[-(1:2)]) stopifnot(!is.null(wb_art_load(bundle,id)))",
    "stopifnot(isTRUE(wb_bundle_validation_report(args[[2L]],level='recipient')$ok))",
    "wb_close_bundle(bundle)", sep = "; ")
  pack_free_log <- tempfile("scwbatac-pack-free-", fileext = ".log")
  pack_free_status <- system2("Rscript", c("-e", shQuote(pack_free_code),
    shQuote(.atac_core_root), shQuote(relocated),
    shQuote(unname(artifact_ids))), stdout = pack_free_log,
    stderr = pack_free_log)
  if (!identical(pack_free_status, 0L)) stop(
    "Pack-free recipient inspection failed: ",
    paste(readLines(pack_free_log, warn = FALSE), collapse = "\n"))
  cat("SCWBATAC_LIVE_ORACLE_OK cells=", length(fixture$cells),
    " motifs=", length(unique(result$run$outputs$deviation_z$key)),
    "\n", sep = "")
  TRUE
}

test_chromvar_cached_oracle <- function() {
  fixture <- .atac_test_fixture()
  on.exit(.atac_test_close(fixture), add = TRUE)
  .atac_cached_match_set(fixture)
  manifest <- system.file("scworkbench", "operation-pack.json",
    package = "scwbATAC")
  scworkbench::wb_operation_pack_activate(manifest)
  set.seed(1701L)
  result <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    fixture$bundle, "scwbatac::chromvar", inputs = list(
      matrix = "ATAC_counts", motif_library = "jaspar2020_core_human",
      motif_match_set = "cached_jaspar", min_motif_peaks = 3L,
      background_iterations = 3L, background_bs = 10L,
      chunk_size = 35L, cell_set = "all"))
  oracle <- readRDS(file.path("tests", "fixtures",
    "old_core_chromvar_oracle.rds"))
  stopifnot(.atac_compare_oracle(result, oracle))
  promoted <- scworkbench::wb_session_promote_current(fixture$bundle,
    result$session, outputs = c("deviation_z", "deviation_raw",
      "motif_summary"), created_by = "phase7-cached-chromvar")
  stopifnot(nrow(promoted$promoted) == 3L,
    all(promoted$promoted$status == "ok"))
  artifact_ids <- stats::setNames(as.character(promoted$promoted$artifact_id),
    as.character(promoted$promoted$output))
  stopifnot(isTRUE(scworkbench::wb_finalize_recipient_bundle(
    fixture$bundle)$ready))
  scworkbench::wb_close_bundle(fixture$bundle)
  fixture$bundle <- scworkbench::wb_open_bundle(fixture$root,
    read_only = TRUE)
  stopifnot(isTRUE(scworkbench::wb_bundle_validation_report(
    fixture$root, level = "recipient")$ok))
  for (name in c("deviation_z", "deviation_raw")) {
    retained <- scworkbench::wb_art_load(fixture$bundle,
      artifact_ids[[name]], expected_kind = "cell_vector_family")
    stopifnot(identical(which(is.na(retained$value)),
      which(is.na(oracle$outputs[[name]]$value))),
      isTRUE(all.equal(as.numeric(retained$value),
        oracle$outputs[[name]]$value, tolerance = 1e-6)))
  }
  cat("SCWBATAC_CACHED_ORACLE_OK cells=", length(fixture$cells),
    " motifs=", length(unique(result$run$outputs$deviation_z$key)),
    "\n", sep = "")
  TRUE
}

if (identical(Sys.getenv("SCWBATAC_RUN_TESTS", ""), "1"))
  stopifnot(test_chromvar_live_oracle(), test_chromvar_cached_oracle())
