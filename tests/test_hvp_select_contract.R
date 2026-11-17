.hvp_test_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
.hvp_core_root <- Sys.getenv("SCWBATAC_CORE_SOURCE", unset = "")
if (nzchar(.hvp_core_root)) {
  devtools::load_all(.hvp_core_root, quiet = TRUE)
  devtools::load_all(.hvp_test_root, quiet = TRUE)
} else {
  library(scworkbench)
  library(scwbATAC)
}

.hvp_fixture <- function() {
  root <- tempfile("scwbatac-hvp-")
  bundle <- scworkbench::wb_init_bundle(root)
  cells <- paste0("cell_", seq_len(8L))
  peaks <- paste0("chr1:", seq.int(1L, 501L, by = 100L), "-",
    seq.int(50L, 550L, by = 100L))
  values <- rbind(
    c(0, 0, 0, 0, 1, 1, 1, 1),
    c(1, 0, 1, 0, 1, 0, 1, 0),
    rep(1, 8),
    c(1, 0, 0, 0, 0, 0, 0, 0),
    c(1, 1, 0, 0, 0, 0, 0, 0),
    c(1, 1, 1, 1, 1, 1, 0, 0))
  dimnames(values) <- list(peaks, cells)
  scworkbench::wb_define_cells(bundle, cells)
  scworkbench::wb_attach_matrix(bundle, "ATAC_counts",
    Matrix::Matrix(values, sparse = TRUE), feature_type = "peak")
  list(root = root, bundle = bundle, cells = cells, peaks = peaks,
    values = values)
}

.hvp_close <- function(fixture) {
  try(scworkbench::wb_close_bundle(fixture$bundle), silent = TRUE)
  unlink(fixture$root, recursive = TRUE, force = TRUE)
}

test_hvp_select_depth_and_detection_contract <- function() {
  fixture <- .hvp_fixture()
  on.exit(.hvp_close(fixture), add = TRUE)
  manifest <- system.file("scworkbench", "operation-pack.json",
    package = "scwbATAC")
  try(scworkbench::wb_operation_pack_unload("scwbatac"), silent = TRUE)
  scworkbench::wb_operation_pack_activate(manifest)
  stopifnot(sum(scworkbench::wb_list_ops() ==
      "scwbatac::hvp_select") == 1L)

  adjusted <- scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), fixture$bundle,
    "scwbatac::hvp_select", inputs = list(
      matrix = "ATAC_counts", n_features = 3L,
      depth_adjustment = TRUE, min_detection_pct = 25,
      max_detection_pct = 75, chunk_size = 3L,
      cell_set = "all", min_cells = 2L))
  stopifnot(identical(sort(names(adjusted$run$outputs)),
    c("fit_summary", "selected_features")))
  table <- adjusted$run$outputs$selected_features$table
  expected_detection <- 100 * rowMeans(fixture$values > 0)
  names(expected_detection) <- fixture$peaks
  observed_detection <- table$detection_pct[
    match(fixture$peaks, table$feature_id)]
  binary <- fixture$values > 0
  log_depth <- log1p(colSums(binary))
  expected_residual_variance <- apply(binary, 1L, function(values) {
    fit <- stats::lm(as.numeric(values) ~ log_depth)
    sum(stats::residuals(fit)^2) / stats::df.residual(fit)
  })
  observed_residual_variance <- table$residual_variance[
    match(fixture$peaks, table$feature_id)]
  stopifnot(isTRUE(all.equal(observed_detection,
      unname(expected_detection), tolerance = 1e-12)),
    isTRUE(all.equal(observed_residual_variance,
      unname(expected_residual_variance), tolerance = 1e-12)),
    sum(table$selected) == 3L,
    all(table$passes_detection_filter ==
      (table$detection_pct >= 25 & table$detection_pct <= 75)),
    all(table$method == "log_depth_residual_variance"),
    all(table$depth_adjustment_status == "applied"),
    all(is.finite(table$score)), all(table$score >= 0),
    all(table$depth_r_squared >= 0 & table$depth_r_squared <= 1))
  summary <- adjusted$run$outputs$fit_summary$table
  stopifnot(isTRUE(summary$depth_adjustment),
    identical(summary$depth_metric, "detected_peaks"),
    summary$n_passing_detection_filter == 4L,
    summary$n_selected_peaks == 3L)

  raw <- scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    fixture$bundle, "scwbatac::hvp_select", inputs = list(
      matrix = "ATAC_counts", n_features = 6L,
      depth_adjustment = FALSE, min_detection_pct = 0,
      max_detection_pct = 100, chunk_size = 2L,
      cell_set = "all", min_cells = 2L))
  raw_table <- raw$run$outputs$selected_features$table
  stopifnot(all(raw_table$method == "raw_binary_variance"),
    isTRUE(all.equal(raw_table$score, raw_table$raw_variance,
      tolerance = 1e-12)), sum(raw_table$selected) == 6L)

  rejected <- try(scworkbench::wb_session_run_op(
    scworkbench::wb_session_new(), fixture$bundle,
    "scwbatac::hvp_select", inputs = list(
      matrix = "ATAC_counts", min_detection_pct = 80,
      max_detection_pct = 20, cell_set = "all", min_cells = 2L)),
    silent = TRUE)
  stopifnot(inherits(rejected, "try-error"),
    grepl("min_detection_pct", as.character(rejected), fixed = TRUE))

  promoted <- scworkbench::wb_session_promote_current(fixture$bundle,
    adjusted$session, outputs = c("selected_features", "fit_summary"),
    created_by = "hvp-select-contract")
  stopifnot(all(promoted$promoted$status == "ok"),
    setequal(promoted$promoted$output,
      c("selected_features", "fit_summary")))
  TRUE
}

if (identical(Sys.getenv("SCWBATAC_RUN_HVP", ""), "1")) {
  test_hvp_select_depth_and_detection_contract()
  cat("PASS test_hvp_select_depth_and_detection_contract\n")
}
