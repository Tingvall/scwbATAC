.atac_hvp_select <- function(readers, writers, params, context) {
  capability <- context$typed_slice_reads$matrix_chunks
  features <- as.character(capability$features)
  cells <- as.character(capability$cells)
  n_cells <- length(cells)
  if (!length(features) || n_cells < 2L)
    .atac_fail("hvp_select requires at least one peak and two cells")

  annotations <- readers$feature_annotations()
  rows <- annotations$rows
  annotation_ids <- vapply(rows, function(row)
    .atac_scalar(row$feature_id), character(1))
  if (!identical(features, annotation_ids))
    .atac_fail("hvp_select feature annotations are not matrix aligned")
  feature_types <- vapply(rows, function(row)
    .atac_scalar(row$values$feature_type), character(1))
  if (any(!nzchar(feature_types)) || any(feature_types != "peak"))
    .atac_fail("hvp_select requires a matrix whose features are all annotated as peaks")

  minimum <- as.numeric(params$min_detection_pct %||% 0)
  maximum <- as.numeric(params$max_detection_pct %||% 100)
  if (!is.finite(minimum) || !is.finite(maximum) || minimum < 0 ||
      maximum > 100 || minimum > maximum)
    .atac_fail("hvp_select requires 0 <= min_detection_pct <= max_detection_pct <= 100")
  adjust_depth <- isTRUE(params$depth_adjustment %||% TRUE)
  if (adjust_depth && n_cells < 3L)
    .atac_fail("hvp_select depth adjustment requires at least three cells")

  detected <- numeric(length(features))
  accessibility_depth <- numeric(n_cells)
  if (adjust_depth) {
    sum_y_log_depth <- numeric(length(features))
    sum_log_depth <- 0
    sum_log_depth_squared <- 0
  }
  offset <- 0L
  chunks <- scworkbench::wb_partition_axis(cells,
    as.integer(params$chunk_size %||% 256L))
  for (chunk_cells in chunks) {
    primitive <- readers$matrix_chunks(features = features,
      cells = chunk_cells)
    if (!identical(as.character(primitive$feature_id), features) ||
        !identical(as.character(primitive$cell_id), chunk_cells))
      .atac_fail("hvp_select typed slice changed the authenticated matrix axes")
    binary <- primitive$values > 0
    chunk_depth <- as.numeric(Matrix::colSums(binary))
    index <- seq.int(offset + 1L, offset + length(chunk_cells))
    accessibility_depth[index] <- chunk_depth
    offset <- offset + length(chunk_cells)
    detected <- detected + as.numeric(Matrix::rowSums(binary))
    if (adjust_depth) {
      log_depth <- log1p(chunk_depth)
      sparse_log_depth <- Matrix::Matrix(log_depth, ncol = 1L,
        sparse = TRUE)
      sum_y_log_depth <- sum_y_log_depth +
        as.numeric(binary %*% sparse_log_depth)
      sum_log_depth <- sum_log_depth + sum(log_depth)
      sum_log_depth_squared <- sum_log_depth_squared + sum(log_depth^2)
    }
  }

  mean_accessibility <- detected / n_cells
  detection_pct <- 100 * mean_accessibility
  total_ss <- pmax(0, detected - detected^2 / n_cells)
  raw_variance <- total_ss / (n_cells - 1L)
  if (adjust_depth) {
    centered_depth_ss <- pmax(0,
      sum_log_depth_squared - sum_log_depth^2 / n_cells)
    centered_cross_product <- sum_y_log_depth -
      detected * sum_log_depth / n_cells
    depth_varies <- is.finite(centered_depth_ss) &&
      centered_depth_ss > sqrt(.Machine$double.eps)
    if (depth_varies) {
      depth_slope <- centered_cross_product / centered_depth_ss
      explained_ss <- pmin(total_ss, pmax(0,
        centered_cross_product^2 / centered_depth_ss))
      residual_ss <- pmax(0, total_ss - explained_ss)
      residual_variance <- residual_ss / (n_cells - 2L)
      depth_r_squared <- ifelse(total_ss > 0, explained_ss / total_ss, 0)
      adjustment_status <- "applied"
    } else {
      depth_slope <- rep(0, length(features))
      residual_variance <- raw_variance
      depth_r_squared <- rep(0, length(features))
      adjustment_status <- "constant_depth"
    }
  } else {
    depth_slope <- rep(0, length(features))
    residual_variance <- raw_variance
    depth_r_squared <- rep(0, length(features))
    adjustment_status <- "not_requested"
  }
  score <- if (adjust_depth) residual_variance else raw_variance
  method <- if (adjust_depth) "log_depth_residual_variance" else
    "raw_binary_variance"
  passes <- is.finite(score) & detection_pct >= minimum &
    detection_pct <= maximum
  if (!any(passes))
    .atac_fail("hvp_select detection filters removed every peak")

  order_index <- order(!passes, -score, -raw_variance, features,
    na.last = TRUE, method = "radix")
  rank <- integer(length(features))
  rank[order_index] <- seq_along(order_index)
  passing_order <- order_index[passes[order_index]]
  selected_index <- passing_order[seq_len(min(length(passing_order),
    as.integer(params$n_features %||% 25000L)))]
  selected <- rep(FALSE, length(features))
  selected[selected_index] <- TRUE
  matrix_id <- as.character(context$input_ids$matrix)[1]

  table <- data.frame(
    feature_id = features,
    rank = as.integer(rank),
    mean_accessibility = as.numeric(mean_accessibility),
    detected_cells = as.integer(detected),
    detection_pct = as.numeric(detection_pct),
    raw_variance = as.numeric(raw_variance),
    residual_variance = as.numeric(residual_variance),
    depth_slope = as.numeric(depth_slope),
    depth_r_squared = as.numeric(depth_r_squared),
    score = as.numeric(score),
    passes_detection_filter = passes,
    selected = selected,
    method = method,
    depth_adjustment_status = adjustment_status,
    matrix_id = matrix_id,
    stringsAsFactors = FALSE)
  table <- table[order(table$rank), , drop = FALSE]
  rownames(table) <- NULL

  summary <- data.frame(
    matrix_id = matrix_id,
    method = method,
    depth_adjustment = adjust_depth,
    depth_metric = "detected_peaks",
    depth_transform = "log1p",
    depth_adjustment_status = adjustment_status,
    n_cells = as.integer(n_cells),
    n_input_peaks = as.integer(length(features)),
    n_passing_detection_filter = as.integer(sum(passes)),
    n_selected_peaks = as.integer(sum(selected)),
    min_detection_pct = minimum,
    max_detection_pct = maximum,
    min_cell_depth = as.numeric(min(accessibility_depth)),
    median_cell_depth = as.numeric(stats::median(accessibility_depth)),
    max_cell_depth = as.numeric(max(accessibility_depth)),
    chunk_size = as.integer(params$chunk_size %||% 256L),
    stringsAsFactors = FALSE)
  meta <- list(matrix_id = matrix_id, selection_kind =
      "highly_variable_peaks", method = method,
    depth_adjustment = adjust_depth, depth_metric = "detected_peaks",
    min_detection_pct = minimum, max_detection_pct = maximum,
    n_selected_features = as.integer(sum(selected)))
  writers$write_terminal_primitive("selected_features",
    scworkbench::wb_prim_feature_table(table,
      table_type = "feature_selection_table", meta = meta))
  writers$write_terminal_primitive("fit_summary",
    scworkbench::wb_prim_obs_table(summary,
      table_type = "hvp_select_summary", entity_kind = "matrix",
      meta = meta))
  invisible(list(status = "ok"))
}
