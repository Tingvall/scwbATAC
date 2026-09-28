.atac_fail <- function(message) stop(paste0("scwbATAC: ", message),
  call. = FALSE)

.atac_scalar <- function(value) if (is.null(value) || !length(value)) "" else
  as.character(value)[1L]

.atac_peak_ranges <- function(snapshot, genome_pkg) {
  rows <- snapshot$peaks
  ids <- vapply(rows, function(x) .atac_scalar(x$feature_id), character(1))
  chr <- vapply(rows, function(x) .atac_scalar(x$chr), character(1))
  start <- vapply(rows, function(x) as.integer(x$start %||% NA_integer_),
    integer(1))
  end <- vapply(rows, function(x) as.integer(x$end %||% NA_integer_),
    integer(1))
  bad <- !nzchar(chr) | is.na(start) | is.na(end)
  if (any(bad)) {
    parsed <- lapply(trimws(ids[bad]), function(id) {
      match <- regmatches(id, regexec(
        "^([^:]+):(\\d+)-(\\d+)$", id, perl = TRUE))[[1L]]
      if (length(match) != 4L) match <- regmatches(id, regexec(
        "^([^-]+)-(\\d+)-(\\d+)$", id, perl = TRUE))[[1L]]
      if (length(match) != 4L) c(NA_character_, NA_character_,
        NA_character_, NA_character_) else match
    })
    chr[bad] <- vapply(parsed, `[[`, character(1), 2L)
    start[bad] <- suppressWarnings(as.integer(vapply(parsed, `[[`,
      character(1), 3L)))
    end[bad] <- suppressWarnings(as.integer(vapply(parsed, `[[`,
      character(1), 4L)))
  }
  if (anyNA(chr) || any(!nzchar(chr)) || anyNA(start) || anyNA(end) ||
      any(start < 1L) || any(end < start))
    .atac_fail("peak coordinates are absent, malformed or unparseable")
  ranges <- GenomicRanges::GRanges(seqnames = chr,
    ranges = IRanges::IRanges(start = start, end = end))
  names(ranges) <- ids
  build <- sub("^.*\\.([A-Za-z0-9]+)$", "\\1", genome_pkg)
  if (!identical(build, genome_pkg)) GenomeInfoDb::genome(ranges) <- build
  ranges
}

.atac_filter_seqlevels <- function(ranges, genome, params) {
  seqs <- as.character(GenomicRanges::seqnames(ranges))
  observed <- unique(seqs)
  keep_levels <- if (identical(params$seqlevels_mode, "standard"))
    intersect(observed, as.character(
      GenomeInfoDb::standardChromosomes(genome))) else observed
  if (!is.null(params$seqlevels_include)) keep_levels <- intersect(
    keep_levels, as.character(params$seqlevels_include))
  if (!is.null(params$seqlevels_exclude)) keep_levels <- setdiff(
    keep_levels, as.character(params$seqlevels_exclude))
  keep <- seqs %in% keep_levels
  if (!any(keep)) .atac_fail("no peaks remain after seqlevel filtering")
  list(keep = keep, kept = unique(seqs[keep]),
    dropped = unique(seqs[!keep]), input = length(seqs),
    retained = sum(keep))
}

.atac_motif_from_json <- function(value) {
  payload <- jsonlite::fromJSON(.atac_scalar(value),
    simplifyVector = TRUE)
  raw <- payload$profileMatrix
  profile <- if (is.matrix(raw)) apply(raw, c(1, 2), as.numeric) else
    do.call(rbind, lapply(raw, as.numeric))
  rownames(profile) <- as.character(payload$rownames %||% c("A", "C",
    "G", "T"))
  TFBSTools::PFMatrix(ID = .atac_scalar(payload$ID),
    name = .atac_scalar(payload$name),
    matrixClass = .atac_scalar(payload$matrixClass),
    strand = .atac_scalar(payload$strand), bg = unlist(payload$bg),
    tags = payload$tags %||% list(), profileMatrix = profile)
}

.atac_motifs <- function(snapshot, selected) {
  rows <- snapshot$motifs
  if (!is.null(selected)) rows <- Filter(function(x)
    .atac_scalar(x$motif_key) %in% as.character(selected), rows)
  if (!length(rows) || any(vapply(rows, function(x)
      !nzchar(.atac_scalar(x$payload_json)), logical(1))))
    .atac_fail("selected motif library contains no usable PWM payloads")
  values <- lapply(rows, function(x) .atac_motif_from_json(x$payload_json))
  names(values) <- vapply(rows, function(x) .atac_scalar(x$motif_key),
    character(1))
  do.call(TFBSTools::PFMatrixList, values)
}

.atac_cached_matches <- function(snapshot, feature_ids, motif_keys,
    library_snapshot, genomic_snapshot) {
  row <- snapshot$match_set
  if (is.null(row) || !identical(.atac_scalar(row$matrix_id),
      .atac_scalar(genomic_snapshot$matrix_id)) ||
      !identical(.atac_scalar(row$motif_library),
        .atac_scalar(library_snapshot$motif_library)) ||
      !identical(.atac_scalar(row$organism),
        .atac_scalar(genomic_snapshot$organism)) ||
      !identical(.atac_scalar(row$genome_build),
        .atac_scalar(genomic_snapshot$genome_build)) ||
      !identical(.atac_scalar(row$genome_provider),
        .atac_scalar(genomic_snapshot$genome_provider)))
    .atac_fail("cached motif-match set is not bound to the selected matrix, library and genome")
  membership <- snapshot$membership
  selected <- if (is.null(motif_keys)) unique(vapply(membership,
    function(x) .atac_scalar(x$motif_key), character(1))) else
    unique(as.character(motif_keys))
  if (!length(selected)) selected <- vapply(library_snapshot$motifs,
    function(x) .atac_scalar(x$motif_key), character(1))
  if (!length(selected)) .atac_fail("cached match set has no motifs")
  i <- match(vapply(membership, function(x) .atac_scalar(x$feature_id),
    character(1)), feature_ids)
  j <- match(vapply(membership, function(x) .atac_scalar(x$motif_key),
    character(1)), selected)
  keep <- !is.na(i) & !is.na(j)
  Matrix::sparseMatrix(i = i[keep], j = j[keep],
    x = rep(TRUE, sum(keep)),
    dims = c(length(feature_ids), length(selected)),
    dimnames = list(feature_ids, selected))
}

.atac_background <- function(object, params) {
  attempt <- function(bias = NULL) {
    args <- list(object = object,
      niterations = as.integer(params$background_iterations),
      w = as.numeric(params$background_w),
      bs = as.integer(params$background_bs))
    if (!is.null(bias)) args$bias <- as.numeric(bias)
    do.call(chromVAR::getBackgroundPeaks, args)
  }
  first <- tryCatch(attempt(), error = identity)
  if (!inherits(first, "error")) return(first)
  message <- conditionMessage(first)
  if (!grepl("chol.default\\(cov\\(norm_mat\\)\\)|leading minor",
      message, perl = TRUE)) stop(first)
  bias <- as.numeric(SummarizedExperiment::rowData(object)$bias)
  if (length(bias) != nrow(object)) stop(first)
  scale <- stats::sd(bias, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) scale <- 1
  rank <- rank(bias, ties.method = "first") / max(1L, length(bias))
  for (candidate in list(bias + rank * scale * 1e-4,
      bias + rank * scale * 1e-3, bias + rank * scale * 1e-2, rank)) {
    result <- tryCatch(attempt(candidate), error = identity)
    if (!inherits(result, "error")) return(result)
    first <- result
  }
  stop(first)
}

.atac_family <- function(value) {
  scworkbench::wb_prim_cell_vector_family(
    key = rep(rownames(value), times = ncol(value)),
    cell_id = rep(colnames(value), each = nrow(value)),
    value = as.numeric(value), value_type = "continuous")
}

.atac_chromvar <- function(readers, writers, params, context) {
  genomic <- readers$genomic_source()
  library <- readers$motif_library()
  cached <- readers$motif_match_set()
  source_cells <- as.character(context$typed_slice_reads$chunks$cells)
  source_features <- as.character(context$typed_slice_reads$chunks$features)
  if (!identical(source_cells,
      as.character(unlist(genomic$cell_id, use.names = FALSE))) ||
      !identical(source_features, vapply(genomic$peaks,
        function(x) .atac_scalar(x$feature_id), character(1))) ||
      !identical(.atac_scalar(genomic$feature_type), "peak") ||
      !nzchar(.atac_scalar(genomic$genome_provider)))
    .atac_fail("authenticated source is not one ordered peak matrix with a genome provider")
  provider <- .atac_scalar(genomic$genome_provider)
  genome <- tryCatch(getExportedValue(provider, provider), error = identity)
  if (inherits(genome, "error") || !inherits(genome, "BSgenome"))
    .atac_fail("core-attested provider did not export a BSgenome singleton")
  provider_builds <- unique(stats::na.omit(as.character(
    GenomeInfoDb::genome(genome))))
  if (length(provider_builds) != 1L || !identical(provider_builds[[1L]],
      .atac_scalar(genomic$genome_build)))
    .atac_fail("installed genome provider build differs from the authenticated matrix build")
  library_row <- library$library
  if (nzchar(.atac_scalar(library_row$organism)) &&
      !identical(.atac_scalar(library_row$organism),
        .atac_scalar(genomic$organism)))
    .atac_fail("motif library organism differs from matrix organism")
  if (nzchar(.atac_scalar(library_row$genome_build)) &&
      !identical(.atac_scalar(library_row$genome_build),
        .atac_scalar(genomic$genome_build)))
    .atac_fail("motif library build differs from matrix genome build")
  ranges <- .atac_peak_ranges(genomic, provider)
  filtered <- .atac_filter_seqlevels(ranges, genome, params)
  ranges <- ranges[filtered$keep]
  features <- source_features[filtered$keep]
  chunks <- split(source_cells,
    ceiling(seq_along(source_cells) / as.integer(params$chunk_size)))
  pieces <- lapply(chunks, function(cells) {
    piece <- readers$chunks(features = source_features,
      cells = as.character(cells))
    if (!identical(as.character(piece$feature_id), source_features) ||
        !identical(as.character(piece$cell_id), as.character(cells)))
      .atac_fail("typed slice changed authenticated feature or cell order")
    methods::as(Matrix::Matrix(piece$values[filtered$keep, , drop = FALSE],
      sparse = TRUE), "dgCMatrix")
  })
  counts <- methods::as(do.call(cbind, pieces), "dgCMatrix")
  dimnames(counts) <- list(features, source_cells)
  materialized_bytes <- as.numeric(utils::object.size(counts))
  nonzero <- Matrix::rowSums(counts) > 0
  if (!any(nonzero)) .atac_fail("selected peak matrix contains no nonzero peaks")
  counts <- counts[nonzero, , drop = FALSE]
  ranges <- ranges[nonzero]
  features <- features[nonzero]
  motifs <- .atac_motifs(library, params$motifs)
  motif_meta <- data.frame(motif_key = names(motifs),
    motif_id = vapply(motifs, TFBSTools::ID, character(1)),
    motif_name = vapply(motifs, TFBSTools::name, character(1)),
    stringsAsFactors = FALSE)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts), rowRanges = ranges,
    colData = S4Vectors::DataFrame(cell_id = source_cells,
      row.names = source_cells))
  se <- chromVAR::addGCBias(se, genome = genome)
  previous <- tryCatch(BiocParallel::registered(), error = identity)
  suppressWarnings(BiocParallel::register(BiocParallel::SerialParam(),
    default = TRUE))
  on.exit(if (is.list(previous) && length(previous)) try(
    suppressWarnings(BiocParallel::register(previous[[1L]], default = TRUE)),
    silent = TRUE), add = TRUE)
  match_source <- if (is.null(params$motif_match_set)) "computed" else
    "registry"
  matches <- if (identical(match_source, "computed"))
    motifmatchr::matchMotifs(pwms = motifs, subject = se,
      genome = genome, p.cutoff = as.numeric(params$motif_p_cutoff),
      w = as.integer(params$motif_w)) else
    .atac_cached_matches(cached, features, params$motifs, library, genomic)
  assay <- if (inherits(matches, "SummarizedExperiment"))
    SummarizedExperiment::assay(matches) else matches
  match_counts <- Matrix::colSums(assay)
  keep <- match_counts >= as.integer(params$min_motif_peaks)
  if (!is.null(params$max_motif_peaks)) keep <- keep &
    match_counts <= as.integer(params$max_motif_peaks)
  if (!any(keep)) .atac_fail("no motifs remain after peak-count filters")
  matches <- matches[, keep, drop = FALSE]
  assay <- if (inherits(matches, "SummarizedExperiment"))
    SummarizedExperiment::assay(matches) else matches
  match_counts <- as.integer(match_counts[keep])
  motif_meta <- motif_meta[match(colnames(assay), motif_meta$motif_key),
    , drop = FALSE]
  background <- .atac_background(se, params)
  deviations <- chromVAR::computeDeviations(se, matches,
    background_peaks = background)
  variability <- chromVAR::computeVariability(deviations,
    bootstrap_error = FALSE)
  z <- SummarizedExperiment::assay(deviations, "z")
  raw <- SummarizedExperiment::assay(deviations, "deviations")
  if (!identical(rownames(z), rownames(raw)) ||
      !identical(colnames(z), colnames(raw)))
    .atac_fail("deviation terminals disagree on their ordered axes")
  data <- as.data.frame(SummarizedExperiment::rowData(deviations))
  data$motif_key <- rownames(z)
  data$matched_peaks <- match_counts[match(data$motif_key,
    colnames(assay))]
  meta <- motif_meta[match(data$motif_key, motif_meta$motif_key),
    , drop = FALSE]
  variability$motif_key <- rownames(variability)
  summary <- data.frame(motif_key = as.character(data$motif_key),
    motif_id = as.character(meta$motif_id),
    motif_name = as.character(meta$motif_name),
    motif_library = .atac_scalar(library$motif_library),
    collection = .atac_scalar(library_row$collection),
    source_id = .atac_scalar(library_row$source_id),
    source_version = .atac_scalar(library_row$source_version),
    genome_pkg = provider,
    genome_build = .atac_scalar(genomic$genome_build),
    organism = .atac_scalar(genomic$organism),
    matched_peaks = as.integer(data$matched_peaks),
    fraction_matches = as.numeric(data$fractionMatches),
    fraction_background_overlap = as.numeric(
      data$fractionBackgroundOverlap),
    variability = as.numeric(variability$variability[match(
      data$motif_key, variability$motif_key)]),
    variability_p_value = as.numeric(variability$p_value[match(
      data$motif_key, variability$motif_key)]),
    variability_p_value_adj = as.numeric(variability$p_value_adj[match(
      data$motif_key, variability$motif_key)]),
    mean_abs_z = rowMeans(abs(z), na.rm = TRUE),
    stringsAsFactors = FALSE)
  writers$declare_runtime_cell_family_keys(rownames(z))
  writers$write_terminal_primitive("deviation_z", .atac_family(z))
  writers$write_terminal_primitive("deviation_raw", .atac_family(raw))
  writers$write_terminal_primitive("motif_summary",
    scworkbench::wb_prim_obs_table(summary,
      table_type = "chromvar_summary", entity_kind = "motif",
      meta = list(matrix_id = .atac_scalar(genomic$matrix_id),
        genome_pkg = provider,
        genome_provider_version = as.character(utils::packageVersion(provider)),
        motif_library = .atac_scalar(library$motif_library),
        motif_match_set = params$motif_match_set,
        motif_match_source = match_source,
        seqlevels_mode = params$seqlevels_mode,
        seqlevels_kept = filtered$kept,
        n_peaks_input = as.integer(filtered$input),
        n_peaks_retained = as.integer(filtered$retained),
        n_peaks_nonzero = as.integer(length(features)),
        n_motifs = as.integer(nrow(summary)),
        selected_materialized_matrix_bytes = materialized_bytes)))
  invisible(list(status = "ok"))
}
