if (!exists(".atac_test_fixture", mode = "function"))
  source(if (file.exists("tests/test_chromvar_pack_contract.R"))
    "tests/test_chromvar_pack_contract.R" else
    "test_chromvar_pack_contract.R", local = FALSE)

.atac_negative_run <- function(fixture, extra = list()) {
  args <- modifyList(list(matrix = "ATAC_counts",
    motif_library = "jaspar2020_core_human",
    min_motif_peaks = 3L, background_iterations = 3L,
    background_bs = 10L, chunk_size = 35L, cell_set = "all"), extra)
  scworkbench::wb_session_run_op(scworkbench::wb_session_new(),
    fixture$bundle, "scwbatac::chromvar", inputs = args)
}

.atac_negative_case <- function(label, mutate = NULL, extra = list(),
    message = NULL) {
  fixture <- .atac_test_fixture()
  on.exit(.atac_test_close(fixture), add = TRUE)
  if (is.function(mutate)) mutate(fixture)
  before <- c(artifacts = nrow(scworkbench::wb_art_list(fixture$bundle)),
    commits = nrow(DBI::dbGetQuery(fixture$bundle$con,
      "SELECT * FROM prov_commits;")))
  failure <- try(.atac_negative_run(fixture, extra), silent = TRUE)
  after <- c(artifacts = nrow(scworkbench::wb_art_list(fixture$bundle)),
    commits = nrow(DBI::dbGetQuery(fixture$bundle$con,
      "SELECT * FROM prov_commits;")))
  if (!inherits(failure, "try-error") || !identical(before, after) ||
      (!is.null(message) && !grepl(message, as.character(failure),
        fixed = TRUE))) stop("negative case failed: ", label,
      " result=", paste(as.character(failure), collapse = "\n"))
  invisible(TRUE)
}

test_chromvar_negative_contract <- function() {
  manifest <- system.file("scworkbench", "operation-pack.json",
    package = "scwbATAC")
  scworkbench::wb_operation_pack_activate(manifest)
  set_provider <- function(provider) function(fixture)
    scworkbench::wb_matrix_set_genome(fixture$bundle, "ATAC_counts",
      organism = "human", genome_build = "hg19",
      genome_provider = provider)
  .atac_negative_case("absent provider", set_provider(NULL),
    message = "genome provider is absent or malformed")
  .atac_negative_case("malformed provider", set_provider("bad/provider"),
    message = "genome provider is absent or malformed")
  .atac_negative_case("missing provider package",
    set_provider("BSgenome.NotInstalled.Unknown"),
    message = "is unavailable")
  if (requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE))
    .atac_negative_case("installed incompatible provider",
      set_provider("BSgenome.Mmusculus.UCSC.mm10"),
      message = "installed genome provider build differs")
  .atac_negative_case("wrong motif-library build", function(fixture)
    DBI::dbExecute(fixture$bundle$con, "UPDATE motif_libraries
      SET genome_build='mm10' WHERE motif_library='jaspar2020_core_human';"),
    message = "motif library build differs")
  .atac_negative_case("unknown motif", extra = list(motifs = "absent"),
    message = "no usable PWM payloads")
  .atac_negative_case("duplicate motif selector", extra = list(
    motifs = c("MA0139.1_CTCF", "MA0139.1_CTCF")))
  .atac_negative_case("wrong cached-match build", function(fixture) {
    motif <- scworkbench::wb_motifs(fixture$bundle,
      motif_library = "jaspar2020_core_human")$motif_key[[1L]]
    scworkbench::wb_motif_match_replace(fixture$bundle,
      motif_match_set = "bad_cached", matrix_id = "ATAC_counts",
      motif_library = "jaspar2020_core_human",
      matches = data.frame(motif_key = motif,
        feature_id = fixture$features[[1L]]),
      organism = "human", genome_build = "mm10",
      genome_provider = "BSgenome.Hsapiens.UCSC.hg19")
  }, extra = list(motif_match_set = "bad_cached"),
    message = "not bound to the selected matrix, library and genome")
  .atac_negative_case("wrong cached-match matrix", function(fixture) {
    other <- Matrix::Matrix(matrix(1L,
      nrow = length(fixture$features), ncol = length(fixture$cells),
      dimnames = list(fixture$features, fixture$cells)), sparse = TRUE)
    scworkbench::wb_attach_matrix(fixture$bundle, "other_peaks", other,
      feature_type = "peak")
    scworkbench::wb_matrix_set_genome(fixture$bundle, "other_peaks",
      organism = "human", genome_build = "hg19",
      genome_provider = "BSgenome.Hsapiens.UCSC.hg19")
    motif <- scworkbench::wb_motifs(fixture$bundle,
      motif_library = "jaspar2020_core_human")$motif_key[[1L]]
    scworkbench::wb_motif_match_replace(fixture$bundle,
      motif_match_set = "foreign_matrix", matrix_id = "other_peaks",
      motif_library = "jaspar2020_core_human",
      matches = data.frame(motif_key = motif,
        feature_id = fixture$features[[1L]]),
      organism = "human", genome_build = "hg19",
      genome_provider = "BSgenome.Hsapiens.UCSC.hg19")
  }, extra = list(motif_match_set = "foreign_matrix"),
    message = "not bound to the selected matrix, library and genome")
  cat("SCWBATAC_NEGATIVES_OK cases=9\n")
  TRUE
}

test_chromvar_stale_source_contract <- function() {
  manifest <- system.file("scworkbench", "operation-pack.json",
    package = "scwbATAC")
  scworkbench::wb_operation_pack_activate(manifest)
  mutations <- list(
    matrix = function(fixture) scworkbench::wb_matrix_set_genome(
      fixture$bundle, "ATAC_counts", organism = "human",
      genome_build = "hg19",
      genome_provider = "BSgenome.Hsapiens.UCSC.hg38"),
    peak = function(fixture) DBI::dbExecute(fixture$bundle$con,
      "UPDATE features SET chr='chr9' WHERE matrix_id='ATAC_counts'
        AND feature_id=?;", params = list(fixture$features[[1L]])),
    motif = function(fixture) DBI::dbExecute(fixture$bundle$con,
      "UPDATE motifs SET motif_name='changed-after-stage'
        WHERE motif_library='jaspar2020_core_human' AND motif_key=?;",
      params = list(scworkbench::wb_motifs(fixture$bundle,
        motif_library = "jaspar2020_core_human")$motif_key[[1L]]))
  )
  for (label in names(mutations)) {
    fixture <- .atac_test_fixture()
    set.seed(1701L)
    staged <- .atac_negative_run(fixture)
    mutations[[label]](fixture)
    before <- c(artifacts = nrow(scworkbench::wb_art_list(fixture$bundle)),
      commits = nrow(DBI::dbGetQuery(fixture$bundle$con,
        "SELECT * FROM prov_commits;")))
    before_files <- sort(list.files(fixture$root, all.files = TRUE,
      recursive = TRUE, no.. = TRUE))
    failure <- try(scworkbench::wb_session_promote_current(fixture$bundle,
      staged$session, outputs = c("deviation_z", "deviation_raw",
        "motif_summary")), silent = TRUE)
    after <- c(artifacts = nrow(scworkbench::wb_art_list(fixture$bundle)),
      commits = nrow(DBI::dbGetQuery(fixture$bundle$con,
        "SELECT * FROM prov_commits;")))
    after_files <- sort(list.files(fixture$root, all.files = TRUE,
      recursive = TRUE, no.. = TRUE))
    if (!inherits(failure, "try-error") || !identical(before, after) ||
        !identical(before_files, after_files))
      stop("stale source mutation was not rejected atomically: ", label,
        " result=", paste(as.character(failure), collapse = "\n"))
    .atac_test_close(fixture)
  }
  fixture <- .atac_test_fixture()
  .atac_cached_match_set(fixture)
  set.seed(1701L)
  staged <- .atac_negative_run(fixture,
    extra = list(motif_match_set = "cached_jaspar"))
  motif <- scworkbench::wb_motifs(fixture$bundle,
    motif_library = "jaspar2020_core_human")$motif_key[[1L]]
  scworkbench::wb_motif_match_replace(fixture$bundle,
    motif_match_set = "cached_jaspar", matrix_id = "ATAC_counts",
    motif_library = "jaspar2020_core_human",
    matches = data.frame(motif_key = motif,
      feature_id = fixture$features[[1L]]),
    organism = "human", genome_build = "hg19",
    genome_provider = "BSgenome.Hsapiens.UCSC.hg19",
    match_params = list(p_cutoff = 5e-5, w = 7L))
  before <- c(artifacts = nrow(scworkbench::wb_art_list(fixture$bundle)),
    commits = nrow(DBI::dbGetQuery(fixture$bundle$con,
      "SELECT * FROM prov_commits;")))
  before_files <- sort(list.files(fixture$root, all.files = TRUE,
    recursive = TRUE, no.. = TRUE))
  failure <- try(scworkbench::wb_session_promote_current(fixture$bundle,
    staged$session, outputs = c("deviation_z", "deviation_raw",
      "motif_summary")), silent = TRUE)
  after <- c(artifacts = nrow(scworkbench::wb_art_list(fixture$bundle)),
    commits = nrow(DBI::dbGetQuery(fixture$bundle$con,
      "SELECT * FROM prov_commits;")))
  after_files <- sort(list.files(fixture$root, all.files = TRUE,
    recursive = TRUE, no.. = TRUE))
  if (!inherits(failure, "try-error") || !identical(before, after) ||
      !identical(before_files, after_files))
    stop("stale cached-match membership was not rejected atomically")
  .atac_test_close(fixture)
  cat("SCWBATAC_STALE_SOURCES_OK cases=4\n")
  TRUE
}

if (identical(Sys.getenv("SCWBATAC_RUN_NEGATIVES", ""), "1"))
  test_chromvar_negative_contract()
if (identical(Sys.getenv("SCWBATAC_RUN_STALE", ""), "1"))
  test_chromvar_stale_source_contract()
