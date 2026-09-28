.atac_seqlevel_core <- Sys.getenv("SCWBATAC_CORE_SOURCE", unset = "")
if (nzchar(.atac_seqlevel_core)) {
  devtools::load_all(.atac_seqlevel_core, quiet = TRUE)
  devtools::load_all(normalizePath(".", mustWork = TRUE), quiet = TRUE)
}

test_chromvar_seqlevels_contract <- function() {
  if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE))
    stop("Phase 7 seqlevel contract requires the installed mm10 provider")
  genome <- getExportedValue("BSgenome.Mmusculus.UCSC.mm10",
    "BSgenome.Mmusculus.UCSC.mm10")
  ranges <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr2", "GL456233.1", "chrM"),
    ranges = IRanges::IRanges(start = c(1, 10, 20, 30), width = 50))
  standard <- scwbATAC:::.atac_filter_seqlevels(ranges, genome,
    list(seqlevels_mode = "standard", seqlevels_include = NULL,
      seqlevels_exclude = NULL))
  stopifnot(identical(standard$keep, c(TRUE, TRUE, FALSE, TRUE)),
    identical(standard$kept, c("chr1", "chr2", "chrM")),
    identical(standard$dropped, "GL456233.1"))
  all_levels <- scwbATAC:::.atac_filter_seqlevels(ranges, genome,
    list(seqlevels_mode = "all", seqlevels_include = NULL,
      seqlevels_exclude = NULL))
  stopifnot(all(all_levels$keep))
  excluded <- scwbATAC:::.atac_filter_seqlevels(ranges, genome,
    list(seqlevels_mode = "standard", seqlevels_include = NULL,
      seqlevels_exclude = "chrM"))
  stopifnot(identical(excluded$keep, c(TRUE, TRUE, FALSE, FALSE)))
  TRUE
}

if (identical(Sys.getenv("SCWBATAC_RUN_SEQLEVELS", ""), "1"))
  test_chromvar_seqlevels_contract()
