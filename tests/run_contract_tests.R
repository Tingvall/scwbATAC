# Source-worktree pack suite.  R CMD check sources this file but runs the
# expensive chromVAR method tests only when explicitly requested.
if (identical(Sys.getenv("SCWBATAC_RUN_SUITE", ""), "1")) {
  cases <- list(
    list(name = "live_and_cached", file = "tests/test_chromvar_pack_contract.R",
      flag = "SCWBATAC_RUN_TESTS", n = 2L),
    list(name = "negatives", file = "tests/test_chromvar_negative_contract.R",
      flag = "SCWBATAC_RUN_NEGATIVES", n = 1L),
    list(name = "stale_sources", file = "tests/test_chromvar_negative_contract.R",
      flag = "SCWBATAC_RUN_STALE", n = 1L),
    list(name = "seqlevels", file = "tests/test_chromvar_seqlevels_contract.R",
      flag = "SCWBATAC_RUN_SEQLEVELS", n = 1L)
  )
  discovered <- sum(vapply(cases, `[[`, integer(1), "n"))
  passed <- failed <- 0L
  for (case in cases) {
    output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
      case$file, env = paste0(case$flag, "=1"), stdout = TRUE,
      stderr = TRUE))
    status <- attr(output, "status")
    if (is.null(status)) status <- 0L
    cat(paste(output, collapse = "\n"), "\n")
    if (identical(as.integer(status), 0L)) {
      passed <- passed + case$n
      cat("PASS ", case$name, " tests=", case$n, "\n", sep = "")
    } else {
      failed <- failed + case$n
      cat("FAIL ", case$name, " tests=", case$n,
        " status=", status, "\n", sep = "")
    }
  }
  cat(sprintf("RESULT discovered=%d executed=%d passed=%d skipped=0 failed=%d\n",
    discovered, discovered, passed, failed))
  if (failed) quit(status = 1L)
}
