`%||%` <- function(x, y) if (is.null(x)) y else x

#' Register the scwbATAC operation pack
#'
#' Activation authenticates the installed ABI-v3 manifest and records the
#' fixed and matrix-selected dynamic R dependencies before dispatch.
#' @export
register_operation_pack <- function() {
  scworkbench::wb_operation_pack_register(
    system.file("scworkbench", "operation-pack.json", package = "scwbATAC"),
    list(`scwbatac::chromvar` = .atac_chromvar))
}
