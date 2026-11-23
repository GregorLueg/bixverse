# shared single cell ATAC test fixtures -----------------------------------------

# The chr22 PBMC fixtures (Signac 1.17.1, Seurat 6.0.0, MACS3, AMULET) live in
# the bixverse-rs repository and are too large to ship with the package. Point
# BIXVERSE_ATAC_FIXTURES at `crates/bixverse-atac/tests/fixtures/pbmc3k_chr22`;
# the ATAC test files are skipped without it.

ATAC_FIX <- Sys.getenv(
  "BIXVERSE_ATAC_FIXTURES",
  path.expand(
    "~/repos/shared/bixverse-rs/crates/bixverse-atac/tests/fixtures/pbmc3k_chr22"
  )
)
ATAC_OK <- file.exists(file.path(ATAC_FIX, "fragments.tsv.gz"))

ATAC_FRAGS <- file.path(ATAC_FIX, "fragments.tsv.gz")
ATAC_GTF <- file.path(ATAC_FIX, "genes.gtf.gz")
ATAC_PEAKS <- file.path(ATAC_FIX, "peaks.bed")
ATAC_BARCODES <- if (ATAC_OK) readLines(file.path(ATAC_FIX, "cells.txt"))

#' Read a fixture table
#'
#' @param rel String. Path relative to the fixture directory.
#'
#' @returns A data.table.
atac_tsv <- function(rel) {
  data.table::fread(file.path(ATAC_FIX, rel), na.strings = c("NaN", "NA"))
}

#' Read a fixture matrix whose first column is a row label
#'
#' @param rel String. Path relative to the fixture directory.
#'
#' @returns A numeric matrix.
atac_dense <- function(rel) {
  as.matrix(atac_tsv(rel)[, -1])
}

#' Read the lines of a fixture file
#'
#' @param rel String. Path relative to the fixture directory.
#'
#' @returns Character vector.
atac_lines <- function(rel) {
  readLines(file.path(ATAC_FIX, rel))
}

#' Multi-modal experiment with synthetic RNA keyed by fixture barcodes
#'
#' @param barcodes Character vector. The cell ids.
#' @param peaks Boolean. Add the fixture peaks as the ATAC store.
#'
#' @returns The `SingleCellsMultiModal`.
atac_mm_experiment <- function(barcodes = ATAC_BARCODES, peaks = TRUE) {
  data <- generate_single_cell_test_data(
    params_sc_synthetic_data(
      n_cells = length(barcodes),
      n_genes = 60L,
      n_batches = 1L
    ),
    seed = 3L
  )
  data$obs$cell_id <- barcodes
  rownames(data$counts) <- barcodes
  dir <- tempfile("atac_mm")
  dir.create(dir)
  mm <- load_r_data(
    SingleCellsMultiModal(dir),
    counts = data$counts,
    obs = data$obs,
    var = data$var,
    sc_qc_param = params_sc_min_quality(
      min_unique_genes = 1L,
      min_lib_size = 1L,
      min_cells = 1L
    ),
    .verbose = FALSE
  )
  if (peaks) {
    mm <- add_atac_fragments_sc(
      mm,
      ATAC_FRAGS,
      peaks = ATAC_PEAKS,
      mode = "fragment",
      .verbose = FALSE
    )
  }

  mm
}

#' ATAC-only experiment over the fixture barcodes and peaks
#'
#' @param barcodes Character vector. The cell barcodes.
#' @param mode String. The count mode.
#'
#' @returns The `SingleCellsAtac`.
atac_only_experiment <- function(
  barcodes = ATAC_BARCODES,
  mode = "fragment"
) {
  dir <- tempfile("atac_only")
  dir.create(dir)
  load_atac_fragments(
    SingleCellsAtac(dir),
    ATAC_FRAGS,
    peaks = ATAC_PEAKS,
    mode = mode,
    barcodes = barcodes,
    .verbose = FALSE
  )
}
