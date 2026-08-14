# shared single cell test fixtures ---------------------------------------------

# Sourced from the top of the single cell test files with
# `source("helper_sc.R", local = TRUE)`. tinytest has no helper convention, but
# it `setwd()`s to the test directory and evaluates every top level expression
# of a file in one environment, so this resolves under `run_test_file()` and
# `test_package()` alike.
#
# No `expect_*` calls in here: tinytest only registers expectations returned at
# the top level of a test file.

## directories -----------------------------------------------------------------

#' Create a temporary directory for an on-disk single cell object
#'
#' @param name String. Directory name, appended to `parent`.
#' @param parent String. Where to create it. Defaults to the session `tempdir()`.
#'
#' @returns The path, guaranteed to exist.
#'
#' @keywords internal
sc_test_dir <- function(name, parent = tempdir()) {
  checkmate::qassert(name, "S1")
  checkmate::qassert(parent, "S1")

  path <- file.path(parent, name)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  stopifnot("Test directory does not exist" = dir.exists(path))

  path
}

#' Remove test directories
#'
#' @param ... Strings. The paths to remove.
#'
#' @returns `NULL`, invisibly.
#'
#' @keywords internal
sc_test_cleanup <- function(...) {
  paths <- c(...)
  checkmate::qassert(paths, "S*")

  unlink(paths, recursive = TRUE, force = TRUE)

  invisible(NULL)
}

## synthetic data --------------------------------------------------------------

#' Standard synthetic single cell fixture
#'
#' @description
#' The 1000 cell x 100 gene default of [generate_single_cell_test_data()] plus
#' the QC thresholds every single cell test file uses and the gene/cell indices
#' that pass them.
#'
#' @param min_lib_size Integer. Minimum library size per cell.
#' @param min_genes_exp Integer. Minimum number of expressed genes per cell.
#' @param min_cells_exp Integer. Minimum number of cells expressing a gene.
#' @param hvg_to_keep Integer. Number of HVGs downstream tests ask for.
#' @param no_pcs Integer. Number of principal components downstream tests ask
#' for.
#' @param syn_data_params List. See [params_sc_synthetic_data()].
#' @param seed Integer. Seed handed to the generator.
#'
#' @returns A list with the thresholds, `counts` / `obs` / `var` from the
#' generator, and the 1-indexed `genes_pass` / `cells_pass`.
#'
#' @keywords internal
sc_test_fixture <- function(
  min_lib_size = 300L,
  min_genes_exp = 45L,
  min_cells_exp = 500L,
  hvg_to_keep = 30L,
  no_pcs = 10L,
  syn_data_params = params_sc_synthetic_data(),
  seed = 42L
) {
  checkmate::qassert(min_lib_size, "I1")
  checkmate::qassert(min_genes_exp, "I1")
  checkmate::qassert(min_cells_exp, "I1")
  checkmate::qassert(hvg_to_keep, "I1")
  checkmate::qassert(no_pcs, "I1")
  checkmate::assertList(syn_data_params, names = "named")
  checkmate::qassert(seed, "I1")

  data <- generate_single_cell_test_data(
    syn_data_params = syn_data_params,
    seed = seed
  )

  genes_pass <- which(
    Matrix::colSums(data$counts != 0) >= min_cells_exp
  )

  cells_pass <- which(
    (Matrix::rowSums(data$counts[, genes_pass]) >= min_lib_size) &
      (Matrix::rowSums(data$counts[, genes_pass] != 0) >= min_genes_exp)
  )

  list(
    min_lib_size = min_lib_size,
    min_genes_exp = min_genes_exp,
    min_cells_exp = min_cells_exp,
    hvg_to_keep = hvg_to_keep,
    no_pcs = no_pcs,
    counts = data$counts,
    obs = data$obs,
    var = data$var,
    genes_pass = genes_pass,
    cells_pass = cells_pass
  )
}

#' Multi-batch synthetic fixture
#'
#' @description
#' The four cell type, three batch data the batch correction tests run on.
#'
#' @param batch_effect_strength String. One of `c("strong", "medium", "weak")`.
#' @param n_cells Integer. Number of cells.
#' @param ... Passed on to [sc_test_fixture()].
#'
#' @returns As [sc_test_fixture()].
#'
#' @keywords internal
sc_batch_fixture <- function(
  batch_effect_strength = c("strong", "medium", "weak"),
  n_cells = 900L,
  ...
) {
  batch_effect_strength <- match.arg(batch_effect_strength)
  checkmate::qassert(n_cells, "I1")

  cell_markers <- list(
    cell_type_1 = list(marker_genes = 0:8L),
    cell_type_2 = list(marker_genes = 9:19L),
    cell_type_3 = list(marker_genes = 20:29L),
    cell_type_4 = list(marker_genes = 30:44L)
  )

  sc_test_fixture(
    syn_data_params = params_sc_synthetic_data(
      n_cells = n_cells,
      marker_genes = cell_markers,
      n_batches = 3L,
      batch_effect_strength = batch_effect_strength
    ),
    ...
  )
}

## objects ---------------------------------------------------------------------

#' QC parameters matching a fixture's thresholds
#'
#' @param fixture List. Output of [sc_test_fixture()].
#' @param ... Overrides passed to [params_sc_min_quality()].
#'
#' @returns The parameter list.
#'
#' @keywords internal
sc_test_qc_params <- function(fixture, ...) {
  checkmate::assertList(fixture, names = "named")

  defaults <- list(
    min_unique_genes = fixture$min_genes_exp,
    min_lib_size = fixture$min_lib_size,
    min_cells = fixture$min_cells_exp
  )
  overrides <- list(...)

  do.call(
    params_sc_min_quality,
    utils::modifyList(defaults, overrides)
  )
}

#' Build a `SingleCells` object off a fixture
#'
#' @param dir String. Directory for the object, see [sc_test_dir()].
#' @param fixture List. Output of [sc_test_fixture()].
#' @param obs data.table. Observation table, defaults to the fixture's. Pass a
#' modified one to add grouping columns.
#' @param counts dgRMatrix. Counts, defaults to the fixture's.
#' @param sc_qc_param List. Defaults to [sc_test_qc_params()] off the fixture.
#' @param streaming Integer. One of `0L`, `1L`, `2L`. The suite runs the
#' in-memory path by default, the package default is `1L`.
#' @param ... Passed on to [load_r_data()].
#'
#' @returns The loaded `SingleCells` object.
#'
#' @keywords internal
sc_test_object <- function(
  dir,
  fixture,
  obs = fixture$obs,
  counts = fixture$counts,
  sc_qc_param = sc_test_qc_params(fixture),
  streaming = 0L,
  ...
) {
  checkmate::assertDirectoryExists(dir)
  checkmate::assertList(fixture, names = "named")
  checkmate::assertDataFrame(obs)
  checkmate::assertClass(counts, "dgRMatrix")
  checkmate::assertList(sc_qc_param, names = "named")
  checkmate::assertChoice(streaming, c(0L, 1L, 2L))

  load_r_data(
    object = SingleCells(dir_data = dir),
    counts = counts,
    obs = obs,
    var = fixture$var,
    sc_qc_param = sc_qc_param,
    streaming = streaming,
    .verbose = FALSE,
    ...
  )
}

#' Run the standard HVG -> PCA -> kNN chain
#'
#' @param object `SingleCells` or `SingleCellsSubset`.
#' @param fixture List. Output of [sc_test_fixture()].
#' @param k Integer. Number of neighbours.
#'
#' @returns The object with HVGs, PCA factors and a kNN/sNN graph cached.
#'
#' @keywords internal
sc_test_prepped <- function(object, fixture, k = 15L) {
  checkmate::assertList(fixture, names = "named")
  checkmate::qassert(k, "I1")

  object <- find_hvg_sc(
    object = object,
    hvg_no = fixture$hvg_to_keep,
    .verbose = FALSE
  )
  object <- calculate_pca_sc(
    object = object,
    no_pcs = fixture$no_pcs,
    .verbose = FALSE
  )

  find_neighbours_sc(
    object,
    neighbours_params = params_sc_neighbours(knn = list(k = k)),
    .verbose = FALSE
  )
}
