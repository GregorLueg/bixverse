# helpers ----------------------------------------------------------------------

## utils -----------------------------------------------------------------------

#' Helper to parse the verbosity
#'
#' @param input Boolean or integer to parse
#'
#' @returns The integer controlling the verbosity
#'
#' @keywords internal
parse_verbosity <- function(input) {
  # checks
  checkmate::qassert(input, c("B1", "I1[0, 2]"))

  as.integer(sum(input))
}

## cell ranger outputs ---------------------------------------------------------

#' Helper to generate cell ranger input parameters
#'
#' @param dir_data String. The directory with the Cell Ranger outputs
#'
#' @return A list based on [bixverse::params_sc_mtx_io()].
#'
#' @export
get_cell_ranger_params <- function(dir_data) {
  # checks
  checkmate::assertDirectory(dir_data)
  assertFileExists(dir_data, c("barcodes.tsv", "genes.tsv", "matrix.mtx"))

  res <- params_sc_mtx_io(
    path_mtx = path.expand(file.path(dir_data, "matrix.mtx")),
    path_obs = path.expand(file.path(dir_data, "barcodes.tsv")),
    path_var = path.expand(file.path(dir_data, "genes.tsv")),
    cells_as_rows = FALSE,
    has_hdr = FALSE
  )

  return(res)
}

## seurat assay to list --------------------------------------------------------

#' Transform Seurat raw counts into a List
#'
#' @param seurat_obj `Seurat` class. The class from which to extract the counts
#' from and transform to a list.
#'
#' @return A list with the following elements
#' \itemize{
#'   \item indptr - Index pointers of the sparse data.
#'   \item indices - Indices of the data.
#'   \item data - The underlying data.
#'   \item format - String that defines if the data is CSR or CSC.
#'   \item nrow - The number of rows.
#'   \item ncol - The number of columns.
#' }
get_seurat_counts_to_list <- function(seurat_obj) {
  # checks
  checkmate::assertClass(seurat_obj, "Seurat")

  assay <- seurat_obj@assays$RNA
  if (methods::.hasSlot(assay, "layers")) {
    raw_counts <- assay@layers$counts
  } else {
    raw_counts <- assay@counts
  }

  res <- list(
    indptr = raw_counts@p,
    indices = raw_counts@i,
    data = raw_counts@x,
    format = "csr",
    nrow = raw_counts@Dim[2],
    ncol = raw_counts@Dim[1]
  )

  return(res)
}

## meta cell matrices ----------------------------------------------------------

#' Create sparse dgRMatrix matrices for raw and norm counts
#'
#' @param meta_cell_data Named list. This contains the indptr, indices and data
#' for both raw counts and norm counts.
#'
#' @returns A list of two items
#' \itemize{
#'  \item raw - Sparse matrix representing the raw counts.
#'  \item norm - Sparse matrix representing the norm counts.
#' }
#'
#' @keywords internal
get_meta_cell_matrices <- function(meta_cell_data) {
  checkmate::assertList(meta_cell_data, names = "named")
  checkmate::assertNames(
    names(meta_cell_data),
    must.include = c(
      "indptr",
      "indices",
      "raw_counts",
      "norm_counts",
      "nrow",
      "ncol"
    )
  )

  dims <- as.integer(c(meta_cell_data$nrow, meta_cell_data$ncol))
  p <- as.integer(meta_cell_data$indptr)
  j <- as.integer(meta_cell_data$indices)

  list(
    raw = Matrix::sparseMatrix(
      p = p,
      j = j,
      x = as.numeric(meta_cell_data$raw_counts),
      dims = dims,
      repr = "R",
      index1 = FALSE
    ),
    norm = Matrix::sparseMatrix(
      p = p,
      j = j,
      x = as.numeric(meta_cell_data$norm_counts),
      dims = dims,
      repr = "R",
      index1 = FALSE
    )
  )
}

## qc --------------------------------------------------------------------------

### per cell outliers ----------------------------------------------------------

#' Use MAD outlier detection on per-cell QC metrics
#'
#' @param metric Numerical vector. The QC metric to check for.
#' @param threshold Numeric. How many MADs in either direction to consider for
#' outlier detection.
#' @param direction String. One of `c("twosided", "below", "above")`. Which
#' directionality to consider
#'
#' @returns A list with:
#' \itemize{
#'  \item outlier - Boolean vector indicating which cell is an outlier
#'  \item metrics - The applied thresholds.
#' }
per_cell_qc_outlier <- function(
  metric,
  threshold = 3,
  direction = c("twosided", "below", "above")
) {
  direction <- match.arg(direction)
  checkmate::qassert(metric, "N+")
  checkmate::qassert(threshold, "N1")

  outliers <- rs_mad_outlier(
    x = metric,
    threshold = threshold,
    direction = direction
  )

  med <- median(metric)
  metrics <- c(median = med, upper_threshold = NA, lower_threshold = NA)

  if (direction %in% c("twosided", "above")) {
    metrics["upper_threshold"] <- med + outliers$threshold
  }
  if (direction %in% c("twosided", "below")) {
    metrics["lower_threshold"] <- med - outliers$threshold
  }

  list(outlier = outliers$outlier, metrics = metrics)
}

### multiple metrics -----------------------------------------------------------

#' Run MAD outlier detection on per-cell QC metrics
#'
#' @param metrics Named list of numeric vectors. Each element is a QC metric
#' to check (e.g. `list(log10_lib_size = log10(lib_size), MT = mt_pct)`).
#' @param cells_to_keep Integer. Which cells were included in the analysis.
#' 0-indices for Rust.
#' @param directions Named character vector mapping metric names to direction.
#' One of `"twosided"`, `"below"`, `"above"`. Defaults to `"twosided"` for
#' all metrics if `NULL`.
#' @param threshold Numeric. Number of MADs to use for outlier detection.
#' @param groups Optional grouping variable. A string.
#'
#' @return An object of class `CellQc` containing:
#' \describe{
#'   \item{metrics}{The input metrics list.}
#'   \item{per_metric}{Named list of per-metric results from
#'     \code{\link{per_cell_qc_outlier}}.}
#'   \item{outlier_mat}{Logical matrix with one column per metric.}
#'   \item{combined}{Logical vector. `TRUE` if a cell is an outlier in any
#'     metric.}
#' }
#'
#' @export
run_cell_qc <- function(
  metrics,
  cells_to_keep,
  directions = NULL,
  threshold = 3,
  groups = NULL
) {
  checkmate::assertList(metrics, types = "numeric", names = "unique")
  checkmate::qassert(cells_to_keep, "I+")
  checkmate::qassert(threshold, "N1")

  n <- length(metrics[[1]])
  if (is.null(groups)) {
    groups <- rep("all", n)
  }
  checkmate::assertAtomic(groups, len = n, any.missing = FALSE)
  groups <- as.character(groups)

  if (is.null(directions)) {
    directions <- setNames(rep("twosided", length(metrics)), names(metrics))
  }

  group_levels <- unique(groups)

  results <- Map(
    function(x, dir) {
      outlier <- logical(n)
      thresholds <- list()
      for (g in group_levels) {
        idx <- which(groups == g)
        res <- per_cell_qc_outlier(
          metric = x[idx],
          threshold = threshold,
          direction = dir
        )
        outlier[idx] <- res$outlier
        thresholds[[g]] <- res$metrics
      }
      list(outlier = outlier, metrics = thresholds)
    },
    metrics,
    directions[names(metrics)]
  )

  outlier_mat <- do.call(cbind, lapply(results, `[[`, "outlier"))
  combined <- rowSums(outlier_mat) > 0

  structure(
    list(
      # transform to cell idx
      cell_idx = cells_to_keep + 1L,
      metrics = metrics,
      groups = groups,
      per_metric = results,
      outlier_mat = outlier_mat,
      combined = combined
    ),
    class = "CellQc"
  )
}

## knn -------------------------------------------------------------------------

### class ----------------------------------------------------------------------

#' Generate a new SingleCellNearestNeighbour from data
#'
#' @param data Numerical matrix. Samplex x features. The embedding matrix from
#' which to generate the kNN data.
#' @param neighbours_params List. Output of [bixverse::params_sc_neighbours()].
#' A list with the following items:
#' \itemize{
#'   \item full_snn - Boolean. Shall the full shared nearest neighbour graph
#'   be generated that generates edges between all cells instead of between
#'   only neighbours. (Not used in this function.)
#'   \item pruning - Numeric. Weights below this threshold will be set to 0 in
#'   the generation of the sNN graph. (Not used in this function.)
#'   \item snn_similarity - String. One of `c("rank", "jaccard")`. Defines how
#'   the weight from the SNN graph is calculated. For details, please see
#'   [bixverse::params_sc_neighbours()]. (Not used in this function.)
#'   \item knn - List of kNN parameters. See [bixverse::params_knn_defaults()]
#'   for available parameters and their defaults.
#' }
#' @param seed Integer. Random seed for reproducibility.
#' @param .validate_index Boolean. Shall an exhaustive search against a subset
#' of cells be run to validate the approximate nearest neighbour index.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns The `SingleCellNearestNeighbour` for downstream usage.
#'
#' @export
generate_sc_knn <- function(
  data,
  neighbours_params = params_sc_neighbours(),
  seed = 42L,
  .validate_index = FALSE,
  .verbose = TRUE
) {
  checkmate::assertMatrix(data, mode = "numeric")
  assertScNeighbours(neighbours_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.validate_index, "B1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  knn_data <- rs_sc_knn_w_dist(
    embd = data,
    knn_params = neighbours_params,
    verbose = parse_verbosity(.verbose),
    validate_index = .validate_index,
    seed = seed
  )

  knn <- new_sc_knn(knn_data = knn_data, used_cells = row.names(data))

  knn
}

### metrics --------------------------------------------------------------------

#' Calculate recall at k and distance ratio
#'
#' @description
#' Helper function to compare the results of two `SingleCellNearestNeighbour`
#' against each other. The first one can serve as a reference (ground truth)
#' and you can compare against the second one.
#'
#' @param ref_knn The reference `SingleCellNearestNeighbour`.
#' @param query_knn The query `SingleCellNearestNeighbour`.
#'
#' @returns A list with:
#' \itemize{
#'  \item matches - The intersecting indices between the reference and query
#'  kNN for each sample. In an ideal match up should be equal to k.
#'  \item distance_ratio - The distance ratio. Calculates
#'  `sum(dist_query) / sum(dist_ref)` per sample. Indicates how much worse the
#'  reference is.
#'  \item final_recall - The final recall across all samples.
#'  \item final_ratio - The final distance ratio across all samples.
#' }
#'
#' @export
calc_knn_metrics <- function(ref_knn, query_knn) {
  # checks
  checkmate::assertClass(ref_knn, "SingleCellNearestNeighbour")
  checkmate::assertClass(query_knn, "SingleCellNearestNeighbour")

  res <- rs_compare_knn(knn_data_a = ref_knn, knn_data_b = query_knn)

  res
}
