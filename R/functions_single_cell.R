# helpers ----------------------------------------------------------------------

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

  # get the counts
  raw_counts <- seurat_obj@assays$RNA@layers$counts

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
  # checks
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

  list(
    raw = new(
      "dgRMatrix",
      p = as.integer(meta_cell_data$indptr),
      j = as.integer(meta_cell_data$indices),
      x = as.numeric(meta_cell_data$raw_counts),
      Dim = as.integer(c(meta_cell_data$nrow, meta_cell_data$ncol))
    ),
    norm = new(
      "dgRMatrix",
      p = as.integer(meta_cell_data$indptr),
      j = as.integer(meta_cell_data$indices),
      x = as.numeric(meta_cell_data$norm_counts),
      Dim = as.integer(c(meta_cell_data$nrow, meta_cell_data$ncol))
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
#'   to check (e.g. `list(log10_lib_size = log10(lib_size), MT = mt_pct)`).
#' @param directions Named character vector mapping metric names to direction.
#'   One of `"twosided"`, `"below"`, `"above"`. Defaults to `"twosided"` for
#'   all metrics if `NULL`.
#' @param threshold Numeric. Number of MADs to use for outlier detection.
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
run_cell_qc <- function(metrics, directions = NULL, threshold = 3) {
  checkmate::assertList(metrics, types = "numeric", names = "unique")
  checkmate::qassert(threshold, "N1")

  if (is.null(directions)) {
    directions <- setNames(rep("twosided", length(metrics)), names(metrics))
  }

  results <- Map(
    function(x, dir) {
      per_cell_qc_outlier(metric = x, threshold = threshold, direction = dir)
    },
    metrics,
    directions[names(metrics)]
  )

  outlier_mat <- do.call(cbind, lapply(results, `[[`, "outlier"))
  combined <- rowSums(outlier_mat) > 0

  structure(
    list(
      metrics = metrics,
      per_metric = results,
      outlier_mat = outlier_mat,
      combined = combined
    ),
    class = "CellQc"
  )
}

#### R primitives --------------------------------------------------------------

#' Print a CellQc object
#'
#' @param x A `CellQc` object.
#' @param ... Ignored.
#'
#' @return Invisible `x`.
#'
#' @export
#'
#' @keywords internal
print.CellQc <- function(x, ...) {
  n_cells <- length(x$combined)
  n_outliers <- sum(x$combined)
  metric_names <- colnames(x$outlier_mat)

  cat(sprintf(
    "CellQc: %d cells, %d outliers (%.1f%%)\n",
    n_cells,
    n_outliers,
    100 * n_outliers / n_cells
  ))
  cat("Metrics:\n")

  for (nm in metric_names) {
    n <- sum(x$outlier_mat[, nm])
    thresholds <- x$per_metric[[nm]]$metrics
    parts <- character()
    if (!is.na(thresholds["lower_threshold"])) {
      parts <- c(parts, sprintf("lower = %.2f", thresholds["lower_threshold"]))
    }
    if (!is.na(thresholds["upper_threshold"])) {
      parts <- c(parts, sprintf("upper = %.2f", thresholds["upper_threshold"]))
    }
    cat(sprintf(
      "  - %s: %d outliers [%s]\n",
      nm,
      n,
      paste(parts, collapse = ", ")
    ))
  }

  invisible(x)
}

#' Plot per-cell QC violin plots from a CellQc object
#'
#' @param x A `CellQc` object.
#' @param qc_df A data.table containing the cell-level data.
#' @param ... Ignored.
#'
#' @return A named list of ggplot objects, one per metric.
#'
#' @export
#'
#' @import ggplot2
#'
#' @keywords internal
plot.CellQc <- function(x, qc_df, ...) {
  outlier_colours <- c("FALSE" = "lightgrey", "TRUE" = "orange")

  plots <- purrr::map(names(x$metrics), function(nm) {
    plot_dt <- data.table::copy(qc_df)[, outlier := x$combined]

    ggplot2::ggplot(plot_dt, ggplot2::aes(y = x$metrics[[nm]], x = "cells")) +
      ggplot2::geom_violin() +
      ggplot2::geom_jitter(
        ggplot2::aes(colour = outlier),
        width = 0.05,
        size = 0.4,
        alpha = 0.5,
        show.legend = FALSE
      ) +
      ggplot2::scale_colour_manual(values = outlier_colours) +
      ggplot2::ylab(nm) +
      ggplot2::xlab("") +
      ggplot2::theme_bw()
  })

  setNames(plots, names(x$metrics))
}
