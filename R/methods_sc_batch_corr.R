# methods related to batch corrections -----------------------------------------

## kBET scores -----------------------------------------------------------------

#' Calculate kBET scores
#'
#' @description
#' This function calculates the k-nearest neighbour batch-effect test (kBET).
#' Briefly, the function leverages a Chi-Square statistic to calculate the
#' differences in batch proportions observed in the neighbourhood of a given
#' cell with the overall batch proportions. If the test is significant for that
#' cell it indicates poor mixing for that cell specifically.
#' Large number of positive tests indicate bad mixing overall. For more details,
#' please see Büttner et al.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param batch_column String. The column with the batch information in the
#' obs data of the class.
#' @param threshold Numeric. Number between 0 and 1. Below this threshold, the
#' test is considered significant. Defaults to `0.05`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A `KbetScores` object with the following elements
#' \itemize{
#'   \item kbet_score - Proportion of significant tests over all cells. 0
#'   indicates perfect mixing, 1 indicates no mixing between batches.
#'   \item significant_tests - Logical vector indicating for which cells the
#'   test was below the threshold.
#'   \item p_values - The p-values from the Chi-Square test.
#'   \item chi_square_stats - Per-cell Chi-Square statistics.
#'   \item mean_chi_square - Mean Chi-Square statistic across all cells.
#'   \item median_chi_square - Median Chi-Square statistic across all cells.
#'   \item threshold - The significance threshold used.
#'   \item n_batches - Number of batches in the data.
#' }
#'
#' @export
#'
#' @references Büttner, et al., Nat. Methods, 2019
#'
#' @examples
#' # kBET rejection rate over three batches
#' sc <- demo_single_cells(
#'   syn_data_params = params_sc_synthetic_data(
#'     n_cells = 600L, n_genes = 50L, n_batches = 3L
#'   )
#' )
#' calculate_kbet_sc(sc, batch_column = "batch_index", .verbose = FALSE)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
calculate_kbet_sc <- S7::new_generic(
  name = "calculate_kbet_sc",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    threshold = 0.05,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calculate_kbet_sc ScOrScSubset
S7::method(calculate_kbet_sc, ScOrScSubset) <- function(
  object,
  batch_column,
  threshold = 0.05,
  .verbose = TRUE
) {
  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(threshold, "N1[0, 1]")
  checkmate::qassert(.verbose, "B1")

  batch_index <- unlist(object[[batch_column]])

  if (!length(levels(factor(batch_index))) > 1) {
    warning("The batch column only has one batch. Returning NULL")
    return(NULL)
  }

  knn_mat <- get_knn_mat(object)

  if (is.null(knn_mat)) {
    warning("No kNN matrix found to calculate kBET. Returning NULL")
    return(NULL)
  }

  n_batches <- length(levels(factor(batch_index)))

  rs_res <- rs_kbet(
    knn_mat = knn_mat,
    batch_vector = as.integer(factor(batch_index)),
    verbose = .verbose
  )

  structure(
    list(
      kbet_score = sum(rs_res$pval <= threshold) / length(rs_res$pval),
      significant_tests = rs_res$pval <= threshold,
      p_values = rs_res$pval,
      chi_square_stats = rs_res$chi_square_stats,
      mean_chi_square = rs_res$mean_chi_square,
      median_chi_square = rs_res$median_chi_square,
      threshold = threshold,
      n_batches = n_batches
    ),
    class = "KbetScores"
  )
}

### kbet print -----------------------------------------------------------------

#' Print method for KbetScores
#'
#' @param x A `KbetScores` object.
#' @param ... Additional arguments (ignored).
#'
#' @export
#'
#' @keywords internal
print.KbetScores <- function(x, ...) {
  n_cells <- length(x$p_values)
  n_sig <- sum(x$significant_tests)
  dof <- x$n_batches - 1

  cat("kBET Scores\n")
  cat(sprintf(
    "  Cells: %d | Batches: %d | Threshold: %.3f\n",
    n_cells,
    x$n_batches,
    x$threshold
  ))
  cat(sprintf(
    "  Rejection rate:      %.4f (%d / %d)\n",
    x$kbet_score,
    n_sig,
    n_cells
  ))
  cat(sprintf(
    "  Mean Chi-Square:     %.4f (expected under H0: %d)\n",
    x$mean_chi_square,
    dof
  ))
  cat(sprintf("  Median Chi-Square:   %.4f\n", x$median_chi_square))

  invisible(x)
}

## batch silhouette width ------------------------------------------------------

#' Calculate batch average silhouette width
#'
#' @description
#' Computes the average silhouette width (ASW) on batch labels using pairwise
#' Euclidean distances in the embedding space. For each cell, the function
#' estimates the mean distance to cells of the same batch (a) and the mean
#' distance to cells of the nearest other batch (b), then computes
#' s = (b - a) / max(a, b).
#'
#' Values near 0 indicate good batch mixing, values near 1 indicate batch
#' separation, and negative values suggest overcorrection. This metric is
#' best suited for embedding-based correction methods (e.g. Harmony, fastMNN).
#' For graph-based methods like BBKNN, consider using
#' [bixverse::calculate_lisi_sc()] instead.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param batch_column String. The column with the batch information in the
#' obs data of the class.
#' @param embd_to_use String. Which embedding to compute the ASW on. One of
#' `c("pca", "harmony", "mnn")`. Defaults to `"pca"`.
#' @param max_cells Integer or `NULL`. If not `NULL`, subsample to this many
#' cells for performance. The pairwise distance computation is O(n^2), so
#' subsampling is recommended for large datasets. Defaults to `5000L`.
#' @param seed Integer. Seed for subsampling reproducibility.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A `BatchSilhouetteScores` object with the following elements
#' \itemize{
#'   \item per_cell - Per-cell silhouette scores in `[-1, 1]`.
#'   \item mean_asw - Mean silhouette width across all cells.
#'   \item median_asw - Median silhouette width across all cells.
#'   \item n_batches - Number of batches in the data.
#'   \item embedding_used - Which embedding the ASW was computed on.
#' }
#'
#' @export
#'
#' @examples
#' # batch silhouette width on the PCA embedding
#' sc <- demo_single_cells(
#'   syn_data_params = params_sc_synthetic_data(
#'     n_cells = 600L, n_genes = 50L, n_batches = 3L
#'   )
#' )
#' calculate_batch_asw_sc(
#'   sc,
#'   batch_column = "batch_index",
#'   .verbose = FALSE
#' )
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
calculate_batch_asw_sc <- S7::new_generic(
  name = "calculate_batch_asw_sc",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    embd_to_use = "pca",
    max_cells = 5000L,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calculate_batch_asw_sc ScOrScSubset
S7::method(calculate_batch_asw_sc, ScOrScSubset) <- function(
  object,
  batch_column,
  embd_to_use = "pca",
  max_cells = 5000L,
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(max_cells, c("I1", "0"))
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  batch_index <- unlist(object[[batch_column]])

  if (!length(levels(factor(batch_index))) > 1) {
    warning("The batch column only has one batch. Returning NULL")
    return(NULL)
  }

  if (!embd_to_use %in% get_available_embeddings(object)) {
    warning("The desired embedding was not found. Returning NULL")
    return(NULL)
  }

  embd <- get_embedding(x = object, embd_name = embd_to_use)
  n_batches <- length(levels(factor(batch_index)))

  rs_res <- rs_batch_silhouette_width(
    embedding = embd,
    batch_vector = as.integer(factor(batch_index)),
    max_cells = max_cells,
    seed = seed,
    verbose = .verbose
  )

  structure(
    list(
      per_cell = rs_res$per_cell,
      mean_asw = rs_res$mean_asw,
      median_asw = rs_res$median_asw,
      n_batches = n_batches,
      embedding_used = embd_to_use
    ),
    class = "BatchSilhouetteScores"
  )
}

### batch silhouette print -----------------------------------------------------

#' Print method for BatchSilhouetteScores
#'
#' @param x A `BatchSilhouetteScores` object.
#' @param ... Additional arguments (ignored).
#'
#' @export
#'
#' @keywords internal
print.BatchSilhouetteScores <- function(x, ...) {
  n_cells <- length(x$per_cell)

  cat("Batch Silhouette Width\n")
  cat(sprintf("  Cells: %d | Batches: %d\n", n_cells, x$n_batches))
  cat(sprintf(
    "  Mean ASW:    %.4f (-1 = strong intermixing, 0 = mixed, 1 = separated)\n",
    x$mean_asw
  ))
  cat(sprintf("  Median ASW:  %.4f\n", x$median_asw))

  invisible(x)
}

## LISI ------------------------------------------------------------------------

#' Calculate LISI scores (iLISI or cLISI)
#'
#' @description
#' Computes the Local Inverse Simpson's Index (LISI) on the kNN graph: the
#' effective number of labels in each cell's neighbourhood. On batch labels
#' this is iLISI, where higher means better mixing. On cell type labels it is
#' cLISI, where lower means cell types stay apart. Unlike kBET, LISI does not
#' compare against global proportions, so it also works on graph-based
#' corrections like BBKNN.
#'
#' The normalised score follows scIB and lands in `[0, 1]`, higher is better
#' for both: iLISI as `(median - 1) / (n - 1)`, cLISI as
#' `(n - median) / (n - 1)`.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param label_column String. The column with the batch or cell type labels in
#' the obs data of the class.
#' @param type String. One of `c("batch", "cell_type")`. Decides which
#' normalised score is reported. Defaults to `"batch"`.
#' @param weighted Boolean. Weight the neighbours with a perplexity-calibrated
#' Gaussian kernel on the kNN distances, as in Korsunsky et al. If `FALSE`,
#' all neighbours count equally. Defaults to `FALSE`.
#' @param perplexity Numeric. Perplexity for the weighted version. Defaults to
#' `30`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A `LisiScores` object with the following elements
#' \itemize{
#'   \item per_cell - Per-cell LISI scores in `[1, n_labels]`.
#'   \item mean_lisi - Mean LISI across all cells.
#'   \item median_lisi - Median LISI across all cells.
#'   \item lisi_norm - The normalised score in `[0, 1]`, higher is better.
#'   \item n_labels - Number of distinct labels.
#'   \item type - `"batch"` (iLISI) or `"cell_type"` (cLISI).
#' }
#'
#' @export
#'
#' @references Korsunsky, et al., Nat. Methods, 2019; Luecken, et al., Nat.
#' Methods, 2022
#'
#' @examples
#' # iLISI over the kNN graph
#' sc <- demo_single_cells(
#'   syn_data_params = params_sc_synthetic_data(
#'     n_cells = 600L, n_genes = 50L, n_batches = 3L
#'   )
#' )
#' calculate_lisi_sc(
#'   sc,
#'   label_column = "batch_index",
#'   .verbose = FALSE
#' )
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
calculate_lisi_sc <- S7::new_generic(
  name = "calculate_lisi_sc",
  dispatch_args = "object",
  fun = function(
    object,
    label_column,
    type = c("batch", "cell_type"),
    weighted = FALSE,
    perplexity = 30,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calculate_lisi_sc ScOrScSubset
S7::method(calculate_lisi_sc, ScOrScSubset) <- function(
  object,
  label_column,
  type = c("batch", "cell_type"),
  weighted = FALSE,
  perplexity = 30,
  .verbose = TRUE
) {
  type <- match.arg(type)
  checkmate::qassert(label_column, "S1")
  checkmate::assertChoice(type, c("batch", "cell_type"))
  checkmate::qassert(weighted, "B1")
  checkmate::qassert(perplexity, "N1(0,)")
  checkmate::qassert(.verbose, "B1")

  labels <- unlist(object[[label_column]])

  if (!length(levels(factor(labels))) > 1) {
    warning("The label column only has one label. Returning NULL")
    return(NULL)
  }

  knn_mat <- get_knn_mat(object)

  if (is.null(knn_mat)) {
    warning("No kNN matrix found to calculate LISI. Returning NULL")
    return(NULL)
  }

  knn_dist <- if (weighted) get_knn_dist(object) else NULL
  if (weighted && is.null(knn_dist)) {
    warning("No kNN distances found for weighted LISI. Returning NULL")
    return(NULL)
  }

  rs_res <- rs_lisi(
    knn_mat = knn_mat,
    knn_dist = knn_dist,
    labels = as.integer(factor(labels)),
    perplexity = perplexity,
    verbose = .verbose
  )

  structure(
    list(
      per_cell = rs_res$per_cell,
      mean_lisi = rs_res$mean_lisi,
      median_lisi = rs_res$median_lisi,
      lisi_norm = if (type == "batch") {
        rs_res$ilisi_norm
      } else {
        rs_res$clisi_norm
      },
      n_labels = rs_res$n_labels,
      type = type
    ),
    class = "LisiScores"
  )
}

### lisi print -----------------------------------------------------------------

#' @export
#'
#' @keywords internal
print.LisiScores <- function(x, ...) {
  n_cells <- length(x$per_cell)
  label <- if (x$type == "batch") "iLISI (batch)" else "cLISI (cell type)"

  cat(sprintf("%s\n", label))
  cat(sprintf("  Cells: %d | Labels: %d\n", n_cells, x$n_labels))
  cat(sprintf("  Mean LISI:    %.4f\n", x$mean_lisi))
  cat(sprintf("  Median LISI:  %.4f\n", x$median_lisi))
  cat(sprintf("  Normalised:   %.4f (0 = worst, 1 = best)\n", x$lisi_norm))

  invisible(x)
}

## PCR -------------------------------------------------------------------------

#' Calculate the principal component regression on batch
#'
#' @description
#' Regresses each embedding dimension on the batch labels and weights the
#' per-dimension R-squared by the variance of that dimension. On its own the
#' number says how much of the embedding variance batch explains. For a
#' corrected embedding, the function also runs it on the uncorrected PCA and
#' reports the scIB comparison `(pre - post) / pre`: 1 means the batch
#' variance is gone, 0 means nothing changed, negative means it got worse.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param batch_column String. The column with the batch information in the
#' obs data of the class.
#' @param embd_to_use String. Which embedding to compute the PCR on. Defaults
#' to `"pca"`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A `PcrScores` object with the following elements
#' \itemize{
#'   \item pcr - Variance-weighted R-squared of batch on `embd_to_use`.
#'   \item pcr_pca - The same on the uncorrected PCA.
#'   \item pcr_comparison - `(pcr_pca - pcr) / pcr_pca`. `NA` if
#'   `embd_to_use = "pca"`.
#'   \item var_explained - Variance per embedding dimension.
#'   \item r_squared - Batch R-squared per embedding dimension.
#'   \item embedding_used - Which embedding the PCR was computed on.
#' }
#'
#' @export
#'
#' @references Büttner, et al., Nat. Methods, 2019; Luecken, et al., Nat.
#' Methods, 2022
#'
#' @examples
#' # share of PCA variance explained by batch
#' sc <- demo_single_cells(
#'   syn_data_params = params_sc_synthetic_data(
#'     n_cells = 600L, n_genes = 50L, n_batches = 3L
#'   )
#' )
#' calculate_pcr_sc(sc, batch_column = "batch_index")
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
calculate_pcr_sc <- S7::new_generic(
  name = "calculate_pcr_sc",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    embd_to_use = "pca",
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calculate_pcr_sc ScOrScSubset
S7::method(calculate_pcr_sc, ScOrScSubset) <- function(
  object,
  batch_column,
  embd_to_use = "pca",
  .verbose = TRUE
) {
  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(.verbose, "B1")

  batch_index <- unlist(object[[batch_column]])

  if (!length(levels(factor(batch_index))) > 1) {
    warning("The batch column only has one batch. Returning NULL")
    return(NULL)
  }

  available <- get_available_embeddings(object)
  if (!all(c("pca", embd_to_use) %in% available)) {
    warning("The PCA or the desired embedding was not found. Returning NULL")
    return(NULL)
  }

  batch_vector <- as.integer(factor(batch_index))
  pca_res <- rs_pcr(
    embedding = get_embedding(x = object, embd_name = "pca"),
    batch_vector = batch_vector
  )
  rs_res <- if (embd_to_use == "pca") {
    pca_res
  } else {
    rs_pcr(
      embedding = get_embedding(x = object, embd_name = embd_to_use),
      batch_vector = batch_vector
    )
  }

  pcr_comparison <- if (embd_to_use == "pca") {
    NA_real_
  } else {
    (pca_res$pcr - rs_res$pcr) / pca_res$pcr
  }

  structure(
    list(
      pcr = rs_res$pcr,
      pcr_pca = pca_res$pcr,
      pcr_comparison = pcr_comparison,
      var_explained = rs_res$var_explained,
      r_squared = rs_res$r_squared,
      embedding_used = embd_to_use
    ),
    class = "PcrScores"
  )
}

### pcr print ------------------------------------------------------------------

#' @export
#'
#' @keywords internal
print.PcrScores <- function(x, ...) {
  cat("Principal Component Regression (batch)\n")
  cat(sprintf(
    "  Embedding: %s | Dimensions: %d\n",
    x$embedding_used,
    length(x$r_squared)
  ))
  cat(sprintf("  PCR:             %.4f\n", x$pcr))
  if (!is.na(x$pcr_comparison)) {
    cat(sprintf("  PCR (PCA):       %.4f\n", x$pcr_pca))
    cat(sprintf(
      "  PCR comparison:  %.4f (1 = batch variance removed)\n",
      x$pcr_comparison
    ))
  }

  invisible(x)
}

## cell type silhouette width --------------------------------------------------

#' Calculate cell type average silhouette width
#'
#' @description
#' Average silhouette width on cell type labels in the embedding, rescaled to
#' `[0, 1]` via `(s + 1) / 2` as in scIB. Higher values mean cell types stay
#' separated after correction. Counterpart to
#' [bixverse::calculate_batch_asw_sc()] on the biology side.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param cell_type_column String. The column with the cell type labels in the
#' obs data of the class.
#' @param embd_to_use String. Which embedding to compute the ASW on. Defaults
#' to `"pca"`.
#' @param max_cells Integer or `NULL`. If not `NULL`, subsample to this many
#' cells for performance. Defaults to `5000L`.
#' @param seed Integer. Seed for subsampling reproducibility.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A `CellTypeAswScores` object with the following elements
#' \itemize{
#'   \item per_cell - Per-cell rescaled silhouette scores in `[0, 1]`.
#'   \item mean_asw - Mean rescaled silhouette width.
#'   \item median_asw - Median rescaled silhouette width.
#'   \item n_cell_types - Number of cell types.
#'   \item embedding_used - Which embedding the ASW was computed on.
#' }
#'
#' @export
#'
#' @references Luecken, et al., Nat. Methods, 2022
#'
#' @examples
#' # cell type silhouette width on the PCA embedding
#' sc <- demo_single_cells(
#'   syn_data_params = params_sc_synthetic_data(
#'     n_cells = 600L, n_genes = 50L, n_batches = 3L
#'   )
#' )
#' calculate_cell_type_asw_sc(
#'   sc,
#'   cell_type_column = "cell_grp",
#'   .verbose = FALSE
#' )
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
calculate_cell_type_asw_sc <- S7::new_generic(
  name = "calculate_cell_type_asw_sc",
  dispatch_args = "object",
  fun = function(
    object,
    cell_type_column,
    embd_to_use = "pca",
    max_cells = 5000L,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calculate_cell_type_asw_sc ScOrScSubset
S7::method(calculate_cell_type_asw_sc, ScOrScSubset) <- function(
  object,
  cell_type_column,
  embd_to_use = "pca",
  max_cells = 5000L,
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::qassert(cell_type_column, "S1")
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(max_cells, c("I1", "0"))
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  labels <- unlist(object[[cell_type_column]])

  if (!length(levels(factor(labels))) > 1) {
    warning("The cell type column only has one cell type. Returning NULL")
    return(NULL)
  }

  if (!embd_to_use %in% get_available_embeddings(object)) {
    warning("The desired embedding was not found. Returning NULL")
    return(NULL)
  }

  rs_res <- rs_cell_type_asw(
    embedding = get_embedding(x = object, embd_name = embd_to_use),
    labels = as.integer(factor(labels)),
    max_cells = max_cells,
    verbose = .verbose,
    seed = seed
  )

  structure(
    list(
      per_cell = rs_res$per_cell,
      mean_asw = rs_res$mean_asw,
      median_asw = rs_res$median_asw,
      n_cell_types = length(levels(factor(labels))),
      embedding_used = embd_to_use
    ),
    class = "CellTypeAswScores"
  )
}

### cell type silhouette print -------------------------------------------------

#' @export
#'
#' @keywords internal
print.CellTypeAswScores <- function(x, ...) {
  cat("Cell Type Silhouette Width (rescaled)\n")
  cat(sprintf(
    "  Cells: %d | Cell types: %d | Embedding: %s\n",
    length(x$per_cell),
    x$n_cell_types,
    x$embedding_used
  ))
  cat(sprintf(
    "  Mean ASW:    %.4f (0.5 = no structure, 1 = separated)\n",
    x$mean_asw
  ))
  cat(sprintf("  Median ASW:  %.4f\n", x$median_asw))

  invisible(x)
}

## graph connectivity ----------------------------------------------------------

#' Calculate the graph connectivity per cell type
#'
#' @description
#' For each cell type, restricts the kNN graph to that cell type and takes the
#' fraction of its cells in the largest connected component. A cell type split
#' across batches after correction falls apart into several components and
#' scores low. 1 means every cell type is one connected piece.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param cell_type_column String. The column with the cell type labels in the
#' obs data of the class.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A `GraphConnectivityScores` object with the following elements
#' \itemize{
#'   \item per_cell_type - Named numeric. Connectivity per cell type.
#'   \item mean_connectivity - Mean connectivity across cell types.
#'   \item median_connectivity - Median connectivity across cell types.
#' }
#'
#' @export
#'
#' @references Luecken, et al., Nat. Methods, 2022
#'
#' @examples
#' # connectivity of each cell type in the kNN graph
#' sc <- demo_single_cells(
#'   syn_data_params = params_sc_synthetic_data(
#'     n_cells = 600L, n_genes = 50L, n_batches = 3L
#'   )
#' )
#' calculate_graph_connectivity_sc(sc, cell_type_column = "cell_grp")
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
calculate_graph_connectivity_sc <- S7::new_generic(
  name = "calculate_graph_connectivity_sc",
  dispatch_args = "object",
  fun = function(
    object,
    cell_type_column,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calculate_graph_connectivity_sc ScOrScSubset
S7::method(calculate_graph_connectivity_sc, ScOrScSubset) <- function(
  object,
  cell_type_column,
  .verbose = TRUE
) {
  checkmate::qassert(cell_type_column, "S1")
  checkmate::qassert(.verbose, "B1")

  labels <- factor(unlist(object[[cell_type_column]]))

  knn_mat <- get_knn_mat(object)

  if (is.null(knn_mat)) {
    warning("No kNN matrix found to calculate connectivity. Returning NULL")
    return(NULL)
  }

  rs_res <- rs_graph_connectivity(
    knn_mat = knn_mat,
    labels = as.integer(labels)
  )

  structure(
    list(
      # Rust reports the labels in first-appearance order, not level order
      per_cell_type = stats::setNames(
        rs_res$per_label,
        levels(labels)[unique(as.integer(labels))]
      ),
      mean_connectivity = rs_res$mean,
      median_connectivity = rs_res$median
    ),
    class = "GraphConnectivityScores"
  )
}

### graph connectivity print ---------------------------------------------------

#' @export
#'
#' @keywords internal
print.GraphConnectivityScores <- function(x, ...) {
  cat("Graph Connectivity\n")
  cat(sprintf("  Cell types: %d\n", length(x$per_cell_type)))
  cat(sprintf(
    "  Mean:    %.4f (1 = every cell type connected)\n",
    x$mean_connectivity
  ))
  cat(sprintf("  Median:  %.4f\n", x$median_connectivity))
  worst <- sort(x$per_cell_type)[seq_len(min(3L, length(x$per_cell_type)))]
  cat(sprintf(
    "  Lowest:  %s\n",
    paste(sprintf("%s (%.3f)", names(worst), worst), collapse = ", ")
  ))

  invisible(x)
}

## summary ---------------------------------------------------------------------

#' Calculate a summary of integration metrics
#'
#' @description
#' Runs the batch mixing and (if cell type labels are given) the biological
#' conservation metrics in one go and returns one row per call, so results
#' across correction methods can be `rbind`-ed into one table. Every column is
#' on `[0, 1]` (PCR comparison can go negative) and higher is better, the
#' scIB convention:
#'
#' Batch mixing:
#' \itemize{
#'   \item kbet_accept - `1 - ` kBET rejection rate.
#'   \item batch_asw - `mean(1 - |s|)` over the per-cell batch silhouettes.
#'   \item ilisi - Normalised iLISI.
#'   \item pcr_comparison - `(pre - post) / pre` of the batch PCR.
#' }
#'
#' Biological conservation:
#' \itemize{
#'   \item clisi - Normalised cLISI.
#'   \item cell_type_asw - Rescaled cell type silhouette width.
#'   \item graph_connectivity - Mean graph connectivity over cell types.
#' }
#'
#' The kNN metrics read the kNN graph currently stored in the object, so
#' recompute the neighbours on the corrected embedding first. Embedding metrics
#' are `NA` if `embd_to_use = NULL` (e.g. BBKNN, which only returns a graph).
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param batch_column String. The column with the batch information in the
#' obs data of the class.
#' @param cell_type_column Optional string. The column with the cell type
#' labels. If `NULL`, the conservation metrics are `NA`.
#' @param embd_to_use Optional string. The embedding for ASW and PCR. Defaults
#' to `"pca"`.
#' @param max_cells Integer or `NULL`. Subsampling for the silhouette widths.
#' Defaults to `5000L`.
#' @param seed Integer. Seed for subsampling reproducibility.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A one-row data.table with the columns `embedding`, `kbet_accept`,
#' `batch_asw`, `ilisi`, `pcr_comparison`, `clisi`, `cell_type_asw` and
#' `graph_connectivity`.
#'
#' @export
#'
#' @references Luecken, et al., Nat. Methods, 2022
#'
#' @examples
#' # all metrics on the uncorrected PCA
#' sc <- demo_single_cells(
#'   syn_data_params = params_sc_synthetic_data(
#'     n_cells = 600L, n_genes = 50L, n_batches = 3L
#'   )
#' )
#' calculate_integration_metrics_sc(
#'   sc,
#'   batch_column = "batch_index",
#'   cell_type_column = "cell_grp",
#'   .verbose = FALSE
#' )
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
calculate_integration_metrics_sc <- S7::new_generic(
  name = "calculate_integration_metrics_sc",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    cell_type_column = NULL,
    embd_to_use = "pca",
    max_cells = 5000L,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calculate_integration_metrics_sc ScOrScSubset
S7::method(calculate_integration_metrics_sc, ScOrScSubset) <- function(
  object,
  batch_column,
  cell_type_column = NULL,
  embd_to_use = "pca",
  max_cells = 5000L,
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(cell_type_column, c("S1", "0"))
  checkmate::qassert(embd_to_use, c("S1", "0"))
  checkmate::qassert(max_cells, c("I1", "0"))
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  # NULL results (missing kNN / embedding) become NA
  or_na <- \(x, expr) if (is.null(x)) NA_real_ else expr(x)

  kbet <- calculate_kbet_sc(
    object,
    batch_column = batch_column,
    .verbose = FALSE
  )
  ilisi <- calculate_lisi_sc(
    object,
    label_column = batch_column,
    type = "batch",
    .verbose = FALSE
  )

  batch_asw <- pcr <- NULL
  if (!is.null(embd_to_use)) {
    batch_asw <- calculate_batch_asw_sc(
      object,
      batch_column = batch_column,
      embd_to_use = embd_to_use,
      max_cells = max_cells,
      seed = seed,
      .verbose = FALSE
    )
    pcr <- calculate_pcr_sc(
      object,
      batch_column = batch_column,
      embd_to_use = embd_to_use,
      .verbose = FALSE
    )
  }

  clisi <- ct_asw <- conn <- NULL
  if (!is.null(cell_type_column)) {
    clisi <- calculate_lisi_sc(
      object,
      label_column = cell_type_column,
      type = "cell_type",
      .verbose = FALSE
    )
    conn <- calculate_graph_connectivity_sc(
      object,
      cell_type_column = cell_type_column,
      .verbose = FALSE
    )
    if (!is.null(embd_to_use)) {
      ct_asw <- calculate_cell_type_asw_sc(
        object,
        cell_type_column = cell_type_column,
        embd_to_use = embd_to_use,
        max_cells = max_cells,
        seed = seed,
        .verbose = FALSE
      )
    }
  }

  data.table::data.table(
    embedding = if (is.null(embd_to_use)) NA_character_ else embd_to_use,
    kbet_accept = or_na(kbet, \(x) 1 - x$kbet_score),
    batch_asw = or_na(batch_asw, \(x) mean(1 - abs(x$per_cell))),
    ilisi = or_na(ilisi, \(x) x$lisi_norm),
    pcr_comparison = or_na(pcr, \(x) x$pcr_comparison),
    clisi = or_na(clisi, \(x) x$lisi_norm),
    cell_type_asw = or_na(ct_asw, \(x) x$mean_asw),
    graph_connectivity = or_na(conn, \(x) x$mean_connectivity)
  )
}

## batch aware hvg -------------------------------------------------------------

#' Identify HVGs (batch aware)
#'
#' @description
#' This is a helper function to identify highly variable genes in a batch-aware
#' manner. At the moment the implementation has only the VST-based version
#' (known as Seurat v3). The other methods will be implemented in the future.
#' This function will calculate the HVG per given experimental batch and you
#' can choose the way to combine them. The choices are union (of Top x HVG per
#' batch), based on the average variance per batch or only take genes that are
#' amongst the Top X HVG in all batches. Important. The function returns
#' 0-indices for the genes!
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param batch_column String. The column name of the batch column in the obs
#' table.
#' @param hvg_no Integer. Number of highly variable genes to include. Defaults
#' to `2000L`.
#' @param gene_comb_method String. One of
#' `c("union", "average", "intersection")`. The method to combine the HVG across
#' the different batches. Defaults to `"union"`.
#' @param hvg_params List, see [bixverse::params_sc_hvg()]. This list contains
#' \itemize{
#'   \item method - Which method to use. One of
#'   `c("vst", "meanvarbin", "dispersion")`
#'   \item loess_span - The span for the loess function to standardise the
#'   variance
#'   \item num_bin - Integer. Not yet implemented.
#'   \item bin_method - String. One of `c("equal_width", "equal_freq")`. Not
#'   implemented yet.
#' }
#' @param streaming Optional Boolean. Shall the data be streamed in. Useful for
#' larger data sets where you wish to avoid loading in the whole data. If
#' `NULL`, will automatically detect.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns This function will return a list with:
#' \itemize{
#'   \item hvg_genes - The gene names of the HVGs.
#'   \item hvg_gene_idx - The (0-index) gene features.
#'   \item batch_hvg_data - data.table with the detailed information of the
#'   variance per batch.
#' }
#'
#' @export
#'
#' @examples
#' # highly variable genes taken as the union over the batches
#' sc <- demo_single_cells(
#'   syn_data_params = params_sc_synthetic_data(
#'     n_cells = 600L, n_genes = 50L, n_batches = 3L
#'   )
#' )
#' hvg <- find_hvg_batch_aware_sc(
#'   sc,
#'   hvg_no = 20L,
#'   batch_column = "batch_index",
#'   gene_comb_method = "union",
#'   .verbose = FALSE
#' )
#' head(hvg$hvg_genes)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
find_hvg_batch_aware_sc <- S7::new_generic(
  name = "find_hvg_batch_aware_sc",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    hvg_no = 2000L,
    gene_comb_method = c("union", "average", "intersection"),
    hvg_params = params_sc_hvg(),
    streaming = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method find_hvg_batch_aware_sc ScOrScSubset
S7::method(find_hvg_batch_aware_sc, ScOrScSubset) <- function(
  object,
  batch_column,
  hvg_no = 2000L,
  gene_comb_method = c("union", "average", "intersection"),
  hvg_params = params_sc_hvg(),
  streaming = NULL,
  .verbose = TRUE
) {
  gene_comb_method <- match.arg(gene_comb_method)

  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(hvg_no, "I1")
  checkmate::assertChoice(
    gene_comb_method,
    c("union", "average", "intersection")
  )
  assertScHvg(hvg_params)
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # the residual path does its own per-group split, at fitting time
  if (hvg_params$method == "residual") {
    stop(paste(
      "find_hvg_batch_aware_sc() does not support method = 'residual'.",
      "Fit per batch instead with",
      "fit_residuals_sc(group_column = ...), then call find_hvg_sc():",
      "the per-group ranking and union happen there."
    ))
  }

  streaming <- auto_streaming(
    n_cells = nrow(object),
    streaming = streaming,
    .verbose = .verbose
  )

  batch_indices <- unlist(object[[batch_column]])
  batch_factor <- factor(batch_indices)
  batch_indices <- as.integer(batch_factor) - 1L

  batch_hvgs <- with(
    hvg_params,
    rs_sc_hvg_batch_aware(
      f_path_gene = get_rust_count_gene_f_path(object),
      hvg_method = method,
      cell_indices = get_cells_to_keep(object),
      batch_labels = batch_indices,
      loess_span = loess_span,
      clip_max = NULL,
      n_bins = num_bin,
      binning = bin_method,
      streaming = streaming,
      verbose = parse_verbosity(.verbose)
    )
  )

  batch_hvgs_dt <- data.table::as.data.table(batch_hvgs)
  batch_hvgs_dt[, batch := levels(batch_factor)[batch + 1L]]

  sort_col <- switch(
    hvg_params$method,
    "vst" = "var_std",
    "dispersion" = "dispersion",
    "meanvarbin" = "dispersion_scaled",
    stop("Unknown HVG method: ", hvg_params$method)
  )

  hvg_gene_idx <- switch(
    gene_comb_method,
    union = {
      batch_hvgs_dt[,
        .SD[order(-score)][1:hvg_no],
        by = batch,
        env = list(score = sort_col)
      ][, unique(gene_idx)]
    },
    average = {
      avg_dt <- batch_hvgs_dt[,
        .(score_avg = mean(score)),
        by = gene_idx,
        env = list(score = sort_col)
      ]
      avg_dt[order(-score_avg)][1:hvg_no, gene_idx]
    },
    intersection = {
      top_per_batch <- batch_hvgs_dt[,
        .(gene_idx = .SD[order(-score)][1:hvg_no, gene_idx]),
        by = batch,
        env = list(score = sort_col)
      ]
      top_per_batch[, .N, by = gene_idx][
        N == uniqueN(batch_hvgs_dt$batch),
        gene_idx
      ]
    }
  )
  hvg_genes <- get_gene_names_from_idx(x = object, gene_idx = hvg_gene_idx)

  list(
    hvg_genes = hvg_genes,
    hvg_gene_idx = hvg_gene_idx,
    hvg_data = batch_hvgs_dt
  )
}

## BBKNN -----------------------------------------------------------------------

#' Run BBKNN
#'
#' @description
#' This function implements the batch-balanced k-nearest neighbour algorithm
#' from Polański, et al. Briefly, the algorithm generates a KNN index on a per
#' batch basis and identifies the neighbours of cells for each individual index.
#' Subsequently, it leverages UMAP connectivity calculations to reduce spurious
#' connections. For more details, please refer to Polański, et al.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param batch_column String. The column with the batch information in the
#' obs data of the class.
#' @param no_neighbours_to_keep Integer. Maximum number of neighbours to keep
#' from the BBKNN algorithm. Due to generating neighbours for each batch, there
#' might be a large number of generated neighbours. This will only keep the top
#' `no_neighbours_to_keep` neighbours. Defaults to `5L`.
#' @param embd_to_use String. The embedding to use. Atm, the only option is
#' `"pca"`.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param bbknn_params A list, please see [bixverse::params_sc_bbknn()]. The
#' list has the following parameters:
#' \itemize{
#'   \item neighbours_within_batch - Integer. Number of neighbours to consider
#'   per batch.
#'   \item set_op_mix_ratio - Numeric. Mixing ratio between union (1.0) and
#'   intersection (0.0).
#'   \item local_connectivity - Numeric. UMAP connectivity computation
#'   parameter, how many nearest neighbours of each cell are assumed to be fully
#'   connected.
#'   \item trim - Optional integer. Trim the neighbours of each cell to these
#'   many top connectivities. May help with population independence and improve
#'   the tidiness of clustering. If `NULL`, it defaults to
#'   `10 * neighbours_within_batch`.
#'   \item knn - List of kNN parameters. See [bixverse::params_knn_defaults()]
#'   for available parameters and their defaults.
#' }
#' @param seed Integer. Random seed.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns The object with added kNN matrix based on BBKNN and the graph based
#' on the returned connectivities of the algorithm.
#'
#' @export
#'
#' @references Polański, et al., Bioinformatics, 2020
#'
#' @examples
#' # batch balanced kNN, replacing the graph the demo object carries
#' sc <- demo_single_cells(
#'   syn_data_params = params_sc_synthetic_data(
#'     n_cells = 600L, n_genes = 50L, n_batches = 3L
#'   )
#' )
#' sc <- bbknn_sc(
#'   remove_knn(sc),
#'   batch_column = "batch_index",
#'   no_neighbours_to_keep = 5L,
#'   .verbose = FALSE
#' )
#' dim(get_knn_mat(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
bbknn_sc <- S7::new_generic(
  name = "bbknn_sc",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    no_neighbours_to_keep = 5L,
    embd_to_use = "pca",
    no_embd_to_use = NULL,
    bbknn_params = params_sc_bbknn(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method bbknn_sc ScOrScSubset
S7::method(bbknn_sc, ScOrScSubset) <- function(
  object,
  batch_column,
  no_neighbours_to_keep = 5L,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  bbknn_params = params_sc_bbknn(),
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::qassert(batch_column, "S1")
  checkmate::assertChoice(embd_to_use, c("pca"))
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  assertScBbknn(bbknn_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # presence probe, not a read: the kNN is about to be overwritten, so a stale
  # one here is not a problem worth signalling about
  if (.sc_has_artefact(object, "knn")) {
    warning("Prior kNN matrix found. Will be overwritten.")
  }

  # hard tier: the kNN built here feeds the sNN and everything downstream
  assert_sc_state(object, artefacts = embd_to_use)

  embd <- switch(embd_to_use, pca = get_pca_factors(object))

  if (is.null(embd)) {
    warning(paste(
      "The desired embedding was not found. Please check the parameters.",
      "Returning NULL."
    ))
    return(NULL)
  }

  if (!is.null(no_embd_to_use)) {
    to_take <- min(c(no_embd_to_use, ncol(embd)))
    embd <- embd[, 1:to_take]
  }

  batch_index <- unlist(object[[batch_column]])
  batch_factor <- factor(batch_index)
  batch_index <- as.integer(batch_factor) - 1L

  if (!length(levels(factor(batch_index))) > 1) {
    warning("The batch column only has one batch. Returning object as is.")
    return(object)
  }

  no_generated_neighbours <- length(levels(factor(batch_index))) *
    bbknn_params$neighbours_within_batch

  if (no_neighbours_to_keep > no_generated_neighbours) {
    warning(paste(
      "The number of desired neighbours cannot be generated with these BBKNN",
      "parameters (too few generated neighbours).",
      "Please adopt neighbours_within_batch accordingly.",
      "Returning all neighbours from BBKNN."
    ))
  }

  if (.verbose) {
    message("Running BBKNN algorithm.")
  }

  bbknn_res <- rs_bbknn(
    embd = embd,
    batch_labels = as.integer(batch_index),
    bbknn_params = bbknn_params,
    seed = seed,
    verbose = parse_verbosity(.verbose)
  )

  knn_data <- {
    no_k <- min(no_neighbours_to_keep, no_generated_neighbours)
    filtered <- rs_bbknn_filtering(
      indptr = bbknn_res$distances$indptr,
      indices = bbknn_res$distances$indices,
      data = bbknn_res$distances$data,
      no_neighbours_to_keep = no_k
    )
    list(
      indices = filtered$indices,
      dist = filtered$dist,
      dist_metric = bbknn_params[["ann_dist"]]
    )
  }

  storage.mode(knn_data$indices) <- "integer"

  used_cells <- get_cell_names(object, filtered = TRUE)
  sc_knn <- new_sc_knn(knn_data = knn_data, used_cells = used_cells)
  object <- set_knn(object, knn = sc_knn, from = "pca")

  if (.verbose) {
    message(paste(
      "Generating graph based on BBKNN connectivities.",
      "Weights will be based on the connectivities",
      "and not shared nearest neighbour calculations."
    ))
  }

  sparse_mat <- Matrix::sparseMatrix(
    i = rep(
      seq_along(bbknn_res$connectivities$indptr[-1]),
      diff(bbknn_res$connectivities$indptr)
    ),
    j = bbknn_res$connectivities$indices + 1,
    x = bbknn_res$connectivities$data,
    dims = c(bbknn_res$connectivities$nrow, bbknn_res$connectivities$ncol),
    index1 = TRUE
  )

  snn_graph <- igraph::graph_from_adjacency_matrix(
    sparse_mat,
    mode = "max",
    weighted = TRUE
  )

  object <- set_snn_graph(object, snn_graph = snn_graph, from = "knn")

  object
}

## fastMNN ---------------------------------------------------------------------

#' Run fastMNN
#'
#' @description
#' This function implements the fast mutual nearest neighbour (MNN) from
#' Haghverdi, et al. This version works on the PCA embedding and generates
#' an embedding only and not a fully corrected count matrix. The function will
#' iterate through the batches, identify the MNN and generate correction vectors
#' and generate a corrected embedding which is added to the function.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param batch_column String. The column with the batch information in the
#' obs data of the class.
#' @param batch_hvg_genes Integer vector. These are the highly variable genes,
#' identified by a batch-aware method. Please refer to
#' [bixverse::find_hvg_batch_aware_sc()] for more details. These genes have to
#' be 0-indexed!
#' @param fastmnn_params A list, please see [bixverse::params_sc_fastmnn()]. The
#' list has the following parameters:
#' \itemize{
#'   \item sigma - Numeric. Bandwidth of the Gaussian smoothing kernel (as
#'   proportion of space radius).
#'   \item cos_norm - Logical. Apply cosine normalisation before computing
#'   distances.
#'   \item var_adj - Logical. Apply variance adjustment to avoid kissing
#'   effects.
#'   \item no_pcs - Integer. Number of PCs to use for MNN calculations.
#'   \item knn - List of kNN parameters. See [bixverse::params_knn_defaults()]
#'   for available parameters and their defaults.
#'   \item pca - List of PCA parameters, see [bixverse::params_sc_pca()]
#'   for available parameters and their defaults.
#' }
#' @param use_precomputed_pca Boolean. Should the PCA in the object be used
#' if found. If you decide to do this, make sure that you have run the PCA
#' on the batch-aware HVG ideally.
#' @param seed Integer. Random seed.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns The object with the added fastMNN embeddings to the object.
#'
#' @export
#'
#' @references Haghverdi, et al., Nat Biotechnol, 2018
#'
#' @examples
#' # fastMNN over batch aware highly variable genes
#' sc <- demo_single_cells(
#'   syn_data_params = params_sc_synthetic_data(
#'     n_cells = 600L, n_genes = 50L, n_batches = 3L
#'   )
#' )
#' hvg <- find_hvg_batch_aware_sc(
#'   sc, hvg_no = 20L, batch_column = "batch_index", .verbose = FALSE
#' )
#' sc <- fast_mnn_sc(
#'   sc,
#'   batch_column = "batch_index",
#'   batch_hvg_genes = hvg$hvg_gene_idx,
#'   fastmnn_params = params_sc_fastmnn(
#'     no_pcs = 10L,
#'     knn = list(k = 5L)
#'   ),
#'   .verbose = FALSE
#' )
#' dim(get_embedding(sc, "mnn"))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
fast_mnn_sc <- S7::new_generic(
  name = "fast_mnn_sc",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    batch_hvg_genes,
    fastmnn_params = params_sc_fastmnn(),
    use_precomputed_pca = FALSE,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method fast_mnn_sc ScOrScSubset
S7::method(fast_mnn_sc, ScOrScSubset) <- function(
  object,
  batch_column,
  batch_hvg_genes,
  fastmnn_params = params_sc_fastmnn(),
  use_precomputed_pca = FALSE,
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(batch_hvg_genes, "I+")
  assertScFastmnn(fastmnn_params)
  checkmate::qassert(use_precomputed_pca, "B1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  batch_indices <- unlist(object[[batch_column]])
  batch_factor <- factor(batch_indices)
  batch_indices <- as.integer(batch_factor) - 1L

  # hard tier, but only on the branch that actually consumes the PCA
  if (use_precomputed_pca) {
    assert_sc_state(object, artefacts = "pca")
  }

  pca_data <- if (use_precomputed_pca && !is.null(get_pca_factors(object))) {
    if (.verbose) {
      message("Using pre-computed PCA found in the object")
    }
    get_pca_factors(object)
  } else {
    NULL
  }

  mnn_embd <- rs_mnn(
    f_path_gene = get_rust_count_gene_f_path(object),
    f_path_cell = get_rust_count_cell_f_path(object),
    cell_indices = get_cells_to_keep(object),
    gene_indices = as.integer(batch_hvg_genes),
    batch_indices = batch_indices,
    mnn_params = fastmnn_params,
    precomputed_pca = pca_data,
    verbose = parse_verbosity(.verbose),
    seed = 42L
  )

  colnames(mnn_embd) <- sprintf("mnn_%s", 1:ncol(mnn_embd))

  set_embedding(
    x = object,
    embd = mnn_embd,
    name = "mnn",
    from = if (is.null(pca_data)) character() else "pca"
  )
}

## harmony ---------------------------------------------------------------------

#' Run Harmony
#'
#' @description
#' A version of Harmony by Korsunsky et al., implemented in Rust. Performs
#' batch correction on PCA embeddings and stores the result as a `"harmony"`
#' embedding in the object.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param batch_column String. Column name in the object containing the primary
#' batch labels.
#' @param additional_batch_columns Optional character vector. Additional batch
#' columns to regress out. If `NULL`, only the primary batch column is used.
#' @param modality String. One of `c("rna", "adt")`. You can only use `"adt"`
#' on `SingleCellsMultiModal` class.
#' @param harmony_params List. Output of [bixverse::params_sc_harmony()].
#' @param seed Integer. For reproducibility.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns The object with a `"harmony"` embedding added. If no PCA embeddings
#' are found, returns the object unchanged with a warning.
#'
#' @export
#'
#' @examples
#' # Harmony correction of the PCA embedding
#' sc <- demo_single_cells(
#'   syn_data_params = params_sc_synthetic_data(
#'     n_cells = 600L, n_genes = 50L, n_batches = 3L
#'   )
#' )
#' sc <- harmony_sc(sc, batch_column = "batch_index", .verbose = FALSE)
#' dim(get_embedding(sc, "harmony"))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
harmony_sc <- S7::new_generic(
  name = "harmony_sc",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    additional_batch_columns = NULL,
    modality = c("rna", "adt"),
    harmony_params = params_sc_harmony(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method harmony_sc ScOrScSubset
S7::method(harmony_sc, ScOrScSubset) <- function(
  object,
  batch_column,
  additional_batch_columns = NULL,
  modality = c("rna", "adt"),
  harmony_params = params_sc_harmony(),
  seed = 42L,
  .verbose = TRUE
) {
  modality <- match.arg(modality)

  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(additional_batch_columns, c("S+", "0"))
  assertScHarmonyParams(harmony_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  if (modality != "rna" && !S7::S7_inherits(object, SingleCellsMultiModal)) {
    stop(sprintf(
      "modality = '%s' is only supported for SingleCellsMultiModal.",
      modality
    ))
  }

  # hard tier: the corrected embedding is written back onto the object, so a
  # stale PCA would silently produce a mis-aligned one
  assert_sc_state(object, artefacts = "pca", modality = modality)

  if (is.null(get_pca_factors(object, modality = modality))) {
    warning("No PCA embeddings found in the object. Returning class as is")
    return(object)
  }
  pca_data <- get_pca_factors(object, modality = modality)

  batch_index_ls <- list()

  batch_indices <- unlist(object[[batch_column]])
  batch_factor <- factor(batch_indices)
  batch_indices <- as.integer(batch_factor) - 1L

  batch_index_ls[[1]] <- batch_indices

  if (!is.null(additional_batch_columns)) {
    for (i in seq_along(additional_batch_columns)) {
      batch_indices_i <- unlist(object[[additional_batch_columns[[i]]]])
      batch_factor_i <- factor(batch_indices_i)
      batch_indices_i <- as.integer(batch_factor_i) - 1L
      batch_index_ls[[i + 1]] <- batch_indices_i
    }
  }

  checkmate::assertTRUE(all(
    purrr::map_dbl(batch_index_ls, length) == nrow(pca_data)
  ))

  if (is.null(harmony_params$k)) {
    harmony_params$k <- as.integer(min(round(nrow(pca_data) / 30), 100L))
    if (.verbose) {
      message(sprintf(
        " Auto-determined number of Harmony clusters: %d",
        harmony_params$k
      ))
    }
  }

  harmony_embd <- rs_harmony(
    pca = pca_data,
    harmony_params = harmony_params,
    batch_labels = batch_index_ls,
    seed = seed,
    verbose = parse_verbosity(.verbose)
  )

  colnames(harmony_embd) <- sprintf("harmony_%s", 1:ncol(harmony_embd))

  set_embedding(
    x = object,
    embd = harmony_embd,
    name = "harmony",
    modality = modality,
    from = "pca"
  )
}

## harmony v2 ------------------------------------------------------------------

#' Run Harmony v2
#'
#' @description
#' A version of Harmony v2 by Patikas et al., 2026, implemented in Rust.
#' Performs batch correction on PCA embeddings and stores the result as a
#' `"harmony_v2"` embedding in the object.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param batch_column String. Column name in the object containing the primary
#' batch labels.
#' @param additional_batch_columns Optional character vector. Additional batch
#' columns to regress out. If `NULL`, only the primary batch column is used.
#' @param modality String. One of `c("rna", "adt")`. You can only use `"adt"`
#' on `SingleCellsMultiModal` class.
#' @param harmony_params List. Output of [bixverse::params_sc_harmony_v2()].
#' @param seed Integer. For reproducibility.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns The object with a `"harmony_v2"` embedding added. If no PCA
#' embeddings are found, returns the object unchanged with a warning.
#'
#' @export
#'
#' @examples
#' # the reimplemented Harmony, writing its own embedding
#' sc <- demo_single_cells(
#'   syn_data_params = params_sc_synthetic_data(
#'     n_cells = 600L, n_genes = 50L, n_batches = 3L
#'   )
#' )
#' sc <- harmony_v2_sc(
#'   sc,
#'   batch_column = "batch_index",
#'   .verbose = FALSE
#' )
#' dim(get_embedding(sc, "harmony_v2"))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
harmony_v2_sc <- S7::new_generic(
  name = "harmony_v2_sc",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    additional_batch_columns = NULL,
    modality = c("rna", "adt"),
    harmony_params = params_sc_harmony_v2(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method harmony_v2_sc ScOrScSubset
S7::method(harmony_v2_sc, ScOrScSubset) <- function(
  object,
  batch_column,
  additional_batch_columns = NULL,
  modality = c("rna", "adt"),
  harmony_params = params_sc_harmony_v2(),
  seed = 42L,
  .verbose = TRUE
) {
  modality <- match.arg(modality)

  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(additional_batch_columns, c("S+", "0"))
  assertScHarmonyParamsV2(harmony_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  if (modality != "rna" && !S7::S7_inherits(object, SingleCellsMultiModal)) {
    stop(sprintf(
      "modality = '%s' is only supported for SingleCellsMultiModal.",
      modality
    ))
  }

  # hard tier: the corrected embedding is written back onto the object, so a
  # stale PCA would silently produce a mis-aligned one
  assert_sc_state(object, artefacts = "pca", modality = modality)

  if (is.null(get_pca_factors(object, modality = modality))) {
    warning("No PCA embeddings found in the object. Returning class as is")
    return(object)
  }
  pca_data <- get_pca_factors(object, modality = modality)

  batch_index_ls <- list()

  batch_indices <- unlist(object[[batch_column]])
  batch_factor <- factor(batch_indices)
  batch_indices <- as.integer(batch_factor) - 1L

  batch_index_ls[[1]] <- batch_indices

  if (!is.null(additional_batch_columns)) {
    for (i in seq_along(additional_batch_columns)) {
      batch_indices_i <- unlist(object[[additional_batch_columns[[i]]]])
      batch_factor_i <- factor(batch_indices_i)
      batch_indices_i <- as.integer(batch_factor_i) - 1L
      batch_index_ls[[i + 1]] <- batch_indices_i
    }
  }

  checkmate::assertTRUE(all(
    purrr::map_dbl(batch_index_ls, length) == nrow(pca_data)
  ))

  if (is.null(harmony_params$k)) {
    harmony_params$k <- as.integer(min(round(nrow(pca_data) / 30), 100L))
    if (.verbose) {
      message(sprintf(
        " Auto-determined number of Harmony clusters: %d",
        harmony_params$k
      ))
    }
  }

  harmony_embd <- rs_harmony_v2(
    pca = pca_data,
    harmony_params = harmony_params,
    batch_labels = batch_index_ls,
    seed = seed,
    verbose = parse_verbosity(.verbose)
  )

  colnames(harmony_embd) <- sprintf("harmony_v2_%s", 1:ncol(harmony_embd))

  set_embedding(
    x = object,
    embd = harmony_embd,
    name = "harmony_v2",
    modality = modality,
    from = "pca"
  )
}

## seurat CCA ------------------------------------------------------------------

#' Run Seurat CCA integration
#'
#' @description
#' This function implements the canonical correlation analysis (CCA) anchor
#' integration from Stuart, et al. For each pair of batches a shared CCA
#' embedding is computed, mutual nearest neighbours in that embedding become
#' anchors, these are filtered in gene space, scored by shared neighbours and
#' finally used to apply a kernel-weighted correction on the union PCA
#' embedding. Batches are merged in the order of their pairwise anchor counts.
#'
#' This port deviates from Seurat in two places. It skips the per-gene
#' `ScaleData` step and works from per-cell standardised log-normalised HVG
#' expression, and it never materialises the `N1 x N2` cross-product (the
#' canonical correlations come from a matrix-free randomised SVD). The
#' correction runs on the embedding, not on full log-expression. In practice
#' the anchor structure comes out close to identical at a fraction of the
#' memory.
#'
#' @param object `SingleCells` class.
#' @param batch_column String. The column with the batch information in the
#' obs data of the class.
#' @param batch_hvg_genes Integer vector. These are the highly variable genes,
#' identified by a batch-aware method. Please refer to
#' [bixverse::find_hvg_batch_aware_sc()] for more details. These genes have to
#' be 0-indexed!
#' @param cca_params A list, please see [bixverse::params_sc_seurat_cca()]. The
#' list has the following parameters:
#' \itemize{
#'   \item num_cc - Integer. Number of canonical correlation dimensions. The
#'   effective rank used is `max(num_cc, dims)`.
#'   \item dims - Integer. Dimensions used for the anchor kNN queries and size
#'   of the returned embedding.
#'   \item k_anchor - Integer. Neighbourhood size for the anchor search.
#'   \item k_filter - Integer. Neighbourhood size for the gene-space filter.
#'   \item k_score - Integer. Neighbourhood size for the anchor scoring.
#'   \item k_weight - Integer. Neighbourhood size for the kernel weights.
#'   \item n_top_features - Integer. Top-loading genes for the gene-space
#'   filter.
#'   \item l2_norm - Logical. L2-normalise the CCA embedding per cell.
#'   \item sd - Numeric. Bandwidth divisor of the Gaussian kernel.
#'   \item knn - List of kNN parameters. See [bixverse::params_knn_defaults()]
#'   for available parameters and their defaults.
#'   \item pca - List of PCA parameters, see [bixverse::params_sc_pca()]
#'   for available parameters and their defaults.
#' }
#' @param use_precomputed_pca Boolean. Should the PCA in the object be used
#' if found. If you decide to do this, make sure that you have run the PCA
#' on the batch-aware HVG ideally. Note that CCA still needs the PCA loadings
#' for the gene-space filter, so this saves less work than it does for fastMNN.
#' @param seed Integer. Random seed.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns The object with the added `"cca"` embedding.
#'
#' @export
#'
#' @references Stuart, et al., Cell, 2019
#'
#' @examples
#' # CCA anchor integration over batch aware highly variable genes
#' sc <- demo_single_cells(
#'   syn_data_params = params_sc_synthetic_data(
#'     n_cells = 600L, n_genes = 50L, n_batches = 3L
#'   )
#' )
#' hvg <- find_hvg_batch_aware_sc(
#'   sc, hvg_no = 20L, batch_column = "batch_index", .verbose = FALSE
#' )
#' sc <- seurat_cca_sc(
#'   sc,
#'   batch_column = "batch_index",
#'   batch_hvg_genes = hvg$hvg_gene_idx,
#'   cca_params = params_sc_seurat_cca(num_cc = 10L, dims = 10L),
#'   .verbose = FALSE
#' )
#' dim(get_embedding(sc, "cca"))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
seurat_cca_sc <- S7::new_generic(
  name = "seurat_cca_sc",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    batch_hvg_genes,
    cca_params = params_sc_seurat_cca(),
    use_precomputed_pca = FALSE,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method seurat_cca_sc SingleCells
#'
#' @export
S7::method(seurat_cca_sc, SingleCells) <- function(
  object,
  batch_column,
  batch_hvg_genes,
  cca_params = params_sc_seurat_cca(),
  use_precomputed_pca = FALSE,
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(batch_hvg_genes, "I+")
  assertScSeuratCca(cca_params)
  checkmate::qassert(use_precomputed_pca, "B1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # function body
  batch_indices <- unlist(object[[batch_column]])
  batch_factor <- factor(batch_indices)
  batch_indices <- as.integer(batch_factor) - 1L

  if (!length(levels(batch_factor)) > 1) {
    warning("The batch column only has one batch. Returning object as is.")
    return(object)
  }

  # hard tier, but only on the branch that actually consumes the PCA
  if (use_precomputed_pca) {
    assert_sc_state(object, artefacts = "pca")
  }

  pca_data <- if (use_precomputed_pca && !is.null(get_pca_factors(object))) {
    if (.verbose) {
      message("Using pre-computed PCA found in the object")
    }
    get_pca_factors(object)
  } else {
    NULL
  }

  cca_embd <- rs_seurat_cca(
    f_path_gene = get_rust_count_gene_f_path(object),
    f_path_cell = get_rust_count_cell_f_path(object),
    cell_indices = get_cells_to_keep(object),
    gene_indices = as.integer(batch_hvg_genes),
    batch_indices = batch_indices,
    precomputed_pca = pca_data,
    cca_params = cca_params,
    verbose = parse_verbosity(.verbose),
    seed = seed
  )

  colnames(cca_embd) <- sprintf("cca_%s", 1:ncol(cca_embd))

  object <- set_embedding(
    x = object,
    embd = cca_embd,
    name = "cca",
    from = if (is.null(pca_data)) character() else "pca"
  )

  return(object)
}

## seurat rPCA -----------------------------------------------------------------

#' Run Seurat rPCA integration
#'
#' @description
#' This function implements the reciprocal PCA (rPCA) anchor integration from
#' Stuart, et al. It runs the same anchor pipeline as
#' [bixverse::seurat_cca_sc()] but builds a cheaper per-pair anchor space: each
#' batch keeps its own PCA basis and the other batch's HVG expression is
#' projected into it. Cross-batch mutual nearest neighbours are then found in
#' these projected bases.
#'
#' rPCA is faster than CCA and corrects less aggressively, which makes it the
#' safer choice when batches share most of their cell types. As in Seurat, no
#' gene-space anchor filter is applied, that step is CCA-only.
#'
#' @param object `SingleCells` class.
#' @param batch_column String. The column with the batch information in the
#' obs data of the class.
#' @param batch_hvg_genes Integer vector. These are the highly variable genes,
#' identified by a batch-aware method. Please refer to
#' [bixverse::find_hvg_batch_aware_sc()] for more details. These genes have to
#' be 0-indexed!
#' @param rpca_params A list, please see [bixverse::params_sc_seurat_rpca()].
#' The list has the following parameters:
#' \itemize{
#'   \item dims - Integer. Dimensions used for the per-batch projections, the
#'   anchor kNN queries and the size of the returned embedding.
#'   \item k_anchor - Integer. Neighbourhood size for the anchor search.
#'   \item k_score - Integer. Neighbourhood size for the anchor scoring.
#'   \item k_weight - Integer. Neighbourhood size for the kernel weights.
#'   \item l2_norm - Logical. L2-normalise the projected embeddings per cell.
#'   \item sd - Numeric. Bandwidth divisor of the Gaussian kernel.
#'   \item knn - List of kNN parameters. See [bixverse::params_knn_defaults()]
#'   for available parameters and their defaults.
#'   \item pca - List of PCA parameters, see [bixverse::params_sc_pca()]
#'   for available parameters and their defaults.
#' }
#' @param use_precomputed_pca Boolean. Should the PCA in the object be used
#' if found. This only applies to the union PCA that gets corrected, the
#' per-batch PCAs are always recomputed because rPCA needs them.
#' @param seed Integer. Random seed.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns The object with the added `"rpca"` embedding.
#'
#' @export
#'
#' @references Stuart, et al., Cell, 2019
#'
#' @examples
#' # reciprocal PCA anchor integration, the cheaper sibling of CCA
#' sc <- demo_single_cells(
#'   syn_data_params = params_sc_synthetic_data(
#'     n_cells = 600L, n_genes = 50L, n_batches = 3L
#'   )
#' )
#' hvg <- find_hvg_batch_aware_sc(
#'   sc, hvg_no = 20L, batch_column = "batch_index", .verbose = FALSE
#' )
#' sc <- seurat_rpca_sc(
#'   sc,
#'   batch_column = "batch_index",
#'   batch_hvg_genes = hvg$hvg_gene_idx,
#'   rpca_params = params_sc_seurat_rpca(dims = 10L),
#'   .verbose = FALSE
#' )
#' dim(get_embedding(sc, "rpca"))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
seurat_rpca_sc <- S7::new_generic(
  name = "seurat_rpca_sc",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    batch_hvg_genes,
    rpca_params = params_sc_seurat_rpca(),
    use_precomputed_pca = FALSE,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method seurat_rpca_sc SingleCells
#'
#' @export
S7::method(seurat_rpca_sc, SingleCells) <- function(
  object,
  batch_column,
  batch_hvg_genes,
  rpca_params = params_sc_seurat_rpca(),
  use_precomputed_pca = FALSE,
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(batch_hvg_genes, "I+")
  assertScSeuratRpca(rpca_params)
  checkmate::qassert(use_precomputed_pca, "B1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # function body
  batch_indices <- unlist(object[[batch_column]])
  batch_factor <- factor(batch_indices)
  batch_indices <- as.integer(batch_factor) - 1L

  if (!length(levels(batch_factor)) > 1) {
    warning("The batch column only has one batch. Returning object as is.")
    return(object)
  }

  # hard tier, but only on the branch that actually consumes the PCA
  if (use_precomputed_pca) {
    assert_sc_state(object, artefacts = "pca")
  }

  pca_data <- if (use_precomputed_pca && !is.null(get_pca_factors(object))) {
    if (.verbose) {
      message("Using pre-computed PCA found in the object")
    }
    get_pca_factors(object)
  } else {
    NULL
  }

  rpca_embd <- rs_seurat_rpca(
    f_path_gene = get_rust_count_gene_f_path(object),
    f_path_cell = get_rust_count_cell_f_path(object),
    cell_indices = get_cells_to_keep(object),
    gene_indices = as.integer(batch_hvg_genes),
    batch_indices = batch_indices,
    precomputed_pca = pca_data,
    rpca_params = rpca_params,
    verbose = parse_verbosity(.verbose),
    seed = seed
  )

  colnames(rpca_embd) <- sprintf("rpca_%s", 1:ncol(rpca_embd))

  object <- set_embedding(
    x = object,
    embd = rpca_embd,
    name = "rpca",
    from = if (is.null(pca_data)) character() else "pca"
  )

  return(object)
}
