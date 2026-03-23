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
#' @param object `SingleCells` class.
#' @param batch_column String. The column with the batch information in the
#' obs data of the class.
#' @param threshold Numeric. Number between 0 and 1. Below this threshold, the
#' test is considered significant. Defaults to `0.05`.
#' @param .verbose Boolean. Controls the verbosity of the function.
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

#' @method calculate_kbet_sc SingleCells
#'
#' @export
S7::method(calculate_kbet_sc, SingleCells) <- function(
  object,
  batch_column,
  threshold = 0.05,
  .verbose = TRUE
) {
  # check
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(threshold, "N1[0, 1]")
  checkmate::qassert(.verbose, "B1")

  # function
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
    batch_vector = as.integer(factor(batch_index))
  )

  res <- structure(
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

  return(res)
}

### kbet print -----------------------------------------------------------------

#' Print method for KbetScores
#'
#' @param x A `KbetScores` object.
#' @param ... Additional arguments (ignored).
#'
#' @export
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
#' [bixverse::calculate_batch_lisi_sc()] instead.
#'
#' @param object `SingleCells` class.
#' @param batch_column String. The column with the batch information in the
#' obs data of the class.
#' @param embd_to_use String. Which embedding to compute the ASW on. One of
#' `c("pca", "harmony", "mnn")`. Defaults to `"pca"`.
#' @param max_cells Integer or `NULL`. If not `NULL`, subsample to this many
#' cells for performance. The pairwise distance computation is O(n^2), so
#' subsampling is recommended for large datasets. Defaults to `5000L`.
#' @param seed Integer. Seed for subsampling reproducibility.
#' @param .verbose Boolean. Controls the verbosity of the function.
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

#' @method calculate_batch_asw_sc SingleCells
#'
#' @export
S7::method(calculate_batch_asw_sc, SingleCells) <- function(
  object,
  batch_column,
  embd_to_use = "pca",
  max_cells = 5000L,
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(batch_column, "S1")
  checkmate::assertChoice(embd_to_use, c("pca", "harmony", "mnn"))
  checkmate::qassert(max_cells, c("I1", "0"))
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  batch_index <- unlist(object[[batch_column]])

  if (!length(levels(factor(batch_index))) > 1) {
    warning("The batch column only has one batch. Returning NULL")
    return(NULL)
  }

  embd <- switch(
    embd_to_use,
    pca = get_pca_factors(object),
    harmony = get_embedding(object, "harmony"),
    mnn = get_embedding(object, "mnn")
  )

  if (is.null(embd)) {
    warning("Embedding not found. Returning NULL")
    return(NULL)
  }

  n_batches <- length(levels(factor(batch_index)))

  rs_res <- rs_batch_silhouette_width(
    embedding = embd,
    batch_vector = as.integer(factor(batch_index)),
    max_cells = max_cells,
    seed = seed
  )

  res <- structure(
    list(
      per_cell = rs_res$per_cell,
      mean_asw = rs_res$mean_asw,
      median_asw = rs_res$median_asw,
      n_batches = n_batches,
      embedding_used = embd_to_use
    ),
    class = "BatchSilhouetteScores"
  )

  return(res)
}

### batch silhouette print -----------------------------------------------------

#' Print method for BatchSilhouetteScores
#'
#' @param x A `BatchSilhouetteScores` object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.BatchSilhouetteScores <- function(x, ...) {
  n_cells <- length(x$per_cell)

  cat("Batch Silhouette Width\n")
  cat(sprintf("  Cells: %d | Batches: %d\n", n_cells, x$n_batches))
  cat(sprintf(
    "  Mean ASW:    %.4f (0 = perfect mixing, 1 = separated)\n",
    x$mean_asw
  ))
  cat(sprintf("  Median ASW:  %.4f\n", x$median_asw))

  invisible(x)
}

## batch LISI ------------------------------------------------------------------

#' Calculate batch LISI scores
#'
#' @description
#' Computes the Local Inverse Simpson's Index (LISI) on batch labels using the
#' kNN graph. LISI measures the effective number of batches represented in each
#' cell's neighbourhood. Under perfect mixing, LISI equals the number of
#' batches. Under no mixing, LISI equals 1. Unlike kBET, LISI does not compare
#' against global batch proportions, making it suitable for graph-based
#' correction methods like BBKNN.
#'
#' @param object `SingleCells` class.
#' @param batch_column String. The column with the batch information in the
#' obs data of the class.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @returns A `BatchLisiScores` object with the following elements
#' \itemize{
#'   \item per_cell - Per-cell LISI scores in `[1, n_batches]`.
#'   \item mean_lisi - Mean LISI across all cells.
#'   \item median_lisi - Median LISI across all cells.
#'   \item n_batches - Number of batches in the data.
#' }
#'
#' @export
#'
#' @references Korsunsky, et al., Nat. Methods, 2019
calculate_batch_lisi_sc <- S7::new_generic(
  name = "calculate_batch_lisi_sc",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calculate_batch_lisi_sc SingleCells
#'
#' @export
S7::method(calculate_batch_lisi_sc, SingleCells) <- function(
  object,
  batch_column,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(.verbose, "B1")

  batch_index <- unlist(object[[batch_column]])

  if (!length(levels(factor(batch_index))) > 1) {
    warning("The batch column only has one batch. Returning NULL")
    return(NULL)
  }

  knn_mat <- get_knn_mat(object)

  if (is.null(knn_mat)) {
    warning("No kNN matrix found to calculate batch LISI. Returning NULL")
    return(NULL)
  }

  n_batches <- length(levels(factor(batch_index)))

  rs_res <- rs_batch_lisi(
    knn_mat = knn_mat,
    batch_vector = as.integer(factor(batch_index))
  )

  res <- structure(
    list(
      per_cell = rs_res$per_cell,
      mean_lisi = rs_res$mean_lisi,
      median_lisi = rs_res$median_lisi,
      n_batches = n_batches
    ),
    class = "BatchLisiScores"
  )

  return(res)
}

### batch lisi print -----------------------------------------------------------

#' Print method for BatchLisiScores
#'
#' @param x A `BatchLisiScores` object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.BatchLisiScores <- function(x, ...) {
  n_cells <- length(x$per_cell)

  cat("Batch LISI Scores\n")
  cat(sprintf("  Cells: %d | Batches: %d\n", n_cells, x$n_batches))
  cat(sprintf(
    "  Mean LISI:    %.4f (1 = no mixing, %d = perfect mixing)\n",
    x$mean_lisi,
    x$n_batches
  ))
  cat(sprintf("  Median LISI:  %.4f\n", x$median_lisi))

  invisible(x)
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
#' amongst the Top X HVG in all batches.
#'
#' @param object `SingleCells` class.
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
#' @param streaming Boolean. Shall the genes be streamed in. Useful for larger
#' data sets where you wish to avoid loading in the whole data. Defaults to
#' `FALSE`.
#' @param .verbose Boolean. Controls verbosity and returns run times.
#'
#' @return This function will return a list with:
#' \itemize{
#'   \item hvg_indices - The indices of the batch-aware HVG.
#'   \item batch_hvg_data - data.table with the detailed information of the
#'   variance per batch.
#' }
#'
#' @export
find_hvg_batch_aware_sc <- S7::new_generic(
  name = "find_hvg_batch_aware_sc",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    hvg_no = 2000L,
    gene_comb_method = c("union", "average", "intersection"),
    hvg_params = params_sc_hvg(),
    streaming = FALSE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method find_hvg_batch_aware_sc SingleCells
#'
#' @export
S7::method(find_hvg_batch_aware_sc, SingleCells) <- function(
  object,
  batch_column,
  hvg_no = 2000L,
  gene_comb_method = c("union", "average", "intersection"),
  hvg_params = params_sc_hvg(),
  streaming = FALSE,
  .verbose = TRUE
) {
  gene_comb_method <- match.arg(gene_comb_method)

  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(hvg_no, "I1")
  checkmate::assertChoice(
    gene_comb_method,
    c("union", "average", "intersection")
  )
  assertScHvg(hvg_params)
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(.verbose, "B1")

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
      streaming = streaming,
      verbose = .verbose
    )
  )

  batch_hvgs_dt <- data.table::as.data.table(batch_hvgs)
  batch_hvgs_dt[, batch := levels(batch_factor)[batch + 1L]]

  hvg_genes <- switch(
    gene_comb_method,
    union = {
      batch_hvgs_dt[, .SD[order(-var_std)][1:hvg_no], by = batch][, unique(
        gene_idx
      )]
    },
    average = {
      avg_dt <- batch_hvgs_dt[, .(var_std_avg = mean(var_std)), by = gene_idx]
      avg_dt[order(-var_std_avg)][1:hvg_no, gene_idx]
    },
    intersection = {
      top_per_batch <- batch_hvgs_dt[,
        .(gene_idx = .SD[order(-var_std)][1:hvg_no, gene_idx]),
        by = batch
      ]
      top_per_batch[, .N, by = gene_idx][
        N == uniqueN(batch_hvgs_dt$batch),
        gene_idx
      ]
    }
  )

  return(list(hvg_genes = hvg_genes, hvg_data = batch_hvgs_dt))
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
#' @param object `SingleCells` class.
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
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @returns The object with added kNN matrix based on BBKNN and the graph based
#' on the returned connectivities of the algorithm.
#'
#' @export
#'
#' @references Polański, et al., Bioinformatics, 2020
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

#' @method bbknn_sc SingleCells
#'
#' @export
S7::method(bbknn_sc, SingleCells) <- function(
  object,
  batch_column,
  no_neighbours_to_keep = 15L,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  bbknn_params = params_sc_bbknn(),
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(batch_column, "S1")
  checkmate::assertChoice(embd_to_use, c("pca"))
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  assertScBbknn(bbknn_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  # function body
  if (!is.null(get_knn_mat(object))) {
    warning("Prior kNN matrix found. Will be overwritten.")
  }

  embd <- switch(embd_to_use, pca = get_pca_factors(object))

  if (is.null(embd)) {
    warning(
      paste(
        "The desired embedding was not found. Please check the parameters.",
        "Returning NULL."
      )
    )
    return(NULL)
  }

  if (!is.null(no_embd_to_use)) {
    to_take <- min(c(no_embd_to_use, ncol(embd)))
    embd <- embd[, 1:to_take]
  }

  batch_index <- unlist(object[[batch_column]])

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
    verbose = .verbose
  )

  # extract kNN indices and distances
  # In bbknn_sc, replace the if/else with:
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
  object <- set_knn(object, knn = sc_knn)

  # build graph from BBKNN connectivities
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

  object <- set_snn_graph(object, snn_graph = snn_graph)

  return(object)
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
#' @param object `SingleCells` class.
#' @param batch_column String. The column with the batch information in the
#' obs data of the class.
#' @param batch_hvg_genes Integer vector. These are the highly variable genes,
#' identified by a batch-aware method. Please refer to
#' [bixverse::find_hvg_batch_aware_sc()] for more details.
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
#'   \item random_svd - Logical. Use randomised SVD.
#'   \item knn - List of kNN parameters. See [bixverse::params_knn_defaults()]
#'   for available parameters and their defaults.
#' }
#' @param use_precomputed_pca Boolean. Should the PCA in the object be used
#' if found. If you decide to do this, make sure that you have run the PCA
#' on the batch-aware HVG ideally.
#' @param seed Integer. Random seed.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @returns The object with the added fastMNN embeddings to the object.
#'
#' @export
#'
#' @references Haghverdi, et al., Nat Biotechnol, 2018
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

#' @method fast_mnn_sc SingleCells
#'
#' @export
S7::method(fast_mnn_sc, SingleCells) <- function(
  object,
  batch_column,
  batch_hvg_genes,
  fastmnn_params = params_sc_fastmnn(),
  use_precomputed_pca = FALSE,
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(batch_hvg_genes, "I+")
  assertScFastmnn(fastmnn_params)
  checkmate::qassert(use_precomputed_pca, "B1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  # function body
  batch_indices <- unlist(object[[batch_column]])
  batch_factor <- factor(batch_indices)
  batch_indices <- as.integer(batch_factor) - 1L

  pca_data <- if (use_precomputed_pca && !is.null(get_pca_factors(object))) {
    if (.verbose) {
      message("Using pre-computed PCA found in the object")
    }
    get_pca_factors(object)
  } else {
    NULL
  }

  mnn_embd <- rs_mnn(
    f_path = get_rust_count_gene_f_path(object),
    cell_indices = get_cells_to_keep(object),
    gene_indices = as.integer(batch_hvg_genes - 1L),
    batch_indices = batch_indices,
    mnn_params = fastmnn_params,
    precomputed_pca = pca_data,
    verbose = .verbose,
    seed = 42L
  )

  colnames(mnn_embd) <- sprintf("mnn_%s", 1:ncol(mnn_embd))

  object <- set_embedding(x = object, embd = mnn_embd, name = "mnn")

  return(object)
}

## harmony ---------------------------------------------------------------------

#' Run Harmony
#'
#' @description
#' A version of Harmony by Korsunsky et al., implemented in Rust. Performs
#' batch correction on PCA embeddings and stores the result as a `"harmony"`
#' embedding in the object.
#'
#' @param object `SingleCells` class.
#' @param batch_column String. Column name in the object containing the primary
#' batch labels.
#' @param additional_batch_columns Optional character vector. Additional batch
#' columns to regress out. If `NULL`, only the primary batch column is used.
#' @param harmony_params List. Output of [bixverse::params_sc_harmony()].
#' @param seed Integer. For reproducibility.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @return The object with a `"harmony"` embedding added. If no PCA embeddings
#' are found, returns the object unchanged with a warning.
#'
#' @export
harmony_sc <- S7::new_generic(
  name = "harmony_sc",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    additional_batch_columns = NULL,
    harmony_params = params_sc_harmony(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method harmony_sc SingleCells
#'
#' @export
S7::method(harmony_sc, SingleCells) <- function(
  object,
  batch_column,
  additional_batch_columns = NULL,
  harmony_params = params_sc_harmony(),
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(additional_batch_columns, c("S+", "0"))
  assertScHarmonyParams(harmony_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  # early return
  if (is.null(get_pca_factors(object))) {
    warning(paste(
      "No PCA embeddings found in the object. Returning class as is"
    ))
    return(object)
  } else {
    pca_data <- get_pca_factors(object)
  }

  # function body
  # main batch
  batch_index_ls <- list()

  batch_indices <- unlist(object[[batch_column]])
  batch_factor <- factor(batch_indices)
  batch_indices <- as.integer(batch_factor) - 1L

  batch_index_ls[[1]] <- batch_indices

  # add optional batch effects to regress out
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
    harmony_params$k <- min(round(nrow(pca_data) / 30), 200L)
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
    verbose = .verbose
  )

  colnames(harmony_embd) <- sprintf("harmony_%s", 1:ncol(harmony_embd))

  object <- set_embedding(x = object, embd = harmony_embd, name = "harmony")

  return(object)
}
