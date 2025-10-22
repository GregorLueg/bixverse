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
#' @param object `single_cell_exp` class.
#' @param batch_column String. The column with the batch information in the
#' obs data of the class.
#' @param threshold Numeric. Number between 0 and 1. Below this threshold, the
#' test is considered significant. Defaults to `0.05`.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @returns A list with the following elements
#' \itemize{
#'   \item kbet_score - Number of significant tests over all cells. 0 indicates
#'   perfect mixing, 1 indicates basically zero mixing between batches.
#'   \item significant_tests - Boolean indicating for which cells the statistic
#'   was below the threshold
#'   \item chisquare_pvals - The p-values of the ChiSquare test.
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


#' @method calculate_kbet_sc single_cell_exp
#'
#' @export
S7::method(calculate_kbet_sc, single_cell_exp) <- function(
  object,
  batch_column,
  threshold = 0.05,
  .verbose = TRUE
) {
  # check
  checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))
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

  rs_res <- rs_kbet(
    knn_mat = knn_mat,
    batch_vector = as.integer(factor(batch_index))
  )

  res <- list(
    kbet_score = sum(rs_res <= threshold) / length(rs_res),
    significant_tests = rs_res <= threshold,
    chisquare_pvals = rs_res
  )

  return(res)
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
#' @param object `single_cell_exp` class.
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

#' @method find_hvg_batch_aware_sc single_cell_exp
#'
#' @export
S7::method(find_hvg_batch_aware_sc, single_cell_exp) <- function(
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
  checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))
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
      f_path_gene = get_rust_count_gene_f_path(sc_object.weak_batch_effect),
      hvg_method = method,
      cell_indices = get_cells_to_keep(sc_object.weak_batch_effect),
      batch_labels = batch_index,
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
      batch_hvgs_dt[, .SD[order(-var_std)][1:no_hvg], by = batch][, unique(
        gene_idx
      )]
    },
    average = {
      avg_dt <- batch_hvgs_dt[, .(var_std_avg = mean(var_std)), by = gene_idx]
      avg_dt[order(-var_std_avg)][1:no_hvg, gene_idx]
    },
    intersection = {
      top_per_batch <- batch_hvgs_dt[,
        .(gene_idx = .SD[order(-var_std)][1:no_hvg, gene_idx]),
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
#' from Polański, et al. Briefly, the algorithm generate a KNN index on a per
#' batch basis and identifies the neighbours of cells for each individual index.
#' Subsequently, it leverages UMAP connectivity calculations to reduce spurious
#' connections. For more details, please refer to Polański, et al.
#'
#' @param object `single_cell_exp` class.
#' @param batch_column String. The column with the batch information in the
#' obs data of the class.
#' @param no_neighbours_to_keep Integer. Maximum number of neighbours to keep
#' from the BBKNN algorithm. Due to generating neighbours for each batch, there
#' might be a large number of generated neighbours. This will only keep the top
#' `no_neighbours_to_keep` neighbours.
#' @param embd_to_use String. The embedding to use. Atm, the only option is
#' `"pca"`.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param sc_bbknn_params A list, please see [bixverse::params_sc_bbknn()]. The
#' list has the following parameters:
#' \itemize{
#'   \item neighbours_within_batch - Integer. Number of neighbours to consider
#'   per batch.
#'   \item knn_method - String. One of `c("annoy", "hnsw")`. Defaults to
#'   `"annoy"`.
#'   \item ann_dist - String. One of `c("cosine", "euclidean")`. The distance
#'   metric to be used for the approximate neighbour search. Defaults to
#'   `"cosine"`.
#'   \item set_op_mix_ratio - Numeric. Mixing ratio between union (1.0) and
#'   intersection (0.0).
#'   \item local_connectivity - Numeric. UMAP connectivity computation
#'   parameter, how many nearest neighbours of each cell are assumed to be fully
#'   connected.
#'   \item annoy_n_trees - Integer. Number of trees to use in the generation of
#'   the Annoy index.
#'   \item search_budget - Integer. Search budget per tree for the `annoy`
#'   algorithm.
#'   \item trim - Optional integer. Trim the neighbours of each cell to these
#'   many top connectivities. May help with population independence and improve
#'   the tidiness of clustering. If `NULL`, it defaults to
#'   `10 * neighbours_within_batch`.
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
    no_neighbours_to_keep = 15L,
    embd_to_use = "pca",
    no_embd_to_use = NULL,
    bbknn_params = params_sc_bbknn(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method bbknn_sc single_cell_exp
#'
#' @export
S7::method(bbknn_sc, single_cell_exp) <- function(
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
  checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))
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
  # early return
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

  knn_mat <- if (no_generated_neighbours <= no_neighbours_to_keep) {
    matrix(data = bbknn_res$distances$indices, nrow = nrow(embd), byrow = TRUE)
  } else {
    rs_bbknn_filtering(
      indptr = bbknn_res$distances$indptr,
      indices = bbknn_res$distances$indices,
      no_neighbours_to_keep = no_neighbours_to_keep
    )
  }

  storage.mode(knn_mat) <- "integer"

  object <- set_knn(object, knn_mat = knn_mat)

  if (.verbose) {
    message(paste(
      "Generating graph based on BBKNN connectivities.",
      "Weights will be based on the connectivities and not shared nearest neighour",
      "calculations."
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
#' an embedding only and not fully corrected count matrix. The function will
#' iterate through the batches, identify the MNN and generate correction vectors
#' and generate a corrected embedding which is added to the function.
#'
#' @param object `single_cell_exp` class.
#' @param batch_column String. The column with the batch information in the
#' obs data of the class.
#' @param batch_hvg_genes Integer vector. These are the highly variable genes,
#' identified by a batch-aware method. Please refer to
#' [bixverse::find_hvg_batch_aware_sc()] for more details.
#' @param fastmnn_params A list, please see [bixverse::params_sc_fastmnn()]. The
#' list has the following parameters:
#' Claude fill this out
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
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method fast_mnn_sc single_cell_exp
#'
#' @export
S7::method(fast_mnn_sc, single_cell_exp) <- function(
  object,
  batch_column,
  batch_hvg_genes,
  fastmnn_params = params_sc_fastmnn(),
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))
  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(batch_hvg_genes, "I+")
  assertScFastmnn(fastmnn_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  # function body
  batch_indices <- unlist(object[[batch_column]])
  batch_factor <- factor(batch_indices)
  batch_indices <- as.integer(batch_factor) - 1L

  mnn_embd <- rs_mnn(
    f_path = get_rust_count_gene_f_path(object),
    cell_indices = get_cells_to_keep(object),
    gene_indices = as.integer(batch_hvg_genes - 1L),
    batch_indices = batch_indices,
    mnn_params = fastmnn_params,
    verbose = TRUE,
    seed = 42L
  )

  object <- set_embedding(x = object, embd = mnn_embd, name = "mnn")

  return(object)
}
