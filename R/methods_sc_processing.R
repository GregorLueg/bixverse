# single cell processing methods -----------------------------------------------

## doublet detection -----------------------------------------------------------

### scrublet -------------------------------------------------------------------

#' Doublet detection with Scrublet
#'
#' @description This function implements the doublet detection from Scrublet,
#' see Wolock, et al. Briefly, arteficial doublets are being generated from
#' the data via random combination of initial cells. Highly variable genes (HVG)
#' are being identified and the observed cells are being projected on a PCA
#' space. Subsequently, the simulated doublets are being projected on the same
#' PCA space given the same HVGs; kNN graphs are being generated and a kNN
#' classifier is used to assign a probability that a given cell in the original
#' data is a doublet. For more details, please check the publication.
#'
#' @param object `SingleCells` class.
#' @param scrublet_params A list with the final scrublet parameters, see
#' [bixverse::params_scrublet()] for full details.
#' @param seed Integer. Random seed.
#' @param streaming Boolean. Shall streaming be used during the HVG
#' calculations. Slower, but less memory usage.
#' @param return_combined_pca Boolean. Shall the PCA of the observed cells and
#' simulated doublets be returned.
#' @param return_pairs Boolean. Shall the pairs be returned.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return A `scrublet_res` class that has with the following items:
#' \itemize{
#'   \item predicted_doublets - Boolean vector indicating which observed cells
#'   predicted as doublets (TRUE = doublet, FALSE = singlet).
#'   \item doublet_scores_obs - Numerical vector with the likelihood of being
#'   a doublet for the observed cells.
#'   \item doublet_scores_sim - Numerical vector with the likelihood of being
#'   a doublet for the simulated cells.
#'   \item doublet_errors_obs - Numerical vector with the standard errors of
#'   the scores for the observed cells.
#'   \item z_scores - Z-scores for the observed cells. Represents:
#'   `score - threshold / error`.
#'   \item threshold - Used threshold.
#'   \item detected_doublet_rate - Fraction of cells that are called as
#'   doublet.
#'   \item detectable_doublet_fraction - Fraction of simulated doublets with
#'   scores above the threshold.
#'   \item overall_doublet_rate - Estimated overall doublet rate. Should roughly
#'   match the expected doublet rate.
#'   \item pca - Optional PCA embeddings across the original cells and simulated
#'   doublets.
#'   \item pair_1 - Optional index of the parent cell 1 of the simulated
#'   doublets.
#'   \item pair_2 - Optional index of the parent cell 2 of the simulated
#'   doublets.
#' }
#'
#' @export
#'
#' @references Wollock, et al., Cell Syst, 2020
scrublet_sc <- S7::new_generic(
  name = "scrublet_sc",
  dispatch_args = "object",
  fun = function(
    object,
    scrublet_params = params_scrublet(),
    seed = 42L,
    streaming = FALSE,
    return_combined_pca = FALSE,
    return_pairs = FALSE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method scrublet_sc SingleCells
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(scrublet_sc, SingleCells) <- function(
  object,
  scrublet_params = params_scrublet(),
  seed = 42L,
  streaming = FALSE,
  return_combined_pca = FALSE,
  return_pairs = FALSE,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScScrublet(scrublet_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(return_combined_pca, "B1")
  checkmate::qassert(return_pairs, "B1")
  checkmate::qassert(.verbose, "B1")

  # function body
  cells_to_keep <- get_cells_to_keep(object)

  if (length(cells_to_keep) >= 100000) {
    message("Setting PCA to sparse default. N_cells greater than 100,000")

    scrublet_params$sparse <- TRUE
  }

  scrublet_res <- rs_sc_scrublet(
    f_path_gene = get_rust_count_gene_f_path(object),
    f_path_cell = get_rust_count_cell_f_path(object),
    cells_to_keep = cells_to_keep,
    scrublet_params = scrublet_params,
    seed = seed,
    verbose = .verbose,
    streaming = streaming,
    return_combined_pca = return_combined_pca,
    return_pairs = return_pairs
  )

  attr(scrublet_res, "cell_indices") <- cells_to_keep
  class(scrublet_res) <- "ScrubletRes"

  return(scrublet_res)
}

### boost ----------------------------------------------------------------------

#' Doublet detection with boosted doublet classification
#'
#' @description This function implements the boosted doublet detection. It
#' generates through several iterations simulated doublets, generate kNN graphs,
#' runs Louvain clustering and assesses how often an observed cells clsuters
#' together with the simulated doublets.
#'
#' @param object `SingleCells` class.
#' @param boost_params A list with the final scrublet parameters, see
#' [bixverse::params_boost()] for full details.
#' @param seed Integer. Random seed.
#' @param streaming Boolean. Shall streaming be used during the HVG
#' calculations. Slower, but less memory usage.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return A `boost_res` class that has with the following items:
#' \itemize{
#'   \item predicted_doublets - Boolean vector indicating which observed cells
#'   predicted as doublets (TRUE = doublet, FALSE = singlet).
#'   \item doublet_scores_obs - Numerical vector with the likelihood of being
#'   a doublet for the observed cells.
#'   \item voting_avg - Numerical vector with the average voting score.
#' }
#'
#' @export
doublet_detection_boost_sc <- S7::new_generic(
  name = "doublet_detection_boost_sc",
  dispatch_args = "object",
  fun = function(
    object,
    boost_params = params_boost(),
    seed = 42L,
    streaming = FALSE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method doublet_detection_boost_sc SingleCells
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(doublet_detection_boost_sc, SingleCells) <- function(
  object,
  boost_params = params_boost(),
  seed = 42L,
  streaming = FALSE,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScBoost(boost_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(.verbose, "B1")

  # function body
  cells_to_keep <- get_cells_to_keep(object)

  if (length(cells_to_keep) >= 100000) {
    message("Setting PCA to sparse default. N_cells greater than 100,000")

    boost_params$sparse <- TRUE
  }

  boost_res <- rs_sc_doublet_detection(
    f_path_gene = get_rust_count_gene_f_path(object),
    f_path_cell = get_rust_count_cell_f_path(object),
    cells_to_keep = cells_to_keep,
    boost_params = boost_params,
    seed = seed,
    verbose = .verbose,
    streaming = streaming
  )

  attr(boost_res, "cell_indices") <- cells_to_keep
  class(boost_res) <- "BoostRes"

  return(boost_res)
}

### scdblfinder ----------------------------------------------------------------

#' Run scDblFinder doublet detection on a SingleCells object
#'
#' @description
#' Cluster-aware doublet detection using engineered features and a
#' gradient-boosted classifier. See Germain et al., F1000Research, 2022.
#'
#' @param object `SingleCells` class.
#' @param scdblfinder_params List. Parameters from
#' [bixverse::params_scdblfinder()].
#' @param return_features Boolean. Shall the features used to train the
#' classifier be returned.
#' @param seed Integer. Seed for reproducibility.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @return An S3 object of class `ScDblFinderRes` containing:
#' \describe{
#'   \item{predicted_doublets}{Logical vector of doublet calls.}
#'   \item{doublet_score}{Numeric vector of classifier probabilities.}
#'   \item{cxds_scores}{Numeric vector of the cxds scores.}
#'   \item{weighted}{Numeric vector of the weighted scores.}
#'   \item{threshold}{The threshold used for calling.}
#'   \item{cluster_labels}{Integer vector of final cluster assignments.}
#'   \item{detected_doublet_rate}{Fraction of cells called as doublets.}
#' }
#' with `cell_indices` stored as an attribute.
#'
#' @export
scdblfinder_sc <- S7::new_generic(
  name = "scdblfinder_sc",
  dispatch_args = "object",
  fun = function(
    object,
    scdblfinder_params = params_scdblfinder(),
    return_features = FALSE,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method scdblfinder_sc SingleCells
#'
#' @export
S7::method(scdblfinder_sc, SingleCells) <- function(
  object,
  scdblfinder_params = params_scdblfinder(),
  return_features = FALSE,
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScDblFinder(scdblfinder_params)
  checkmate::qassert(return_features, "B1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  cell_indices <- get_cells_to_keep(object)

  if (length(cell_indices) >= 100000) {
    message("Setting PCA to sparse default. N_cells greater than 100,000")

    scdblfinder_params$sparse <- TRUE
  }

  res <- rs_sc_scdblfinder(
    f_path_gene = get_rust_count_gene_f_path(object),
    f_path_cell = get_rust_count_cell_f_path(object),
    cell_indices = cell_indices,
    params = scdblfinder_params,
    return_features = return_features,
    seed = seed,
    verbose = .verbose,
    debug = FALSE
  )

  features <- if (return_features) {
    feature_matrix <- res$features$feature_mat
    colnames(feature_matrix) <- res$features$feature_names
    rownames(feature_matrix) <- get_cell_names(object, filtered = TRUE)

    feature_matrix
  } else {
    NULL
  }

  structure(
    list(
      predicted_doublets = res$predicted_doublets,
      doublet_score = res$doublet_scores,
      cxds_scores = res$cxds_scores,
      weighted = res$weighted,
      threshold = res$threshold,
      cluster_labels = res$cluster_labels,
      detected_doublet_rate = res$detected_doublet_rate,
      features = features
    ),
    cell_indices = cell_indices,
    class = "ScDblFinderRes"
  )
}

## gene proportions ------------------------------------------------------------

### top n genes ----------------------------------------------------------------

#' Calculate the proportions of reads for the Top N genes
#'
#' @description
#' This is a helper function that calculates proportions of reads to the Top N
#' genes by expression in a given cell. High values here can indicate low
#' complexity, quality cells. The values will be automatically added to the
#' obs table.
#'
#' @param object `SingleCells` class.
#' @param top_n_vals Integer. The Top N thresholds to test.
#' @param streaming Boolean. Shall the cells be streamed in. Useful for larger
#' data sets where you wish to avoid loading in the whole data. Default to
#' `FALSE`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return It will add the columns based on the names in the `gene_set_list` to
#' the obs table.
#'
#' @export
top_genes_perc_sc <- S7::new_generic(
  name = "top_genes_perc_sc",
  dispatch_args = "object",
  fun = function(
    object,
    top_n_vals = c(25L, 50L, 100L),
    streaming = FALSE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method top_genes_perc_sc SingleCells
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(top_genes_perc_sc, SingleCells) <- function(
  object,
  top_n_vals = c(25L, 50L, 100L),
  streaming = FALSE,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(top_n_vals, "I+")
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(.verbose, "B1")

  # function
  rs_results <- rs_sc_get_top_genes_perc(
    f_path_cell = get_rust_count_cell_f_path(object),
    top_n_vals = top_n_vals,
    cell_indices = get_cells_to_keep(object),
    streaming = streaming,
    verbose = .verbose
  )

  names(rs_results) <- sprintf("top_%i_genes_percentage", top_n_vals)

  class(rs_results) <- "sc_proportion_res"
  attr(rs_results, "cell_indices") <- get_cells_to_keep(object)

  duckdb_con <- get_sc_duckdb(object)

  duckdb_con$join_data_obs(get_obs_data(rs_results))

  return(object)
}

### gene set version -----------------------------------------------------------

#' Calculate the proportions of reads for specific gene sets
#'
#' @description
#' This is a helper function that calculates proportions of reads belonging to
#' given gene sets. This can be used for example for the calculation of
#' percentage mitochondrial reads per cell. These will be automatically added
#' to the obs table
#'
#' @param object `SingleCells` class.
#' @param gene_set_list A named list with each element containing the gene
#' identifiers of that set. These should be the same as
#' `get_gene_names(object)`!
#' @param streaming Boolean. Shall the cells be streamed in. Useful for larger
#' data sets where you wish to avoid loading in the whole data. Default to
#' `FALSE`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return It will add the columns based on the names in the `gene_set_list` to
#' the obs table.
#'
#' @export
gene_set_proportions_sc <- S7::new_generic(
  name = "gene_set_proportions_sc",
  dispatch_args = "object",
  fun = function(
    object,
    gene_set_list,
    streaming = FALSE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)


#' @method gene_set_proportions_sc SingleCells
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(gene_set_proportions_sc, SingleCells) <- function(
  object,
  gene_set_list,
  streaming = FALSE,
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::SingleCells")
  checkmate::assertList(gene_set_list, names = "named", types = "character")
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(.verbose, "B1")

  gene_set_list_tidy <- purrr::map(gene_set_list, \(g) {
    get_gene_indices(object, gene_ids = g, rust_index = TRUE)
  })
  names(gene_set_list_tidy) <- names(gene_set_list)

  rs_results <- rs_sc_get_gene_set_perc(
    f_path_cell = get_rust_count_cell_f_path(object),
    cell_indices = get_cells_to_keep(object),
    gene_set_idx = gene_set_list_tidy,
    streaming = streaming,
    verbose = .verbose
  )

  class(rs_results) <- "sc_proportion_res"
  attr(rs_results, "cell_indices") <- get_cells_to_keep(object)

  duckdb_con <- get_sc_duckdb(object)

  duckdb_con$join_data_obs(get_obs_data(rs_results))

  return(object)
}

## hvg -------------------------------------------------------------------------

# generic found in R/base_generics_sc.R

#' @method find_hvg_sc SingleCells
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(find_hvg_sc, SingleCells) <- function(
  object,
  hvg_no = 2000L,
  hvg_params = params_sc_hvg(),
  streaming = FALSE,
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::SingleCells")
  checkmate::qassert(hvg_no, "I1")
  assertScHvg(hvg_params)
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(.verbose, "B1")

  if (length(get_cells_to_keep(object)) == 0) {
    warning(paste(
      "You need to set the cells to keep with set_cells_to_keep().",
      "Returning class as is."
    ))
    return(object)
  }

  res <- with(
    hvg_params,
    rs_sc_hvg(
      f_path_gene = get_rust_count_gene_f_path(object),
      hvg_method = method,
      cell_indices = get_cells_to_keep(object),
      loess_span = loess_span,
      n_bins = num_bin,
      binning = bin_method,
      clip_max = NULL,
      streaming = streaming,
      verbose = .verbose
    )
  )

  object <- set_sc_new_var_cols(object = object, data_list = res)

  hvg <- switch(
    hvg_params$method,
    "vst" = order(res$var_std, decreasing = TRUE)[1:hvg_no],
    "dispersion" = order(res$dispersion, decreasing = TRUE)[1:hvg_no],
    "meanvarbin" = order(res$dispersion_scaled, decreasing = TRUE)[1:hvg_no],
    stop("Unknown HVG method: ", hvg_params$method)
  )

  object <- set_hvg(object, hvg = hvg)

  return(object)
}

## dimension reduction and knn/snn ---------------------------------------------

### pca ------------------------------------------------------------------------

# generic found in R/base_generics_sc.R

#' @method calculate_pca_sc SingleCells
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(calculate_pca_sc, SingleCells) <- function(
  object,
  no_pcs,
  randomised_svd = TRUE,
  sparse_svd = FALSE,
  hvg = NULL,
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::SingleCells")
  checkmate::qassert(no_pcs, "I1")
  checkmate::qassert(randomised_svd, "B1")
  checkmate::qassert(sparse_svd, "B1")
  checkmate::qassert(hvg, c("I+", "0"))
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  if ((length(get_hvg(object)) == 0) && is.null(hvg)) {
    warning(paste(
      "No HVGs identified in the object nor provided.",
      "Please run find_hvg_sc() or provide the indices of the HVG",
      "Returning object as is."
    ))
    return(object)
  }

  selected_hvg <- if (!is.null(hvg)) {
    if (.verbose) {
      message(
        paste(
          "HVGs provided.",
          "Will use these ones and set the internal HVG to the provided genes."
        )
      )
    }
    object <- set_hvg(object, hvg) # this one deals with zero/one indexing internally
    hvg - 1L
  } else {
    get_hvg(object)
  }

  # swap to sparse SVD for large data sets
  n_cells <- length(get_cells_to_keep(object))

  if (n_cells > 500000 & !sparse_svd) {
    message(paste(
      "More than 500,000 cells with sparse SVD = FALSE",
      "Setting sparse SVD to TRUE to avoid high memory pressure."
    ))

    sparse_svd <- TRUE
  }

  # dense path
  if (!sparse_svd) {
    if (.verbose) {
      message(
        sprintf(
          "Using dense SVD solving on scaled data on %i HVG.",
          length(selected_hvg)
        )
      )
    }
    zeallot::`%<-%`(
      c(pca_factors, pca_loadings, singular_values, scaled),
      rs_sc_pca(
        f_path_gene = get_rust_count_gene_f_path(object),
        no_pcs = no_pcs,
        random_svd = randomised_svd,
        cell_indices = get_cells_to_keep(object),
        gene_indices = selected_hvg,
        seed = seed,
        return_scaled = FALSE,
        verbose = .verbose
      )
    )

    object <- set_pca_factors(object, pca_factors)
    object <- set_pca_loadings(object, pca_loadings)
    object <- set_pca_singular_vals(object, singular_values[1:no_pcs])

    return(object)
  } else {
    if (.verbose) {
      message(
        sprintf(
          "Using sparse SVD solving on scaled data on %i HVG.",
          length(selected_hvg)
        )
      )
    }

    zeallot::`%<-%`(
      c(sparse_pca_factors, sparse_pca_loadings, sparse_pca_eigenvals),
      rs_sc_pca_sparse(
        f_path_gene = bixverse:::get_rust_count_gene_f_path(object),
        no_pcs = no_pcs,
        random_svd = randomised_svd,
        cell_indices = get_cells_to_keep(object),
        gene_indices = selected_hvg,
        seed = seed,
        verbose = .verbose
      )
    )

    object <- set_pca_factors(object, sparse_pca_factors)
    object <- set_pca_loadings(object, sparse_pca_loadings)
    object <- set_pca_singular_vals(
      object,
      sparse_pca_eigenvals[1:no_pcs]
    )

    return(object)
  }
}

### neighbours -----------------------------------------------------------------

# generic found in R/base_generics_sc.R
# method shared across SingleCells and MetaCells

#' @export
S7::method(find_neighbours_sc, ScOrMc) <- function(
  object,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  neighbours_params = params_sc_neighbours(),
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) || S7::S7_inherits(object, MetaCells)
  )
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  assertScNeighbours(neighbours_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  if (!embd_to_use %in% get_available_embeddings(object)) {
    warning("The desired embedding was not found. Returning class as is.")
    return(object)
  }

  embd <- get_embedding(x = object, embd_name = embd_to_use)
  if (!is.null(no_embd_to_use)) {
    embd <- embd[, 1:min(no_embd_to_use, ncol(embd))]
  }

  if (.verbose) {
    message(sprintf(
      "Generating kNN data with %s method.",
      neighbours_params$knn_algorithm
    ))
  }
  knn_data <- generate_sc_knn(
    data = embd,
    neighbours_params = neighbours_params,
    seed = seed,
    .verbose = .verbose
  )
  object <- set_knn(object, knn_data)

  if (.verbose) {
    message(sprintf(
      "Generating sNN graph (full: %s).",
      neighbours_params$full_snn
    ))
  }
  snn_graph_rs <- with(
    neighbours_params,
    rs_sc_snn(
      knn_mat = get_knn_mat(knn_data),
      snn_method = snn_similarity,
      pruning = pruning,
      limited_graph = !full_snn,
      verbose = .verbose
    )
  )

  if (.verbose) {
    message("Transforming sNN data to igraph.")
  }
  snn_g <- igraph::make_empty_graph(n = nrow(embd), directed = FALSE)
  snn_g <- igraph::add_edges(
    snn_g,
    snn_graph_rs$edges,
    attr = list(weight = snn_graph_rs$weights)
  )

  set_snn_graph(object, snn_graph = snn_g)
}

### clustering -----------------------------------------------------------------

# generic found in R/base_generics_sc.R
# method shared across SingleCells and MetaCells

S7::method(find_clusters_sc, ScOrMc) <- function(
  object,
  res = 1,
  name = "leiden_clustering"
) {
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) || S7::S7_inherits(object, MetaCells)
  )
  checkmate::qassert(res, "N1")
  checkmate::qassert(name, "S1")

  snn_graph <- get_snn_graph(object)
  if (is.null(snn_graph)) {
    warning(
      "No sNN graph found. Did you run find_neighbours_sc(). Returning class as is."
    )
    return(object)
  }

  leiden_clusters <- igraph::cluster_leiden(
    snn_graph,
    objective_function = "modularity",
    resolution = res
  )

  object[[name]] <- leiden_clusters$membership
  object
}

### knn with distances ---------------------------------------------------------

#' Generate the KNN data with distances
#'
#' @description
#' This function will generate the kNNs based on a given embedding. Available
#' algorithms are:
#' \itemize{
#'   \item `hnsw` - Hierarchical Navigable Small World. A graph-based
#'   approximate nearest neighbour search algorithm; works well on large data
#'   sets. A benign race condition is leveraged during index build, making the
#'   build non-deterministic. Bigger impact on smaller data sets.
#'   \item `ivf` - Inverted file index. Uses first k-means clustering to
#'   identify Voronoi cells and leverages these during querying. Works well
#'   on large data sets with high dimensionality.
#'   \item `nndescent` - Nearest neighbour descent. Similar to `PyNNDescent`,
#'   uses a first index to initialise the graph. Good all-rounder.
#'   \item `annoy` - Approximate nearest neighbours Oh Yeah. Tree-based index,
#'   used across different R single cell packages (Seurat, SCE). This version
#'   is purely memory-based.
#'   \item `exhaustive` - An exhaustive, flat index. On smaller data sets often
#'   faster than the approximate nearest neighbour search algorithms.
#' }
#'
#' @param object `SingleCells` class.
#' @param embd_to_use String. The embedding to use. Whichever you chose, it
#' needs to be part of the object.
#' @param cells_to_use String. Optional cell names to include in the generation
#' of the kNN graph. If `NULL` all (filtered) cells in the object will be used.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param neighbours_params List. Output of [bixverse::params_sc_neighbours()].
#' A list with the following items:
#' \itemize{
#'   \item full_snn - Boolean. Shall the full shared nearest neighbour graph
#'   be generated that generates edges between all cells instead of between
#'   only neighbours. Not used for this function.
#'   \item pruning - Numeric. Weights below this threshold will be set to 0 in
#'   the generation of the sNN graph. Not used for this function.
#'   \item snn_similarity - String. One of `c("rank", "jaccard")`. Defines how
#'   the weight from the SNN graph is calculated. For details, please see
#'   [bixverse::params_sc_neighbours()]. Not used for this function.
#'   \item knn - List of kNN parameters. See [bixverse::params_knn_defaults()]
#'   for available parameters and their defaults.
#' }
#' @param seed Integer. For reproducibility.
#' @param .validate_index Boolean. Shall an exhaustive search against a subset
#' of cells be run to validate the approximate nearest neighbour index.
#' @param .verbose Boolean. Controls verbosity and returns run times.
#'
#' @return Initialised `sc_knn` with the kNN data.
#'
#' @export
generate_knn_sc <- S7::new_generic(
  name = "generate_knn_sc",
  dispatch_args = "object",
  fun = function(
    object,
    embd_to_use = "pca",
    cells_to_use = NULL,
    no_embd_to_use = NULL,
    neighbours_params = params_sc_neighbours(),
    seed = 42L,
    .validate_index = TRUE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method generate_knn_sc SingleCells
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(generate_knn_sc, SingleCells) <- function(
  object,
  embd_to_use = "pca",
  cells_to_use = NULL,
  no_embd_to_use = NULL,
  neighbours_params = params_sc_neighbours(),
  seed = 42L,
  .validate_index = TRUE,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  checkmate::qassert(cells_to_use, c("S+", "0"))
  assertScNeighbours(neighbours_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.validate_index, "B1")
  checkmate::qassert(.verbose, "B1")

  # function body
  embd <- get_embedding(x = object, embd_name = embd_to_use)

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

  if (!is.null(cells_to_use)) {
    embd <- embd[which(rownames(embd) %in% cells_to_use), ]
  }

  knn_data <- rs_sc_knn_w_dist(
    embd = embd,
    knn_params = neighbours_params,
    verbose = .verbose,
    validate_index = .validate_index,
    seed = seed
  )

  knn_obj <- new_sc_knn(knn_data = knn_data, used_cells = row.names(embd))

  knn_obj
}
