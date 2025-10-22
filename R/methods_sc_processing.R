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
#' @param object `single_cell_exp` class.
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

#' @method scrublet_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(scrublet_sc, single_cell_exp) <- function(
  object,
  scrublet_params = params_scrublet(),
  seed = 42L,
  streaming = FALSE,
  return_combined_pca = FALSE,
  return_pairs = FALSE,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))
  assertScScrublet(scrublet_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(return_combined_pca, "B1")
  checkmate::qassert(return_pairs, "B1")
  checkmate::qassert(.verbose, "B1")

  # function body
  cells_to_keep <- get_cells_to_keep(object)

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
  class(scrublet_res) <- "scrublet_res"

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
#' @param object `single_cell_exp` class.
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

#' @method doublet_detection_boost_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(doublet_detection_boost_sc, single_cell_exp) <- function(
  object,
  boost_params = params_boost(),
  seed = 42L,
  streaming = FALSE,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))
  assertScBoost(boost_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(.verbose, "B1")

  # function body
  cells_to_keep <- get_cells_to_keep(object)

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
  class(boost_res) <- "boost_res"

  return(boost_res)
}

## gene proportions ------------------------------------------------------------

#' Calculate the proportions of reads for specific gene sets
#'
#' @description
#' This is a helper function that calculates proportions of reads belonging to
#' given gene sets. This can be used for example for the calculation of
#' percentage mitochondrial reads per cell. These will be automatically added
#' to the obs table
#'
#' @param object `single_cell_exp` class.
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

#' @method gene_set_proportions_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(gene_set_proportions_sc, single_cell_exp) <- function(
  object,
  gene_set_list,
  streaming = FALSE,
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::single_cell_exp")
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

#' Identify HVGs
#'
#' @description
#' This is a helper function to identify highly variable genes. At the moment
#' the implementation has only the VST-based version (known as Seurat v3). The
#' other methods will be implemented in the future.
#'
#' @param object `single_cell_exp` class.
#' @param hvg_no Integer. Number of highly variable genes to include. Defaults
#' to `2000L`.
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
#' @return It will add the mean, var, var_exp, var_std of each gene to the
#' the var table.
#'
#' @export
find_hvg_sc <- S7::new_generic(
  name = "find_hvg_sc",
  dispatch_args = "object",
  fun = function(
    object,
    hvg_no = 2000L,
    hvg_params = params_sc_hvg(),
    streaming = FALSE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method find_hvg_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(find_hvg_sc, single_cell_exp) <- function(
  object,
  hvg_no = 2000L,
  hvg_params = params_sc_hvg(),
  streaming = FALSE,
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::single_cell_exp")
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
      clip_max = NULL,
      streaming = streaming,
      verbose = .verbose
    )
  )

  object <- set_sc_new_var_cols(object = object, data_list = res)

  hvg <- order(res$var_std, decreasing = TRUE)[1:hvg_no]

  object <- set_hvg(object, hvg = hvg)

  return(object)
}

## dimension reduction and knn/snn ---------------------------------------------

### pca ------------------------------------------------------------------------

#' Run PCA for single cell
#'
#' @description
#' This function will run PCA (option of full SVD and randomised SVD for now)
#' on the detected highly variable genes.
#'
#' @param object `single_cell_exp` class.
#' @param no_pcs Integer. Number of PCs to calculate.
#' @param randomised_svd Boolean. Shall randomised SVD be used. Faster, but
#' less precise.
#' @param seed Integer. Controls reproducibility. Only relevant if
#' `randomised_svd = TRUE`.
#' @param .verbose Boolean. Controls verbosity and returns run times.
#'
#' @return The function will add the PCA factors, loadings and singular values
#' to the object cache in memory.
#'
#' @export
calculate_pca_sc <- S7::new_generic(
  name = "calculate_pca_sc",
  dispatch_args = "object",
  fun = function(
    object,
    no_pcs,
    randomised_svd = TRUE,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calculate_pca_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(calculate_pca_sc, single_cell_exp) <- function(
  object,
  no_pcs,
  randomised_svd = TRUE,
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::qassert(no_pcs, "I1")
  checkmate::qassert(randomised_svd, "B1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  if (length(get_hvg(object)) == 0) {
    warning(paste(
      "No HVGs identified in this object. Did you run find_hvg_sc()?",
      "Returning object as is."
    ))
    return(object)
  }

  zeallot::`%<-%`(
    c(pca_factors, pca_loadings, singular_values, scaled),
    rs_sc_pca(
      f_path_gene = get_rust_count_gene_f_path(object),
      no_pcs = no_pcs,
      random_svd = randomised_svd,
      cell_indices = get_cells_to_keep(object),
      gene_indices = get_hvg(object),
      seed = seed,
      return_scaled = FALSE,
      verbose = .verbose
    )
  )

  object <- set_pca_factors(object, pca_factors)
  object <- set_pca_loadings(object, pca_loadings)
  object <- set_pca_singular_vals(object, singular_values[1:no_pcs])

  return(object)
}

### neighbours -----------------------------------------------------------------

#' Find the neighbours for single cell.
#'
#' @description
#' This function will generate the kNNs based on a given embedding (atm,
#' only option is PCA). Two different algorithms are implemented with different
#' speed and accuracy to approximate the nearest neighbours. `"annoy"` is more
#' rapid and based on the `Approximate Nearest Neigbours Oh Yeah` algorithm,
#' whereas `"hnsw"` implements a `Hierarchical Navigatable Small Worlds` vector
#' search that is slower, but more precise. Subsequently, the kNN data will
#' be used to generate an sNN igraph for clustering methods.
#'
#' @param object `single_cell_exp` class.
#' @param embd_to_use String. The embedding to use. Whichever you chose, it
#' needs to be part of the object.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param neighbours_params List. Output of [bixverse::params_sc_neighbours()].
#' A list with the following items:
#' \itemize{
#'   \item k - Integer. Number of neighbours to identify.
#'   \item n_trees -  Integer. Number of trees to use for the `annoy` algorithm.
#'   The higher, the longer the algorithm takes, but the more precise the
#'   approximated nearest neighbours.
#'   \item search_budget - Integer. Search budget per tree for the `annoy`
#'   algorithm. The higher, the longer the algorithm takes, but the more precise
#'   the approximated nearest neighbours.
#'   \item knn_algorithm - String. One of `c("annoy", "hnsw")`. `"hnsw"` takes
#'   longer, is more precise and more memory friendly. `"annoy"` is faster, less
#'   precise and will take more memory.
#'   \item ann_dist - String. One of `c("cosine", "euclidean")`.
#'   \item full_snn - Boolean. Shall the sNN graph be generated across all
#'   cells (standard in the `bluster` package.) Defaults to `FALSE`.
#'   \item pruning - Value below which the weight in the sNN graph is set to 0.
#'   \item snn_similarity - String. One of `c("rank", "jaccard")`. Defines how
#'   the weight form the SNN graph is calculated. For details, please see
#'   [bixverse::params_sc_neighbours()].
#' }
#' @param seed Integer. For reproducibility.
#' @param .verbose Boolean. Controls verbosity and returns run times.
#'
#' @return The object with added KNN matrix.
#'
#' @export
find_neighbours_sc <- S7::new_generic(
  name = "find_neighbours_sc",
  dispatch_args = "object",
  fun = function(
    object,
    embd_to_use = "pca",
    no_embd_to_use = NULL,
    neighbours_params = params_sc_neighbours(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method find_neighbours_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(find_neighbours_sc, single_cell_exp) <- function(
  object,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  neighbours_params = params_sc_neighbours(),
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  assertScNeighbours(neighbours_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  # get the embedding
  checkmate::assertTRUE(embd_to_use %in% get_available_embeddings(object))
  embd <- get_embedding(x = object, embd_name = embd_to_use)

  if (!is.null(no_embd_to_use)) {
    to_take <- min(c(no_embd_to_use, ncol(embd)))
    embd <- embd[, 1:to_take]
  }

  if (.verbose) {
    message(sprintf(
      "Generating kNN data with %s method.",
      neighbours_params$knn_algorithm
    ))
  }

  knn_data <- with(
    neighbours_params,
    rs_sc_knn(
      embd = embd,
      no_neighbours = k,
      seed = seed,
      n_trees = n_trees,
      search_budget = search_budget,
      verbose = .verbose,
      algorithm_type = knn_algorithm,
      ann_dist = ann_dist
    )
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
      knn_data,
      snn_method = snn_similarity,
      pruning = pruning,
      limited_graph = !full_snn,
      verbose = .verbose
    )
  )

  if (.verbose) {
    message("Transforming sNN data to igraph.")
  }

  snn_g <- igraph::make_graph(snn_graph_rs$edges + 1, directed = FALSE)
  igraph::E(snn_g)$weight <- snn_graph_rs$weights

  object <- set_snn_graph(object, snn_graph = snn_g)

  return(object)
}

### clustering -----------------------------------------------------------------

#' Graph-based clustering of cells on the sNN graph
#'
#' @description
#' This function will apply Leiden clustering on the sNN graph with the
#' given resolution and add a column to the obs table.
#'
#' @param object `single_cell_exp` class.
#' @param res Numeric. The resolution parameter for [igraph::cluster_leiden()].
#' @param name String. The name to add to the obs table in the DuckDB.
#'
#' @return The object with added clustering in the obs table.
#'
#' @export
find_clusters_sc <- S7::new_generic(
  name = "find_clusters_sc",
  dispatch_args = "object",
  fun = function(
    object,
    res = 1,
    name = "leiden_clustering"
  ) {
    S7::S7_dispatch()
  }
)

#' @method find_clusters_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(find_clusters_sc, single_cell_exp) <- function(
  object,
  res = 1,
  name = "leiden_clustering"
) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::qassert(res, "N1")
  checkmate::qassert(name, "S1")

  snn_graph <- get_snn_graph(object)

  if (is.null(snn_graph)) {
    warning(paste(
      "No sNN graph found. Did you run find_neighbours_sc()",
      "Returning class as is."
    ))
    return(object)
  }

  leiden_clusters <- igraph::cluster_leiden(
    snn_graph,
    objective_function = "modularity",
    resolution = res
  )

  duckdb_con <- get_sc_duckdb(object)

  new_data <- data.table::data.table(
    cell_idx = get_cells_to_keep(object) + 1, # needs to be 1-indexed
    new_data = leiden_clusters$membership
  )
  data.table::setnames(new_data, "new_data", name)

  duckdb_con$join_data_obs(new_data = new_data)

  return(object)
}
