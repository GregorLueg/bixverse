# generics shared between single and meta cells --------------------------------

# contains generics that are shared between single cells, meta cells (and in
# the future also spatial transcriptomics).

## single/meta cells -----------------------------------------------------------

### obs table ------------------------------------------------------------------

#' Getter the obs table
#'
#' @param object `SingleCells`, `MetaCells` class.
#' @param indices Optional integer vector. The integer positions of the cells
#' to return.
#' @param cols Optional string vector. The columns from the obs table to return.
#' @param filtered Boolean. Whether to return all cells or filtered to `to_keep`
#' cells. Not relevant for `MetaCells`.
#'
#' @return The obs table
#'
#' @export
get_sc_obs <- S7::new_generic(
  name = "get_sc_obs",
  dispatch_args = "object",
  fun = function(
    object,
    indices = NULL,
    cols = NULL,
    filtered = FALSE
  ) {
    S7::S7_dispatch()
  }
)

### var table ------------------------------------------------------------------

#' Getter the var table
#'
#' @param object `SingleCells`, `MetaCells` class.
#' @param indices Optional integer vector. The integer positions of the genes
#' to return.
#' @param cols Optional string vector. The columns from the var table to return.
#'
#' @return The vars table
#'
#' @export
get_sc_var <- S7::new_generic(
  name = "get_sc_var",
  dispatch_args = "object",
  fun = function(
    object,
    indices = NULL,
    cols = NULL
  ) {
    S7::S7_dispatch()
  }
)

### counts ---------------------------------------------------------------------

#' Getter the counts
#'
#' @param object `SingleCells`, `MetaCells` class.
#' @param assay String. Which slot to return. One of `c("raw", "norm")`.
#' Defaults to `"raw"`.
#' @param return_format String. One of `c("cell", "gene")`. Return data in
#' cell-centric compressed format (CSR) or gene-centric compressed format (CSC).
#' Defaults to `"cell"`. Not relevant for `MetaCells`.
#' @param cell_indices Optional cell indices.
#' @param gene_indices Optional gene indices.
#' @param use_cells_to_keep Boolean. Shall cells to keep be found in the class,
#' shall the counts be reduced to these. Not relevant for `MetaCells`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return The counts table
#'
#' @export
get_sc_counts <- S7::new_generic(
  name = "get_sc_counts",
  dispatch_args = "object",
  fun = function(
    object,
    assay = c("raw", "norm"),
    return_format = c("cell", "gene"),
    cell_indices = NULL,
    gene_indices = NULL,
    use_cells_to_keep = TRUE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

### available features ---------------------------------------------------------

#' Returns the available features for single cell applications
#'
#' @description
#' Returns a data.table with available features in the obs table and in the
#' count matrices.
#'
#' @param object `SingleCells`, `MetaCells` class.
#'
#' @return A data.table with available features.
#'
#' @export
get_sc_available_features <- S7::new_generic(
  name = "get_sc_available_features",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)

### rename columns -------------------------------------------------------------

#' Rename columns in the obs or var table
#'
#' @description
#' Renames the columns in the obs or var table of single cell-related classes.
#'
#' @param object `SingleCells`, `MetaCells` class.
#' @param table String. One of `c("obs", "var")`. In which of the tables to
#' rename the columns.
#' @param old Character vector. The old column names.
#' @param new Character vector. The new column names.
#'
#' @return Invisible self
#'
#' @export
setnames_sc <- S7::new_generic(
  name = "setnames_sc",
  dispatch_args = "object",
  fun = function(
    object,
    table = c("obs", "var"),
    old,
    new
  ) {
    S7::S7_dispatch()
  }
)

### ScMap ----------------------------------------------------------------------

#### setters -------------------------------------------------------------------

#' Set gene mapping for ScMap object
#'
#' @param x An object to set gene mapping for
#' @param gene_map Named integer indicating indices and names of the genes
#'
#' @export
#'
#' @keywords internal
set_gene_mapping <- function(x, gene_map) {
  UseMethod("set_gene_mapping")
}

#' Set cell mapping for ScMap object
#'
#' @param x An object to set cell mapping for
#' @param cell_map Named integer indicating indices and names of the cells
#'
#' @export
#'
#' @keywords internal
set_cell_mapping <- function(x, cell_map) {
  UseMethod("set_cell_mapping")
}

#' Set cells to keep for ScMap object
#'
#' @param x An object to set cells to keep for
#' @param cells_to_keep String or integer. The names or indices of the cells
#' to keep in downstream analysis.
#'
#' @export
#'
#' @keywords internal
set_cells_to_keep <- function(x, cells_to_keep) {
  UseMethod("set_cells_to_keep")
}

#' Set the HVG genes
#'
#' @param x An object to set the HVGs for
#' @param hvg String or integer. The names or indices of the highly variable
#' genes.
#'
#' @export
#'
#' @keywords internal
set_hvg <- function(x, hvg) {
  UseMethod("set_hvg")
}

#### getters -------------------------------------------------------------------

#' Get the HVG
#'
#' @param x An object to get HVG from.
#'
#' @export
get_hvg <- function(x) {
  UseMethod("get_hvg")
}

#' Get the gene names
#'
#' @param x An object to get the gene names from.
#'
#' @export
get_gene_names <- function(x) {
  UseMethod("get_gene_names")
}

#' Get the cell names
#'
#' @param x An object to get the cell names from.
#' @param filtered Boolean. Shall, if found only the cell names of the
#' `cells_to_keep` be returned (see [bixverse::set_cells_to_keep()]. Defaults
#' to `FALSE`
#'
#' @export
get_cell_names <- function(x, filtered = FALSE) {
  UseMethod("get_cell_names")
}

#' Get the index position for a gene
#'
#' @param x An object to get the gene index from.
#' @param gene_ids String vector. The gene ids to search for.
#' @param rust_index Bool. Shall rust-based indexing be returned.
#'
#' @export
get_gene_indices <- function(x, gene_ids, rust_index) {
  UseMethod("get_gene_indices")
}

#' Get the index position for a gene
#'
#' @param x An object to get the gene index from.
#' @param cell_ids String vector. The cell ids to search for.
#' @param rust_index Bool. Shall rust-based indexing be returned.
#'
#' @export
get_cell_indices <- function(x, cell_ids, rust_index) {
  UseMethod("get_cell_indices")
}

#' Get the cells to keep
#'
#' @param x An object to get the gene index from.
#'
#' @export
get_cells_to_keep <- function(x) {
  UseMethod("get_cells_to_keep")
}

#' Get the gene names based on the gene idx
#'
#' @param x An object to get the gene names from.
#' @param gene_idx Integer. The original gene indices for which to return
#' the gene names.
#' @param rust_based Boolean. Is it Rust-based, i.e., 0-index or R-based, i.e.,
#' 1-indexed.
#'
#' @export
get_gene_names_from_idx <- function(x, gene_idx, rust_based = TRUE) {
  UseMethod("get_gene_names_from_idx")
}

### ScCache --------------------------------------------------------------------

#### setters -------------------------------------------------------------------

#' Set/add PCA factors
#'
#' @param x An object to add the PCA factors for.
#' @param pca_factor Numerical matrix. The matrix with the PCA factors.
#'
#' @export
#'
#' @keywords internal
set_pca_factors <- function(x, pca_factor) {
  UseMethod("set_pca_factors")
}

#' Set/add PCA loadings
#'
#' @param x An object to add the PCA loadings for.
#' @param pca_loading Numerical matrix. The Matrix with the PCA loadings.
#'
#' @export
#'
#' @keywords internal
set_pca_loadings <- function(x, pca_loading) {
  UseMethod("set_pca_loadings")
}

#' Set/add PCA singular values
#'
#' @param x An object to add the singular values for.
#' @param singular_vals Numerical vector. The singular values.
#'
#' @export
#'
#' @keywords internal
set_pca_singular_vals <- function(x, singular_vals) {
  UseMethod("set_pca_singular_vals")
}

#' Add additional embeddings to the class
#'
#' @param x An object to add the singular values for.
#' @param embd Numerical matrix representing the additional embedding.
#' @param name String. Name of the embedding.
#'
#' @export
#'
#' @keywords internal
set_embedding <- function(x, embd, name) {
  UseMethod("set_embedding")
}

#' Set/add KNN
#'
#' @param x An object to add the KNN data to
#' @param knn `SingleCellNearestNeighbour` class to add to the classes.
#'
#' @export
#'
#' @keywords internal
set_knn <- function(x, knn) {
  UseMethod("set_knn")
}

#' Set/add KNN
#'
#' @param x An object to add the KNN data to
#' @param snn_graph Igraph. The sNN graph for subsequent clustering.
#'
#' @export
#'
#' @keywords internal
set_snn_graph <- function(x, snn_graph) {
  UseMethod("set_snn_graph")
}

#' Remove the KNN data
#'
#' @param x An object from which to remove the kNN data
#'
#' @export
#'
#' @keywords internal
remove_knn <- function(x) {
  UseMethod("remove_knn")
}

#' Remove the sNN graph
#'
#' @param x An object from which to remove the sNN graph
#'
#' @export
#'
#' @keywords internal
remove_snn_graph <- function(x) {
  UseMethod("remove_snn_graph")
}

#### getters -------------------------------------------------------------------

#' Get the PCA factors
#'
#' @param x An object to get PCA factors from.
#'
#' @returns The PCA factors from the object (if found).
#'
#' @export
get_pca_factors <- function(x) {
  UseMethod("get_pca_factors")
}

#' Get the PCA loadings
#'
#' @param x An object to get PCA loadings from.
#'
#' @returns The PCA feature loadings from the object (if found).
#'
#' @export
get_pca_loadings <- function(x) {
  UseMethod("get_pca_loadings")
}

#' Get the PCA singular values
#'
#' @param x An object to get PCA singular values from.
#'
#' @returns The PCA singular values from the object (if found).
#'
#' @export
get_pca_singular_val <- function(x) {
  UseMethod("get_pca_singular_val")
}

#' Get the embedding from the cache
#'
#' @description
#' General wrapper function that can be used to pull out any embedding stored
#' in the `ScCache`.
#'
#' @param x An object to get embedding from
#' @param embd_name String. The name of the embedding to return. The function
#' will throw an error if the embedding does not exist.
#'
#' @returns Get the specified embeddings from the object (if found).
#'
#' @export
get_embedding <- function(x, embd_name) {
  UseMethod("get_embedding")
}

#' Get the available embeddings from the cache
#'
#' @description
#' Returns the available embedding names from the cache.
#'
#' @param x An object to get embedding from
#'
#' @return Get the names of the available embeddings.
#'
#' @export
get_available_embeddings <- function(x) {
  UseMethod("get_available_embeddings")
}

#' Get the sNN graph
#'
#' @param x An object to get the sNN graph from.
#'
#' @returns The igraph that has the shared nearest neighbours.
#'
#' @export
get_snn_graph <- function(x) {
  UseMethod("get_snn_graph")
}

#' Get the KNN object
#'
#' @param x An object to get the KNN class from.
#'
#' @returns The `SingleCellNearestNeighbour` object.
#'
#' @export
get_knn_obj <- function(x) {
  UseMethod("get_knn_obj")
}

### others ---------------------------------------------------------------------

#### knn getter methods --------------------------------------------------------

#' Get the KNN matrix
#'
#' @param x An object to get the kNN matrix from.
#'
#' @export
get_knn_mat <- function(x) {
  UseMethod("get_knn_mat")
}

#' Get the KNN distance measures
#'
#' @param x An object to get the kNN distances from.
#'
#' @export
get_knn_dist <- function(x) {
  UseMethod("get_knn_dist")
}

### methods --------------------------------------------------------------------

#### hvg -----------------------------------------------------------------------

#' Identify HVGs
#'
#' @description
#' This is a helper function to identify highly variable genes for `SingleCells`
#' (using the Rust-based streaming of data) or `MetaCells`.
#'
#' @param object `SingleCells`, `MetaCells` class.
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
#' `FALSE`. Not used for `MetaCells`.
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

#### pca -----------------------------------------------------------------------

#' Run PCA for single cell
#'
#' @description
#' This function will run PCA on the detected highly variable genes. You can
#' use randomised SVD for speed and there is an option for sparse SVD for very
#' large data sets to avoid memory pressure.
#'
#' @param object `SingleCells`, `MetaCells` class.
#' @param no_pcs Integer. Number of PCs to calculate.
#' @param randomised_svd Boolean. Shall randomised SVD be used. Faster, but
#' less precise.
#' @param sparse_svd Boolean. Shall sparse solvers be used that do not do
#' scaling. If set to yes, in the case of `random_svd = FALSE`, Lanczos
#' iterations are used to solve the sparse SVD. With `random_svd = TRUE`, the
#' sparse initial matrix is multiplied with the random matrix, yielding a
#' much smaller dense matrix that does not increase the memory pressure
#' massively. Not used for `MetaCells`.
#' @param hvg Optional integer. If you want to provide your own HVG genes.
#' Otherwise, the function will default to what is found in
#' [bixverse::get_hvg()]. Please provide 1-indexed genes here! If you provide
#' these, the internal HVG will be overwritten.
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
    sparse_svd = FALSE,
    hvg = NULL,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#### nearest neighbours --------------------------------------------------------

#' Find the neighbours for single cell.
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
#' Subsequently, the kNN graph will be additionally transformed into a shared
#' nearest neighbour graph for clustering methods.
#'
#' @param object `SingleCells`, `MetaCells` class.
#' @param embd_to_use String. The embedding to use. Whichever you chose, it
#' needs to be part of the object.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param neighbours_params List. Output of [bixverse::params_sc_neighbours()].
#' A list with the following items:
#' \itemize{
#'   \item full_snn - Boolean. Shall the full shared nearest neighbour graph
#'   be generated that generates edges between all cells instead of between
#'   only neighbours.
#'   \item pruning - Numeric. Weights below this threshold will be set to 0 in
#'   the generation of the sNN graph.
#'   \item snn_similarity - String. One of `c("rank", "jaccard")`. Defines how
#'   the weight from the SNN graph is calculated. For details, please see
#'   [bixverse::params_sc_neighbours()].
#'   \item knn - List of kNN parameters. See [bixverse::params_knn_defaults()]
#'   for available parameters and their defaults.
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

#### clustering ----------------------------------------------------------------

#' Graph-based clustering of cells on the sNN graph
#'
#' @description
#' This function will apply Leiden clustering on the sNN graph with the
#' given resolution and add a column to the obs table.
#'
#' @param object `SingleCells`, `MetaCells` class.
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
