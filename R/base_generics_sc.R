# generics shared between single and meta cells --------------------------------

# contains generics that are shared between single cells, meta cells (and in
# the future also spatial transcriptomics).

## single/meta cells -----------------------------------------------------------

### obs table ------------------------------------------------------------------

#' Getter the obs table
#'
#' @param object `SingleCells`, `MetaCells`, `SingleCellsMultiModal` class.
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
#' @param object `SingleCells`, `MetaCells`, `SingleCellsMultiModal` class.
#' @param indices Optional integer vector. The integer positions of the genes
#' to return.
#' @param cols Optional string vector. The columns from the var table to return.
#' @param modality String. The modality to return. One of `c("rna", "adt")`.
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
    cols = NULL,
    modality = c("rna", "adt")
  ) {
    S7::S7_dispatch()
  }
)

### counts ---------------------------------------------------------------------

#' Getter the counts
#'
#' @param object `SingleCells`, `MetaCells`, `SingleCellsMultiModal` class.
#' @param assay String. Which slot to return. One of `c("raw", "norm")`.
#' Defaults to `"raw"`.
#' @param return_format String. One of `c("cell", "gene")`. Return data in
#' cell-centric compressed format (CSR) or gene-centric compressed format (CSC).
#' Defaults to `"cell"`. Not relevant for `MetaCells`.
#' @param cell_indices Optional cell indices.
#' @param gene_indices Optional gene indices.
#' @param modality String. The modality to return. One of `c("rna", "adt")`.
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
    modality = c("rna", "adt"),
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
#' @param object `SingleCells`, `MetaCells` (or potentially other) class.
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
#' @param object `SingleCells`, `MetaCells` (or potentially other) class.
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

### drop columns ---------------------------------------------------------------

#' Drop columns from the obs or var table
#'
#' @description
#' Drops the named columns from the obs or var table of single cell-related
#' classes. Protected identifier and bookkeeping columns (`cell_idx`,
#' `cell_id`, `to_keep` for obs; `gene_idx`, `gene_id` for var) are refused
#' with a warning. Columns that do not exist trigger a warning and are
#' skipped.
#'
#' @param object `SingleCells` (or other compatible) class.
#' @param table String. One of `c("obs", "var")`.
#' @param cols Character vector. Column names to drop.
#'
#' @return Invisible self.
#'
#' @export
drop_cols_sc <- S7::new_generic(
  name = "drop_cols_sc",
  dispatch_args = "object",
  fun = function(
    object,
    table = c("obs", "var"),
    cols
  ) {
    S7::S7_dispatch()
  }
)

### ScMap ----------------------------------------------------------------------

#### setters -------------------------------------------------------------------

#' @title Set gene mapping
#'
#' @description Set a gene mapping for a given object. This is used for the
#' single cell-related classes with streaming from disk.
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

#' @title Set cell mapping
#'
#' @description Set a cell mapping for a given object. This is used for the
#' single cell-related classes with streaming from disk.
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

#' @title Set cells to keep
#'
#' @description Set the cells to keep. This is used for the single cell-related
#' classes with streaming from disk and tells subsequent (Rust) methods which
#' cells to include.
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

#' @title Reset the cells to keep
#'
#' @description
#' Restores every cell found in the binary count file and wipes the cache,
#' taking the object back to a pristine state. Filtering only ever flips a
#' `to_keep` flag in the DuckDB, so nothing was deleted and nothing is lost by
#' resetting.
#'
#' Wiping the cache is not optional: a PCA or a kNN computed on a filtered
#' subset does not describe the full cell set, and keeping it would recreate
#' exactly the mismatch this guards against. With `force = FALSE` you are asked
#' to confirm before that happens.
#'
#' @param object `SingleCells` or `SingleCellsMultiModal` class.
#' @param force Boolean. Skip the confirmation prompt. Defaults to `FALSE`, in
#' which case an interactive session asks before wiping the cache and a
#' non-interactive one errors, because there is no one there to ask.
#'
#' @returns The object with every cell restored and an empty cache. Unchanged
#' if the confirmation was declined.
#'
#' @export
reset_cells_to_keep <- S7::new_generic(
  name = "reset_cells_to_keep",
  dispatch_args = "object",
  fun = function(
    object,
    force = FALSE
  ) {
    S7::S7_dispatch()
  }
)

#' @title Set the HVG genes
#'
#' @description
#' Stores within the class the index positions of the HVG. This is used for
#' the single cell-related classes and methods.
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

#' @title Get the HVG
#'
#' @description
#' Returns the HVG indices. Pending class type this are 1-based (for R) or
#' 0-based for Rust.
#'
#' @param x An object to get HVG from.
#'
#' @returns Indices of the stored HVG genes.
#'
#' @export
get_hvg <- function(x) {
  UseMethod("get_hvg")
}

#' @title Get the gene names
#'
#' @description
#' Get the main gene names (for example symbols or Ensembl identifiers).
#'
#' @param x An object to get the gene names from.
#'
#' @return The primary gene identifiers stored in the class.
#'
#' @export
get_gene_names <- function(x) {
  UseMethod("get_gene_names")
}

#' @title Get the cell names
#'
#' @description
#' Returns the cell names (usually barcodes).
#'
#' @param x An object to get the cell names from.
#' @param filtered Boolean. Shall, if found only the cell names of the
#' `cells_to_keep` be returned (see [bixverse::set_cells_to_keep()]. Defaults
#' to `FALSE`
#'
#' @return The cell names (barcodes)
#'
#' @export
get_cell_names <- function(x, filtered = FALSE) {
  UseMethod("get_cell_names")
}

#' @title Get the index position for a gene
#'
#' @description
#' Returns the index for a given gene based on the internal gene mapping. This
#' is used for the single cell-related classes and methods.
#'
#' @param x An object to get the gene index from.
#' @param gene_ids String vector. The gene ids to search for.
#' @param rust_index Bool. Shall Rust-based indexing be returned.
#'
#' @returns The indices of the genes
#'
#' @export
get_gene_indices <- function(x, gene_ids, rust_index) {
  UseMethod("get_gene_indices")
}

#' @title Get the index position for a gene
#'
#' @description
#' Returns the index for a given gene based on the internal gene mapping. This
#' is used for the single cell-related classes and methods.
#'
#' @param x An object to get the gene index from.
#' @param cell_ids String vector. The cell ids to search for.
#' @param rust_index Bool. Shall rust-based indexing be returned.
#'
#' @returns The indices of the cells
#'
#' @export
get_cell_indices <- function(x, cell_ids, rust_index) {
  UseMethod("get_cell_indices")
}

#' @title Get the cell idx (R-based) and cell names
#'
#' @description
#' Returns the cell indices (R-based) and the cell names (usually barcodes)
#' from the object for further downstream usage.
#'
#' @param x An object to get the cell info from
#' @param filtered Boolean. If `TRUE`, only the cells to keep will be returned.
#'
#' @returns A named vector with elements -> cell_idx, names -> cell_names.
#'
#' @export
get_cell_info <- function(x, filtered = TRUE) {
  UseMethod("get_cell_info")
}

#' @title Get the cells to keep
#'
#' @param x An object from which to get the cells to keep from. These are
#' 0-indexed.
#'
#' @returns Integer vector with 0-indices of the cells to keep.
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
#' @param ... Other parameters.
#'
#' @export
#'
#' @keywords internal
set_pca_factors <- function(x, pca_factor, ...) {
  UseMethod("set_pca_factors")
}

#' Set/add PCA loadings
#'
#' @param x An object to add the PCA loadings for.
#' @param pca_loading Numerical matrix. The Matrix with the PCA loadings.
#' @param ... Other parameters.
#'
#' @export
#'
#' @keywords internal
set_pca_loadings <- function(x, pca_loading, ...) {
  UseMethod("set_pca_loadings")
}

#' Set/add PCA singular values
#'
#' @param x An object to add the singular values for.
#' @param singular_vals Numerical vector. The singular values.
#' @param ... Other parameters.
#'
#' @export
#'
#' @keywords internal
set_pca_singular_vals <- function(x, singular_vals, ...) {
  UseMethod("set_pca_singular_vals")
}

#' Add additional embeddings to the class
#'
#' @param x An object to add the singular values for.
#' @param embd Numerical matrix representing the additional embedding.
#' @param name String. Name of the embedding.
#' @param ... Other parameters.
#'
#' @export
#'
#' @keywords internal
set_embedding <- function(x, embd, name, ...) {
  UseMethod("set_embedding")
}

#' Set/add KNN
#'
#' @param x An object to add the KNN data to
#' @param knn `SingleCellNearestNeighbour` class to add to the classes.
#' @param ... Other parameters.
#'
#' @export
#'
#' @keywords internal
set_knn <- function(x, knn, ...) {
  UseMethod("set_knn")
}

#' Set/add KNN
#'
#' @param x An object to add the KNN data to.
#' @param snn_graph Igraph. The sNN graph for subsequent clustering.
#' @param ... Other parameters.
#'
#' @export
#'
#' @keywords internal
set_snn_graph <- function(x, snn_graph, ...) {
  UseMethod("set_snn_graph")
}

#' Remove the KNN data
#'
#' @param x An object from which to remove the kNN data.
#' @param ... Other parameters.
#'
#' @export
#'
#' @keywords internal
remove_knn <- function(x, ...) {
  UseMethod("remove_knn")
}

#' Remove the sNN graph
#'
#' @param x An object from which to remove the sNN graph
#' @param ... Other parameters.
#'
#' @export
#'
#' @keywords internal
remove_snn_graph <- function(x, ...) {
  UseMethod("remove_snn_graph")
}

#' Set/add the MAGIC imputed layer
#'
#' @param x An object to add the imputed layer to.
#' @param magic `ScMagic` class with the imputed counts.
#' @param ... Other parameters.
#'
#' @export
#'
#' @keywords internal
set_magic <- function(x, magic, ...) {
  UseMethod("set_magic")
}

#' Remove the MAGIC imputed layer
#'
#' @param x An object from which to remove the imputed layer.
#' @param ... Other parameters.
#'
#' @export
#'
#' @keywords internal
remove_magic <- function(x, ...) {
  UseMethod("remove_magic")
}

#### getters -------------------------------------------------------------------

#' @title Get the PCA factors
#'
#' @description
#' Returns the PCA factors (sample-based scores). This function is used for the
#' single cell-related classes and methods.
#'
#' @param x An object to get PCA factors from.
#' @param ... Other parameters.
#'
#' @returns The PCA factors from the object (if found).
#'
#' @export
get_pca_factors <- function(x, ...) {
  UseMethod("get_pca_factors")
}

#' @title Get the PCA loadings
#'
#' @description
#' Returns the PCA loadings (feature-based scores). This function is used for
#' the single cell-related classes and methods.
#'
#' @param x An object to get PCA loadings from.
#' @param ... Other parameters.
#'
#' @returns The PCA feature loadings from the object (if found).
#'
#' @export
get_pca_loadings <- function(x, ...) {
  UseMethod("get_pca_loadings")
}

#' @title Get the PCA singular values
#'
#' @description
#' Returns the PCA singular values (can be useful to assess cumulative variance
#' explained). This function is used for the single cell-related classes and
#' methods.
#'
#' @param x An object to get PCA singular values from.
#' @param ... Other parameters.
#'
#' @returns The PCA singular values from the object (if found).
#'
#' @export
get_pca_singular_val <- function(x, ...) {
  UseMethod("get_pca_singular_val")
}

#' @title Get the embedding
#'
#' @description
#' General wrapper function that can be used to pull out any embedding stored
#' in the class. This function is used for the single cell-related classes and
#' methods.
#'
#' @param x An object to get embedding from
#' @param embd_name String. The name of the embedding to return. The function
#' will throw an error if the embedding does not exist.
#' @param ... Other parameters.
#'
#' @returns Get the specified embeddings from the object (if found).
#'
#' @export
get_embedding <- function(x, embd_name, ...) {
  UseMethod("get_embedding")
}

#' Get the available embeddings
#'
#' @description
#' Returns the available embedding as names from the class. This function is
#' used for the single cell-related classes and methods.
#'
#' @param x An object to get embedding from
#' @param ... Other parameters.
#'
#' @return Get the names of the available embeddings.
#'
#' @export
get_available_embeddings <- function(x, ...) {
  UseMethod("get_available_embeddings")
}

#' Get the sNN graph
#'
#' @description
#' Returns the shared nearest neighbour graph from the object. This function is
#' used for the single cell-related classes and methods.
#'
#' @param x An object to get the sNN graph from.
#' @param ... Other parameters.
#'
#' @returns The igraph that has the shared nearest neighbours.
#'
#' @export
get_snn_graph <- function(x, ...) {
  UseMethod("get_snn_graph")
}

#' Get the KNN object
#'
#' @description
#' Returns the `SingleCellNearestNeighbour` from the object. This function is
#' used for the single cell-related classes and methods.
#'
#' @param x An object to get the KNN class from.
#' @param ... Other parameters.
#'
#' @returns The `SingleCellNearestNeighbour` object.
#'
#' @export
get_knn_obj <- function(x, ...) {
  UseMethod("get_knn_obj")
}

#' Get the MAGIC imputed layer
#'
#' @description
#' Returns the `ScMagic` layer written by [bixverse::run_magic_sc()]. This
#' function is used for the single cell-related classes and methods.
#'
#' @param x An object to get the imputed layer from.
#' @param ... Other parameters.
#'
#' @returns The `ScMagic` object, or `NULL` when nothing was imputed.
#'
#' @export
get_magic <- function(x, ...) {
  UseMethod("get_magic")
}

### others ---------------------------------------------------------------------

#### obs -----------------------------------------------------------------------

#' Get the ready obs data from various sub method
#'
#' @description
#' Helper method that creates data.tables with cell indices which were used
#' in the given analysis + the values that are to be added to the obs table
#' in the DuckDB.
#'
#' @param x An object to get the data from.
#' @param columns Optional string. For some of the functions you can decide
#' to only extract specific columns.
#' @param ... Other parameters
#'
#' @returns Returns a data.table with a cell_idx column for the cells included
#' in the analysis and additional columns to be added to the obs table.
#'
#' @export
get_data <- function(x, columns = NULL, ...) {
  UseMethod("get_data")
}

#### knn getter methods --------------------------------------------------------

#' Get the KNN matrix
#'
#' @description
#' Getter for an integer matrix of samples x neighbours.
#'
#' @param x An object to get the kNN matrix from.
#' @param ... Other parameters.
#'
#' @export
get_knn_mat <- function(x, ...) {
  UseMethod("get_knn_mat")
}

#' Get the KNN distance
#'
#' @description
#' Getter for an integer matrix of samples x distances. Useful in combination
#' with [get_knn_mat()].
#'
#' @param x An object to get the kNN distances from.
#' @param ... Other parameters.
#'
#' @export
get_knn_dist <- function(x, ...) {
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
#' @param object `SingleCells`, `MetaCells` (or potentially other) class.
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
#' @param streaming Optional Boolean. Shall the data be streamed in. Useful for
#' larger data sets where you wish to avoid loading in the whole data. If
#' `NULL`, will automatically detect. Not used for `MetaCells`.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
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
    streaming = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' Identify HVGs without mutating object state
#'
#' @description
#' Like [find_hvg_sc()] but does not mutate `object`. Returns a data.table
#' with per-gene HVG statistics plus `is_hvg`/`hvg_rank` for the top `hvg_no`
#' genes. Useful for computing HVGs on a subset of cells (e.g. a specific
#' cell type) for downstream methods like NMF, without overwriting the HVGs
#' stored on the object.
#'
#' @param object `SingleCells` or `MetaCells` class.
#' @param cell_ids Optional character. Cell ids (or meta cell ids) to restrict
#' the HVG calculation to. If `NULL`, uses [get_cells_to_keep()] for
#' `SingleCells` and all meta cells for `MetaCells`.
#' @param hvg_no Integer. Number of top HVGs to flag. Defaults to `3000L`.
#' @param hvg_params List, see [params_sc_hvg()].
#' @param streaming Optional Boolean. Stream the data. Ignored for `MetaCells`.
#' @param .verbose Boolean or integer. Verbosity.
#'
#' @return data.table with `gene_idx`, `gene_id`, the HVG statistics returned
#' by the Rust HVG function, an `is_hvg` boolean and an `hvg_rank` integer
#' (`NA` for non-HVGs).
#'
#' @export
get_hvg_data_sc <- S7::new_generic(
  name = "get_hvg_data_sc",
  dispatch_args = "object",
  fun = function(
    object,
    cell_ids = NULL,
    hvg_no = 3000L,
    hvg_params = params_sc_hvg(),
    streaming = NULL,
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
#' @param object `SingleCells`, `MetaCells` (or potentially other) class.
#' @param no_pcs Integer. Number of PCs to calculate.
#' @param pca_params Named list. Controls the parameters to be used for the
#' PCA calculation which is single cell-specific, see [params_sc_pca()]
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
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
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
    pca_params = params_sc_pca(),
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
#'   \item `kmknn` - An exact kNN search that leverages k-means clustering under
#'   the hood to prune out data points. The default setting.
#'   \item `exhaustive` - An exhaustive, flat index. On smaller data sets often
#'   faster than the approximate nearest neighbour search algorithms.
#'   \item `hnsw` - Hierarchical Navigable Small World. A graph-based
#'   approximate nearest neighbour search algorithm; works well on large data
#'   sets. A benign race condition is leveraged during index build, making the
#'   build non-deterministic. Bigger impact on smaller data sets.
#'   \item `nndescent` - Nearest neighbour descent. Leverages concepts from
#'   `PyNNDescent` and works well on very large data sets similar to `hnsw`.
#'   \item `ivf` - Inverted file index. Uses first k-means clustering to
#'   identify Voronoi cells and leverages these during querying. Works well
#'   on large data sets with high dimensionality and when you need to return
#'   large number of neighbours.
#'   \item `annoy` - Approximate nearest neighbours Oh Yeah. Tree-based index,
#'   used across different R single cell packages (Seurat, SCE). This version
#'   is purely memory-based.
#' }
#' Subsequently, the kNN graph will be additionally transformed into a shared
#' nearest neighbour graph for clustering methods.
#'
#' @param object `SingleCells`, `MetaCells` (or potentially other) class.
#' @param embd_to_use String. The embedding to use. Whichever you chose, it
#' needs to be part of the object.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param modality String. One of `c("rna", "adt")`. You can only use `"adt"`
#' on `SingleCellsMultiModal` class.
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
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
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
    modality = c("rna", "adt"),
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
#' @param object `SingleCells`, `MetaCells` (or potentially other) class.
#' @param cluster_algorithm String. One of `c("leiden", "louvain")`.
#' @param res Numeric. The resolution parameter for [igraph::cluster_leiden()]
#' or [igraph::cluster_louvain()].
#' @param name String. The name to add to the obs table in the DuckDB.
#' @param modality String. On which modality to run the UMAP. One of
#' `c("rna", "adt", "wnn")`. The two latter options are only available for
#' multi-modal versions with the added data.
#' @param seed Integer. For reproducibility.
#'
#' @return The object with added clustering in the obs table.
#'
#' @export
find_clusters_sc <- S7::new_generic(
  name = "find_clusters_sc",
  dispatch_args = "object",
  fun = function(
    object,
    cluster_algorithm = c("leiden", "louvain"),
    res = 1.0,
    name = "leiden_clustering",
    modality = c("rna", "adt", "wnn"),
    seed = 42L
  ) {
    S7::S7_dispatch()
  }
)

#### auc -----------------------------------------------------------------------

#' Calculate AUC scores (akin to AUCell)
#'
#' @description
#' Calculates an AUC-type score akin to AUCell across the gene sets, see Aibar
#' et al. Three statistics are on offer, all consuming the same within-cell
#' ranking but weighting it differently. `"wilcox"` (the default) is the AUC
#' derived from the Mann-Whitney U statistic over the full ranking; its null
#' sits at 0.5 for any gene set size, which makes it a good fit for pathway
#' activity. `"recovery"` is the recovery-curve AUC under a rank cutoff, i.e.
#' the actual AUCell statistic of Aibar, et al., and is top-heavy: only genes
#' inside the top `max_rank` of the cell contribute. `"ap"` is average
#' precision, the most top-heavy of the three, but its null tracks the gene set
#' prevalence, so raw values are not comparable across gene sets of different
#' size unless `standardise` is on. Data can be streamed in chunks of 50k cells
#' per or loaded in in one go.
#'
#' @param object `SingleCells`, `MetaCells` (or potentially other) class.
#' @param gs_list Named list. The elements have the gene identifiers of the
#' respective gene sets.
#' @param aucell_params List with the AUCell parameters, see
#' [bixverse::params_sc_aucell()] with the following elements:
#' \itemize{
#'   \item auc_type - String. Which statistic to calculate. One of
#'   `c("recovery", "wilcox", "ap")`. `"recovery"` is the SCENIC one.
#'   \item max_rank - Optional numeric. Rank cutoff for `"recovery"`. If `NULL`,
#'   the top 5% of the gene universe is used. Ignored by the other statistics.
#'   \item standardise - Boolean. Shall each gene set's scores be z-scored
#'   across the cells.
#' }
#' @param streaming Optional Boolean. Shall the data be streamed in. Useful for
#' larger data sets where you wish to avoid loading in the whole data. If
#' `NULL`, will automatically detect. Ignored when applied to `MetaCells`.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @return AUCell results in form of a matrix that is cells x gene sets or as
#' `ScMatrixRes` pending the input.
#'
#' @export
#'
#' @references Aibar, et al., Nat Methods, 2017
aucell_sc <- S7::new_generic(
  name = "aucell_sc",
  dispatch_args = "object",
  fun = function(
    object,
    gs_list,
    aucell_params = params_sc_aucell(),
    streaming = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#### scenic --------------------------------------------------------------------

##### gene filtering -----------------------------------------------------------

#' Filter genes for SCENIC GRN inference
#'
#' @description
#' Filters genes by minimum total counts and minimum expressed-cell fraction
#' using the SCENIC inclusion criteria. Returns a character vector of gene
#' identifiers passing both filters.
#'
#' @param object `SingleCells`, `MetaCells` (or potentially other) class.
#' @param scenic_params List. SCENIC parameters, see
#' [bixverse::params_scenic()]. Only `min_counts` and `min_cells` are used
#' by this function.
#' @param cells_to_take Optional string vector. Cell identifiers to restrict
#' to. If `NULL`, defaults to all filtered cells in the class.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns A character vector of gene identifiers passing the SCENIC
#' inclusion criteria.
#'
#' @export
scenic_gene_filter_sc <- S7::new_generic(
  name = "scenic_gene_filter_sc",
  dispatch_args = "object",
  fun = function(
    object,
    scenic_params = params_scenic(),
    cells_to_take = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

##### grn ----------------------------------------------------------------------

#' Run SCENIC GRN inference
#'
#' @description
#' Runs SCENIC GRN inference on the provided genes using the specified
#' transcription factors as predictors. Returns a `ScenicGrn` object
#' containing the TF-gene importance matrix for further processing.
#'
#' @param object `SingleCells`, `MetaCells` (or potentially other) class.
#' @param tf_ids Character vector. Transcription factor gene identifiers to
#' use as predictors. Must be a subset of gene identifiers present in the
#' object.
#' @param scenic_params List. SCENIC parameters, see
#' [bixverse::params_scenic()].
#' @param genes_to_take Optional character vector. Target gene identifiers.
#' If `NULL`, genes are selected automatically via
#' [bixverse::scenic_gene_filter_sc()] using the `min_counts` and `min_cells`
#' thresholds in `scenic_params`.
#' @param cells_to_take Optional string vector. Cell identifiers to restrict
#' to. If `NULL`, defaults to all filtered cells in the class.
#' @param streaming Optional Boolean. Shall the data be streamed in. Useful for
#' larger data sets where you wish to avoid loading in the whole data. If
#' `NULL`, will automatically detect. Ignored when applied to `MetaCells`.
#' @param random_seed Integer. Used for reproducibility. Defaults to `42L`.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns A `ScenicGrn` object.
#'
#' @details
#' TF identifiers that are not found in the object's gene list are silently
#' dropped with a warning indicating how many were removed. TF indices are
#' intersected with the target gene indices so that TFs not passing the gene
#' filter are excluded from the predictor set but remain as potential targets
#' if present in `genes_to_take`. You have the option to generate the TF-gene
#' importance values with three distinct methods. For the `random_forest` and
#' the `extratrees` version, a batching strategy is applied in the default
#' settings. Correlated genes are identified and clustered together via
#' k-means clustering on the feature loadings of the PCA. These are then
#' divided into batches of `gene_batch_size` and the regression learners
#' are leveraging multi-target regression to fit all genes in the batch in one
#' go. This massively accelerates the algorithm and the importance values per
#' gene-TF pair are calculated then individually. Due to the batching by
#' similar gene, the signal dilution is limited. If you wish to run the
#' traditional approach, you can set gene_batch_size to `1L` or use the
#' `grnboost2` learner that can only fit one gene at a given time.
#'
#' @export
scenic_grn_sc <- S7::new_generic(
  name = "scenic_grn_sc",
  dispatch_args = "object",
  fun = function(
    object,
    tf_ids,
    scenic_params = params_scenic(),
    genes_to_take = NULL,
    cells_to_take = NULL,
    streaming = NULL,
    random_seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#### nmf -----------------------------------------------------------------------

#' Run single-run NMF on single cell or meta cell data
#'
#' @description
#' Runs a single HALS NMF on a chosen subset of cells and genes. For
#' `SingleCells`, the counts are streamed from disk via the Rust binary
#' files; for `MetaCells`, the in-memory sparse counts are used.
#'
#' @param object `SingleCells` or `MetaCells` class.
#' @param k Integer. Number of latent factors to return.
#' @param cell_ids Optional character. Cell ids (or meta cell ids) to restrict
#' the NMF to. If `NULL`, uses [get_cells_to_keep()] for `SingleCells` and all
#' meta cells for `MetaCells`.
#' @param gene_ids Optional character. Gene ids to restrict the NMF to. If
#' `NULL`, uses [get_hvg()] on the object.
#' @param preprocessing String. One of `c("none", "sd", "sqrt_sd")`.
#' @param use_second_layer Boolean. If `TRUE`, runs NMF on the normalised
#' counts (recommended); if `FALSE`, on the raw counts.
#' @param nmf_hals_params List, see [params_nmf_hals()].
#' @param seed Integer. Random seed for initialisation.
#' @param .verbose Boolean or integer. Verbosity.
#'
#' @returns An `NmfResult` object.
#'
#' @export
nmf_sc <- S7::new_generic(
  name = "nmf_sc",
  dispatch_args = "object",
  fun = function(
    object,
    k,
    cell_ids = NULL,
    gene_ids = NULL,
    preprocessing = "none",
    use_second_layer = TRUE,
    nmf_hals_params = params_nmf_hals(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' Run stabilised (multi-run) NMF on single cell or meta cell data
#'
#' @description
#' Runs `n_runs` HALS NMF with random initialisations seeded by `seed + i`.
#' The `nmf_init` field in `nmf_hals_params` is ignored; random init is
#' always used.
#'
#' @inheritParams nmf_sc
#' @param n_runs Integer. Number of random restarts.
#'
#' @returns A `StabilisedNmfResult` object.
#'
#' @export
stabilised_nmf_sc <- S7::new_generic(
  name = "stabilised_nmf_sc",
  dispatch_args = "object",
  fun = function(
    object,
    k,
    cell_ids = NULL,
    gene_ids = NULL,
    preprocessing = "none",
    use_second_layer = TRUE,
    nmf_hals_params = params_nmf_hals(),
    n_runs = 30L,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' Run consensus NMF on single cell or meta cell data
#'
#' @description
#' Runs `n_runs` HALS NMF restarts, pools their components, drops unstable ones
#' by local density, k-means clusters the survivors into `k` groups and refits
#' the partner factor against the per-cluster median. This is cNMF: the answer
#' is the structure the restarts agree on, and the mean silhouette of those
#' clusters (`stability`) says how much they agreed.
#'
#' Prefer this over [stabilised_nmf_sc()], which picks the lowest-loss restart
#' and tells you nothing about whether that restart is reproducible.
#'
#' @details
#' Use [nmf_k_sweep_sc()] first if you do not already know `k`.
#'
#' Memory is the thing to watch. The restarts are dense and all live at once:
#' roughly `n_cells * k * n_runs` plus `n_runs * k * n_genes` floats on top of
#' the counts. At 200k cells, `k = 20` and `n_runs = 50` that is a few hundred
#' megabytes before anything else, which is exactly the regime where running on
#' `MetaCells` instead is the honest answer.
#'
#' The density filter is the part that bites. If it leaves fewer than `k`
#' components the fit errors rather than returning something degenerate. With
#' few restarts the filter is jumpy, so either raise `n_runs` or set
#' `density_threshold = 2` in [params_nmf_consensus()] to switch it off.
#'
#' @param object `SingleCells` or `MetaCells` class.
#' @param k Integer. Number of latent factors. Must be at least 2.
#' @param cell_ids Optional character. Cell ids (or meta cell ids) to restrict
#' the NMF to. If `NULL`, uses [get_cells_to_keep()] for `SingleCells` and all
#' meta cells for `MetaCells`.
#' @param gene_ids Optional character. Gene ids to restrict the NMF to. If
#' `NULL`, uses [get_hvg()] on the object.
#' @param preprocessing String. One of `c("none", "sd", "sqrt_sd")`.
#' @param use_second_layer Boolean. If `TRUE`, runs NMF on the normalised
#' counts (recommended); if `FALSE`, on the raw counts.
#' @param nmf_hals_params List, see [params_nmf_hals()]. The `nmf_init` field is
#' ignored, restarts always use random initialisation.
#' @param nmf_consensus_params List, see [params_nmf_consensus()].
#' @param n_runs Integer. Number of random restarts. Must be at least 2.
#' @param seed Integer. Base random seed. Restart `i` uses `seed + i`, and the
#' k-means step is seeded from it too.
#' @param .verbose Boolean or integer. Verbosity.
#'
#' @returns A `ConsensusNmfResult` object.
#'
#' @references Kotliar et al., eLife, 2019
#'
#' @export
consensus_nmf_sc <- S7::new_generic(
  name = "consensus_nmf_sc",
  dispatch_args = "object",
  fun = function(
    object,
    k,
    cell_ids = NULL,
    gene_ids = NULL,
    preprocessing = "none",
    use_second_layer = TRUE,
    nmf_hals_params = params_nmf_hals(),
    nmf_consensus_params = params_nmf_consensus(),
    n_runs = 30L,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' Sweep k for consensus NMF on single cell or meta cell data
#'
#' @description
#' Runs the consensus clustering step across a range of `k` and reports
#' stability against reconstruction error, without keeping any factors. Pick
#' the `k` where stability is still high and the error curve has not yet
#' flattened out, then run [consensus_nmf_sc()] there.
#'
#' @details
#' This is a diagnostic, so it leaves the object alone and hands the result back
#' directly. `plot()` on it gives you the two curves.
#'
#' Cost is `length(k_range) * n_runs` full NMF fits. On the `SingleCells` path
#' the counts are read from disk once and reused across every `k`, but the fits
#' themselves are not free, so keep both modest on a first pass.
#'
#' @inheritParams consensus_nmf_sc
#' @param k_range Integer vector. The ranks to evaluate. Every entry must be at
#' least 2.
#'
#' @returns An `NmfKSweepResult`, which is a data.table with one row per `k`.
#'
#' @references Kotliar et al., eLife, 2019
#'
#' @export
nmf_k_sweep_sc <- S7::new_generic(
  name = "nmf_k_sweep_sc",
  dispatch_args = "object",
  fun = function(
    object,
    k_range,
    cell_ids = NULL,
    gene_ids = NULL,
    preprocessing = "none",
    use_second_layer = TRUE,
    nmf_hals_params = params_nmf_hals(),
    nmf_consensus_params = params_nmf_consensus(),
    n_runs = 30L,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#### i/o -----------------------------------------------------------------------

#' Save memory-bound data to disk
#'
#' @description
#' Helper function that stores the memory-bound data to disk for checkpointing
#' or when you close the session for quick recovery of prior work. You have the
#' option to save as `".rds"` or `".qs2"` (you need to have the package `"qs2"`
#' installed for this option!).
#'
#' @param object `SingleCells`, `MetaCells` or `SingleCellsMultiModal` class.
#' @param type String. One of `c("qs2", "rds")`. Defines which binary format to
#' use. Will default to `"qs2"` for speed.
#'
#' @returns The object with added information on the data on disk.
#'
#' @export
save_sc_exp_to_disk <- S7::new_generic(
  name = "save_sc_exp_to_disk",
  dispatch_args = "object",
  fun = function(
    object,
    type = c("qs2", "rds")
  ) {
    S7::S7_dispatch()
  }
)

#' Load an existing SingleCells from disk
#'
#' @description
#' Helper function that can load the parameters to access the on-disk stored
#' data into the class.
#'
#' @param object `SingleCells`, `MetaCells` or `SingleCellsMultiModal` class.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns The object with added information on the data on disk.
#'
#' @export
load_existing <- S7::new_generic(
  name = "load_existing",
  dispatch_args = "object",
  fun = function(
    object,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#### plotting ------------------------------------------------------------------

#' Extract grouped gene statistics for dot plots
#'
#' @description
#' Extracts per-group mean expression and percentage of expressing cells for a
#' set of genes. Returns a long-format data.table suitable for dot plots.
#'
#' @param object A single cell class.
#' @param features Character vector. Gene IDs to extract.
#' @param grouping_variable String. Column name in the obs table to group by.
#' @param scale_exp Boolean. Whether to min-max scale mean expression per gene.
#' @param modality String. One of `c("rna", "adt")`. ADT is only available for
#' `SingleCellsMultiModal`.
#'
#' @return A data.table with columns: gene, group, mean_exp, scaled_exp, pct_exp.
#'
#' @export
extract_dot_plot_data <- S7::new_generic(
  name = "extract_dot_plot_data",
  dispatch_args = "object",
  fun = function(
    object,
    features,
    grouping_variable,
    scale_exp = TRUE,
    modality = c("rna", "adt")
  ) {
    S7::S7_dispatch()
  }
)

#' Extract normalised gene expression for plotting
#'
#' @description
#' Extracts dense normalised (log1p) expression values for a set of genes,
#' optionally with additional observation metadata columns.
#'
#' @param object A single cell class.
#' @param features Character vector. Gene IDs to extract.
#' @param obs_cols Optional character vector. Column names from the obs table
#' to include.
#' @param scale Boolean. Whether to z-score the expression values.
#' @param clip Optional numeric. If `scale = TRUE`, clip z-scores to
#' `[-clip, clip]`.
#' @param modality String. One of `c("rna", "adt")`. ADT is only available for
#' `SingleCellsMultiModal`.
#' @param layer String. One of `c("norm", "magic")`. With `"magic"` the values
#' come from the imputed layer [bixverse::run_magic_sc()] wrote, which only
#' holds the genes it was asked for. Imputation inflates gene-gene correlation,
#' so this is for looking at things, not for measuring them. Note that
#' [bixverse::extract_dot_plot_data()] deliberately has no such argument:
#' group means of imputed values are exactly the quantity MAGIC manufactures.
#'
#' @return A data.table with a `cell_id` column, one column per gene, and
#' any requested obs columns.
#'
#' @export
extract_gene_expression <- S7::new_generic(
  name = "extract_gene_expression",
  dispatch_args = "object",
  fun = function(
    object,
    features,
    obs_cols = NULL,
    scale = FALSE,
    clip = NULL,
    modality = c("rna", "adt"),
    layer = c("norm", "magic")
  ) {
    S7::S7_dispatch()
  }
)

#### hotspot -------------------------------------------------------------------

#' Calculate the local auto-correlation of a gene
#'
#' @description
#' This method implements the HotSpot approach (see DeTomaso, et al.) to
#' calculate the auto-correlation of a given gene in the kNN graph based on
#' the chosen embedding. This can be used to identify genes that have strong
#' local correlations and vary across the kNN graph.
#'
#' @param object `SingleCells` or `MetaCells` class.
#' @param embd_to_use String. The embedding to use. Defaults to `"pca"`.
#' @param use_knn Boolean. Shall the internal kNN be used. If set to yes, you
#' need to ensure consistency. If you provide `cells_to_take`, the function
#' will regenerate the kNN graph with these cells.
#' @param hotspot_params List with hotspot parameters, see
#' [bixverse::params_sc_hotspot()] with the following elements:
#' \itemize{
#'   \item model - String. Which of the available models to use for the
#'   gene expression. Choices are one of `c("danb", "normal", "bernoulli")`.
#'   \item normalise - Boolean. Shall the data be normalised.
#'   \item weighted_graph - Boolean. Shall the Gaussian kernel be applied to
#'   the neighbour distances. If `FALSE`, every retained edge weighs one.
#'   \item neighborhood_factor - Float. Kernel width for `weighted_graph`.
#'   \item knn - List of kNN parameters. See [bixverse::params_knn_defaults()]
#'   for available parameters and their defaults.
#' }
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param cells_to_take Optional string vector. If you want to only use
#' selected cells. If `NULL` will default to all cells_to_keep in the class.
#' @param genes_to_take Optional string vector. If you wish to limit the
#' search to a subset of genes. If `NULL` will default to all genes in the
#' class.
#' @param streaming Optional Boolean. Shall the data be streamed in. Useful for
#' larger data sets where you wish to avoid loading in the whole data. If
#' `NULL`, will automatically detect. Ignored for `MetaCells`, which are held
#' in memory.
#' @param random_seed Integer. Used for reproducibility.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns A data.table with the auto-correlations on a per gene basis and
#' various statistics.
#'
#' @details
#' Should a gene not be found in sufficient cells, the gene will be
#' automatically filtered out from the results. This can occur for example
#' if you have filtered out the cells that contain a given gene. The underlying
#' genes are still available, but the cells that might contain them are not
#' included.
#'
#' Whether the neighbour distances need squaring before the kernel sees them
#' follows from the metric. With `use_knn = TRUE` it is taken from the metric
#' stored on the cached kNN graph, otherwise from `ann_dist` in
#' `hotspot_params`.
#'
#' @export
hotspot_autocor_sc <- S7::new_generic(
  name = "hotspot_autocor_sc",
  dispatch_args = "object",
  fun = function(
    object,
    embd_to_use = "pca",
    use_knn = TRUE,
    hotspot_params = params_sc_hotspot(),
    no_embd_to_use = NULL,
    cells_to_take = NULL,
    genes_to_take = NULL,
    streaming = NULL,
    random_seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' Calculate the local pairwise gene-gene correlation
#'
#' @description
#' This method implements the HotSpot approach (see DeTomaso, et al.) to
#' calculate the local gene-gene correlations and their Z-scores.
#'
#' @param object `SingleCells` or `MetaCells` class.
#' @param embd_to_use String. The embedding to use. Defaults to `"pca"`.
#' @param use_knn Boolean. Shall the internal kNN be used. If set to yes, you
#' need to ensure consistency. If you provide `cells_to_take`, the function
#' will regenerate the kNN graph with these cells.
#' @param hotspot_params List with hotspot parameters, see
#' [bixverse::params_sc_hotspot()] with the following elements:
#' \itemize{
#'   \item model - String. Which of the available models to use for the
#'   gene expression. Choices are one of `c("danb", "normal", "bernoulli")`.
#'   \item normalise - Boolean. Shall the data be normalised.
#'   \item weighted_graph - Boolean. Shall the Gaussian kernel be applied to
#'   the neighbour distances. If `FALSE`, every retained edge weighs one.
#'   \item neighborhood_factor - Float. Kernel width for `weighted_graph`.
#'   \item knn - List of kNN parameters. See [bixverse::params_knn_defaults()]
#'   for available parameters and their defaults.
#' }
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param cells_to_take Optional string vector. If you want to only use
#' selected cells. If `NULL` will default to all cells_to_keep in the class.
#' @param genes_to_take Optional string vector. If you wish to limit the
#' search to a subset of genes. If `NULL` will default to all genes in the
#' class.
#' @param streaming Optional Boolean. Shall the data be streamed in. Useful for
#' larger data sets where you wish to avoid loading in the whole data. If
#' `NULL`, will automatically detect. Ignored for `MetaCells`, which are held
#' in memory.
#' @param working_mem_gb Numeric. Approximate working memory (GB) the streaming
#'  pair path may use for resident gene panels. Ignored when `streaming` is
#'  `FALSE`. Larger values mean fewer disk re-reads. Note this excludes the two
#' dense N_genes x N_genes output matrices, which scale with `genes_to_use`.
#' Defaults to `4` (4 GB of memory allocated).
#' @param random_seed Integer. Used for reproducibility.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns A `sc_hotspot` class that can be used for subsequent analysis.
#'
#' @details
#' Should a gene not be found in sufficient cells, the pairs with this gene
#' will be set to 0. Please ensure prior to running the function that you
#' are only calculating gene-gene auto-correlations that occur in sufficient
#' cells.
#'
#' Whether the neighbour distances need squaring before the kernel sees them
#' follows from the metric. With `use_knn = TRUE` it is taken from the metric
#' stored on the cached kNN graph, otherwise from `ann_dist` in
#' `hotspot_params`.
#'
#' @export
hotspot_gene_cor_sc <- S7::new_generic(
  name = "hotspot_gene_cor_sc",
  dispatch_args = "object",
  fun = function(
    object,
    embd_to_use = "pca",
    use_knn = TRUE,
    hotspot_params = params_sc_hotspot(),
    no_embd_to_use = NULL,
    cells_to_take = NULL,
    genes_to_take = NULL,
    streaming = NULL,
    working_mem_gb = 4,
    random_seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#### vision -------------------------------------------------------------------

#' Calculate VISION scores
#'
#' @description
#' Calculates an VISION-type scores for pathways based on DeTomaso, et al.
#' Compared to other score types, you can also calculate delta-type scores
#' between positive and negative gene indices, think epithelial vs mesenchymal
#' gene signature, etc.
#'
#' @param object `SingleCells`, `MetaCells` (or potentially other) class.
#' @param gs_list Named nested list. The elements have the gene identifiers of
#' the respective gene sets and have the option to have a `"pos"` and `"neg"`
#' gene sets. The names need to be part of the variables of the object.
#' @param streaming Optional Boolean. Shall the data be streamed in. Useful for
#' larger data sets where you wish to avoid loading in the whole data. If
#' `NULL`, will automatically detect. Ignored when applied to `MetaCells`.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns The VISION scores in form of a matrix that is cells x gene sets
#' or as `ScMatrixRes` pending the input.
#'
#' @references DeTomaso, et al., Nat. Commun., 2019
#'
#' @export
vision_sc <- S7::new_generic(
  name = "vision_sc",
  dispatch_args = "object",
  fun = function(
    object,
    gs_list,
    streaming = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' Calculate VISION scores (with auto-correlation scores)
#'
#' @description
#' Calculates VISION-type scores for pathways based on DeTomaso, et al.
#' Compared to other score types, you can also calculate delta-type scores
#' between positive and negative gene indices, think epithelial vs mesenchymal
#' gene signature, etc. Additionally, this function also calculates the auto-
#' correlation values, answering the question if a given signature shows non-
#' random enrichment on the kNN graph. The kNN graph (and distance measures)
#' will be generated on-the-fly based on the embedding you wish to use.
#'
#' @param object `SingleCells`, `MetaCells` (or potentially other) class.
#' @param gs_list Named nested list. The elements have the gene identifiers of
#' the respective gene sets and have the option to have a `"pos"` and `"neg"`
#' gene sets. The names need to be part of the variables of the object.
#' @param vision_params List with vision parameters, see
#' [bixverse::params_sc_vision()] with the following elements:
#' \itemize{
#'   \item n_perm - Integer. Number of random permutations
#'   \item n_cluster - Integer. Number of random clusters to generate to
#'   associate each set with.
#'   \item knn - List of kNN parameters. See [bixverse::params_knn_defaults()]
#'   for available parameters and their defaults.
#' }
#' @param embd_to_use String. The embedding to use. Whichever you chose, it
#' needs to be part of the object.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param use_knn Boolean. Shall the internal kNN be used. If set to yes, you
#' need to ensure consistency.
#' @param random_seed Integer. The random seed.
#' @param streaming Optional Boolean. Shall the data be streamed in. Useful for
#' larger data sets where you wish to avoid loading in the whole data. If
#' `NULL`, will automatically detect. Ignored when applied to `MetaCells`.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item vision_matrix - Matrix of cells x signatures with the VISION
#'   pathway scores as values.
#'   \item auto_cor_dt - data.table with the auto-correlation results per
#'   gene set, i.e., `auto_cor` (1 - Gaery's C), `p_val` and `fdr`.
#' }
#'
#' @references DeTomaso, et al., Nat. Commun., 2019
#'
#' @export
vision_w_autocor_sc <- S7::new_generic(
  name = "vision_w_autocor_sc",
  dispatch_args = "object",
  fun = function(
    object,
    gs_list,
    embd_to_use,
    no_embd_to_use = NULL,
    use_knn = TRUE,
    vision_params = params_sc_vision(),
    streaming = NULL,
    random_seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)
