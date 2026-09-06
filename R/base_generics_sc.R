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
#' @returns The obs table
#'
#' @export
#'
#' @examples
#' # the obs table, restricted to a few columns
#' sc <- demo_single_cells(prepped = FALSE)
#' head(get_sc_obs(sc, cols = c("cell_id", "cell_grp", "lib_size")), 3)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#' @returns The vars table
#'
#' @export
#'
#' @examples
#' # the per gene statistics the HVG step wrote
#' sc <- demo_single_cells()
#' head(get_sc_var(sc, cols = c("gene_id", "mean", "var_std")), 3)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#' @returns The counts table
#'
#' @export
#'
#' @examples
#' # raw counts for the first ten genes, cell-centric (CSR)
#' sc <- demo_single_cells(prepped = FALSE)
#' counts <- get_sc_counts(sc, gene_indices = 1:10, .verbose = FALSE)
#' dim(counts)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#' @returns A data.table with available features.
#'
#' @export
#'
#' @examples
#' # what can be queried from the obs table and the counts
#' sc <- demo_single_cells(prepped = FALSE)
#' head(get_sc_available_features(sc), 3)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#' @returns Invisible self
#'
#' @export
#'
#' @examples
#' # rename a column in the obs table
#' sc <- demo_single_cells(prepped = FALSE)
#' sc <- setnames_sc(sc, table = "obs", old = "cell_grp", new = "cell_type")
#' head(get_sc_obs(sc)$cell_type, 3)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#' @returns Invisible self.
#'
#' @export
#'
#' @examples
#' # drop a column that is no longer needed
#' sc <- demo_single_cells(prepped = FALSE)
#' sc <- drop_cols_sc(sc, table = "obs", cols = "batch_index")
#' names(get_sc_obs(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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

#' Set gene mapping
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
#'
#' @examples
#' # the mapping is normally written during ingestion
#' sc <- demo_single_cells(prepped = FALSE)
#' genes <- get_gene_names(sc)
#' sc <- set_gene_mapping(sc, stats::setNames(seq_along(genes), genes))
#' head(get_gene_names(sc), 3)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
set_gene_mapping <- function(x, gene_map) {
  UseMethod("set_gene_mapping")
}

#' Set cell mapping
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
#'
#' @examples
#' # the mapping is normally written during ingestion
#' sc <- demo_single_cells(prepped = FALSE)
#' cells <- get_cell_names(sc)
#' sc <- set_cell_mapping(sc, stats::setNames(seq_along(cells), cells))
#' head(get_cell_names(sc), 3)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
set_cell_mapping <- function(x, cell_map) {
  UseMethod("set_cell_mapping")
}

#' Set cells to keep
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
#'
#' @examples
#' # restrict everything downstream to the first 100 cells
#' sc <- demo_single_cells(prepped = FALSE)
#' sc <- set_cells_to_keep(sc, get_cell_names(sc)[1:100])
#' length(get_cells_to_keep(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
set_cells_to_keep <- function(x, cells_to_keep) {
  UseMethod("set_cells_to_keep")
}

#' Reset the cells to keep
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
#'
#' @examples
#' # a filter taken back off, cache wiped with it
#' sc <- demo_single_cells(prepped = FALSE)
#' sc <- set_cells_to_keep(sc, get_cell_names(sc)[1:100])
#' sc <- reset_cells_to_keep(sc, force = TRUE)
#' length(get_cells_to_keep(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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

#' Set the HVG genes
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
#'
#' @examples
#' # HVGs picked by hand instead of by find_hvg_sc()
#' sc <- demo_single_cells(prepped = FALSE)
#' sc <- set_hvg(sc, get_gene_names(sc)[1:20])
#' head(get_hvg(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
set_hvg <- function(x, hvg) {
  UseMethod("set_hvg")
}

#### getters -------------------------------------------------------------------

#' Get the HVG
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
#'
#' @examples
#' # stored 0-based, so map them back through the gene names
#' sc <- demo_single_cells()
#' get_gene_names_from_idx(sc, head(get_hvg(sc), 3))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
get_hvg <- function(x) {
  UseMethod("get_hvg")
}

#' Get the gene names
#'
#' @description
#' Get the main gene names (for example symbols or Ensembl identifiers).
#'
#' @param x An object to get the gene names from.
#'
#' @returns The primary gene identifiers stored in the class.
#'
#' @export
#'
#' @examples
#' # the primary gene identifiers held by the object
#' sc <- demo_single_cells(prepped = FALSE)
#' head(get_gene_names(sc), 3)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
get_gene_names <- function(x) {
  UseMethod("get_gene_names")
}

#' Get the cell names
#'
#' @description
#' Returns the cell names (usually barcodes).
#'
#' @param x An object to get the cell names from.
#' @param filtered Boolean. Shall, if found only the cell names of the
#' `cells_to_keep` be returned (see [bixverse::set_cells_to_keep()]. Defaults
#' to `FALSE`
#'
#' @returns The cell names (barcodes)
#'
#' @export
#'
#' @examples
#' # barcodes of the cells that passed quality control
#' sc <- demo_single_cells(prepped = FALSE)
#' head(get_cell_names(sc, filtered = TRUE), 3)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
get_cell_names <- function(x, filtered = FALSE) {
  UseMethod("get_cell_names")
}

#' Get the index position for a gene
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
#'
#' @examples
#' # R-based positions of two genes
#' sc <- demo_single_cells(prepped = FALSE)
#' get_gene_indices(sc, c("gene_01", "gene_02"), rust_index = FALSE)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
get_gene_indices <- function(x, gene_ids, rust_index) {
  UseMethod("get_gene_indices")
}

#' Get the index position for a gene
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
#'
#' @examples
#' # R-based positions of two cells
#' sc <- demo_single_cells(prepped = FALSE)
#' get_cell_indices(sc, c("cell_001", "cell_002"), rust_index = FALSE)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
get_cell_indices <- function(x, cell_ids, rust_index) {
  UseMethod("get_cell_indices")
}

#' Get the cell idx (R-based) and cell names
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
#'
#' @examples
#' # cell indices carrying the barcodes as names
#' sc <- demo_single_cells(prepped = FALSE)
#' head(get_cell_info(sc), 3)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
get_cell_info <- function(x, filtered = TRUE) {
  UseMethod("get_cell_info")
}

#' Get the cells to keep
#'
#' @description
#' Returns the indices of the cells that survived quality control. These are
#' stored 0-indexed for Rust, so add one before using them in R.
#'
#' @param x An object from which to get the cells to keep from. These are
#' 0-indexed.
#'
#' @returns Integer vector with 0-indices of the cells to keep.
#'
#' @export
#'
#' @examples
#' # 0-based for Rust, so add one before indexing in R
#' sc <- demo_single_cells(prepped = FALSE)
#' head(get_cells_to_keep(sc) + 1, 3)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # Rust indices translated back into gene identifiers
#' sc <- demo_single_cells(prepped = FALSE)
#' get_gene_names_from_idx(sc, gene_idx = 0:2, rust_based = TRUE)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # an externally computed embedding pushed into the cache
#' sc <- demo_single_cells(prepped = FALSE)
#' sc <- set_pca_factors(sc, matrix(stats::rnorm(dim(sc)[1] * 2), ncol = 2))
#' dim(get_pca_factors(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # only the first two loading vectors kept
#' sc <- demo_single_cells()
#' sc <- set_pca_loadings(sc, get_pca_loadings(sc)[, 1:2])
#' dim(get_pca_loadings(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # singular values from a decomposition done elsewhere
#' sc <- demo_single_cells(prepped = FALSE)
#' sc <- set_pca_singular_vals(sc, c(4.1, 2.3))
#' get_pca_singular_val(sc)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # the first two PCs stored as an embedding in their own right
#' sc <- demo_single_cells()
#' sc <- set_embedding(sc, get_pca_factors(sc)[, 1:2], name = "pca_2d")
#' get_available_embeddings(sc)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # a kNN built outside the object put into the cache
#' sc <- demo_single_cells()
#' knn <- generate_knn_sc(sc, .validate_index = FALSE, .verbose = FALSE)
#' sc <- set_knn(sc, knn)
#' dim(get_knn_mat(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # the sNN graph taken out and put back
#' sc <- demo_single_cells()
#' snn <- get_snn_graph(sc)
#' sc <- set_snn_graph(remove_snn_graph(sc), snn)
#' igraph::vcount(get_snn_graph(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # drop the cached kNN, for example before rebuilding it
#' sc <- demo_single_cells()
#' sc <- remove_knn(sc)
#' is.null(get_knn_obj(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # drop the cached sNN graph
#' sc <- demo_single_cells()
#' sc <- remove_snn_graph(sc)
#' is.null(get_snn_graph(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
remove_snn_graph <- function(x, ...) {
  UseMethod("remove_snn_graph")
}

#' Set/add the MAGIC imputed layer
#'
#' @param x An object to add the imputed layer to.
#' @param magic `ScMagic` class with the imputed counts.
#' @param ... Other parameters.
#'
#' @returns The object with the imputed layer attached.
#'
#' @export
#'
#' @examples
#' # the imputed layer taken out and put back
#' sc <- demo_single_cells()
#' sc <- run_magic_sc(sc, features = get_gene_names(sc)[1:5], .verbose = FALSE)
#' magic <- get_magic(sc)
#' sc <- set_magic(remove_magic(sc), magic)
#' dim(get_magic(sc)$data)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
set_magic <- function(x, magic, ...) {
  UseMethod("set_magic")
}

#' Remove the MAGIC imputed layer
#'
#' @param x An object from which to remove the imputed layer.
#' @param ... Other parameters.
#'
#' @returns The object with the imputed layer dropped.
#'
#' @export
#'
#' @examples
#' # drop the imputed layer again
#' sc <- demo_single_cells()
#' sc <- run_magic_sc(sc, features = get_gene_names(sc)[1:5], .verbose = FALSE)
#' sc <- remove_magic(sc)
#' is.null(get_magic(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
remove_magic <- function(x, ...) {
  UseMethod("remove_magic")
}

#### getters -------------------------------------------------------------------

#' Get the PCA factors
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
#'
#' @examples
#' # cells x PCs
#' sc <- demo_single_cells()
#' dim(get_pca_factors(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
get_pca_factors <- function(x, ...) {
  UseMethod("get_pca_factors")
}

#' Get the PCA loadings
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
#'
#' @examples
#' # HVGs x PCs
#' sc <- demo_single_cells()
#' dim(get_pca_loadings(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
get_pca_loadings <- function(x, ...) {
  UseMethod("get_pca_loadings")
}

#' Get the PCA singular values
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
#'
#' @examples
#' # singular values, largest first
#' sc <- demo_single_cells()
#' head(get_pca_singular_val(sc), 3)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
get_pca_singular_val <- function(x, ...) {
  UseMethod("get_pca_singular_val")
}

#' Get the embedding
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
#'
#' @examples
#' # any stored embedding, by name
#' sc <- demo_single_cells()
#' dim(get_embedding(sc, "pca"))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#' @returns Get the names of the available embeddings.
#'
#' @export
#'
#' @examples
#' # what is in the cache to plot against
#' sc <- demo_single_cells()
#' get_available_embeddings(sc)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # the sNN graph the clustering methods run on
#' sc <- demo_single_cells()
#' igraph::vcount(get_snn_graph(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # the cached kNN, cells x neighbours
#' sc <- demo_single_cells()
#' dim(get_knn_mat(get_knn_obj(sc)))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # the imputed layer, only the genes MAGIC was asked for
#' sc <- demo_single_cells()
#' sc <- run_magic_sc(sc, features = get_gene_names(sc)[1:5], .verbose = FALSE)
#' dim(get_magic(sc)$data)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # the cell indices and cluster memberships a run produced
#' sc <- demo_single_cells()
#' res <- fast_cluster_sc(sc, resolutions = 1.0, .verbose = FALSE)
#' head(get_data(res), 3)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # cells x neighbours
#' sc <- demo_single_cells()
#' dim(get_knn_mat(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # cells x neighbour distances
#' sc <- demo_single_cells()
#' dim(get_knn_dist(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#' @returns It will add the mean, var, var_exp, var_std of each gene to the
#' the var table.
#'
#' @export
#'
#' @examples
#' # the twenty most variable genes by the vst method
#' sc <- demo_single_cells(prepped = FALSE)
#' sc <- find_hvg_sc(sc, hvg_no = 20L, .verbose = FALSE)
#' length(get_hvg(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#' @returns data.table with `gene_idx`, `gene_id`, the HVG statistics returned
#' by the Rust HVG function, an `is_hvg` boolean and an `hvg_rank` integer
#' (`NA` for non-HVGs).
#'
#' @export
#'
#' @examples
#' # HVG statistics without touching the object
#' sc <- demo_single_cells(prepped = FALSE)
#' dt <- get_hvg_data_sc(sc, hvg_no = 20L, .verbose = FALSE)
#' head(dt[(is_hvg), c("gene_id", "hvg_rank")], 3)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#' @returns The function will add the PCA factors, loadings and singular values
#' to the object cache in memory.
#'
#' @export
#'
#' @examples
#' # PCA on the highly variable genes
#' sc <- demo_single_cells(prepped = FALSE)
#' sc <- find_hvg_sc(sc, hvg_no = 30L, .verbose = FALSE)
#' sc <- calculate_pca_sc(sc, no_pcs = 10L, .verbose = FALSE)
#' dim(get_pca_factors(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'   Set `extract_knn = TRUE` in the kNN parameters to hand back the descent
#'   graph directly instead of beam searching it. That drops the query pass
#'   altogether, so it is much faster, but recall goes down a little.
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
#' @returns The object with added KNN matrix.
#'
#' @export
#'
#' @examples
#' # kNN and the sNN graph on top of the PCA
#' sc <- demo_single_cells(prepped = FALSE)
#' sc <- find_hvg_sc(sc, hvg_no = 30L, .verbose = FALSE)
#' sc <- calculate_pca_sc(sc, no_pcs = 10L, .verbose = FALSE)
#' sc <- find_neighbours_sc(
#'   sc,
#'   neighbours_params = params_sc_neighbours(knn = list(k = 15L)),
#'   .verbose = FALSE
#' )
#' dim(get_knn_mat(sc))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#' @returns The object with added clustering in the obs table.
#'
#' @export
#'
#' @examples
#' # Leiden on the cached sNN graph
#' sc <- demo_single_cells()
#' sc <- find_clusters_sc(sc, res = 1.0, name = "clusters")
#' table(get_sc_obs(sc)$clusters)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#' ranking but weighting it differently. `"recovery"` (default) is the
#' recovery-curve AUC under a rank cutoff, i.e. the AUCell statistic of Aibar,
#' et al., and is top-heavy: only genes inside the top `max_rank` of the cell
#' contribute. `"wilcox"` is the AUC derived from the Mann-Whitney U statistic
#' over the full ranking; its null sits at 0.5 for any gene set size, which
#' makes it a good fit for pathway activity.  `"ap"` is average precision, the
#' most top-heavy of the three, but its null tracks the gene set prevalence, so
#' raw values are not comparable across gene sets of different size unless
#' `standardise` is on.
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
#' @returns AUCell results in form of a matrix that is cells x gene sets or as
#' `ScMatrixRes` pending the input.
#'
#' @export
#'
#' @references Aibar, et al., Nat Methods, 2017
#'
#' @examples
#' # recovery curve AUC for two marker programmes
#' sc <- demo_single_cells()
#' gs_list <- list(
#'   type_1 = get_gene_names(sc)[1:10],
#'   type_2 = get_gene_names(sc)[11:20]
#' )
#' res <- aucell_sc(sc, gs_list = gs_list, .verbose = FALSE)
#' dim(res)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # genes clearing the SCENIC count and prevalence thresholds
#' sc <- demo_single_cells()
#' genes <- scenic_gene_filter_sc(
#'   sc,
#'   scenic_params = params_scenic(min_counts = 100L, min_cells = 0.05),
#'   .verbose = FALSE
#' )
#' length(genes)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # TF to gene importances from a small random forest
#' sc <- demo_single_cells()
#' res <- scenic_grn_sc(
#'   sc,
#'   tf_ids = get_gene_names(sc)[1:5],
#'   scenic_params = params_scenic(
#'     min_counts = 100L,
#'     learner_params = list(n_trees = 20L)
#'   ),
#'   .verbose = FALSE
#' )
#' res
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # three factors on the highly variable genes
#' sc <- demo_single_cells()
#' res <- nmf_sc(sc, k = 3L, .verbose = FALSE)
#' res
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # five random restarts, the best one reported
#' sc <- demo_single_cells()
#' res <- stabilised_nmf_sc(sc, k = 3L, n_runs = 5L, .verbose = FALSE)
#' res
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # ten restarts pooled, density filter off on data this small
#' sc <- demo_single_cells()
#' res <- consensus_nmf_sc(
#'   sc,
#'   k = 3L,
#'   n_runs = 10L,
#'   nmf_consensus_params = params_nmf_consensus(density_threshold = 2),
#'   .verbose = FALSE
#' )
#' res
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @param k_range Integer vector. The ranks to evaluate. Every entry must be at
#' least 2.
#'
#' @returns An `NmfKSweepResult`, which is a data.table with one row per `k`.
#'
#' @references Kotliar et al., eLife, 2019
#'
#' @export
#'
#' @examples
#' # stability against reconstruction error across three ranks
#' sc <- demo_single_cells()
#' res <- nmf_k_sweep_sc(
#'   sc,
#'   k_range = 2:4,
#'   n_runs = 5L,
#'   nmf_consensus_params = params_nmf_consensus(density_threshold = 2),
#'   .verbose = FALSE
#' )
#' res
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#' @returns `NULL`, invisibly. Called for the side effect of writing the
#' in-memory maps and caches next to the counts. It does not return the
#' object, so do not assign the result.
#'
#' @export
#'
#' @examples
#' # checkpoint the in-memory map and cache next to the counts
#' sc <- demo_single_cells()
#' save_sc_exp_to_disk(sc, type = "rds")
#' list.files(sc@dir_data)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # a fresh handle over a directory written earlier
#' sc <- demo_single_cells(prepped = FALSE)
#' dir <- sc@dir_data
#' sc <- load_existing(SingleCells(dir_data = dir), .verbose = FALSE)
#' dim(sc)
#'
#' unlink(dir, recursive = TRUE, force = TRUE)
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
#' @returns A data.table with columns: gene, group, mean_exp, scaled_exp and
#' pct_exp.
#'
#' @export
#'
#' @examples
#' # mean expression and expressing fraction per cell group
#' sc <- demo_single_cells()
#' dt <- extract_dot_plot_data(
#'   sc,
#'   features = get_gene_names(sc)[1:5],
#'   grouping_variable = "cell_grp"
#' )
#' head(dt, 3)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#' @returns A data.table with a `cell_id` column, one column per gene, and
#' any requested obs columns.
#'
#' @export
#'
#' @examples
#' # normalised expression of three genes with a cell annotation
#' sc <- demo_single_cells()
#' dt <- extract_gene_expression(
#'   sc,
#'   features = get_gene_names(sc)[1:3],
#'   obs_cols = "cell_grp"
#' )
#' head(dt, 3)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # local auto-correlation of every gene on the cached kNN graph
#' sc <- demo_single_cells()
#' res <- hotspot_autocor_sc(sc, .verbose = FALSE)
#' head(res, 3)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#'
#' @examples
#' # local gene-gene correlations over a subset of the genes
#' sc <- demo_single_cells()
#' res <- hotspot_gene_cor_sc(
#'   sc,
#'   genes_to_take = get_gene_names(sc)[1:20],
#'   .verbose = FALSE
#' )
#' res
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#' @param gs_list Named nested list. Every element must itself be a list with
#' at least a `"pos"` element holding the gene identifiers, and optionally a
#' `"neg"` one. A bare character vector is not accepted. The gene identifiers
#' need to be part of the variables of the object.
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
#'
#' @examples
#' # a signed signature alongside a plain one
#' sc <- demo_single_cells()
#' gs_list <- list(
#'   programme_a = list(
#'     pos = get_gene_names(sc)[1:10],
#'     neg = get_gene_names(sc)[11:20]
#'   ),
#'   programme_b = list(pos = get_gene_names(sc)[21:30])
#' )
#' res <- vision_sc(sc, gs_list = gs_list, .verbose = FALSE)
#' dim(res)
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
#' @param gs_list Named nested list. Every element must itself be a list with
#' at least a `"pos"` element holding the gene identifiers, and optionally a
#' `"neg"` one. A bare character vector is not accepted. The gene identifiers
#' need to be part of the variables of the object.
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
#' @returns A list with the following elements:
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
#'
#' @examples
#' # scores plus whether they sit non-randomly on the kNN graph
#' sc <- demo_single_cells()
#' gs_list <- list(
#'   programme_a = list(pos = get_gene_names(sc)[1:10]),
#'   programme_b = list(pos = get_gene_names(sc)[21:30])
#' )
#' res <- vision_w_autocor_sc(
#'   sc,
#'   gs_list = gs_list,
#'   embd_to_use = "pca",
#'   vision_params = params_sc_vision(n_perm = 50L, n_cluster = 3L),
#'   .verbose = FALSE
#' )
#' res$auto_cor_dt
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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

#### dialogue ------------------------------------------------------------------

#' Find multicellular programmes with DIALOGUE
#'
#' @description
#' DIALOGUE looks for programmes of cell-type-specific genes whose activity
#' covaries across the samples that several cell types share. Where
#' co-expression modules ask what varies together *within* a cell type, this
#' asks what varies together *between* them, sample by sample.
#'
#' @details
#' The algorithm works in three stages. Each cell type's features are  collapsed
#' to one row per sample and put through a sparse multi-CCA, giving  every
#' programme a weight vector per cell type plus a provisional gene signature.
#' Then, for every ordered pair of cell types and every candidate gene, a mixed
#' model asks whether a cell's own programme score tracks the partner's
#' expression of that gene in the same sample. Finally the partners are
#' meta-analysed and the scores refit onto the surviving genes by non-negative
#' least squares.
#'
#' `features` is mandatory and there is no default. DIALOGUE does not compute
#' it: it is whatever low-dimensional description of each cell type you trust,
#' and the method is only as good as that choice. Do *not* hand it a slice of a
#' global PCA. Those components mostly carry between-cell-type identity, which
#' is near-constant inside one cell type and gets dropped by the ANOVA filter,
#' leaving the decomposition to work off whatever is left. Run a PCA per cell
#' type instead, which for `SingleCells` means [SingleCellsSubset()] followed by
#' [calculate_pca_sc()] on each subset.
#'
#' The method is unforgiving about study design, and the failure modes are
#' errors rather than bad answers. It needs at least two cell types, at least
#' five samples present in *every* cell type, and enough cells per sample per
#' cell type to clear `abn_c` in [params_dialogue_pmd()]. A handful of samples
#' with thousands of cells each is the regime it was built for; many samples
#' with a dozen cells each is not.
#'
#' @param object `SingleCells`, `SingleCellsSubset` or `MetaCells` class.
#' @param cell_type_col String. Column in the obs table holding the cell type
#' labels.
#' @param sample_col String. Column in the obs table holding the sample labels.
#' The random intercept in stage two is over these. For `MetaCells` the meta
#' cells must have been built *within* samples, otherwise the level is not
#' well-defined.
#' @param features Named list of numeric matrices, one per cell type. Names must
#' match the levels in `cell_type_col`, row names must cover that cell type's
#' cells, and each needs at least two columns. Rows are matched by name, not by
#' position.
#' @param quality_col Optional string. Column in the obs table to use as the
#' cell quality covariate. If `NULL`, defaults to the z-scored log library size.
#' @param gene_ids Optional character. Genes to consider when building
#' signatures. If `NULL`, uses [get_hvg()] on the object.
#' @param pmd_params List, see [params_dialogue_pmd()].
#' @param hlm_params List, see [params_dialogue_hlm()].
#' @param refine_params List, see [params_dialogue_refine()].
#' @param .verbose Boolean or integer. Verbosity.
#'
#' @returns A `DialogueResult` object.
#'
#' @references Jerby-Arnon & Regev, Nature Biotechnology, 2022
#'
#' @export
#'
#' @examples
#' # a planted multicellular programme recovered across three cell types
#' data <- generate_dialogue_test_data()
#' dir <- tempfile("bixverse_dlg")
#' dir.create(dir)
#' object <- load_r_data(
#'   SingleCells(dir_data = dir),
#'   counts = data$counts,
#'   obs = data$obs,
#'   var = data$var,
#'   sc_qc_param = params_sc_min_quality(
#'     min_unique_genes = 10L,
#'     min_lib_size = 50L,
#'     min_cells = 10L
#'   ),
#'   .verbose = FALSE
#' )
#' res <- dialogue_sc(
#'   object,
#'   cell_type_col = "cell_grp",
#'   sample_col = "sample_id",
#'   features = data$features,
#'   gene_ids = data$var$gene_id,
#'   pmd_params = params_dialogue_pmd(k = 2L, n_permutations = 20L),
#'   .verbose = FALSE
#' )
#' res
#'
#' unlink(dir, recursive = TRUE, force = TRUE)
dialogue_sc <- S7::new_generic(
  name = "dialogue_sc",
  dispatch_args = "object",
  fun = function(
    object,
    cell_type_col,
    sample_col,
    features,
    quality_col = NULL,
    gene_ids = NULL,
    pmd_params = params_dialogue_pmd(),
    hlm_params = params_dialogue_hlm(),
    refine_params = params_dialogue_refine(),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)
