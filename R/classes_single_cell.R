# s3 ---------------------------------------------------------------------------

# the S3 classes serve as helpers to store key indices, cell names and
# gene names. additionally, there is a version that stores the in-memory held
# data for the class

## mapping class ---------------------------------------------------------------

#' Helper function to construct relevant maps
#'
#' @description
#' Helper class that contains various mappings and makes it easier to use
#' getters/setters with
#'
#' @return Generates an empty version of the `sc_mapper` class.
#'
#' @export
new_sc_mapper <- function() {
  sc_mapper <- list(
    gene_mapping = NULL,
    cell_mapping = NULL,
    cells_to_keep_idx = NULL,
    hvg_gene_indices = NULL
  )

  class(sc_mapper) <- "sc_mapper"

  return(sc_mapper)
}

## cache class -----------------------------------------------------------------

#' Helper function to hold relevant cached data
#'
#' @description
#' Helper class that contains various data that is held in memory and not on
#' disk.
#'
#' @return Generates an empty version of the `sc_cache` class.
#'
#' @export
new_sc_cache <- function() {
  sc_cache <- list(
    pca_factors = NULL,
    pca_loadings = NULL,
    pca_singular_vals = NULL,
    other_embeddings = list(),
    knn_matrix = NULL,
    snn_graph = NULL
  )

  class(sc_cache) <- "sc_cache"

  return(sc_cache)
}

## setters ---------------------------------------------------------------------

### generics -------------------------------------------------------------------

#### sc_mapper -----------------------------------------------------------------

#' Set gene mapping for sc_mapper object
#'
#' @param x An object to set gene mapping for
#' @param gene_map Named integer indicating indices and names of the genes
#'
#' @export
set_gene_mapping <- function(x, gene_map) {
  UseMethod("set_gene_mapping")
}

#' Set cell mapping for sc_mapper object
#'
#' @param x An object to set cell mapping for
#' @param cell_map Named integer indicating indices and names of the cells
#'
#' @export
set_cell_mapping <- function(x, cell_map) {
  UseMethod("set_cell_mapping")
}

#' Set cells to keep for sc_mapper object
#'
#' @param x An object to set cells to keep for
#' @param cells_to_keep String or integer. The names or indices of the cells
#' to keep in downstream analysis.
#'
#' @export
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
set_hvg <- function(x, hvg) {
  UseMethod("set_hvg")
}

#### sc_cache ------------------------------------------------------------------

#' Set/add PCA factors
#'
#' @param x An object to add the PCA factors for.
#' @param pca_factor Numerical matrix. The matrix with the PCA factors.
#'
#' @export
set_pca_factors <- function(x, pca_factor) {
  UseMethod("set_pca_factors")
}

#' Set/add PCA loadings
#'
#' @param x An object to add the PCA loadings for.
#' @param pca_loading Numerical matrix. The Matrix with the PCA loadings.
#'
#' @export
set_pca_loadings <- function(x, pca_loading) {
  UseMethod("set_pca_loadings")
}

#' Set/add PCA singular values
#'
#' @param x An object to add the singular values for.
#' @param singular_vals Numerical vector. The singular values.
#'
#' @export
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
set_embedding <- function(x, embd, name) {
  UseMethod("set_embedding")
}

#' Set/add KNN
#'
#' @param x An object to add the KNN data to
#' @param knn_mat Numerical matrix. The matrix with the KNN data
#'
#' @export
set_knn <- function(x, knn_mat) {
  UseMethod("set_knn")
}

#' Set/add KNN
#'
#' @param x An object to add the KNN data to
#' @param snn_graph Igraph. The sNN graph for subsequent clustering.
#'
#' @export
set_snn_graph <- function(x, snn_graph) {
  UseMethod("set_snn_graph")
}

#' Remove the KNN data
#'
#' @param x An object from which to remove the kNN data
#'
#' @export
remove_knn <- function(x) {
  UseMethod("remove_knn")
}

#' Remove the sNN graph
#'
#' @param x An object from which to remove the sNN graph
#'
#' @export
remove_snn_graph <- function(x) {
  UseMethod("remove_knn")
}

### methods --------------------------------------------------------------------

#### sc_mapper -----------------------------------------------------------------

#' @rdname set_gene_mapping
#'
#' @export
set_gene_mapping.sc_mapper <- function(x, gene_map) {
  # checks
  checkmate::assertClass(x, "sc_mapper")
  checkmate::qassert(gene_map, "N+")
  checkmate::assertNamed(gene_map)

  x[["gene_mapping"]] <- gene_map

  return(x)
}

#' @rdname set_cell_mapping
#'
#' @export
set_cell_mapping.sc_mapper <- function(x, cell_map) {
  # checks
  checkmate::assertClass(x, "sc_mapper")
  checkmate::qassert(cell_map, "N+")
  checkmate::assertNamed(cell_map)

  x[["cell_mapping"]] <- cell_map

  return(x)
}

#' @rdname set_cells_to_keep
#'
#' @export
set_cells_to_keep.sc_mapper <- function(x, cells_to_keep) {
  # checks
  checkmate::assertClass(x, "sc_mapper")
  checkmate::qassert(cells_to_keep, c("N+", "S+"))

  # transform to Rust 0-based indexing
  res <- if (is.numeric(cells_to_keep)) {
    cells_to_keep - 1
  } else {
    x[["cell_mapping"]][cells_to_keep] - 1
  }

  x[["cells_to_keep_idx"]] <- res

  return(x)
}

#' @rdname set_hvg
#'
#' @export
set_hvg.sc_mapper <- function(x, hvg) {
  # checks
  checkmate::assertClass(x, "sc_mapper")
  checkmate::qassert(hvg, c("N+", "S+"))

  # transform to Rust 0-based indexing
  res <- if (is.numeric(hvg)) {
    hvg - 1
  } else {
    x[["gene_mapping"]][hvg] - 1
  }

  x[["hvg_gene_indices"]] <- res

  return(x)
}

#### sc_cache ------------------------------------------------------------------

#' @rdname set_pca_factors
#'
#' @export
set_pca_factors.sc_cache <- function(x, pca_factor) {
  # checks
  checkmate::assertClass(x, "sc_cache")
  checkmate::assertMatrix(pca_factor, mode = "numeric")

  x[["pca_factors"]] <- pca_factor

  return(x)
}

#' @rdname set_pca_loadings
#'
#' @export
set_pca_loadings.sc_cache <- function(x, pca_loading) {
  # checks
  checkmate::assertClass(x, "sc_cache")
  checkmate::assertMatrix(pca_loading, mode = "numeric")

  x[["pca_loadings"]] <- pca_loading

  return(x)
}

#' @rdname set_pca_singular_vals
#'
#' @export
set_pca_singular_vals.sc_cache <- function(x, singular_vals) {
  # checks
  checkmate::assertClass(x, "sc_cache")
  checkmate::qassert(singular_vals, "N+")

  x[["pca_singular_vals"]] <- singular_vals

  return(x)
}

#' @rdname set_embedding
#'
#' @export
set_embedding.sc_cache <- function(x, embd, name) {
  # checks
  checkmate::assertClass(x, "sc_cache")
  checkmate::assertMatrix(embd, mode = "numeric")
  checkmate::qassert(name, "S1")

  x[["other_embeddings"]][[name]] <- embd

  return(x)
}

#' @rdname set_knn
#'
#' @export
set_knn.sc_cache <- function(x, knn_mat) {
  # checks
  checkmate::assertClass(x, "sc_cache")
  checkmate::assertMatrix(knn_mat, mode = "numeric")

  x[["knn_matrix"]] <- knn_mat

  return(x)
}

#' @rdname set_snn_graph
#'
#' @export
set_snn_graph.sc_cache <- function(x, snn_graph) {
  # checks
  checkmate::assertClass(x, "sc_cache")
  checkmate::assertClass(snn_graph, "igraph")

  x[["snn_graph"]] <- snn_graph

  return(x)
}

#' @rdname remove_knn
#'
#' @export
remove_knn.sc_cache <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_cache")

  x[["knn_matrix"]] <- NULL

  return(x)
}

#' @rdname remove_snn_graph
#'
#' @export
remove_snn_graph.sc_cache <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_cache")

  x[["snn_graph"]] <- NULL

  return(x)
}

## getters ---------------------------------------------------------------------

### generics -------------------------------------------------------------------

#### sc_mapper -----------------------------------------------------------------

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

#' Get the HVG
#'
#' @param x An object to get HVG from.
#'
#' @export
get_hvg <- function(x) {
  UseMethod("get_hvg")
}

#### sc_cache ------------------------------------------------------------------

#' Get the PCA factors
#'
#' @param x An object to get PCA factors from.
#'
#' @export
get_pca_factors <- function(x) {
  UseMethod("get_pca_factors")
}

#' Get the PCA loadings
#'
#' @param x An object to get PCA loadings from.
#'
#' @export
get_pca_loadings <- function(x) {
  UseMethod("get_pca_loadings")
}

#' Get the PCA singular values
#'
#' @param x An object to get PCA singular values from.
#'
#' @export
get_pca_singular_val <- function(x) {
  UseMethod("get_pca_singular_val")
}

#' Get the embedding from the cache
#'
#' @description
#' General wrapper function that can be used to pull out any embedding stored
#' in the `sc_cache`.
#'
#' @param x An object to get embedding from
#' @param embd_name String. The name of the embedding to return. The function
#' will throw an error if the embedding does not exist.
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
#' @export
get_available_embeddings <- function(x) {
  UseMethod("get_available_embeddings")
}

#' Get the sNN graph
#'
#' @param x An object to get the sNN graph from.
#'
#' @export
get_snn_graph <- function(x) {
  UseMethod("get_snn_graph")
}

### methods --------------------------------------------------------------------

#### sc_mapper -----------------------------------------------------------------

#' @rdname get_gene_names
#'
#' @export
get_gene_names.sc_mapper <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_mapper")

  return(names(x[["gene_mapping"]]))
}

#' @rdname get_gene_indices
#'
#' @export
get_gene_indices.sc_mapper <- function(x, gene_ids, rust_index) {
  # checks
  checkmate::assertClass(x, "sc_mapper")
  checkmate::qassert(gene_ids, "S+")
  checkmate::qassert(rust_index, "B1")

  gene_map <- x$gene_mapping
  indices <- which(names(gene_map) %in% gene_ids)
  if (rust_index) {
    indices <- indices - 1
  }

  return(as.integer(indices))
}

#' @rdname get_cell_indices
#'
#' @export
get_cell_indices.sc_mapper <- function(x, cell_ids, rust_index) {
  # checks
  checkmate::assertClass(x, "sc_mapper")
  checkmate::qassert(cell_ids, "S+")
  checkmate::qassert(rust_index, "B1")

  cell_map <- x$cell_mapping
  indices <- which(names(cell_map) %in% cell_ids)
  if (rust_index) {
    indices <- indices - 1
  }

  return(as.integer(indices))
}

#' @rdname get_cells_to_keep
#'
#' @export
get_cells_to_keep.sc_mapper <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_mapper")

  return(as.integer(x[["cells_to_keep_idx"]]))
}

#' @rdname get_cell_names
#'
#' @export
get_cell_names.sc_mapper <- function(x, filtered = FALSE) {
  # checks
  checkmate::assertClass(x, "sc_mapper")
  checkmate::qassert(filtered, "B1")

  cell_names <- names(x[["cell_mapping"]])

  if (filtered) {
    cells_to_keep <- get_cells_to_keep(x)
    if (length(cells_to_keep) > 0) {
      cell_names <- cell_names[cells_to_keep + 1]
    }
  }

  return(cell_names)
}

#' @rdname get_hvg
#'
#' @export
get_hvg.sc_mapper <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_mapper")

  hvg_indices <- as.integer(x[["hvg_gene_indices"]])
  if (length(hvg_indices) == 0) {
    warning("No highly variable features found in the class.")
  }

  return(hvg_indices)
}

#### sc_cache ------------------------------------------------------------------

#' @rdname get_pca_factors
#'
#' @export
get_pca_factors.sc_cache <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_cache")

  res <- x[["pca_factors"]]

  return(res)
}

#' @rdname get_pca_loadings
#'
#' @export
get_pca_loadings.sc_cache <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_cache")

  res <- x[["pca_loadings"]]

  return(res)
}

#' @rdname get_pca_singular_val
#'
#' @export
get_pca_singular_val.sc_cache <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_cache")

  res <- x[["pca_singular_vals"]]

  return(res)
}

#' @rdname get_embedding
#'
#' @export
get_embedding.sc_cache <- function(x, embd_name) {
  possible_embedding_names <- c(names(x[["other_embeddings"]]), "pca")

  checkmate::assertClass(x, "sc_cache")
  checkmate::assertChoice(embd_name, possible_embedding_names)

  embd <- if (embd_name == "pca") {
    get_pca_factors(x)
  } else {
    x[["other_embeddings"]][[embd_name]]
  }

  if (is.null(embd)) {
    stop(paste(
      "The sought embedding was not found.",
      "Please verify which functions you have run."
    ))
  }

  return(embd)
}

#' @rdname get_available_embeddings
#'
#' @export
get_available_embeddings.sc_cache <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_cache")

  pca_available <- ifelse(!is.null(x[["pca_factors"]]), "pca", NULL)
  other_embeddings <- names(x[["other_embeddings"]])

  c(pca_available, other_embeddings)
}

#' @rdname get_knn_mat
#'
#' @export
get_knn_mat.sc_cache <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_cache")

  return(x[["knn_matrix"]])
}

#' @rdname get_snn_graph
#'
#' @export
get_snn_graph.sc_cache <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_cache")

  return(x[["snn_graph"]])
}

# s7 ---------------------------------------------------------------------------

# main S7 class for single cell experiments.

## single cell class -----------------------------------------------------------

#' @title bixverse single cell class
#'
#' @description
#' This is the `bixverse`-based single cell class. Under the hood it uses a
#' DuckDB for obs and vars storing, and a Rust-based binarised file format to
#' store the raw and normalised counts. In both cases, the idea is not to hold
#' any data that is not needed at a given point of time in memory, but leverage
#' speedy on-disk computations and streaming engines powered by Rust and DuckDB
#' to run the analysis.
#'
#' @param dir_data String. This is the directory in which the experimental files
#' will be
#'
#' @section Properties:
#' \describe{
#'   \item{db_connection}{This contains an R6 class with DuckDB pointers and
#'   wrappers to interact with the table-like data for this experiment.}
#'   \item{count_connection}{This contains an R6-like environment that points
#'   to Rust functions that can work on the counts more specifically.}
#'   \item{dir_data}{Path to the directory in which the data will be saved on
#'   disk.}
#'   \item{sc_cache}{Class with cached data. Contains less memory-heavy objects
#'   such as embeddings, kNN information or sNN graphs.}
#'   \item{sc_map}{Class containing various mapping information such as HVG
#'   indices, cells to keep, etc.}
#'   \item{dims}{Dimensions of the original data.}
#' }
#'
#' @return Returns the `single_cell_exp` class for further operations.
#'
#' @export
single_cell_exp <- S7::new_class(
  name = "single_cell_exp",
  properties = list(
    db_connection = S7::class_any,
    count_connection = S7::class_any,
    dir_data = S7::class_character,
    sc_cache = S7::class_any,
    sc_map = S7::class_any,
    dims = S7::class_integer
  ),
  constructor = function(dir_data) {
    # nightly_feature()
    # checks
    checkmate::assertDirectoryExists(dir_data)

    # generate the Rust pointer
    count_connection <- SingeCellCountData$new(
      f_path_cells = file.path(dir_data, "counts_cells.bin"),
      f_path_genes = file.path(dir_data, "counts_genes.bin")
    )

    # generate the DuckDB connector
    db_connection <- single_cell_duckdb_con$new(
      db_dir = dir_data
    )

    S7::new_object(
      S7::S7_object(),
      db_connection = db_connection,
      count_connection = count_connection,
      dir_data = dir_data,
      sc_cache = new_sc_cache(),
      sc_map = new_sc_mapper(),
      dims = c(0L, 0L)
    )
  }
)

### getters --------------------------------------------------------------------

#### env getters ---------------------------------------------------------------

#' Getter for the single cell DuckDB connection
#'
#' @param object `single_cell_exp` class.
#'
#' @return The DuckDB connector
#'
#' @export
get_sc_duckdb <- S7::new_generic(
  name = "get_sc_duckdb",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_sc_duckdb single_cell_exp
#'
#' @export
S7::method(get_sc_duckdb, single_cell_exp) <- function(object) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")

  return(S7::prop(object, "db_connection"))
}

#' Getter for the single cell Rust pointer
#'
#' @param object `single_cell_exp` class.
#'
#' @return The Rust structure
#'
#' @export
get_sc_rust_ptr <- S7::new_generic(
  name = "get_sc_rust_ptr",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_sc_rust_ptr single_cell_exp
#'
#' @export
S7::method(get_sc_rust_ptr, single_cell_exp) <- function(object) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")

  return(S7::prop(object, "count_connection"))
}

#' Get the file path to the counts in gene format
#'
#' @param object `single_cell_exp` class.
#'
#' @return The path to the `counts_genes.bin`
get_rust_count_gene_f_path <- S7::new_generic(
  name = "get_rust_count_gene_f_path",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_rust_count_gene_f_path single_cell_exp
#'
#' @export
S7::method(get_rust_count_gene_f_path, single_cell_exp) <- function(object) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")

  f_path <- file.path(S7::prop(object, "dir_data"), "counts_genes.bin")

  return(f_path)
}

#' Get the file path to the counts in cell format
#'
#' @param object `single_cell_exp` class.
#'
#' @return The path to the `counts_cells.bin`
get_rust_count_cell_f_path <- S7::new_generic(
  name = "get_rust_count_cell_f_path",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_rust_count_cell_f_path single_cell_exp
#'
#' @export
S7::method(get_rust_count_cell_f_path, single_cell_exp) <- function(object) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")

  f_path <- file.path(S7::prop(object, "dir_data"), "counts_cells.bin")

  return(f_path)
}

#### duckdb getters ------------------------------------------------------------

#' @method get_sc_obs single_cell_exp
#'
#' @export
S7::method(get_sc_obs, single_cell_exp) <- function(
  object,
  indices = NULL,
  cols = NULL,
  filtered = FALSE
) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::qassert(indices, c("0", "I+"))
  checkmate::qassert(filtered, "B1")
  duckdb_con <- get_sc_duckdb(object)

  obs_table <- duckdb_con$get_obs_table(
    indices = indices,
    cols = cols,
    filtered = filtered
  )

  obs_table

  return(obs_table)
}

#' @method get_sc_var single_cell_exp
#'
#' @export
S7::method(get_sc_var, single_cell_exp) <- function(
  object,
  indices = NULL,
  cols = NULL
) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::qassert(indices, c("0", "I+"))

  duckdb_con <- get_sc_duckdb(object)

  var_table <- duckdb_con$get_vars_table(indices = indices, cols = cols)

  return(var_table)
}

#' @method `[[` single_cell_exp
#'
#' @export
S7::method(`[[`, single_cell_exp) <- function(x, i, ...) {
  if (missing(i)) {
    i <- NULL
  }

  if (checkmate::qtest(i, "S+")) {
    get_sc_obs(x, cols = i, filtered = TRUE)
  } else if (checkmate::qtest(i, "I+")) {
    get_sc_obs(x, indices = i, filtered = TRUE)
  } else if (checkmate::qtest(i, "0")) {
    get_sc_obs(x, filtered = TRUE)
  } else {
    stop("Invalid type")
  }
}

### map getter -----------------------------------------------------------------

#' Getter for the different maps in the object
#'
#' @param object `single_cell_exp` class.
#'
#' @returns Returns the maps within the class.
#'
#' @export
get_sc_map <- S7::new_generic(
  name = "get_sc_maps",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_sc_var single_cell_exp
#'
#' @export
S7::method(get_sc_map, single_cell_exp) <- function(
  object
) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")

  res <- S7::prop(object, "sc_map")

  return(res)
}

### sc cache getter ------------------------------------------------------------

#' Getter the memory-stored data from the class
#'
#' @param object `single_cell_exp` class.
#'
#' @returns Returns the sc_cache class from the object with all the
#' memory-stored data.
#'
#' @export
get_sc_cache <- S7::new_generic(
  name = "get_sc_cache",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_sc_cache single_cell_exp
#'
#' @export
S7::method(get_sc_cache, single_cell_exp) <- function(
  object
) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")

  res <- S7::prop(object, "sc_cache")

  return(res)
}

#### count getters -------------------------------------------------------------

#' @method get_sc_counts single_cell_exp
#'
#' @export
S7::method(get_sc_counts, single_cell_exp) <- function(
  object,
  assay = c("raw", "norm"),
  return_format = c("cell", "gene"),
  cell_indices = NULL,
  gene_indices = NULL,
  use_cells_to_keep = TRUE,
  .verbose = TRUE
) {
  assay <- match.arg(assay)
  return_format <- match.arg(return_format)

  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::assertChoice(assay, c("raw", "norm"))
  checkmate::assertChoice(return_format, c("cell", "gene"))
  checkmate::qassert(cell_indices, c("0", "I+"))
  checkmate::qassert(gene_indices, c("0", "I+"))
  checkmate::qassert(.verbose, "B1")

  requireNamespace("Matrix", quietly = TRUE)

  rust_con <- get_sc_rust_ptr(object)

  sc_map <- get_sc_map(object)

  cells_to_keep <- get_cells_to_keep(object)

  if (use_cells_to_keep) {
    if (!is.null(cell_indices) & length(cells_to_keep) > 0) {
      # deal with the case that the cells to keep where specified
      cell_indices <- as.integer(intersect(cell_indices, cells_to_keep + 1))
    } else if (length(cells_to_keep) > 0) {
      cell_indices <- as.integer(cells_to_keep + 1)
    }
  }

  # Get raw data from Rust
  count_data <- get_counts_from_rust(
    rust_con = rust_con,
    assay = assay,
    return_format = return_format,
    cell_indices = cell_indices,
    gene_indices = gene_indices,
    .verbose = .verbose
  )

  # Create sparse matrix
  count_data <- create_sparse_matrix(
    count_data = count_data,
    return_format = return_format
  )

  # Set names and subset if needed
  count_data <- finalise_matrix(
    matrix = count_data,
    return_format = return_format,
    cell_indices = cell_indices,
    gene_indices = gene_indices,
    sc_map = sc_map
  )

  return(count_data)
}

#' @method `[` single_cell_exp
#'
#' @export
S7::method(`[`, single_cell_exp) <- function(
  x,
  i,
  j,
  ...,
  assay = c("raw", "norm"),
  return_format = c("cell", "gene"),
  use_cells_to_keep = TRUE,
  drop = TRUE
) {
  if (missing(i)) {
    i <- NULL
  }
  if (missing(j)) {
    j <- NULL
  }

  assay <- match.arg(assay)
  return_format <- match.arg(return_format)

  # asserts
  checkmate::qassert(i, c("I+", "0"))
  checkmate::qassert(j, c("I+", "0"))
  checkmate::assertChoice(assay, c("raw", "norm"))
  checkmate::assertChoice(return_format, c("cell", "gene"))

  get_sc_counts(
    object = x,
    assay = assay,
    return_format = return_format,
    cell_indices = i,
    gene_indices = j,
    use_cells_to_keep = use_cells_to_keep,
    .verbose = FALSE
  )
}

##### helpers ------------------------------------------------------------------

#' Helper function to get counts from Rust
#'
#' @param rust_con `SingeCellCountData` class. The connector to Rust.
#' @param assay String. One of `c("raw", "norm")`.
#' @param return_format String. One of `c("cell", "gene")`.
#' @param cell_indices Optional integer vector. The index positions of cells to
#' return.
#' @param gene_indices Optional integer vector. The index positions of genes to
#' return.
#' @param .verbose Boolean. Controls verbosity
#'
#' @return Returns a list with:
#' \itemize{
#'  \item indptr - The index pointers of the compressed sparse format.
#'  \item indices - The indices of the data.
#'  \item data - The underlying data.
#'  \item no_cells - Number of cells.
#'  \item no_genes - Number of genes.
#' }
get_counts_from_rust <- function(
  rust_con,
  assay,
  return_format,
  cell_indices,
  gene_indices,
  .verbose = TRUE
) {
  # checks
  checkmate::assertClass(rust_con, "SingeCellCountData")
  checkmate::assertChoice(assay, c("raw", "norm"))
  checkmate::assertChoice(return_format, c("cell", "gene"))
  checkmate::qassert(cell_indices, c("0", "I+"))
  checkmate::qassert(gene_indices, c("0", "I+"))
  checkmate::qassert(.verbose, "B1")

  res <- if (return_format == "cell") {
    if (is.null(cell_indices)) {
      rust_con$return_full_mat(
        assay = assay,
        cell_based = TRUE,
        verbose = .verbose
      )
    } else {
      rust_con$get_cells_by_indices(indices = cell_indices, assay = assay)
    }
  } else {
    # gene
    if (is.null(gene_indices)) {
      rust_con$return_full_mat(
        assay = assay,
        cell_based = FALSE,
        verbose = .verbose
      )
    } else {
      rust_con$get_genes_by_indices(indices = gene_indices, assay = assay)
    }
  }

  return(res)
}

#' Helper function transform Rust counts into sparse matrices
#'
#' @param count_data A list. Output of [bixverse::get_counts_from_rust()].
#' @param return_format String. One of `c("cell", "gene")`.
#'
#' @return The sparse matrix in CSR or CSC format, pending the choice in
#' `return_format`
create_sparse_matrix <- function(count_data, return_format) {
  # checks
  checkmate::assertList(count_data)
  checkmate::assertNames(
    names(count_data),
    must.include = c("indices", "indptr", "data", "no_cells", "no_genes")
  )
  checkmate::assertChoice(return_format, c("cell", "gene"))

  matrix_class <- if (return_format == "cell") "dgRMatrix" else "dgCMatrix"
  index_slot <- if (return_format == "cell") "j" else "i"

  sparse_mat <- with(count_data, {
    args <- list(
      p = as.integer(indptr),
      x = as.numeric(data),
      Dim = as.integer(c(no_cells, no_genes))
    )
    args[[index_slot]] <- as.integer(indices)

    do.call(new, c(matrix_class, args))
  })

  return(sparse_mat)
}

#' Finalise the count matrix
#'
#' @param matrix The sparse matrix to finalise
#' @param return_format String. One of `c("cell", "gene")`.
#' @param cell_indices Optional integer. The cell indices to return.
#' @param gene_indices Optional integer. The gene indices to return.
#' @param sc_map A `sc_mapper` class. Contains various mapping information.
#'
#' @return The finalised matrix.
finalise_matrix <- function(
  matrix,
  return_format,
  cell_indices,
  gene_indices,
  sc_map
) {
  checkmate::assert(
    checkmate::checkClass(matrix, "dgRMatrix"),
    checkmate::checkClass(matrix, "dgCMatrix")
  )
  checkmate::assertChoice(return_format, c("cell", "gene"))
  checkmate::qassert(cell_indices, c("0", "I+"))
  checkmate::qassert(gene_indices, c("0", "I+"))
  checkmate::assertClass(sc_map, "sc_mapper")

  gene_names <- get_gene_names(sc_map)
  cell_names <- get_cell_names(sc_map)

  if (return_format == "cell") {
    rownames(matrix) <- if (is.null(cell_indices)) {
      cell_names
    } else {
      cell_names[cell_indices]
    }

    if (is.null(gene_indices)) {
      colnames(matrix) <- gene_names
    } else {
      matrix <- matrix[, gene_indices]
      colnames(matrix) <- gene_names[gene_indices]
    }
  } else {
    colnames(matrix) <- if (is.null(gene_indices)) {
      gene_names
    } else {
      gene_names[gene_indices]
    }

    # Set row names (cells) and subset if needed
    if (is.null(cell_indices)) {
      rownames(matrix) <- cell_names
    } else {
      matrix <- matrix[cell_indices, ]
      rownames(matrix) <- cell_names[cell_indices]
    }
  }

  matrix
}

#### sc map --------------------------------------------------------------------

#' @name get_cell_names.single_cell_exp
#'
#' @title Get the cell names from a `single_cell_exp`.
#'
#' @rdname get_cell_names
#'
#' @method get_cell_names single_cell_exp
S7::method(get_cell_names, single_cell_exp) <- function(
  x,
  filtered = FALSE
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")

  # add the data using the S3 method
  cell_names <- get_cell_names(
    x = S7::prop(x, "sc_map"),
    filtered = filtered
  )

  return(cell_names)
}

#' @name get_gene_names.single_cell_exp
#'
#' @title Get the gene names from a `single_cell_exp`.
#'
#' @rdname get_gene_names
#'
#' @method get_gene_names single_cell_exp
S7::method(get_gene_names, single_cell_exp) <- function(
  x
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")

  # add the data using the S3 method
  gene_names <- get_gene_names(
    x = S7::prop(x, "sc_map")
  )

  return(gene_names)
}

#' @name get_cells_to_keep.single_cell_exp
#'
#' @title Get the cells to keep from a `single_cell_exp`.
#'
#' @rdname get_cells_to_keep
#'
#' @method get_cells_to_keep single_cell_exp
S7::method(get_cells_to_keep, single_cell_exp) <- function(
  x
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")

  # forward to S3
  res <- get_cells_to_keep(
    x = S7::prop(x, "sc_map")
  )
  # special case that this has not been set. Return all cell indices then
  if (length(res) == 0) {
    no_cells <- S7::prop(x, "dims")[1]
    res <- seq_len(no_cells) - 1 # 0 index for Rust
  }

  return(as.integer(res))
}

#' @name get_gene_indices.single_cell_exp
#'
#' @title Get the gene indices from a `single_cell_exp`.
#'
#' @rdname get_gene_indices
#'
#' @method get_gene_indices single_cell_exp
S7::method(get_gene_indices, single_cell_exp) <- function(
  x,
  gene_ids,
  rust_index
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")
  checkmate::qassert(gene_ids, "S+")
  checkmate::qassert(rust_index, "B1")

  # forward to S3
  res <- get_gene_indices(
    x = S7::prop(x, "sc_map"),
    gene_ids = gene_ids,
    rust_index = rust_index
  )

  return(res)
}


#' @name get_cell_indices.single_cell_exp
#'
#' @title Set the gene mapping for a `single_cell_exp` class.
#'
#' @rdname get_cell_indices
#'
#' @method get_cell_indices single_cell_exp
S7::method(get_cell_indices, single_cell_exp) <- function(
  x,
  cell_ids,
  rust_index
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")
  checkmate::qassert(cell_ids, "S+")
  checkmate::qassert(rust_index, "B1")

  # add the data using the S3 method
  res <- get_cell_indices(
    x = S7::prop(x, "sc_map"),
    cell_ids = cell_ids,
    rust_index = rust_index
  )

  return(res)
}

#' @name get_hvg.single_cell_exp
#'
#' @title Get the highly variable gene indices from a `single_cell_exp`.
#'
#' @rdname get_hvg
#'
#' @method get_hvg single_cell_exp
S7::method(get_hvg, single_cell_exp) <- function(
  x
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")

  # forward to S3
  res <- get_hvg(
    x = S7::prop(x, "sc_map")
  )

  return(res)
}

#### sc_cache ------------------------------------------------------------------

#' @name get_pca_factors.single_cell_exp
#'
#' @title Get the PCA factors from a `single_cell_exp`.
#'
#' @rdname get_pca_factors
#'
#' @method get_pca_factors single_cell_exp
S7::method(get_pca_factors, single_cell_exp) <- function(
  x
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")

  # forward to S3
  res <- get_pca_factors(
    x = S7::prop(x, "sc_cache")
  )

  if (is.null(res)) {
    return(NULL)
  }

  rownames(res) <- get_cell_names(x, filtered = TRUE)
  colnames(res) <- sprintf("PC_%i", 1:ncol(res))

  return(res)
}

#' @name get_pca_loadings.single_cell_exp
#'
#' @title Get the PCA loadings from a `single_cell_exp`.
#'
#' @rdname get_pca_loadings
#'
#' @method get_pca_loadings single_cell_exp
S7::method(get_pca_loadings, single_cell_exp) <- function(
  x
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")

  # forward to S3
  res <- get_pca_loadings(
    x = S7::prop(x, "sc_cache")
  )

  colnames(res) <- sprintf("PC_%i", 1:ncol(res))
  rownames(res) <- get_gene_names(x)[get_hvg(x) + 1]

  return(res)
}

#' @name get_pca_singular_val.single_cell_exp
#'
#' @title Get the PCA singular values from a `single_cell_exp`.
#'
#' @rdname get_pca_singular_val
#'
#' @method get_pca_singular_val single_cell_exp
S7::method(get_pca_singular_val, single_cell_exp) <- function(
  x
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")

  # forward to S3
  res <- get_pca_singular_val(
    x = S7::prop(x, "sc_cache")
  )

  return(res)
}

#' @name get_embedding.single_cell_exp
#'
#' @title Get the embeddings from a `single_cell_exp`.
#'
#' @rdname get_embedding
#'
#' @method get_embedding single_cell_exp
S7::method(get_embedding, single_cell_exp) <- function(
  x,
  embd_name
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")
  checkmate::qassert(embd_name, "S1")

  # forward to S3
  res <- get_embedding(
    x = S7::prop(x, "sc_cache"),
    embd_name = embd_name
  )

  rownames(res) <- get_cell_names(x, filtered = TRUE)

  return(res)
}

#' @name get_available_embeddings.single_cell_exp
#'
#' @title Get the embeddings from a `single_cell_exp`.
#'
#' @rdname get_available_embeddings
#'
#' @method get_available_embeddings single_cell_exp
S7::method(get_available_embeddings, single_cell_exp) <- function(
  x
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")

  # forward to S3
  res <- get_available_embeddings(
    x = S7::prop(x, "sc_cache")
  )

  return(res)
}

#' @name get_knn_mat.single_cell_exp
#'
#' @title Get the KNN matrix from a `single_cell_exp`.
#'
#' @rdname get_knn_mat
#'
#' @method get_knn_mat single_cell_exp
S7::method(get_knn_mat, single_cell_exp) <- function(
  x
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")

  # forward to S3
  res <- get_knn_mat(
    x = S7::prop(x, "sc_cache")
  )

  return(res)
}

#' @name get_snn_graph.single_cell_exp
#'
#' @title Get the SNN graph from a `single_cell_exp`.
#'
#' @rdname get_snn_graph
#'
#' @method get_snn_graph single_cell_exp
S7::method(get_snn_graph, single_cell_exp) <- function(
  x
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")

  # forward to S3
  res <- get_snn_graph(
    x = S7::prop(x, "sc_cache")
  )

  return(res)
}

### setters --------------------------------------------------------------------

#### obs -----------------------------------------------------------------------

#' Add a new column to the obs table
#'
#' @param object `bixverse::single_cell_exp` class.
#' @param col_name String. The name of the column to add.
#' @param new_data Atomic vector. The data to add to the column. Needs to be
#' of same length as [bixverse::get_cells_to_keep()] and have the same order.
#'
#' @return The class with updated obs table in the DuckDB
#'
#' @export
set_sc_new_obs_col <- S7::new_generic(
  name = "set_sc_new_obs_col",
  dispatch_args = "object",
  fun = function(
    object,
    col_name,
    new_data
  ) {
    S7::S7_dispatch()
  }
)

#' @method set_sc_new_obs_col single_cell_exp
#'
#' @export
S7::method(set_sc_new_obs_col, single_cell_exp) <- function(
  object,
  col_name,
  new_data
) {
  cells_to_keep <- (get_cells_to_keep(object) + 1)

  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::qassert(col_name, "S1")
  checkmate::qassert(new_data, sprintf("a%i", length(cells_to_keep)))

  duckdb_con <- get_sc_duckdb(object)

  new_data <- data.table::data.table(new_data)
  data.table::setnames(new_data, "new_data", col_name)
  new_data[, cell_idx := sort(cells_to_keep)]

  duckdb_con$join_data_obs(new_data)

  return(object)
}

#' Add multiple new columns to the obs table
#'
#' @param object `bixverse::single_cell_exp` class.
#' @param new_data Named list. The names will be the column names and the
#' elements will be added to the obs table. Needs to be of same length as
#' [bixverse::get_cells_to_keep()] and have the same order!
#'
#' @return The class with updated obs table in the DuckDB
#'
#' @export
set_sc_new_obs_col_multiple <- S7::new_generic(
  name = "set_sc_new_obs_col_multiple",
  dispatch_args = "object",
  fun = function(
    object,
    new_data
  ) {
    S7::S7_dispatch()
  }
)

#' @method set_sc_new_obs_col_multiple single_cell_exp
#'
#' @export
S7::method(set_sc_new_obs_col_multiple, single_cell_exp) <- function(
  object,
  new_data
) {
  cells_to_keep <- (get_cells_to_keep(object) + 1)

  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::assertList(new_data, names = "named", types = "atomic")
  checkmate::assertTRUE(all(
    purrr::map_dbl(new_data, length) == length(cells_to_keep)
  ))

  duckdb_con <- get_sc_duckdb(object)

  new_data <- data.table::as.data.table(new_data)
  new_data[, cell_idx := sort(cells_to_keep)]

  duckdb_con$join_data_obs(new_data)

  return(object)
}

#' @method `[[<-` single_cell_exp
#'
#' @export
S7::method(`[[<-`, single_cell_exp) <- function(x, i, ..., value) {
  checkmate::assertClass(x, "bixverse::single_cell_exp")
  checkmate::qassert(i, "S+")

  if (length(i) == 1) {
    checkmate::qassert(value, "a")
  } else {
    checkmate::assertList(value, names = "named", types = "atomic")
  }

  if (length(i) == 1) {
    x <- set_sc_new_obs_col(object = x, col_name = i, new_data = value)
  } else {
    names(value) <- i
    x <- set_sc_new_obs_col_multiple(object = x, new_data = value)
  }

  return(x)
}


#' Add a new column to the var table
#'
#' @param object `bixverse::single_cell_exp` class.
#' @param data_list Named list with the data to add.
#'
#' @return The class with updated var table in the DuckDB
#'
#' @export
set_sc_new_var_cols <- S7::new_generic(
  name = "set_sc_new_var_cols",
  dispatch_args = "object",
  fun = function(
    object,
    data_list
  ) {
    S7::S7_dispatch()
  }
)

#' @method set_sc_new_var_cols single_cell_exp
#'
#' @export
S7::method(set_sc_new_var_cols, single_cell_exp) <- function(
  object,
  data_list
) {
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::assertList(data_list, types = "atomic", names = "named")

  new_data <- data.table::as.data.table(data_list)

  duckdb_con <- get_sc_duckdb(object)
  duckdb_con$add_data_var(new_data = new_data)

  return(object)
}

#### sc map --------------------------------------------------------------------

#' @name set_gene_mapping.single_cell_exp
#'
#' @title Set the gene mapping for a `single_cell_exp` class.
#'
#' @rdname set_gene_mapping
#'
#' @method set_gene_mapping single_cell_exp
S7::method(set_gene_mapping, single_cell_exp) <- function(
  x,
  gene_map
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")
  checkmate::qassert(gene_map, "N+")
  checkmate::assertNamed(gene_map, "named")

  # add the data using the S3 method
  S7::prop(x, "sc_map") <- set_gene_mapping(
    x = S7::prop(x, "sc_map"),
    gene_map = gene_map
  )

  return(x)
}

#' @name set_cell_mapping.single_cell_exp
#'
#' @title Set the cell mapping for a `single_cell_exp` class.
#'
#' @rdname set_cell_mapping
#'
#' @method set_cell_mapping single_cell_exp
S7::method(set_cell_mapping, single_cell_exp) <- function(
  x,
  cell_map
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")
  checkmate::qassert(cell_map, "N+")
  checkmate::assertNamed(cell_map, "named")

  # add the data using the S3 method
  S7::prop(x, "sc_map") <- set_cell_mapping(
    x = S7::prop(x, "sc_map"),
    cell_map = cell_map
  )

  return(x)
}

#' @name set_cells_to_keep.single_cell_exp
#'
#' @title Set the cell mapping for a `single_cell_exp` class.
#'
#' @rdname set_cells_to_keep
#'
#' @method set_cells_to_keep single_cell_exp
S7::method(set_cells_to_keep, single_cell_exp) <- function(
  x,
  cells_to_keep
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")
  checkmate::qassert(cells_to_keep, c("I+", "S+"))

  # add the data using the S3 method
  S7::prop(x, "sc_map") <- set_cells_to_keep(
    x = S7::prop(x, "sc_map"),
    cells_to_keep = cells_to_keep
  )

  # remove the cells from the obs table
  duckdb_con <- get_sc_duckdb(x)

  duckdb_con$set_cells_to_keep(
    cell_idx_to_keep = as.integer(get_cells_to_keep(x) + 1)
  )

  return(x)
}

#' @name set_hvg.single_cell_exp
#'
#' @title Set the highly variable genes for a `single_cell_exp` class.
#'
#' @rdname set_hvg
#'
#' @method set_hvg single_cell_exp
S7::method(set_hvg, single_cell_exp) <- function(
  x,
  hvg
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")
  checkmate::qassert(hvg, c("I+", "S+"))

  # add the data using the S3 method
  S7::prop(x, "sc_map") <- set_hvg(
    x = S7::prop(x, "sc_map"),
    hvg = hvg
  )

  return(x)
}

#### sc cache ------------------------------------------------------------------

#' @name set_pca_factors.single_cell_exp
#'
#' @title Set the PCA factors for a `single_cell_exp` class.
#'
#' @rdname set_pca_factors
#'
#' @method set_pca_factors single_cell_exp
S7::method(set_pca_factors, single_cell_exp) <- function(
  x,
  pca_factor
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")
  checkmate::assertMatrix(pca_factor, mode = "numeric")

  # add the data using the S3 method
  S7::prop(x, "sc_cache") <- set_pca_factors(
    x = S7::prop(x, "sc_cache"),
    pca_factor = pca_factor
  )

  return(x)
}

#' @name set_pca_loadings.single_cell_exp
#'
#' @title Set the PCA factors for a `single_cell_exp` class.
#'
#' @rdname set_pca_loadings
#'
#' @method set_pca_loadings single_cell_exp
S7::method(set_pca_loadings, single_cell_exp) <- function(
  x,
  pca_loading
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")
  checkmate::assertMatrix(pca_loading, mode = "numeric")

  # add the data using the S3 method
  S7::prop(x, "sc_cache") <- set_pca_loadings(
    x = S7::prop(x, "sc_cache"),
    pca_loading = pca_loading
  )

  return(x)
}

#' @name set_pca_singular_vals.single_cell_exp
#'
#' @title Set the PCA singular values for a `single_cell_exp` class.
#'
#' @rdname set_pca_singular_vals
#'
#' @method set_pca_singular_vals single_cell_exp
S7::method(set_pca_singular_vals, single_cell_exp) <- function(
  x,
  singular_vals
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")
  checkmate::qassert(singular_vals, "N+")

  # add the data using the S3 method
  S7::prop(x, "sc_cache") <- set_pca_singular_vals(
    x = S7::prop(x, "sc_cache"),
    singular_vals = singular_vals
  )

  return(x)
}

#' @name set_embedding.single_cell_exp
#'
#' @title Set additional embeddings to `single_cell_exp`.
#'
#' @rdname set_embedding
#'
#' @method set_embedding single_cell_exp
S7::method(set_embedding, single_cell_exp) <- function(
  x,
  embd,
  name
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")
  checkmate::assertMatrix(embd, mode = "numeric")
  checkmate::qassert(name, "S1")

  # add the data using the S3 method
  S7::prop(x, "sc_cache") <- set_embedding(
    x = S7::prop(x, "sc_cache"),
    embd = embd,
    name = name
  )

  return(x)
}

#' @name set_knn.single_cell_exp
#'
#' @title Set the KNN matrix for a `single_cell_exp` class.
#'
#' @rdname set_knn
#'
#' @method set_knn single_cell_exp
S7::method(set_knn, single_cell_exp) <- function(
  x,
  knn_mat
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")
  checkmate::assertMatrix(knn_mat, mode = "numeric")

  # add the data using the S3 method
  S7::prop(x, "sc_cache") <- set_knn(
    x = S7::prop(x, "sc_cache"),
    knn_mat = knn_mat
  )

  return(x)
}


#' @name set_snn_graph.single_cell_exp
#'
#' @title Set the sNN graph for a `single_cell_exp` class.
#'
#' @rdname set_snn_graph
#'
#' @method set_snn_graph single_cell_exp
S7::method(set_snn_graph, single_cell_exp) <- function(
  x,
  snn_graph
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")
  checkmate::assertClass(snn_graph, "igraph")

  # add the data using the S3 method
  S7::prop(x, "sc_cache") <- set_snn_graph(
    x = S7::prop(x, "sc_cache"),
    snn_graph = snn_graph
  )

  return(x)
}

#' @name remove_knn.single_cell_exp
#'
#' @title Remove the KNN matrix from a `single_cell_exp` class.
#'
#' @rdname remove_knn
#'
#' @method remove_knn single_cell_exp
S7::method(remove_knn, single_cell_exp) <- function(
  x
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")

  # add the data using the S3 method
  S7::prop(x, "sc_cache") <- remove_knn(
    x = S7::prop(x, "sc_cache")
  )

  return(x)
}


#' @name remove_snn_graph.single_cell_exp
#'
#' @title Remove the sNN graph from a `single_cell_exp` class.
#'
#' @rdname remove_snn_graph
#'
#' @method remove_snn_graph single_cell_exp
S7::method(remove_snn_graph, single_cell_exp) <- function(
  x
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")

  # add the data using the S3 method
  S7::prop(x, "sc_cache") <- remove_snn_graph(
    x = S7::prop(x, "sc_cache")
  )

  return(x)
}
