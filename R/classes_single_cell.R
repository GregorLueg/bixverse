# s3 ---------------------------------------------------------------------------

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
    other_embeddings = list(),
    knn_matrix = NULL,
    snn_matrix = NULL
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
set_cell_to_keep <- function(x, cells_to_keep) {
  UseMethod("set_cell_to_keep")
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

### methods --------------------------------------------------------------------

#### sc_mapper -----------------------------------------------------------------

#' Set gene mapping method for sc_mapper
#'
#' @param x An `sc_mapper` object
#' @param gene_map Named integer vector for gene mapping
#'
#' @return Updated `sc_mapper` object with gene mapping set
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

#' Set cell mapping method for sc_mapper
#'
#' @param x An `sc_mapper` object
#' @param cell_map Named integer vector for cell mapping
#'
#' @return Updated `sc_mapper` object with cell mapping set
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

#' Set cells to keep for sc_mapper
#'
#' @param x An `sc_mapper` object
#' @param cells_to_keep String or integer. The names or indices of the cells
#' to keep in downstream analysis.
#'
#' @return Updated `sc_mapper` object with column index mapping set
#'
#' @export
set_cell_to_keep.sc_mapper <- function(x, cells_to_keep) {
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

#' Set the HVG for sc_mapper
#'
#' @param x An `sc_mapper` object
#' @param hvg String or integer. The names or indices of the highly variable
#' genes.
#'
#' @return Updated `sc_mapper` object with column index mapping set
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

#' Set PCA factors for sc_cache
#'
#' @param x An `sc_cache` object
#' @param pca_factor Numerical matrix. The matrix with the PCA factors.
#'
#' @return Updated `sc_cache` object with added PCA factors.
#'
#' @export
set_pca_factors.sc_cache <- function(x, pca_factor) {
  # checks
  checkmate::assertClass(x, "sc_cache")
  checkmate::assertMatrix(pca_factor, mode = "numeric")

  x[["pca_factors"]] <- pca_factor

  return(x)
}

#' Set PCA loadings for sc_mapper
#'
#' @param x An `sc_cache` object
#' @param pca_loading Numerical matrix. The matrix with the PCA factors.
#'
#' @return Updated `sc_cache` object with added PCA loadings.
#'
#' @export
set_pca_loadings.sc_cache <- function(x, pca_loading) {
  # checks
  checkmate::assertClass(x, "sc_cache")
  checkmate::assertMatrix(pca_loading, mode = "numeric")

  x[["pca_loadings"]] <- pca_loading

  return(x)
}

## getters ---------------------------------------------------------------------

### generics -------------------------------------------------------------------

#### sc_mapper -----------------------------------------------------------------

#' Get the gene names
#'
#' @param x An object to get the gene names for
#'
#' @export
get_gene_names <- function(x) {
  UseMethod("get_gene_names")
}

#' Get the cell names
#'
#' @param x An object to get the cell names for
#'
#' @export
get_cell_names <- function(x) {
  UseMethod("get_cell_names")
}

#' Get the index position for a gene
#'
#' @param x An object to get the gene index for.
#' @param gene_ids String vector. The gene ids to search for.
#' @param rust_index Bool. Shall rust-based indexing be returned.
get_gene_indices <- function(x, gene_ids, rust_index) {
  UseMethod("get_gene_indices")
}

#' Get the cells to keep
#'
#' @param x An object to get the gene index for.
get_cells_to_keep <- function(x) {
  UseMethod("get_cells_to_keep")
}

#' Get the HVG
#'
#' @param x An object to get HVG for.
get_hvg <- function(x) {
  UseMethod("get_hvg")
}

#### sc_cache ------------------------------------------------------------------

#' Get the PCA factors
#'
#' @param x An object to get PCA factors for.
get_pca_factors <- function(x) {
  UseMethod("get_pca_factors")
}

#' Get the PCA loadings
#'
#' @param x An object to get PCA loadings for.
get_pca_loadings <- function(x) {
  UseMethod("get_pca_loadings")
}

### methods --------------------------------------------------------------------

#### sc_mapper -----------------------------------------------------------------

#' Get the gene names
#'
#' @param x An `sc_mapper` object
#'
#' @return The gene names
#'
#' @export
get_gene_names.sc_mapper <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_mapper")

  return(names(x[["gene_mapping"]]))
}

#' Get the cell names
#'
#' @param x An `sc_mapper` object
#' @param gene_ids String vector. The gene ids to search for.
#' @param rust_index Bool. Shall rust-based indexing be returned.
#'
#' @return The gene indices
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

#' Get the gene indices
#'
#' @param x An `sc_mapper` object
#'
#' @return The cell names
#'
#' @export
get_cells_to_keep.sc_mapper <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_mapper")

  return(as.integer(x[["cells_to_keep_idx"]]))
}

#' Get the indices of the cells to keep
#'
#' @param x An `sc_mapper` object
#'
#' @return The cell indices (0-based for Rust)
#'
#' @export
get_cell_names.sc_mapper <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_mapper")

  return(names(x[["cell_mapping"]]))
}

#' Get the indices of the HVG
#'
#' @param x An `sc_mapper` object
#'
#' @return The gene indices (0-based for Rust) of the HVG
#'
#' @export
get_hvg.sc_mapper <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_mapper")

  return(as.integer(x[["hvg_gene_indices"]]))
}

#### sc_cache ------------------------------------------------------------------

#' Get the PCA factors
#'
#' @param x An `sc_cache` object
#'
#' @return The PCA factor matrix.
#'
#' @export
get_pca_factors.sc_cache <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_cache")

  return(x[["pca_factors"]])
}

#' Get the PCA loadings
#'
#' @param x An `sc_cache` object
#'
#' @return The PCA loading matrix.
#'
#' @export
get_pca_loadings.sc_cache <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_cache")

  return(x[["pca_loadings"]])
}

# s7 ---------------------------------------------------------------------------

## single cell class -----------------------------------------------------------

#' @title bixverse single cell class (nightly!)
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
#'   \item{sc_cache}{List with cached data. Future feature, nothing implemented
#'   yet.}
#'   \item{dims}{Dimensions of the original data.}
#'   \item{index_maps}{A list of two named numerics that contains cell id to
#'   cell idx and gene id to gene idx info.}
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
    nightly_feature()
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

#' Getter the obs table
#'
#' @param object `single_cell_exp` class.
#' @param indices Optional integer vector. The integer positions of the cells
#' to return.
#' @param cols Optional string vector. The columns from the obs table to return.
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
    cols = NULL
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_sc_obs single_cell_exp
#'
#' @export
S7::method(get_sc_obs, single_cell_exp) <- function(
  object,
  indices = NULL,
  cols = NULL
) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::qassert(indices, c("0", "I+"))

  duckdb_con <- get_sc_duckdb(object)

  obs_table <- duckdb_con$get_obs_table(indices = indices, cols = cols)

  return(obs_table)
}

#' Getter the var table
#'
#' @param object `single_cell_exp` class.
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
    get_sc_obs(x, cols = i)
  } else if (checkmate::qtest(i, "I+")) {
    get_sc_obs(x, indices = i)
  } else if (checkmate::qtest(i, "0")) {
    get_sc_obs(x)
  } else {
    stop("Invalid type")
  }
}

### map getters ----------------------------------------------------------------

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

#### count getters -------------------------------------------------------------

#' Getter the counts
#'
#' @param object `single_cell_exp` class.
#' @param assay String. Which slot to return. One of `c("raw", "norm")`.
#' Defaults to `"raw"`.
#' @param return_format String. One of `c("cell", "gene")`. Return data in
#' cell-centric compressed format (CSR) or gene-centric compressed format (CSC).
#' Defaults to `"cell"`.
#' @param cell_indices Optional cell indices.
#' @param gene_indices Optional gene indices.
#'
#' @return The obs table
#'
#' @export
get_sc_counts <- S7::new_generic(
  name = "get_sc_obs",
  dispatch_args = "object",
  fun = function(
    object,
    assay = c("raw", "norm"),
    return_format = c("cell", "gene"),
    cell_indices = NULL,
    gene_indices = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_sc_counts single_cell_exp
#'
#' @export
S7::method(get_sc_counts, single_cell_exp) <- function(
  object,
  assay = c("raw", "norm"),
  return_format = c("cell", "gene"),
  cell_indices = NULL,
  gene_indices = NULL,
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
#' @title Get cell names for `single_cell_exp`
#'
#' @method get_cell_names single_cell_exp
S7::method(get_cell_names, single_cell_exp) <- function(
  x
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")

  # add the data using the S3 method
  cell_names <- get_cell_names(
    x = S7::prop(x, "sc_map")
  )

  return(cell_names)
}

#' @name get_gene_names.single_cell_exp
#'
#' @title Get gene names for `single_cell_exp`
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
#' @title Get gene indices for `single_cell_exp`
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

  return(res)
}

#' @name get_gene_indices.single_cell_exp
#'
#' @title Get gene indices for `single_cell_exp`
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

#' @name get_hvg.single_cell_exp
#'
#' @title Get HVG gene indices for `single_cell_exp`
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
#' @title Get the PCA factors for `single_cell_exp`
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

  return(res)
}

#' @name get_pca_loadings.single_cell_exp
#'
#' @title Get the PCA loadings for `single_cell_exp`
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

  return(res)
}

### setters --------------------------------------------------------------------

#### obs -----------------------------------------------------------------------

#' Add a new column to the obs table
#'
#' @param object `bixverse::single_cell_exp` class.
#' @param col_name String. The name of the column to add.
#' @param new_data Atomic vector. The data to add to the column. Needs to be
#' of same length as nrow obs table.
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
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::qassert(col_name, "S1")
  checkmate::qassert(new_data, "a")

  duckdb_con <- get_sc_duckdb(object)

  new_data <- data.table::data.table(new_data)
  data.table::setnames(new_data, "new_data", col_name)

  duckdb_con$add_data_obs(new_data = new_data)

  return(object)
}

#' @method `[[<-` single_cell_exp
#'
#' @export
S7::method(`[[<-`, single_cell_exp) <- function(x, i, value) {
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
    for (j in seq_along(i)) {
      col_name_j <- names(value)[j]
      new_data_j <- value[[j]]
      x <- set_sc_new_obs_col(
        object = x,
        col_name = col_name_j,
        new_data = new_data_j
      )
    }
  }

  return(x)
}

#' Add a new column to the obs table
#'
#' @param object `bixverse::single_cell_exp` class.
#' @param data_list Named list with the data to add.
#'
#' @return The class with updated obs table in the DuckDB
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

  duckdb_con <- get_sc_duckdb(object)

  for (i in seq_along(data_list)) {
    new_data_i <- data.table::data.table(new_data = data_list[[i]])
    data.table::setnames(new_data_i, "new_data", names(data_list)[i])
    duckdb_con$add_data_var(new_data = new_data_i)
  }

  return(object)
}

#### sc map --------------------------------------------------------------------

#' @name set_gene_mapping.single_cell_exp
#'
#' @title Set gene mapping method for `single_cell_exp`
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
#' @title Set cell mapping method for `single_cell_exp`
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

#' @name set_cell_to_keep.single_cell_exp
#'
#' @title Set cell method for `single_cell_exp`
#'
#' @method set_cell_to_keep single_cell_exp
S7::method(set_cell_to_keep, single_cell_exp) <- function(
  x,
  cells_to_keep
) {
  # checks
  checkmate::assertClass(x, "bixverse::single_cell_exp")
  checkmate::qassert(cells_to_keep, c("I+", "S+"))

  # add the data using the S3 method
  S7::prop(x, "sc_map") <- set_cell_to_keep(
    x = S7::prop(x, "sc_map"),
    cells_to_keep = cells_to_keep
  )

  return(x)
}

#' @name set_hvg.single_cell_exp
#'
#' @title Set HVG genes for `single_cell_exp`
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
#' @title Set PCA factors method for `single_cell_exp`
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
#' @title Set PCA factors method for `single_cell_exp`
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
