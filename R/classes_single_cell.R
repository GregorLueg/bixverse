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
#' @return Generates an empty version of the `ScMap` class.
#'
#' @export
#'
#' @keywords internal
new_sc_mapper <- function() {
  sc_mapper <- list(
    gene_mapping = NULL,
    cell_mapping = NULL,
    cells_to_keep_idx = NULL,
    hvg_gene_indices = NULL
  )

  class(sc_mapper) <- "ScMap"

  return(sc_mapper)
}

## cache class -----------------------------------------------------------------

#' Helper function to hold relevant cached data
#'
#' @description
#' Helper class that contains various data that is held in memory and not on
#' disk.
#'
#' @return Generates an empty version of the `ScCache` class.
#'
#' @export
#'
#' @keywords internal
new_sc_cache <- function() {
  sc_cache <- list(
    pca_factors = NULL,
    pca_loadings = NULL,
    pca_singular_vals = NULL,
    other_embeddings = list(),
    knn = NULL,
    snn_graph = NULL
  )

  class(sc_cache) <- "ScCache"

  return(sc_cache)
}

## setters ---------------------------------------------------------------------

### methods --------------------------------------------------------------------

#### ScMap ---------------------------------------------------------------------

#' @rdname set_gene_mapping
#'
#' @export
set_gene_mapping.ScMap <- function(x, gene_map) {
  # checks
  checkmate::assertClass(x, "ScMap")
  checkmate::qassert(gene_map, "N+")
  checkmate::assertNamed(gene_map)

  x[["gene_mapping"]] <- gene_map

  return(x)
}

#' @rdname set_cell_mapping
#'
#' @export
set_cell_mapping.ScMap <- function(x, cell_map) {
  # checks
  checkmate::assertClass(x, "ScMap")
  checkmate::qassert(cell_map, "N+")
  checkmate::assertNamed(cell_map)

  x[["cell_mapping"]] <- cell_map

  return(x)
}

#' @rdname set_cells_to_keep
#'
#' @export
set_cells_to_keep.ScMap <- function(x, cells_to_keep) {
  # checks
  checkmate::assertClass(x, "ScMap")
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
set_hvg.ScMap <- function(x, hvg) {
  # checks
  checkmate::assertClass(x, "ScMap")
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

#### ScCache -------------------------------------------------------------------

#' @rdname set_pca_factors
#'
#' @export
set_pca_factors.ScCache <- function(x, pca_factor) {
  # checks
  checkmate::assertClass(x, "ScCache")
  checkmate::assertMatrix(pca_factor, mode = "numeric")

  x[["pca_factors"]] <- pca_factor

  return(x)
}

#' @rdname set_pca_loadings
#'
#' @export
set_pca_loadings.ScCache <- function(x, pca_loading) {
  # checks
  checkmate::assertClass(x, "ScCache")
  checkmate::assertMatrix(pca_loading, mode = "numeric")

  x[["pca_loadings"]] <- pca_loading

  return(x)
}

#' @rdname set_pca_singular_vals
#'
#' @export
set_pca_singular_vals.ScCache <- function(x, singular_vals) {
  # checks
  checkmate::assertClass(x, "ScCache")
  checkmate::qassert(singular_vals, "N+")

  x[["pca_singular_vals"]] <- singular_vals

  return(x)
}

#' @rdname set_embedding
#'
#' @export
set_embedding.ScCache <- function(x, embd, name) {
  # checks
  checkmate::assertClass(x, "ScCache")
  checkmate::assertMatrix(embd, mode = "numeric")
  checkmate::qassert(name, "S1")

  x[["other_embeddings"]][[name]] <- embd

  return(x)
}

#' @rdname set_knn
#'
#' @export
set_knn.ScCache <- function(x, knn) {
  # checks
  checkmate::assertClass(x, "ScCache")
  checkmate::assertClass(knn, "SingleCellNearestNeighbour")

  x[["knn"]] <- knn

  return(x)
}

#' @rdname set_snn_graph
#'
#' @export
set_snn_graph.ScCache <- function(x, snn_graph) {
  # checks
  checkmate::assertClass(x, "ScCache")
  checkmate::assertClass(snn_graph, "igraph")

  x[["snn_graph"]] <- snn_graph

  return(x)
}

#' @rdname remove_knn
#'
#' @export
remove_knn.ScCache <- function(x) {
  # checks
  checkmate::assertClass(x, "ScCache")

  x[["knn"]] <- NULL

  return(x)
}

#' @rdname remove_snn_graph
#'
#' @export
remove_snn_graph.ScCache <- function(x) {
  # checks
  checkmate::assertClass(x, "ScCache")

  x[["snn_graph"]] <- NULL

  return(x)
}

## getters ---------------------------------------------------------------------

### generics -------------------------------------------------------------------

#### ScMap ---------------------------------------------------------------------

### methods --------------------------------------------------------------------

#### ScMap ---------------------------------------------------------------------

#' @rdname get_gene_names
#'
#' @export
get_gene_names.ScMap <- function(x) {
  # checks
  checkmate::assertClass(x, "ScMap")

  return(names(x[["gene_mapping"]]))
}

#' @rdname get_gene_indices
#'
#' @export
get_gene_indices.ScMap <- function(x, gene_ids, rust_index) {
  checkmate::assertClass(x, "ScMap")
  checkmate::qassert(gene_ids, "S+")
  checkmate::qassert(rust_index, "B1")

  gene_map <- x$gene_mapping
  indices <- match(gene_ids, names(gene_map))

  missing <- is.na(indices)
  if (any(missing)) {
    warning(sprintf(
      "With the provided gene_ids a total of %i could not be matched.",
      sum(missing)
    ))
    indices <- indices[!missing]
  }

  if (length(indices) == 0) {
    stop(
      "The gene indices have length 0. Please double check provided parameters!"
    )
  }

  if (rust_index) {
    indices <- indices - 1L
  }

  return(as.integer(indices))
}

#' @rdname get_cell_indices
#'
#' @export
get_cell_indices.ScMap <- function(x, cell_ids, rust_index) {
  checkmate::assertClass(x, "ScMap")
  checkmate::qassert(cell_ids, "S+")
  checkmate::qassert(rust_index, "B1")

  cell_map <- x$cell_mapping
  indices <- match(cell_ids, names(cell_map))

  missing <- is.na(indices)
  if (any(missing)) {
    warning(sprintf(
      "With the provided cell_ids a total of %i could not be matched.",
      sum(missing)
    ))
    indices <- indices[!missing]
  }

  if (length(indices) == 0) {
    stop(
      "The cell indices have length 0. Please double check provided parameters!"
    )
  }

  if (rust_index) {
    indices <- indices - 1L
  }

  return(as.integer(indices))
}

#' @rdname get_cells_to_keep
#'
#' @export
get_cells_to_keep.ScMap <- function(x) {
  # checks
  checkmate::assertClass(x, "ScMap")

  return(as.integer(x[["cells_to_keep_idx"]]))
}

#' @rdname get_cell_names
#'
#' @export
get_cell_names.ScMap <- function(x, filtered = FALSE) {
  # checks
  checkmate::assertClass(x, "ScMap")
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
get_hvg.ScMap <- function(x) {
  # checks
  checkmate::assertClass(x, "ScMap")

  hvg_indices <- as.integer(x[["hvg_gene_indices"]])
  if (length(hvg_indices) == 0) {
    warning("No highly variable features found in the class.")
  }

  return(hvg_indices)
}

#' @rdname get_gene_names_from_idx
#'
#' @export
get_gene_names_from_idx.ScMap <- function(x, gene_idx, rust_based = TRUE) {
  # checks
  checkmate::assertClass(x, "ScMap")
  checkmate::qassert(gene_idx, "I+")
  checkmate::qassert(rust_based, "B1")

  # body
  gene_map <- x[["gene_mapping"]]

  gene_map_rev <- names(gene_map)
  names(gene_map_rev) <- if (rust_based) {
    as.integer(gene_map - 1)
  } else {
    as.integer(gene_map)
  }

  return(gene_map_rev[as.character(gene_idx)])
}

#### ScCache -------------------------------------------------------------------

#' @rdname get_pca_factors
#'
#' @export
get_pca_factors.ScCache <- function(x) {
  # checks
  checkmate::assertClass(x, "ScCache")

  res <- x[["pca_factors"]]

  return(res)
}

#' @rdname get_pca_loadings
#'
#' @export
get_pca_loadings.ScCache <- function(x) {
  # checks
  checkmate::assertClass(x, "ScCache")

  res <- x[["pca_loadings"]]

  return(res)
}

#' @rdname get_pca_singular_val
#'
#' @export
get_pca_singular_val.ScCache <- function(x) {
  # checks
  checkmate::assertClass(x, "ScCache")

  res <- x[["pca_singular_vals"]]

  return(res)
}

#' @rdname get_embedding
#'
#' @export
get_embedding.ScCache <- function(x, embd_name) {
  possible_embedding_names <- c(names(x[["other_embeddings"]]), "pca")

  checkmate::assertClass(x, "ScCache")
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
get_available_embeddings.ScCache <- function(x) {
  checkmate::assertClass(x, "ScCache")

  pca_available <- if (!is.null(x[["pca_factors"]])) "pca" else character(0)
  other_embeddings <- names(x[["other_embeddings"]])

  res <- c(pca_available, other_embeddings)
  if (length(res) == 0L) {
    return("")
  }
  res
}

#' @rdname get_knn_mat
#'
#' @export
get_knn_mat.ScCache <- function(x) {
  # checks
  checkmate::assertClass(x, "ScCache")
  knn <- x[["knn"]]
  if (is.null(knn)) {
    return(NULL)
  }
  get_knn_mat(knn)
}

#' @rdname get_knn_dist
#'
#' @export
get_knn_dist.ScCache <- function(x) {
  # checks
  checkmate::assertClass(x, "ScCache")
  knn <- x[["knn"]]
  if (is.null(knn)) {
    return(NULL)
  }
  get_knn_dist(knn)
}

#' @rdname get_knn_obj
#'
#' @export
get_knn_obj.ScCache <- function(x) {
  # checks
  checkmate::assertClass(x, "ScCache")

  x[["knn"]]
}

#' @rdname get_snn_graph
#'
#' @export
get_snn_graph.ScCache <- function(x) {
  # checks
  checkmate::assertClass(x, "ScCache")

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
#' @return Returns the `SingleCells` class for further operations.
#'
#' @export
SingleCells <- S7::new_class(
  name = "SingleCells",
  properties = list(
    db_connection = S7::class_any,
    count_connection = S7::class_any,
    dir_data = S7::class_character,
    sc_cache = S7::class_any,
    sc_map = S7::class_any,
    dims = S7::class_integer
  ),
  constructor = function(dir_data) {
    # checks
    checkmate::assertDirectoryExists(dir_data)

    dir_data <- path.expand(dir_data)

    # generate the Rust pointer
    count_connection <- SingeCellCountData$new(
      f_path_cells = file.path(dir_data, "counts_cells.bin"),
      f_path_genes = file.path(dir_data, "counts_genes.bin")
    )

    # generate the DuckDB connector
    db_connection <- SingleCellDuckDB$new(
      db_dir = dir_data
    )

    S7::new_object(
      S7::S7_object(),
      db_connection = db_connection,
      count_connection = count_connection,
      dir_data = path.expand(dir_data),
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
#' @param object `SingleCells` class.
#'
#' @return The DuckDB connector
#'
#' @export
#'
#' @keywords internal
get_sc_duckdb <- S7::new_generic(
  name = "get_sc_duckdb",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_sc_duckdb SingleCells
#'
#' @export
S7::method(get_sc_duckdb, SingleCells) <- function(object) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))

  return(S7::prop(object, "db_connection"))
}

#' Getter for the single cell Rust pointer
#'
#' @param object `SingleCells` class.
#'
#' @return The Rust structure
#'
#' @export
#'
#' @keywords internal
get_sc_rust_ptr <- S7::new_generic(
  name = "get_sc_rust_ptr",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_sc_rust_ptr SingleCells
#'
#' @export
S7::method(get_sc_rust_ptr, SingleCells) <- function(object) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))

  return(S7::prop(object, "count_connection"))
}

#' Get the file path to the counts in gene format
#'
#' @param object `SingleCells` class.
#'
#' @return The path to the `counts_genes.bin`
#'
#' @keywords internal
get_rust_count_gene_f_path <- S7::new_generic(
  name = "get_rust_count_gene_f_path",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_rust_count_gene_f_path SingleCells
#'
#' @export
S7::method(get_rust_count_gene_f_path, SingleCells) <- function(object) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))

  f_path <- file.path(S7::prop(object, "dir_data"), "counts_genes.bin")

  return(f_path)
}

#' Get the file path to the counts in cell format
#'
#' @param object `SingleCells` class.
#'
#' @return The path to the `counts_cells.bin`
#'
#' @keywords internal
get_rust_count_cell_f_path <- S7::new_generic(
  name = "get_rust_count_cell_f_path",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_rust_count_cell_f_path SingleCells
#'
#' @export
S7::method(get_rust_count_cell_f_path, SingleCells) <- function(object) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))

  f_path <- file.path(S7::prop(object, "dir_data"), "counts_cells.bin")

  return(f_path)
}

#### duckdb getters ------------------------------------------------------------

#' @method get_sc_obs SingleCells
#'
#' @export
S7::method(get_sc_obs, SingleCells) <- function(
  object,
  indices = NULL,
  cols = NULL,
  filtered = FALSE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
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

#' @method get_sc_var SingleCells
#'
#' @export
S7::method(get_sc_var, SingleCells) <- function(
  object,
  indices = NULL,
  cols = NULL
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(indices, c("0", "I+"))

  duckdb_con <- get_sc_duckdb(object)

  var_table <- duckdb_con$get_vars_table(indices = indices, cols = cols)

  return(var_table)
}

#' @method `[[` SingleCells
#'
#' @export
S7::method(`[[`, SingleCells) <- function(x, i, ...) {
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
#' @param object `SingleCells` class.
#'
#' @returns Returns the maps within the class.
#'
#' @export
#'
#' @keywords internal
get_sc_map <- S7::new_generic(
  name = "get_sc_maps",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_sc_var SingleCells
#'
#' @export
S7::method(get_sc_map, SingleCells) <- function(
  object
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))

  res <- S7::prop(object, "sc_map")

  return(res)
}

### sc cache getter ------------------------------------------------------------

#' Getter the memory-stored data from the class
#'
#' @param object `SingleCells` class.
#'
#' @returns Returns the ScCache class from the object with all the
#' memory-stored data.
#'
#' @export
#'
#' @keywords internal
get_sc_cache <- S7::new_generic(
  name = "get_sc_cache",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_sc_cache SingleCells
#'
#' @export
S7::method(get_sc_cache, SingleCells) <- function(
  object
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))

  res <- S7::prop(object, "sc_cache")

  return(res)
}

#### count getters -------------------------------------------------------------

#' @method get_sc_counts SingleCells
#'
#' @export
S7::method(get_sc_counts, SingleCells) <- function(
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
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
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

#' @method `[` SingleCells
#'
#' @export
S7::method(`[`, SingleCells) <- function(
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
  # transform cell ids and gene ids to indices
  if (checkmate::qtest(i, "S+")) {
    i <- get_cell_indices(x = x, cell_ids = i, rust_index = FALSE)
  }
  if (checkmate::qtest(j, "S+")) {
    j <- get_gene_indices(x = x, gene_ids = j, rust_index = FALSE)
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
#'
#' @keywords internal
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
#'
#' @keywords internal
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
#' @param sc_map A `ScMap` class. Contains various mapping information.
#'
#' @return The finalised matrix.
#'
#' @keywords internal
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
  checkmate::assertClass(sc_map, "ScMap")

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

#### ScMap ---------------------------------------------------------------------

#' @name get_cell_names.SingleCells
#'
#' @title Get the cell names from a `SingleCells`.
#'
#' @rdname get_cell_names
#'
#' @method get_cell_names SingleCells
S7::method(get_cell_names, SingleCells) <- function(
  x,
  filtered = FALSE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))

  # add the data using the S3 method
  cell_names <- get_cell_names(
    x = S7::prop(x, "sc_map"),
    filtered = filtered
  )

  return(cell_names)
}

#' @name get_gene_names.SingleCells
#'
#' @title Get the gene names from a `SingleCells`.
#'
#' @rdname get_gene_names
#'
#' @method get_gene_names SingleCells
S7::method(get_gene_names, SingleCells) <- function(
  x
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))

  # add the data using the S3 method
  gene_names <- get_gene_names(
    x = S7::prop(x, "sc_map")
  )

  return(gene_names)
}

#' @name get_cells_to_keep.SingleCells
#'
#' @title Get the cells to keep from a `SingleCells`.
#'
#' @rdname get_cells_to_keep
#'
#' @method get_cells_to_keep SingleCells
S7::method(get_cells_to_keep, SingleCells) <- function(
  x
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))

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

#' @name get_gene_indices.SingleCells
#'
#' @title Get the gene indices from a `SingleCells`.
#'
#' @rdname get_gene_indices
#'
#' @method get_gene_indices SingleCells
S7::method(get_gene_indices, SingleCells) <- function(
  x,
  gene_ids,
  rust_index
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))
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


#' @name get_cell_indices.SingleCells
#'
#' @title Set the gene mapping for a `SingleCells` class.
#'
#' @rdname get_cell_indices
#'
#' @method get_cell_indices SingleCells
S7::method(get_cell_indices, SingleCells) <- function(
  x,
  cell_ids,
  rust_index
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))
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

#' @name get_hvg.SingleCells
#'
#' @title Get the highly variable gene indices from a `SingleCells`.
#'
#' @rdname get_hvg
#'
#' @method get_hvg SingleCells
S7::method(get_hvg, SingleCells) <- function(
  x
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))

  # forward to S3
  res <- get_hvg(
    x = S7::prop(x, "sc_map")
  )

  return(res)
}

#' @name get_gene_names_from_idx.SingleCells
#'
#' @title Get the highly variable gene indices from a `SingleCells`.
#'
#' @rdname get_gene_names_from_idx
#'
#' @method get_gene_names_from_idx SingleCells
S7::method(get_gene_names_from_idx, SingleCells) <- function(
  x,
  gene_idx,
  rust_based = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))
  checkmate::qassert(gene_idx, "I+")
  checkmate::qassert(rust_based, "B1")

  # forward to S3
  res <- get_gene_names_from_idx(
    x = S7::prop(x, "sc_map"),
    gene_idx = gene_idx,
    rust_based = rust_based
  )

  return(res)
}

#### ScCache -------------------------------------------------------------------

#' @name get_pca_factors.SingleCells
#'
#' @title Get the PCA factors from a `SingleCells`.
#'
#' @rdname get_pca_factors
#'
#' @method get_pca_factors SingleCells
S7::method(get_pca_factors, SingleCells) <- function(
  x
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))

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

#' @name get_pca_loadings.SingleCells
#'
#' @title Get the PCA loadings from a `SingleCells`.
#'
#' @rdname get_pca_loadings
#'
#' @method get_pca_loadings SingleCells
S7::method(get_pca_loadings, SingleCells) <- function(
  x
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))

  # forward to S3
  res <- get_pca_loadings(
    x = S7::prop(x, "sc_cache")
  )

  colnames(res) <- sprintf("PC_%i", 1:ncol(res))
  rownames(res) <- get_gene_names(x)[get_hvg(x) + 1]

  return(res)
}

#' @name get_pca_singular_val.SingleCells
#'
#' @title Get the PCA singular values from a `SingleCells`.
#'
#' @rdname get_pca_singular_val
#'
#' @method get_pca_singular_val SingleCells
S7::method(get_pca_singular_val, SingleCells) <- function(
  x
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))

  # forward to S3
  res <- get_pca_singular_val(
    x = S7::prop(x, "sc_cache")
  )

  return(res)
}

#' @name get_embedding.SingleCells
#'
#' @title Get the embeddings from a `SingleCells`.
#'
#' @rdname get_embedding
#'
#' @method get_embedding SingleCells
S7::method(get_embedding, SingleCells) <- function(
  x,
  embd_name
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))
  checkmate::qassert(embd_name, "S1")

  # forward to S3
  res <- get_embedding(
    x = S7::prop(x, "sc_cache"),
    embd_name = embd_name
  )

  rownames(res) <- get_cell_names(x, filtered = TRUE)

  return(res)
}

#' @name get_available_embeddings.SingleCells
#'
#' @title Get the embeddings from a `SingleCells`.
#'
#' @rdname get_available_embeddings
#'
#' @method get_available_embeddings SingleCells
S7::method(get_available_embeddings, SingleCells) <- function(
  x
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))

  # forward to S3
  res <- get_available_embeddings(
    x = S7::prop(x, "sc_cache")
  )

  return(res)
}

#' @name get_knn_mat.SingleCells
#'
#' @title Get the KNN matrix from a `SingleCells`.
#'
#' @rdname get_knn_mat
#'
#' @method get_knn_mat SingleCells
S7::method(get_knn_mat, SingleCells) <- function(
  x
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))

  # forward to S3
  res <- get_knn_mat(
    x = S7::prop(x, "sc_cache")
  )

  return(res)
}

#' @name get_knn_dist.SingleCells
#'
#' @title Get the KNN distances from a `SingleCells`.
#'
#' @rdname get_knn_dist
#'
#' @method get_knn_dist SingleCells
S7::method(get_knn_dist, SingleCells) <- function(
  x
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))

  # forward to S3
  res <- get_knn_dist(
    x = S7::prop(x, "sc_cache")
  )

  return(res)
}

#' @name get_knn_obj.SingleCells
#'
#' @title Get the KNN object from a `SingleCells`.
#'
#' @rdname get_knn_obj
#'
#' @method get_knn_obj SingleCells
S7::method(get_knn_obj, SingleCells) <- function(
  x
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))

  # forward to S3
  res <- get_knn_obj(
    x = S7::prop(x, "sc_cache")
  )

  return(res)
}

#' @name get_snn_graph.SingleCells
#'
#' @title Get the SNN graph from a `SingleCells`.
#'
#' @rdname get_snn_graph
#'
#' @method get_snn_graph SingleCells
S7::method(get_snn_graph, SingleCells) <- function(
  x
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))

  # forward to S3
  res <- get_snn_graph(
    x = S7::prop(x, "sc_cache")
  )

  return(res)
}

#### available variables -------------------------------------------------------

#' @method get_sc_available_features SingleCells
#'
#' @export
S7::method(get_sc_available_features, SingleCells) <- function(
  object
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))

  # get the duckdb and obs features
  duckdb_con <- get_sc_duckdb(object)
  obs_cols <- duckdb_con$get_obs_cols()
  obs_features <- data.table::data.table(
    feature_name = obs_cols,
    origin = "obs"
  )

  feature_cols <- get_gene_names(object)
  counts_features <- data.table::data.table(
    feature_name = feature_cols,
    origin = "counts"
  )

  final_features <- data.table::rbindlist(list(obs_features, counts_features))

  return(final_features)
}

### setters --------------------------------------------------------------------

#### obs -----------------------------------------------------------------------

#' @method setnames_sc SingleCells
#'
#' @export
S7::method(setnames_sc, SingleCells) <- function(
  object,
  table = c("obs", "var"),
  old,
  new
) {
  table <- match.arg(table)

  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::assertChoice(table, c("obs", "var"))
  checkmate::qassert(old, "S+")
  checkmate::qassert(new, "S+")
  checkmate::assertTRUE(length(old) == length(new))

  duckdb_con <- get_sc_duckdb(object)

  for (i in seq_along(old)) {
    old_i <- old[[i]]
    new_i <- new[[i]]
    duckdb_con$rename_column(table = table, old = old_i, new = new_i)
  }

  invisible(object)
}

#' Add a new column to the obs table
#'
#' @param object `bixverse::SingleCells` class.
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

#' @method set_sc_new_obs_col SingleCells
#'
#' @export
S7::method(set_sc_new_obs_col, SingleCells) <- function(
  object,
  col_name,
  new_data
) {
  cells_to_keep <- (get_cells_to_keep(object) + 1)

  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
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
#' @param object `bixverse::SingleCells` class.
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

#' @method set_sc_new_obs_col_multiple SingleCells
#'
#' @export
S7::method(set_sc_new_obs_col_multiple, SingleCells) <- function(
  object,
  new_data
) {
  cells_to_keep <- (get_cells_to_keep(object) + 1)

  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
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

#' @method `[[<-` SingleCells
#'
#' @export
S7::method(`[[<-`, SingleCells) <- function(x, i, ..., value) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))
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
#' @param object `bixverse::SingleCells` class.
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

#' @method set_sc_new_var_cols SingleCells
#'
#' @export
S7::method(set_sc_new_var_cols, SingleCells) <- function(
  object,
  data_list
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::assertList(data_list, types = "atomic", names = "named")

  new_data <- data.table::as.data.table(data_list)

  duckdb_con <- get_sc_duckdb(object)
  duckdb_con$add_data_var(new_data = new_data)

  return(object)
}

#### sc map --------------------------------------------------------------------

#' @name set_gene_mapping.SingleCells
#'
#' @title Set the gene mapping for a `SingleCells` class.
#'
#' @rdname set_gene_mapping
#'
#' @method set_gene_mapping SingleCells
S7::method(set_gene_mapping, SingleCells) <- function(
  x,
  gene_map
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))
  checkmate::qassert(gene_map, "N+")
  checkmate::assertNamed(gene_map, "named")

  # add the data using the S3 method
  S7::prop(x, "sc_map") <- set_gene_mapping(
    x = S7::prop(x, "sc_map"),
    gene_map = gene_map
  )

  return(x)
}

#' @name set_cell_mapping.SingleCells
#'
#' @title Set the cell mapping for a `SingleCells` class.
#'
#' @rdname set_cell_mapping
#'
#' @method set_cell_mapping SingleCells
S7::method(set_cell_mapping, SingleCells) <- function(
  x,
  cell_map
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))
  checkmate::qassert(cell_map, "N+")
  checkmate::assertNamed(cell_map, "named")

  # add the data using the S3 method
  S7::prop(x, "sc_map") <- set_cell_mapping(
    x = S7::prop(x, "sc_map"),
    cell_map = cell_map
  )

  return(x)
}

#' @name set_cells_to_keep.SingleCells
#'
#' @title Set the cell mapping for a `SingleCells` class.
#'
#' @rdname set_cells_to_keep
#'
#' @method set_cells_to_keep SingleCells
S7::method(set_cells_to_keep, SingleCells) <- function(
  x,
  cells_to_keep
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))
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

#' @name set_hvg.SingleCells
#'
#' @title Set the highly variable genes for a `SingleCells` class.
#'
#' @rdname set_hvg
#'
#' @method set_hvg SingleCells
S7::method(set_hvg, SingleCells) <- function(
  x,
  hvg
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))
  checkmate::qassert(hvg, c("I+", "S+"))

  # add the data using the S3 method
  S7::prop(x, "sc_map") <- set_hvg(
    x = S7::prop(x, "sc_map"),
    hvg = hvg
  )

  return(x)
}

#### sc cache ------------------------------------------------------------------

#' @name set_pca_factors.SingleCells
#'
#' @title Set the PCA factors for a `SingleCells` class.
#'
#' @rdname set_pca_factors
#'
#' @method set_pca_factors SingleCells
S7::method(set_pca_factors, SingleCells) <- function(
  x,
  pca_factor
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))
  checkmate::assertMatrix(pca_factor, mode = "numeric")

  # add the data using the S3 method
  S7::prop(x, "sc_cache") <- set_pca_factors(
    x = S7::prop(x, "sc_cache"),
    pca_factor = pca_factor
  )

  return(x)
}

#' @name set_pca_loadings.SingleCells
#'
#' @title Set the PCA factors for a `SingleCells` class.
#'
#' @rdname set_pca_loadings
#'
#' @method set_pca_loadings SingleCells
S7::method(set_pca_loadings, SingleCells) <- function(
  x,
  pca_loading
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))
  checkmate::assertMatrix(pca_loading, mode = "numeric")

  # add the data using the S3 method
  S7::prop(x, "sc_cache") <- set_pca_loadings(
    x = S7::prop(x, "sc_cache"),
    pca_loading = pca_loading
  )

  return(x)
}

#' @name set_pca_singular_vals.SingleCells
#'
#' @title Set the PCA singular values for a `SingleCells` class.
#'
#' @rdname set_pca_singular_vals
#'
#' @method set_pca_singular_vals SingleCells
S7::method(set_pca_singular_vals, SingleCells) <- function(
  x,
  singular_vals
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))
  checkmate::qassert(singular_vals, "N+")

  # add the data using the S3 method
  S7::prop(x, "sc_cache") <- set_pca_singular_vals(
    x = S7::prop(x, "sc_cache"),
    singular_vals = singular_vals
  )

  return(x)
}

#' @name set_embedding.SingleCells
#'
#' @title Set additional embeddings to `SingleCells`.
#'
#' @rdname set_embedding
#'
#' @method set_embedding SingleCells
S7::method(set_embedding, SingleCells) <- function(
  x,
  embd,
  name
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))
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

#' @name set_knn.SingleCells
#'
#' @title Set the KNN matrix for a `SingleCells` class.
#'
#' @rdname set_knn
#'
#' @method set_knn SingleCells
S7::method(set_knn, SingleCells) <- function(
  x,
  knn
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))
  checkmate::assertClass(knn, "SingleCellNearestNeighbour")

  # add the data using the S3 method
  S7::prop(x, "sc_cache") <- set_knn(
    x = S7::prop(x, "sc_cache"),
    knn = knn
  )

  return(x)
}


#' @name set_snn_graph.SingleCells
#'
#' @title Set the sNN graph for a `SingleCells` class.
#'
#' @rdname set_snn_graph
#'
#' @method set_snn_graph SingleCells
S7::method(set_snn_graph, SingleCells) <- function(
  x,
  snn_graph
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))
  checkmate::assertClass(snn_graph, "igraph")

  # add the data using the S3 method
  S7::prop(x, "sc_cache") <- set_snn_graph(
    x = S7::prop(x, "sc_cache"),
    snn_graph = snn_graph
  )

  return(x)
}

#' @name remove_knn.SingleCells
#'
#' @title Remove the KNN matrix from a `SingleCells` class.
#'
#' @rdname remove_knn
#'
#' @method remove_knn SingleCells
S7::method(remove_knn, SingleCells) <- function(
  x
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))

  # add the data using the S3 method
  S7::prop(x, "sc_cache") <- remove_knn(
    x = S7::prop(x, "sc_cache")
  )

  return(x)
}


#' @name remove_snn_graph.SingleCells
#'
#' @title Remove the sNN graph from a `SingleCells` class.
#'
#' @rdname remove_snn_graph
#'
#' @method remove_snn_graph SingleCells
S7::method(remove_snn_graph, SingleCells) <- function(
  x
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))

  # add the data using the S3 method
  S7::prop(x, "sc_cache") <- remove_snn_graph(
    x = S7::prop(x, "sc_cache")
  )

  return(x)
}

### generic / primitives -------------------------------------------------------

#' @name print.SingleCells
#' @title print Method for SingleCells object
#'
#' @description
#' Print a SingleCells object.
#'
#' @param x An object of class `SingleCells`.
#' @param ... Additional arguments (currently not used).
#'
#' @returns Invisibly returns `x`.
#'
#' @method print SingleCells
#'
#' @keywords internal
S7::method(print, SingleCells) <- function(x, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCells))

  dims <- S7::prop(x, "dims")
  no_cells_to_keep <- length(get_cells_to_keep(x))
  sc_map <- S7::prop(x, "sc_map")
  sc_cache <- S7::prop(x, "sc_cache")

  hvg_calculated <- !is.null(sc_map[["hvg_gene_indices"]])
  pca_calculated <- !is.null(sc_cache[["pca_factors"]])
  knn_generated <- !is.null(sc_cache[["knn"]])
  snn_generated <- !is.null(sc_cache[["snn_graph"]])

  other_embeddings <- names(sc_cache[["other_embeddings"]])
  other_embeddings_str <- if (length(other_embeddings) > 0) {
    paste(other_embeddings, collapse = ", ")
  } else {
    "none"
  }

  cat(
    "Single cell experiment (Single Cells).\n",
    sprintf("  No cells (original): %i\n", dims[1]),
    sprintf("   To keep n: %i\n", no_cells_to_keep),
    sprintf("  No genes: %i\n", dims[2]),
    sprintf("  HVG calculated: %s\n", hvg_calculated),
    sprintf("  PCA calculated: %s\n", pca_calculated),
    sprintf("  Other embeddings: %s\n", other_embeddings_str),
    sprintf("  KNN generated: %s\n", knn_generated),
    sprintf("  SNN generated: %s\n", snn_generated),
    sep = ""
  )

  invisible(x)
}

#' @name dim.SingleCells
#' @title dim Method for SingleCells object
#'
#' @description
#' Returns the dimensions of a SingleCells object.
#'
#' @param x An object of class `SingleCells`.
#'
#' @returns An integer vector of length 2 with the number of cells and genes.
#'
#' @method dim SingleCells
#'
#' @keywords internal
S7::method(dim, SingleCells) <- function(x) {
  S7::prop(x, "dims")
}
