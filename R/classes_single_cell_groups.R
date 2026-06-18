# single cell with grouping class --------------------------------------------------------------

## s7 object -------------------------------------------------------------------

#' @title bixverse single cell with grouping class
#'
#' @description
#' This is the `bixverse`-based single cell class with grouping variable.
#' It is generated internally when requesting grouping based meta-cell generating
#' functions on top of the main single cell class.
#' The class is derived from the original single cell class, pointing
#' to the original count matrix. To prevent issues with state management
#' of the single cell object, the obs and var are kept in memory.
#'
#' @param dir_data String. This is the directory in which the experimental files
#' will be saved
#' @param obs_data data.table with the original obs table.
#' @param var_data data.table with the variable/feature informations.
#' @param grouping String describing the grouping variable.
#'
#' @section Properties:
#' \describe{
#'   \item{count_connection}{This contains an R6-like environment that points
#'   to Rust functions that can work on the RNAseq counts more specifically.}
#'   \item{obs_table}{The single cell observation table.}
#'   \item{var_table}{The single cell variable table.}
#'   \item{grouping_variable}{The grouping variable.}
#'   \item{sc_cache}{Class with embeddings, kNN/sNN graphs, etc. Shared with
#'   [bixverse::SingleCells()].}
#'   \item{dims}{Dimensions of the original data.}
#'   \item{other_data}{List that contains additional data and results, such
#'   as for example the hvg indices}
#' }
#'
#' @return Returns the `SingleCellsGrps` class for further operations.
#'
#' @export
SingleCellsGrps <- S7::new_class(
  name = "SingleCellsGrps",
  properties = list(
    count_connection = S7::class_any,
    obs_table = S7::class_data.frame,
    var_table = S7::class_data.frame,
    groups = S7::class_character,
    sc_cache = S7::class_any,
    sc_map = S7::class_any,
    dims = S7::class_integer,
    other_data = S7::class_list
  ),
  constructor = function(
    sc_object,
    grouping_variable,
    other_data = list()
  ) {
    # checks
    checkmate::assertClass(sc_object, "bixverse::SingleCells")

    dir_data <- S7::prop(sc_object, "dir_data")
    sc_map <- S7::prop(sc_object, "sc_map")
    obs_table <- get_sc_obs(sc_object)
    var_table <- get_sc_var(sc_object)

    checkmate::assertString(grouping_variable)
    checkmate::assertNames(
      x = colnames(obs_table),
      must.include = grouping_variable
    )
    # function body
    # generate the Rust pointer
    count_connection <- SingleCellCountData$new(
      f_path_cells = file.path(dir_data, "counts_cells.bin"),
      f_path_genes = file.path(dir_data, "counts_genes.bin")
    )
    #
    groups = unique(sort(obs_table[[grouping_variable]]))
    #
    cells_by_group <- split(
      obs_table[, .(cell_idx, cell_id)],
      obs_table[[grouping_variable]]
    )
    cells_by_group <- purrr::map(
      cells_by_group,
      ~ {
        setNames(.x$cell_idx, .x$cell_id)
      }
    )
    sc_map$cell_mapping <- cells_by_group

    S7::new_object(
      S7::S7_object(),
      count_connection = count_connection,
      obs_table = obs_table,
      var_table = var_table,
      groups = groups,
      sc_cache = new_sc_cache(),
      sc_map = sc_map,
      dims = c(nrow(obs_table), nrow(var_table)),
      other_data = other_data
    )
  }
)

## primitives ------------------------------------------------------------------

#' @name print.SingleCellsGrps
#'
#' @title print Method for SingleCellsGrps object
#'
#' @description
#' Print a SingleCellsGrps object.
#'
#' @param x An object of class `SingleCellsGrps`.
#' @param ... Additional arguments (currently not used).
#'
#' @returns Invisibly returns `x`.
#'
#' @method print SingleCellsGrps
#'
#' @keywords internal
S7::method(print, SingleCellsGrps) <- function(x, ...) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))

  # various params
  dims <- S7::prop(x, "dims")
  ngroups <- S7::prop(x, "groups")
  sc_cache <- S7::prop(x, "sc_cache")
  other_data <- S7::prop(x, "other_data")
  hvg_calculated <- !is.null(other_data[["hvg"]])
  pca_calculated <- !is.null(sc_cache[["pca_factors"]])
  knn_generated <- !is.null(sc_cache[["knn"]])
  snn_generated <- !is.null(sc_cache[["snn_graph"]])
  other_embeddings <- names(sc_cache[["other_embeddings"]])
  other_embeddings_str <- if (length(other_embeddings) > 0) {
    paste(other_embeddings, collapse = ", ")
  } else {
    "none"
  }

  # print
  cat(
    "Single cell experiment (Groups).\n",
    sprintf("  No cells: %i\n", dims[1]),
    sprintf("  No genes: %i\n", dims[2]),
    sprintf("  No of groups: %i\n", length(ngroups)),
    sprintf("    Group names: %s \n", paste(ngroups, collapse = ", ")),
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
#'
#' @title dim Method for SingleCells object
#'
#' @description
#' Returns the dimensions of a SingleCells object. (Only taking into
#' consideration the cells to keep).
#'
#' @param x An object of class `SingleCells`.
#'
#' @returns An integer vector of length 2 with the number of cells and genes.
#'
#' @method dim SingleCells
#'
#' @keywords internal
S7::method(dim, SingleCellsGrps) <- function(x) {
  n_cells <- nrow(S7::prop(x, "obs_table"))
  n_genes <- nrow(S7::prop(x, "var_table"))
  c(n_cells, n_genes)
}

#' @name head.SingleCellsGrps
#'
#' @title head Method for SingleCellsGrps object
#'
#' @description
#' Returns the first `n` rows of the obs table from a `SingleCellsGrps` object.
#'
#' @param x An object of class `SingleCellsGrps`.
#' @param n Integer. Number of rows to return. Defaults to `6L`.
#' @param ... Additional arguments (currently not used).
#'
#' @returns A data.table with the first `n` rows of the obs table.
#'
#' @method head SingleCellsGrps
#'
#' @keywords internal
S7::method(head, SingleCellsGrps) <- function(x, n = 6L, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))
  checkmate::qassert(n, "I1[1,)")
  get_sc_obs(x, indices = seq_len(n), filtered = FALSE)
}

## getters ---------------------------------------------------------------------

### obs ------------------------------------------------------------------------

#' @method get_sc_obs SingleCellsGrps
#'
#' @export
S7::method(get_sc_obs, SingleCellsGrps) <- function(
  object,
  indices = NULL,
  cols = NULL,
  filtered = FALSE
) {
  checkmate::assertClass(object, "bixverse::SingleCellsGrps")
  checkmate::qassert(indices, c("0", "I+"))
  checkmate::qassert(filtered, "B1")

  obs_table <- data.table::copy(S7::prop(object, "obs_table"))

  if (!is.null(indices)) {
    obs_table <- obs_table[indices, ]
  }

  if (!is.null(cols)) {
    obs_table <- obs_table[, cols, with = FALSE]
  }

  return(obs_table)
}

#' @method `[[` SingleCellsGrps
#'
#' @export
S7::method(`[[`, SingleCellsGrps) <- function(x, i, ...) {
  if (missing(i)) {
    i <- NULL
  }
  if (checkmate::qtest(i, "S+")) {
    get_sc_obs(x, cols = i, filtered = FALSE)
  } else if (checkmate::qtest(i, "I+")) {
    get_sc_obs(x, indices = i, filtered = FALSE)
  } else if (checkmate::qtest(i, "0")) {
    get_sc_obs(x, filtered = FALSE)
  } else {
    stop("Invalid type")
  }
}

### vars -----------------------------------------------------------------------

#' @method get_sc_var SingleCellsGrps
#'
#' @export
S7::method(get_sc_var, SingleCellsGrps) <- function(
  object,
  indices = NULL,
  cols = NULL,
  modality = c("rna", "adt")
) {
  modality <- match.arg(modality)

  checkmate::assertClass(object, "bixverse::SingleCellsGrps")
  checkmate::qassert(indices, c("0", "I+"))

  if (modality != "rna") {
    stop(
      paste(
        "SingleCellsGrps only supports modality = 'rna'.",
        "Use SingleCellsMultiModal for ADT."
      )
    )
  }

  var_table <- object@var_table

  if (!is.null(indices)) {
    var_table <- var_table[indices, ]
  }

  if (!is.null(cols)) {
    var_table <- var_table[, cols, with = FALSE]
  }

  return(var_table)
}

### counts ---------------------------------------------------------------------

#' @method get_sc_counts SingleCellsGrps
#'
#' @export
S7::method(get_sc_counts, SingleCellsGrps) <- function(
  object,
  group = NULL,
  assay = c("raw", "norm"),
  return_format = c("cell", "gene"),
  cell_indices = NULL,
  gene_indices = NULL,
  use_cells_to_keep = TRUE,
  modality = c("rna", "adt"),
  .verbose = TRUE
) {
  assay <- match.arg(assay)
  return_format <- match.arg(return_format)
  modality <- match.arg(modality)
  if (modality != "rna") {
    stop(
      paste(
        "SingleCellsGrps only supports modality = 'rna'.",
        "Use SingleCellsMultiModal for ADT."
      )
    )
  }

  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsGrps))
  checkmate::assertChoice(group, S7::prop(object, "groups"))
  checkmate::assertChoice(assay, c("raw", "norm"))
  checkmate::assertChoice(return_format, c("cell", "gene"))
  checkmate::qassert(cell_indices, c("0", "I+"))
  checkmate::qassert(gene_indices, c("0", "I+"))
  checkmate::qassert(.verbose, "B1")

  requireNamespace("Matrix", quietly = TRUE)

  rust_con <- get_sc_rust_ptr(object)

  sc_map <- get_sc_map(object)
  ## tmp fix, need to think about a better implementation
  sc_map$cell_mapping <- sc_map$cell_mapping[[group]]

  cells_to_keep <- intersect(get_cells_to_keep(object), sc_map$cell_mapping)

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

### map getter -----------------------------------------------------------------

#' @method get_sc_var SingleCellsGrps
#'
#' @export
S7::method(get_sc_map, SingleCellsGrps) <- function(
  object
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsGrps))

  res <- S7::prop(object, "sc_map")

  return(res)
}

#' @name get_cells_to_keep.SingleCellsGrps
#'
#' @rdname get_cells_to_keep
#'
#' @method get_cells_to_keep SingleCellsGrps
S7::method(get_cells_to_keep, SingleCellsGrps) <- function(
  x
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))

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

#### env getters ---------------------------------------------------------------

#' @method get_sc_rust_ptr SingleCellsGrps
#'
#' @export
S7::method(get_sc_rust_ptr, SingleCellsGrps) <- function(object) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsGrps))

  return(S7::prop(object, "count_connection"))
}

#' @method get_rust_count_gene_f_path SingleCellsGrps
#'
#' @export
S7::method(get_rust_count_gene_f_path, SingleCellsGrps) <- function(object) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsGrps))

  f_path <- file.path(S7::prop(object, "dir_data"), "counts_genes.bin")

  return(f_path)
}


#' @method get_rust_count_cell_f_path SingleCellsGrps
#'
#' @export
S7::method(get_rust_count_cell_f_path, SingleCellsGrps) <- function(object) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsGrps))

  f_path <- file.path(S7::prop(object, "dir_data"), "counts_cells.bin")

  return(f_path)
}

#' @method `[` SingleCellsGrps
#'
#' @export
S7::method(`[`, SingleCellsGrps) <- function(
  x,
  i,
  j,
  ...,
  assay = c("raw", "norm"),
  drop = TRUE
) {
  if (missing(i)) {
    i <- NULL
  }
  if (missing(j)) {
    j <- NULL
  }
  assay <- match.arg(assay)

  # transform (meta)cell ids and gene ids to indices
  if (checkmate::qtest(i, "S+")) {
    i <- get_cell_indices(x = x, cell_ids = i, rust_index = FALSE)
  }
  if (checkmate::qtest(j, "S+")) {
    j <- get_gene_indices(x = x, gene_ids = j, rust_index = FALSE)
  }

  checkmate::qassert(i, c("I+", "0"))
  checkmate::qassert(j, c("I+", "0"))
  checkmate::assertChoice(assay, c("raw", "norm"))

  get_sc_counts(
    object = x,
    assay = assay,
    cell_indices = i,
    gene_indices = j
  )
}


## setters ---------------------------------------------------------------------

### obs table ------------------------------------------------------------------

#' @method `[[<-` SingleCellsGrps
#'
#' @export
S7::method(`[[<-`, SingleCellsGrps) <- function(x, i, ..., value) {
  checkmate::assertClass(x, "bixverse::SingleCellsGrps")
  checkmate::qassert(i, "S+")

  if (length(i) == 1) {
    checkmate::qassert(value, "a")
    S7::prop(x, "obs_table")[, (i) := value]
  } else {
    checkmate::assertList(value, names = "named", types = "atomic")
    S7::prop(x, "obs_table")[, (i) := value]
  }

  return(x)
}

### others ---------------------------------------------------------------------

#### setters -------------------------------------------------------------------

#' @name set_hvg.SingleCellsGrps
#'
#' @rdname set_hvg
#'
#' @method set_hvg SingleCellsGrps
S7::method(set_hvg, SingleCellsGrps) <- function(
  x,
  hvg
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))
  checkmate::qassert(hvg, c("I+", "S+"))

  var_table <- S7::prop(x, "var_table")

  res <- if (is.numeric(hvg)) {
    as.integer(hvg)
  } else {
    idx <- match(hvg, var_table$gene_id)
    missing <- is.na(idx)
    if (any(missing)) {
      warning(sprintf(
        "With the provided hvg gene_ids a total of %i could not be matched.",
        sum(missing)
      ))
      idx <- idx[!missing]
    }
    if (length(idx) == 0) {
      stop(
        "The HVG indices have length 0. Please double check provided parameters!"
      )
    }
    as.integer(idx)
  }

  other_data <- S7::prop(x, "other_data")
  other_data[["hvg"]] <- res
  S7::prop(x, "other_data") <- other_data

  return(x)
}

#### getters -------------------------------------------------------------------

#' @name get_hvg.SingleCellsGrps
#'
#' @rdname get_hvg
#'
#' @method get_hvg SingleCellsGrps
S7::method(get_hvg, SingleCellsGrps) <- function(
  x
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))

  hvg_indices <- as.integer(S7::prop(x, "other_data")[["hvg"]])
  if (length(hvg_indices) == 0) {
    warning("No highly variable features found in the class.")
  }

  return(hvg_indices)
}

#' @name get_cell_indices.SingleCellsGrps
#'
#' @rdname get_cell_indices
#'
#' @method get_cell_indices SingleCellsGrps
S7::method(get_cell_indices, SingleCellsGrps) <- function(
  x,
  cell_ids,
  rust_index
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))
  checkmate::qassert(cell_ids, "S+")
  checkmate::qassert(rust_index, "B1")

  indices <- match(cell_ids, x@obs_table[["meta_cell_id"]])

  missing <- is.na(indices)
  if (any(missing)) {
    warning(sprintf(
      "With the provided meta_cell_ids a total of %i could not be matched.",
      sum(missing)
    ))
    indices <- indices[!missing]
  }

  return(indices)
}

#' @name get_gene_indices.SingleCellsGrps
#'
#' @rdname get_gene_indices
#'
#' @method get_gene_indices SingleCellsGrps
S7::method(get_gene_indices, SingleCellsGrps) <- function(
  x,
  gene_ids,
  rust_index
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))
  checkmate::qassert(gene_ids, "S+")
  checkmate::qassert(rust_index, "B1")

  indices <- match(gene_ids, x@var_table[["gene_id"]])

  missing <- is.na(indices)
  if (any(missing)) {
    warning(sprintf(
      "With the provided gene_ids a total of %i could not be matched.",
      sum(missing)
    ))
    indices <- indices[!missing]
  }

  if (rust_index) {
    indices <- indices - 1L
  }

  if (length(indices) == 0) {
    stop(
      "The gene indices have length 0. Please double check provided parameters!"
    )
  }

  return(indices)
}

### sc cache -------------------------------------------------------------------

#### setters -------------------------------------------------------------------

#' @name set_pca_factors.SingleCellsGrps
#'
#' @rdname set_pca_factors
#'
#' @method set_pca_factors SingleCellsGrps
S7::method(set_pca_factors, SingleCellsGrps) <- function(
  x,
  pca_factor,
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))
  checkmate::assertMatrix(pca_factor, mode = "numeric")

  S7::prop(x, "sc_cache") <- set_pca_factors(
    x = S7::prop(x, "sc_cache"),
    pca_factor = pca_factor
  )

  return(x)
}

#' @name set_pca_loadings.SingleCellsGrps
#'
#' @rdname set_pca_loadings
#'
#' @method set_pca_loadings SingleCellsGrps
S7::method(set_pca_loadings, SingleCellsGrps) <- function(
  x,
  pca_loading,
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))
  checkmate::assertMatrix(pca_loading, mode = "numeric")

  S7::prop(x, "sc_cache") <- set_pca_loadings(
    x = S7::prop(x, "sc_cache"),
    pca_loading = pca_loading
  )

  return(x)
}

#' @name set_pca_singular_vals.SingleCellsGrps
#'
#' @rdname set_pca_singular_vals
#'
#' @method set_pca_singular_vals SingleCellsGrps
S7::method(set_pca_singular_vals, SingleCellsGrps) <- function(
  x,
  singular_vals,
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))
  checkmate::qassert(singular_vals, "N+")

  S7::prop(x, "sc_cache") <- set_pca_singular_vals(
    x = S7::prop(x, "sc_cache"),
    singular_vals = singular_vals
  )

  return(x)
}

#' @name set_embedding.SingleCellsGrps
#'
#' @rdname set_embedding
#'
#' @method set_embedding SingleCellsGrps
S7::method(set_embedding, SingleCellsGrps) <- function(
  x,
  embd,
  name,
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))
  checkmate::assertMatrix(embd, mode = "numeric")
  checkmate::qassert(name, "S1")

  S7::prop(x, "sc_cache") <- set_embedding(
    x = S7::prop(x, "sc_cache"),
    embd = embd,
    name = name
  )

  return(x)
}

#' @name set_knn.SingleCellsGrps
#'
#' @rdname set_knn
#'
#' @method set_knn SingleCellsGrps
S7::method(set_knn, SingleCellsGrps) <- function(
  x,
  knn,
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))
  checkmate::assertClass(knn, "SingleCellNearestNeighbour")

  S7::prop(x, "sc_cache") <- set_knn(
    x = S7::prop(x, "sc_cache"),
    knn = knn
  )

  return(x)
}

#' @name set_snn_graph.SingleCellsGrps
#'
#' @rdname set_snn_graph
#'
#' @method set_snn_graph SingleCellsGrps
S7::method(set_snn_graph, SingleCellsGrps) <- function(
  x,
  snn_graph,
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))
  checkmate::assertClass(snn_graph, "igraph")

  S7::prop(x, "sc_cache") <- set_snn_graph(
    x = S7::prop(x, "sc_cache"),
    snn_graph = snn_graph
  )

  return(x)
}

#' @name remove_knn.SingleCellsGrps
#'
#' @rdname remove_knn
#'
#' @method remove_knn SingleCellsGrps
S7::method(remove_knn, SingleCellsGrps) <- function(
  x,
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))

  S7::prop(x, "sc_cache") <- remove_knn(
    x = S7::prop(x, "sc_cache")
  )

  return(x)
}

#' @name remove_snn_graph.SingleCellsGrps
#'
#' @rdname remove_snn_graph
#'
#' @method remove_snn_graph SingleCellsGrps
S7::method(remove_snn_graph, SingleCellsGrps) <- function(
  x,
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))

  S7::prop(x, "sc_cache") <- remove_snn_graph(
    x = S7::prop(x, "sc_cache")
  )

  return(x)
}

#### getters -------------------------------------------------------------------

#' @name get_pca_factors.SingleCellsGrps
#'
#' @rdname get_pca_factors
#'
#' @method get_pca_factors SingleCellsGrps
S7::method(get_pca_factors, SingleCellsGrps) <- function(
  x,
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))

  res <- get_pca_factors(
    x = S7::prop(x, "sc_cache")
  )

  if (is.null(res)) {
    return(NULL)
  }

  rownames(res) <- S7::prop(x, "obs_table")$meta_cell_id
  colnames(res) <- sprintf("PC_%i", 1:ncol(res))

  return(res)
}

#' @name get_pca_loadings.SingleCellsGrps
#'
#' @rdname get_pca_loadings
#'
#' @method get_pca_loadings SingleCellsGrps
S7::method(get_pca_loadings, SingleCellsGrps) <- function(
  x,
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))

  res <- get_pca_loadings(
    x = S7::prop(x, "sc_cache")
  )

  if (is.null(res)) {
    return(NULL)
  }

  colnames(res) <- sprintf("PC_%i", 1:ncol(res))
  rownames(res) <- S7::prop(x, "var_table")$gene_id[get_hvg(x)]

  return(res)
}

#' @name get_pca_singular_val.SingleCellsGrps
#'
#' @rdname get_pca_singular_val
#'
#' @method get_pca_singular_val SingleCellsGrps
S7::method(get_pca_singular_val, SingleCellsGrps) <- function(
  x,
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))

  res <- get_pca_singular_val(
    x = S7::prop(x, "sc_cache")
  )

  return(res)
}

#' @name get_embedding.SingleCellsGrps
#'
#' @rdname get_embedding
#'
#' @method get_embedding SingleCellsGrps
S7::method(get_embedding, SingleCellsGrps) <- function(
  x,
  embd_name,
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))
  checkmate::qassert(embd_name, "S1")

  res <- get_embedding(
    x = S7::prop(x, "sc_cache"),
    embd_name = embd_name
  )

  rownames(res) <- S7::prop(x, "obs_table")$meta_cell_id

  return(res)
}

#' @name get_available_embeddings.SingleCellsGrps
#'
#' @rdname get_available_embeddings
#'
#' @method get_available_embeddings SingleCellsGrps
S7::method(get_available_embeddings, SingleCellsGrps) <- function(
  x,
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))

  res <- get_available_embeddings(
    x = S7::prop(x, "sc_cache")
  )

  return(res)
}

#' @name get_knn_mat.SingleCellsGrps
#'
#' @rdname get_knn_mat
#'
#' @method get_knn_mat SingleCellsGrps
S7::method(get_knn_mat, SingleCellsGrps) <- function(
  x,
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))

  res <- get_knn_mat(
    x = S7::prop(x, "sc_cache")
  )

  return(res)
}

#' @name get_knn_dist.SingleCellsGrps
#'
#' @rdname get_knn_dist
#'
#' @method get_knn_dist SingleCellsGrps
S7::method(get_knn_dist, SingleCellsGrps) <- function(
  x,
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))

  res <- get_knn_dist(
    x = S7::prop(x, "sc_cache")
  )

  return(res)
}

#' @name get_knn_obj.SingleCellsGrps
#'
#' @rdname get_knn_obj
#'
#' @method get_knn_obj SingleCellsGrps
S7::method(get_knn_obj, SingleCellsGrps) <- function(
  x,
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))

  res <- get_knn_obj(
    x = S7::prop(x, "sc_cache")
  )

  return(res)
}

#' @name get_snn_graph.SingleCellsGrps
#'
#' @rdname get_snn_graph
#'
#' @method get_snn_graph SingleCellsGrps
S7::method(get_snn_graph, SingleCellsGrps) <- function(
  x,
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsGrps))

  res <- get_snn_graph(
    x = S7::prop(x, "sc_cache")
  )

  return(res)
}
