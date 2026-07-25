# single cell subset class -----------------------------------------------------

## s7 object -------------------------------------------------------------------

#' @title bixverse single cell subset class
#'
#' @description
#' Subset view onto a [bixverse::SingleCells()] object, restricted to cells
#' belonging to a single level of a grouping variable. The Rust count
#' connection is shared with the parent (no data copy). `obs_table` and
#' `var_table` are held in memory; `sc_map` is rebuilt to point only at the
#' subset cells but stays in the original index space so Rust calls remain
#' valid without further translation.
#'
#' @section Properties:
#' \describe{
#'   \item{count_connection}{Shared Rust pointer to the on-disk counts.}
#'   \item{dir_data}{Directory holding the binary count files.}
#'   \item{obs_table}{Subset obs (rows for the chosen group only).
#'   `cell_idx` keeps the original 1-indexed position in the parent.}
#'   \item{var_table}{Variable/feature table (unchanged from parent).}
#'   \item{grouping_column}{Column in obs used to define the subset.}
#'   \item{group}{Value of `grouping_column` represented by this subset.}
#'   \item{sc_cache}{Fresh `ScCache` for subset-specific PCA, kNN, sNN,
#'   embeddings.}
#'   \item{sc_map}{`ScMap` restricted to the subset cells. `cell_mapping`
#'   stays 1-indexed and `cells_to_keep_idx` stays 0-indexed, both in the
#'   original parent index space.}
#'   \item{subset_to_original}{Integer vector. 1-indexed original cell
#'   positions, in subset row order. `subset_to_original[i]` is the parent
#'   position of subset row `i`.}
#'   \item{dims}{`c(n_cells_subset, n_genes)`.}
#' }
#'
#' @return A `SingleCellsSubset` object.
#'
#' @export
SingleCellsSubset <- S7::new_class(
  name = "SingleCellsSubset",
  properties = list(
    count_connection = S7::class_any,
    dir_data = S7::class_character,
    obs_table = S7::class_data.frame,
    var_table = S7::class_data.frame,
    grouping_column = S7::class_character,
    group = S7::class_character,
    sc_cache = S7::class_any,
    sc_map = S7::class_any,
    subset_to_original = S7::class_integer,
    dims = S7::class_integer
  ),
  constructor = function(
    sc_object,
    grouping_column,
    group
  ) {
    checkmate::assertClass(sc_object, "bixverse::SingleCells")
    checkmate::assertString(grouping_column)
    checkmate::assertString(group)

    dir_data <- S7::prop(sc_object, "dir_data")
    sc_map <- S7::prop(sc_object, "sc_map")
    count_connection <- sc_object@count_connection

    obs_table <- get_sc_obs(sc_object, filtered = TRUE)
    checkmate::assertNames(
      x = colnames(obs_table),
      must.include = grouping_column
    )
    obs_table <- obs_table[get(grouping_column) == group, ]

    if (nrow(obs_table) == 0L) {
      stop(sprintf(
        "No cells found for group '%s' in column '%s'.",
        group,
        grouping_column
      ))
    }

    var_table <- get_sc_var(sc_object)

    # 1-indexed original positions per subset row. Mirrors obs_table$cell_idx
    # but exposed as a first-class property so translation does not depend on
    # an obs column name.
    subset_to_original <- as.integer(obs_table$cell_idx)

    if (!is.null(sc_map$cells_to_keep_idx)) {
      # cells_to_keep_idx is 0-indexed; compare to 0-indexed subset positions
      sc_map$cells_to_keep_idx <- sc_map$cells_to_keep_idx[
        sc_map$cells_to_keep_idx %in% (subset_to_original - 1L)
      ]
    }
    # HVG is subset-specific
    sc_map$hvg_gene_indices <- NULL

    S7::new_object(
      S7::S7_object(),
      count_connection = count_connection,
      dir_data = dir_data,
      obs_table = obs_table,
      var_table = var_table,
      grouping_column = grouping_column,
      group = group,
      sc_cache = new_sc_cache(),
      sc_map = sc_map,
      subset_to_original = subset_to_original,
      dims = c(nrow(obs_table), nrow(var_table))
    )
  }
)

## primitives ------------------------------------------------------------------

#' @name print.SingleCellsSubset
#'
#' @title print Method for SingleCellsSubset object
#'
#' @param x An object of class `SingleCellsSubset`.
#' @param ... Additional arguments (currently not used).
#'
#' @returns Invisibly returns `x`.
#'
#' @method print SingleCellsSubset
#'
#' @keywords internal
S7::method(print, SingleCellsSubset) <- function(x, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))

  dims <- S7::prop(x, "dims")
  group <- S7::prop(x, "group")
  grouping_column <- S7::prop(x, "grouping_column")
  sc_cache <- S7::prop(x, "sc_cache")
  sc_map <- S7::prop(x, "sc_map")
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
    "Single cell experiment (subset).\n",
    sprintf("  No cells: %i\n", dims[1]),
    sprintf("  No genes: %i\n", dims[2]),
    sprintf("  Group: %s = %s\n", grouping_column, group),
    sprintf("  HVG calculated: %s\n", hvg_calculated),
    sprintf("  PCA calculated: %s\n", pca_calculated),
    sprintf("  Other embeddings: %s\n", other_embeddings_str),
    sprintf("  KNN generated: %s\n", knn_generated),
    sprintf("  SNN generated: %s\n", snn_generated),
    sep = ""
  )
  invisible(x)
}

#' @name dim.SingleCellsSubset
#'
#' @title dim Method for SingleCellsSubset object
#'
#' @param x An object of class `SingleCellsSubset`.
#'
#' @returns An integer vector of length 2 with the number of cells and genes.
#'
#' @method dim SingleCellsSubset
#'
#' @keywords internal
S7::method(dim, SingleCellsSubset) <- function(x) {
  c(nrow(S7::prop(x, "obs_table")), nrow(S7::prop(x, "var_table")))
}

#' @name head.SingleCellsSubset
#'
#' @title head Method for SingleCellsSubset object
#'
#' @param x An object of class `SingleCellsSubset`.
#' @param n Integer. Number of rows to return. Defaults to `6L`.
#' @param ... Additional arguments (currently not used).
#'
#' @returns A data.table with the first `n` rows of the obs table.
#'
#' @method head SingleCellsSubset
#'
#' @keywords internal
S7::method(head, SingleCellsSubset) <- function(x, n = 6L, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  checkmate::qassert(n, "I1[1,)")
  n <- min(n, nrow(S7::prop(x, "obs_table")))
  get_sc_obs(x, indices = seq_len(n))
}

## getters ---------------------------------------------------------------------

### obs ------------------------------------------------------------------------

#' @method get_sc_obs SingleCellsSubset
#'
#' @export
S7::method(get_sc_obs, SingleCellsSubset) <- function(
  object,
  indices = NULL,
  cols = NULL,
  filtered = FALSE
) {
  checkmate::assertClass(object, "bixverse::SingleCellsSubset")
  checkmate::qassert(indices, c("0", "I+"))
  checkmate::qassert(filtered, "B1")

  # `filtered` kept for parent API parity; the subset is always filtered.
  obs_table <- data.table::copy(S7::prop(object, "obs_table"))
  if (!is.null(indices)) {
    obs_table <- obs_table[indices, ]
  }
  if (!is.null(cols)) {
    obs_table <- obs_table[, cols, with = FALSE]
  }
  obs_table
}

#' @method `[[` SingleCellsSubset
#'
#' @export
S7::method(`[[`, SingleCellsSubset) <- function(x, i, ...) {
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

### vars -----------------------------------------------------------------------

#' @method get_sc_var SingleCellsSubset
#'
#' @export
S7::method(get_sc_var, SingleCellsSubset) <- function(
  object,
  indices = NULL,
  cols = NULL,
  modality = c("rna", "adt")
) {
  modality <- match.arg(modality)

  checkmate::assertClass(object, "bixverse::SingleCellsSubset")
  checkmate::qassert(indices, c("0", "I+"))

  if (modality != "rna") {
    stop(
      "SingleCellsSubset only supports modality = 'rna'."
    )
  }

  var_table <- object@var_table
  if (!is.null(indices)) {
    var_table <- var_table[indices, ]
  }
  if (!is.null(cols)) {
    var_table <- var_table[, cols, with = FALSE]
  }
  var_table
}

### counts ---------------------------------------------------------------------

#' @method get_sc_counts SingleCellsSubset
#'
#' @export
S7::method(get_sc_counts, SingleCellsSubset) <- function(
  object,
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
    stop("SingleCellsSubset only supports modality = 'rna'.")
  }

  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsSubset))
  checkmate::assertChoice(assay, c("raw", "norm"))
  checkmate::assertChoice(return_format, c("cell", "gene"))
  checkmate::qassert(cell_indices, c("0", "I+"))
  checkmate::qassert(gene_indices, c("0", "I+"))
  checkmate::qassert(.verbose, "B1")

  requireNamespace("Matrix", quietly = TRUE)

  rust_con <- get_sc_rust_ptr(object)
  sc_map <- get_sc_map(object)
  subset_to_original <- S7::prop(object, "subset_to_original")

  # `cell_indices` arrives as 1-indexed SUBSET positions. Translate to
  # 1-indexed ORIGINAL positions before any Rust interaction.
  original_cell_indices <- if (!is.null(cell_indices)) {
    if (any(cell_indices < 1L | cell_indices > length(subset_to_original))) {
      stop("cell_indices out of bounds for the subset.")
    }
    as.integer(subset_to_original[cell_indices])
  } else {
    NULL
  }

  cells_to_keep <- get_cells_to_keep(object) # 0-indexed original space

  if (use_cells_to_keep) {
    if (!is.null(original_cell_indices) && length(cells_to_keep) > 0L) {
      original_cell_indices <- as.integer(
        intersect(original_cell_indices, cells_to_keep + 1L)
      )
    } else if (length(cells_to_keep) > 0L) {
      original_cell_indices <- as.integer(cells_to_keep + 1L)
    }
  }

  count_data <- get_counts_from_rust(
    rust_con = rust_con,
    assay = assay,
    return_format = return_format,
    cell_indices = original_cell_indices,
    gene_indices = gene_indices,
    .verbose = .verbose
  )

  count_data <- create_sparse_matrix(
    count_data = count_data,
    return_format = return_format
  )

  count_data <- finalise_matrix(
    matrix = count_data,
    return_format = return_format,
    cell_indices = original_cell_indices,
    gene_indices = gene_indices,
    sc_map = sc_map
  )

  count_data
}

### map getter -----------------------------------------------------------------

#' @method get_sc_map SingleCellsSubset
#'
#' @export
S7::method(get_sc_map, SingleCellsSubset) <- function(object) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsSubset))
  S7::prop(object, "sc_map")
}

#' @name get_cells_to_keep.SingleCellsSubset
#'
#' @rdname get_cells_to_keep
#'
#' @method get_cells_to_keep SingleCellsSubset
S7::method(get_cells_to_keep, SingleCellsSubset) <- function(x) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  res <- get_cells_to_keep(x = S7::prop(x, "sc_map"))
  if (length(res) == 0) {
    res <- S7::prop(x, "subset_to_original") - 1L
  }
  as.integer(res)
}

#' @name get_cell_names.SingleCellsSubset
#'
#' @rdname get_cell_names
#'
#' @method get_cell_names SingleCellsSubset
S7::method(get_cell_names, SingleCellsSubset) <- function(x, filtered = FALSE) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  # subset obs is already filtered to the group; `filtered` accepted for parity
  S7::prop(x, "obs_table")$cell_id
}

#' @name get_gene_names_from_idx.SingleCellsSubset
#'
#' @rdname get_gene_names_from_idx
#'
#' @method get_gene_names_from_idx SingleCellsSubset
S7::method(get_gene_names_from_idx, SingleCellsSubset) <- function(
  x,
  gene_idx,
  rust_based = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  idx <- if (rust_based) gene_idx + 1L else gene_idx
  S7::prop(x, "var_table")$gene_id[idx]
}

#### env getters ---------------------------------------------------------------

#' @method get_sc_rust_ptr SingleCellsSubset
#'
#' @export
S7::method(get_sc_rust_ptr, SingleCellsSubset) <- function(object) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsSubset))
  S7::prop(object, "count_connection")
}

#' @method get_rust_count_gene_f_path SingleCellsSubset
#'
#' @export
S7::method(get_rust_count_gene_f_path, SingleCellsSubset) <- function(object) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsSubset))
  file.path(S7::prop(object, "dir_data"), "counts_genes.bin")
}

#' @method get_rust_count_cell_f_path SingleCellsSubset
#'
#' @export
S7::method(get_rust_count_cell_f_path, SingleCellsSubset) <- function(object) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsSubset))
  file.path(S7::prop(object, "dir_data"), "counts_cells.bin")
}

#' @method `[` SingleCellsSubset
#'
#' @export
S7::method(`[`, SingleCellsSubset) <- function(
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
    gene_indices = j,
    ...
  )
}

## setters ---------------------------------------------------------------------

### obs table ------------------------------------------------------------------

#' @method `[[<-` SingleCellsSubset
#'
#' @export
S7::method(`[[<-`, SingleCellsSubset) <- function(x, i, ..., value) {
  checkmate::assertClass(x, "bixverse::SingleCellsSubset")
  checkmate::qassert(i, "S+")

  if (length(i) == 1) {
    checkmate::qassert(value, "a")
    S7::prop(x, "obs_table")[, (i) := value]
  } else {
    checkmate::assertList(value, names = "named", types = "atomic")
    S7::prop(x, "obs_table")[, (i) := value]
  }
  x
}

### hvg ------------------------------------------------------------------------

#' @name set_hvg.SingleCellsSubset
#'
#' @rdname set_hvg
#'
#' @method set_hvg SingleCellsSubset
S7::method(set_hvg, SingleCellsSubset) <- function(x, hvg) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  checkmate::qassert(hvg, c("I+", "S+"))

  var_table <- S7::prop(x, "var_table")

  # Resolve to 1-indexed first
  idx_1 <- if (is.numeric(hvg)) {
    as.integer(hvg)
  } else {
    m <- match(hvg, var_table$gene_id)
    missing <- is.na(m)
    if (any(missing)) {
      warning(sprintf(
        "With the provided hvg gene_ids a total of %i could not be matched.",
        sum(missing)
      ))
      m <- m[!missing]
    }
    if (length(m) == 0L) {
      stop(
        "The HVG indices have length 0. Please double check provided parameters!"
      )
    }
    as.integer(m)
  }

  # Store 0-indexed (parent SingleCells convention)
  sc_map <- S7::prop(x, "sc_map")
  sc_map$hvg_gene_indices <- idx_1 - 1L
  S7::prop(x, "sc_map") <- sc_map
  x
}

#' @name get_hvg.SingleCellsSubset
#'
#' @rdname get_hvg
#'
#' @method get_hvg SingleCellsSubset
S7::method(get_hvg, SingleCellsSubset) <- function(x) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  hvg <- as.integer(S7::prop(x, "sc_map")[["hvg_gene_indices"]])
  if (length(hvg) == 0L) {
    warning("No highly variable features found in the class.")
  }
  hvg
}

### cell / gene index lookup ---------------------------------------------------

#' @name get_cell_indices.SingleCellsSubset
#'
#' @rdname get_cell_indices
#'
#' @method get_cell_indices SingleCellsSubset
S7::method(get_cell_indices, SingleCellsSubset) <- function(
  x,
  cell_ids,
  rust_index
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  checkmate::qassert(cell_ids, "S+")
  checkmate::qassert(rust_index, "B1")

  obs_table <- S7::prop(x, "obs_table")
  subset_pos <- match(cell_ids, obs_table[["cell_id"]])

  missing <- is.na(subset_pos)
  if (any(missing)) {
    warning(sprintf(
      "With the provided cell_ids a total of %i could not be matched.",
      sum(missing)
    ))
    subset_pos <- subset_pos[!missing]
  }
  if (length(subset_pos) == 0L) {
    stop("The cell indices have length 0.")
  }

  if (rust_index) {
    # Translate subset position -> 0-indexed ORIGINAL position
    subset_to_original <- S7::prop(x, "subset_to_original")
    as.integer(subset_to_original[subset_pos] - 1L)
  } else {
    # Return 1-indexed SUBSET position (what `[` and `[[` expect)
    as.integer(subset_pos)
  }
}

#' @name get_gene_indices.SingleCellsSubset
#'
#' @rdname get_gene_indices
#'
#' @method get_gene_indices SingleCellsSubset
S7::method(get_gene_indices, SingleCellsSubset) <- function(
  x,
  gene_ids,
  rust_index
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  checkmate::qassert(gene_ids, "S+")
  checkmate::qassert(rust_index, "B1")

  # Gene space is unchanged in the subset, so no translation needed.
  indices <- match(gene_ids, x@var_table[["gene_id"]])

  missing <- is.na(indices)
  if (any(missing)) {
    warning(sprintf(
      "With the provided gene_ids a total of %i could not be matched.",
      sum(missing)
    ))
    indices <- indices[!missing]
  }
  if (length(indices) == 0L) {
    stop(
      "The gene indices have length 0. Please double check provided parameters!"
    )
  }

  if (rust_index) as.integer(indices - 1L) else as.integer(indices)
}

### sc cache -------------------------------------------------------------------

#### setters -------------------------------------------------------------------

#' @name set_pca_factors.SingleCellsSubset
#'
#' @rdname set_pca_factors
#'
#' @method set_pca_factors SingleCellsSubset
S7::method(set_pca_factors, SingleCellsSubset) <- function(
  x,
  pca_factor,
  ...
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  checkmate::assertMatrix(pca_factor, mode = "numeric")
  S7::prop(x, "sc_cache") <- set_pca_factors(
    x = S7::prop(x, "sc_cache"),
    pca_factor = pca_factor
  )
  x
}

#' @name set_pca_loadings.SingleCellsSubset
#'
#' @rdname set_pca_loadings
#'
#' @method set_pca_loadings SingleCellsSubset
S7::method(set_pca_loadings, SingleCellsSubset) <- function(
  x,
  pca_loading,
  ...
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  checkmate::assertMatrix(pca_loading, mode = "numeric")
  S7::prop(x, "sc_cache") <- set_pca_loadings(
    x = S7::prop(x, "sc_cache"),
    pca_loading = pca_loading
  )
  x
}

#' @name set_pca_singular_vals.SingleCellsSubset
#'
#' @rdname set_pca_singular_vals
#'
#' @method set_pca_singular_vals SingleCellsSubset
S7::method(set_pca_singular_vals, SingleCellsSubset) <- function(
  x,
  singular_vals,
  ...
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  checkmate::qassert(singular_vals, "N+")
  S7::prop(x, "sc_cache") <- set_pca_singular_vals(
    x = S7::prop(x, "sc_cache"),
    singular_vals = singular_vals
  )
  x
}

#' @name set_embedding.SingleCellsSubset
#'
#' @rdname set_embedding
#'
#' @method set_embedding SingleCellsSubset
S7::method(set_embedding, SingleCellsSubset) <- function(
  x,
  embd,
  name,
  ...
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  checkmate::assertMatrix(embd, mode = "numeric")
  checkmate::qassert(name, "S1")
  S7::prop(x, "sc_cache") <- set_embedding(
    x = S7::prop(x, "sc_cache"),
    embd = embd,
    name = name
  )
  x
}

#' @name set_knn.SingleCellsSubset
#'
#' @rdname set_knn
#'
#' @method set_knn SingleCellsSubset
S7::method(set_knn, SingleCellsSubset) <- function(x, knn, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  checkmate::assertClass(knn, "SingleCellNearestNeighbour")
  S7::prop(x, "sc_cache") <- set_knn(
    x = S7::prop(x, "sc_cache"),
    knn = knn
  )
  x
}

#' @name set_snn_graph.SingleCellsSubset
#'
#' @rdname set_snn_graph
#'
#' @method set_snn_graph SingleCellsSubset
S7::method(set_snn_graph, SingleCellsSubset) <- function(x, snn_graph, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  checkmate::assertClass(snn_graph, "igraph")
  S7::prop(x, "sc_cache") <- set_snn_graph(
    x = S7::prop(x, "sc_cache"),
    snn_graph = snn_graph
  )
  x
}

#' @name remove_knn.SingleCellsSubset
#'
#' @rdname remove_knn
#'
#' @method remove_knn SingleCellsSubset
S7::method(remove_knn, SingleCellsSubset) <- function(x, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  S7::prop(x, "sc_cache") <- remove_knn(x = S7::prop(x, "sc_cache"))
  x
}

#' @name remove_snn_graph.SingleCellsSubset
#'
#' @rdname remove_snn_graph
#'
#' @method remove_snn_graph SingleCellsSubset
S7::method(remove_snn_graph, SingleCellsSubset) <- function(x, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  S7::prop(x, "sc_cache") <- remove_snn_graph(x = S7::prop(x, "sc_cache"))
  x
}

#### getters -------------------------------------------------------------------

#' @name get_pca_factors.SingleCellsSubset
#'
#' @rdname get_pca_factors
#'
#' @method get_pca_factors SingleCellsSubset
S7::method(get_pca_factors, SingleCellsSubset) <- function(x, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  res <- get_pca_factors(x = S7::prop(x, "sc_cache"))
  if (is.null(res)) {
    return(NULL)
  }
  # TODO: confirm cell id column name
  rownames(res) <- S7::prop(x, "obs_table")$cell_id
  colnames(res) <- sprintf("PC_%i", seq_len(ncol(res)))
  res
}

#' @name get_pca_loadings.SingleCellsSubset
#'
#' @rdname get_pca_loadings
#'
#' @method get_pca_loadings SingleCellsSubset
S7::method(get_pca_loadings, SingleCellsSubset) <- function(x, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  res <- get_pca_loadings(x = S7::prop(x, "sc_cache"))
  if (is.null(res)) {
    return(NULL)
  }
  colnames(res) <- sprintf("PC_%i", seq_len(ncol(res)))
  # get_hvg returns 0-indexed; +1 for R subsetting
  hvg <- get_hvg(x) + 1L
  rownames(res) <- S7::prop(x, "var_table")$gene_id[hvg]
  res
}

#' @name get_pca_singular_val.SingleCellsSubset
#'
#' @rdname get_pca_singular_val
#'
#' @method get_pca_singular_val SingleCellsSubset
S7::method(get_pca_singular_val, SingleCellsSubset) <- function(x, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  get_pca_singular_val(x = S7::prop(x, "sc_cache"))
}

#' @name get_embedding.SingleCellsSubset
#'
#' @rdname get_embedding
#'
#' @method get_embedding SingleCellsSubset
S7::method(get_embedding, SingleCellsSubset) <- function(x, embd_name, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  checkmate::qassert(embd_name, "S1")
  res <- get_embedding(x = S7::prop(x, "sc_cache"), embd_name = embd_name)
  if (is.null(res)) {
    return(NULL)
  }
  # TODO: confirm cell id column name
  rownames(res) <- S7::prop(x, "obs_table")$cell_id
  res
}

#' @name get_available_embeddings.SingleCellsSubset
#'
#' @rdname get_available_embeddings
#'
#' @method get_available_embeddings SingleCellsSubset
S7::method(get_available_embeddings, SingleCellsSubset) <- function(x, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  get_available_embeddings(x = S7::prop(x, "sc_cache"))
}

#' @name get_knn_mat.SingleCellsSubset
#'
#' @rdname get_knn_mat
#'
#' @method get_knn_mat SingleCellsSubset
S7::method(get_knn_mat, SingleCellsSubset) <- function(x, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  get_knn_mat(x = S7::prop(x, "sc_cache"))
}

#' @name get_knn_dist.SingleCellsSubset
#'
#' @rdname get_knn_dist
#'
#' @method get_knn_dist SingleCellsSubset
S7::method(get_knn_dist, SingleCellsSubset) <- function(x, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  get_knn_dist(x = S7::prop(x, "sc_cache"))
}

#' @name get_knn_obj.SingleCellsSubset
#'
#' @rdname get_knn_obj
#'
#' @method get_knn_obj SingleCellsSubset
S7::method(get_knn_obj, SingleCellsSubset) <- function(x, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  get_knn_obj(x = S7::prop(x, "sc_cache"))
}

#' @name get_snn_graph.SingleCellsSubset
#'
#' @rdname get_snn_graph
#'
#' @method get_snn_graph SingleCellsSubset
S7::method(get_snn_graph, SingleCellsSubset) <- function(x, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  get_snn_graph(x = S7::prop(x, "sc_cache"))
}
