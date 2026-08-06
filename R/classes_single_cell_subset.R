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
#' @param sc_object A [bixverse::SingleCells()] object to subset.
#' @param grouping_column String. Column in the obs table that defines the
#' grouping.
#' @param group String. Level of `grouping_column` to retain.
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

#' @noRd
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

#' @noRd
S7::method(dim, SingleCellsSubset) <- function(x) {
  c(nrow(S7::prop(x, "obs_table")), nrow(S7::prop(x, "var_table")))
}

#' @noRd
S7::method(head, SingleCellsSubset) <- function(x, ..., n = 6L) {
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

#' @name get_gene_names.SingleCellsSubset
#'
#' @rdname get_gene_names
#'
#' @method get_gene_names SingleCellsSubset
S7::method(get_gene_names, SingleCellsSubset) <- function(x) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsSubset))
  # gene space is unchanged in the subset
  S7::prop(x, "var_table")$gene_id
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

# Registered as a plain S3 method instead of via `S7::method()`: S7 registers
# base-generic methods dynamically with the function object, which makes
# `R CMD check`'s replacement-function check choke on any `<-` generic.
#' @exportS3Method "[[<-" "bixverse::SingleCellsSubset"
#'
#' @noRd
`[[<-.bixverse::SingleCellsSubset` <- function(x, i, ..., value) {
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

## write-back ------------------------------------------------------------------

#' Merge obs columns from subsets back into the parent object
#'
#' @description
#' Takes columns computed on one or more [bixverse::SingleCellsSubset()]
#' objects and writes them into the parent's obs table in the DuckDB, joining
#' on `cell_idx`. The join is a left join, so every parent cell that is not
#' part of any of the provided subsets ends up as `NA`.
#'
#' Pass a single subset or a list of them, e.g. the output of
#' [apply_pipeline_per_group()]. All subsets are written in a single join, and
#' the function refuses to run if two subsets claim the same cell.
#'
#' @param object `SingleCells` class. The parent the subsets were derived from.
#' @param subsets A `SingleCellsSubset` or a list of them. Every element must
#' originate from `object`; this is verified against the cell names.
#' @param cols Optional character vector. Columns in the subset obs tables to
#' merge. If `NULL`, every column that is present in all subsets but absent
#' from the parent obs table is taken, which after a pipeline run is exactly
#' the set of newly generated columns.
#' @param new_names Optional character vector. Names to give the merged columns
#' in the parent obs table. Same length as `cols`, which must be given
#' explicitly if this is used.
#' @param prefix_values Boolean. Prefix every merged value with the subset's
#' group, i.e. `"<group>_<value>"`. Coerces the columns to character, so this
#' only makes sense for discrete labels. Needed when merging a list of subsets
#' into a shared column, since sub-cluster `1` of one group and sub-cluster `1`
#' of another are otherwise indistinguishable.
#' @param overwrite Boolean. Allow merged columns to replace existing columns
#' of the same name in the parent obs table. Defaults to `FALSE`, in which case
#' a name clash is an error.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @returns The parent `SingleCells` object with the merged columns added to
#' the obs table in the DuckDB.
#'
#' @seealso [add_sc_new_obs()] for the equivalent on result objects that carry
#' their own `cell_idx`, such as [fast_cluster_sc()].
#'
#' @export
merge_subset_obs <- function(
  object,
  subsets,
  cols = NULL,
  new_names = NULL,
  prefix_values = FALSE,
  overwrite = FALSE,
  .verbose = TRUE
) {
  # checks
  checkmate::assertClass(object, "bixverse::SingleCells")
  if (S7::S7_inherits(subsets, SingleCellsSubset)) {
    subsets <- list(subsets)
  }
  checkmate::assertList(subsets, min.len = 1L)
  purrr::walk(
    subsets,
    \(x) checkmate::assertClass(x, "bixverse::SingleCellsSubset")
  )
  checkmate::qassert(cols, c("0", "S+"))
  checkmate::qassert(new_names, c("0", "S+"))
  checkmate::qassert(prefix_values, "B1")
  checkmate::qassert(overwrite, "B1")
  checkmate::qassert(.verbose, "B1")

  parent_names <- get_cell_names(object, filtered = FALSE)
  parent_cols <- get_sc_duckdb(object)$get_obs_cols()

  # cell_idx is a 1-indexed position into the parent's full obs table. Verify
  # that assumption holds before writing anything into the DuckDB.
  purrr::iwalk(subsets, \(sub, i) {
    idx <- S7::prop(sub, "subset_to_original")
    if (max(idx) > length(parent_names) || min(idx) < 1L) {
      stop(sprintf("Subset %s has cell indices outside of the parent.", i))
    }
    if (!identical(get_cell_names(sub), parent_names[idx])) {
      stop(sprintf(
        "Subset %s does not originate from the provided parent object.",
        i
      ))
    }
  })

  shared_cols <- Reduce(
    intersect,
    purrr::map(subsets, \(x) names(S7::prop(x, "obs_table")))
  )

  if (is.null(cols)) {
    if (!is.null(new_names)) {
      stop("`cols` needs to be given explicitly when using `new_names`.")
    }
    cols <- setdiff(shared_cols, c(parent_cols, "cell_idx"))
    if (length(cols) == 0) {
      stop("No new obs columns found in the subset(s). Nothing to merge.")
    }
  } else {
    if ("cell_idx" %in% cols) {
      stop("`cell_idx` is the join key and cannot be merged.")
    }
    missing_cols <- setdiff(cols, shared_cols)
    if (length(missing_cols) > 0) {
      stop(sprintf(
        "Column(s) %s missing from at least one subset obs table.",
        toString(sprintf("`%s`", missing_cols))
      ))
    }
  }

  if (is.null(new_names)) {
    new_names <- cols
  }
  checkmate::assertTRUE(length(new_names) == length(cols))

  clashes <- intersect(new_names, parent_cols)
  if (length(clashes) > 0 && !overwrite) {
    stop(sprintf(
      "Column(s) %s already exist in the parent obs. Set overwrite = TRUE.",
      toString(sprintf("`%s`", clashes))
    ))
  }

  to_merge <- purrr::map(subsets, \(sub) {
    dt <- get_sc_obs(sub, cols = c("cell_idx", cols))
    if (prefix_values) {
      grp <- S7::prop(sub, "group")
      dt[,
        (cols) := lapply(.SD, \(v) sprintf("%s_%s", grp, v)),
        .SDcols = cols
      ]
    }
    data.table::setnames(dt, cols, new_names)
    dt
  })

  if (length(to_merge) > 1 && !prefix_values) {
    purrr::walk(new_names, \(col) {
      seen <- unlist(purrr::map(
        to_merge,
        \(dt) unique(as.character(dt[[col]]))
      ))
      shared <- unique(seen[duplicated(seen)])
      shared <- shared[!is.na(shared)]
      if (length(shared) > 0) {
        warning(sprintf(
          paste(
            "Value(s) %s of column `%s` appear in more than one subset.",
            "Set prefix_values = TRUE to disambiguate them."
          ),
          toString(sprintf(
            "`%s`",
            shared[seq_len(min(5L, length(shared)))]
          )),
          col
        ))
      }
    })
  }

  merged <- data.table::rbindlist(to_merge, use.names = TRUE)

  if (anyDuplicated(merged$cell_idx) > 0) {
    stop("Duplicated cells across the subsets. Cells can be merged only once.")
  }

  merged <- .add_is_obs_attr(merged)

  get_sc_duckdb(object)$join_data_obs(merged)

  if (.verbose) {
    message(sprintf(
      "Merged %i column(s) for %i cells into obs (%i cells left as NA).",
      length(new_names),
      nrow(merged),
      length(parent_names) - nrow(merged)
    ))
  }

  return(object)
}
