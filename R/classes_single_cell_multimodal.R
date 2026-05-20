# s3 ---------------------------------------------------------------------------

## adt counts ------------------------------------------------------------------

#' Generates a new `ADTCounts` class
#'
#' @description
#' This function generates a new `ADTCounts` class which uses CLR normalisation
#' under the hood.
#'
#' @param raw_counts Numeric matrix. The raw ADT counts.
#' @param cell_info Named integer vector. Output of [get_cell_info()]. Defines
#' as elements the cell indices (R-based) and as names the barcodes.
#' @param clean_clr_counts Boolean. Shall the per-protein 1st percentile be
#' removed from the CLR normalised counts.
#' @param percentile Numeric. The percentile to remove to reduce background
#' effects.
#'
#' @returns `ADTCounts` that contains the raw and normalised ADT counts.
#'
#' @export
new_adt_counts_clr <- function(
  raw_counts,
  cell_info,
  clean_clr_counts = TRUE,
  percentile = 0.01
) {
  # checks
  checkmate::assertMatrix(
    raw_counts,
    mode = "numeric",
    row.names = "named",
    col.names = "named"
  )
  checkmate::assertInteger(cell_info, names = "named")
  checkmate::qassert(clean_clr_counts, "B1")
  checkmate::qassert(percentile, "N1(0,1)")

  raw_counts <- raw_counts[names(cell_info), ]
  norm_counts <- rs_adt_clr(counts = raw_counts)
  colnames(norm_counts) <- colnames(raw_counts)
  rownames(norm_counts) <- rownames(raw_counts)

  if (clean_clr_counts) {
    for (j in seq_len(ncol(norm_counts))) {
      norm_counts[, j] <- norm_counts[, j] -
        quantile(norm_counts[, j], percentile)
    }
  }

  structure(
    list(
      raw_counts = raw_counts,
      norm_counts = norm_counts,
      cell_info = cell_info,
      other_data = list(),
      type = "CLR"
    ),
    class = "ADTCounts"
  )
}

#' Generates a new `ADTCounts` class via DSB normalisation
#'
#' @description
#' This function generates a new `ADTCounts` class using DSB normalisation from
#' Mulè et al., instead of CLR. When `empty_drops` is provided, the per-protein
#' background is estimated from empty droplets. Without `empty_drops`, the
#' function falls back to a 2-component k-means on the log-transformed cell
#' counts ("ModelNegativeADTnorm").
#'
#' @param raw_counts Numeric matrix. Cells x proteins matrix of raw ADT counts
#' with cell barcodes as row names and protein names as column names.
#' @param cell_info Named integer vector. Output of [get_cell_info()]. Defines
#' as elements the cell indices (R-based) and as names the barcodes.
#' @param empty_drops Optional numeric matrix. Cells x proteins matrix of
#' empty-droplet ADT counts. If provided, used to estimate per-protein ambient
#' background.
#' @param isotype_names Optional string vector. Column names in `raw_counts`
#' identifying isotype control proteins. Required when
#' `dsb_params$use_isotype_controls = TRUE`.
#' @param dsb_params List. Output of [params_sc_dsb()] with DSB parameters.
#' @param scale_factor String. One of `c("standardise", "mean_subtract")`. Only
#' used when `empty_drops` is provided.
#' @param seed Integer. Random seed for k-means initialisation.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns `ADTCounts` containing the raw counts and the DSB-normalised
#' counts.
#'
#' @export
#'
#' @references Mulè et al., Nat Commun, 2022
new_adt_counts_dsb <- function(
  raw_counts,
  cell_info,
  empty_drops = NULL,
  isotype_names = NULL,
  dsb_params = params_sc_dsb(),
  scale_factor = c("standardise", "mean_subtract"),
  seed = 42L,
  .verbose = TRUE
) {
  scale_factor <- match.arg(scale_factor)

  # checks
  checkmate::assertMatrix(
    raw_counts,
    mode = "numeric",
    row.names = "named",
    col.names = "named"
  )
  checkmate::assertInteger(cell_info, names = "named")
  if (!is.null(empty_drops)) {
    checkmate::assertMatrix(
      empty_drops,
      mode = "numeric",
      col.names = "named"
    )
    checkmate::assertTRUE(ncol(empty_drops) == ncol(raw_counts))
    checkmate::assertTRUE(
      all(colnames(empty_drops) == colnames(raw_counts))
    )
  }
  checkmate::qassert(isotype_names, c("S+", "0"))
  if (!is.null(isotype_names)) {
    checkmate::assertSubset(isotype_names, colnames(raw_counts))
  }
  if (isTRUE(dsb_params$use_isotype_controls) && is.null(isotype_names)) {
    stop(
      "dsb_params$use_isotype_controls = TRUE but isotype_names is NULL."
    )
  }
  assertScDsbParams(dsb_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  raw_counts <- raw_counts[names(cell_info), ]

  # resolve isotype names to 0-based column indices
  isotype_indices <- if (is.null(isotype_names)) {
    integer(0)
  } else {
    as.integer(match(isotype_names, colnames(raw_counts)) - 1L)
  }

  dsb_res <- rs_dsb(
    raw_counts = raw_counts,
    background_counts = empty_drops,
    isotype_indices = isotype_indices,
    dsb_params = dsb_params,
    scale_factor = scale_factor,
    seed = seed,
    verbose = parse_verbosity(.verbose)
  )

  norm_counts <- dsb_res$norm_counts
  colnames(norm_counts) <- colnames(raw_counts)
  rownames(norm_counts) <- rownames(raw_counts)

  other_data <- list(
    protein_background_mean = dsb_res$protein_background_mean,
    protein_background_sd = dsb_res$protein_background_sd,
    technical_component = dsb_res$technical_component,
    cellwise_background_mean = dsb_res$cellwise_background_mean,
    isotype_names = isotype_names,
    dsb_params = dsb_params,
    scale_factor = scale_factor
  )

  structure(
    list(
      raw_counts = raw_counts,
      norm_counts = norm_counts,
      cell_info = cell_info,
      other_data = other_data,
      type = "DSB"
    ),
    class = "ADTCounts"
  )
}

### methods --------------------------------------------------------------------

#### primitives ----------------------------------------------------------------

#' @export
print.ADTCounts <- function(x, ...) {
  cat("ADTCounts\n")
  cat("  Cells:    ", nrow(x$raw_counts), "\n")
  cat("  Proteins: ", ncol(x$raw_counts), "\n")
  cat("  Type:     ", x$type, "\n")
  invisible(x)
}

#### getters -------------------------------------------------------------------

#' Get the ADT feature info
#'
#' @description
#' Returns the number of cells expressing a given ADT
#'
#' @param x The object from which to get the ADT feature information
#'
#' @returns A data.table with the ADT feature information
#'
#' @export
get_adt_feature_info <- function(x) {
  UseMethod("get_adt_feature_info")
}

#' @rdname get_adt_feature_info
#'
#' @export
get_adt_feature_info.ADTCounts <- function(x) {
  # checks
  checkmate::assertClass(x, "ADTCounts")

  # body
  raw_counts <- x$raw_counts

  res <- data.table::data.table(
    feature_idx = seq_along(colnames(raw_counts)),
    feature_id = colnames(raw_counts),
    nnz = colSums(raw_counts > 0)
  )

  return(res)
}

#' Get the ADT sample info
#'
#' @description
#' Returns sample information in terms of ADT capture, number of features, etc.
#'
#' @param x The object from which to get the ADT sample information
#'
#' @returns A data.table with the ADT sample information
#'
#' @export
get_adt_sample_info <- function(x) {
  UseMethod("get_adt_sample_info")
}

#' @rdname get_adt_sample_info
#'
#' @export
get_adt_sample_info.ADTCounts <- function(x) {
  # checks
  checkmate::assertClass(x, "ADTCounts")

  # body
  res <- data.table::data.table(
    cell_idx = x$cell_info,
    adt_nnz = rowSums(x$raw_counts != 0),
    adt_lib_size = rowSums(x$raw_counts)
  )

  .add_is_obs_attr(res)

  return(res)
}


# s7 ---------------------------------------------------------------------------

## helpers ---------------------------------------------------------------------

#' Helper function to generate ADT vars
#'
#' @param adt_counts Numeric matrix. Of shape samples x features with the ADT
#' raw counts.
#'
#' @returns A data.table with the feature_idx and feature_id of the ADT data.
#'
#' @keywords internal
.generate_adt_var <- function(adt_counts) {
  # checks
  checkmate::assertMatrix(
    adt_counts,
    mode = "numeric",
    row.names = "named",
    col.names = "named"
  )

  data.table::data.table(
    feature_idx = seq_along(colnames(adt_counts)),
    feature_id = colnames(adt_counts)
  )
}

#' Helper to get the cache slot
#'
#' @param modality String. One of `c("rna", "adt")`.
#'
#' @returns The name of the cache.
#'
#' @keywords internal
.cache_slot_from_modality <- function(modality) {
  modality <- match.arg(modality, c("rna", "adt"))

  switch(modality, rna = "sc_cache", adt = "adt_cache")
}

#' Helper to get the ADT feature names
#'
#' @param x The class from which to extract the ADT feature names.
#'
#' @returns The feature names if available
#'
#' @keywords internal
.adt_feature_names <- function(x) {
  adt <- S7::prop(x, "adt_counts")
  if (is.null(adt)) {
    return(character(0))
  }
  colnames(adt$raw_counts)
}

## single cell class (multi modal) ---------------------------------------------

#' @title bixverse SingleCells (multi modal) class
#'
#' @description
#' This is the `bixverse`-based SingleCells class for multiple modalities. Under
#' the hood it uses a DuckDB for obs and vars storing, and a Rust-based
#' binarised file format to store the raw and normalised counts for single cell
#' RNAseq. In both cases, the idea is not to hold any data that is not needed at
#' a given point of time in memory, but leverage speedy on-disk computations and
#' streaming engines powered by Rust and DuckDB to run the analysis.
#'
#' @param dir_data String. This is the directory in which the experimental files
#' will be
#'
#' @section Properties:
#' \describe{
#'   \item{db_connection}{This contains an R6 class with DuckDB pointers and
#'   wrappers to interact with the table-like data for this experiment.}
#'   \item{count_connection}{This contains an R6-like environment that points
#'   to Rust functions that can work on the RNAseq counts more specifically.}
#'   \item{adt_counts}{...}
#'   \item{peak_connection}{Future feature: to store ATAC Seq counts in the
#'   future.}
#'   \item{dir_data}{Path to the directory in which the data will be saved on
#'   disk.}
#'   \item{sc_cache}{Class with cached data. Contains less memory-heavy objects
#'   such as embeddings, kNN information or sNN graphs for the single cell
#'   RNAseq.}
#'   \item{adt_cache}{Class with cached data. Contains less memory-heavy objects
#'   such as embeddings, kNN information or sNN graphs for the single cell
#'   Antibody-Derived Tags.}
#'   \item{atac_cache}{Class with cached data. Contains less memory-heavy objects
#'   such as embeddings, kNN information or sNN graphs for the single cell
#'   chromatin accessability.}
#'   \item{sc_map}{Class containing various mapping information such as HVG
#'   indices, cells to keep, etc.}
#'   \item{other_data}{List that contains additional data and results, such
#'   as for example the WNN graph.}
#'   \item{dims}{Dimensions of the original data.}
#' }
#'
#' @return Returns the `SingleCellsMultiModal` class for further operations.
#'
#' @export
SingleCellsMultiModal <- S7::new_class(
  name = "SingleCellsMultiModal",
  parent = SingleCells,
  properties = list(
    db_connection = S7::class_any,
    count_connection = S7::class_any,
    adt_counts = S7::class_any,
    peak_connection = S7::class_any,
    dir_data = S7::class_character,
    sc_cache = S7::class_any,
    adt_cache = S7::class_any,
    atac_cache = S7::class_any,
    sc_map = S7::class_any,
    other_data = S7::class_list,
    dims = S7::class_integer
  ),
  constructor = function(dir_data) {
    # checks
    checkmate::assertDirectoryExists(dir_data)

    dir_data <- path.expand(dir_data)

    # generate the Rust pointer
    count_connection <- SingleCellCountData$new(
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
      adt_counts = NULL,
      peak_connection = NULL,
      dir_data = path.expand(dir_data),
      sc_cache = new_sc_cache(),
      adt_cache = new_sc_cache(),
      atac_cache = new_sc_cache(),
      sc_map = new_sc_mapper(),
      other_data = list(),
      dims = c(0L, 0L)
    )
  }
)

## primitives ------------------------------------------------------------------

#' @name print.SingleCellsMultiModal
#'
#' @title print Method for SingleCellsMultiModal object
#'
#' @description
#' Print a SingleCellsMultiModal object.
#'
#' @param x An object of class `SingleCellsMultiModal`.
#' @param ... Additional arguments (currently not used).
#'
#' @returns Invisibly returns `x`.
#'
#' @method print SingleCellsMultiModal
#'
#' @keywords internal
S7::method(print, SingleCellsMultiModal) <- function(x, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsMultiModal))

  dims <- S7::prop(x, "dims")
  no_cells_to_keep <- length(get_cells_to_keep(x))
  sc_map <- S7::prop(x, "sc_map")
  sc_cache <- S7::prop(x, "sc_cache")
  adt_cache <- S7::prop(x, "adt_cache")
  adt_counts <- S7::prop(x, "adt_counts")

  # RNA status
  rna_hvg <- !is.null(sc_map[["hvg_gene_indices"]])
  rna_pca <- !is.null(sc_cache[["pca_factors"]])
  rna_knn <- !is.null(sc_cache[["knn"]])
  rna_snn <- !is.null(sc_cache[["snn_graph"]])
  rna_other <- names(sc_cache[["other_embeddings"]])
  rna_other_str <- if (length(rna_other) > 0) {
    paste(rna_other, collapse = ", ")
  } else {
    "none"
  }

  # ADT status
  adt_present <- !is.null(adt_counts)
  adt_n_features <- if (adt_present) ncol(adt_counts$raw_counts) else 0L
  adt_pca <- !is.null(adt_cache[["pca_factors"]])
  adt_knn <- !is.null(adt_cache[["knn"]])
  adt_snn <- !is.null(adt_cache[["snn_graph"]])
  adt_other <- names(adt_cache[["other_embeddings"]])
  adt_other_str <- if (length(adt_other) > 0) {
    paste(adt_other, collapse = ", ")
  } else {
    "none"
  }

  cat(
    "Single cell experiment (Multi-modal).\n",
    sprintf("  No cells (original): %i\n", dims[1]),
    sprintf("   To keep n: %i\n", no_cells_to_keep),
    "  RNA:\n",
    sprintf("    No genes: %i\n", dims[2]),
    sprintf("    HVG calculated: %s\n", rna_hvg),
    sprintf("    PCA calculated: %s\n", rna_pca),
    sprintf("    Other embeddings: %s\n", rna_other_str),
    sprintf("    KNN generated: %s\n", rna_knn),
    sprintf("    SNN generated: %s\n", rna_snn),
    "  ADT:\n",
    sprintf("    Present: %s\n", adt_present),
    sprintf("    No features: %i\n", adt_n_features),
    sprintf("    PCA calculated: %s\n", adt_pca),
    sprintf("    Other embeddings: %s\n", adt_other_str),
    sprintf("    KNN generated: %s\n", adt_knn),
    sprintf("    SNN generated: %s\n", adt_snn),
    "  ATAC: not yet implemented\n",
    sep = ""
  )

  invisible(x)
}

## methods ---------------------------------------------------------------------

### getters --------------------------------------------------------------------

#' @method get_sc_var SingleCellsMultiModal
#'
#' @export
S7::method(get_sc_var, SingleCellsMultiModal) <- function(
  object,
  indices = NULL,
  cols = NULL,
  modality = c("rna", "adt")
) {
  modality <- match.arg(modality)

  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(indices, c("0", "I+"))
  checkmate::assertChoice(modality, c("rna", "adt"))

  duckdb_con <- get_sc_duckdb(object)

  var_table <- switch(
    modality,
    rna = duckdb_con$get_vars_table(indices = indices, cols = cols),
    adt = duckdb_con$get_vars_adt_table(indices = indices, cols = cols)
  )

  return(var_table)
}

### adt ------------------------------------------------------------------------

#### add function --------------------------------------------------------------

#' Add ADT counts to `SingleCellsMultiModal`
#'
#' @description
#' This method allows you to add ADT counts to a `SingleCellsMultiModal`.
#' Assumes the transcriptomics data has already been ingested and the cells
#' to keep are known. The method subsets the ADT counts to the kept cells,
#' applies the requested normalisation (CLR or DSB), populates the `"var_adt"`
#' table in the DuckDB, and attaches an `ADTCounts` to the class.
#'
#' @param object `SingleCellsMultiModal` class.
#' @param adt_counts Numeric matrix. Cells x features matrix of raw ADT counts.
#' @param method String. One of `c("clr", "dsb")`. Normalisation method.
#' @param ... Additional arguments forwarded to the normalisation constructor.
#' For `method = "clr"`: `clean_clr_counts`, `percentile`. For `method = "dsb"`:
#' `empty_drops`, `isotype_names`, `dsb_params`, `scale_factor`, `seed`,
#' `verbose`. See [new_adt_counts_clr()] and [new_adt_counts_dsb()].
#'
#' @returns Returns a `SingleCellsMultiModal` with the ADT data added.
#'
#' @export
add_adt_counts_sc <- S7::new_generic(
  name = "add_adt_counts_sc",
  dispatch_args = "object",
  fun = function(
    object,
    adt_counts,
    method = c("clr", "dsb"),
    ...
  ) {
    S7::S7_dispatch()
  }
)

#' @method add_adt_counts_sc SingleCellsMultiModal
#'
#' @export
S7::method(add_adt_counts_sc, SingleCellsMultiModal) <- function(
  object,
  adt_counts,
  method = c("clr", "dsb"),
  ...
) {
  method <- match.arg(method)

  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsMultiModal))
  checkmate::assertMatrix(
    adt_counts,
    mode = "numeric",
    row.names = "named",
    col.names = "named"
  )

  cell_info <- get_cell_info(object, filtered = TRUE)
  checkmate::assertSubset(names(cell_info), rownames(adt_counts))
  checkmate::assertTRUE(length(cell_info) >= 1)

  adt_obj <- switch(
    method,
    clr = new_adt_counts_clr(
      raw_counts = adt_counts,
      cell_info = cell_info,
      ...
    ),
    dsb = new_adt_counts_dsb(
      raw_counts = adt_counts,
      cell_info = cell_info,
      ...
    )
  )

  adt_vars <- get_adt_feature_info(x = adt_obj)
  adt_obs <- get_adt_sample_info(x = adt_obj)

  duckdb_con <- get_sc_duckdb(object)
  duckdb_con$populate_var_adt_from_data.table(var_dt = adt_vars)

  object <- add_sc_new_obs(object, adt_obs)

  S7::prop(object, "adt_counts") <- adt_obj

  return(object)
}

#### counts --------------------------------------------------------------------

S7::method(get_sc_counts, SingleCellsMultiModal) <- function(
  object,
  assay = c("raw", "norm"),
  return_format = c("cell", "gene"),
  cell_indices = NULL,
  gene_indices = NULL,
  use_cells_to_keep = TRUE,
  modality = c("rna", "adt"),
  .verbose = TRUE
) {
  modality <- match.arg(modality)
  assay <- match.arg(assay)
  return_format <- match.arg(return_format)

  if (modality == "rna") {
    rna_method <- S7::method(get_sc_counts, SingleCells)
    return(rna_method(
      object = object,
      assay = assay,
      return_format = return_format,
      cell_indices = cell_indices,
      gene_indices = gene_indices,
      use_cells_to_keep = use_cells_to_keep,
      modality = "rna",
      .verbose = .verbose
    ))
  }

  # ADT path
  checkmate::qassert(cell_indices, c("0", "I+"))
  checkmate::qassert(gene_indices, c("0", "I+"))

  adt <- S7::prop(object, "adt_counts")
  if (is.null(adt)) {
    stop("No ADT counts in this object. Add them with add_adt_counts_sc().")
  }

  mat <- if (assay == "raw") adt$raw_counts else adt$norm_counts

  cells_to_keep_1idx <- get_cells_to_keep(object) + 1L

  if (!is.null(cell_indices)) {
    adt_rows <- match(cell_indices, cells_to_keep_1idx)
    if (any(is.na(adt_rows))) {
      warning(sprintf(
        "%i cell_indices not in the kept-cell set; dropping.",
        sum(is.na(adt_rows))
      ))
      adt_rows <- adt_rows[!is.na(adt_rows)]
    }
    if (length(adt_rows) == 0L) {
      stop("No valid cell_indices for the ADT matrix.")
    }
    mat <- mat[adt_rows, , drop = FALSE]
  } else {
    mat <- mat[cells_to_keep_1idx, , drop = FALSE]
  }

  if (!is.null(gene_indices)) {
    mat <- mat[, gene_indices, drop = FALSE]
  }

  mat
}

S7::method(`[`, SingleCellsMultiModal) <- function(
  x,
  i,
  j,
  ...,
  assay = c("raw", "norm"),
  return_format = c("cell", "gene"),
  use_cells_to_keep = TRUE,
  modality = c("rna", "adt"),
  drop = TRUE
) {
  modality <- match.arg(modality)
  if (missing(i)) {
    i <- NULL
  }
  if (missing(j)) {
    j <- NULL
  }

  if (checkmate::qtest(i, "S+")) {
    i <- get_cell_indices(x = x, cell_ids = i, rust_index = FALSE)
  }

  if (checkmate::qtest(j, "S+")) {
    if (modality == "rna") {
      j <- get_gene_indices(x = x, gene_ids = j, rust_index = FALSE)
    } else {
      feature_names <- .adt_feature_names(x)
      j_match <- match(j, feature_names)
      if (any(is.na(j_match))) {
        warning(sprintf(
          "%i ADT feature names not found; dropping.",
          sum(is.na(j_match))
        ))
        j_match <- j_match[!is.na(j_match)]
      }
      j <- as.integer(j_match)
    }
  }

  assay <- match.arg(assay)
  return_format <- match.arg(return_format)
  checkmate::qassert(i, c("I+", "0"))
  checkmate::qassert(j, c("I+", "0"))

  get_sc_counts(
    object = x,
    assay = assay,
    return_format = return_format,
    cell_indices = i,
    gene_indices = j,
    use_cells_to_keep = use_cells_to_keep,
    modality = modality,
    .verbose = FALSE
  )
}

#### setters -------------------------------------------------------------------

S7::method(set_pca_factors, SingleCellsMultiModal) <- function(
  x,
  pca_factor,
  modality = c("rna", "adt"),
  ...
) {
  slot <- .cache_slot_from_modality(modality)
  checkmate::assertMatrix(pca_factor, mode = "numeric")
  S7::prop(x, slot) <- set_pca_factors(S7::prop(x, slot), pca_factor)
  x
}

S7::method(set_pca_loadings, SingleCellsMultiModal) <- function(
  x,
  pca_loading,
  modality = c("rna", "adt"),
  ...
) {
  slot <- .cache_slot_from_modality(modality)
  checkmate::assertMatrix(pca_loading, mode = "numeric")
  S7::prop(x, slot) <- set_pca_loadings(S7::prop(x, slot), pca_loading)
  x
}

S7::method(set_pca_singular_vals, SingleCellsMultiModal) <- function(
  x,
  singular_vals,
  modality = c("rna", "adt"),
  ...
) {
  slot <- .cache_slot_from_modality(modality)
  checkmate::qassert(singular_vals, "N+")
  S7::prop(x, slot) <- set_pca_singular_vals(S7::prop(x, slot), singular_vals)
  x
}

S7::method(set_embedding, SingleCellsMultiModal) <- function(
  x,
  embd,
  name,
  modality = c("rna", "adt", "other"),
  ...
) {
  modality <- match.arg(modality)
  checkmate::assertMatrix(embd, mode = "numeric")
  checkmate::qassert(name, "S1")

  # wnn-derived embeddings live under other_data[["wnn"]][["embeddings"]]
  if (modality == "other") {
    other_data <- S7::prop(x, "other_data")
    other_data[[name]] <- embd
    S7::prop(x, "other_data") <- other_data
    return(x)
  }

  slot <- .cache_slot_from_modality(modality)
  S7::prop(x, slot) <- set_embedding(
    S7::prop(x, slot),
    embd = embd,
    name = name
  )
  x
}

S7::method(set_knn, SingleCellsMultiModal) <- function(
  x,
  knn,
  modality = c("rna", "adt"),
  ...
) {
  slot <- .cache_slot_from_modality(modality)
  checkmate::assertClass(knn, "SingleCellNearestNeighbour")
  S7::prop(x, slot) <- set_knn(S7::prop(x, slot), knn)
  x
}

S7::method(set_snn_graph, SingleCellsMultiModal) <- function(
  x,
  snn_graph,
  modality = c("rna", "adt"),
  ...
) {
  slot <- .cache_slot_from_modality(modality)
  checkmate::assertClass(snn_graph, "igraph")
  S7::prop(x, slot) <- set_snn_graph(S7::prop(x, slot), snn_graph)
  x
}

S7::method(remove_knn, SingleCellsMultiModal) <- function(
  x,
  modality = c("rna", "adt"),
  ...
) {
  slot <- .cache_slot_from_modality(modality)
  S7::prop(x, slot) <- remove_knn(S7::prop(x, slot))
  x
}

S7::method(remove_snn_graph, SingleCellsMultiModal) <- function(
  x,
  modality = c("rna", "adt"),
  ...
) {
  slot <- .cache_slot_from_modality(modality)
  S7::prop(x, slot) <- remove_snn_graph(S7::prop(x, slot))
  x
}

#### getters -------------------------------------------------------------------

S7::method(get_pca_factors, SingleCellsMultiModal) <- function(
  x,
  modality = c("rna", "adt"),
  ...
) {
  slot <- .cache_slot_from_modality(modality)
  res <- get_pca_factors(S7::prop(x, slot))
  if (is.null(res)) {
    return(NULL)
  }
  rownames(res) <- get_cell_names(x, filtered = TRUE)
  colnames(res) <- sprintf("PC_%i", seq_len(ncol(res)))
  res
}

S7::method(get_pca_loadings, SingleCellsMultiModal) <- function(
  x,
  modality = c("rna", "adt"),
  ...
) {
  modality <- match.arg(modality)
  slot <- .cache_slot_from_modality(modality)
  res <- get_pca_loadings(S7::prop(x, slot))
  if (is.null(res)) {
    return(NULL)
  }
  colnames(res) <- sprintf("PC_%i", seq_len(ncol(res)))
  rownames(res) <- if (modality == "rna") {
    get_gene_names(x)[get_hvg(x) + 1]
  } else {
    .adt_feature_names(x)
  }
  res
}

S7::method(get_pca_singular_val, SingleCellsMultiModal) <- function(
  x,
  modality = c("rna", "adt"),
  ...
) {
  slot <- .cache_slot_from_modality(modality)
  get_pca_singular_val(S7::prop(x, slot))
}

S7::method(get_embedding, SingleCellsMultiModal) <- function(
  x,
  embd_name,
  modality = c("rna", "adt", "other"),
  ...
) {
  modality <- match.arg(modality)

  if (modality == "other") {
    res <- S7::prop(x, "other_data")[[embd_name]]
    if (is.null(res)) {
      warning("The desired embedding is not available.")
    }
    return(res)
  }

  slot <- .cache_slot_from_modality(modality)
  checkmate::qassert(embd_name, "S1")
  res <- get_embedding(S7::prop(x, slot), embd_name = embd_name)
  rownames(res) <- get_cell_names(x, filtered = TRUE)
  res
}

S7::method(get_available_embeddings, SingleCellsMultiModal) <- function(
  x,
  modality = c("rna", "adt"),
  ...
) {
  slot <- .cache_slot_from_modality(modality)
  get_available_embeddings(S7::prop(x, slot))
}

S7::method(get_knn_mat, SingleCellsMultiModal) <- function(
  x,
  modality = c("rna", "adt"),
  ...
) {
  slot <- .cache_slot_from_modality(modality)
  get_knn_mat(S7::prop(x, slot))
}

S7::method(get_knn_dist, SingleCellsMultiModal) <- function(
  x,
  modality = c("rna", "adt"),
  ...
) {
  slot <- .cache_slot_from_modality(modality)
  get_knn_dist(S7::prop(x, slot))
}

S7::method(get_knn_obj, SingleCellsMultiModal) <- function(
  x,
  modality = c("rna", "adt", "wnn"),
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsMultiModal))
  checkmate::assertChoice(modality, c("rna", "adt", "wnn"))

  # special case WNN
  if (modality == "wnn") {
    wnn_data <- S7::prop(x, "other_data")[["wnn"]]

    if (is.null(wnn_data)) {
      warning(paste(
        "No WNN data found. Did you run generate_wnn_graph_sc()",
        "Returning NULL."
      ))
      return(NULL)
    }

    knn <- wnn_data[["wnn"]][["knn"]]

    return(knn)
  }

  slot <- .cache_slot_from_modality(modality)
  get_knn_obj(S7::prop(x, slot))
}

S7::method(get_snn_graph, SingleCellsMultiModal) <- function(
  x,
  modality = c("rna", "adt", "wnn"),
  ...
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsMultiModal))
  checkmate::assertChoice(modality, c("rna", "adt", "wnn"))

  # special case WNN
  if (modality == "wnn") {
    wnn_data <- S7::prop(x, "other_data")[["wnn"]]

    if (is.null(wnn_data)) {
      warning(paste(
        "No WNN data found. Did you run generate_wnn_graph_sc()",
        "Returning NULL."
      ))
      return(NULL)
    }

    snn <- wnn_data[["wnn"]][["snn"]]

    return(snn)
  }

  slot <- .cache_slot_from_modality(modality)
  get_snn_graph(S7::prop(x, slot))
}
