# single cell atac -------------------------------------------------------------

## globals ---------------------------------------------------------------------

# scale of the TF-IDF term frequencies, as Signac and muon use it
.ATAC_TF_IDF_SCALE <- 1e4

# names of the mitochondrial genome across annotations
.ATAC_MITO_CHROMS <- c("chrM", "MT", "chrMT", "M")

# warn when more of the requested cells than this have no fragment
.ATAC_MISSING_CELLS_WARN <- 0.01

# warn when more of the fragments than this lie on chromosomes without targets
.ATAC_OFF_TARGET_WARN <- 0.5

# what one count of an ATAC store is
.ATAC_COUNT_MODES <- c("insertion", "paired-insertion", "fragment")

# store metadata key of the count mode
.ATAC_COUNT_MODE_KEY <- "count_mode"

# the LSI fields of `AtacInfo`, which only the memory file carries
.ATAC_LSI_FIELDS <- c(
  "depth_cor",
  "lsi_dropped",
  "lsi_features",
  "lsi_loadings",
  "lsi_singular_values"
)

## atac info -------------------------------------------------------------------

#' Generate a new `AtacInfo` object
#'
#' @description
#' What an ATAC store holds and how far it was processed. The kind, count mode
#' and TF-IDF flag are read from the store; the LSI fields come from
#' [bixverse::calculate_lsi_sc()] and are saved with the object.
#'
#' @param kind String. One of `c("peak", "tile")`.
#' @param count_mode String. One of `c("insertion", "paired-insertion",
#' "fragment")`. What one count is: up to two insertions per fragment, one per
#' fragment and feature for the latter two.
#' @param tf_idf Boolean. Whether the store is TF-IDF normalised.
#'
#' @returns An `AtacInfo` object, a list with
#' \itemize{
#'   \item kind - `"peak"` or `"tile"`.
#'   \item count_mode - What one count is.
#'   \item tf_idf - Whether the store is TF-IDF normalised.
#'   \item depth_cor - Correlation of each LSI component with the library
#'   size, all components including any dropped one. `NULL` before the LSI.
#'   \item lsi_dropped - 0-based LSI components left out of the `"lsi"`
#'   embedding.
#'   \item lsi_features - 0-based feature rows the LSI was fitted on.
#'   \item lsi_loadings - Features x components loadings of the LSI.
#'   \item lsi_singular_values - Singular values of the LSI.
#' }
#'
#' @keywords internal
new_atac_info <- function(kind, count_mode, tf_idf = FALSE) {
  checkmate::assertChoice(kind, c("peak", "tile"))
  checkmate::assertChoice(count_mode, .ATAC_COUNT_MODES)
  checkmate::qassert(tf_idf, "B1")

  structure(
    list(
      kind = kind,
      count_mode = count_mode,
      tf_idf = tf_idf,
      depth_cor = NULL,
      lsi_dropped = integer(),
      lsi_features = NULL,
      lsi_loadings = NULL,
      lsi_singular_values = NULL
    ),
    class = "AtacInfo"
  )
}

#' @export
print.AtacInfo <- function(x, ...) {
  cat(
    "AtacInfo\n",
    sprintf("  kind: %s\n", x$kind),
    sprintf("  count mode: %s\n", x$count_mode),
    sprintf("  TF-IDF: %s\n", x$tf_idf),
    sprintf("  LSI: %s\n", !is.null(x$depth_cor)),
    sep = ""
  )
  invisible(x)
}

#' Read the ATAC state a store records
#'
#' @description
#' Rebuilds the kind, count mode and TF-IDF flag from the store's header and
#' metadata. Stores written before the count mode was recorded are taken as
#' insertion counts, with a warning.
#'
#' @param store_dir String. The directory holding the store's two files.
#' @param lsi Optional list with the LSI fields of an earlier `AtacInfo`, as
#' saved in the memory file.
#'
#' @returns An `AtacInfo` object.
#'
#' @keywords internal
.read_atac_info <- function(store_dir, lsi = NULL) {
  checkmate::assertDirectoryExists(store_dir)
  checkmate::assertList(lsi, null.ok = TRUE)

  cells <- file.path(store_dir, "counts_cells.bin")
  kind <- rs_store_feature_kind(cells)
  if (!kind %in% c("peak", "tile")) {
    stop(sprintf(
      "%s holds %s features, not peaks or tiles.",
      store_dir,
      kind
    ))
  }
  mode <- rs_store_metadata(cells)[.ATAC_COUNT_MODE_KEY]
  if (is.na(mode) || !mode %in% .ATAC_COUNT_MODES) {
    # 10x counts and the add_atac_fragments_sc() default are insertions
    warning(sprintf(
      paste(
        "The ATAC store in %s records no count mode; assuming 'insertion'.",
        "Re-count it to record one."
      ),
      store_dir
    ))
    mode <- "insertion"
  }

  info <- new_atac_info(
    kind = kind,
    count_mode = unname(mode),
    tf_idf = rs_store_norm_kind(cells) == "tf_idf"
  )
  for (field in intersect(names(lsi), .ATAC_LSI_FIELDS)) {
    info[field] <- list(lsi[[field]])
  }

  info
}

#' Record the count mode in both files of an ATAC store
#'
#' @param store_dir String. The directory holding the store's two files.
#' @param count_mode String. One of `c("insertion", "paired-insertion",
#' "fragment")`.
#'
#' @returns Invisible `NULL`.
#'
#' @keywords internal
.write_count_mode <- function(store_dir, count_mode) {
  checkmate::assertDirectoryExists(store_dir)
  checkmate::assertChoice(count_mode, .ATAC_COUNT_MODES)

  rs_set_store_metadata(
    f_path_cells = file.path(store_dir, "counts_cells.bin"),
    f_path_genes = file.path(store_dir, "counts_genes.bin"),
    keys = .ATAC_COUNT_MODE_KEY,
    values = count_mode
  )

  invisible(NULL)
}

#' The LSI fields of an `AtacInfo`, for the memory file
#'
#' @param info `AtacInfo` or `NULL`.
#'
#' @returns A list with the LSI fields, or `NULL`.
#'
#' @keywords internal
.atac_lsi_payload <- function(info) {
  if (is.null(info)) {
    return(NULL)
  }
  checkmate::assertClass(info, "AtacInfo")

  info[.ATAC_LSI_FIELDS]
}

## class -----------------------------------------------------------------------

#' @title bixverse single cell ATAC class
#'
#' @description
#' A single cell ATAC experiment without RNA: a peak or genome tile store,
#' with the same on-disk layout as [bixverse::SingleCells()]. The var table
#' holds the features as `gene_id` (`chr:start-end`, 0-based half-open) with
#' `chrom`, `start` and `end`, so the usual store machinery reads it unchanged.
#' Create one with [bixverse::load_atac_fragments()] or
#' [bixverse::load_atac_tenx()], reopen one with
#' [bixverse::load_existing()].
#'
#' `modality = "atac"` and the default `"rna"` both name this object's own
#' store and cache, so code written for the ATAC side of a
#' [bixverse::SingleCellsMultiModal()] runs unchanged.
#'
#' @param dir_data String. The directory with the experimental files.
#'
#' @section Properties:
#' \describe{
#'   \item{db_connection}{R6 class with the DuckDB connection.}
#'   \item{count_connection}{Rust connection to the peak or tile store.}
#'   \item{dir_data}{Path to the directory with the data.}
#'   \item{sc_cache}{The ATAC cache: LSI and other embeddings, kNN, graph.}
#'   \item{sc_map}{Cell and feature mappings, the cell filter.}
#'   \item{dims}{Dimensions of the store, cells x features.}
#'   \item{atac_info}{`AtacInfo` with the store state and the LSI fit.}
#' }
#'
#' @returns Returns the `SingleCellsAtac` class for further operations.
#'
#' @export
#'
#' @examples
#' # an empty handle; fill it with load_atac_fragments() or load_atac_tenx()
#' dir <- tempfile("bixverse_atac")
#' dir.create(dir)
#' object <- SingleCellsAtac(dir_data = dir)
#'
#' unlink(dir, recursive = TRUE, force = TRUE)
SingleCellsAtac <- S7::new_class(
  name = "SingleCellsAtac",
  parent = SingleCells,
  properties = list(
    atac_info = S7::class_any
  ),
  constructor = function(dir_data) {
    checkmate::assertDirectoryExists(dir_data)

    S7::new_object(
      SingleCells(dir_data = dir_data),
      atac_info = NULL
    )
  }
)

## primitives ------------------------------------------------------------------

#' @noRd
S7::method(print, SingleCellsAtac) <- function(x, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsAtac))

  dims <- S7::prop(x, "dims")
  sc_cache <- S7::prop(x, "sc_cache")
  info <- S7::prop(x, "atac_info")
  other <- names(sc_cache[["other_embeddings"]])
  ga <- file.exists(file.path(
    S7::prop(x, "dir_data"),
    "gene_activity",
    "counts_cells.bin"
  ))

  cat(
    "Single cell experiment (ATAC).\n",
    sprintf("  No cells (original): %i\n", dims[1]),
    sprintf("   To keep n: %i\n", length(get_cells_to_keep(x))),
    sprintf("  No features: %i\n", dims[2]),
    sprintf(
      "  Kind: %s (%s counts)\n",
      info$kind %||% "unknown",
      info$count_mode %||% "unknown"
    ),
    sprintf("  TF-IDF: %s\n", isTRUE(info$tf_idf)),
    sprintf(
      "  Embeddings: %s\n",
      if (length(other)) paste(other, collapse = ", ") else "none"
    ),
    sprintf("  KNN generated: %s\n", !is.null(sc_cache[["knn"]])),
    sprintf("  Graph generated: %s\n", !is.null(sc_cache[["connectivities"]])),
    sprintf("  Gene activity: %s\n", ga),
    sprintf("  Stale artefacts: %s\n", .print_stale_str(x)),
    sep = ""
  )

  invisible(x)
}

## host helpers ----------------------------------------------------------------

# `SingleCellsAtac` and the ATAC side of `SingleCellsMultiModal` share every
# ATAC function. These helpers hold the few places where the two differ.

#' Directory of the ATAC store of an object
#'
#' @param x `SingleCellsAtac` or `SingleCellsMultiModal`.
#'
#' @returns String. The experiment directory for `SingleCellsAtac`, its `atac`
#' subdirectory for `SingleCellsMultiModal`.
#'
#' @keywords internal
.atac_store_dir <- function(x) {
  dir_data <- S7::prop(x, "dir_data")
  if (S7::S7_inherits(x, SingleCellsAtac)) {
    return(dir_data)
  }
  checkmate::assertTRUE(S7::S7_inherits(x, SingleCellsMultiModal))

  file.path(dir_data, "atac")
}

#' Paths of the two files of the ATAC store
#'
#' @param x `SingleCellsAtac` or `SingleCellsMultiModal`.
#'
#' @returns A list with `cells` and `genes`, the cell- and feature-based files.
#'
#' @keywords internal
.atac_paths <- function(x) {
  store_dir <- .atac_store_dir(x)

  list(
    cells = file.path(store_dir, "counts_cells.bin"),
    genes = file.path(store_dir, "counts_genes.bin")
  )
}

#' The modality name of the ATAC cache of an object
#'
#' @param x `SingleCellsAtac` or `SingleCellsMultiModal`.
#'
#' @returns String. `"rna"` for `SingleCellsAtac`, whose own cache is the ATAC
#' one, `"atac"` otherwise.
#'
#' @keywords internal
.atac_modality <- function(x) {
  if (S7::S7_inherits(x, SingleCellsAtac)) "rna" else "atac"
}

#' The obs column with the ATAC library size
#'
#' @param x `SingleCellsAtac` or `SingleCellsMultiModal`.
#'
#' @returns String.
#'
#' @keywords internal
.atac_lib_col <- function(x) {
  if (S7::S7_inherits(x, SingleCellsAtac)) "lib_size" else "atac_lib_size"
}

#' Resolve a modality against the object it is used on
#'
#' @description
#' On a `SingleCellsAtac` `"atac"` names the object's own cache, i.e.
#' `"rna"`. Any modality other than `"rna"` on an object without other
#' modalities is an error.
#'
#' @param x The single cell object.
#' @param modality String. The requested modality.
#'
#' @returns String. The modality to use on `x`.
#'
#' @keywords internal
.resolve_modality <- function(x, modality) {
  checkmate::qassert(modality, "S1")
  if (S7::S7_inherits(x, SingleCellsAtac)) {
    if (!modality %in% c("rna", "atac")) {
      stop(sprintf("A SingleCellsAtac has no '%s' modality.", modality))
    }
    return("rna")
  }
  if (modality != "rna" && !S7::S7_inherits(x, SingleCellsMultiModal)) {
    stop(sprintf(
      paste(
        "modality = '%s' needs a SingleCellsMultiModal (or a SingleCellsAtac",
        "for 'atac')."
      ),
      modality
    ))
  }

  modality
}

#' The ATAC state of an object
#'
#' @param x `SingleCellsAtac` or `SingleCellsMultiModal`.
#'
#' @returns The `AtacInfo`. Errors when no ATAC store was added.
#'
#' @keywords internal
.require_atac <- function(x) {
  info <- S7::prop(x, "atac_info")
  if (is.null(info)) {
    stop(paste(
      "No ATAC data in this object. Add it with add_atac_tenx_sc() or",
      "add_atac_fragments_sc()."
    ))
  }

  info
}

#' The ATAC feature table
#'
#' @param x `SingleCellsAtac` or `SingleCellsMultiModal`.
#' @param cols Optional character vector of columns, in the
#' `feature_idx` / `feature_id` naming.
#'
#' @returns A data.table keyed `feature_idx` / `feature_id` (1-based index).
#'
#' @keywords internal
.atac_var <- function(x, cols = NULL) {
  checkmate::qassert(cols, c("0", "S+"))
  duckdb_con <- get_sc_duckdb(x)

  if (S7::S7_inherits(x, SingleCellsMultiModal)) {
    return(duckdb_con$get_vars_atac_table(cols = cols))
  }

  rename <- c(feature_idx = "gene_idx", feature_id = "gene_id")
  if (!is.null(cols)) {
    cols <- ifelse(cols %in% names(rename), rename[cols], cols)
  }
  var <- duckdb_con$get_vars_table(cols = cols)
  old <- intersect(rename, names(var))
  data.table::setnames(var, old, names(rename)[match(old, rename)])

  var
}

#' Add or replace columns of the ATAC feature table
#'
#' @param x `SingleCellsAtac` or `SingleCellsMultiModal`.
#' @param new_data data.table with one row per feature, in store order.
#'
#' @returns Invisible `NULL`.
#'
#' @keywords internal
.add_atac_var <- function(x, new_data) {
  checkmate::assertDataTable(new_data)
  duckdb_con <- get_sc_duckdb(x)

  if (S7::S7_inherits(x, SingleCellsAtac)) {
    duckdb_con$add_data_var(new_data = data.table::copy(new_data))
    return(invisible(NULL))
  }

  var <- duckdb_con$get_vars_atac_table()
  checkmate::assertTRUE(nrow(var) == nrow(new_data))
  drop <- intersect(c("feature_idx", names(new_data)), names(var))
  if (length(drop)) {
    var[, (drop) := NULL]
  }
  duckdb_con$populate_var_atac_from_data.table(cbind(var, new_data))

  invisible(NULL)
}

#' Copy an object's cell filter into its default gene activity store
#'
#' @description
#' The gene activity store at `<dir>/gene_activity` shares the cell axis of
#' the experiment; without this, a filter set after it was written would not
#' reach it.
#'
#' @param x The single cell object owning the gene activity.
#'
#' @returns Invisible `NULL`.
#'
#' @keywords internal
.sync_gene_activity_filter <- function(x) {
  ga_dir <- file.path(S7::prop(x, "dir_data"), "gene_activity")
  if (!file.exists(file.path(ga_dir, "counts_cells.bin"))) {
    return(invisible(NULL))
  }
  ga_db <- SingleCellDuckDB$new(db_dir = ga_dir)
  ga_db$set_cells_to_keep(cell_idx_to_keep = get_cells_to_keep(x) + 1L)

  invisible(NULL)
}

#' Stamp the cell IDs of a freshly written store into its headers
#'
#' @description
#' Only their hash is stored. Stores that share a cell axis, such as a gene
#' activity store and its source, carry the same one, which is how
#' [bixverse::get_gene_activity()] detects a re-ingest.
#'
#' @param dir_data String. The directory with the store and its DuckDB.
#'
#' @returns Invisible `NULL`.
#'
#' @keywords internal
.stamp_store_cell_ids <- function(dir_data) {
  checkmate::assertDirectoryExists(dir_data)

  duckdb_con <- SingleCellDuckDB$new(db_dir = dir_data)
  rs_set_cell_ids(
    f_path_cells = file.path(dir_data, "counts_cells.bin"),
    f_path_genes = file.path(dir_data, "counts_genes.bin"),
    ids = as.character(names(duckdb_con$get_obs_index_map()))
  )

  invisible(NULL)
}
