# spatial i/o, h5ad -----------------------------------------------------------

# Separate from `methods_sp_io.R` on purpose. The counts, obs and var of an
# h5ad are the single cell reader's job and `load_h5ad()` already does it, so
# the only thing this file adds is `obsm/spatial` and whatever survived in
# `uns/spatial`. The reuse is a function call, not a fork.
#
# Index spaces, in the order they matter here:
#
# | thing                        | base | space                              |
# |------------------------------|------|------------------------------------|
# | rows of `obsm/spatial`       | 1    | original on-disk obs order         |
# | DuckDB `cell_idx`            | 1    | global, after QC filtering         |
# | `get_spot_indices_for_exp()` | 0    | global                             |
# | coordinates                  | n/a  | full-res pixels, (x, y)            |

## helpers ---------------------------------------------------------------------

#' Read the full obs identifier vector of an h5ad
#'
#' @description
#' The identifiers of **every** row on disk, before any QC filtering. Needed to
#' line `obsm/spatial` up with the obs table [bixverse::load_h5ad()] leaves
#' behind, which holds only the rows that survived.
#'
#' Mirrors the precedence [bixverse::SingleCellDuckDB()] uses when it populates
#' `cell_id`: the `_index` first, then a caller-named column.
#'
#' @param h5_path String. Path to the h5ad file.
#' @param cell_id_col Optional string. Obs column holding the identifiers, for
#' files that carry no `_index`.
#'
#' @return Character vector, one entry per row of the on-disk obs table.
#'
#' @keywords internal
.sp_h5ad_obs_ids <- function(h5_path, cell_id_col = NULL) {
  checkmate::assertFileExists(h5_path)
  checkmate::qassert(cell_id_col, c("0", "S1"))
  on.exit(tryCatch(rhdf5::h5closeAll(), error = function(e) invisible()))

  h5_content <- rhdf5::h5ls(h5_path) |> data.table::setDT()
  idx <- .resolve_h5_index(h5_path, "/obs", h5_content)$idx

  if (is.null(idx) && !is.null(cell_id_col)) {
    idx <- as.vector(rhdf5::h5read(h5_path, paste0("/obs/", cell_id_col)))
  }

  if (is.null(idx)) {
    stop(paste(
      "The h5ad carries no `_index` in /obs, so the spot coordinates cannot",
      "be matched back onto the loaded cells. Pass `cell_id_col` naming the",
      "obs column that holds the barcodes."
    ))
  }

  return(as.character(idx))
}

#' Map the loaded obs rows back onto the rows of `obsm/spatial`
#'
#' @description
#' This is the load-bearing join of the h5ad reader, and it is the same trap
#' [bixverse::load_visium()] has with `tissue_positions.csv`.
#' [bixverse::load_h5ad()] drops cells that fail QC, so the obs table is a
#' subset of the file in file order while `obsm/spatial` still has every row.
#' Zipping the two by position silently shifts every coordinate from the first
#' dropped cell onwards and nothing downstream notices.
#'
#' @param obs_ids Character vector. `cell_id` of the loaded rows, in ascending
#' `cell_idx` order.
#' @param all_ids Character vector. Identifiers of every row on disk.
#' @param n_coords Integer. Rows in `obsm/spatial`.
#'
#' @return Integer vector of 1-based positions into `obsm/spatial`, one per
#' loaded obs row.
#'
#' @keywords internal
.sp_h5ad_coord_rows <- function(obs_ids, all_ids, n_coords) {
  checkmate::qassert(obs_ids, "S+")
  checkmate::qassert(all_ids, "S+")
  checkmate::qassert(n_coords, "I1[1,)")

  if (length(all_ids) != n_coords) {
    stop(sprintf(
      paste(
        "The h5ad has %i obs rows but `obsm/spatial` has %i. The file is",
        "inconsistent and the coordinates cannot be trusted."
      ),
      length(all_ids),
      n_coords
    ))
  }

  if (anyDuplicated(all_ids) > 0L) {
    stop(paste(
      "The obs identifiers of the h5ad are not unique, so the coordinates",
      "cannot be matched onto the loaded cells without ambiguity."
    ))
  }

  rows <- match(obs_ids, all_ids)
  if (anyNA(rows)) {
    stop(sprintf(
      "%i loaded cell(s) have no row in the h5ad obs index.",
      sum(is.na(rows))
    ))
  }

  return(as.integer(rows))
}

#' Turn the crate's scale factor pair into a named list
#'
#' @param spatial List. Output of [bixverse::rs_sp_read_h5ad_spatial()].
#'
#' @return Named list of single numerics. Empty when the file shipped none.
#'
#' @keywords internal
.sp_h5ad_scale_factors <- function(spatial) {
  values <- as.numeric(spatial$scale_factor_values)
  if (length(values) == 0L) {
    return(list())
  }
  stats::setNames(as.list(values), as.character(spatial$scale_factor_names))
}

#' Which graph layouts a sample can actually use
#'
#' @description
#' `obsm/spatial` gives pixel coordinates but not the lattice indices, and
#' those only turn up in `obs` when the file happened to keep them. `"knn"` and
#' `"radius"` work off the coordinates alone and are therefore always
#' available; `"hex"` and `"square"` are exact table look-ups over
#' `array_row`/`array_col` and are not.
#'
#' @param has_array_indices Boolean. Whether the obs table carries both.
#'
#' @return Character vector of layout names usable with
#' [bixverse::params_sp_graph()].
#'
#' @keywords internal
.sp_available_layouts <- function(has_array_indices) {
  checkmate::qassert(has_array_indices, "B1")

  if (has_array_indices) {
    c("hex", "square", "knn", "radius")
  } else {
    c("knn", "radius")
  }
}

#' Say what the reader concluded about the coordinates
#'
#' @description
#' The orientation and the available layouts are the two things a caller cannot
#' see from the returned object but will trip over later, so they get said at
#' load time.
#'
#' @param spatial List. Output of [bixverse::rs_sp_read_h5ad_spatial()].
#' @param exp_id String. The experiment identifier.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return Invisibly `NULL`. Warns when the orientation was not settled by the
#' file.
#'
#' @keywords internal
.sp_h5ad_report <- function(spatial, exp_id, .verbose = TRUE) {
  layouts <- .sp_available_layouts(isTRUE(spatial$has_array_indices))

  if (.verbose) {
    message(sprintf(
      " Coordinates read as (%s), settled by: %s.",
      switch(
        spatial$orientation,
        xy = "column 0 = x",
        yx = "column 0 = y"
      ),
      switch(
        spatial$evidence,
        image_tissue = "the spots land on the tissue in the shipped image",
        image_frame = "only this order fits inside the shipped image",
        obs_pixel_columns = "the obs pxl_*_in_fullres columns",
        assumed = "nothing, the parameter was taken"
      )
    ))
    message(sprintf(
      " Graph layouts available for '%s': %s.",
      exp_id,
      toString(layouts)
    ))
    if (length(spatial$obsm_keys) > 1L) {
      message(sprintf(
        " obsm also holds %s, which were not read.",
        toString(setdiff(spatial$obsm_keys, "spatial"))
      ))
    }
  }

  if (identical(spatial$evidence, "assumed")) {
    warning(sprintf(
      paste(
        "Nothing in the h5ad settles the column order of `obsm/spatial`, so",
        "it was assumed to be (%s). Every spatial statistic is unaffected by",
        "this (a swap is a reflection and preserves every distance), but",
        "`image_features_sp()` would cut tiles from the transpose. Set",
        "`assume_orientation` in `params_sp_h5ad_io()` if you know better."
      ),
      spatial$orientation
    ))
  }

  if (!isTRUE(spatial$has_array_indices)) {
    warning(sprintf(
      paste(
        "The obs table of '%s' has no array_row/array_col, so the 'hex' and",
        "'square' graph layouts are unavailable. Use 'knn' or 'radius'."
      ),
      exp_id
    ))
  }

  invisible(NULL)
}

## methods ---------------------------------------------------------------------

#' Load an h5ad with spatial coordinates into `SpatialSpot`
#'
#' @description
#' Reads an h5ad carrying `obsm/spatial` into a [bixverse::SpatialSpot()].
#' Counts, obs and var go through [bixverse::load_h5ad()] unchanged; this adds
#' the coordinates, the scale factors that survived and the sample
#' registration.
#'
#' Written for the file academics actually share: coordinates yes, images
#' usually stripped, scale factors sometimes. A sample with nothing but
#' coordinates registers fine and runs graph construction, Moran's I, SPARK-X
#' and neighbourhood enrichment. Only [bixverse::image_features_sp()] needs
#' more, and it says so at call time.
#'
#' @param object `SpatialSpot` class.
#' @param h5_path String. Path to the `.h5ad` file. Must carry
#' `obsm/spatial`; a bare `filtered_feature_bc_matrix.h5` does not and errors
#' with a message saying so.
#' @param sp_io_param List. Output of [bixverse::params_sp_h5ad_io()]. A list
#' with the following elements:
#' \itemize{
#'   \item library_id - String or `NULL`. Which `uns/spatial` library to read.
#'   \item assume_orientation - String. Column order to fall back on.
#'   \item in_tissue_only - Boolean. Keep only `in_tissue == 1` spots.
#'   \item technology - String. Recorded on the registered sample.
#' }
#' @param sc_qc_param List. Output of [bixverse::params_sc_min_quality()]. A
#' list with the following elements:
#' \itemize{
#'   \item min_unique_genes - Integer. Minimum genes detected in a spot.
#'   \item min_lib_size - Integer. Minimum library size in a spot.
#'   \item min_cells - Integer. Minimum spots a gene must be detected in.
#'   \item target_size - Float. Target size to normalise to.
#' }
#' @param exp_id String. The experiment identifier written into the obs table
#' and used as the key of the resulting `SpatialSample`.
#' @param streaming Integer. CSR-to-CSC conversion mode. `0L` -> in-memory,
#' `1L` -> light streaming, `2L` -> heavy streaming. Defaults to `1L`.
#' @param raw_count_slot Where raw counts live. `"auto"` detects per file via
#' [bixverse::detect_raw_count_slot()]; otherwise one of `"X"`, `"raw.X"`,
#' `"layers.counts"`.
#' @param cell_id_col Optional string. Obs column holding the barcodes, for
#' files that carry no `_index`.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @section Coordinate order:
#' The column order of `obsm/spatial` is not fixed by the AnnData spec and the
#' two academic collections this reader was built against were checked file by
#' file rather than assumed. It gets resolved in this order:
#'
#' \enumerate{
#'   \item Whether the spots land on the tissue in the image the file ships.
#'   \item Whether only one order fits inside that image at all.
#'   \item The `obs` `pxl_row_in_fullres` / `pxl_col_in_fullres` labels. Last,
#'   because `scanpy.read_visium` swaps those two names relative to Space
#'   Ranger, so they mean the opposite in any file that went through it.
#'   \item `assume_orientation` out of [bixverse::params_sp_h5ad_io()].
#' }
#'
#' What was used is reported at load time and a warning is raised when nothing
#' in the file settled it. Note that a swap of `x` and `y` is a reflection and
#' preserves every pairwise distance, so the graph, Moran's I, SPARK-X and
#' neighbourhood enrichment are all numerically identical either way. The one
#' thing that cares is [bixverse::image_features_sp()].
#'
#' @section Graph layouts:
#' `"knn"` and `"radius"` always work. `"hex"` and `"square"` are exact
#' look-ups over `array_row`/`array_col` and need those columns in `obs`, which
#' `obsm/spatial` does not supply. The reader warns at load time when they are
#' missing rather than letting [bixverse::build_spatial_graph_sp()] fail later.
#'
#' @section One sample per object:
#' Like [bixverse::load_visium()], this rewrites the whole store, so it errors
#' when the object already carries a different sample. Re-running with the same
#' `exp_id` is fine.
#'
#' @return The class with updated shape information, a populated DuckDB and one
#' registered `SpatialSample`.
#'
#' @export
load_spatial_h5ad <- S7::new_generic(
  name = "load_spatial_h5ad",
  dispatch_args = "object",
  fun = function(
    object,
    h5_path,
    sp_io_param = params_sp_h5ad_io(),
    sc_qc_param = params_sc_min_quality(),
    exp_id,
    streaming = 1L,
    raw_count_slot = c("auto", "X", "raw.X", "layers.counts"),
    cell_id_col = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_spatial_h5ad SpatialSpot
#'
#' @export
S7::method(load_spatial_h5ad, SpatialSpot) <- function(
  object,
  h5_path,
  sp_io_param = params_sp_h5ad_io(),
  sc_qc_param = params_sc_min_quality(),
  exp_id,
  streaming = 1L,
  raw_count_slot = c("auto", "X", "raw.X", "layers.counts"),
  cell_id_col = NULL,
  .verbose = TRUE
) {
  raw_count_slot <- match.arg(raw_count_slot)

  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SpatialSpot))
  checkmate::assertFileExists(h5_path)
  assertSpH5adIo(sp_io_param)
  assertScMinQC(sc_qc_param)
  checkmate::qassert(exp_id, "S1")
  checkmate::qassert(streaming, "I1")
  checkmate::assertTRUE(streaming %in% c(0L, 1L, 2L))
  checkmate::assertChoice(
    raw_count_slot,
    c("auto", "X", "raw.X", "layers.counts")
  )
  checkmate::qassert(cell_id_col, c("0", "S1"))
  checkmate::qassert(.verbose, "B1")

  # Same guard as `load_visium()`: the ingest below rewrites counts, obs and
  # var while `add_spatial_sample()` merges, so a second `exp_id` would
  # register a sample whose data has just been overwritten.
  registered <- get_sample_ids(object)
  if (length(registered) > 0L && !identical(registered, exp_id)) {
    stop(sprintf(
      paste(
        "The object already holds sample(s) %s. `load_spatial_h5ad()`",
        "rewrites the whole store, so it cannot add '%s' on top."
      ),
      toString(sprintf("'%s'", registered)),
      exp_id
    ))
  }

  h5_path <- path.expand(h5_path)

  # The spatial extras go first so a file without coordinates fails before
  # anything has been written to disk.
  if (.verbose) {
    message(sprintf("Reading the spatial extras of '%s'.", exp_id))
  }
  spatial <- rs_sp_read_h5ad_spatial(
    h5_path = h5_path,
    library_id = sp_io_param$library_id,
    orientation = sp_io_param$assume_orientation
  )
  .sp_h5ad_report(spatial, exp_id = exp_id, .verbose = .verbose)

  all_ids <- .sp_h5ad_obs_ids(h5_path, cell_id_col = cell_id_col)

  object <- load_h5ad(
    object,
    h5_path = h5_path,
    sc_qc_param = sc_qc_param,
    streaming = streaming,
    raw_count_slot = raw_count_slot,
    cell_id_col = cell_id_col,
    .verbose = .verbose
  )

  duckdb_con <- get_sc_duckdb(object)
  obs_keys <- duckdb_con$get_obs_table(cols = c("cell_idx", "cell_id"))
  data.table::setorderv(obs_keys, "cell_idx")

  coord_rows <- .sp_h5ad_coord_rows(
    obs_ids = as.character(obs_keys[["cell_id"]]),
    all_ids = all_ids,
    n_coords = nrow(spatial$coords)
  )

  if (.verbose) {
    message(sprintf(
      "Attaching %i of %i spot coordinates.",
      length(coord_rows),
      nrow(spatial$coords)
    ))
  }

  # Written under the Space Ranger names because that is what the default
  # `coord_cols` of `new_spatial_sample()` looks for and what `load_visium()`
  # writes. A file that already carried them gets them replaced by the
  # orientation-resolved pair, which is the point.
  coord_dt <- data.table::data.table(
    cell_idx = as.integer(obs_keys[["cell_idx"]]),
    pxl_col_in_fullres = as.numeric(spatial$coords[coord_rows, 1L]),
    pxl_row_in_fullres = as.numeric(spatial$coords[coord_rows, 2L]),
    exp_id = exp_id
  )
  duckdb_con$join_data_obs(new_data = coord_dt)

  # one row, purely to learn the column names
  has_in_tissue <- "in_tissue" %in%
    names(duckdb_con$get_obs_table(indices = 1L))
  if (sp_io_param$in_tissue_only && has_in_tissue) {
    object <- .apply_in_tissue_filter(
      object = object,
      in_tissue_only = TRUE,
      .verbose = .verbose
    )
  } else if (sp_io_param$in_tissue_only && .verbose) {
    message(" No `in_tissue` column in the obs table, keeping every spot.")
  }

  if (.verbose) {
    message("Registering the spatial sample.")
  }

  kept_exp_ids <- duckdb_con$get_obs_table(
    cols = "exp_id",
    filtered = TRUE
  )[["exp_id"]]

  sample <- new_spatial_sample(
    exp_id = exp_id,
    technology = sp_io_param$technology,
    n_spots = sum(kept_exp_ids == exp_id),
    array_cols = if (isTRUE(spatial$has_array_indices)) {
      c(row = "array_row", col = "array_col")
    } else {
      character(0)
    },
    scale_factors = .sp_h5ad_scale_factors(spatial),
    image_paths = list()
  )

  object <- add_spatial_sample(object, sample)

  return(object)
}
