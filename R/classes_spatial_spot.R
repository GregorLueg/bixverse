# SpatialSpot S7 and sp_cache S3 -----------------------------------------------

# SpatialSpot inherits from SingleCells, adding two properties:
# - `samples`:  named list (keyed by exp_id) of SpatialSample objects.
# - `sp_cache`: a new S3 cache for per-sample spatial state. Sibling of
#               the inherited `sc_cache`, not a child class. This follows
#               the SingleCellsMultiModal pattern (separate cache per
#               modality / per concern) rather than overriding sc_cache.

# s3 ---------------------------------------------------------------------------

## sp_cache --------------------------------------------------------------------

#' Construct an empty SpCache
#'
#' @description
#' Helper that builds the per-sample spatial cache used by
#' [bixverse::SpatialSpot()]. Each slot is a named list keyed by the
#' experiment identifier (`exp_id`).
#'
#' @return An empty `SpCache` S3 object.
#'
#' @export
#'
#' @keywords internal
new_sp_cache <- function() {
  sp_cache <- list(
    per_sample_spatial_graph = list(),
    per_sample_pca = list(),
    per_sample_morans_i = list(),
    per_sample_sparkx = list(),
    per_sample_nhood_enrichment = list()
  )

  class(sp_cache) <- "SpCache"

  return(sp_cache)
}

## setters (SpCache) -----------------------------------------------------------

#' @rdname set_per_sample_spatial_graph
#'
#' @export
set_per_sample_spatial_graph.SpCache <- function(x, exp_id, graph, ...) {
  checkmate::assertClass(x, "SpCache")
  checkmate::qassert(exp_id, "S1")
  # Accept any CSC sparse matrix (`dgCMatrix`, `dsCMatrix`, etc.) - the
  # plan calls for CSR storage, but Matrix's CSC base class covers both
  # the general and symmetric variants returned by `Matrix::Matrix(..., sparse=TRUE)`.
  checkmate::assertClass(graph, "CsparseMatrix")

  x[["per_sample_spatial_graph"]][[exp_id]] <- graph

  return(x)
}

#' @rdname set_per_sample_pca
#'
#' @export
set_per_sample_pca.SpCache <- function(x, exp_id, pca, ...) {
  checkmate::assertClass(x, "SpCache")
  checkmate::qassert(exp_id, "S1")
  checkmate::assertList(pca, names = "named")
  checkmate::assertNames(
    names(pca),
    must.include = c("factors", "loadings", "singular_vals")
  )

  x[["per_sample_pca"]][[exp_id]] <- pca

  return(x)
}

#' @rdname set_per_sample_morans_i
#'
#' @export
set_per_sample_morans_i.SpCache <- function(x, exp_id, morans_i, ...) {
  checkmate::assertClass(x, "SpCache")
  checkmate::qassert(exp_id, "S1")
  checkmate::assertDataTable(morans_i)

  x[["per_sample_morans_i"]][[exp_id]] <- morans_i

  return(x)
}

#' @rdname set_per_sample_sparkx
#'
#' @export
set_per_sample_sparkx.SpCache <- function(x, exp_id, sparkx, ...) {
  checkmate::assertClass(x, "SpCache")
  checkmate::qassert(exp_id, "S1")
  checkmate::assertDataTable(sparkx)

  x[["per_sample_sparkx"]][[exp_id]] <- sparkx

  return(x)
}

#' @rdname set_per_sample_nhood_enrichment
#'
#' @export
set_per_sample_nhood_enrichment.SpCache <- function(
  x,
  exp_id,
  label_col,
  result,
  ...
) {
  checkmate::assertClass(x, "SpCache")
  checkmate::qassert(exp_id, "S1")
  checkmate::qassert(label_col, "S1")

  if (is.null(x[["per_sample_nhood_enrichment"]][[exp_id]])) {
    x[["per_sample_nhood_enrichment"]][[exp_id]] <- list()
  }
  x[["per_sample_nhood_enrichment"]][[exp_id]][[label_col]] <- result

  return(x)
}

## getters (SpCache) -----------------------------------------------------------

#' @rdname get_per_sample_spatial_graph
#'
#' @export
get_per_sample_spatial_graph.SpCache <- function(x, exp_id, ...) {
  checkmate::assertClass(x, "SpCache")
  checkmate::qassert(exp_id, "S1")

  return(x[["per_sample_spatial_graph"]][[exp_id]])
}

#' @rdname get_per_sample_pca
#'
#' @export
get_per_sample_pca.SpCache <- function(x, exp_id, ...) {
  checkmate::assertClass(x, "SpCache")
  checkmate::qassert(exp_id, "S1")

  return(x[["per_sample_pca"]][[exp_id]])
}

#' @rdname get_per_sample_morans_i
#'
#' @export
get_per_sample_morans_i.SpCache <- function(x, exp_id, ...) {
  checkmate::assertClass(x, "SpCache")
  checkmate::qassert(exp_id, "S1")

  return(x[["per_sample_morans_i"]][[exp_id]])
}

#' @rdname get_per_sample_sparkx
#'
#' @export
get_per_sample_sparkx.SpCache <- function(x, exp_id, ...) {
  checkmate::assertClass(x, "SpCache")
  checkmate::qassert(exp_id, "S1")

  return(x[["per_sample_sparkx"]][[exp_id]])
}

#' @rdname get_per_sample_nhood_enrichment
#'
#' @export
get_per_sample_nhood_enrichment.SpCache <- function(
  x,
  exp_id,
  label_col,
  ...
) {
  checkmate::assertClass(x, "SpCache")
  checkmate::qassert(exp_id, "S1")
  checkmate::qassert(label_col, "S1")

  return(x[["per_sample_nhood_enrichment"]][[exp_id]][[label_col]])
}

## primitives ------------------------------------------------------------------

#' Print method for SpCache objects
#'
#' @param x A `SpCache` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.SpCache <- function(x, ...) {
  ids <- unique(unlist(lapply(x, names)))
  cat(
    "SpCache (per-sample spatial cache)\n",
    sprintf("  Samples with cached state: %i\n", length(ids)),
    sprintf(
      "  IDs: %s\n",
      if (length(ids) > 0) paste(ids, collapse = ", ") else "none"
    ),
    sep = ""
  )
  invisible(x)
}

# s7 ---------------------------------------------------------------------------

## SpatialSpot ----------------------------------------------------------------

#' @title bixverse SpatialSpot class
#'
#' @description
#' Spatial transcriptomics class for spot-based assays (Visium, Xenium).
#' Inherits from [bixverse::SingleCells()] - the on-disk DuckDB + Rust
#' binary count storage and all global single-cell methods (HVG, PCA,
#' KNN, SNN, embeddings) work unchanged. Adds per-experiment spatial
#' metadata in `samples` and per-experiment spatial analysis results in
#' `sp_cache`.
#'
#' The obs table on disk is expected to include an `exp_id` column to
#' distinguish spots from different slides/sections.
#'
#' @param dir_data String. Path to the directory holding the SpatialSpot
#' on-disk files (DuckDB plus the two `counts_*.bin` Rust files).
#'
#' @section Properties:
#' \describe{
#'   \item{db_connection}{R6 [bixverse::SpatialSpotDuckDB()] instance.}
#'   \item{count_connection}{R6-like environment pointing to the Rust
#'   sparse count storage.}
#'   \item{dir_data}{Path to the on-disk data directory.}
#'   \item{sc_cache}{Inherited `ScCache` for global single-cell-style
#'   state (PCA, KNN, SNN, embeddings).}
#'   \item{sc_map}{Inherited `ScMap` with gene/cell mappings and HVG
#'   indices.}
#'   \item{dims}{Two-element integer with `c(no_spots, no_genes)`.}
#'   \item{samples}{Named list keyed by `exp_id` holding one
#'   [bixverse::new_spatial_sample()] per slide/section.}
#'   \item{sp_cache}{[bixverse::new_sp_cache()] holding per-sample
#'   spatial analysis results.}
#' }
#'
#' @return Returns the `SpatialSpot` class for further operations.
#'
#' @export
SpatialSpot <- S7::new_class(
  name = "SpatialSpot",
  parent = SingleCells,
  properties = list(
    db_connection = S7::class_any,
    count_connection = S7::class_any,
    dir_data = S7::class_character,
    sc_cache = S7::class_any,
    sc_map = S7::class_any,
    dims = S7::class_integer,
    samples = S7::class_list,
    sp_cache = S7::class_any
  ),
  constructor = function(dir_data) {
    # checks
    checkmate::assertDirectoryExists(dir_data)
    dir_data <- path.expand(dir_data)

    # Rust pointer for counts
    count_connection <- SingleCellCountData$new(
      f_path_cells = file.path(dir_data, "counts_cells.bin"),
      f_path_genes = file.path(dir_data, "counts_genes.bin")
    )

    # DuckDB connector (spatial subclass)
    db_connection <- SpatialSpotDuckDB$new(db_dir = dir_data)

    S7::new_object(
      S7::S7_object(),
      db_connection = db_connection,
      count_connection = count_connection,
      dir_data = dir_data,
      sc_cache = new_sc_cache(),
      sc_map = new_sc_mapper(),
      dims = c(0L, 0L),
      samples = list(),
      sp_cache = new_sp_cache()
    )
  }
)

## getters ---------------------------------------------------------------------

### samples --------------------------------------------------------------------

#' @method get_samples SpatialSpot
#'
#' @export
S7::method(get_samples, SpatialSpot) <- function(object) {
  checkmate::assertTRUE(S7::S7_inherits(object, SpatialSpot))
  return(S7::prop(object, "samples"))
}

#' @method get_sample SpatialSpot
#'
#' @export
S7::method(get_sample, SpatialSpot) <- function(object, exp_id) {
  checkmate::assertTRUE(S7::S7_inherits(object, SpatialSpot))
  checkmate::qassert(exp_id, "S1")

  samples <- S7::prop(object, "samples")
  if (!exp_id %in% names(samples)) {
    stop(sprintf("exp_id '%s' not found in samples.", exp_id))
  }
  return(samples[[exp_id]])
}

#' @method get_sample_ids SpatialSpot
#'
#' @export
S7::method(get_sample_ids, SpatialSpot) <- function(object) {
  checkmate::assertTRUE(S7::S7_inherits(object, SpatialSpot))
  return(names(S7::prop(object, "samples")))
}

### sp_cache -------------------------------------------------------------------

#' @method get_sp_cache SpatialSpot
#'
#' @export
S7::method(get_sp_cache, SpatialSpot) <- function(object) {
  checkmate::assertTRUE(S7::S7_inherits(object, SpatialSpot))
  return(S7::prop(object, "sp_cache"))
}

### spot indices --------------------------------------------------------------

#' @method get_spot_indices_for_exp SpatialSpot
#'
#' @export
S7::method(get_spot_indices_for_exp, SpatialSpot) <- function(
  object,
  exp_id,
  filtered = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SpatialSpot))
  checkmate::qassert(exp_id, "S1")
  checkmate::qassert(filtered, "B1")

  duckdb_con <- get_sc_duckdb(object)
  return(duckdb_con$get_spot_indices_for_exp(
    exp_id = exp_id,
    filtered = filtered
  ))
}

### spatial coords -------------------------------------------------------------

#' @method get_spatial_coords SpatialSpot
#'
#' @export
S7::method(get_spatial_coords, SpatialSpot) <- function(
  object,
  exp_id,
  filtered = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SpatialSpot))
  checkmate::qassert(exp_id, "S1")
  checkmate::qassert(filtered, "B1")

  sample <- get_sample(object, exp_id)
  coord_cols <- sample$coord_cols

  duckdb_con <- get_sc_duckdb(object)
  # `get_spot_indices_for_exp` returns 0-based indices in original obs
  # space; DuckDB obs uses 1-based `cell_idx`.
  indices_rust <- duckdb_con$get_spot_indices_for_exp(
    exp_id = exp_id,
    filtered = filtered
  )
  indices_r <- as.integer(indices_rust + 1L)

  obs_sub <- duckdb_con$get_obs_table(
    indices = indices_r,
    cols = c(coord_cols[["x"]], coord_cols[["y"]]),
    filtered = FALSE
  )

  coords <- cbind(
    x = as.numeric(obs_sub[[coord_cols[["x"]]]]),
    y = as.numeric(obs_sub[[coord_cols[["y"]]]])
  )

  return(coords)
}

### images ---------------------------------------------------------------------

#' @method get_image SpatialSpot
#'
#' @export
S7::method(get_image, SpatialSpot) <- function(
  object,
  exp_id,
  resolution = c("lowres", "hires", "fullres")
) {
  resolution <- match.arg(resolution)

  checkmate::assertTRUE(S7::S7_inherits(object, SpatialSpot))
  checkmate::qassert(exp_id, "S1")

  sample <- get_sample(object, exp_id)
  rel_path <- sample$image_paths[[resolution]]
  if (is.null(rel_path)) {
    stop(sprintf(
      "No image at resolution '%s' registered for exp_id '%s'.",
      resolution,
      exp_id
    ))
  }
  abs_path <- file.path(S7::prop(object, "dir_data"), rel_path)
  checkmate::assertFileExists(abs_path)

  ext <- tolower(tools::file_ext(abs_path))
  img <- switch(
    ext,
    png = {
      if (!requireNamespace("png", quietly = TRUE)) {
        stop("Reading PNG images requires the 'png' package.")
      }
      png::readPNG(abs_path)
    },
    jpg = ,
    jpeg = {
      if (!requireNamespace("jpeg", quietly = TRUE)) {
        stop("Reading JPEG images requires the 'jpeg' package.")
      }
      jpeg::readJPEG(abs_path)
    },
    tif = ,
    tiff = stop(
      "TIFF images are not yet supported. Provide a PNG or JPEG."
    ),
    stop(sprintf("Unsupported image extension '%s'.", ext))
  )

  return(img)
}

## sp_cache forwarding (SpatialSpot -> SpCache) --------------------------------

#' @rdname set_per_sample_spatial_graph
#'
#' @export
set_per_sample_spatial_graph.SpatialSpot <- function(
  x,
  exp_id,
  graph,
  ...
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpot))
  S7::prop(x, "sp_cache") <- set_per_sample_spatial_graph(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id,
    graph = graph
  )
  return(x)
}

#' @rdname set_per_sample_pca
#'
#' @export
set_per_sample_pca.SpatialSpot <- function(x, exp_id, pca, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpot))
  S7::prop(x, "sp_cache") <- set_per_sample_pca(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id,
    pca = pca
  )
  return(x)
}

#' @rdname set_per_sample_morans_i
#'
#' @export
set_per_sample_morans_i.SpatialSpot <- function(x, exp_id, morans_i, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpot))
  S7::prop(x, "sp_cache") <- set_per_sample_morans_i(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id,
    morans_i = morans_i
  )
  return(x)
}

#' @rdname set_per_sample_sparkx
#'
#' @export
set_per_sample_sparkx.SpatialSpot <- function(x, exp_id, sparkx, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpot))
  S7::prop(x, "sp_cache") <- set_per_sample_sparkx(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id,
    sparkx = sparkx
  )
  return(x)
}

#' @rdname set_per_sample_nhood_enrichment
#'
#' @export
set_per_sample_nhood_enrichment.SpatialSpot <- function(
  x,
  exp_id,
  label_col,
  result,
  ...
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpot))
  S7::prop(x, "sp_cache") <- set_per_sample_nhood_enrichment(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id,
    label_col = label_col,
    result = result
  )
  return(x)
}

#' @rdname get_per_sample_spatial_graph
#'
#' @export
get_per_sample_spatial_graph.SpatialSpot <- function(x, exp_id, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpot))
  return(get_per_sample_spatial_graph(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id
  ))
}

#' @rdname get_per_sample_pca
#'
#' @export
get_per_sample_pca.SpatialSpot <- function(x, exp_id, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpot))
  return(get_per_sample_pca(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id
  ))
}

#' @rdname get_per_sample_morans_i
#'
#' @export
get_per_sample_morans_i.SpatialSpot <- function(x, exp_id, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpot))
  return(get_per_sample_morans_i(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id
  ))
}

#' @rdname get_per_sample_sparkx
#'
#' @export
get_per_sample_sparkx.SpatialSpot <- function(x, exp_id, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpot))
  return(get_per_sample_sparkx(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id
  ))
}

#' @rdname get_per_sample_nhood_enrichment
#'
#' @export
get_per_sample_nhood_enrichment.SpatialSpot <- function(
  x,
  exp_id,
  label_col,
  ...
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpot))
  return(get_per_sample_nhood_enrichment(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id,
    label_col = label_col
  ))
}

## sample registration --------------------------------------------------------

#' @method add_spatial_sample SpatialSpot
#'
#' @export
S7::method(add_spatial_sample, SpatialSpot) <- function(object, sample) {
  checkmate::assertTRUE(S7::S7_inherits(object, SpatialSpot))
  assertSpatialSample(sample)

  samples <- S7::prop(object, "samples")
  samples[[sample$exp_id]] <- sample
  S7::prop(object, "samples") <- samples

  return(object)
}

## primitives ------------------------------------------------------------------

#' @name print.SpatialSpot
#'
#' @title print Method for SpatialSpot object
#'
#' @description
#' Print method that summarises the inherited single-cell-style state
#' plus a per-experiment breakdown of spatial state.
#'
#' @param x A `SpatialSpot` object.
#' @param ... Unused.
#'
#' @returns Invisibly returns `x`.
#'
#' @method print SpatialSpot
#'
#' @keywords internal
S7::method(print, SpatialSpot) <- function(x, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpot))

  dims <- S7::prop(x, "dims")
  no_cells_to_keep <- length(get_cells_to_keep(x))
  sc_map <- S7::prop(x, "sc_map")
  sc_cache <- S7::prop(x, "sc_cache")
  sp_cache <- S7::prop(x, "sp_cache")
  samples <- S7::prop(x, "samples")

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
    "Spatial transcriptomics experiment (SpatialSpot).\n",
    sprintf("  No spots (original): %i\n", dims[1]),
    sprintf("   To keep n: %i\n", no_cells_to_keep),
    sprintf("  No genes: %i\n", dims[2]),
    sprintf("  HVG calculated: %s\n", hvg_calculated),
    sprintf("  PCA (global) calculated: %s\n", pca_calculated),
    sprintf("  Other embeddings: %s\n", other_embeddings_str),
    sprintf("  KNN generated: %s\n", knn_generated),
    sprintf("  SNN generated: %s\n", snn_generated),
    sprintf("  No experiments: %i\n", length(samples)),
    sep = ""
  )

  if (length(samples) > 0) {
    for (id in names(samples)) {
      s <- samples[[id]]
      has_graph <- !is.null(sp_cache$per_sample_spatial_graph[[id]])
      has_pca <- !is.null(sp_cache$per_sample_pca[[id]])
      n_imgs <- length(s$image_paths)
      cat(sprintf(
        paste0(
          "   - %s (%s, %i spots, %i image(s), ",
          "graph: %s, per-sample PCA: %s)\n"
        ),
        id,
        s$technology,
        s$n_spots,
        n_imgs,
        has_graph,
        has_pca
      ))
    }
  }

  invisible(x)
}
