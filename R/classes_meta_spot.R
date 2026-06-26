# MetaSpot S7 ------------------------------------------------------------------

# MetaSpot inherits from MetaCells. Adds spatial structure at the
# metaspot level: per-sample centroids and per-sample CSR adjacency
# graphs at metaspot resolution, plus the SpatialSample metadata copied
# from the parent SpatialSpot.
#
# Unlike SpatialSpot, these aren't held in a cache: metaspot identities
# and positions only exist because aggregation already ran, so they are
# core data and live as first-class properties.

# s7 ---------------------------------------------------------------------------

#' @title bixverse MetaSpot class
#'
#' @description
#' Metaspot (aggregated SpatialSpot) class. Inherits from
#' [bixverse::MetaCells()] - all in-memory obs/var/data primitives work
#' unchanged. Adds per-sample centroids and adjacency graphs at metaspot
#' resolution.
#'
#' The `obs_table` must contain an `exp_id` column identifying which
#' slide each metaspot belongs to.
#'
#' @param meta_cell_data Named list. Output of meta-spot generation Rust
#' functions. Must contain `aggregated` and `assignments`, same shape as
#' [bixverse::MetaCells()].
#' @param var_data data.table with the variable/feature information.
#' Must contain `gene_idx` and `gene_id`.
#' @param samples Named list keyed by `exp_id` of
#' [bixverse::new_spatial_sample()] objects (typically copied from the
#' parent SpatialSpot).
#' @param per_sample_coords Named list keyed by `exp_id` of two-column
#' numeric matrices with the metaspot centroids for each sample.
#' @param per_sample_graph Named list keyed by `exp_id` of `dgCMatrix`
#' (CSR-style) adjacency matrices at metaspot resolution. Optional;
#' defaults to an empty list.
#' @param aggregation_method String describing how metaspots were
#' generated. E.g. `"hex_bin_16um"`, `"superspot_walktrap"`.
#'
#' @section Properties:
#' \describe{
#'   \item{obs_table}{The metaspot observation table. Must include
#'   `exp_id`.}
#'   \item{var_table}{The metaspot variable table.}
#'   \item{data}{Named list with `raw` and `norm` aggregated count
#'   matrices.}
#'   \item{sc_cache}{Inherited `ScCache` for global analysis state.}
#'   \item{original_assignment}{Inherited assignment metadata.}
#'   \item{dims}{Two-element integer with `c(no_metaspots, no_genes)`.}
#'   \item{other_data}{Inherited list for additional data.}
#'   \item{meta_cell_method}{Inherited origin string.}
#'   \item{samples}{Named list of [bixverse::new_spatial_sample()]
#'   keyed by `exp_id`.}
#'   \item{per_sample_coords}{Named list of n_metaspot x 2 coord
#'   matrices keyed by `exp_id`.}
#'   \item{per_sample_graph}{Named list of `dgCMatrix` adjacency
#'   matrices at metaspot resolution, keyed by `exp_id`.}
#'   \item{aggregation_method}{String identifying the aggregation
#'   algorithm used.}
#' }
#'
#' @return Returns the `MetaSpot` class for further operations.
#'
#' @export
MetaSpot <- S7::new_class(
  name = "MetaSpot",
  parent = MetaCells,
  properties = list(
    obs_table = S7::class_data.frame,
    var_table = S7::class_data.frame,
    data = S7::class_list,
    sc_cache = S7::class_any,
    original_assignment = S7::class_list,
    dims = S7::class_integer,
    other_data = S7::class_list,
    meta_cell_method = S7::class_character,
    samples = S7::class_list,
    per_sample_coords = S7::class_list,
    per_sample_graph = S7::class_list,
    aggregation_method = S7::class_character
  ),
  constructor = function(
    meta_cell_data,
    var_data,
    samples,
    per_sample_coords,
    per_sample_graph = list(),
    aggregation_method
  ) {
    # checks
    checkmate::assertList(meta_cell_data)
    checkmate::assertNames(
      names(meta_cell_data),
      must.include = c("aggregated", "assignments")
    )
    checkmate::assertDataTable(var_data)
    checkmate::assertNames(
      names(var_data),
      must.include = c("gene_idx", "gene_id")
    )
    checkmate::assertList(samples, names = "named")
    for (s in samples) assertSpatialSample(s)
    checkmate::assertList(per_sample_coords, names = "named")
    checkmate::assertList(per_sample_graph, names = "named")
    checkmate::qassert(aggregation_method, "S1")

    # build the metaspot-level obs as MetaCells does, then attach exp_id
    n_digits <- nchar(as.character(meta_cell_data$assignments$n_metacells))
    format_str <- sprintf("meta_spot_%%0%dd", n_digits)

    obs_data <- data.table::data.table(
      meta_cell_idx = 1:meta_cell_data$assignments$n_metacells,
      meta_cell_id = sprintf(
        format_str,
        1:meta_cell_data$assignments$n_metacells
      ),
      no_originating_cells = purrr::map_dbl(
        meta_cell_data$assignments$metacells,
        length
      ),
      original_cell_idx = meta_cell_data$assignments$metacells
    )

    if (!is.null(meta_cell_data$assignments$exp_id)) {
      obs_data[, exp_id := meta_cell_data$assignments$exp_id]
    } else {
      stop(
        paste(
          "MetaSpot obs requires an 'exp_id' column.",
          "meta_cell_data$assignments$exp_id was NULL."
        )
      )
    }

    c(raw_counts, norm_counts) %<-%
      get_meta_cell_matrices(meta_cell_data$aggregated)

    rownames(raw_counts) <- rownames(norm_counts) <- obs_data$meta_cell_id
    colnames(raw_counts) <- colnames(norm_counts) <- var_data$gene_id

    # consistency: every per_sample_coords key must be a known exp_id
    sample_ids <- names(samples)
    bad_coord_keys <- setdiff(names(per_sample_coords), sample_ids)
    if (length(bad_coord_keys) > 0) {
      stop(sprintf(
        "per_sample_coords has keys not present in samples: %s",
        paste(bad_coord_keys, collapse = ", ")
      ))
    }
    bad_graph_keys <- setdiff(names(per_sample_graph), sample_ids)
    if (length(bad_graph_keys) > 0) {
      stop(sprintf(
        "per_sample_graph has keys not present in samples: %s",
        paste(bad_graph_keys, collapse = ", ")
      ))
    }

    for (id in names(per_sample_coords)) {
      m <- per_sample_coords[[id]]
      checkmate::assertMatrix(m, mode = "numeric", ncols = 2)
    }
    for (id in names(per_sample_graph)) {
      checkmate::assertClass(per_sample_graph[[id]], "CsparseMatrix")
    }

    S7::new_object(
      S7::S7_object(),
      obs_table = obs_data,
      var_table = var_data,
      data = list(raw = raw_counts, norm = norm_counts),
      sc_cache = new_sc_cache(),
      original_assignment = meta_cell_data$assignments[c(
        "assignments",
        "n_cells",
        "n_metacells",
        "n_unassigned"
      )],
      dims = as.integer(c(nrow(raw_counts), ncol(norm_counts))),
      other_data = list(),
      meta_cell_method = aggregation_method,
      samples = samples,
      per_sample_coords = per_sample_coords,
      per_sample_graph = per_sample_graph,
      aggregation_method = aggregation_method
    )
  }
)

## getters --------------------------------------------------------------------

#' @method get_samples MetaSpot
#'
#' @export
S7::method(get_samples, MetaSpot) <- function(object) {
  checkmate::assertTRUE(S7::S7_inherits(object, MetaSpot))
  return(S7::prop(object, "samples"))
}

#' @method get_sample_ids MetaSpot
#'
#' @export
S7::method(get_sample_ids, MetaSpot) <- function(object) {
  checkmate::assertTRUE(S7::S7_inherits(object, MetaSpot))
  return(names(S7::prop(object, "samples")))
}

#' @method get_spatial_coords MetaSpot
#'
#' @export
S7::method(get_spatial_coords, MetaSpot) <- function(
  object,
  exp_id,
  filtered = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, MetaSpot))
  checkmate::qassert(exp_id, "S1")
  # `filtered` accepted for signature parity with SpatialSpot, but
  # MetaSpot has no to_keep filter at this stage - all metaspots are
  # kept.
  checkmate::qassert(filtered, "B1")

  coords <- S7::prop(object, "per_sample_coords")
  if (!exp_id %in% names(coords)) {
    stop(sprintf("No coords found for exp_id '%s'.", exp_id))
  }
  return(coords[[exp_id]])
}

#' @method get_per_sample_graph MetaSpot
#'
#' @export
S7::method(get_per_sample_graph, MetaSpot) <- function(object, exp_id) {
  checkmate::assertTRUE(S7::S7_inherits(object, MetaSpot))
  checkmate::qassert(exp_id, "S1")

  graphs <- S7::prop(object, "per_sample_graph")
  if (!exp_id %in% names(graphs)) {
    stop(sprintf("No graph found for exp_id '%s'.", exp_id))
  }
  return(graphs[[exp_id]])
}

## primitives ------------------------------------------------------------------

#' @name print.MetaSpot
#'
#' @title print Method for MetaSpot object
#'
#' @description
#' Prints the inherited MetaCells summary plus a per-sample spatial
#' breakdown.
#'
#' @param x A `MetaSpot` object.
#' @param ... Unused.
#'
#' @returns Invisibly returns `x`.
#'
#' @method print MetaSpot
#'
#' @keywords internal
S7::method(print, MetaSpot) <- function(x, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, MetaSpot))

  dims <- S7::prop(x, "dims")
  original_assignment <- S7::prop(x, "original_assignment")
  sc_cache <- S7::prop(x, "sc_cache")
  samples <- S7::prop(x, "samples")
  per_sample_graph <- S7::prop(x, "per_sample_graph")
  aggregation_method <- S7::prop(x, "aggregation_method")

  pca_calculated <- !is.null(sc_cache[["pca_factors"]])
  knn_generated <- !is.null(sc_cache[["knn"]])
  snn_generated <- !is.null(sc_cache[["snn_graph"]])

  cat(
    "Spatial transcriptomics experiment (MetaSpot).\n",
    sprintf("  Aggregation method: %s\n", aggregation_method),
    sprintf("  No metaspots: %i\n", dims[1]),
    sprintf("  No genes: %i\n", dims[2]),
    sprintf("  No original spots: %i\n", original_assignment$n_cells),
    sprintf(
      "  No unassigned spots: %i\n",
      original_assignment$n_unassigned
    ),
    sprintf("  PCA calculated: %s\n", pca_calculated),
    sprintf("  KNN generated: %s\n", knn_generated),
    sprintf("  SNN generated: %s\n", snn_generated),
    sprintf("  No experiments: %i\n", length(samples)),
    sep = ""
  )

  if (length(samples) > 0) {
    for (id in names(samples)) {
      s <- samples[[id]]
      has_graph <- !is.null(per_sample_graph[[id]])
      cat(sprintf(
        "   - %s (%s, %i original spots, graph: %s)\n",
        id,
        s$technology,
        s$n_spots,
        has_graph
      ))
    }
  }

  invisible(x)
}
