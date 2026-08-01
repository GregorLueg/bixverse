# spatial analysis methods -----------------------------------------------------

# Every method here is registered against the SpOrSpSubset union, so a
# SpatialSpot and a SpatialSpotSubset behave identically. The only difference is
# how `exp_id` resolves: a subset carries exactly one.
#
# Index spaces, in the order they cross the boundary:
#
# | thing                       | base | space                               |
# |-----------------------------|------|-------------------------------------|
# | DuckDB `cell_idx`           | 1    | global                              |
# | `get_spot_indices_for_exp()`| 0    | global                              |
# | `spots_to_keep` into Rust   | 0    | global                              |
# | `neighbours` / `weights`    | 0    | local, positional vs spots_to_keep  |
# | `labels`                    | 0    | local                               |
# | coordinates                 | n/a  | full-res pixels, (x, y)             |

## helpers ---------------------------------------------------------------------

### adjacency <> sparse matrix -------------------------------------------------

#' Pack a spatial adjacency pair into a sparse matrix
#'
#' @description
#' The `SpCache` stores the spatial graph as a sparse matrix because that is
#' what a user wants to poke at. Entry `(i, j)` is the weight of the directed
#' edge from spot `i` to spot `j`, exactly as [bixverse::rs_sp_build_graph()]
#' returned it, i.e. before any symmetrisation. Every downstream statistic
#' symmetrises internally.
#'
#' Weights of exactly zero are not representable this way. No layout produces
#' them (binary weights are one, distance-based weights are positive), but a
#' gaussian bandwidth small enough to underflow would silently lose edges.
#'
#' @param indices List of integer vectors. Local 0-based neighbour indices.
#' @param weights List of numeric vectors. Aligned with `indices`.
#' @param n_spots Integer. Number of spots.
#'
#' @return A `dgCMatrix` of dimension `n_spots x n_spots`.
#'
#' @keywords internal
.sp_adjacency_to_sparse <- function(indices, weights, n_spots) {
  checkmate::assertList(indices)
  checkmate::assertList(weights)
  checkmate::qassert(n_spots, "I1[1,)")
  checkmate::assertTRUE(length(indices) == n_spots)

  n_edges <- lengths(indices)
  if (sum(n_edges) == 0L) {
    stop("The spatial graph has no edges.")
  }

  Matrix::sparseMatrix(
    i = rep.int(seq_len(n_spots), n_edges),
    j = as.integer(unlist(indices, use.names = FALSE)) + 1L,
    x = as.numeric(unlist(weights, use.names = FALSE)),
    dims = c(n_spots, n_spots),
    repr = "C"
  )
}

#' Unpack a sparse matrix back into a spatial adjacency pair
#'
#' @description
#' Inverse of [.sp_adjacency_to_sparse()]. Goes through the transpose so the
#' row slices are contiguous in the CSC layout rather than a per-row scan.
#'
#' @param graph A `dgCMatrix` as stored in the `SpCache`.
#'
#' @return A list with `indices` (local 0-based) and `weights`.
#'
#' @keywords internal
.sp_sparse_to_adjacency <- function(graph) {
  checkmate::assertClass(graph, "CsparseMatrix")

  # CSC of t(graph) is CSR of graph
  transposed <- Matrix::t(graph)
  ptr <- transposed@p
  n_spots <- ncol(transposed)

  slices <- lapply(seq_len(n_spots), \(i) {
    if (ptr[i] == ptr[i + 1L]) integer(0) else seq.int(ptr[i] + 1L, ptr[i + 1L])
  })

  list(
    indices = lapply(slices, \(s) as.integer(transposed@i[s])),
    weights = lapply(slices, \(s) as.numeric(transposed@x[s]))
  )
}

### resolution -----------------------------------------------------------------

#' Assert the object is one of the two spatial classes
#'
#' @description
#' `S7::S7_inherits()` does not take a union, so the two members get checked
#' one at a time.
#'
#' @param object The object to check.
#'
#' @return Invisibly `TRUE`, otherwise it errors.
#'
#' @keywords internal
.assert_sp_object <- function(object) {
  if (
    !S7::S7_inherits(object, SpatialSpot) &&
      !S7::S7_inherits(object, SpatialSpotSubset)
  ) {
    stop("The object must be a `SpatialSpot` or a `SpatialSpotSubset`.")
  }
  invisible(TRUE)
}

#' Resolve which experiments a spatial method should run over
#'
#' @param object A [bixverse::SpatialSpot()] or
#' [bixverse::SpatialSpotSubset()].
#' @param exp_id Optional string. `NULL` means every registered sample.
#'
#' @return Character vector of experiment identifiers.
#'
#' @keywords internal
.sp_resolve_exp_ids <- function(object, exp_id) {
  checkmate::qassert(exp_id, c("0", "S1"))

  available <- get_sample_ids(object)
  if (length(available) == 0L) {
    stop("No spatial sample registered on the object.")
  }
  if (is.null(exp_id)) {
    return(available)
  }
  checkmate::assertChoice(exp_id, available)

  return(exp_id)
}

#' Resolve a gene selection to on-disk 0-based gene indices
#'
#' @param object A [bixverse::SpatialSpot()] or
#' [bixverse::SpatialSpotSubset()].
#' @param genes `NULL` for every gene, a character vector of gene identifiers,
#' or an integer vector of 1-indexed positions in the var table.
#'
#' @return A list with `gene_indices` (0-based, for Rust) and `gene_ids`.
#'
#' @keywords internal
.sp_resolve_genes <- function(object, genes) {
  gene_names <- get_gene_names(object)

  idx <- if (is.null(genes)) {
    seq_along(gene_names)
  } else if (checkmate::qtest(genes, "S+")) {
    get_gene_indices(x = object, gene_ids = genes, rust_index = FALSE)
  } else if (checkmate::qtest(genes, "I+")) {
    checkmate::assertIntegerish(
      genes,
      lower = 1L,
      upper = length(gene_names),
      any.missing = FALSE
    )
    as.integer(genes)
  } else {
    stop("`genes` must be NULL, a character vector or an integer vector.")
  }

  list(
    gene_indices = as.integer(idx - 1L),
    gene_ids = gene_names[idx]
  )
}

#' Resolve the graph parameters for one sample
#'
#' @description
#' The lattice layouts need the array coordinates, which live in the obs table
#' rather than in the parameter list. They get injected here so
#' [bixverse::params_sp_graph()] stays a pure parameter bundle.
#'
#' @param graph_params List. Output of [bixverse::params_sp_graph()].
#' @param obs_dt `data.table`. Obs rows of this sample, in spot order.
#' @param sample A [bixverse::new_spatial_sample()].
#'
#' @return The parameter list Rust expects.
#'
#' @keywords internal
.sp_graph_params_for_sample <- function(graph_params, obs_dt, sample) {
  if (!graph_params$layout %in% c("hex", "square")) {
    return(graph_params)
  }

  array_cols <- sample$array_cols
  if (length(array_cols) != 2L) {
    stop(sprintf(
      "The '%s' layout needs array coordinates, which sample '%s' has none.",
      graph_params$layout,
      sample$exp_id
    ))
  }
  checkmate::assertNames(
    colnames(obs_dt),
    must.include = unname(array_cols)
  )

  graph_params$array_row <- as.integer(obs_dt[[array_cols[["row"]]]])
  graph_params$array_col <- as.integer(obs_dt[[array_cols[["col"]]]])

  return(graph_params)
}

#' Pull the obs rows of one experiment in spot order
#'
#' @description
#' Ascending `cell_idx`, matching what [bixverse::get_spot_indices_for_exp()]
#' and [bixverse::get_spatial_coords()] return.
#'
#' @param object A [bixverse::SpatialSpot()] or
#' [bixverse::SpatialSpotSubset()].
#' @param exp_id String. The experiment identifier.
#' @param cols Optional character vector of columns to pull.
#'
#' @return A `data.table` sorted by `cell_idx`.
#'
#' @keywords internal
.sp_obs_for_exp <- function(object, exp_id, cols = NULL) {
  if (!is.null(cols)) {
    cols <- unique(c("cell_idx", "exp_id", cols))
  }
  obs <- get_sc_obs(object, cols = cols, filtered = TRUE)
  # The mask is built outside the `[` call on purpose. The whole `i` expression
  # is evaluated with the columns in scope, and `obs` has a column named
  # `exp_id` while this function has an argument of the same name. Inside `i`
  # the column wins on both sides, so `obs[get("exp_id") == exp_id, ]` compares
  # the column with itself and quietly keeps every sample.
  keep <- obs[["exp_id"]] == exp_id
  obs <- obs[keep, ]
  data.table::setorderv(obs, "cell_idx")

  return(obs)
}

#' Fetch the cached spatial graph or fail with a pointed message
#'
#' @description
#' The graph is stored in the local index space: row `i` is the `i`-th spot of
#' `spot_idx` as it stood when the graph was built. Nothing in the matrix
#' records which spots those were, so a `to_keep` change to a different set of
#' the same size would pair expression, labels and adjacency against three
#' different spots and still return finite, plausible numbers. The hash stored
#' next to the graph is what turns that into an error.
#'
#' @param object A [bixverse::SpatialSpot()] or
#' [bixverse::SpatialSpotSubset()].
#' @param exp_id String. The experiment identifier.
#' @param spot_idx Integer vector. The 0-based global spot indices the caller
#' is about to run against.
#'
#' @return A list with `indices` and `weights`.
#'
#' @keywords internal
.sp_require_graph <- function(object, exp_id, spot_idx) {
  checkmate::qassert(spot_idx, "I+")

  graph <- get_per_sample_spatial_graph(object, exp_id = exp_id)
  if (is.null(graph)) {
    stop(sprintf(
      "No spatial graph for exp_id '%s'. Run build_spatial_graph_sp() first.",
      exp_id
    ))
  }
  if (nrow(graph) != length(spot_idx)) {
    stop(sprintf(
      "The cached graph of '%s' has %i spots, the object %i.",
      exp_id,
      nrow(graph),
      length(spot_idx)
    ))
  }

  stored_key <- get_sp_cache(object)[["per_sample_graph_spots"]][[exp_id]]
  if (is.null(stored_key)) {
    warning(sprintf(
      paste(
        "The cached graph of '%s' records no spot identity, so only its spot",
        "count can be checked. Rebuild it with build_spatial_graph_sp()."
      ),
      exp_id
    ))
  } else if (!identical(stored_key, rlang::hash(as.integer(spot_idx)))) {
    stop(sprintf(
      paste(
        "The cached graph of '%s' was built over a different set of spots.",
        "Re-run build_spatial_graph_sp() after changing `to_keep`."
      ),
      exp_id
    ))
  }

  return(.sp_sparse_to_adjacency(graph))
}

#' Map a registered image resolution onto its scale factor
#'
#' @param sample A [bixverse::new_spatial_sample()].
#' @param resolution String. One of `c("lowres", "hires", "cytassist",
#' "fullres")`.
#'
#' @return Float. The scale from full-resolution pixels into that image.
#'
#' @keywords internal
.sp_image_scalef <- function(sample, resolution) {
  if (resolution == "fullres") {
    return(1.0)
  }

  key <- switch(
    resolution,
    lowres = "tissue_lowres_scalef",
    hires = "tissue_hires_scalef",
    cytassist = "cytassist_scalef"
  )
  scalef <- sample$scale_factors[[key]]

  if (is.null(scalef) || !is.finite(scalef) || scalef <= 0) {
    stop(sprintf(
      "Sample '%s' has no usable '%s' for resolution '%s'.",
      sample$exp_id,
      key,
      resolution
    ))
  }

  return(as.numeric(scalef))
}

## methods ---------------------------------------------------------------------

### spatial graph --------------------------------------------------------------

S7::method(build_spatial_graph_sp, SpOrSpSubset) <- function(
  object,
  exp_id = NULL,
  graph_params = params_sp_graph(),
  knn_params = NULL,
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  .assert_sp_object(object)
  assertSpGraph(graph_params)
  if (!is.null(knn_params)) {
    assertScKnn(knn_params)
  }
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  exp_ids <- .sp_resolve_exp_ids(object, exp_id)
  verbosity <- parse_verbosity(.verbose)

  for (id in exp_ids) {
    sample <- get_sample(object, id)
    coords <- get_spatial_coords(object, exp_id = id, filtered = TRUE)
    obs <- .sp_obs_for_exp(object, id)
    spots_to_keep <- get_spot_indices_for_exp(
      object,
      exp_id = id,
      filtered = TRUE
    )

    if (nrow(obs) != nrow(coords) || length(spots_to_keep) != nrow(coords)) {
      stop(sprintf(
        "Got %i obs rows and %i spot indices for %i coordinates of '%s'.",
        nrow(obs),
        length(spots_to_keep),
        nrow(coords),
        id
      ))
    }

    rust_params <- .sp_graph_params_for_sample(graph_params, obs, sample)

    if (verbosity > 0L) {
      message(sprintf(
        "Building the '%s' spatial graph for '%s' (%i spots).",
        graph_params$layout,
        id,
        nrow(coords)
      ))
    }

    res <- rs_sp_build_graph(
      coords = coords,
      graph_params = rust_params,
      knn_params = knn_params,
      seed = seed,
      verbose = verbosity
    )

    if (verbosity > 0L) {
      # `degree` is the row sum of the symmetrised graph W + t(W), not the
      # neighbour count. On a binary lattice it is exactly twice it.
      message(sprintf(
        paste(
          "  %i directed edges, median %.3g neighbours per spot,",
          "median symmetrised degree %.3g."
        ),
        sum(lengths(res$indices)),
        stats::median(lengths(res$indices)),
        stats::median(res$degree)
      ))
    }

    object <- set_per_sample_spatial_graph(
      x = object,
      exp_id = id,
      graph = .sp_adjacency_to_sparse(
        indices = res$indices,
        weights = res$weights,
        n_spots = nrow(coords)
      ),
      spot_idx = as.integer(spots_to_keep)
    )
  }

  return(object)
}

### moran's i ------------------------------------------------------------------

S7::method(morans_i_sp, SpOrSpSubset) <- function(
  object,
  exp_id = NULL,
  genes = NULL,
  svg_params = params_sp_svg(),
  streaming = TRUE,
  .verbose = TRUE
) {
  # checks
  .assert_sp_object(object)
  assertSpSvg(svg_params)
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  exp_ids <- .sp_resolve_exp_ids(object, exp_id)
  verbosity <- parse_verbosity(.verbose)
  gene_sel <- .sp_resolve_genes(object, genes)
  f_path_gene <- get_rust_count_gene_f_path(object)

  for (id in exp_ids) {
    spots_to_keep <- get_spot_indices_for_exp(
      object,
      exp_id = id,
      filtered = TRUE
    )
    adjacency <- .sp_require_graph(object, id, as.integer(spots_to_keep))

    if (verbosity > 0L) {
      message(sprintf(
        "Moran's I for '%s' over %i genes.",
        id,
        length(gene_sel$gene_indices)
      ))
    }

    res <- rs_sp_morans_i(
      f_path_gene = f_path_gene,
      spots_to_keep = as.integer(spots_to_keep),
      neighbours = adjacency$indices,
      weights = adjacency$weights,
      gene_indices = gene_sel$gene_indices,
      svg_params = svg_params,
      streaming = streaming,
      verbose = verbosity
    )

    object <- set_per_sample_morans_i(
      x = object,
      exp_id = id,
      morans_i = new_sp_morans_res(
        morans_res = res,
        gene_ids = gene_sel$gene_ids,
        exp_id = id
      )
    )
  }

  return(object)
}

### sparkx ---------------------------------------------------------------------

S7::method(sparkx_sp, SpOrSpSubset) <- function(
  object,
  exp_id = NULL,
  genes = NULL,
  sparkx_params = params_sp_sparkx(),
  seed = 42L,
  streaming = TRUE,
  .verbose = TRUE
) {
  # checks
  .assert_sp_object(object)
  assertSpSparkx(sparkx_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  exp_ids <- .sp_resolve_exp_ids(object, exp_id)
  verbosity <- parse_verbosity(.verbose)
  gene_sel <- .sp_resolve_genes(object, genes)
  f_path_gene <- get_rust_count_gene_f_path(object)

  for (id in exp_ids) {
    coords <- get_spatial_coords(object, exp_id = id, filtered = TRUE)
    spots_to_keep <- get_spot_indices_for_exp(
      object,
      exp_id = id,
      filtered = TRUE
    )

    if (nrow(coords) != length(spots_to_keep)) {
      stop(sprintf(
        "Got %i coordinate rows for %i spots of exp_id '%s'.",
        nrow(coords),
        length(spots_to_keep),
        id
      ))
    }

    if (verbosity > 0L) {
      message(sprintf(
        "SPARK-X for '%s' over %i genes.",
        id,
        length(gene_sel$gene_indices)
      ))
    }

    res <- rs_sp_sparkx(
      f_path_gene = f_path_gene,
      spots_to_keep = as.integer(spots_to_keep),
      coords = coords,
      gene_indices = gene_sel$gene_indices,
      sparkx_params = sparkx_params,
      seed = seed,
      streaming = streaming,
      verbose = verbosity
    )

    object <- set_per_sample_sparkx(
      x = object,
      exp_id = id,
      sparkx = new_sp_sparkx_res(
        sparkx_res = res,
        gene_ids = gene_sel$gene_ids,
        exp_id = id
      )
    )
  }

  return(object)
}

### neighbourhood enrichment ---------------------------------------------------

S7::method(nhood_enrichment_sp, SpOrSpSubset) <- function(
  object,
  label_col,
  exp_id = NULL,
  nhood_params = params_sp_nhood(),
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  .assert_sp_object(object)
  checkmate::qassert(label_col, "S1")
  assertSpNhood(nhood_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  exp_ids <- .sp_resolve_exp_ids(object, exp_id)
  verbosity <- parse_verbosity(.verbose)

  for (id in exp_ids) {
    spots_to_keep <- get_spot_indices_for_exp(
      object,
      exp_id = id,
      filtered = TRUE
    )
    adjacency <- .sp_require_graph(object, id, as.integer(spots_to_keep))
    obs <- .sp_obs_for_exp(object, id, cols = label_col)

    if (nrow(obs) != length(adjacency$indices)) {
      stop(sprintf(
        "The cached graph of '%s' has %i spots, the obs table %i.",
        id,
        length(adjacency$indices),
        nrow(obs)
      ))
    }

    labels_raw <- obs[[label_col]]
    if (anyNA(labels_raw)) {
      stop(sprintf(
        "Column '%s' has %i missing value(s) in exp_id '%s'.",
        label_col,
        sum(is.na(labels_raw)),
        id
      ))
    }

    # 0-based label codes, local to this sample
    label_factor <- factor(as.character(labels_raw))
    label_levels <- levels(label_factor)
    if (length(label_levels) < 2L) {
      stop(sprintf(
        "Column '%s' has a single level in exp_id '%s'.",
        label_col,
        id
      ))
    }

    if (verbosity > 0L) {
      message(sprintf(
        "Neighbourhood enrichment for '%s' over %i levels of '%s'.",
        id,
        length(label_levels),
        label_col
      ))
    }

    res <- rs_sp_nhood_enrichment(
      indices = adjacency$indices,
      weights = adjacency$weights,
      labels = as.integer(label_factor) - 1L,
      label_levels = label_levels,
      nhood_params = nhood_params,
      seed = seed,
      verbose = verbosity
    )

    object <- set_per_sample_nhood_enrichment(
      x = object,
      exp_id = id,
      label_col = label_col,
      result = new_sp_nhood_res(
        nhood_res = res,
        label_col = label_col,
        exp_id = id,
        n_perm = nhood_params$n_perm
      )
    )
  }

  return(object)
}

### image features -------------------------------------------------------------

S7::method(image_features_sp, SpOrSpSubset) <- function(
  object,
  exp_id = NULL,
  resolution = c("hires", "lowres", "cytassist", "fullres"),
  image_params = params_sp_image(),
  .verbose = TRUE
) {
  resolution <- match.arg(resolution)

  # checks
  .assert_sp_object(object)
  checkmate::assertChoice(
    resolution,
    c("hires", "lowres", "cytassist", "fullres")
  )
  assertSpImage(image_params)
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  exp_ids <- .sp_resolve_exp_ids(object, exp_id)
  verbosity <- parse_verbosity(.verbose)
  dir_data <- S7::prop(object, "dir_data")

  for (id in exp_ids) {
    sample <- get_sample(object, id)
    rel_path <- sample$image_paths[[resolution]]
    if (is.null(rel_path)) {
      stop(sprintf(
        "No image at resolution '%s' registered for exp_id '%s'.",
        resolution,
        id
      ))
    }
    image_path <- file.path(dir_data, rel_path)
    checkmate::assertFileExists(image_path)

    spot_diameter <- sample$scale_factors[["spot_diameter_fullres"]]
    if (is.null(spot_diameter) || !is.finite(spot_diameter)) {
      stop(sprintf(
        "Sample '%s' has no usable 'spot_diameter_fullres'.",
        id
      ))
    }

    coords <- get_spatial_coords(object, exp_id = id, filtered = TRUE)
    # 0-based global spot indices, back to the 1-indexed `cell_idx` the obs
    # table uses, so the features can be joined onto it later.
    spot_cell_idx <- as.integer(
      get_spot_indices_for_exp(object, exp_id = id, filtered = TRUE) + 1L
    )

    if (verbosity > 0L) {
      message(sprintf(
        "Cutting %i tiles out of the '%s' image of '%s'.",
        nrow(coords),
        resolution,
        id
      ))
    }

    res <- rs_sp_image_features(
      image_path = image_path,
      image_scalef = .sp_image_scalef(sample, resolution),
      coords = coords,
      spot_diameter_fullres = as.numeric(spot_diameter),
      image_params = image_params,
      verbose = verbosity
    )

    object <- set_per_sample_image_features(
      x = object,
      exp_id = id,
      resolution = resolution,
      features = new_sp_image_features_res(
        image_res = res,
        spot_cell_idx = spot_cell_idx,
        exp_id = id,
        resolution = resolution
      )
    )
  }

  return(object)
}
