# single cell processing methods -----------------------------------------------

## helpers ---------------------------------------------------------------------

### remapping communities ------------------------------------------------------

#' Remap communities by size
#'
#' @description
#' Remap community assignment by size
#'
#' @param labels Integer vector that represents the community assignment.
#'
#' @returns The remapped communities by size
#'
#' @keywords internal
.remap_communities_by_size <- function(labels) {
  counts <- sort(table(labels), decreasing = TRUE)
  remap <- setNames(seq_along(counts) - 1L, names(counts))
  unname(remap[as.character(labels)])
}

### doublet detection helpers --------------------------------------------------

#' Convert verbose argument to integer level
#'
#' @param .verbose `TRUE`, `FALSE`, or an integer.
#'
#' @return Integer verbosity level.
#'
#' @keywords internal
.verbose_level <- function(.verbose) {
  if (isTRUE(.verbose)) {
    return(1L)
  }
  if (isFALSE(.verbose)) {
    return(0L)
  }
  as.integer(.verbose)
}

#' Reduce verbosity by one level
#'
#' @param .verbose `TRUE`, `FALSE`, or an integer.
#'
#' @return Integer verbosity level, floored at 0.
#'
#' @keywords internal
.demote_verbosity <- function(.verbose) {
  max(0L, .verbose_level(.verbose) - 1L)
}

#' Assert that a grouping column exists in the obs table
#'
#' @param object A `SingleCells` object.
#' @param group_by Character. Column name to validate.
#'
#' @return Invisibly `NULL`. Stops with an informative message if the column
#' is absent.
#'
#' @keywords internal
.assert_group_by <- function(object, group_by) {
  duckdb_con <- get_sc_duckdb(object)
  obs_cols <- duckdb_con$get_obs_cols()
  if (!(group_by %in% obs_cols)) {
    stop(sprintf(
      "Column '%s' not found in the obs table. Available: %s",
      group_by,
      paste(obs_cols, collapse = ", ")
    ))
  }
}

#' Split cell indices into per-group lists
#'
#' Filters the obs table to `cells_to_use` and splits the resulting (0-indexed)
#' cell indices by the levels of `group_by`.
#'
#' @param object A `SingleCells` object.
#' @param group_by Character. Column name to group by.
#' @param cells_to_use Integer vector of 0-indexed cell indices.
#'
#' @return Named list of 0-indexed integer vectors, one element per group.
#'
#' @keywords internal
.split_cells_by_group <- function(object, group_by, cells_to_use) {
  # cells_to_use is rust 0-indexed; obs cell_idx is R 1-based
  obs <- get_sc_obs(
    object,
    cols = c("cell_idx", group_by),
    filtered = TRUE
  )

  obs <- obs[obs$cell_idx %in% (cells_to_use + 1L), ]

  group_vec <- as.character(obs[[group_by]])
  if (any(is.na(group_vec))) {
    stop(sprintf(
      "Column '%s' contains %d NA value(s); cannot group on it.",
      group_by,
      sum(is.na(group_vec))
    ))
  }

  split(as.integer(obs$cell_idx - 1L), group_vec)
}

#' Validate that all groups meet minimum cell count requirements
#'
#' Stops if any group has fewer than 50 cells. Warns if any group has fewer
#' than 500 cells.
#'
#' @param groups Named list of cell index vectors as returned by
#' [.split_cells_by_group()].
#'
#' @return Invisibly `NULL`.
#'
#' @keywords internal
.validate_group_sizes <- function(groups) {
  sizes <- lengths(groups)
  too_small <- sizes < 50L
  if (any(too_small)) {
    stop(sprintf(
      "Group(s) with fewer than 50 cells (hard minimum): %s",
      paste(
        sprintf("'%s' (n=%d)", names(sizes)[too_small], sizes[too_small]),
        collapse = ", "
      )
    ))
  }
  low <- sizes < 500L
  if (any(low)) {
    warning(sprintf(
      "Group(s) with fewer than 500 cells (results may be unreliable): %s",
      paste(
        sprintf("'%s' (n=%d)", names(sizes)[low], sizes[low]),
        collapse = ", "
      )
    ))
  }
  invisible(NULL)
}

#' Run a function over each group with optional progress reporting
#'
#' @param groups Named list of cell index vectors.
#' @param per_group_fn Function with signature `(cells, name, inner_verbose)`.
#' @param .verbose Verbosity level passed in from the calling function.
#' One level is demoted before being forwarded to `per_group_fn`.
#' @param label Character. Label shown on the progress bar.
#'
#' @return Named list of results, one element per group.
#'
#' @keywords internal
.run_per_group <- function(groups, per_group_fn, .verbose, label) {
  show_progress <- .verbose_level(.verbose) > 0L
  inner_verbose <- .demote_verbosity(.verbose)

  if (show_progress) {
    cli::cli_progress_bar(label, total = length(groups))
  }

  results <- vector("list", length(groups))
  names(results) <- names(groups)

  for (i in seq_along(groups)) {
    g <- names(groups)[i]
    if (show_progress) {
      cli::cli_progress_update(status = g)
    }
    results[[g]] <- per_group_fn(groups[[g]], g, inner_verbose)
  }

  if (show_progress) {
    cli::cli_progress_done()
  }

  results
}

## per-class concatenation -----------------------------------------------------

#' Extract and concatenate a named field from a list of group results
#'
#' @param group_results Named list of per-group result objects.
#' @param field Character. Name of the field to extract from each element.
#'
#' @return Unnamed vector of concatenated values.
#'
#' @keywords internal
.concat_per_cell <- function(group_results, field) {
  unlist(lapply(group_results, `[[`, field), use.names = FALSE)
}

#' Concatenate per-group Scrublet results into a single `ScrubletRes` object
#'
#' Cell-level vectors are re-ordered by the original cell index. Per-group
#' scalar summaries (threshold, rates) are retained as named vectors.
#'
#' @param group_results Named list of `ScrubletRes` objects.
#' @param group_by Character. Name of the grouping column; stored as an
#' attribute on the returned object.
#' @param return_combined_pca Logical. If `TRUE`, PCA results are included
#' per group.
#' @param return_pairs Logical. If `TRUE`, simulated doublet pair indices are
#' included per group.
#'
#' @return A `ScrubletRes` object with attributes `cell_indices`, `grouped`,
#' and `group_by_col`.
#'
#' @keywords internal
.concat_scrublet <- function(
  group_results,
  group_by,
  return_combined_pca,
  return_pairs
) {
  all_cell_idx <- unlist(
    lapply(group_results, function(r) attr(r, "cell_indices")),
    use.names = FALSE
  )
  ord <- order(all_cell_idx)

  group_labels <- unlist(
    Map(
      function(name, r) rep(name, length(r$predicted_doublets)),
      names(group_results),
      group_results
    ),
    use.names = FALSE
  )

  res <- list(
    predicted_doublets = .concat_per_cell(group_results, "predicted_doublets")[
      ord
    ],
    doublet_scores_obs = .concat_per_cell(group_results, "doublet_scores_obs")[
      ord
    ],
    doublet_errors_obs = .concat_per_cell(group_results, "doublet_errors_obs")[
      ord
    ],
    z_scores = .concat_per_cell(group_results, "z_scores")[ord],
    threshold = vapply(group_results, `[[`, numeric(1), "threshold"),
    detected_doublet_rate = vapply(
      group_results,
      `[[`,
      numeric(1),
      "detected_doublet_rate"
    ),
    detectable_doublet_fraction = vapply(
      group_results,
      `[[`,
      numeric(1),
      "detectable_doublet_fraction"
    ),
    overall_doublet_rate = vapply(
      group_results,
      `[[`,
      numeric(1),
      "overall_doublet_rate"
    ),
    doublet_scores_sim = lapply(group_results, `[[`, "doublet_scores_sim"),
    cell_groups = group_labels[ord]
  )

  if (return_combined_pca) {
    res$pca <- lapply(group_results, `[[`, "pca")
  }
  if (return_pairs) {
    res$pair_1 <- lapply(group_results, `[[`, "pair_1")
    res$pair_2 <- lapply(group_results, `[[`, "pair_2")
  }

  attr(res, "cell_indices") <- sort(all_cell_idx)
  attr(res, "grouped") <- TRUE
  attr(res, "group_by_col") <- group_by
  class(res) <- "ScrubletRes"
  res
}

#' Concatenate per-group boost results into a single `BoostRes` object
#'
#' Cell-level vectors are re-ordered by the original cell index.
#'
#' @param group_results Named list of `BoostRes` objects.
#' @param group_by Character. Name of the grouping column; stored as an
#' attribute on the returned object.
#'
#' @return A `BoostRes` object with attributes `cell_indices`, `grouped`,
#' and `group_by_col`.
#'
#' @keywords internal
.concat_boost <- function(group_results, group_by) {
  all_cell_idx <- unlist(
    lapply(group_results, function(r) attr(r, "cell_indices")),
    use.names = FALSE
  )
  ord <- order(all_cell_idx)

  group_labels <- unlist(
    Map(
      function(name, r) rep(name, length(r$doublet)),
      names(group_results),
      group_results
    ),
    use.names = FALSE
  )

  res <- list(
    doublet = .concat_per_cell(group_results, "doublet")[ord],
    doublet_score = .concat_per_cell(group_results, "doublet_score")[ord],
    voting_avg = .concat_per_cell(group_results, "voting_avg")[ord],
    cell_groups = group_labels[ord]
  )

  attr(res, "cell_indices") <- sort(all_cell_idx)
  attr(res, "grouped") <- TRUE
  attr(res, "group_by_col") <- group_by
  class(res) <- "BoostRes"
  res
}

#' Concatenate per-group scDblFinder results into a single `ScDblFinderRes`
#' object
#'
#' Cell-level vectors are re-ordered by the original cell index. Cluster labels
#' are prefixed with the group name to avoid collisions across groups. If
#' features were requested, per-group feature matrices are row-bound and
#' similarly reordered.
#'
#' @param group_results Named list of `ScDblFinderRes` objects.
#' @param group_by Character. Name of the grouping column; stored as an
#' attribute on the returned object.
#'
#' @return A `ScDblFinderRes` object with attributes `cell_indices`, `grouped`,
#' and `group_by_col`.
#'
#' @keywords internal
.concat_scdblfinder <- function(group_results, group_by) {
  all_cell_idx <- unlist(
    lapply(group_results, function(r) attr(r, "cell_indices")),
    use.names = FALSE
  )
  ord <- order(all_cell_idx)

  cluster_labels <- unlist(
    Map(
      function(name, r) sprintf("%s_%s", name, as.character(r$cluster_labels)),
      names(group_results),
      group_results
    ),
    use.names = FALSE
  )

  group_labels <- unlist(
    Map(
      function(name, r) rep(name, length(r$predicted_doublets)),
      names(group_results),
      group_results
    ),
    use.names = FALSE
  )

  has_features <- !is.null(group_results[[1]]$features)
  features <- if (has_features) {
    mat <- do.call(rbind, lapply(group_results, `[[`, "features"))
    mat[ord, , drop = FALSE]
  } else {
    NULL
  }

  structure(
    list(
      predicted_doublets = .concat_per_cell(
        group_results,
        "predicted_doublets"
      )[ord],
      doublet_score = .concat_per_cell(group_results, "doublet_score")[ord],
      cxds_scores = .concat_per_cell(group_results, "cxds_scores")[ord],
      weighted = .concat_per_cell(group_results, "weighted")[ord],
      threshold = vapply(group_results, `[[`, numeric(1), "threshold"),
      cluster_labels = cluster_labels[ord],
      detected_doublet_rate = vapply(
        group_results,
        `[[`,
        numeric(1),
        "detected_doublet_rate"
      ),
      features = features,
      cell_groups = group_labels[ord]
    ),
    cell_indices = sort(all_cell_idx),
    grouped = TRUE,
    group_by_col = group_by,
    class = "ScDblFinderRes"
  )
}

### dispatchers ----------------------------------------------------------------

#' Run Scrublet doublet detection on a set of cells
#'
#' @param object A `SingleCells` object.
#' @param cells_to_use Integer vector of 0-indexed cell indices.
#' @param scrublet_params List of Scrublet parameters from [params_scrublet()].
#' @param seed Integer. Random seed.
#' @param streaming Optional logical. Whether to stream the count data. If
#' `NULL`, resolved automatically via [auto_streaming()].
#' @param return_combined_pca Logical. If `TRUE`, PCA embeddings for observed
#' cells and simulated doublets are included in the result.
#' @param return_pairs Logical. If `TRUE`, parent cell indices of simulated
#' doublets are included in the result.
#' @param .verbose Logical or integer. Verbosity level.
#'
#' @return A `ScrubletRes` object with `cell_indices` set as an attribute.
#'
#' @keywords internal
.scrublet_run <- function(
  object,
  cells_to_use,
  scrublet_params,
  seed,
  streaming,
  return_combined_pca,
  return_pairs,
  .verbose
) {
  streaming <- auto_streaming(
    n_cells = length(cells_to_use),
    streaming = streaming,
    .verbose = .verbose
  )

  scrublet_res <- rs_sc_scrublet(
    f_path_gene = get_rust_count_gene_f_path(object),
    f_path_cell = get_rust_count_cell_f_path(object),
    cells_to_keep = cells_to_use,
    scrublet_params = scrublet_params,
    seed = seed,
    verbose = parse_verbosity(.verbose),
    streaming = streaming,
    return_combined_pca = return_combined_pca,
    return_pairs = return_pairs
  )

  attr(scrublet_res, "cell_indices") <- cells_to_use
  class(scrublet_res) <- "ScrubletRes"
  scrublet_res
}

#' Run boosted doublet detection on a set of cells
#'
#' @param object A `SingleCells` object.
#' @param cells_to_use Integer vector of 0-indexed cell indices.
#' @param boost_params List of boost parameters from [params_boost()].
#' @param seed Integer. Random seed.
#' @param streaming Optional logical. Whether to stream the count data. If
#' `NULL`, resolved automatically via [auto_streaming()].
#' @param .verbose Logical or integer. Verbosity level.
#'
#' @return A `BoostRes` object with `cell_indices` set as an attribute.
#'
#' @keywords internal
.boost_run <- function(
  object,
  cells_to_use,
  boost_params,
  seed,
  streaming,
  .verbose
) {
  streaming <- auto_streaming(
    n_cells = length(cells_to_use),
    streaming = streaming,
    .verbose = .verbose
  )

  if (boost_params$fast_cluster & is.null(boost_params$n_centroids)) {
    message(paste(
      "Fast clustering activated without any n_centroids set.",
      "Setting n_centroids to sqrt(N) * 2"
    ))
    boost_params$n_centroids <- as.integer(sqrt(length(cells_to_use)) * 4)
  }

  boost_res <- rs_sc_doublet_detection(
    f_path_gene = get_rust_count_gene_f_path(object),
    f_path_cell = get_rust_count_cell_f_path(object),
    cells_to_keep = cells_to_use,
    boost_params = boost_params,
    seed = seed,
    verbose = parse_verbosity(.verbose),
    streaming = streaming
  )

  attr(boost_res, "cell_indices") <- cells_to_use
  class(boost_res) <- "BoostRes"
  boost_res
}

#' Run scDblFinder doublet detection on a set of cells
#'
#' @param object A `SingleCells` object.
#' @param cells_to_use Integer vector of 0-indexed cell indices.
#' @param scdblfinder_params List of scDblFinder parameters from
#' [params_scdblfinder()].
#' @param return_features Logical. If `TRUE`, the classifier feature matrix is
#' included in the result, with cells as rows and feature names as column names.
#' @param streaming Optional logical. Whether to stream the count data. If
#' `NULL`, resolved automatically via [auto_streaming()].
#' @param seed Integer. Random seed.
#' @param .verbose Logical or integer. Verbosity level.
#'
#' @return A `ScDblFinderRes` object with `cell_indices` set as an attribute.
#'
#' @keywords internal
.scdblfinder_run <- function(
  object,
  cells_to_use,
  scdblfinder_params,
  return_features,
  streaming,
  seed,
  .verbose
) {
  streaming <- auto_streaming(
    n_cells = length(cells_to_use),
    streaming = streaming,
    .verbose = .verbose
  )

  if (
    scdblfinder_params$fast_cluster & is.null(scdblfinder_params$n_centroids)
  ) {
    message(paste(
      "Fast clustering activated without any n_centroids set.",
      "Setting n_centroids to sqrt(N) * 2"
    ))
    scdblfinder_params$n_centroids <- as.integer(sqrt(length(cells_to_use)) * 4)
  }

  res <- rs_sc_scdblfinder(
    f_path_gene = get_rust_count_gene_f_path(object),
    f_path_cell = get_rust_count_cell_f_path(object),
    cell_indices = cells_to_use,
    params = scdblfinder_params,
    return_features = return_features,
    streaming = streaming,
    seed = seed,
    verbose = parse_verbosity(.verbose)
  )

  features <- if (return_features) {
    feature_matrix <- res$features$feature_mat
    colnames(feature_matrix) <- res$features$feature_names
    rownames(feature_matrix) <- get_cell_names(object, filtered = TRUE)[
      cells_to_use + 1L
    ]
    feature_matrix
  } else {
    NULL
  }

  structure(
    list(
      predicted_doublets = res$predicted_doublets,
      doublet_score = res$doublet_scores,
      cxds_scores = res$cxds_scores,
      weighted = res$weighted,
      threshold = res$threshold,
      cluster_labels = res$cluster_labels,
      detected_doublet_rate = res$detected_doublet_rate,
      features = features
    ),
    cell_indices = cells_to_use,
    class = "ScDblFinderRes"
  )
}

## doublet detection -----------------------------------------------------------

### scrublet -------------------------------------------------------------------

#' Doublet detection with Scrublet
#'
#' @description This function implements the doublet detection from Scrublet,
#' see Wolock, et al. Briefly, arteficial doublets are being generated from
#' the data via random combination of initial cells. Highly variable genes (HVG)
#' are being identified and the observed cells are being projected on a PCA
#' space. Subsequently, the simulated doublets are being projected on the same
#' PCA space given the same HVGs; kNN graphs are being generated and a kNN
#' classifier is used to assign a probability that a given cell in the original
#' data is a doublet. For more details, please check the publication.
#'
#' @param object `SingleCells` class.
#' @param scrublet_params A list with the final scrublet parameters, see
#' [bixverse::params_scrublet()] for full details.
#' @param seed Integer. Random seed.
#' @param streaming Optional Boolean. Shall the data be streamed in. Useful for
#' larger data sets where you wish to avoid loading in the whole data. If
#' `NULL`, will automatically detect.
#' @param cells_to_use Optional string. Names of the cells to use for the
#' generation of the Scrublet. Useful when you wish to run doublet detection
#' on individual batches within your data. The object returned will be
#' specifically using these cells.
#' @param group_by Optional grouping variable. Useful if you want to run the
#' method on a per-sample basis.
#' @param return_combined_pca Boolean. Shall the PCA of the observed cells and
#' simulated doublets be returned.
#' @param return_pairs Boolean. Shall the pairs be returned.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @return A `scrublet_res` class that has with the following items:
#' \itemize{
#'   \item predicted_doublets - Boolean vector indicating which observed cells
#'   predicted as doublets (TRUE = doublet, FALSE = singlet).
#'   \item doublet_scores_obs - Numerical vector with the likelihood of being
#'   a doublet for the observed cells.
#'   \item doublet_scores_sim - Numerical vector with the likelihood of being
#'   a doublet for the simulated cells.
#'   \item doublet_errors_obs - Numerical vector with the standard errors of
#'   the scores for the observed cells.
#'   \item z_scores - Z-scores for the observed cells. Represents:
#'   `score - threshold / error`.
#'   \item threshold - Used threshold.
#'   \item detected_doublet_rate - Fraction of cells that are called as
#'   doublet.
#'   \item detectable_doublet_fraction - Fraction of simulated doublets with
#'   scores above the threshold.
#'   \item overall_doublet_rate - Estimated overall doublet rate. Should roughly
#'   match the expected doublet rate.
#'   \item pca - Optional PCA embeddings across the original cells and simulated
#'   doublets.
#'   \item pair_1 - Optional index of the parent cell 1 of the simulated
#'   doublets.
#'   \item pair_2 - Optional index of the parent cell 2 of the simulated
#'   doublets.
#' }
#'
#' @export
#'
#' @references Wollock, et al., Cell Syst, 2020
scrublet_sc <- S7::new_generic(
  name = "scrublet_sc",
  dispatch_args = "object",
  fun = function(
    object,
    scrublet_params = params_scrublet(),
    seed = 42L,
    streaming = NULL,
    cells_to_use = NULL,
    group_by = NULL,
    return_combined_pca = FALSE,
    return_pairs = FALSE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method scrublet_sc SingleCells
S7::method(scrublet_sc, SingleCells) <- function(
  object,
  scrublet_params = params_scrublet(),
  seed = 42L,
  streaming = NULL,
  cells_to_use = NULL,
  group_by = NULL,
  return_combined_pca = FALSE,
  return_pairs = FALSE,
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScScrublet(scrublet_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(group_by, c("S1", "0"))
  checkmate::qassert(return_combined_pca, "B1")
  checkmate::qassert(return_pairs, "B1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  cells_to_use <- if (!is.null(cells_to_use)) {
    get_cell_indices(object, cell_ids = cells_to_use, rust_index = TRUE)
  } else {
    get_cells_to_keep(object)
  }

  if (is.null(group_by)) {
    return(.scrublet_run(
      object = object,
      cells_to_use = cells_to_use,
      scrublet_params = scrublet_params,
      seed = seed,
      streaming = streaming,
      return_combined_pca = return_combined_pca,
      return_pairs = return_pairs,
      .verbose = .verbose
    ))
  }

  .assert_group_by(object, group_by)
  groups <- .split_cells_by_group(object, group_by, cells_to_use)
  .validate_group_sizes(groups)

  group_results <- .run_per_group(
    groups = groups,
    per_group_fn = function(cells, name, inner_v) {
      .scrublet_run(
        object = object,
        cells_to_use = cells,
        scrublet_params = scrublet_params,
        seed = seed,
        streaming = streaming,
        return_combined_pca = return_combined_pca,
        return_pairs = return_pairs,
        .verbose = inner_v
      )
    },
    .verbose = .verbose,
    label = "Running Scrublet per group"
  )

  .concat_scrublet(group_results, group_by, return_combined_pca, return_pairs)
}

### boost ----------------------------------------------------------------------

#' Doublet detection with boosted doublet classification
#'
#' @description This function implements the boosted doublet detection. It
#' generates through several iterations simulated doublets, generate kNN graphs,
#' runs Louvain clustering and assesses how often an observed cells clsuters
#' together with the simulated doublets.
#'
#' @param object `SingleCells` class.
#' @param boost_params A list with the final scrublet parameters, see
#' [bixverse::params_boost()] for full details.
#' @param cells_to_use Optional string. Names of the cells to use for the
#' run of the boosted doublet detection. Useful when you wish to run doublet
#' detection on individual batches within your data. The object returned will be
#' specifically using these cells.
#' @param group_by Optional grouping variable. Useful if you want to run the
#' method on a per-sample basis.
#' @param seed Integer. Random seed.
#' @param streaming Optional Boolean. Shall the data be streamed in. Useful for
#' larger data sets where you wish to avoid loading in the whole data. If
#' `NULL`, will automatically detect.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @return A `boost_res` class that has with the following items:
#' \itemize{
#'   \item predicted_doublets - Boolean vector indicating which observed cells
#'   predicted as doublets (TRUE = doublet, FALSE = singlet).
#'   \item doublet_scores_obs - Numerical vector with the likelihood of being
#'   a doublet for the observed cells.
#'   \item voting_avg - Numerical vector with the average voting score.
#' }
#'
#' @export
doublet_detection_boost_sc <- S7::new_generic(
  name = "doublet_detection_boost_sc",
  dispatch_args = "object",
  fun = function(
    object,
    boost_params = params_boost(),
    cells_to_use = NULL,
    group_by = NULL,
    seed = 42L,
    streaming = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method doublet_detection_boost_sc SingleCells
S7::method(doublet_detection_boost_sc, SingleCells) <- function(
  object,
  boost_params = params_boost(),
  cells_to_use = NULL,
  group_by = NULL,
  seed = 42L,
  streaming = NULL,
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScBoost(boost_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(group_by, c("S1", "0"))
  checkmate::qassert(.verbose, c("B1", "N1[0,2]"))

  cells_to_use <- if (!is.null(cells_to_use)) {
    get_cell_indices(object, cell_ids = cells_to_use, rust_index = TRUE)
  } else {
    get_cells_to_keep(object)
  }

  if (is.null(group_by)) {
    return(.boost_run(
      object = object,
      cells_to_use = cells_to_use,
      boost_params = boost_params,
      seed = seed,
      streaming = streaming,
      .verbose = .verbose
    ))
  }

  .assert_group_by(object, group_by)
  groups <- .split_cells_by_group(object, group_by, cells_to_use)
  .validate_group_sizes(groups)

  group_results <- .run_per_group(
    groups = groups,
    per_group_fn = function(cells, name, inner_v) {
      .boost_run(
        object = object,
        cells_to_use = cells,
        boost_params = boost_params,
        seed = seed,
        streaming = streaming,
        .verbose = inner_v
      )
    },
    .verbose = .verbose,
    label = "Running boost per group"
  )

  .concat_boost(group_results, group_by)
}

### scdblfinder ----------------------------------------------------------------

#' Run scDblFinder doublet detection on a SingleCells object
#'
#' @description
#' Cluster-aware doublet detection using engineered features and a
#' gradient-boosted classifier. See Germain et al., F1000Research, 2022.
#'
#' @param object `SingleCells` class.
#' @param scdblfinder_params List. Parameters from
#' [bixverse::params_scdblfinder()].
#' @param cells_to_use Optional string. Names of the cells to use for the
#' run of the boosted doublet detection. Useful when you wish to run doublet
#' detection on individual batches within your data. The object returned will be
#' specifically using these cells.
#' @param group_by Optional grouping variable. Useful if you want to run the
#' method on a per-sample basis.
#' @param streaming Optional Boolean. Shall the data be streamed in. Useful for
#' larger data sets where you wish to avoid loading in the whole data. If
#' `NULL`, will automatically detect.
#' @param return_features Boolean. Shall the features used to train the
#' classifier be returned.
#' @param seed Integer. Seed for reproducibility.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @return An S3 object of class `ScDblFinderRes` containing:
#' \describe{
#'   \item{predicted_doublets}{Logical vector of doublet calls.}
#'   \item{doublet_score}{Numeric vector of classifier probabilities.}
#'   \item{cxds_scores}{Numeric vector of the cxds scores.}
#'   \item{weighted}{Numeric vector of the weighted scores.}
#'   \item{threshold}{The threshold used for calling.}
#'   \item{cluster_labels}{Integer vector of final cluster assignments.}
#'   \item{detected_doublet_rate}{Fraction of cells called as doublets.}
#' }
#' with `cell_indices` stored as an attribute.
#'
#' @export
scdblfinder_sc <- S7::new_generic(
  name = "scdblfinder_sc",
  dispatch_args = "object",
  fun = function(
    object,
    scdblfinder_params = params_scdblfinder(),
    return_features = FALSE,
    cells_to_use = NULL,
    group_by = NULL,
    streaming = NULL,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method scdblfinder_sc SingleCells
#'
#' @export
S7::method(scdblfinder_sc, SingleCells) <- function(
  object,
  scdblfinder_params = params_scdblfinder(),
  return_features = FALSE,
  cells_to_use = NULL,
  group_by = NULL,
  streaming = NULL,
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScDblFinder(scdblfinder_params)
  checkmate::qassert(return_features, "B1")
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(group_by, c("S1", "0"))
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  cells_to_use <- if (!is.null(cells_to_use)) {
    get_cell_indices(object, cell_ids = cells_to_use, rust_index = TRUE)
  } else {
    get_cells_to_keep(object)
  }

  if (is.null(group_by)) {
    return(.scdblfinder_run(
      object = object,
      cells_to_use = cells_to_use,
      scdblfinder_params = scdblfinder_params,
      return_features = return_features,
      streaming = streaming,
      seed = seed,
      .verbose = .verbose
    ))
  }

  .assert_group_by(object, group_by)
  groups <- .split_cells_by_group(object, group_by, cells_to_use)
  .validate_group_sizes(groups)

  group_results <- .run_per_group(
    groups = groups,
    per_group_fn = function(cells, name, inner_v) {
      .scdblfinder_run(
        object = object,
        cells_to_use = cells,
        scdblfinder_params = scdblfinder_params,
        return_features = return_features,
        streaming = streaming,
        seed = seed,
        .verbose = inner_v
      )
    },
    .verbose = .verbose,
    label = "Running scDblFinder per group"
  )

  .concat_scdblfinder(group_results, group_by)
}

## gene proportions ------------------------------------------------------------

### top n genes ----------------------------------------------------------------

#' Calculate the proportions of reads for the Top N genes
#'
#' @description
#' This is a helper function that calculates proportions of reads to the Top N
#' genes by expression in a given cell. High values here can indicate low
#' complexity, quality cells. The values will be automatically added to the
#' obs table.
#'
#' @param object `SingleCells` class.
#' @param top_n_vals Integer. The Top N thresholds to test.
#' @param streaming Optional Boolean. Shall the data be streamed in. Useful for
#' larger data sets where you wish to avoid loading in the whole data. If
#' `NULL`, will automatically detect.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @return It will add the columns based on the names in the `gene_set_list` to
#' the obs table.
#'
#' @export
top_genes_perc_sc <- S7::new_generic(
  name = "top_genes_perc_sc",
  dispatch_args = "object",
  fun = function(
    object,
    top_n_vals = c(25L, 50L, 100L),
    streaming = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method top_genes_perc_sc SingleCells
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(top_genes_perc_sc, SingleCells) <- function(
  object,
  top_n_vals = c(25L, 50L, 100L),
  streaming = NULL,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(top_n_vals, "I+")
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  streaming <- auto_streaming(
    n_cells = nrow(object),
    streaming = streaming,
    .verbose = .verbose
  )

  # function
  rs_results <- rs_sc_get_top_genes_perc(
    f_path_cell = get_rust_count_cell_f_path(object),
    top_n_vals = top_n_vals,
    cell_indices = get_cells_to_keep(object),
    streaming = streaming,
    verbose = parse_verbosity(.verbose)
  )

  names(rs_results) <- sprintf("top_%i_genes_percentage", top_n_vals)

  res <- new_sc_list(res = rs_results, cell_indices = get_cells_to_keep(object))

  duckdb_con <- get_sc_duckdb(object)

  duckdb_con$join_data_obs(get_data(res))

  return(object)
}

### gene set version -----------------------------------------------------------

#' Calculate the proportions of reads for specific gene sets
#'
#' @description
#' This is a helper function that calculates proportions of reads belonging to
#' given gene sets. This can be used for example for the calculation of
#' percentage mitochondrial reads per cell. These will be automatically added
#' to the obs table
#'
#' @param object `SingleCells` class.
#' @param gene_set_list A named list with each element containing the gene
#' identifiers of that set. These should be the same as
#' `get_gene_names(object)`!
#' @param streaming Optional Boolean. Shall the data be streamed in. Useful for
#' larger data sets where you wish to avoid loading in the whole data. If
#' `NULL`, will automatically detect.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @return It will add the columns based on the names in the `gene_set_list` to
#' the obs table.
#'
#' @export
gene_set_proportions_sc <- S7::new_generic(
  name = "gene_set_proportions_sc",
  dispatch_args = "object",
  fun = function(
    object,
    gene_set_list,
    streaming = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)


#' @method gene_set_proportions_sc SingleCells
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(gene_set_proportions_sc, SingleCells) <- function(
  object,
  gene_set_list,
  streaming = NULL,
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::SingleCells")
  checkmate::assertList(gene_set_list, names = "named", types = "character")
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  streaming <- auto_streaming(
    n_cells = nrow(object),
    streaming = streaming,
    .verbose = .verbose
  )

  gene_set_list_tidy <- purrr::map(gene_set_list, \(g) {
    get_gene_indices(object, gene_ids = g, rust_index = TRUE)
  })
  names(gene_set_list_tidy) <- names(gene_set_list)

  rs_results <- rs_sc_get_gene_set_perc(
    f_path_cell = get_rust_count_cell_f_path(object),
    cell_indices = get_cells_to_keep(object),
    gene_set_idx = gene_set_list_tidy,
    streaming = streaming,
    verbose = parse_verbosity(.verbose)
  )

  res <- new_sc_list(res = rs_results, cell_indices = get_cells_to_keep(object))

  duckdb_con <- get_sc_duckdb(object)

  duckdb_con$join_data_obs(get_data(res))

  return(object)
}

## hvg -------------------------------------------------------------------------

### manipulating object state --------------------------------------------------

# generic found in R/base_generics_sc.R

#' @method find_hvg_sc SingleCells
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(find_hvg_sc, SingleCells) <- function(
  object,
  hvg_no = 2000L,
  hvg_params = params_sc_hvg(),
  streaming = NULL,
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::SingleCells")
  checkmate::qassert(hvg_no, "I1")
  assertScHvg(hvg_params)
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(.verbose, c("B1", "I1[0, 2]"))

  streaming <- auto_streaming(
    n_cells = nrow(object),
    streaming = streaming,
    .verbose = .verbose
  )

  res <- with(
    hvg_params,
    rs_sc_hvg(
      f_path_gene = get_rust_count_gene_f_path(object),
      hvg_method = method,
      cell_indices = get_cells_to_keep(object),
      loess_span = loess_span,
      n_bins = num_bin,
      binning = bin_method,
      clip_max = NULL,
      streaming = streaming,
      verbose = parse_verbosity(.verbose)
    )
  )

  object <- set_sc_new_var_cols(object = object, data_list = res)

  hvg <- switch(
    hvg_params$method,
    "vst" = order(res$var_std, decreasing = TRUE)[1:hvg_no],
    "dispersion" = order(res$dispersion, decreasing = TRUE)[1:hvg_no],
    "meanvarbin" = order(res$dispersion_scaled, decreasing = TRUE)[1:hvg_no],
    stop("Unknown HVG method: ", hvg_params$method)
  )

  object <- set_hvg(object, hvg = hvg)

  return(object)
}

### without manipulating object state ------------------------------------------

#' @method get_hvg_data_sc SingleCells
#'
#' @export
S7::method(get_hvg_data_sc, SingleCells) <- function(
  object,
  cell_ids = NULL,
  hvg_no = 3000L,
  hvg_params = params_sc_hvg(),
  streaming = NULL,
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(cell_ids, c("0", "S+"))
  checkmate::qassert(hvg_no, "I1")
  assertScHvg(hvg_params)
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  cell_indices <- if (is.null(cell_ids)) {
    get_cells_to_keep(object)
  } else {
    get_cell_indices(object, cell_ids = cell_ids, rust_index = TRUE)
  }

  streaming <- auto_streaming(
    n_cells = length(cell_indices),
    streaming = streaming,
    .verbose = .verbose
  )

  res <- with(
    hvg_params,
    rs_sc_hvg(
      f_path_gene = get_rust_count_gene_f_path(object),
      hvg_method = method,
      cell_indices = cell_indices,
      loess_span = loess_span,
      n_bins = num_bin,
      binning = bin_method,
      clip_max = NULL,
      streaming = streaming,
      verbose = parse_verbosity(.verbose)
    )
  )

  var_table <- get_sc_var(object, cols = c("gene_idx", "gene_id"))

  build_hvg_table(
    var_table = var_table,
    res = res,
    hvg_no = hvg_no,
    hvg_method = hvg_params$method
  )
}

## dimension reduction and knn/snn ---------------------------------------------

### pca ------------------------------------------------------------------------

# generic found in R/base_generics_sc.R

#' @method calculate_pca_sc SingleCells
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(calculate_pca_sc, SingleCells) <- function(
  object,
  no_pcs,
  pca_params = params_sc_pca(),
  sparse_svd = FALSE,
  hvg = NULL,
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::SingleCells")
  checkmate::qassert(no_pcs, "I1")
  assertScPca(pca_params)
  checkmate::qassert(sparse_svd, "B1")
  checkmate::qassert(hvg, c("I+", "0"))
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  if ((length(get_hvg(object)) == 0) && is.null(hvg)) {
    warning(paste(
      "No HVGs identified in the object nor provided.",
      "Please run find_hvg_sc() or provide the indices of the HVG",
      "Returning object as is."
    ))
    return(object)
  }

  selected_hvg <- if (!is.null(hvg)) {
    if (.verbose) {
      message(
        paste(
          "HVGs provided.",
          "Will use these ones and set the internal HVG to the provided genes."
        )
      )
    }
    object <- set_hvg(object, hvg) # this one deals with zero/one indexing internally
    hvg - 1L
  } else {
    get_hvg(object)
  }

  # swap to sparse SVD for large data sets
  n_cells <- length(get_cells_to_keep(object))

  if (n_cells > 500000 & !sparse_svd) {
    message(paste(
      "More than 500,000 cells with sparse SVD = FALSE",
      "Setting sparse SVD to TRUE to avoid high memory pressure."
    ))

    sparse_svd <- TRUE
  }

  # dense path
  if (!sparse_svd) {
    if (.verbose) {
      message(
        sprintf(
          "Using dense SVD solving on scaled data on %i HVG.",
          length(selected_hvg)
        )
      )
    }
    zeallot::`%<-%`(
      c(pca_factors, pca_loadings, singular_values, scaled),
      rs_sc_pca(
        f_path_gene = get_rust_count_gene_f_path(object),
        f_path_cell = get_rust_count_cell_f_path(object),
        no_pcs = no_pcs,
        pca_params = pca_params,
        cell_indices = get_cells_to_keep(object),
        gene_indices = selected_hvg,
        seed = seed,
        return_scaled = FALSE,
        verbose = parse_verbosity(.verbose)
      )
    )

    object <- set_pca_factors(object, pca_factors)
    object <- set_pca_loadings(object, pca_loadings)
    object <- set_pca_singular_vals(object, singular_values[1:no_pcs])

    return(object)
  } else {
    if (.verbose) {
      message(
        sprintf(
          "Using sparse SVD solving on scaled data on %i HVG.",
          length(selected_hvg)
        )
      )
    }

    zeallot::`%<-%`(
      c(sparse_pca_factors, sparse_pca_loadings, sparse_pca_eigenvals),
      rs_sc_pca_sparse(
        f_path_gene = get_rust_count_gene_f_path(object),
        f_path_cell = get_rust_count_cell_f_path(object),
        no_pcs = no_pcs,
        pca_params = pca_params,
        cell_indices = get_cells_to_keep(object),
        gene_indices = selected_hvg,
        seed = seed,
        verbose = parse_verbosity(.verbose)
      )
    )

    object <- set_pca_factors(object, sparse_pca_factors)
    object <- set_pca_loadings(object, sparse_pca_loadings)
    object <- set_pca_singular_vals(
      object,
      sparse_pca_eigenvals[1:no_pcs]
    )

    return(object)
  }
}

### neighbours -----------------------------------------------------------------

# generic found in R/base_generics_sc.R
# method shared across SingleCells and MetaCells

S7::method(find_neighbours_sc, ScOrMc) <- function(
  object,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  modality = c("rna", "adt"),
  neighbours_params = params_sc_neighbours(),
  seed = 42L,
  .verbose = TRUE
) {
  modality <- match.arg(modality)

  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) ||
      S7::S7_inherits(object, MetaCells) ||
      S7::S7_inherits(object, SingleCellsSubset)
  )
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  assertScNeighbours(neighbours_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  if (modality != "rna" && !S7::S7_inherits(object, SingleCellsMultiModal)) {
    stop(sprintf(
      "modality = '%s' is only supported for SingleCellsMultiModal.",
      modality
    ))
  }

  if (!embd_to_use %in% get_available_embeddings(object, modality = modality)) {
    warning("The desired embedding was not found. Returning class as is.")
    return(object)
  }

  embd <- get_embedding(
    x = object,
    embd_name = embd_to_use,
    modality = modality
  )
  if (!is.null(no_embd_to_use)) {
    embd <- embd[, 1:min(no_embd_to_use, ncol(embd))]
  }

  if (.verbose) {
    message(sprintf(
      "Generating kNN data with %s method.",
      neighbours_params$knn_algorithm
    ))
  }
  knn_data <- generate_sc_knn(
    data = embd,
    neighbours_params = neighbours_params,
    seed = seed,
    .verbose = .verbose
  )
  object <- set_knn(object, knn_data, modality = modality)

  if (.verbose) {
    message(sprintf(
      "Generating sNN graph (full: %s).",
      neighbours_params$full_snn
    ))
  }
  snn_graph_rs <- with(
    neighbours_params,
    rs_sc_snn(
      knn_mat = get_knn_mat(knn_data),
      snn_method = snn_similarity,
      pruning = pruning,
      limited_graph = !full_snn,
      verbose = parse_verbosity(.verbose)
    )
  )

  if (.verbose) {
    message("Transforming sNN data to igraph.")
  }
  snn_g <- igraph::make_empty_graph(n = nrow(embd), directed = FALSE)
  snn_g <- igraph::add_edges(
    snn_g,
    snn_graph_rs$edges,
    attr = list(weight = snn_graph_rs$weights)
  )

  object <- set_snn_graph(object, snn_graph = snn_g, modality = modality)

  return(object)
}

### clustering -----------------------------------------------------------------

# generic found in R/base_generics_sc.R
# method shared across SingleCells and MetaCells

S7::method(find_clusters_sc, ScOrMc) <- function(
  object,
  cluster_algorithm = c("leiden", "louvain"),
  res = 1.0,
  name = "leiden_clustering",
  modality = c("rna", "adt", "wnn"),
  seed = 42L
) {
  cluster_algorithm <- match.arg(cluster_algorithm)
  modality <- match.arg(modality)

  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) ||
      S7::S7_inherits(object, MetaCells) ||
      S7::S7_inherits(object, SingleCellsSubset)
  )
  checkmate::qassert(res, "N1")
  checkmate::qassert(name, "S1")
  checkmate::assertChoice(cluster_algorithm, c("leiden", "louvain"))
  checkmate::qassert(seed, "I1")

  if (modality != "rna" && !S7::S7_inherits(object, SingleCellsMultiModal)) {
    stop(sprintf(
      "modality = '%s' is only supported for SingleCellsMultiModal.",
      modality
    ))
  }

  snn_graph <- get_snn_graph(object, modality = modality)
  if (is.null(snn_graph)) {
    warning(
      paste(
        "No sNN graph found. Did you run find_neighbours_sc().",
        "Returning class as is."
      )
    )
    return(object)
  }

  # set a global seed to ensure more reproducibility here
  set.seed(seed)
  clusters <- switch(
    cluster_algorithm,
    leiden = igraph::cluster_leiden(
      graph = snn_graph,
      objective_function = "modularity",
      resolution = res
    ),
    louvain = igraph::cluster_louvain(graph = snn_graph, resolution = res)
  )

  object[[name]] <- .remap_communities_by_size(clusters$membership)

  object
}

### fast clustering ------------------------------------------------------------

#' Run fast Louvain clustering on a SingleCells object
#'
#' @description
#' Runs k-means on the chosen embedding, builds a kNN graph on the centroids,
#' applies Louvain clustering and propagates memberships back to the cells.
#' Optionally runs a grid over multiple seeds and returns stability statistics.
#'
#' @param object `SingleCells` class.
#' @param embd_to_use String. Embedding name. Defaults to `"pca"`.
#' @param no_embd_to_use Optional integer. Number of dimensions to keep.
#' @param resolutions Numeric vector. Louvain resolutions.
#' @param km_type String. One of `c("kmeans", "minibatch")`. The former runs
#' standard k-means, the latter a mini-batch version that can be useful for
#' large data sets.
#' @param n_centroids Optional integer. Number of k-means centroids. Defaults
#' to `sqrt(n_cells)` Rust-side if `NULL`.
#' @param fc_params List. Output of [params_sc_fast_cluster()].
#' @param snn Boolean. Convert kNN to sNN.
#' @param return_kmeans Boolean. Return k-means assignments and centroids.
#' @param grid_search Boolean. Run multi-seed grid version.
#' @param no_seeds Integer. Number of additional seeds (only used when
#' `grid_search = TRUE`).
#' @param seed Integer. Reproducibility.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns `SingleCellFastClusters` S3 object with:
#' \describe{
#'   \item{memberships}{data.table with `cell_idx` and one column per
#'   resolution (`res_<value>`).}
#'   \item{stats}{data.table of grid statistics, or `NULL`.}
#'   \item{k_means_cluster}{Integer vector of k-means assignments, or `NULL`.}
#'   \item{centroids}{Numeric matrix of centroids, or `NULL`.}
#'   \item{resolutions}{Resolutions used.}
#' }
#' with `cell_indices` stored as an attribute (0-indexed).
#'
#' @export
fast_cluster_sc <- S7::new_generic(
  name = "fast_cluster_sc",
  dispatch_args = "object",
  fun = function(
    object,
    embd_to_use = "pca",
    no_embd_to_use = NULL,
    resolutions = c(2.0, 1.0, 0.5),
    km_type = c("kmeans", "minibatch"),
    n_centroids = NULL,
    fc_params = params_sc_fast_cluster(),
    snn = TRUE,
    return_kmeans = FALSE,
    grid_search = FALSE,
    no_seeds = 10L,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method fast_cluster_sc SingleCells
#'
#' @export
S7::method(fast_cluster_sc, SingleCells) <- function(
  object,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  resolutions = c(2.0, 1.0, 0.5),
  km_type = c("kmeans", "minibatch"),
  n_centroids = NULL,
  fc_params = params_sc_fast_cluster(),
  snn = TRUE,
  return_kmeans = FALSE,
  grid_search = FALSE,
  no_seeds = 10L,
  seed = 42L,
  .verbose = TRUE
) {
  km_type <- match.arg(km_type)

  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  checkmate::qassert(resolutions, "N+")
  checkmate::qassert(n_centroids, c("I1", "0"))
  checkmate::qassert(snn, "B1")
  checkmate::qassert(return_kmeans, "B1")
  checkmate::qassert(grid_search, "B1")
  checkmate::qassert(no_seeds, "I1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  if (!embd_to_use %in% get_available_embeddings(object)) {
    stop(sprintf("Embedding '%s' was not found.", embd_to_use))
  }

  embd <- get_embedding(x = object, embd_name = embd_to_use)
  if (!is.null(no_embd_to_use)) {
    embd <- embd[, 1:min(no_embd_to_use, ncol(embd))]
  }

  cells_to_use <- get_cells_to_keep(object)

  if (grid_search) {
    res <- rs_fast_cluster_sc_grid(
      embd = embd,
      km_type = km_type,
      resolutions = resolutions,
      n_centroids = n_centroids,
      fc_params = fc_params,
      snn = snn,
      return_kmeans = return_kmeans,
      no_seeds = no_seeds,
      seed = seed,
      verbose = parse_verbosity(.verbose)
    )
    memberships <- res$membership$memberships
    stats <- data.table::as.data.table(res$membership$stats)
    stats[, resolution := resolutions]
    data.table::setcolorder(stats, "resolution")
  } else {
    res <- rs_fast_cluster_sc(
      embd = embd,
      km_type = km_type,
      resolutions = resolutions,
      n_centroids = n_centroids,
      fc_params = fc_params,
      snn = snn,
      return_kmeans = return_kmeans,
      seed = seed,
      verbose = parse_verbosity(.verbose)
    )
    memberships <- res$membership
    stats <- NULL
  }

  membership_dt <- data.table::data.table(cell_idx = cells_to_use + 1L)
  membership_dt[, (paste0("res_", resolutions)) := memberships]

  structure(
    list(
      memberships = membership_dt,
      stats = stats,
      k_means_cluster = res$k_means_cluster,
      centroids = res$centroids,
      resolutions = resolutions
    ),
    cell_indices = cells_to_use,
    class = "SingleCellFastClusters"
  )
}

### generate knn class ---------------------------------------------------------

#' Generate a `SingleCellNearestNeighbour` from a single cell class
#'
#' @description
#' This function will generate the kNNs based on a given embedding. Available
#' algorithms are:
#' \itemize{
#'   \item `kmknn` - An exact kNN search that leverages k-means clustering under
#'   the hood to prune out data points. The default setting.
#'   \item `exhaustive` - An exhaustive, flat index. On smaller data sets often
#'   faster than the approximate nearest neighbour search algorithms.
#'   \item `hnsw` - Hierarchical Navigable Small World. A graph-based
#'   approximate nearest neighbour search algorithm; works well on large data
#'   sets. A benign race condition is leveraged during index build, making the
#'   build non-deterministic. Bigger impact on smaller data sets.
#'   \item `nndescent` - Nearest neighbour descent. Leverages concepts from
#'   `PyNNDescent` and works well on very large data sets similar to `hnsw`.
#'   \item `ivf` - Inverted file index. Uses first k-means clustering to
#'   identify Voronoi cells and leverages these during querying. Works well
#'   on large data sets with high dimensionality and when you need to return
#'   large number of neighbours.
#'   \item `annoy` - Approximate nearest neighbours Oh Yeah. Tree-based index,
#'   used across different R single cell packages (Seurat, SCE). This version
#'   is purely memory-based.
#' }
#'
#' @param object `SingleCells` class.
#' @param embd_to_use String. The embedding to use. Whichever you chose, it
#' needs to be part of the object.
#' @param cells_to_use String. Optional cell names to include in the generation
#' of the kNN graph. If `NULL` all (filtered) cells in the object will be used.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param neighbours_params List. Output of [bixverse::params_sc_neighbours()].
#' A list with the following items:
#' \itemize{
#'   \item full_snn - Boolean. Shall the full shared nearest neighbour graph
#'   be generated that generates edges between all cells instead of between
#'   only neighbours. Not used for this function.
#'   \item pruning - Numeric. Weights below this threshold will be set to 0 in
#'   the generation of the sNN graph. Not used for this function.
#'   \item snn_similarity - String. One of `c("rank", "jaccard")`. Defines how
#'   the weight from the SNN graph is calculated. For details, please see
#'   [bixverse::params_sc_neighbours()]. Not used for this function.
#'   \item knn - List of kNN parameters. See [bixverse::params_knn_defaults()]
#'   for available parameters and their defaults.
#' }
#' @param seed Integer. For reproducibility.
#' @param .validate_index Boolean. Shall an exhaustive search against a subset
#' of cells be run to validate the approximate nearest neighbour index.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @return Initialised `sc_knn` with the kNN data.
#'
#' @export
generate_knn_sc <- S7::new_generic(
  name = "generate_knn_sc",
  dispatch_args = "object",
  fun = function(
    object,
    embd_to_use = "pca",
    cells_to_use = NULL,
    no_embd_to_use = NULL,
    neighbours_params = params_sc_neighbours(),
    seed = 42L,
    .validate_index = TRUE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method generate_knn_sc SingleCells
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(generate_knn_sc, SingleCells) <- function(
  object,
  embd_to_use = "pca",
  cells_to_use = NULL,
  no_embd_to_use = NULL,
  neighbours_params = params_sc_neighbours(),
  seed = 42L,
  .validate_index = TRUE,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  checkmate::qassert(cells_to_use, c("S+", "0"))
  assertScNeighbours(neighbours_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.validate_index, "B1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # function body
  embd <- get_embedding(x = object, embd_name = embd_to_use)

  # early return
  if (is.null(embd)) {
    warning(
      paste(
        "The desired embedding was not found. Please check the parameters.",
        "Returning NULL."
      )
    )

    return(NULL)
  }

  if (!is.null(cells_to_use)) {
    embd <- embd[which(rownames(embd) %in% cells_to_use), ]
  }

  knn_data <- rs_sc_knn_w_dist(
    embd = embd,
    knn_params = neighbours_params,
    verbose = parse_verbosity(.verbose),
    validate_index = .validate_index,
    seed = seed
  )

  knn_obj <- new_sc_knn(knn_data = knn_data, used_cells = row.names(embd))

  knn_obj
}
