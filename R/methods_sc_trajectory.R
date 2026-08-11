# ------------------------------------------------------------------------------
# trajectory inference for single cells. both Palantir and PAGA only consume a
# kNN graph, so the methods are shared across all of the single cell classes via
# the ScOrMc union.
# ------------------------------------------------------------------------------

# helpers ----------------------------------------------------------------------

#' Guard against non-RNA modalities on single-modality classes
#'
#' @param object One of the single cell classes.
#' @param modality String. The requested modality.
#'
#' @returns `NULL`, invisibly. Errors if the modality is unavailable.
#'
#' @keywords internal
.assert_modality_available <- function(object, modality) {
  if (modality != "rna" && !S7::S7_inherits(object, SingleCellsMultiModal)) {
    stop(sprintf(
      "modality = '%s' is only supported for SingleCellsMultiModal.",
      modality
    ))
  }

  invisible(NULL)
}

#' Resolve cell names to kNN row positions
#'
#' @description
#' The kNN indices are positions within the cells the graph was built over, not
#' global cell indices, so cell names are resolved against `used_cells`.
#'
#' @param cell_ids Character vector. The cell names to resolve.
#' @param used_cells Character vector. The cells the kNN graph was built over.
#' @param arg_name String. Name of the argument, used in the error message.
#'
#' @returns Integer vector with the 0-indexed positions for Rust.
#'
#' @keywords internal
.cells_to_knn_idx <- function(cell_ids, used_cells, arg_name) {
  # checks
  checkmate::qassert(cell_ids, "S+")
  checkmate::qassert(used_cells, "S+")
  checkmate::qassert(arg_name, "S1")

  indices <- match(cell_ids, used_cells)

  if (anyNA(indices)) {
    stop(sprintf(
      "`%s`: %i cell(s) not found in the kNN graph, e.g. '%s'.",
      arg_name,
      sum(is.na(indices)),
      cell_ids[which(is.na(indices))[1]]
    ))
  }

  as.integer(indices - 1L)
}

# palantir ---------------------------------------------------------------------

#' Run Palantir trajectory inference
#'
#' @description
#' Palantir models differentiation as a Markov chain over a diffusion map of the
#' cells. It returns a pseudotime per cell, the probability of each cell
#' reaching each terminal state, and the differentiation entropy derived from
#' those probabilities. For details, please refer to Setty, et al.
#'
#' The kNN graph stored on the object feeds the diffusion kernel, so run
#' [bixverse::find_neighbours_sc()] first. Whether the stored distances are
#' treated as squared is derived from the distance metric of that graph. The
#' geodesics themselves are measured over a second kNN graph that Palantir
#' builds internally on the multiscale space, controlled by the `knn` element of
#' [bixverse::params_sc_palantir()].
#'
#' For `SingleCellsMultiModal` the `"wnn"` graph works, but its distances are
#' kernel-derived and not a metric, so the diffusion map on top is a heuristic
#' on a heuristic. You will get a warning saying as much.
#'
#' @param object One of `SingleCells`, `SingleCellsSubset`, `MetaCells` or
#' `SingleCellsMultiModal`.
#' @param early_cell String. Name of the cell to start the trajectory from. Must
#' be one of the cells the kNN graph was built over.
#' @param terminal_states Optional character vector. Names of the terminal state
#' cells. If `NULL`, they are detected from the waypoint Markov chain.
#' @param modality String. One of `c("rna", "adt", "wnn")`. Which kNN graph to
#' run over. Anything but `"rna"` requires a `SingleCellsMultiModal` object.
#' @param palantir_params List. See [bixverse::params_sc_palantir()].
#' @param seed Integer. For reproducibility.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns A `PalantirRes` S3 object with:
#' \itemize{
#'   \item pseudotime - data.table with `cell_id`, `pseudotime` (min-max scaled
#'   to `[0, 1]`; the start cell is not pinned to 0, and a start cell far from 0
#'   means the refinement disagreed with the anchor) and `entropy`.
#'   \item branch_probs - Numeric matrix of cells x terminal states with the
#'   fate probabilities. Rows need not sum to one, as sub-threshold values are
#'   zeroed without renormalisation.
#'   \item terminal_states - Character vector with the terminal state cell
#'   names. Sets the column order of `branch_probs`.
#'   \item waypoints - Character vector with the waypoint cell names. The first
#'   element is the start cell.
#'   \item start_cell - String. The start cell that was actually used.
#'   \item multiscale - Numeric matrix of cells x components with the multiscale
#'   diffusion components.
#'   \item run_info - List with `iterations`, `converged`, `eigen_converged`,
#'   `eigen_residual`, `repair_edges`, `stranded_waypoints`, `n_waypoints` and
#'   `modality`. A `eigen_converged` of `FALSE` means the diffusion eigensolve
#'   ran out of restarts and the embedding is under-resolved.
#' }
#'
#' @references Setty, et al., Nat. Biotechnol., 2019.
#'
#' @export
run_palantir_sc <- S7::new_generic(
  name = "run_palantir_sc",
  dispatch_args = "object",
  fun = function(
    object,
    early_cell,
    terminal_states = NULL,
    modality = c("rna", "adt", "wnn"),
    palantir_params = params_sc_palantir(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

# method shared across SingleCells, SingleCellsSubset and MetaCells

S7::method(run_palantir_sc, ScOrMc) <- function(
  object,
  early_cell,
  terminal_states = NULL,
  modality = c("rna", "adt", "wnn"),
  palantir_params = params_sc_palantir(),
  seed = 42L,
  .verbose = TRUE
) {
  modality <- match.arg(modality)

  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) ||
      S7::S7_inherits(object, MetaCells) ||
      S7::S7_inherits(object, SingleCellsSubset)
  )
  checkmate::qassert(early_cell, "S1")
  checkmate::qassert(terminal_states, c("S+", "0"))
  assertScPalantirParams(palantir_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  .assert_modality_available(object, modality)

  # hard tier: the kNN indices are handed to Rust as is
  assert_sc_state(object, artefacts = "knn", modality = modality)

  knn_data <- get_knn_obj(object, modality = modality)

  if (is.null(knn_data)) {
    warning(
      paste(
        "No kNN data found. Did you run find_neighbours_sc().",
        "Returning NULL."
      )
    )
    return(NULL)
  }

  if (!knn_data$dist_metric %in% c("euclidean", "cosine")) {
    warning(sprintf(
      paste(
        "The kNN graph uses '%s', which is not a metric.",
        "The diffusion map is built on it regardless, but interpret the",
        "pseudotime with care."
      ),
      knn_data$dist_metric
    ))
  }

  early_cell_idx <- .cells_to_knn_idx(
    early_cell,
    knn_data$used_cells,
    "early_cell"
  )

  terminal_state_idx <- if (is.null(terminal_states)) {
    NULL
  } else {
    .cells_to_knn_idx(terminal_states, knn_data$used_cells, "terminal_states")
  }

  if (.verbose) {
    message(sprintf(
      "Running Palantir over the %s kNN graph (%i cells).",
      modality,
      knn_data$n
    ))
  }

  palantir_res <- rs_palantir(
    knn_data = knn_data,
    palantir_params = palantir_params,
    early_cell = early_cell_idx,
    terminal_states = terminal_state_idx,
    seed = seed,
    verbose = parse_verbosity(.verbose)
  )

  new_palantir_res(
    rs_res = palantir_res,
    used_cells = knn_data$used_cells,
    modality = modality
  )
}

# paga -------------------------------------------------------------------------

#' Run PAGA graph abstraction
#'
#' @description
#' PAGA abstracts the single cell kNN graph into a cluster-level graph, where
#' the edge weight between two clusters measures how much more connected they
#' are than expected under a random null model. It also returns the maximum
#' spanning forest of that graph, which is the backbone typically plotted. For
#' details, please refer to Wolf, et al.
#'
#' Only the kNN indices are used, no distances and no expression data, so the
#' method is safe on any of the available graphs, including the WNN one. Run
#' [bixverse::find_neighbours_sc()] and a clustering first.
#'
#' @param object One of `SingleCells`, `SingleCellsSubset`, `MetaCells` or
#' `SingleCellsMultiModal`.
#' @param cluster_col String. The obs column holding the cluster assignments.
#' Empty factor levels are retained.
#' @param modality String. One of `c("rna", "adt", "wnn")`. Which kNN graph to
#' run over. Anything but `"rna"` requires a `SingleCellsMultiModal` object.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns A `PagaRes` S3 object with:
#' \itemize{
#'   \item connectivities - Symmetric sparse matrix of clusters x clusters with
#'   a zero diagonal. Values lie in `(0, 1]`.
#'   \item connectivities_tree - Maximum spanning forest of `connectivities`,
#'   carrying the original connectivity values on the retained edges.
#'   \item sizes - data.table with `cluster` and `n_cells`.
#'   \item params - List with `cluster_col` and `modality`.
#' }
#'
#' @references Wolf, et al., Genome Biol., 2019.
#'
#' @export
run_paga_sc <- S7::new_generic(
  name = "run_paga_sc",
  dispatch_args = "object",
  fun = function(
    object,
    cluster_col,
    modality = c("rna", "adt", "wnn"),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

# method shared across SingleCells, SingleCellsSubset and MetaCells

S7::method(run_paga_sc, ScOrMc) <- function(
  object,
  cluster_col,
  modality = c("rna", "adt", "wnn"),
  .verbose = TRUE
) {
  modality <- match.arg(modality)

  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) ||
      S7::S7_inherits(object, MetaCells) ||
      S7::S7_inherits(object, SingleCellsSubset)
  )
  checkmate::qassert(cluster_col, "S1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  .assert_modality_available(object, modality)

  # hard tier: the kNN indices are handed to Rust as is
  assert_sc_state(object, artefacts = "knn", modality = modality)

  knn_data <- get_knn_obj(object, modality = modality)

  if (is.null(knn_data)) {
    warning(
      paste(
        "No kNN data found. Did you run find_neighbours_sc().",
        "Returning NULL."
      )
    )
    return(NULL)
  }

  clusters <- as.factor(unlist(object[[cluster_col]], use.names = FALSE))

  # the kNN rows and the (filtered) obs rows are aligned positionally, which is
  # also how clusterings get written back in find_clusters_sc()
  if (length(clusters) != knn_data$n) {
    stop(sprintf(
      paste(
        "`%s` has %i entries but the kNN graph covers %i cells.",
        "Regenerate the neighbours or the clustering so they align."
      ),
      cluster_col,
      length(clusters),
      knn_data$n
    ))
  }

  if (.verbose) {
    message(sprintf(
      "Running PAGA over the %s kNN graph (%i cells, %i clusters).",
      modality,
      knn_data$n,
      nlevels(clusters)
    ))
  }

  paga_res <- rs_paga(
    knn_mat = get_knn_mat(knn_data),
    partitions = as.integer(clusters) - 1L,
    n_partitions = nlevels(clusters)
  )

  new_paga_res(
    rs_res = paga_res,
    cluster_levels = levels(clusters),
    cluster_col = cluster_col,
    modality = modality
  )
}
