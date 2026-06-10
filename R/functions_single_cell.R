# helpers ----------------------------------------------------------------------

## consts ----------------------------------------------------------------------

# Auto-streaming fires up with more than 100,000 cells
.N_CELLS_STREAMING_THRESHOLD <- 1e5

## utils -----------------------------------------------------------------------

#' Helper to parse the verbosity
#'
#' @param input Boolean or integer to parse
#'
#' @returns The integer controlling the verbosity
#'
#' @keywords internal
parse_verbosity <- function(input) {
  # checks
  checkmate::qassert(input, c("B1", "I1[0, 2]"))

  as.integer(sum(input))
}

#' Helper to set streaming to TRUE on large data sets
#'
#' @param n_cells Integer. Number of cells in the data set.
#' @param streaming Optional boolean. Parsed from the main function.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns Boolean if streaming should be used.
#'
#' @keywords internal
auto_streaming <- function(n_cells, streaming = NULL, .verbose = TRUE) {
  # checks
  checkmate::qassert(n_cells, "I1")
  checkmate::qassert(streaming, c("0", "B1"))

  res <- ifelse(
    test = is.logical(streaming),
    yes = streaming,
    no = n_cells > .N_CELLS_STREAMING_THRESHOLD
  )

  if (res && .verbose) {
    message("Setting streaming for large data to TRUE.")
  }

  return(res)
}

## cell ranger outputs ---------------------------------------------------------

#' Helper to generate cell ranger input parameters
#'
#' @param dir_data String. The directory with the Cell Ranger outputs
#'
#' @return A list based on [bixverse::params_sc_mtx_io()].
#'
#' @export
get_cell_ranger_params <- function(dir_data) {
  # checks
  checkmate::assertDirectory(dir_data)
  assertFileExists(dir_data, c("barcodes.tsv", "genes.tsv", "matrix.mtx"))

  res <- params_sc_mtx_io(
    path_mtx = path.expand(file.path(dir_data, "matrix.mtx")),
    path_obs = path.expand(file.path(dir_data, "barcodes.tsv")),
    path_var = path.expand(file.path(dir_data, "genes.tsv")),
    cells_as_rows = FALSE,
    has_hdr = FALSE
  )

  return(res)
}


## seurat assay to list --------------------------------------------------------

#' Transform Seurat raw counts into a List
#'
#' @param seurat_obj `Seurat` class. The class from which to extract the counts
#' from and transform to a list.
#'
#' @return A list with the following elements
#' \itemize{
#'   \item indptr - Index pointers of the sparse data.
#'   \item indices - Indices of the data.
#'   \item data - The underlying data.
#'   \item format - String that defines if the data is CSR or CSC.
#'   \item nrow - The number of rows.
#'   \item ncol - The number of columns.
#' }
get_seurat_counts_to_list <- function(seurat_obj) {
  # checks
  checkmate::assertClass(seurat_obj, "Seurat")

  assay <- seurat_obj@assays$RNA
  if (methods::.hasSlot(assay, "layers")) {
    raw_counts <- assay@layers$counts
  } else {
    raw_counts <- assay@counts
  }

  res <- list(
    indptr = raw_counts@p,
    indices = raw_counts@i,
    data = raw_counts@x,
    format = "csr",
    nrow = raw_counts@Dim[2],
    ncol = raw_counts@Dim[1]
  )

  return(res)
}

## meta cell matrices ----------------------------------------------------------

#' Create sparse dgRMatrix matrices for raw and norm counts
#'
#' @param meta_cell_data Named list. This contains the indptr, indices and data
#' for both raw counts and norm counts.
#'
#' @returns A list of two items
#' \itemize{
#'  \item raw - Sparse matrix representing the raw counts.
#'  \item norm - Sparse matrix representing the norm counts.
#' }
#'
#' @keywords internal
get_meta_cell_matrices <- function(meta_cell_data) {
  checkmate::assertList(meta_cell_data, names = "named")
  checkmate::assertNames(
    names(meta_cell_data),
    must.include = c(
      "indptr",
      "indices",
      "raw_counts",
      "norm_counts",
      "nrow",
      "ncol"
    )
  )

  dims <- as.integer(c(meta_cell_data$nrow, meta_cell_data$ncol))
  p <- as.integer(meta_cell_data$indptr)
  j <- as.integer(meta_cell_data$indices)

  list(
    raw = Matrix::sparseMatrix(
      p = p,
      j = j,
      x = as.numeric(meta_cell_data$raw_counts),
      dims = dims,
      repr = "R",
      index1 = FALSE
    ),
    norm = Matrix::sparseMatrix(
      p = p,
      j = j,
      x = as.numeric(meta_cell_data$norm_counts),
      dims = dims,
      repr = "R",
      index1 = FALSE
    )
  )
}

## qc --------------------------------------------------------------------------

### per cell outliers ----------------------------------------------------------

#' Use MAD outlier detection on per-cell QC metrics
#'
#' @param metric Numerical vector. The QC metric to check for.
#' @param threshold Numeric. How many MADs in either direction to consider for
#' outlier detection.
#' @param direction String. One of `c("twosided", "below", "above")`. Which
#' directionality to consider
#'
#' @returns A list with:
#' \itemize{
#'  \item outlier - Boolean vector indicating which cell is an outlier
#'  \item metrics - The applied thresholds.
#' }
per_cell_qc_outlier <- function(
  metric,
  threshold = 3,
  direction = c("twosided", "below", "above")
) {
  direction <- match.arg(direction)
  checkmate::qassert(metric, "N+")
  checkmate::qassert(threshold, "N1")

  outliers <- rs_mad_outlier(
    x = metric,
    threshold = threshold,
    direction = direction
  )

  med <- median(metric)
  metrics <- c(median = med, upper_threshold = NA, lower_threshold = NA)

  if (direction %in% c("twosided", "above")) {
    metrics["upper_threshold"] <- med + outliers$threshold
  }
  if (direction %in% c("twosided", "below")) {
    metrics["lower_threshold"] <- med - outliers$threshold
  }

  list(outlier = outliers$outlier, metrics = metrics)
}

### per group outliers ---------------------------------------------------------

#' MAD outlier detection on per-group QC metrics
#'
#' Aggregates each metric to per-group medians and applies `per_cell_qc_outlier`
#' to those medians, flagging whole groups (e.g. donors, samples) whose median
#' deviates by more than `threshold` MADs. Directionality is respected per
#' metric.
#'
#' @param metrics Named list of numeric vectors. The QC metrics.
#' @param groups Character vector. Group label per observation.
#' @param directions Named character vector mapping metric names to direction.
#' One of `"twosided"`, `"below"`, `"above"`.
#' @param threshold Numeric. Number of MADs to use for outlier detection.
#'
#' @return A `data.table` with one row per metric/group and columns `metric`,
#' `group`, `group_median`, `lower_threshold`, `upper_threshold`, `is_outlier`.
#'
#' @keywords internal
per_group_qc_outlier <- function(metrics, groups, directions, threshold = 3) {
  checkmate::assertList(metrics, types = "numeric", names = "unique")
  checkmate::qassert(threshold, "N1")

  stats <- lapply(names(metrics), function(nm) {
    dt <- data.table::data.table(value = metrics[[nm]], group = groups)
    medians <- dt[, .(group_median = median(value)), by = group]
    res <- per_cell_qc_outlier(
      metric = medians$group_median,
      threshold = threshold,
      direction = directions[[nm]]
    )
    data.table::data.table(
      metric = nm,
      group = medians$group,
      group_median = medians$group_median,
      lower_threshold = res$metrics[["lower_threshold"]],
      upper_threshold = res$metrics[["upper_threshold"]],
      is_outlier = res$outlier
    )
  })

  data.table::rbindlist(stats)
}

### multiple metrics -----------------------------------------------------------

#' Run MAD outlier detection on per-cell QC metrics
#'
#' @param metrics Named list of numeric vectors. Each element is a QC metric
#' to check (e.g. `list(log10_lib_size = log10(lib_size), MT = mt_pct)`).
#' @param cells_to_keep Integer. Which cells were included in the analysis.
#' 0-indices for Rust.
#' @param directions Named character vector mapping metric names to direction.
#' One of `"twosided"`, `"below"`, `"above"`. Defaults to `"twosided"` for
#' all metrics if `NULL`.
#' @param threshold Numeric. Number of MADs to use for outlier detection.
#' @param groups Optional grouping variable. An atomic vector of length equal to
#' the metrics. Per-group outlier detection only runs when more than one group
#' is present.
#'
#' @return An object of class `CellQc` containing:
#' \describe{
#'   \item{cell_idx}{Integer vector of 1-indexed cell positions.}
#'   \item{metrics}{The input metrics list.}
#'   \item{groups}{Character vector of group labels.}
#'   \item{per_metric}{Named list of per-metric results from
#'     \code{\link{per_cell_qc_outlier}}.}
#'   \item{outlier_mat}{Logical matrix with one column per metric.}
#'   \item{combined}{Logical vector. `TRUE` if a cell is an outlier in any
#'     metric.}
#'   \item{per_group_stats}{A `data.table` of per-group outlier statistics
#'     (see \code{\link{per_group_qc_outlier}}), or `NULL` for a single group.}
#' }
#'
#' @export
run_cell_qc <- function(
  metrics,
  cells_to_keep,
  directions = NULL,
  threshold = 3,
  groups = NULL
) {
  checkmate::assertList(metrics, types = "numeric", names = "unique")
  checkmate::qassert(cells_to_keep, "I+")
  checkmate::qassert(threshold, "N1")

  n <- length(metrics[[1]])
  if (is.null(groups)) {
    groups <- rep("all", n)
  }
  checkmate::assertAtomic(groups, len = n, any.missing = FALSE)
  groups <- as.character(groups)

  if (is.null(directions)) {
    directions <- setNames(rep("twosided", length(metrics)), names(metrics))
  }

  checkmate::assertNames(names(directions), must.include = names(metrics))

  directions <- directions[names(metrics)]
  group_levels <- unique(groups)

  results <- Map(
    function(x, dir) {
      outlier <- logical(n)
      thresholds <- list()
      for (g in group_levels) {
        idx <- which(groups == g)
        res <- per_cell_qc_outlier(
          metric = x[idx],
          threshold = threshold,
          direction = dir
        )
        outlier[idx] <- res$outlier
        thresholds[[g]] <- res$metrics
      }
      list(outlier = outlier, metrics = thresholds)
    },
    metrics,
    directions
  )

  outlier_mat <- do.call(cbind, lapply(results, `[[`, "outlier"))
  combined <- rowSums(outlier_mat) > 0

  per_group_stats <- if (length(group_levels) > 1L) {
    per_group_qc_outlier(metrics, groups, directions, threshold)
  } else {
    NULL
  }

  structure(
    list(
      cell_idx = cells_to_keep + 1L,
      metrics = metrics,
      groups = groups,
      per_metric = results,
      outlier_mat = outlier_mat,
      combined = combined,
      per_group_stats = per_group_stats
    ),
    class = "CellQc"
  )
}

## knn -------------------------------------------------------------------------

### class ----------------------------------------------------------------------

#' Generate a new SingleCellNearestNeighbour from data
#'
#' @param data Numerical matrix. Samplex x features. The embedding matrix from
#' which to generate the kNN data.
#' @param neighbours_params List. Output of [bixverse::params_sc_neighbours()].
#' A list with the following items:
#' \itemize{
#'   \item full_snn - Boolean. Shall the full shared nearest neighbour graph
#'   be generated that generates edges between all cells instead of between
#'   only neighbours. (Not used in this function.)
#'   \item pruning - Numeric. Weights below this threshold will be set to 0 in
#'   the generation of the sNN graph. (Not used in this function.)
#'   \item snn_similarity - String. One of `c("rank", "jaccard")`. Defines how
#'   the weight from the SNN graph is calculated. For details, please see
#'   [bixverse::params_sc_neighbours()]. (Not used in this function.)
#'   \item knn - List of kNN parameters. See [bixverse::params_knn_defaults()]
#'   for available parameters and their defaults.
#' }
#' @param seed Integer. Random seed for reproducibility.
#' @param .validate_index Boolean. Shall an exhaustive search against a subset
#' of cells be run to validate the approximate nearest neighbour index.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns The `SingleCellNearestNeighbour` for downstream usage.
#'
#' @export
generate_sc_knn <- function(
  data,
  neighbours_params = params_sc_neighbours(),
  seed = 42L,
  .validate_index = FALSE,
  .verbose = TRUE
) {
  checkmate::assertMatrix(data, mode = "numeric")
  assertScNeighbours(neighbours_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.validate_index, "B1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  knn_data <- rs_sc_knn_w_dist(
    embd = data,
    knn_params = neighbours_params,
    verbose = parse_verbosity(.verbose),
    validate_index = .validate_index,
    seed = seed
  )

  knn <- new_sc_knn(knn_data = knn_data, used_cells = row.names(data))

  knn
}

### metrics --------------------------------------------------------------------

#' Calculate recall at k and distance ratio
#'
#' @description
#' Helper function to compare the results of two `SingleCellNearestNeighbour`
#' against each other. The first one can serve as a reference (ground truth)
#' and you can compare against the second one.
#'
#' @param ref_knn The reference `SingleCellNearestNeighbour`.
#' @param query_knn The query `SingleCellNearestNeighbour`.
#'
#' @returns A list with:
#' \itemize{
#'  \item matches - The intersecting indices between the reference and query
#'  kNN for each sample. In an ideal match up should be equal to k.
#'  \item distance_ratio - The distance ratio. Calculates
#'  `sum(dist_query) / sum(dist_ref)` per sample. Indicates how much worse the
#'  reference is.
#'  \item final_recall - The final recall across all samples.
#'  \item final_ratio - The final distance ratio across all samples.
#' }
#'
#' @export
calc_knn_metrics <- function(ref_knn, query_knn) {
  # checks
  checkmate::assertClass(ref_knn, "SingleCellNearestNeighbour")
  checkmate::assertClass(query_knn, "SingleCellNearestNeighbour")

  res <- rs_compare_knn(knn_data_a = ref_knn, knn_data_b = query_knn)

  res
}

## annotation helpers ----------------------------------------------------------

#' Helper function to prepare cell markers
#'
#' @description
#' This function is a helper to generate a list of cell type to cell marker
#' annotations that can be fed into for example `[REF]`.
#'
#' @param obj A single cell class, i.e. one of `SingleCells` or
#' `SingleCellsMultiModal`.
#' @param marker_df A data.table with the columns `"cell_type"` and `"gene_id"`.
#' You need to ensure that the gene id type matches the obj gene id.
#'
#' @returns A list of cell type to marker associations ready for subsequent
#' usage. Genes not found in the object will be automatically removed.
#'
#' @export
prepare_cell_markers <- function(obj, marker_df) {
  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(obj, SingleCells)
  )
  checkmate::assertDataTable(marker_df)
  checkmate::assertNames(
    names(marker_df),
    must.include = c("gene_id", "cell_type")
  )

  marker_df[, gene_idx := get_sc_map(obj)$gene_mapping[marker_df$gene_id] - 1L]

  res <- marker_df[
    !is.na(gene_idx),
    .(
      data = list(list(
        cell_type = cell_type[1],
        positive_indices = as.integer(gene_idx),
        negative_indices = NULL
      ))
    ),
    by = cell_type
  ]

  res_ls <- res$data

  names(res_ls) <- res$cell_type

  return(res_ls)
}

## plotting extraction ---------------------------------------------------------

### helpers --------------------------------------------------------------------

#' Parse an optional modality suffix from a feature id
#'
#' @param feature String. Feature id, optionally suffixed with `_rna` or
#' `_adt`.
#' @param default_modality String. Modality used when no suffix is present.
#'
#' @return A list with `id` and `modality`.
#'
#' @keywords internal
.parse_feature_modality <- function(feature, default_modality) {
  if (grepl("_rna$", feature)) {
    list(id = sub("_rna$", "", feature), modality = "rna")
  } else if (grepl("_adt$", feature)) {
    list(id = sub("_adt$", "", feature), modality = "adt")
  } else {
    list(id = feature, modality = default_modality)
  }
}

### extractors -----------------------------------------------------------------

#' Extract embedding coordinates for plotting
#'
#' @description
#' Pulls an embedding into a long data.table with standardised coordinate
#' columns (`dim_1`, `dim_2`, ...) and, optionally, observation metadata. The
#' embedding name is stored as an `embedding` attribute for axis labelling.
#'
#' @param object A single cell class.
#' @param embedding String. Name of the embedding (e.g. `"umap"`, `"pca"`).
#' @param obs_cols Optional character vector. Obs columns to attach.
#' @param ... Additional arguments forwarded to [get_embedding()] (e.g.
#' `modality`).
#'
#' @return A data.table with `cell_id`, `dim_*` columns and any requested obs
#' columns.
#'
#' @export
extract_embedding_data <- function(object, embedding, obs_cols = NULL, ...) {
  checkmate::qassert(embedding, "S1")
  checkmate::qassert(obs_cols, c("0", "S+"))

  embd <- get_embedding(object, embedding, ...)
  dt <- data.table::as.data.table(embd, keep.rownames = "cell_id")
  data.table::setnames(
    dt,
    setdiff(names(dt), "cell_id"),
    sprintf("dim_%i", seq_len(ncol(embd)))
  )

  if (!is.null(obs_cols)) {
    obs_dt <- object[[obs_cols]]
    for (col in names(obs_dt)) {
      data.table::set(dt, j = col, value = obs_dt[[col]])
    }
  }

  data.table::setattr(dt, "embedding", embedding)

  dt
}

#' Extract per-cell expression mapped onto an embedding
#'
#' @description
#' Combines [extract_gene_expression()] with [extract_embedding_data()] and
#' melts to long format, ready for faceted feature plots. The expression
#' source and the embedding source are chosen independently via
#' `expr_modality` and `embd_modality`, so you can colour an embedding from
#' one modality by expression from another (e.g. RNA expression on an
#' ADT-derived UMAP, or either modality on a WNN embedding). All sources key
#' on the same kept-cell barcodes, so the merge stays aligned regardless of
#' the chosen combination.
#'
#' @param object A single cell class.
#' @param features Character vector. Gene/feature IDs to extract, taken from
#' `expr_modality`.
#' @param embedding String. Name of the embedding.
#' @param scale Boolean. Whether to z-score the expression values.
#' @param clip Optional numeric. Clip z-scores if `scale = TRUE`.
#' @param obs_cols Optional character vector. Obs columns to attach.
#' @param expr_modality String. Modality the expression is pulled from. One of
#' `c("rna", "adt")`.
#' @param embd_modality String. Modality the embedding is pulled from. One of
#' `c("rna", "adt", "wnn")`. Use `"wnn"` for WNN-derived embeddings.
#' @param ... Additional arguments forwarded to [extract_embedding_data()] and
#' onward to [get_embedding()]. Do not pass `modality` here; the embedding
#' modality is set via `embd_modality` and passing it again will error.
#'
#' @return A long data.table with `cell_id`, `dim_*`, `gene` and `expression`.
#'
#' @export
extract_feature_plot_data <- function(
  object,
  features,
  embedding,
  scale = FALSE,
  clip = NULL,
  obs_col = NULL,
  expr_modality = c("rna", "adt"),
  embd_modality = c("rna", "adt", "wnn"),
  ...
) {
  expr_modality <- match.arg(expr_modality)
  embd_modality <- match.arg(embd_modality)

  expr <- extract_gene_expression(
    object = object,
    features = features,
    scale = scale,
    clip = clip,
    modality = expr_modality
  )
  embd <- extract_embedding_data(
    object,
    embedding = embedding,
    modality = embd_modality,
    obs_col = obs_col,
    ...
  )

  dt <- merge(expr, embd, by = "cell_id")
  dim_cols <- grep("^dim_", names(dt), value = TRUE)
  feature_cols <- setdiff(names(expr), "cell_id")
  add_cols <- setdiff(names(dt), c(dim_cols, feature_cols, "cell_id"))

  long <- data.table::melt(
    dt,
    id.vars = c("cell_id", dim_cols, add_cols),
    measure.vars = feature_cols,
    variable.name = "gene",
    value.name = "expression"
  )
  data.table::setattr(long, "embedding", embedding)

  long
}

#' Extract per-cell expression grouped for violin plots
#'
#' @description
#' Combines [extract_gene_expression()] with a grouping obs column and melts to
#' long format, ready for stacked (one gene per row) violin plots.
#'
#' @param object A single cell class.
#' @param features Character vector. Gene IDs to extract.
#' @param grouping_variable String. Obs column to group by.
#' @param scale Boolean. Whether to z-score the expression values.
#' @param clip Optional numeric. Clip z-scores if `scale = TRUE`.
#' @param modality String. One of `c("rna", "adt")`.
#'
#' @return A long data.table with `cell_id`, `group`, `gene` and `expression`.
#' `gene` is an ordered factor following `features`.
#'
#' @export
extract_gene_violin_data <- function(
  object,
  features,
  grouping_variable,
  scale = FALSE,
  clip = NULL,
  modality = c("rna", "adt")
) {
  modality <- match.arg(modality)
  checkmate::qassert(grouping_variable, "S1")

  expr <- extract_gene_expression(
    object = object,
    features = features,
    obs_cols = grouping_variable,
    scale = scale,
    clip = clip,
    modality = modality
  )

  feature_cols <- setdiff(names(expr), c("cell_id", grouping_variable))

  long <- data.table::melt(
    expr,
    id.vars = c("cell_id", grouping_variable),
    measure.vars = feature_cols,
    variable.name = "gene",
    value.name = "expression"
  )

  data.table::setnames(long, grouping_variable, "group")
  long[, group := as.factor(group)]
  long[, gene := factor(gene, levels = feature_cols)]

  long
}

#' Extract a pair of features for scatter / hex plots
#'
#' @description
#' Extracts two features into a wide data.table with `feature_1` and
#' `feature_2` value columns, ready for a scatter or hex plot. Each feature may
#' carry a `_rna` or `_adt` suffix to choose its modality independently (e.g.
#' `"ENSG00000167286_rna"` against `"CD3_adt"`); features without a suffix fall
#' back to `modality`. For `SingleCells` / `MetaCells` only RNA exists, so an
#' `_adt` feature there errors via [extract_gene_expression()].
#'
#' @param object A single cell class.
#' @param feature_1 String. First feature, optionally `_rna` / `_adt` suffixed.
#' @param feature_2 String. Second feature, optionally `_rna` / `_adt` suffixed.
#' @param obs_cols Optional character vector. Obs columns to attach (e.g. to
#' colour the scatter).
#' @param scale Boolean. Whether to z-score the expression values per feature.
#' @param clip Optional numeric. Clip z-scores if `scale = TRUE`.
#' @param modality String. Fallback modality for unsuffixed features. One of
#' `c("rna", "adt")`.
#'
#' @return A data.table with `cell_id`, `feature_1`, `feature_2` and any
#' requested obs columns. The original feature labels are stored in a
#' `features` attribute as `c(feature_1, feature_2)`.
#'
#' @export
extract_feature_pair <- function(
  object,
  feature_1,
  feature_2,
  obs_cols = NULL,
  scale = FALSE,
  clip = NULL,
  modality = c("rna", "adt")
) {
  modality <- match.arg(modality)
  checkmate::qassert(feature_1, "S1")
  checkmate::qassert(feature_2, "S1")
  checkmate::qassert(obs_cols, c("0", "S+"))

  f1 <- .parse_feature_modality(feature_1, modality)
  f2 <- .parse_feature_modality(feature_2, modality)

  # obs ride along on the first extract (kept-order there), then key on cell_id
  expr_1 <- extract_gene_expression(
    object = object,
    features = f1$id,
    obs_cols = obs_cols,
    scale = scale,
    clip = clip,
    modality = f1$modality
  )
  expr_2 <- extract_gene_expression(
    object = object,
    features = f2$id,
    scale = scale,
    clip = clip,
    modality = f2$modality
  )

  data.table::setnames(expr_1, f1$id, "feature_1")
  v2 <- data.table::data.table(
    cell_id = expr_2$cell_id,
    feature_2 = expr_2[[f2$id]]
  )

  dt <- merge(expr_1, v2, by = "cell_id")
  data.table::setcolorder(dt, c("cell_id", "feature_1", "feature_2"))
  data.table::setattr(dt, "features", c(feature_1, feature_2))

  dt
}
