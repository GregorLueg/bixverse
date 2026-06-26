# ------------------------------------------------------------------------------
# this file contains additional, smaller s3 classes relevant for single cell
# analyses. general generics are found on the top of the file, otherwise you
# have per class the class generation and specific methods, getters, setters.
# ------------------------------------------------------------------------------

# additional single cell classes and methods -----------------------------------

## helpers ---------------------------------------------------------------------

#' Add "is_obs" attribute
#'
#' @param x data.table. The data.table to which to add the attribute.
#'
#' @returns `x` with added attribute
#'
#' @keywords internal
.add_is_obs_attr <- function(x) {
  data.table::setattr(x, "is_obs", TRUE)

  return(x)
}

## kNN class -------------------------------------------------------------------

#' Helper function to generate kNN data with distances
#'
#' @description
#' Wrapper class that stores kNN data for subsequent usage in various other
#' functions, methods and helpers.
#'
#' @param knn_data Named list with the following items:
#' \itemize{
#'   \item indices - Integer matrix containing the indices of the nearest
#'   neighbours (0-indexed).
#'   \item dist - Numerical matrix containing the distances to the nearest
#'   neighbours.
#'   \item dist_metric - String. Distance metric used.
#' }
#' @param used_cells Character vector. The cells used to generate the kNN graph
#' with the distances.
#'
#' @return Generates the `SingleCellNearestNeighbour` class.
#'
#' @export
new_sc_knn <- function(knn_data, used_cells) {
  # checks
  checkmate::assertList(knn_data)
  checkmate::assertNames(
    names(knn_data),
    must.include = c("indices", "dist", "dist_metric")
  )
  checkmate::qassert(used_cells, "S+")

  k <- ncol(knn_data$indices)
  n <- nrow(knn_data$indices)

  sc_knn <- list(
    indices = knn_data$indices,
    dist = knn_data$dist,
    k = k,
    n = n,
    dist_metric = knn_data$dist_metric,
    used_cells = used_cells
  )

  class(sc_knn) <- "SingleCellNearestNeighbour"

  return(sc_knn)
}

### primitives -----------------------------------------------------------------

#' @export
print.SingleCellNearestNeighbour <- function(x, ...) {
  cat(
    sprintf("SingleCellNearestNeighbour: %i cells, k = %i\n", x$n, x$k),
    sprintf("  Distance metric: %s\n", x$dist_metric),
    sprintf(
      "  Index range: [%i, %i]\n",
      min(x$indices, na.rm = TRUE),
      max(x$indices, na.rm = TRUE)
    ),
    sprintf(
      "  Distance range: [%.4f, %.4f]\n",
      min(x$dist, na.rm = TRUE),
      max(x$dist, na.rm = TRUE)
    ),
    sep = ""
  )
  invisible(x)
}

### methods --------------------------------------------------------------------

#### getters -------------------------------------------------------------------

#' @rdname get_knn_mat
#'
#' @export
get_knn_mat.SingleCellNearestNeighbour <- function(x) {
  # checks
  checkmate::assertClass(x, "SingleCellNearestNeighbour")

  return(x[["indices"]])
}

#' @rdname get_knn_dist
#'
#' @export
get_knn_dist.SingleCellNearestNeighbour <- function(x) {
  # checks
  checkmate::assertClass(x, "SingleCellNearestNeighbour")

  return(x[["dist"]])
}

#### converstion ---------------------------------------------------------------

#' @export
print.SingleCellNearestNeighbour <- function(x, ...) {
  cat(
    sprintf("SingleCellNearestNeighbour: %i cells, k = %i\n", x$n, x$k),
    sprintf("  Distance metric: %s\n", x$dist_metric),
    sprintf(
      "  Index range: [%i, %i]\n",
      min(x$indices, na.rm = TRUE),
      max(x$indices, na.rm = TRUE)
    ),
    sprintf(
      "  Distance range: [%.4f, %.4f]\n",
      min(x$dist, na.rm = TRUE),
      max(x$dist, na.rm = TRUE)
    ),
    sep = ""
  )
  invisible(x)
}

#' Convert SingleCellNearestNeighbour to manifoldsR NearestNeighbours
#'
#' @description
#' Converts the bixverse single cell kNN class into the manifoldsR
#' `NearestNeighbours` class for use in UMAP, tSNE, etc. Note that the
#' bixverse class stores 0-indexed indices whereas manifoldsR expects
#' 1-indexed, so the conversion handles that.
#'
#' @param x `SingleCellNearestNeighbour` class.
#'
#' @return A `NearestNeighbours` class compatible with manifoldsR.
#'
#' @export
sc_knn_to_nearest_neighbours <- function(x) {
  checkmate::assertClass(x, "SingleCellNearestNeighbour")

  manifoldsR::new_nearest_neighbour(
    indices = as.vector(t(x$indices)) + 1L,
    dist = as.vector(t(x$dist)),
    k = x$k,
    n = x$n
  )
}

## single cell qc results ------------------------------------------------------

### R primitives ---------------------------------------------------------------

#' Print a CellQc object
#'
#' @param x A `CellQc` object.
#' @param ... Ignored.
#'
#' @return Invisible `x`.
#'
#' @export
#'
#' @keywords internal
print.CellQc <- function(x, ...) {
  n_cells <- length(x$combined)
  n_outliers <- sum(x$combined)
  metric_names <- colnames(x$outlier_mat)
  group_levels <- unique(x$groups)
  grouped <- !(length(group_levels) == 1L && group_levels == "all")

  cat(sprintf(
    "CellQc: %d cells, %d outliers (%.1f%%)\n",
    n_cells,
    n_outliers,
    100 * n_outliers / n_cells
  ))
  if (grouped) {
    cat(sprintf("Groups: %d\n", length(group_levels)))
  }
  cat("Metrics:\n")

  for (nm in metric_names) {
    n <- sum(x$outlier_mat[, nm])
    cat(sprintf("  - %s: %d outliers\n", nm, n))
    for (g in group_levels) {
      thresholds <- x$per_metric[[nm]]$metrics[[g]]
      parts <- character()
      if (!is.na(thresholds["lower_threshold"])) {
        parts <- c(
          parts,
          sprintf("lower = %.2f", thresholds["lower_threshold"])
        )
      }
      if (!is.na(thresholds["upper_threshold"])) {
        parts <- c(
          parts,
          sprintf("upper = %.2f", thresholds["upper_threshold"])
        )
      }
      prefix <- if (grouped) sprintf("    [%s] ", g) else "    "
      cat(sprintf("%s%s\n", prefix, paste(parts, collapse = ", ")))
    }
  }

  if (!is.null(x$per_group_stats)) {
    outlier_groups <- x$per_group_stats[(is_outlier)]
    cat(sprintf("Group outliers: %d\n", nrow(outlier_groups)))
    for (i in seq_len(nrow(outlier_groups))) {
      cat(sprintf(
        "  - [%s] %s\n",
        outlier_groups$metric[i],
        outlier_groups$group[i]
      ))
    }
  }

  invisible(x)
}

### getters --------------------------------------------------------------------

#' @rdname get_data
#'
#' @export
get_data.CellQc <- function(x, ...) {
  # checks
  checkmate::assertClass(x, "CellQc")

  # function body
  obs_dt <- data.table::as.data.table(
    x$metrics
  )[, `:=`(
    cell_idx = x$cell_idx,
    grp = x$groups,
    global_outlier = apply(x$outlier_mat, 1, FUN = function(x) {
      any(x)
    })
  )]

  outlier <- data.table::as.data.table(x$outlier_mat)
  colnames(outlier) <- sprintf("%s_is_outlier", colnames(outlier))

  obs_dt <- cbind(obs_dt, outlier)

  obs_dt <- .add_is_obs_attr(obs_dt)

  return(obs_dt)
}


## list results ----------------------------------------------------------------

#' Generate an ScListRes instance
#'
#' @param res Named list. The results from a single cell analysis in list form.
#' @param cell_indices Integer. The cells that were included in the analysis.
#'
#' @returns `ScListRes` class for subsequant usage.
#'
#' @keywords internal
#'
#' @export
new_sc_list <- function(res, cell_indices) {
  # checks
  checkmate::assertList(res, names = "named")
  checkmate::assertInteger(cell_indices, len = length(res[[1]]))

  structure(
    res,
    class = "ScListRes",
    cell_indices = cell_indices
  )
}

### getters --------------------------------------------------------------------

#' @rdname get_data
#'
#' @export
get_data.ScListRes <- function(x, ...) {
  # checks
  checkmate::assertClass(x, "ScListRes")

  # function body
  obs_dt <- data.table::as.data.table(unclass(x))
  obs_dt[, cell_idx := (attr(x, "cell_indices") + 1)] # was zero indexed

  obs_dt <- .add_is_obs_attr(obs_dt)

  return(obs_dt)
}

## matrix results --------------------------------------------------------------

#' Generate an ScMatrixRes instance
#'
#' @param res Numeric matrix. Of shape samples x features.
#' @param cell_indices Integer. The cells that were included in the analysis.
#'
#' @returns `ScMatrixRes` class for subsequant usage.
#'
#' @keywords internal
#'
#' @export
new_sc_matrix <- function(res, cell_indices) {
  # checks
  checkmate::assertMatrix(
    res,
    mode = "numeric",
    row.names = "named",
    col.names = "named"
  )
  checkmate::assertInteger(cell_indices, len = nrow(res))

  structure(
    res,
    class = "ScMatrixRes",
    cell_indices = cell_indices
  )
}

### getters --------------------------------------------------------------------

#' @rdname get_data
#'
#' @export
get_data.ScMatrixRes <- function(x, columns = NULL, ...) {
  # checks
  checkmate::assertClass(x, "ScMatrixRes")
  checkmate::qassert(columns, c("0", "S+"))

  if (!is.null(columns)) {
    x <- x[, columns]
  }

  # function body
  obs_dt <- data.table::as.data.table(unclass(x))
  obs_dt[, cell_idx := (attr(x, "cell_indices") + 1)] # was zero indexed

  obs_dt <- .add_is_obs_attr(obs_dt)

  return(obs_dt)
}

## scrublet --------------------------------------------------------------------

### helpers --------------------------------------------------------------------

#' Update doublet calls and summary statistics for a new Scrublet threshold
#'
#' Recomputes `predicted_doublets`, `z_scores`, and the three rate fields
#' (`detected_doublet_rate`, `detectable_doublet_fraction`,
#' `overall_doublet_rate`) for either the full result or a single group.
#'
#' @param scrublet_res A `ScrubletRes` object.
#' @param threshold Numeric. The new score threshold to apply.
#' @param sample_name Character or `NULL`. If `NULL`, updates apply to the
#' whole object (ungrouped). Otherwise, updates are scoped to the named group.
#' @param .verbose Logical. If `TRUE`, prints updated rate summaries to the
#' console.
#'
#' @return The `ScrubletRes` object with updated doublet calls and summary
#' statistics.
#'
#' @keywords internal
.update_scrublet_threshold <- function(
  scrublet_res,
  threshold,
  sample_name,
  .verbose
) {
  if (is.null(sample_name)) {
    cell_mask <- rep(TRUE, length(scrublet_res$doublet_scores_obs))
    sim_scores <- scrublet_res$doublet_scores_sim
  } else {
    cell_mask <- scrublet_res$cell_groups == sample_name
    sim_scores <- scrublet_res$doublet_scores_sim[[sample_name]]
  }

  obs_scores <- scrublet_res$doublet_scores_obs[cell_mask]
  obs_errors <- scrublet_res$doublet_errors_obs[cell_mask]

  new_preds <- obs_scores > threshold
  new_z <- (obs_scores - threshold) / obs_errors

  scrublet_res$predicted_doublets[cell_mask] <- new_preds
  scrublet_res$z_scores[cell_mask] <- new_z

  detected_rate <- sum(new_preds) / length(new_preds)
  detectable_fraction <- sum(sim_scores > threshold) / length(sim_scores)
  overall_rate <- if (detectable_fraction > 0.01) {
    detected_rate / detectable_fraction
  } else {
    0.0
  }

  if (is.null(sample_name)) {
    scrublet_res$threshold <- threshold
    scrublet_res$detected_doublet_rate <- detected_rate
    scrublet_res$detectable_doublet_fraction <- detectable_fraction
    scrublet_res$overall_doublet_rate <- overall_rate
  } else {
    scrublet_res$threshold[sample_name] <- threshold
    scrublet_res$detected_doublet_rate[sample_name] <- detected_rate
    scrublet_res$detectable_doublet_fraction[sample_name] <- detectable_fraction
    scrublet_res$overall_doublet_rate[sample_name] <- overall_rate
  }

  if (.verbose) {
    prefix <- if (is.null(sample_name)) "" else sprintf("[%s] ", sample_name)
    cat(sprintf(
      "%sDetected doublet rate = %.1f%%\n",
      prefix,
      100 * detected_rate
    ))
    cat(sprintf(
      "%sEstimated detectable doublet fraction = %.1f%%\n",
      prefix,
      100 * detectable_fraction
    ))
    cat(sprintf("%sOverall doublet rate:\n", prefix))
    cat(sprintf("%s  Estimated = %.1f%%\n", prefix, 100 * overall_rate))
  }

  scrublet_res
}

### histogram plotting ---------------------------------------------------------

#' Plot Scrublet score distributions
#'
#' @description
#' Plots histograms of doublet scores for observed transcriptomes and simulated
#' doublets side by side. The threshold is shown as a dashed red vertical line.
#' For grouped results, plots a single group at a time.
#'
#' @param x A `ScrubletRes` object.
#' @param break_number Integer. Number of breaks to use in the histograms.
#' @param for_sample Optional character. For grouped results, the name of the
#' group to plot. Defaults to the first group if `NULL`.
#' @param ... Additional arguments (unused; required by the S3 generic).
#'
#' @return A `patchwork` object with two `ggplot2` histograms.
#'
#' @export
#'
#' @keywords internal
plot.ScrubletRes <- function(x, break_number = 31L, for_sample = NULL, ...) {
  checkmate::assertClass(x, "ScrubletRes")
  checkmate::qassert(break_number, "I1")
  checkmate::qassert(for_sample, c("S1", "0"))

  is_grouped <- isTRUE(attr(x, "grouped"))

  if (is_grouped) {
    groups <- unique(x$cell_groups)
    sample_name <- if (is.null(for_sample)) groups[1] else for_sample

    if (!(sample_name %in% groups)) {
      stop(sprintf(
        "Sample '%s' not found. Available: %s",
        sample_name,
        paste(groups, collapse = ", ")
      ))
    }

    obs_scores <- x$doublet_scores_obs[x$cell_groups == sample_name]
    sim_scores <- x$doublet_scores_sim[[sample_name]]
    thr <- x$threshold[[sample_name]]
    title_suffix <- sprintf(" [%s]", sample_name)
  } else {
    obs_scores <- x$doublet_scores_obs
    sim_scores <- x$doublet_scores_sim
    thr <- x$threshold
    title_suffix <- ""
  }

  obs_plot <- ggplot2::ggplot(
    data.frame(score = obs_scores),
    ggplot2::aes(x = score)
  ) +
    ggplot2::geom_histogram(
      breaks = seq(0, 1, length.out = break_number),
      fill = "grey",
      colour = "darkgrey",
      ggplot2::aes(y = ggplot2::after_stat(density))
    ) +
    ggplot2::geom_vline(
      xintercept = thr,
      linewidth = 0.5,
      colour = "red",
      linetype = "dashed"
    ) +
    ggplot2::labs(
      title = paste0("Observed transcriptomes", title_suffix),
      x = "Doublet score",
      y = "Probability density"
    ) +
    ggplot2::theme_bw()

  sim_plot <- ggplot2::ggplot(
    data.frame(score = sim_scores),
    ggplot2::aes(x = score)
  ) +
    ggplot2::geom_histogram(
      breaks = seq(0, 1, length.out = break_number),
      fill = "grey",
      colour = "darkgrey",
      ggplot2::aes(y = ggplot2::after_stat(density))
    ) +
    ggplot2::geom_vline(
      xintercept = thr,
      linewidth = 0.5,
      colour = "red",
      linetype = "dashed"
    ) +
    ggplot2::labs(
      title = paste0("Simulated doublets", title_suffix),
      x = "Doublet score",
      y = "Probability density"
    ) +
    ggplot2::theme_bw()

  patchwork::wrap_plots(obs_plot, sim_plot, ncol = 2)
}

### readjusting thresholds -----------------------------------------------------

#' Manually readjust Scrublet doublet call thresholds
#'
#' @description
#' Updates doublet calls and associated summary statistics in a `ScrubletRes`
#' object using a user-supplied threshold. Intended for use after visual
#' inspection of the score histograms via [plot.ScrubletRes()].
#'
#' @param scrublet_res A `ScrubletRes` object.
#' @param threshold Numeric in `[0, 1]`. The new threshold to apply.
#' @param for_sample Optional character. For grouped results, the name of the
#' group to update. Defaults to the first group if `NULL`. Ignored for
#' ungrouped results.
#' @param .verbose Logical. If `TRUE`, prints updated rate summaries to the
#' console.
#'
#' @return The `ScrubletRes` object with updated `predicted_doublets`,
#' `z_scores`, `threshold`, `detected_doublet_rate`,
#' `detectable_doublet_fraction`, and `overall_doublet_rate`.
#'
#' @export
call_doublets_manual <- function(
  scrublet_res,
  threshold,
  for_sample = NULL,
  .verbose = TRUE
) {
  checkmate::assertClass(scrublet_res, "ScrubletRes")
  checkmate::qassert(threshold, "N1[0,1]")
  checkmate::qassert(for_sample, c("S1", "0"))
  checkmate::qassert(.verbose, "B1")

  is_grouped <- isTRUE(attr(scrublet_res, "grouped"))

  if (!is_grouped) {
    return(.update_scrublet_threshold(scrublet_res, threshold, NULL, .verbose))
  }

  groups <- unique(scrublet_res$cell_groups)
  sample_name <- if (is.null(for_sample)) groups[1] else for_sample

  if (!(sample_name %in% groups)) {
    stop(sprintf(
      "Sample '%s' not found. Available: %s",
      sample_name,
      paste(groups, collapse = ", ")
    ))
  }

  .update_scrublet_threshold(scrublet_res, threshold, sample_name, .verbose)
}

### obs data -------------------------------------------------------------------

#' @rdname get_data
#'
#' @export
get_data.ScrubletRes <- function(x, ...) {
  # checks
  checkmate::assertClass(x, "ScrubletRes")

  # function body
  obs_dt <- data.table::data.table(
    doublet = x$predicted_doublets,
    doublet_score = x$doublet_scores_obs
  )
  obs_dt[, cell_idx := (attr(x, "cell_indices") + 1)] # was zero indexed

  obs_dt <- .add_is_obs_attr(obs_dt)

  return(obs_dt)
}

### print ----------------------------------------------------------------------

#' Print a ScrubletRes object
#'
#' @param x A `ScrubletRes` object.
#' @param ... Ignored.
#'
#' @return Invisible `x`.
#'
#' @export
#'
#' @keywords internal
print.ScrubletRes <- function(x, ...) {
  n_cells <- length(x$predicted_doublets)
  n_doublets <- sum(x$predicted_doublets)
  is_grouped <- isTRUE(attr(x, "grouped"))

  header <- if (is_grouped) {
    sprintf(
      "ScrubletRes (grouped by '%s', %d groups): %d cells, %d doublets (%.1f%%)\n",
      attr(x, "group_by_col"),
      length(unique(x$cell_groups)),
      n_cells,
      n_doublets,
      100 * n_doublets / n_cells
    )
  } else {
    sprintf(
      "ScrubletRes: %d cells, %d doublets (%.1f%%)\n",
      n_cells,
      n_doublets,
      100 * n_doublets / n_cells
    )
  }
  cat(header)

  if (is_grouped) {
    cat("Per-group:\n")
    for (g in names(x$threshold)) {
      mask <- x$cell_groups == g
      n_g <- sum(mask)
      d_g <- sum(x$predicted_doublets[mask])
      cat(sprintf(
        "  - %s: %d cells, %d doublets (%.1f%%), threshold = %.4f\n",
        g,
        n_g,
        d_g,
        100 * d_g / n_g,
        x$threshold[[g]]
      ))
    }
  } else {
    cat(sprintf("  Threshold:              %.4f\n", x$threshold))
    cat(sprintf(
      "  Detected doublet rate:  %.1f%%\n",
      100 * x$detected_doublet_rate
    ))
    cat(sprintf(
      "  Detectable fraction:    %.1f%%\n",
      100 * x$detectable_doublet_fraction
    ))
    cat(sprintf(
      "  Overall doublet rate:   %.1f%%\n",
      100 * x$overall_doublet_rate
    ))
    cat(sprintf(
      "  Simulated doublets:     %d\n",
      length(x$doublet_scores_sim)
    ))
  }

  extras <- character()
  if (!is.null(x$pca)) {
    extras <- c(extras, "PCA")
  }
  if (!is.null(x$pair_1)) {
    extras <- c(extras, "doublet pairs")
  }
  if (length(extras) > 0) {
    cat(sprintf(
      "  Optional data:          %s\n",
      paste(extras, collapse = ", ")
    ))
  }

  invisible(x)
}

## boost -----------------------------------------------------------------------

### obs data -------------------------------------------------------------------

#' @rdname get_data
#'
#' @export
get_data.BoostRes <- function(x, ...) {
  # checks
  checkmate::assertClass(x, "BoostRes")

  # function body
  obs_dt <- data.table::as.data.table(
    unclass(x)
  )
  obs_dt[, cell_idx := (attr(x, "cell_indices") + 1)] # was zero indexed

  obs_dt <- .add_is_obs_attr(obs_dt)

  return(obs_dt)
}

### print ----------------------------------------------------------------------

#' Print a BoostRes object
#'
#' @param x A `BoostRes` object.
#' @param ... Ignored.
#'
#' @return Invisible `x`.
#'
#' @export
#'
#' @keywords internal
print.BoostRes <- function(x, ...) {
  n_cells <- length(x$doublet)
  n_doublets <- sum(x$doublet)
  is_grouped <- isTRUE(attr(x, "grouped"))

  if (is_grouped) {
    cat(sprintf(
      "BoostRes (grouped by '%s', %d groups): %d cells, %d doublets (%.1f%%)\n",
      attr(x, "group_by_col"),
      length(unique(x$cell_groups)),
      n_cells,
      n_doublets,
      100 * n_doublets / n_cells
    ))
    cat("Per-group:\n")
    for (g in unique(x$cell_groups)) {
      mask <- x$cell_groups == g
      n_g <- sum(mask)
      d_g <- sum(x$doublet[mask])
      cat(sprintf(
        "  - %s: %d cells, %d doublets (%.1f%%)\n",
        g,
        n_g,
        d_g,
        100 * d_g / n_g
      ))
    }
  } else {
    cat(sprintf(
      "BoostRes: %d cells, %d doublets (%.1f%%)\n",
      n_cells,
      n_doublets,
      100 * n_doublets / n_cells
    ))
    cat(sprintf(
      "  Score range: [%.4f, %.4f]\n",
      min(x$doublet_score),
      max(x$doublet_score)
    ))
  }

  invisible(x)
}

## scdblfinder -----------------------------------------------------------------

### obs data -------------------------------------------------------------------

#' @rdname get_data
#'
#' @export
get_data.ScDblFinderRes <- function(x, ...) {
  # checks
  checkmate::assertClass(x, "ScDblFinderRes")

  # results
  obs_dt <- data.table::as.data.table(x[c(
    "predicted_doublets",
    "doublet_score",
    "cxds_scores",
    "weighted",
    "cluster_labels"
  )])
  obs_dt[, cell_idx := (attr(x, "cell_indices") + 1)]

  obs_dt <- .add_is_obs_attr(obs_dt)

  return(obs_dt)
}

### getters --------------------------------------------------------------------

#' Get the feature matrix used for the classifier
#'
#' @param x An object to get the feature matrix from. This will only include
#' the values of the observed cells.
#'
#' @export
get_feature_mat <- function(x) {
  UseMethod("get_feature_mat")
}

#' @rdname get_feature_mat
#'
#' @export
get_feature_mat.ScDblFinderRes <- function(x, ...) {
  checkmate::assertClass(x, "ScDblFinderRes")

  if (is.null(x$features)) {
    warning(paste(
      "You did not extract the features during run of the function.",
      "Returning NULL."
    ))
  }

  return(x$features)
}

#' Get either the cxds or weighted scores
#'
#' @param x An object to get the weighted or cxds scores from.
#'
#' @export
get_scores <- function(x, ..., score_type = c("weighted", "cxds_scores")) {
  UseMethod("get_scores")
}

#' @rdname get_scores
#'
#' @export
get_scores.ScDblFinderRes <- function(
  x,
  ...,
  score_type = c("weighted", "cxds_scores")
) {
  score_type <- match.arg(score_type)

  checkmate::assertClass(x, "ScDblFinderRes")
  checkmate::assertChoice(score_type, c("weighted", "cxds_scores"))

  x[[score_type]]
}

### print ----------------------------------------------------------------------

#' Print a ScDblFinderRes object
#'
#' @param x A `ScDblFinderRes` object.
#' @param ... Ignored.
#'
#' @return Invisible `x`.
#'
#' @export
#'
#' @keywords internal
print.ScDblFinderRes <- function(x, ...) {
  n_cells <- length(x$predicted_doublets)
  n_doublets <- sum(x$predicted_doublets)
  n_clusters <- length(unique(x$cluster_labels))
  features_extracted <- !is.null(x$features)
  is_grouped <- isTRUE(attr(x, "grouped"))

  if (is_grouped) {
    cat(sprintf(
      "ScDblFinderRes (grouped by '%s', %d groups): %d cells, %d doublets (%.1f%%)\n",
      attr(x, "group_by_col"),
      length(unique(x$cell_groups)),
      n_cells,
      n_doublets,
      100 * n_doublets / n_cells
    ))
    cat("Per-group:\n")
    for (g in names(x$threshold)) {
      mask <- x$cell_groups == g
      n_g <- sum(mask)
      d_g <- sum(x$predicted_doublets[mask])
      cat(sprintf(
        "  - %s: %d cells, %d doublets (%.1f%%), threshold = %.4f\n",
        g,
        n_g,
        d_g,
        100 * d_g / n_g,
        x$threshold[[g]]
      ))
    }
    cat(sprintf("  Total clusters:    %d (prefixed)\n", n_clusters))
    cat(sprintf("  Features available: %s\n", features_extracted))
  } else {
    cat(sprintf(
      "ScDblFinderRes: %d cells, %d doublets (%.1f%%)\n",
      n_cells,
      n_doublets,
      100 * n_doublets / n_cells
    ))
    cat(sprintf("  Threshold:        %.4f\n", x$threshold))
    cat(sprintf(
      "  Score range:      [%.4f, %.4f]\n",
      min(x$doublet_score),
      max(x$doublet_score)
    ))
    cat(sprintf("  Final clusters:   %d\n", n_clusters))
    cat(sprintf("  Features available: %s\n", features_extracted))
  }

  invisible(x)
}

## hotspot gene <> gene results ------------------------------------------------

#' Helper function to generate HotSpot data
#'
#' @description
#' Wrapper class that stores the HotSpot
#'
#' @param hotspot_res Named list with the following items:
#' \itemize{
#'   \item z - Pairwise gene-gene matrix with the Z-scores.
#'   \item cor - Pairwise gene-gene matrix with the local correlation
#'   coefficients.
#' }
#' @param used_genes Character vector. The used genes to generate the matrices.
#' @param used_cells Character vector. The used cells to generate the matrices.
#'
#' @return Generates the `sc_hotspot` class.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @keywords internal
new_sc_hotspot_res <- function(hotspot_res, used_genes, used_cells) {
  # checks
  checkmate::assertList(hotspot_res, types = "matrix")
  checkmate::qassert(used_genes, "S+")
  checkmate::qassert(used_cells, "S+")

  # class generation
  sc_hotspot <- list(
    z = {
      hotspot_res$z %>% `rownames<-`(used_genes) %>% `colnames<-`(used_genes)
    },
    cor = {
      hotspot_res$cor %>% `rownames<-`(used_genes) %>% `colnames<-`(used_genes)
    },
    params = list(used_cells = used_cells),
    module_membership = NULL
  )

  class(sc_hotspot) <- "Hotspot"

  return(sc_hotspot)
}

### functions ------------------------------------------------------------------

#### primitives ----------------------------------------------------------------

#' @method print Hotspot
#'
#' @export
#'
#' @keywords internal
print.Hotspot <- function(x, ...) {
  n_genes <- nrow(x$z)
  n_cells <- length(x$params$used_cells)
  has_modules <- !is.null(x$module_membership)

  cat("Hotspot gene-gene local correlation results\n")
  cat(sprintf("  Genes: %d\n", n_genes))
  cat(sprintf("  Cells: %d\n", n_cells))

  if (has_modules) {
    membership <- x$module_membership
    n_assigned <- sum(!is.na(membership$cluster_member))
    n_modules <- length(unique(na.omit(membership$cluster_member)))
    cat(sprintf(
      "  Modules: %d (%d genes assigned, %d unassigned)\n",
      n_modules,
      n_assigned,
      n_genes - n_assigned
    ))
  } else {
    cat("  Modules: not yet computed (see generate_hotspot_membership)\n")
  }

  invisible(x)
}

#### plotting ------------------------------------------------------------------

#' Plot the Hotspot Z-score matrix
#'
#' @description
#' Produces a heatmap of the pairwise gene-gene Z-score matrix. If module
#' membership has been computed via [generate_hotspot_membership()], genes
#' are ordered by module and unassigned genes are excluded. If no membership
#' is available, all genes are shown.
#'
#' @param x A `Hotspot` object.
#' @param max_genes Integer. Maximum number of genes to plot. If the number
#' of genes exceeds this, a random subsample is drawn (stratified by module
#' if membership is available). Set to `NULL` to disable. Default 500.
#' @param seed Integer. Seed for reproducible subsampling.
#' @param ... Further arguments (currently unused).
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @method plot Hotspot
#'
#' @export
#'
#' @keywords internal
plot.Hotspot <- function(x, max_genes = 500L, seed = 42L, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting.")
  }
  checkmate::assertClass(x, "Hotspot")
  checkmate::qassert(max_genes, "I1")
  checkmate::qassert(seed, "I1")

  z_mat <- x$z

  if (!is.null(x$module_membership)) {
    membership <- data.table::copy(x$module_membership)
    membership <- membership[!is.na(cluster_member)]
    data.table::setorder(membership, cluster_member)

    if (!is.null(max_genes) && nrow(membership) > max_genes) {
      set.seed(seed)
      membership <- membership[,
        .SD[sample(.N, min(.N, ceiling(max_genes * .N / nrow(membership))))],
        by = cluster_member
      ]
      data.table::setorder(membership, cluster_member)
    }

    gene_order <- membership$gene_id
  } else {
    gene_order <- rownames(z_mat)

    if (!is.null(max_genes) && length(gene_order) > max_genes) {
      set.seed(seed)
      gene_order <- sample(gene_order, max_genes)
    }
  }

  z_mat <- z_mat[gene_order, gene_order]

  plot_dt <- data.table::as.data.table(
    as.table(z_mat),
    keep.rownames = FALSE
  )
  data.table::setnames(plot_dt, c("gene_1", "gene_2", "z_score"))

  plot_dt[, gene_1 := factor(gene_1, levels = gene_order)]
  plot_dt[, gene_2 := factor(gene_2, levels = rev(gene_order))]

  cap <- quantile(abs(plot_dt$z_score), 0.99, na.rm = TRUE)
  plot_dt[, z_capped := pmax(pmin(z_score, cap), -cap)]

  ggplot2::ggplot(plot_dt, ggplot2::aes(x = gene_1, y = gene_2)) +
    ggplot2::geom_tile(ggplot2::aes(fill = z_capped)) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "Z-score"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL, y = NULL)
}

#### getters -------------------------------------------------------------------

#' @method get_params Hotspot
#'
#' @export
S7::method(get_params, S7::new_S3_class("Hotspot")) <-
  function(object, to_json = FALSE, pretty_json = FALSE) {
    get_params.Hotspot(
      object = object,
      to_json = to_json,
      pretty_json = pretty_json
    )
  }


#' @rdname get_params
#'
#' @export
get_params.Hotspot <- function(
  object,
  to_json = FALSE,
  pretty_json = FALSE
) {
  # Checks
  checkmate::assertClass(
    object,
    "Hotspot"
  )
  checkmate::qassert(to_json, "B1")
  checkmate::qassert(pretty_json, "B1")

  to_ret <- object[["params"]]
  if (to_json) {
    to_ret <- jsonlite::toJSON(to_ret)
  }
  if (to_json && pretty_json) {
    to_ret <- jsonlite::prettify(to_ret)
  }

  return(to_ret)
}

#' Get the hotspot gene membership table
#'
#' @param x The object from which to retrieve the hotspot gene membership
#'
#' @export
get_hotspot_membership <- function(x) {
  UseMethod("get_hotspot_membership")
}


#' @rdname get_hotspot_membership
#'
#' @export
get_hotspot_membership.Hotspot <- function(
  x
) {
  # checks
  checkmate::assertClass(x, "Hotspot")

  x[["module_membership"]]
}

#### methods -------------------------------------------------------------------

#' Identify hotspot gene clusters
#'
#' @param x An object to generate the hotspot gene clusters for.
#' @param fdr_threshold Numeric. The maximum FDR for a given gene-gene local
#' correlation to be included.
#' @param min_size Integer. Minimum cluster size.
#'
#' @export
generate_hotspot_membership <- function(
  x,
  fdr_threshold = 0.05,
  min_size = 10L
) {
  UseMethod("generate_hotspot_membership")
}

#' @rdname generate_hotspot_membership
#'
#' @export
generate_hotspot_membership.Hotspot <- function(
  x,
  fdr_threshold = 0.05,
  min_size = 10L
) {
  # checks
  checkmate::assertClass(x, "Hotspot")
  checkmate::qassert(fdr_threshold, "N1[0, 1]")
  checkmate::qassert(min_size, "I1")

  gene_membership <- data.table::data.table(
    gene_id = rownames(x$z),
    cluster_member = rs_hotspot_cluster_genes(
      x$z,
      fdr_threshold = 0.05,
      10L
    )
  )

  x[["module_membership"]] <- gene_membership

  x
}

## miloR -----------------------------------------------------------------------

#' Generate a new miloR result
#'
#' @description
#' This is a helper class that wraps the miloR results together. The general
#' idea of the approach is to use the kNN graph generated in the single cell
#' data, generate representative neighbourhoods and calculate differential
#' abundances within these neighbourhoods. For further details on the method,
#' please refer to Dann, et al.
#'
#' @param nhoods Sparse dgCMatrix with cells x neighbourhoods.
#' @param sample_counts Integer matrix. Represents neighbourhoods x cells from
#' sample of interest.
#' @param spatial_dist Numeric. The spatial distance between the neighbourhoods
#' to calculate the spatial FDR.
#' @param params Named list. The parameters that were used to generate these
#' results.
#'
#' @returns An `miloR` class that contains the provided data and has
#' subsequent methods to calculate differential abundance statistics.
#'
#' @references Dann, et al., Nat Biotechnol, 2022
#'
#' @keywords internal
new_sc_miloR_res <- function(nhoods, sample_counts, spatial_dist, params) {
  # checks
  checkmate::checkClass(nhoods, "dgCMatrix")
  checkmate::checkMatrix(
    sample_counts,
    row.names = "named",
    col.names = "named",
    mode = "integer"
  )
  checkmate::qassert(spatial_dist, "N+")
  checkmate::assertList(params, names = "named")

  # function body
  sc_milor <- list(
    nhoods = nhoods,
    sample_counts = sample_counts,
    spatial_dist = spatial_dist,
    nhoods_info = NULL,
    model = NULL,
    params = params
  )

  class(sc_milor) <- "miloR"

  return(sc_milor)
}

### getters --------------------------------------------------------------------

#' @method get_params miloR
#'
#' @export
S7::method(get_params, S7::new_S3_class("miloR")) <-
  function(object, to_json = FALSE, pretty_json = FALSE) {
    get_params.miloR(
      object = object,
      to_json = to_json,
      pretty_json = pretty_json
    )
  }

#' @rdname get_params
#'
#' @export
get_params.miloR <- function(
  object,
  to_json = FALSE,
  pretty_json = FALSE
) {
  # Checks
  checkmate::assertClass(
    object,
    "miloR"
  )
  checkmate::qassert(to_json, "B1")
  checkmate::qassert(pretty_json, "B1")

  to_ret <- object[["params"]]
  if (to_json) {
    to_ret <- jsonlite::toJSON(to_ret)
  }
  if (to_json && pretty_json) {
    to_ret <- jsonlite::prettify(to_ret)
  }

  return(to_ret)
}

#' Get the differential abundance results
#'
#' @param x An object from which to get the differential abundance results
#' from.
#'
#' @returns The differential abundance results stored in the object if found.
#'
#' @export
get_differential_abundance_res <- function(x) {
  UseMethod("get_differential_abundance_res")
}

#' @rdname get_differential_abundance_res
#'
#' @export
get_differential_abundance_res.miloR <- function(
  x
) {
  # checks
  checkmate::assertClass(x, "miloR")

  res <- x[["nhoods_info"]]

  if (is.null(res)) {
    warning(paste(
      "No differential abundance results found in x.",
      "Did you run test_nhoods()?",
      "Returning NULL."
    ))
  }

  res
}

#' Get the fitted model
#'
#' @param x An object from which to get the differential abundance results
#' from.
#'
#' @returns The model object, please refer to [edgeR::glmQLFTest()].
#'
#' @export
get_model_fit <- function(x) {
  UseMethod("get_model_fit")
}

#' @rdname get_model_fit
#'
#' @export
get_model_fit.miloR <- function(
  x
) {
  # checks
  checkmate::assertClass(x, "miloR")

  res <- x[["model"]]

  if (is.null(res)) {
    warning(paste(
      "No DGEGLM results found in x.",
      "Did you run test_nhoods()?",
      "Returning NULL."
    ))
  }

  res
}

#' Get the index cells
#'
#' @param x An object from which get the index cells.
#'
#' @returns The indices of the cells in the neighbourhood.
#'
#' @export
get_index_cells <- function(x) {
  UseMethod("get_index_cells")
}

#' @rdname get_index_cells
#'
#' @export
get_index_cells.miloR <- function(x) {
  # checks
  checkmate::assertClass(x, "miloR")

  x[["params"]]$index_cell
}

### methods --------------------------------------------------------------------

#### helpers -------------------------------------------------------------------

#' Spatial FDR correction for neighbourhoods
#'
#' @param nhoods Sparse matrix of cells x neighbourhoods
#' @param pvalues Numeric vector. The p-values.
#' @param weighting String. Weighting scheme, one of
#' `c("k-distance", "graph-overlap")`
#' @param kth_distances Numeric vector. The k-th nearest neighbour distances.
#' Must be supplied if `weighting == "k-distance"`.
#'
#' @return Vector of spatially-corrected FDR values
#'
#' @noRd
spatial_fdr_correction <- function(
  nhoods,
  pvalues,
  weighting = c("k-distance", "graph-overlap"),
  kth_distances = NULL
) {
  weighting <- match.arg(weighting)

  # checks
  checkmate::checkClass(nhoods, "dgCMatrix")
  checkmate::qassert(pvalues, "n+")
  checkmate::assertChoice(weighting, c("k-distance", "graph-overlap"))
  checkmate::qassert(kth_distances, c("0", "N+"))

  # handle NAs
  haspval <- !is.na(pvalues)
  if (!all(haspval)) {
    pvalues <- pvalues[haspval]
  }

  # weights
  if (weighting == "k-distance") {
    if (is.null(kth_distances)) {
      stop("k-distance weighting requires kth.distances")
    }
    t.connect <- kth_distances[haspval]
  } else if (weighting == "graph-overlap") {
    intersect_mat <- Matrix::crossprod(nhoods)
    diag(intersect_mat) <- 0
    t.connect <- Matrix::rowSums(intersect_mat)
  }

  w <- 1 / t.connect
  w[is.infinite(w)] <- 1

  o <- order(pvalues)
  pvalues <- pvalues[o]
  w <- w[o]

  adjp <- numeric(length(o))
  adjp[o] <- rev(cummin(rev(sum(w) * pvalues / cumsum(w))))
  adjp <- pmin(adjp, 1)

  # put NA's back
  if (!all(haspval)) {
    refp <- rep(NA_real_, length(haspval))
    refp[haspval] <- adjp
    adjp <- refp
  }

  adjp
}

#### neighbourhood testing -----------------------------------------------------

#' Test neighbourhoods for differential abundance
#'
#' @description
#' Performs differential abundance testing on single-cell neighbourhoods using
#' edgeR's quasi-likelihood negative binomial framework. The function fits a
#' generalised linear model to neighbourhood cell counts, tests for differential
#' abundance between conditions, and applies spatial FDR correction to account
#' for overlapping neighbourhoods. This implementation follows the approach
#' described in Dann et al., using graph-based neighbourhoods to identify
#' regions of significant compositional changes in single-cell data.
#'
#' @param x `miloR` object for which to run the differential abundance
#' analysis.
#' @param design Formula for the experimental design
#' @param design_df data.frame. Contains the metadata to be used for the
#' generation of the model matrix.
#' @param coef Optional string/integer. For more complex experimental designs,
#' you can specify which coefficient to test. If NULL, tests the last
#' coefficient in the design matrix (typically the main effect of interest).
#' @param norm_method String. Normalisation method to use. One of
#' `c("TMM", "RLE", "logMS")`. Defaults to TMM (trimmed mean of M-values).
#' @param min_mean Numeric. Minimum mean count threshold for filtering
#' neighbourhoods. Neighbourhoods with mean counts below this value are excluded.
#' Defaults to 0 (no filtering).
#' @param robust Logical. If TRUE, uses robust estimation of the quasi-likelihood
#' dispersion. Recommended for datasets with potential outliers. Defaults to TRUE.
#' @param fdr_weighting String. Spatial FDR weighting scheme. One of
#' `c("k-distance", "graph-overlap", "none")`. k-distance uses the distance to
#' the k-th nearest neighbour, graph-overlap uses neighbourhood overlap counts.
#' Defaults to k-distance.
#'
#' @return The `miloR` object with added model and results from the
#' differential abundance analysis.
#'
#' @references Dann et al., 2022, Nat Biotechnol
#'
#' @export
test_nhoods <- function(
  x,
  design,
  design_df,
  coef = NULL,
  norm_method = c("TMM", "RLE", "logMS"),
  min_mean = 0,
  robust = TRUE,
  fdr_weighting = c("k-distance", "graph-overlap", "none")
) {
  UseMethod("test_nhoods")
}

#' @rdname test_nhoods
#'
#' @export
test_nhoods.miloR <- function(
  x,
  design,
  design_df,
  coef = NULL,
  norm_method = c("TMM", "RLE", "logMS"),
  min_mean = 0,
  robust = TRUE,
  fdr_weighting = c("k-distance", "graph-overlap", "none")
) {
  norm_method <- match.arg(norm_method)
  fdr_weighting <- match.arg(fdr_weighting)

  # checks
  checkmate::assertClass(x, "miloR")
  checkmate::assertFormula(design)
  checkmate::assertDataFrame(design_df, row.names = "named")
  checkmate::qassert(coef, c("0", "S1", "X1"))
  checkmate::qassert(min_mean, "N1")
  checkmate::qassert(robust, "B1")
  checkmate::assertChoice(norm_method, c("TMM", "RLE", "logMS"))
  checkmate::assertChoice(
    fdr_weighting,
    c("k-distance", "graph-overlap", "none")
  )

  mm <- stats::model.matrix(design, data = design_df)

  if (ncol(x$sample_counts) != nrow(mm)) {
    stop(
      "Design matrix (",
      nrow(mm),
      ") and sample counts (",
      ncol(x$sample_counts),
      ") dimensions don't match"
    )
  }

  if (any(colnames(x$sample_counts) != rownames(mm))) {
    if (!all(colnames(x$sample_counts) %in% rownames(mm))) {
      stop("Sample names in counts and design matrix don't match")
    }
    warning("Reordering design matrix to match sample counts")
    mm <- mm[colnames(x$sample_counts), ]
  }

  keep_nh <- if (min_mean > 0) {
    rowMeans(x$sample_counts) >= min_mean
  } else {
    rep(TRUE, nrow(x$sample_counts))
  }

  dge <- edgeR::DGEList(
    counts = x$sample_counts[keep_nh, , drop = FALSE],
    lib.size = colSums(x$sample_counts)
  )

  if (norm_method %in% c("TMM", "RLE")) {
    dge <- edgeR::calcNormFactors(dge, method = norm_method)
  }

  dge <- edgeR::estimateDisp(dge, mm)
  fit <- edgeR::glmQLFit(dge, mm, robust = robust, legacy = TRUE)

  if (is.null(coef)) {
    coef <- ncol(mm)
  }

  res <- edgeR::topTags(
    edgeR::glmQLFTest(fit, coef = coef),
    sort.by = "none",
    n = Inf
  ) %>%
    as.data.frame()

  res$Nhood <- which(keep_nh)

  if (fdr_weighting != "none") {
    spatial_fdr <- spatial_fdr_correction(
      nhoods = x$nhoods[, keep_nh, drop = FALSE],
      pvalues = res$PValue,
      weighting = fdr_weighting,
      kth_distances = x$spatial_dist[keep_nh]
    )
    res$SpatialFDR <- spatial_fdr
  } else {
    res$SpatialFDR <- NA_real_
  }

  res <- data.table::setDT(res)
  data.table::setcolorder(
    res,
    c("Nhood", "logFC", "logCPM", "F", "PValue", "FDR", "SpatialFDR")
  )

  x[["nhoods_info"]] <- res
  x[["model"]] <- fit

  return(x)
}

#### neighbourhood annotations -------------------------------------------------

#' Add neighbourhood info on majority cell type
#'
#' @description
#' This function adds cell type composition information to the `nhoods_info`
#' slot within the `miloR` object. For each neighbourhood, it calculates
#' the proportion of the majority cell type and identifies which cell type
#' is most abundant. This is useful for annotating differential abundance
#' results with the cellular composition of each neighbourhood.
#'
#' @param x `miloR` object on which to tag on additional neighbourhood
#' information.
#' @param cell_info Character vector. Represents the cell type annotations
#' you wish to add to the different neighbourhoods. Must be the same length
#' as the number of cells (rows) in the nhoods matrix.
#'
#' @return Modified `miloR` object with updated `nhoods_info` containing
#' `majority_celltype` and `majority_prop` columns.
#'
#' @export
add_nhoods_info <- function(
  x,
  cell_info
) {
  UseMethod("add_nhoods_info")
}

#' @rdname add_nhoods_info
#'
#' @export
add_nhoods_info.miloR <- function(x, cell_info) {
  # checks
  checkmate::assertClass(x, "miloR")
  checkmate::qassert(cell_info, "S+")

  # early return
  if (is.null(x[["nhoods_info"]])) {
    warning(paste(
      "No neighbourhood information found.",
      "Did you run test_nhoods() on the object?",
      "Returning object as is."
    ))
    return(x)
  }

  # check dimensions
  if (length(cell_info) != nrow(x$nhoods)) {
    stop(
      "Length of cell_info (",
      length(cell_info),
      ") doesn't match number of cells in nhoods (",
      nrow(x$nhoods),
      ")"
    )
  }

  nhood_indices <- Matrix::which(x$nhoods != 0, arr.ind = TRUE)

  celltype_counts <- table(
    celltype = cell_info[nhood_indices[, 1]],
    nhood = nhood_indices[, 2]
  ) %>%
    t() %>%
    unclass() %>%
    as.matrix()

  # calculate majority cell type and proportion for each neighbourhood
  nhood_info <- data.table::data.table(
    Nhood = as.integer(rownames(celltype_counts))
  )

  nhood_info[, `:=`(
    majority_celltype = colnames(celltype_counts)[apply(
      celltype_counts,
      1,
      which.max
    )],
    majority_prop = apply(celltype_counts, 1, max) / rowSums(celltype_counts)
  )]

  x[["nhoods_info"]] <- x[["nhoods_info"]][
    nhood_info,
    on = "Nhood",
    `:=`(
      majority_celltype = i.majority_celltype,
      majority_prop = i.majority_prop
    )
  ]

  x
}

## scenic grns -----------------------------------------------------------------

#' Constructor for SCENIC GRN results
#'
#' @description
#' Stores the TF-gene importance matrix and associated metadata from a SCENIC
#' GRN inference run. Intended as input to downstream regulon generation
#' functions.
#'
#' @param importance_matrix Matrix. Importance matrix of shape
#' `(n_genes, n_tfs)` with gene identifiers as rownames and TF identifiers
#' as colnames.
#' @param gene_ids Character vector. Gene identifiers corresponding to rows.
#' @param tf_ids Character vector. TF identifiers corresponding to columns.
#' @param params List. The full SCENIC parameters used for the run.
#'
#' @return An object of class `ScenicGrn`.
#'
#' @export
#'
#' @keywords internal
new_scenic_grn <- function(
  importance_matrix,
  gene_ids,
  tf_ids,
  params
) {
  checkmate::assertMatrix(importance_matrix, mode = "numeric")
  checkmate::qassert(gene_ids, "S+")
  checkmate::qassert(tf_ids, "S+")
  checkmate::assertList(params)
  checkmate::assertTRUE(nrow(importance_matrix) == length(gene_ids))
  checkmate::assertTRUE(ncol(importance_matrix) == length(tf_ids))

  scenic_grn <- list(
    importance_matrix = importance_matrix,
    gene_ids = gene_ids,
    tf_ids = tf_ids,
    params = params,
    tf_to_gene_results = data.table(),
    cis_targets_results = data.table()
  )

  class(scenic_grn) <- "ScenicGrn"
  return(scenic_grn)
}

### primitives -----------------------------------------------------------------

#### transforms ----------------------------------------------------------------

#' @export
as.matrix.ScenicGrn <- function(x, ...) {
  x$importance_matrix
}

#### prints --------------------------------------------------------------------

#' Print a ScenicGrn object
#'
#' @param x A `ScenicGrn` object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @returns Invisibly returns `x`.
#'
#' @export
#'
#' @keywords internal
print.ScenicGrn <- function(x, ...) {
  tf_to_gene_generated <- nrow(x$tf_to_gene_results) > 0
  cis_targets_results_generated <- nrow(x$cis_targets_results) > 0

  cat("ScenicGrn (GRN results)\n")
  cat("  No genes:                ", length(x$gene_ids), "\n")
  cat("  No TFs:                  ", length(x$tf_ids), "\n")
  cat("  Applied learner:         ", x$params$learner_type, "\n")
  cat("  TF to gene generated:    ", tf_to_gene_generated, "\n")
  cat("  CisTarget res generated: ", cis_targets_results_generated, "\n")

  invisible(x)
}

### getters --------------------------------------------------------------------

#' @method get_params ScenicGrn
#'
#' @export
S7::method(get_params, S7::new_S3_class("ScenicGrn")) <-
  function(object, to_json = FALSE, pretty_json = FALSE) {
    get_params.ScenicGrn(
      object = object,
      to_json = to_json,
      pretty_json = pretty_json
    )
  }

#' @rdname get_params
#'
#' @export
get_params.ScenicGrn <- function(
  object,
  to_json = FALSE,
  pretty_json = FALSE
) {
  # Checks
  checkmate::assertClass(
    object,
    "ScenicGrn"
  )
  checkmate::qassert(to_json, "B1")
  checkmate::qassert(pretty_json, "B1")

  to_ret <- object[["params"]]
  if (to_json) {
    to_ret <- jsonlite::toJSON(to_ret)
  }
  if (to_json && pretty_json) {
    to_ret <- jsonlite::prettify(to_ret)
  }

  return(to_ret)
}

#' Extract the TF to gene data from the ScenicGrn object
#'
#' @param x `ScenicGrn` object from which to extract the TF to gene data.table.
#'
#' @returns data.table with TF to gene information
#'
#' @export
get_tf_to_gene <- function(x) {
  UseMethod("get_tf_to_gene")
}

#' @rdname get_tf_to_gene
#'
#' @export
get_tf_to_gene.ScenicGrn <- function(x) {
  # checks
  checkmate::assertClass(x, "ScenicGrn")

  # return a proper copy
  tf_to_gene <- data.table::copy(x[["tf_to_gene_results"]])

  if (nrow(tf_to_gene) == 0) {
    warning(paste(
      "You did not seem to have run identify_tf_to_genes().",
      "Returning empty data.table"
    ))
  }

  return(tf_to_gene)
}

#' Extract the TF to gene data from the ScenicGrn object
#'
#' @param x `ScenicGrn` object from which to extract the TF to gene data.table.
#'
#' @returns data.table with TF to gene information
#'
#' @export
get_cistarget_res <- function(x) {
  UseMethod("get_cistarget_res")
}

#' @rdname get_cistarget_res
#'
#' @export
get_cistarget_res.ScenicGrn <- function(x) {
  # checks
  checkmate::assertClass(x, "ScenicGrn")

  # return a proper copy
  cis_targets_results <- data.table::copy(x[["cis_targets_results"]])

  if (nrow(cis_targets_results) == 0) {
    warning(paste(
      "You did not seem to have run identify_tf_to_genes().",
      "Returning empty data.table"
    ))
  }

  return(cis_targets_results)
}

### generate tf to gene --------------------------------------------------------

### generate tf to gene --------------------------------------------------------
#' Identify the TF to gene regulation
#'
#' @description
#' Generates the TF to gene associations from the importance matrix produced
#' by the tree-based regression models. Two filtering strategies are available:
#'
#' \itemize{
#'   \item **Threshold** (`method = "threshold"`): For each gene, computes
#'     mean + `n_sd` * SD of the importance scores across all TFs and retains
#'     only pairs exceeding that threshold. This adapts to the per-gene
#'     importance distribution and is less sensitive to differences between
#'     learners.
#'   \item **Top-k** (`method = "top_k"`): Selects the top `k_tfs` TFs per
#'     gene and/or the top `k_genes` genes per TF. Both margins can be combined
#'     (union). At least one of `k_tfs` or `k_genes` must be provided.
#' }
#' Both methods accept an optional `min_importance` floor.
#'
#' @param x `ScenicGrn` object.
#' @param method Character. Either `"top_k"` or `"threshold"`.
#' @param k_tfs Optional integer. Top TFs per gene (only used when
#' `method = "top_k"`).
#' @param k_genes Optional integer. Top genes per TF (only used when
#' `method = "top_k"`).
#' @param n_sd Numeric. Number of standard deviations above the per-gene mean
#' to use as the threshold (only used when `method = "threshold"`). Default
#' is `2`.
#' @param min_importance Optional numeric in \[0, 1\]. Absolute minimum
#'   importance score for inclusion.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns The `ScenicGrn` object with a `tf_to_gene_results` data.table
#'   added.
#'
#' @export
identify_tf_to_genes <- function(
  x,
  method = c("threshold", "top_k"),
  k_tfs = NULL,
  k_genes = NULL,
  n_sd = 2,
  min_importance = NULL,
  .verbose = TRUE
) {
  UseMethod("identify_tf_to_genes")
}

#' @rdname identify_tf_to_genes
#'
#' @export
identify_tf_to_genes.ScenicGrn <- function(
  x,
  method = c("threshold", "top_k"),
  k_tfs = NULL,
  k_genes = NULL,
  n_sd = 2,
  min_importance = NULL,
  .verbose = TRUE
) {
  method <- match.arg(method)

  checkmate::assertClass(x, "ScenicGrn")
  checkmate::qassert(n_sd, "N1(0,)")
  checkmate::qassert(min_importance, c("N1[0, 1]", "0"))
  checkmate::qassert(.verbose, "B1")

  gene_tf_imp <- x$importance_matrix

  if (method == "top_k") {
    checkmate::qassert(k_tfs, c("I1", "0"))
    checkmate::qassert(k_genes, c("I1", "0"))
    if (is.null(k_tfs) && is.null(k_genes)) {
      stop("k_tfs and k_genes cannot both be NULL when method = 'top_k'.")
    }

    if (.verbose) {
      message("Extracting TF to gene associations via top-k filtering.")
    }

    parts <- list()
    if (!is.null(k_tfs)) {
      parts[[length(parts) + 1L]] <- data.table::as.data.table(
        rs_top_k_targets(gene_tf_imp, k_tfs, 1L, min_importance)
      )
    }
    if (!is.null(k_genes)) {
      parts[[length(parts) + 1L]] <- data.table::as.data.table(
        rs_top_k_targets(gene_tf_imp, k_genes, 2L, min_importance)
      )
    }
    tf_gene_dt <- data.table::rbindlist(parts)
  } else {
    if (.verbose) {
      message(
        sprintf(
          "Extracting TF to gene associations via per-gene threshold (mean + %.1f * SD).",
          n_sd
        )
      )
    }

    tf_gene_dt <- data.table::as.data.table(
      rs_importance_threshold(
        matrix = gene_tf_imp,
        n_sd = n_sd,
        min_value = min_importance
      )
    )
  }

  tf_gene_dt <- unique(tf_gene_dt, by = c("tf", "gene"))
  x$tf_to_gene_results <- tf_gene_dt
  return(x)
}

### tf to gene correlation -----------------------------------------------------

#### helpers -------------------------------------------------------------------

#' Get TF gene correlations for SingleCells
#'
#' @param tf_to_gene data.table. The TF to gene data.table
#' @param object `SingleCells` class.
#' @param spearman Boolean. Shall the Spearman correlation be used.
#'
#' @returns Adds a `pairwise_cor` column to the tf_to_gene data.table
#'
#' @keywords internal
.tf_gene_cor_sc <- function(tf_to_gene, object, spearman) {
  # checks
  checkmate::assertDataTable(tf_to_gene)
  checkmate::assertNames(
    names(tf_to_gene),
    must.include = c("tf", "gene", "importance")
  )
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(spearman, "B1")

  indices_1 <- get_sc_map(object)$gene_mapping[tf_to_gene$tf] - 1L
  indices_2 <- get_sc_map(object)$gene_mapping[tf_to_gene$gene] - 1L

  pairwise_cors <- rs_pairwise_gene_cors(
    f_path = get_rust_count_gene_f_path(object),
    gene_indices_1 = as.integer(indices_1),
    gene_indices_2 = as.integer(indices_2),
    cells_to_keep = as.integer(get_cells_to_keep(object)),
    spearman = FALSE
  )

  tf_to_gene[, pairwise_cor := pairwise_cors]

  return(tf_to_gene)
}

#' Get TF gene correlations for MetaCells
#'
#' @param tf_to_gene data.table. The TF to gene data.table
#' @param object `MetaCells` class.
#' @param spearman Boolean. Shall the Spearman correlation be used.
#'
#' @returns Adds a `pairwise_cor` column to the tf_to_gene data.table
#'
#' @keywords internal
.tf_gene_cor_mc <- function(tf_to_gene, object, spearman) {
  # checks
  checkmate::assertDataTable(tf_to_gene)
  checkmate::assertNames(
    names(tf_to_gene),
    must.include = c("tf", "gene", "importance")
  )
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::qassert(spearman, "B1")

  indices_1 <- get_gene_indices(
    x = object,
    gene_ids = tf_to_gene$tf,
    rust_index = TRUE
  )
  indices_2 <- get_gene_indices(
    x = object,
    gene_ids = tf_to_gene$gene,
    rust_index = TRUE
  )

  count_list <- mc_counts_to_list(
    object = object,
    assay = "norm"
  )

  pairwise_cors <- rs_pairwise_gene_cors_mc(
    sparse_data = count_list,
    gene_indices_1 = indices_1,
    gene_indices_2 = indices_2,
    spearman = spearman
  )

  tf_to_gene[, pairwise_cor := pairwise_cors]

  return(tf_to_gene)
}

#### main function -------------------------------------------------------------

#' Generate TF to gene correlations
#'
#' @description
#' This function will calculate the correlations between the identified TF to
#' gene pairs. You need to have run [identify_tf_to_genes()]!
#'
#' @param x `ScenicGrn` object for which to generate the TF to gene
#' associations.
#' @param object `SingleCells` or `MetaCells` object that was used to generate
#' the original GRNs.
#' @param cor_filter Optional float. If you wish to filter out TF genes below
#' a certain correlation. If `NULL` all genes will be kept.
#' @param remove_self Boolean. Shall self loops (where TF controls its own
#' expression) be removed. Defaults to `TRUE`.
#' @param spearman Boolean. Shall Spearman correlation be used. Defaults to
#' `TRUE`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns Adds the correlations coefficients between to the TF to gene.
#'
#' @export
tf_to_genes_correlations <- function(
  x,
  object,
  cor_filter = NULL,
  remove_self = TRUE,
  spearman = TRUE,
  .verbose = TRUE
) {
  UseMethod("tf_to_genes_correlations")
}

#' @rdname tf_to_genes_correlations
#'
#' @export
tf_to_genes_correlations.ScenicGrn <- function(
  x,
  object,
  cor_filter = NULL,
  remove_self = TRUE,
  spearman = TRUE,
  .verbose = TRUE
) {
  # checks
  checkmate::assertClass(x, "ScenicGrn")
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) || S7::S7_inherits(object, MetaCells)
  )
  checkmate::qassert(cor_filter, c("0", "N1"))
  checkmate::qassert(remove_self, "B1")
  checkmate::qassert(spearman, "B1")
  checkmate::qassert(.verbose, "B1")

  # early return
  if (nrow(x$tf_to_gene_results) == 0) {
    warning(paste(
      "No TF to gene pairs found. Returning class as is.",
      "Did you run identify_tf_to_genes()?"
    ))
    return(x)
  }

  tf_to_gene <- get_tf_to_gene(x)

  if (.verbose) {
    message("Calculating the pairwise correlations between the TFs and genes")
  }

  tf_to_gene <- if (S7::S7_inherits(object, SingleCells)) {
    .tf_gene_cor_sc(tf_to_gene, object, spearman)
  } else {
    .tf_gene_cor_mc(tf_to_gene, object, spearman)
  }

  if (!is.null(cor_filter)) {
    if (.verbose) {
      message(sprintf(
        "Removing TF <> gene pairs with cors <= %.3f",
        cor_filter
      ))
    }
    tf_to_gene <- tf_to_gene[pairwise_cor >= cor_filter]
  }

  if (remove_self) {
    if (.verbose) {
      message("Removing self loops (TF controlling its own expression")
    }
    tf_to_gene <- tf_to_gene[tf != gene]
  }

  x[["tf_to_gene_results"]] <- tf_to_gene

  return(x)
}

### scenic cistarget -----------------------------------------------------------

#' Run the SCENIC motif enrichment
#'
#' @description
#' This function will run the motif enrichment based on RCisTarget (or in this
#' case the internal implementation via [run_cistarget()]). You need to provide
#' the expected rankings and TF annotations (for details, please see
#' [run_cistarget()]). Briefly, this function will run CisTarget and add the
#' results to the `ScenicGrn` object and add additionally a column `"in_motif"`,
#' for a given TF to gene set to say if it was part of the motifs associated
#' with this TF (or not). You have the option to limit this to only the
#' the high confidence TFs (default), or also include the low confidence TFs
#' (i.e., links from TF to motif that are less certain).
#'
#' @param x `ScenicGrn` object for which to generate the TF to gene
#' associations.
#' @param motif_rankings Integer matrix. Motif rankings for genes. Row names are
#' gene identifiers, column names are motif identifiers. Lower values indicate
#' higher regulatory potential.
#' @param annot_data data.table. Motif annotation database mapping motifs to
#' transcription factors. Must contain columns: `motif`, `TF`, and
#' `annotationSource`.
#' @param cis_target_params List. Output of [bixverse::params_cistarget()]:
#' \itemize{
#'   \item{auc_threshold - Numeric. Proportion of genes to use for AUC
#'   threshold calculation. Default 0.05 means top 5 percent of genes.}
#'   \item{nes_threshold - Numeric. Normalised Enrichment Score threshold for
#'   determining significant motifs. Default is 3.0.}
#'   \item{rcc_method - Character. Recovery curve calculation method: "approx"
#'   (approximate, faster) or "icistarget" (exact, slower).}
#'   \item{high_conf_cats - Character vector. Annotation categories considered
#'   high confidence (e.g., "directAnnotation", "inferredBy_Orthology").}
#'   \item{low_conf_cats - Character vector. Annotation categories considered
#'   lower confidence (e.g., "inferredBy_MotifSimilarity").}
#' }
#' @param gene_id_to_symbol Named character vector. Mapping from gene
#' identifiers used internally (e.g., Ensembl IDs) to the identifiers used
#' in the motif rankings (e.g., HGNC symbols). Names are internal IDs,
#' values are ranking IDs. If `NULL` (default), no mapping is applied and
#' the internal IDs are assumed to match the ranking rownames.
#' @param only_high_conf_tf Boolean. Shall only the high confidence TF to
#' motif association be used. Defaults to `TRUE`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns Adds a data.table with the first tf to gene results to the class.
#'
#' @export
tf_to_genes_motif_enrichment <- function(
  x,
  motif_rankings,
  annot_data,
  cis_target_params = params_cistarget(),
  gene_id_to_symbol = NULL,
  only_high_conf_tf = TRUE,
  .verbose = TRUE
) {
  UseMethod("tf_to_genes_motif_enrichment")
}

#' @rdname tf_to_genes_motif_enrichment
#'
#' @export
tf_to_genes_motif_enrichment.ScenicGrn <- function(
  x,
  motif_rankings,
  annot_data,
  cis_target_params = params_cistarget(),
  gene_id_to_symbol = NULL,
  only_high_conf_tf = TRUE,
  .verbose = TRUE
) {
  checkmate::assertClass(x, "ScenicGrn")
  checkmate::assertMatrix(
    motif_rankings,
    mode = "integer",
    row.names = "named",
    col.names = "named"
  )
  checkmate::assertDataTable(annot_data)
  checkmate::assertNames(
    names(annot_data),
    must.include = c("motif", "TF", "annotationSource")
  )
  assertCistargetParams(cis_target_params)
  checkmate::qassert(only_high_conf_tf, "B1")
  checkmate::qassert(.verbose, "B1")
  checkmate::assertCharacter(
    gene_id_to_symbol,
    names = "named",
    null.ok = TRUE
  )

  tf_gene_dt <- get_tf_to_gene(x)

  if (!is.null(gene_id_to_symbol)) {
    mapped_tf <- gene_id_to_symbol[tf_gene_dt$tf]
    mapped_gene <- gene_id_to_symbol[tf_gene_dt$gene]
    unmapped <- is.na(mapped_tf) | is.na(mapped_gene)
    if (any(unmapped) && .verbose) {
      message(sprintf(
        "  Dropping %d TF-gene pairs with no symbol mapping.",
        sum(unmapped)
      ))
    }
    tf_gene_dt_mapped <- tf_gene_dt[!unmapped]
    tf_gene_dt_mapped[, `:=`(
      tf_symbol = gene_id_to_symbol[tf],
      gene_symbol = gene_id_to_symbol[gene]
    )]
    tf_gene_lists <- split(
      tf_gene_dt_mapped$gene_symbol,
      tf_gene_dt_mapped$tf_symbol
    )
  } else {
    tf_gene_lists <- split(tf_gene_dt$gene, tf_gene_dt$tf)
    tf_gene_dt_mapped <- tf_gene_dt
  }

  if (.verbose) {
    message(paste(
      "Running CisTarget to associate TF to gene pairs",
      "with motif enrichment info."
    ))
  }

  cis_res <- run_cistarget(
    gs_list = tf_gene_lists,
    rankings = motif_rankings,
    annot_data = annot_data,
    cis_target_params = cis_target_params
  )
  x[["cis_targets_results"]] <- cis_res

  # helper: explode a semicolon-separated TF column, keeping only rows where
  # gs_name matches one of the TFs, then explode leading edge genes
  # build reverse mapping once if needed
  if (!is.null(gene_id_to_symbol)) {
    symbol_to_id <- setNames(names(gene_id_to_symbol), gene_id_to_symbol)
  }

  explode_leading_edge <- function(dt, tf_col) {
    dt <- dt[!is.na(get(tf_col))]
    if (nrow(dt) == 0L) {
      return(data.table::data.table(tf = character(), le_gene = character()))
    }
    le <- dt[,
      .(single_tf = unlist(strsplit(get(tf_col), ";"))),
      by = .(gs_name, leading_edge_genes)
    ][
      single_tf == gs_name,
      .(le_gene = unlist(strsplit(leading_edge_genes, ";"))),
      by = .(tf = gs_name)
    ]

    # map back to internal IDs if a mapping was provided
    if (!is.null(gene_id_to_symbol)) {
      le[, tf := symbol_to_id[tf]]
      le[, le_gene := symbol_to_id[le_gene]]
      le <- le[!is.na(tf) & !is.na(le_gene)]
    }

    le
  }

  le_high <- explode_leading_edge(cis_res, "TF_highConf")

  tf_gene_dt[, in_leading_edge := FALSE]
  tf_gene_dt[le_high, in_leading_edge := TRUE, on = .(tf, gene = le_gene)]

  if (!only_high_conf_tf) {
    if (.verbose) {
      message(" Adding also low confidence TF to motif info")
    }

    le_low <- explode_leading_edge(cis_res, "TF_lowConf")
    tf_gene_dt[le_low, in_leading_edge := TRUE, on = .(tf, gene = le_gene)]
  }

  data.table::setorder(tf_gene_dt, tf)
  x$tf_to_gene_results <- tf_gene_dt
  return(x)
}

## fast clusters ---------------------------------------------------------------

### getters --------------------------------------------------------------------

#' @rdname get_data
#'
#' @export
get_data.SingleCellFastClusters <- function(x, ...) {
  checkmate::assertClass(x, "SingleCellFastClusters")
  res <- data.table::copy(x$memberships)

  res <- .add_is_obs_attr(res)

  res
}

#' Get k-means centroids from a fast cluster result
#'
#' @param x `SingleCellFastClusters` object.
#'
#' @export
get_centroids <- function(x) {
  UseMethod("get_centroids")
}

#' @rdname get_centroids
#'
#' @export
get_centroids.SingleCellFastClusters <- function(x) {
  checkmate::assertClass(x, "SingleCellFastClusters")
  if (is.null(x$centroids)) {
    warning(paste(
      "No centroids stored. Did you set return_kmeans = TRUE?",
      "Returning NULL."
    ))
  }
  x$centroids
}

#' Get k-means cluster assignments from a fast cluster result
#'
#' @param x `SingleCellFastClusters` object.
#'
#' @export
get_kmeans_clusters <- function(x) {
  UseMethod("get_kmeans_clusters")
}

#' @rdname get_kmeans_clusters
#'
#' @export
get_kmeans_clusters.SingleCellFastClusters <- function(x) {
  checkmate::assertClass(x, "SingleCellFastClusters")
  if (is.null(x$k_means_cluster)) {
    warning(paste(
      "No k-means clusters stored. Did you set return_kmeans = TRUE?",
      "Returning NULL."
    ))
  }
  x$k_means_cluster
}

### primitives -----------------------------------------------------------------

#' @export
#'
#' @keywords internal
print.SingleCellFastClusters <- function(x, ...) {
  cat(sprintf(
    "SingleCellFastClusters: %d cells, %d resolutions\n",
    nrow(x$memberships),
    length(x$resolutions)
  ))
  cat(sprintf("  Resolutions: %s\n", paste(x$resolutions, collapse = ", ")))
  cat(sprintf("  Grid stats stored: %s\n", !is.null(x$stats)))
  cat(sprintf("  k-means stored:    %s\n", !is.null(x$centroids)))

  invisible(x)
}

## sctypes ---------------------------------------------------------------------

### getters --------------------------------------------------------------------

#' Get the ScType score matrix
#'
#' @param x `ScTypeResults` object.
#'
#' @returns A numeric matrix of cells x cell types.
#'
#' @export
get_scores <- function(x) {
  UseMethod("get_scores")
}

#' @rdname get_scores
#'
#' @export
get_scores.ScTypeResults <- function(x) {
  checkmate::assertClass(x, "ScTypeResults")
  m <- matrix(
    x$scores,
    nrow = x$n_cells,
    ncol = x$n_cell_types,
    byrow = TRUE
  )
  colnames(m) <- x$cell_types
  m
}

### cluster scoring ------------------------------------------------------------

#' Score clusters based on ScType
#'
#' @param x `ScTypeResults` object.
#' @param cluster_labels Integer vector. Cluster assignment, of length of the
#' scored cells.
#'
#' @returns A `data.table` with cluster_id, cell_type, scores and n_cells.
#'
#' @export
score_clusters <- function(x, cluster_labels) {
  UseMethod("score_clusters")
}

#' @rdname score_clusters
#'
#' @export
score_clusters.ScTypeResults <- function(x, cluster_labels) {
  checkmate::assertClass(x, "ScTypeResults")
  checkmate::assertIntegerish(cluster_labels, len = x$n_cells)

  res <- rs_sc_type_cluster_assignment(
    sc_type_res = x,
    cluster_labels = as.integer(cluster_labels)
  )

  res <- data.table::setDT(res)

  return(res)
}

### primitives -----------------------------------------------------------------

#' @export
#'
#' @keywords internal
print.ScTypeResults <- function(x, ...) {
  cat(sprintf(
    "ScTypeResults: %d cells, %d cell types\n",
    x$n_cells,
    x$n_cell_types
  ))
  cat(sprintf("  Cell types: %s\n", paste(x$cell_types, collapse = ", ")))

  invisible(x)
}

## nmf result classes ----------------------------------------------------------

### NmfResult ------------------------------------------------------------------

#' Constructor for single-run NMF results
#'
#' @description
#' Stores the W (genes x components) and H (components x cells) matrices from
#' a single HALS NMF run, plus convergence metadata and provenance information.
#'
#' @param nmf_res Named list. Output of [rs_nmf_single_sc()] or
#' [rs_nmf_single_mc()]. Must contain `w`, `h`, `final_loss`, `n_iter`,
#' `converged`.
#' @param gene_ids Character vector. Gene identifiers for the rows of W.
#' @param cell_ids Character vector. Cell (or meta cell) identifiers for the
#' columns of H.
#' @param cell_indices Integer vector. 0-indexed cell positions used in the
#' run (Rust-side indices, kept for re-running / cross-referencing).
#' @param source_class String. One of `c("SingleCells", "MetaCells")`.
#' @param params List. The full set of parameters used for the run.
#'
#' @returns An object of class `NmfResult`.
#'
#' @export
#'
#' @keywords internal
new_nmf_result <- function(
  nmf_res,
  gene_ids,
  cell_ids,
  cell_indices,
  source_class,
  params
) {
  checkmate::assertList(nmf_res)
  checkmate::assertNames(
    names(nmf_res),
    must.include = c("w", "h", "final_loss", "n_iter", "converged")
  )
  checkmate::qassert(gene_ids, "S+")
  checkmate::qassert(cell_ids, "S+")
  checkmate::qassert(cell_indices, "I+")
  checkmate::assertChoice(source_class, c("SingleCells", "MetaCells"))
  checkmate::assertList(params)

  w <- t(nmf_res$h)
  h <- t(nmf_res$w)

  k <- ncol(w)
  comp_names <- sprintf("comp_%02d", seq_len(k))

  rownames(w) <- gene_ids
  colnames(w) <- comp_names
  rownames(h) <- comp_names
  colnames(h) <- cell_ids

  res <- list(
    w = w,
    h = h,
    gene_ids = gene_ids,
    cell_ids = cell_ids,
    cell_indices = cell_indices,
    final_loss = nmf_res$final_loss,
    n_iter = nmf_res$n_iter,
    converged = nmf_res$converged,
    source_class = source_class,
    params = params
  )

  class(res) <- "NmfResult"
  res
}

#### primitives ----------------------------------------------------------------

#' @export
#'
#' @keywords internal
print.NmfResult <- function(x, ...) {
  cat("NmfResult (single-run HALS NMF)\n")
  cat(sprintf("  Source class:     %s\n", x$source_class))
  cat(sprintf("  No genes:         %d\n", nrow(x$w)))
  cat(sprintf("  No cells:         %d\n", ncol(x$h)))
  cat(sprintf("  No components:    %d\n", ncol(x$w)))
  cat(sprintf("  Final loss:       %.4g\n", x$final_loss))
  cat(sprintf("  Iterations:       %d\n", x$n_iter))
  cat(sprintf("  Converged:        %s\n", x$converged))
  cat(sprintf("  Preprocessing:    %s\n", x$params$preprocessing))
  invisible(x)
}

#' @export
#'
#' @keywords internal
dim.NmfResult <- function(x) {
  c(nrow(x$w), ncol(x$h), ncol(x$w))
}

#### transforms ----------------------------------------------------------------

#' @export
as.matrix.NmfResult <- function(x, which = c("w", "h"), ...) {
  which <- match.arg(which)
  x[[which]]
}

#### getters -------------------------------------------------------------------

#' Get the W (gene loadings) matrix
#'
#' @param x An object holding NMF results.
#'
#' @export
get_w <- function(x) {
  UseMethod("get_w")
}

#' @rdname get_w
#'
#' @export
get_w.NmfResult <- function(x) {
  checkmate::assertClass(x, "NmfResult")
  x$w
}

#' Get the H (cell activations) matrix
#'
#' @param x An object holding NMF results.
#'
#' @export
get_h <- function(x) {
  UseMethod("get_h")
}

#' @rdname get_h
#'
#' @export
get_h.NmfResult <- function(x) {
  checkmate::assertClass(x, "NmfResult")
  x$h
}

#' @method get_params NmfResult
#'
#' @export
S7::method(get_params, S7::new_S3_class("NmfResult")) <-
  function(object, to_json = FALSE, pretty_json = FALSE) {
    get_params.NmfResult(
      object = object,
      to_json = to_json,
      pretty_json = pretty_json
    )
  }

#' @rdname get_params
#'
#' @export
get_params.NmfResult <- function(
  object,
  to_json = FALSE,
  pretty_json = FALSE
) {
  checkmate::assertClass(object, "NmfResult")
  checkmate::qassert(to_json, "B1")
  checkmate::qassert(pretty_json, "B1")

  to_ret <- object[["params"]]
  if (to_json) {
    to_ret <- jsonlite::toJSON(to_ret)
  }
  if (to_json && pretty_json) {
    to_ret <- jsonlite::prettify(to_ret)
  }
  to_ret
}

#' @rdname get_data
#'
#' @export
get_data.NmfResult <- function(x, ...) {
  checkmate::assertClass(x, "NmfResult")

  obs_dt <- data.table::as.data.table(t(x$h))
  obs_dt[, cell_idx := (x$cell_indices + 1L)]

  obs_dt <- .add_is_obs_attr(obs_dt)
  obs_dt
}

### StabilisedNmfResult --------------------------------------------------------

#' Constructor for stabilised (multi-run) NMF results
#'
#' @description
#' Stores the column-bound W matrices across runs, per-run H matrices,
#' per-run losses and convergence flags, plus the index of the best run.
#'
#' @param nmf_res Named list. Output of [rs_nmf_multi_sc()] or
#' [rs_nmf_multi_mc()]. Must contain `w_all`, `h_per_run`, `losses`,
#' `converged`, `best_idx`.
#' @param gene_ids Character vector. Gene identifiers for the rows of W.
#' @param cell_ids Character vector. Cell identifiers for the columns of H.
#' @param cell_indices Integer vector. 0-indexed cell positions used in the run.
#' @param source_class String. One of `c("SingleCells", "MetaCells")`.
#' @param params List. The full set of parameters used for the run.
#'
#' @returns An object of class `StabilisedNmfResult`.
#'
#' @export
#'
#' @keywords internal
new_stabilised_nmf_result <- function(
  nmf_res,
  gene_ids,
  cell_ids,
  cell_indices,
  source_class,
  params
) {
  checkmate::assertList(nmf_res)
  checkmate::assertNames(
    names(nmf_res),
    must.include = c("w_all", "h_per_run", "losses", "converged", "best_idx")
  )
  checkmate::qassert(gene_ids, "S+")
  checkmate::qassert(cell_ids, "S+")
  checkmate::qassert(cell_indices, "I+")
  checkmate::assertChoice(source_class, c("SingleCells", "MetaCells"))
  checkmate::assertList(params)

  k <- params$k
  n_runs <- length(nmf_res$h_per_run)
  comp_names <- sprintf("comp_%02d", seq_len(k))

  # Each h_per_run entry is (k x n_genes) from Rust; transpose to (n_genes x k)
  # and column-bind into w_all-style matrix.
  w_all <- do.call(cbind, lapply(nmf_res$h_per_run, t))

  # Each w in the original w_all column block is the cell-loading matrix;
  # transpose into per-run H matrices of shape (k x n_cells).
  h_per_run <- lapply(seq_len(n_runs), function(i) {
    col_start <- (i - 1L) * k + 1L
    col_end <- i * k
    t(nmf_res$w_all[, col_start:col_end, drop = FALSE])
  })

  w_all_colnames <- as.vector(outer(
    comp_names,
    sprintf("run_%02d", seq_len(n_runs)),
    function(comp, run) paste(run, comp, sep = ".")
  ))
  rownames(w_all) <- gene_ids
  colnames(w_all) <- w_all_colnames

  h_per_run <- lapply(h_per_run, function(h) {
    rownames(h) <- comp_names
    colnames(h) <- cell_ids
    h
  })
  names(h_per_run) <- sprintf("run_%02d", seq_len(n_runs))

  res <- list(
    w_all = w_all,
    h_per_run = h_per_run,
    gene_ids = gene_ids,
    cell_ids = cell_ids,
    cell_indices = cell_indices,
    losses = nmf_res$losses,
    converged = nmf_res$converged,
    best_idx = nmf_res$best_idx,
    source_class = source_class,
    params = params
  )

  class(res) <- "StabilisedNmfResult"
  res
}

#### primitives ----------------------------------------------------------------

#' @export
#'
#' @keywords internal
print.StabilisedNmfResult <- function(x, ...) {
  n_runs <- length(x$h_per_run)
  cat("StabilisedNmfResult (multi-run HALS NMF)\n")
  cat(sprintf("  Source class:     %s\n", x$source_class))
  cat(sprintf("  No genes:         %d\n", nrow(x$w_all)))
  cat(sprintf("  No cells:         %d\n", ncol(x$h_per_run[[1]])))
  cat(sprintf("  No components:    %d\n", x$params$k))
  cat(sprintf("  No runs:          %d\n", n_runs))
  cat(sprintf("  No converged:     %d / %d\n", sum(x$converged), n_runs))
  cat(sprintf(
    "  Loss range:       [%.4g, %.4g]\n",
    min(x$losses),
    max(x$losses)
  ))
  cat(sprintf(
    "  Best run:         %d (loss = %.4g)\n",
    x$best_idx,
    x$losses[x$best_idx]
  ))
  cat(sprintf("  Preprocessing:    %s\n", x$params$preprocessing))
  invisible(x)
}

#### getters -------------------------------------------------------------------

#' @rdname get_w
#'
#' @export
get_w.StabilisedNmfResult <- function(x) {
  checkmate::assertClass(x, "StabilisedNmfResult")
  x$w_all
}

#' @rdname get_h
#'
#' @export
get_h.StabilisedNmfResult <- function(x) {
  checkmate::assertClass(x, "StabilisedNmfResult")
  x$h_per_run
}

#' Get the best run from a stabilised NMF result
#'
#' @description
#' Extracts the run with the lowest reconstruction loss and returns it as a
#' single-run `NmfResult` for downstream methods.
#'
#' @param x `StabilisedNmfResult` object.
#'
#' @returns An `NmfResult` containing the W/H of the best run.
#'
#' @export
get_best_run <- function(x) {
  UseMethod("get_best_run")
}

#' @rdname get_best_run
#'
#' @export
get_best_run.StabilisedNmfResult <- function(x) {
  checkmate::assertClass(x, "StabilisedNmfResult")

  k <- x$params$k
  best <- x$best_idx
  col_start <- (best - 1L) * k + 1L
  col_end <- best * k

  w_best <- x$w_all[, col_start:col_end, drop = FALSE]
  h_best <- x$h_per_run[[best]]

  comp_names <- sprintf("comp_%02d", seq_len(k))
  colnames(w_best) <- comp_names
  rownames(h_best) <- comp_names

  res <- list(
    w = w_best,
    h = h_best,
    gene_ids = x$gene_ids,
    cell_ids = x$cell_ids,
    cell_indices = x$cell_indices,
    final_loss = x$losses[best],
    n_iter = NA_integer_,
    converged = x$converged[best],
    source_class = x$source_class,
    params = x$params
  )

  class(res) <- "NmfResult"
  res
}

#' @method get_params StabilisedNmfResult
#'
#' @export
S7::method(get_params, S7::new_S3_class("StabilisedNmfResult")) <-
  function(object, to_json = FALSE, pretty_json = FALSE) {
    get_params.StabilisedNmfResult(
      object = object,
      to_json = to_json,
      pretty_json = pretty_json
    )
  }

#' @rdname get_params
#'
#' @export
get_params.StabilisedNmfResult <- function(
  object,
  to_json = FALSE,
  pretty_json = FALSE
) {
  checkmate::assertClass(object, "StabilisedNmfResult")
  checkmate::qassert(to_json, "B1")
  checkmate::qassert(pretty_json, "B1")

  to_ret <- object[["params"]]
  if (to_json) {
    to_ret <- jsonlite::toJSON(to_ret)
  }
  if (to_json && pretty_json) {
    to_ret <- jsonlite::prettify(to_ret)
  }
  to_ret
}
