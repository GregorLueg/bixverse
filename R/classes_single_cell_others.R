# additional single cell classes and methods -----------------------------------

# this file contains additional, smaller s3 classes relevant for single cell
# analyses. general generics are found on the top of the file, otherwise you
# have per class the class generation and specific methods, getters, setters.

## general generics ------------------------------------------------------------

#' Get the ready obs data from various sub method
#'
#' @description
#' Helper method that creates data.tables with cell indices which were used
#' in the given analysis + the values that are to be added to the obs table
#' in the DuckDB.
#'
#' @param x An object to set gene mapping for
#' @param ... Other parameters
#'
#' @returns Returns a data.table with a cell_idx column for the cells included
#' in the analysis and additional columns to be added to the obs table.
#'
#' @export
get_obs_data <- function(x, ...) {
  UseMethod("get_obs_data")
}

# methods ----------------------------------------------------------------------

## gene proportion analysis ----------------------------------------------------

#' @rdname get_obs_data
#'
#' @export
get_obs_data.sc_proportion_res <- function(x, ...) {
  # checks
  checkmate::assertClass(x, "sc_proportion_res")

  # function body
  obs_dt <- data.table::as.data.table(unclass(x))
  obs_dt[, cell_idx := (attr(x, "cell_indices") + 1)] # was zero indexed

  return(obs_dt)
}

# additional S3 classes --------------------------------------------------------

# contains generics/methods for additional S3 classes related to single cell
# analysis

## scrublet --------------------------------------------------------------------

### histogram plotting ---------------------------------------------------------

#' @export
plot.ScrubletRes <- function(
  x,
  break_number = 31L,
  ...
) {
  # checks
  checkmate::assertClass(x, "ScrubletRes")
  checkmate::qassert(break_number, "I1")

  # plotting
  obs_plot <- ggplot2::ggplot(
    data.frame(score = x$doublet_scores_obs),
    ggplot2::aes(x = score)
  ) +
    ggplot2::geom_histogram(
      breaks = seq(0, 1, length.out = break_number),
      fill = "grey",
      colour = "darkgrey",
      ggplot2::aes(y = ggplot2::after_stat(density))
    ) +
    ggplot2::geom_vline(
      xintercept = x$threshold,
      linewidth = 0.5,
      colour = "red",
      linetype = "dashed"
    ) +
    ggplot2::labs(
      title = "Observed transcriptomes",
      x = "Doublet score",
      y = "Probability density"
    ) +
    ggplot2::theme_bw()

  sim_plot <- ggplot2::ggplot(
    data.frame(score = x$doublet_scores_sim),
    ggplot2::aes(x = score)
  ) +
    ggplot2::geom_histogram(
      breaks = seq(0, 1, length.out = break_number),
      fill = "grey",
      colour = "darkgrey",
      ggplot2::aes(y = ggplot2::after_stat(density))
    ) +
    ggplot2::geom_vline(
      xintercept = x$threshold,
      linewidth = 0.5,
      colour = "red",
      linetype = "dashed"
    ) +
    ggplot2::labs(
      title = "Simulated doublets",
      x = "Doublet score",
      y = "Probability density"
    ) +
    ggplot2::theme_bw()

  patchwork::wrap_plots(obs_plot, sim_plot, ncol = 2)
}

### readjusting thresholds -----------------------------------------------------

#' Helper function to manually readjust Scrublet thresholds
#'
#' @description
#' Can update the Scrublet thresholding after manual introspection of the
#' histogram.
#'
#' @param scrublet_res `ScrubletRes` result class.
#' @param threshold Numeric. The new threshold to use.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @returns `ScrubletRes` class with updated doublet calls based on the new
#' threshold.
#'
#' @export
call_doublets_manual <- function(scrublet_res, threshold, .verbose = TRUE) {
  # checks
  checkmate::assertClass(scrublet_res, "ScrubletRes")
  checkmate::qassert(threshold, "N1[0,1]")
  checkmate::qassert(.verbose, "B1")

  # function
  predicted_doublets <- scrublet_res$doublet_scores_obs > threshold

  z_scores <- (scrublet_res$doublet_scores_obs - threshold) /
    scrublet_res$doublet_errors_obs

  detected_doublet_rate <- sum(predicted_doublets) /
    length(predicted_doublets)

  detectable_doublet_fraction <- sum(
    scrublet_res$doublet_scores_sim > threshold
  ) /
    length(scrublet_res$doublet_scores_sim)

  overall_doublet_rate <- if (detectable_doublet_fraction > 0.01) {
    detected_doublet_rate / detectable_doublet_fraction
  } else {
    0.0
  }

  if (.verbose) {
    cat(sprintf(
      "Detected doublet rate = %.1f%%\n",
      100 * detected_doublet_rate
    ))
    cat(sprintf(
      "Estimated detectable doublet fraction = %.1f%%\n",
      100 * detectable_doublet_fraction
    ))
    cat(sprintf("Overall doublet rate:\n"))
    cat(sprintf("  Estimated = %.1f%%\n", 100 * overall_doublet_rate))
  }

  scrublet_res$predicted_doublets <- predicted_doublets
  scrublet_res$z_scores <- z_scores
  scrublet_res$threshold <- threshold
  scrublet_res$detected_doublet_rate <- detected_doublet_rate
  scrublet_res$detectable_doublet_fraction <- detectable_doublet_fraction
  scrublet_res$overall_doublet_rate <- overall_doublet_rate

  scrublet_res
}

### obs data -------------------------------------------------------------------

#' @rdname get_obs_data
#'
#' @export
get_obs_data.ScrubletRes <- function(x, ...) {
  # checks
  checkmate::assertClass(x, "ScrubletRes")

  # function body
  obs_dt <- data.table::data.table(
    doublet = x$predicted_doublets,
    doublet_score = x$doublet_scores_obs
  )
  obs_dt[, cell_idx := (attr(x, "cell_indices") + 1)] # was zero indexed

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

  cat(sprintf(
    "ScrubletRes: %d cells, %d doublets (%.1f%%)\n",
    n_cells,
    n_doublets,
    100 * n_doublets / n_cells
  ))
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
  cat(sprintf("  Simulated doublets:     %d\n", length(x$doublet_scores_sim)))

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

#' @rdname get_obs_data
#'
#' @export
get_obs_data.BoostRes <- function(x, ...) {
  # checks
  checkmate::assertClass(x, "BoostRes")

  # function body
  obs_dt <- data.table::as.data.table(
    unclass(x)
  )
  obs_dt[, cell_idx := (attr(x, "cell_indices") + 1)] # was zero indexed

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

  invisible(x)
}

## scdblfinder -----------------------------------------------------------------

### obs data -------------------------------------------------------------------

#' @rdname get_obs_data
#'
#' @export
get_obs_data.ScDblFinderRes <- function(x, ...) {
  checkmate::assertClass(x, "ScDblFinderRes")
  obs_dt <- data.table::as.data.table(unclass(x))
  obs_dt[, cell_idx := (attr(x, "cell_indices") + 1)]
  return(obs_dt)
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

  invisible(x)
}

## kNN with distances ----------------------------------------------------------

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

  sc_knn <- list(
    indices = knn_data$indices,
    dist = knn_data$dist,
    dist_metric = knn_data$dist_metric,
    used_cells = used_cells
  )

  class(sc_knn) <- "SingleCellNearestNeighbour"

  return(sc_knn)
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
    module_memership = NULL
  )

  class(sc_hotspot) <- "Hotspot"

  return(sc_hotspot)
}

### functions ------------------------------------------------------------------

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

  x[["module_memership"]]
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

  x[["module_memership"]] <- gene_membership

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

#' Identify the TF to gene regulation
#'
#' @description
#' The function will generate the TF to gene associations that can be further
#' subfiltered subsequently.
#'
#' @details
#' You have three options to extract the information
#' based on the importance scores generated by the tree-based regression model.
#' \itemize{
#'  \item `k_tfs` - This defines how many of the Top TFs per given gene you wish
#'  to include in a given analysis. If you assume that there are 5 TFs
#'  controlling your gene, you would set this to `5L`. For each gene, the Top 5
#'  TFs are identified and subsequently the TF -> c(gene_1, gene_2, gene_3)
#'  associations are generated.
#'  \item k_genes - This defines how many target genes per given TF you wish
#'  to include. If you assume that a TF has a 100 potential target genes, you
#'  would set this to `100L`. For each TF, the Top 100 potential downstream
#'  targets are being included in this case.
#'  \item min_importance - If you want to set a threshold for the minimum
#'  importance score to filter out noisy genes.
#' }
#' You can provide all three parameters at once, in this case you will get a
#' union of the TF -> gene, gene <- TF approach, filtered by min_importance.
#' This is the first step and you can subsequently filter by correlation of
#' TF to target gene and motif enrichment for a given TF.
#'
#' @param x `ScenicGrn` object for which to generate the TF to gene
#' associations.
#' @param k_tfs Optional integer. How many TFs per given gene you want to
#' include.
#' @param k_genes Optional integer. How many genes you want to include
#' downstream of each TF. Warning. `k_tfs = NULL` and `k_genes = NULL` does
#' not work.
#' @param min_importance Optional float. If you want a minimum importance
#' score for including a TF to gene association.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns Adds a data.table with the first tf to gene results to the class.
#'
#' @export
identify_tf_to_genes <- function(
  x,
  k_tfs,
  k_genes,
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
  k_tfs,
  k_genes,
  min_importance = NULL,
  .verbose = TRUE
) {
  # checks
  checkmate::assertClass(x, "ScenicGrn")
  checkmate::qassert(k_tfs, c("I1", "0"))
  checkmate::qassert(k_genes, c("I1", "0"))
  checkmate::qassert(min_importance, c("N1[0, 1]", "0"))
  checkmate::qassert(.verbose, "B1")
  if (is.null(k_tfs) & is.null(k_genes)) {
    stop("k_tfs and k_genes cannot be both NULL.")
  }

  if (.verbose) {
    message(
      "Extracting TF to gene associations based on the importance values."
    )
  }

  gene_tf_imp <- x$importance_matrix
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
  tf_gene_dt <- unique(tf_gene_dt, by = c("tf", "gene"))

  x$tf_to_gene_results <- tf_gene_dt
  return(x)
}

### tf to gene correlation -----------------------------------------------------

#' Generate TF to gene correlations
#'
#' @description
#' This function will calculate the correlations between the identified TF to
#' gene pairs. You need to have run [identify_tf_to_genes()]!
#'
#' @param x `ScenicGrn` object for which to generate the TF to gene
#' associations.
#' @param object `SingleCells` object that was used to generate the original
#' GRNs.
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
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
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

  indices_1 <- get_sc_map(object)$gene_mapping[tf_to_gene$tf] - 1L
  indices_2 <- get_sc_map(object)$gene_mapping[tf_to_gene$gene] - 1L

  pairwise_cors <- rs_pairwise_gene_cors(
    f_path = get_rust_count_gene_f_path(object),
    gene_indices_1 = indices_1,
    gene_indices_2 = indices_2,
    cells_to_keep = get_cells_to_keep(object),
    spearman = FALSE
  )

  tf_to_gene[, pairwise_cor := pairwise_cors]

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

  tf_gene_dt <- get_tf_to_gene(x)
  tf_gene_lists <- split(tf_gene_dt$gene, tf_gene_dt$tf)

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
  explode_leading_edge <- function(dt, tf_col) {
    dt <- dt[!is.na(get(tf_col))]
    if (nrow(dt) == 0L) {
      return(data.table::data.table(tf = character(), le_gene = character()))
    }
    dt[,
      .(single_tf = unlist(strsplit(get(tf_col), ";"))),
      by = .(gs_name, leading_edge_genes)
    ][
      single_tf == gs_name,
      .(le_gene = unlist(strsplit(leading_edge_genes, ";"))),
      by = .(tf = gs_name)
    ]
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
