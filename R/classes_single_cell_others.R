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
plot.scrublet_res <- function(
  x,
  break_number = 31L,
  ...
) {
  # checks
  checkmate::assertClass(x, "scrublet_res")
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
#' @param scrublet_res `scrublet_res` result class.
#' @param threshold Numeric. The new threshold to use.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @returns `scrublet_res` class with updated doublet calls based on the new
#' threshold.
#'
#' @export
call_doublets_manual <- function(scrublet_res, threshold, .verbose = TRUE) {
  # checks
  checkmate::assertClass(scrublet_res, "scrublet_res")
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
get_obs_data.scrublet_res <- function(x, ...) {
  # checks
  checkmate::assertClass(x, "scrublet_res")

  # function body
  obs_dt <- data.table::data.table(
    doublet = x$predicted_doublets,
    doublet_score = x$doublet_scores_obs
  )
  obs_dt[, cell_idx := (attr(x, "cell_indices") + 1)] # was zero indexed

  return(obs_dt)
}

## boost -----------------------------------------------------------------------

#' @rdname get_obs_data
#'
#' @export
get_obs_data.boost_res <- function(x, ...) {
  # checks
  checkmate::assertClass(x, "boost_res")

  # function body
  obs_dt <- data.table::as.data.table(
    unclass(x)
  )
  obs_dt[, cell_idx := (attr(x, "cell_indices") + 1)] # was zero indexed

  return(obs_dt)
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
#' @return Generates the `sc_knn` class.
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

  class(sc_knn) <- "sc_knn"

  return(sc_knn)
}

### methods --------------------------------------------------------------------

#### getters -------------------------------------------------------------------

#' @rdname get_knn_mat
#'
#' @export
get_knn_mat.sc_knn <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_knn")

  return(x[["indices"]])
}

#' @rdname get_knn_dist
#'
#' @export
get_knn_dist.sc_knn <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_knn")

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

  class(sc_hotspot) <- "sc_hotspot"

  return(sc_hotspot)
}

### functions ------------------------------------------------------------------

#### getters -------------------------------------------------------------------

#' @method get_params sc_hotspot
#'
#' @export
S7::method(get_params, S7::new_S3_class("sc_hotspot")) <-
  function(object, to_json = FALSE, pretty_json = FALSE) {
    get_params.sc_hotspot(
      object = object,
      to_json = to_json,
      pretty_json = pretty_json
    )
  }


#' @rdname get_params
#'
#' @export
get_params.sc_hotspot <- function(
  object,
  to_json = FALSE,
  pretty_json = FALSE
) {
  # Checks
  checkmate::assertClass(
    object,
    "sc_hotspot"
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
get_hotspot_membership.sc_hotspot <- function(
  x
) {
  # checks
  checkmate::assertClass(x, "sc_hotspot")

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
set_hotspot_membership <- function(x, fdr_threshold = 0.05, min_size = 10L) {
  UseMethod("set_hotspot_membership")
}

#' @rdname set_hotspot_membership
#'
#' @export
set_hotspot_membership.sc_hotspot <- function(
  x,
  fdr_threshold = 0.05,
  min_size = 10L
) {
  # checks
  checkmate::assertClass(x, "sc_hotspot")
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
#' @returns An `sc_miloR` class that contains the provided data and has
#' subsequent methods to calculate differential abundance statistics.
#'
#' @references Dann, et al., Nat Biotechnol, 2022
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

  class(sc_milor) <- "sc_miloR"

  return(sc_milor)
}

### getters --------------------------------------------------------------------

#' @method get_params sc_miloR
#'
#' @export
S7::method(get_params, S7::new_S3_class("sc_miloR")) <-
  function(object, to_json = FALSE, pretty_json = FALSE) {
    get_params.sc_miloR(
      object = object,
      to_json = to_json,
      pretty_json = pretty_json
    )
  }

#' @rdname get_params
#'
#' @export
get_params.sc_miloR <- function(
  object,
  to_json = FALSE,
  pretty_json = FALSE
) {
  # Checks
  checkmate::assertClass(
    object,
    "sc_miloR"
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
get_differential_abundance_res.sc_miloR <- function(
  x
) {
  # checks
  checkmate::assertClass(x, "sc_miloR")

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
get_model_fit.sc_miloR <- function(
  x
) {
  # checks
  checkmate::assertClass(x, "sc_miloR")

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
get_index_cells.sc_miloR <- function(x) {
  # checks
  checkmate::assertClass(x, "sc_miloR")

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
#' @param x `sc_miloR` object for which to run the differential abundance
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
#' @return The `sc_miloR` object with added model and results from the
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
test_nhoods.sc_miloR <- function(
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
  checkmate::assertClass(x, "sc_miloR")
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
#' slot within the `sc_miloR` object. For each neighbourhood, it calculates
#' the proportion of the majority cell type and identifies which cell type
#' is most abundant. This is useful for annotating differential abundance
#' results with the cellular composition of each neighbourhood.
#'
#' @param x `sc_miloR` object on which to tag on additional neighbourhood
#' information.
#' @param cell_info Character vector. Represents the cell type annotations
#' you wish to add to the different neighbourhoods. Must be the same length
#' as the number of cells (rows) in the nhoods matrix.
#'
#' @return Modified `sc_miloR` object with updated `nhoods_info` containing
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
add_nhoods_info.sc_miloR <- function(x, cell_info) {
  # checks
  checkmate::assertClass(x, "sc_miloR")
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
