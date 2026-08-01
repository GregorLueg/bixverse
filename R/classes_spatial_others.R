# ------------------------------------------------------------------------------
# small S3 result classes for the spatial methods. One constructor per method,
# plus the primitives that make the result readable at the prompt.
#
# Moran's I and SPARK-X return data.tables with an extra class in front. The
# SpCache setters assert a data.table and users want one anyway; the class only
# buys a print method that does not dump 19k rows and a diagnostic plot.
# ------------------------------------------------------------------------------

# spatial result classes -------------------------------------------------------

## moran's i -------------------------------------------------------------------

#' Helper function to generate Moran's I results
#'
#' @description
#' Wraps the [bixverse::rs_sp_morans_i()] output into a `data.table` with the
#' gene identifiers joined back on. The null expectation and variance are
#' scalars across genes, so they ride along as attributes rather than as two
#' constant columns.
#'
#' @param morans_res Named list. Output of [bixverse::rs_sp_morans_i()].
#' @param gene_ids Character vector. Gene identifiers for the tested genes, in
#' the order of `morans_res$gene_idx`.
#' @param exp_id String. The experiment the result belongs to.
#'
#' @return An `SpMoransRes` object, which is a `data.table` with columns
#' `gene_id`, `gene_idx`, `morans_i`, `z`, `pval` and `fdr`.
#'
#' @export
#'
#' @keywords internal
new_sp_morans_res <- function(morans_res, gene_ids, exp_id) {
  # checks
  checkmate::assertList(morans_res)
  checkmate::assertNames(
    names(morans_res),
    must.include = c("gene_idx", "i", "e_i", "var_i", "z", "pval", "fdr")
  )
  checkmate::qassert(gene_ids, "S+")
  checkmate::qassert(exp_id, "S1")
  checkmate::assertTRUE(length(gene_ids) == length(morans_res$gene_idx))

  res <- data.table::data.table(
    gene_id = gene_ids,
    gene_idx = as.integer(morans_res$gene_idx),
    morans_i = as.numeric(morans_res$i),
    z = as.numeric(morans_res$z),
    pval = as.numeric(morans_res$pval),
    fdr = as.numeric(morans_res$fdr)
  )

  data.table::setattr(res, "e_i", as.numeric(morans_res$e_i))
  data.table::setattr(res, "var_i", as.numeric(morans_res$var_i))
  data.table::setattr(res, "exp_id", exp_id)
  data.table::setattr(res, "class", c("SpMoransRes", class(res)))

  return(res)
}

### primitives -----------------------------------------------------------------

#' @method print SpMoransRes
#'
#' @export
#'
#' @keywords internal
print.SpMoransRes <- function(x, ...) {
  cat("Moran's I spatially variable genes\n")
  cat(sprintf("  Experiment: %s\n", attr(x, "exp_id")))
  cat(sprintf("  Genes tested: %i\n", nrow(x)))
  cat(sprintf("  Constant genes (NaN): %i\n", sum(is.nan(x$morans_i))))
  cat(sprintf(
    "  E[I]: %.3g, Var[I]: %.3g\n",
    attr(x, "e_i"),
    attr(x, "var_i")
  ))
  cat(sprintf(
    "  Significant at 5%% FDR: %i\n",
    sum(x$fdr < 0.05, na.rm = TRUE)
  ))

  top <- data.table::data.table(
    gene_id = x$gene_id,
    morans_i = x$morans_i,
    fdr = x$fdr
  )
  data.table::setorderv(top, "fdr")
  cat("  Top genes:\n")
  print(utils::head(as.data.frame(top), 5L))

  invisible(x)
}

#' Plot the Moran's I results
#'
#' @description
#' Moran's I against the significance of the Cliff-Ord test, with the null
#' expectation marked. Genes to the right of the line and above the FDR cut are
#' the spatially structured ones.
#'
#' @param x An `SpMoransRes` object.
#' @param fdr_cutoff Float. FDR threshold used to colour the points.
#' @param ... Further arguments (currently unused).
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @method plot SpMoransRes
#'
#' @export
#'
#' @keywords internal
plot.SpMoransRes <- function(x, fdr_cutoff = 0.05, ...) {
  checkmate::assertClass(x, "SpMoransRes")
  checkmate::qassert(fdr_cutoff, "N1(0,1]")

  plot_dt <- data.table::data.table(
    morans_i = x$morans_i,
    neg_log_fdr = -log10(pmax(x$fdr, .Machine$double.xmin)),
    significant = !is.na(x$fdr) & x$fdr < fdr_cutoff
  )
  plot_dt <- plot_dt[is.finite(morans_i)]

  ggplot2::ggplot(
    plot_dt,
    ggplot2::aes(x = morans_i, y = neg_log_fdr, colour = significant)
  ) +
    ggplot2::geom_point(size = 0.6, alpha = 0.6) +
    ggplot2::geom_vline(
      xintercept = attr(x, "e_i"),
      linetype = "dashed",
      colour = "grey40"
    ) +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = "grey70", "TRUE" = "#B2182B"),
      name = sprintf("FDR < %.2g", fdr_cutoff)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Moran's I",
      y = expression(-log[10] ~ "FDR"),
      title = sprintf("Moran's I - %s", attr(x, "exp_id"))
    )
}

## sparkx ----------------------------------------------------------------------

#' Helper function to generate SPARK-X results
#'
#' @description
#' Wraps the [bixverse::rs_sp_sparkx()] output into a `data.table`. The
#' per-kernel statistic and p-value matrices ride along as attributes; they are
#' genes x kernels and only interesting when a single kernel is under
#' suspicion.
#'
#' @param sparkx_res Named list. Output of [bixverse::rs_sp_sparkx()].
#' @param gene_ids Character vector. Gene identifiers for the tested genes, in
#' the order of `sparkx_res$gene_idx`.
#' @param exp_id String. The experiment the result belongs to.
#'
#' @return An `SpSparkxRes` object, which is a `data.table` with columns
#' `gene_id`, `gene_idx`, `pval` and `fdr`.
#'
#' @export
#'
#' @keywords internal
new_sp_sparkx_res <- function(sparkx_res, gene_ids, exp_id) {
  # checks
  checkmate::assertList(sparkx_res)
  checkmate::assertNames(
    names(sparkx_res),
    must.include = c(
      "gene_idx",
      "kernels",
      "stat_per_kernel",
      "pval_per_kernel",
      "pval_combined",
      "fdr"
    )
  )
  checkmate::qassert(gene_ids, "S+")
  checkmate::qassert(exp_id, "S1")
  checkmate::assertTRUE(length(gene_ids) == length(sparkx_res$gene_idx))

  res <- data.table::data.table(
    gene_id = gene_ids,
    gene_idx = as.integer(sparkx_res$gene_idx),
    pval = as.numeric(sparkx_res$pval_combined),
    fdr = as.numeric(sparkx_res$fdr)
  )

  stat_per_kernel <- sparkx_res$stat_per_kernel
  pval_per_kernel <- sparkx_res$pval_per_kernel
  dimnames(stat_per_kernel) <- list(gene_ids, sparkx_res$kernels)
  dimnames(pval_per_kernel) <- list(gene_ids, sparkx_res$kernels)

  data.table::setattr(res, "kernels", sparkx_res$kernels)
  data.table::setattr(res, "stat_per_kernel", stat_per_kernel)
  data.table::setattr(res, "pval_per_kernel", pval_per_kernel)
  data.table::setattr(res, "exp_id", exp_id)
  data.table::setattr(res, "class", c("SpSparkxRes", class(res)))

  return(res)
}

### primitives -----------------------------------------------------------------

#' @method print SpSparkxRes
#'
#' @export
#'
#' @keywords internal
print.SpSparkxRes <- function(x, ...) {
  cat("SPARK-X spatially variable genes\n")
  cat(sprintf("  Experiment: %s\n", attr(x, "exp_id")))
  cat(sprintf("  Genes tested: %i\n", nrow(x)))
  cat(sprintf(
    "  Kernels: %i (%s)\n",
    length(attr(x, "kernels")),
    toString(utils::head(attr(x, "kernels"), 3L))
  ))
  cat(sprintf(
    "  Significant at 5%% FDR: %i\n",
    sum(x$fdr < 0.05, na.rm = TRUE)
  ))

  top <- data.table::data.table(
    gene_id = x$gene_id,
    pval = x$pval,
    fdr = x$fdr
  )
  data.table::setorderv(top, "fdr")
  cat("  Top genes:\n")
  print(utils::head(as.data.frame(top), 5L))

  invisible(x)
}

#' Plot the SPARK-X results
#'
#' @description
#' Histogram of the Cauchy-combined p-values. A flat bar at every bin with a
#' spike at zero is what a working test on real tissue looks like; a hump in the
#' middle means the null is off.
#'
#' @param x An `SpSparkxRes` object.
#' @param bins Integer. Number of histogram bins.
#' @param ... Further arguments (currently unused).
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @method plot SpSparkxRes
#'
#' @export
#'
#' @keywords internal
plot.SpSparkxRes <- function(x, bins = 50L, ...) {
  checkmate::assertClass(x, "SpSparkxRes")
  checkmate::qassert(bins, "I1[2,)")

  plot_dt <- data.table::data.table(pval = x$pval)
  plot_dt <- plot_dt[is.finite(pval)]

  ggplot2::ggplot(plot_dt, ggplot2::aes(x = pval)) +
    ggplot2::geom_histogram(bins = bins, fill = "#2166AC", colour = NA) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Cauchy-combined p-value",
      y = "Genes",
      title = sprintf("SPARK-X - %s", attr(x, "exp_id"))
    )
}

## neighbourhood enrichment ----------------------------------------------------

#' Helper function to generate neighbourhood enrichment results
#'
#' @description
#' Wraps the [bixverse::rs_sp_nhood_enrichment()] output. All four matrices are
#' K x K over the label levels, with matching row and column order.
#'
#' @param nhood_res Named list. Output of [bixverse::rs_sp_nhood_enrichment()].
#' @param label_col String. The obs column the labels came from.
#' @param exp_id String. The experiment the result belongs to.
#' @param n_perm Integer. Number of permutations that were run.
#'
#' @return An `SpNhoodRes` object, a list with `label_levels`, `observed`,
#' `perm_mean`, `perm_sd`, `z_scores` and `params`.
#'
#' @export
#'
#' @keywords internal
new_sp_nhood_res <- function(nhood_res, label_col, exp_id, n_perm) {
  # checks
  checkmate::assertList(nhood_res)
  checkmate::assertNames(
    names(nhood_res),
    must.include = c(
      "label_levels",
      "observed",
      "perm_mean",
      "perm_sd",
      "z_scores"
    )
  )
  checkmate::qassert(label_col, "S1")
  checkmate::qassert(exp_id, "S1")
  checkmate::qassert(n_perm, "I1[1,)")

  levels <- as.character(nhood_res$label_levels)
  name_it <- \(m) {
    dimnames(m) <- list(levels, levels)
    m
  }

  res <- list(
    label_levels = levels,
    observed = name_it(nhood_res$observed),
    perm_mean = name_it(nhood_res$perm_mean),
    perm_sd = name_it(nhood_res$perm_sd),
    z_scores = name_it(nhood_res$z_scores),
    params = list(
      label_col = label_col,
      exp_id = exp_id,
      n_perm = as.integer(n_perm)
    )
  )

  class(res) <- "SpNhoodRes"

  return(res)
}

### primitives -----------------------------------------------------------------

#' @method print SpNhoodRes
#'
#' @export
#'
#' @keywords internal
print.SpNhoodRes <- function(x, ...) {
  z <- x$z_scores
  diag_z <- diag(z)
  off <- z[upper.tri(z)]

  cat("Neighbourhood enrichment\n")
  cat(sprintf("  Experiment: %s\n", x$params$exp_id))
  cat(sprintf("  Label column: %s\n", x$params$label_col))
  cat(sprintf(
    "  Levels: %i (%s)\n",
    length(x$label_levels),
    toString(
      utils::head(x$label_levels, 5L)
    )
  ))
  cat(sprintf("  Permutations: %i\n", x$params$n_perm))
  cat(sprintf(
    "  Z range: %.2f to %.2f\n",
    min(z, na.rm = TRUE),
    max(z, na.rm = TRUE)
  ))
  cat(sprintf(
    "  Self-neighbouring (diagonal) median Z: %.2f\n",
    stats::median(diag_z, na.rm = TRUE)
  ))
  cat(sprintf(
    "  Cross-label (off-diagonal) median Z: %.2f\n",
    stats::median(off, na.rm = TRUE)
  ))

  invisible(x)
}

#' Plot the neighbourhood enrichment Z-scores
#'
#' @description
#' Heatmap of the K x K Z-score matrix. Red is more co-occurrence than the
#' permuted null gives, blue is less.
#'
#' @param x An `SpNhoodRes` object.
#' @param ... Further arguments (currently unused).
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @method plot SpNhoodRes
#'
#' @export
#'
#' @keywords internal
plot.SpNhoodRes <- function(x, ...) {
  checkmate::assertClass(x, "SpNhoodRes")

  levels <- x$label_levels
  plot_dt <- data.table::as.data.table(as.table(x$z_scores))
  data.table::setnames(plot_dt, c("label_1", "label_2", "z"))
  plot_dt[, label_1 := factor(label_1, levels = levels)]
  plot_dt[, label_2 := factor(label_2, levels = rev(levels))]

  cap <- max(abs(plot_dt$z), na.rm = TRUE)

  ggplot2::ggplot(plot_dt, ggplot2::aes(x = label_1, y = label_2)) +
    ggplot2::geom_tile(ggplot2::aes(fill = z)) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-cap, cap),
      name = "Z-score"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = sprintf(
        "Neighbourhood enrichment - %s (%s)",
        x$params$exp_id,
        x$params$label_col
      )
    )
}

## image features --------------------------------------------------------------

#' Helper function to generate histology image feature results
#'
#' @description
#' Wraps the [bixverse::rs_sp_image_features()] output. Spots whose tile falls
#' off the image are dropped inside Rust, so `cell_idx` can be shorter than the
#' sample.
#'
#' @param image_res Named list. Output of [bixverse::rs_sp_image_features()].
#' @param spot_cell_idx Integer vector. 1-indexed `cell_idx` of every spot that
#' went in, in coordinate row order.
#' @param exp_id String. The experiment the result belongs to.
#' @param resolution String. The image resolution the tiles were cut from.
#'
#' @return An `SpImageFeatures` object, a list with `values` (spots x
#' features), `feature_names`, `cell_idx` and `params`.
#'
#' @export
#'
#' @keywords internal
new_sp_image_features_res <- function(
  image_res,
  spot_cell_idx,
  exp_id,
  resolution
) {
  # checks
  checkmate::assertList(image_res)
  checkmate::assertNames(
    names(image_res),
    must.include = c("spot_idx", "feature_names", "values")
  )
  checkmate::qassert(spot_cell_idx, "I+")
  checkmate::qassert(exp_id, "S1")
  checkmate::qassert(resolution, "S1")

  # `spot_idx` is 0-based positional into the coordinate matrix, so it maps
  # straight onto the cell_idx vector that produced those coordinates.
  kept <- as.integer(image_res$spot_idx) + 1L
  if (any(kept < 1L) || any(kept > length(spot_cell_idx))) {
    stop("Image features returned spot indices outside the coordinate matrix.")
  }
  cell_idx <- as.integer(spot_cell_idx[kept])

  values <- image_res$values
  dimnames(values) <- list(as.character(cell_idx), image_res$feature_names)

  res <- list(
    values = values,
    feature_names = as.character(image_res$feature_names),
    cell_idx = cell_idx,
    params = list(
      exp_id = exp_id,
      resolution = resolution,
      n_dropped = length(spot_cell_idx) - length(cell_idx)
    )
  )

  class(res) <- "SpImageFeatures"

  return(res)
}

### primitives -----------------------------------------------------------------

#' @method print SpImageFeatures
#'
#' @export
#'
#' @keywords internal
print.SpImageFeatures <- function(x, ...) {
  cat("Histology image features\n")
  cat(sprintf("  Experiment: %s\n", x$params$exp_id))
  cat(sprintf("  Resolution: %s\n", x$params$resolution))
  cat(sprintf("  Spots: %i\n", length(x$cell_idx)))
  cat(sprintf("  Spots dropped (tile off image): %i\n", x$params$n_dropped))
  cat(sprintf(
    "  Features: %i (%s, ...)\n",
    length(x$feature_names),
    toString(utils::head(x$feature_names, 3L))
  ))

  invisible(x)
}

### getters --------------------------------------------------------------------

#' @rdname get_data
#'
#' @method get_data SpImageFeatures
#'
#' @export
get_data.SpImageFeatures <- function(x, columns = NULL, ...) {
  checkmate::assertClass(x, "SpImageFeatures")
  checkmate::qassert(columns, c("0", "S+"))

  cols <- if (is.null(columns)) x$feature_names else columns
  checkmate::assertSubset(cols, x$feature_names)

  res <- data.table::as.data.table(x$values[, cols, drop = FALSE])
  res[, cell_idx := x$cell_idx]
  data.table::setcolorder(res, "cell_idx")

  return(.add_is_obs_attr(res))
}
