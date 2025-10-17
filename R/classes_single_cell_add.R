# additional single cell classes and methods -----------------------------------

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
