## plots -----------------------------------------------------------------------

### bulk dge -------------------------------------------------------------------

#' @title Helper plot function of distribution of genes by samples
#'
#' @param samples data.table with sample information with nb_detected_genes and
#' a column specifying the cohort.
#' @param group_col String specifying the column with cohort information
#'
#' @return ggplot object, i.e. boxplot with number of genes by cohort
plot_preprocessing_genes <- function(samples, group_col) {
  # checks
  checkmate::assertDataFrame(samples)
  checkmate::assertNames(
    colnames(samples),
    must.include = c(group_col, "nb_detected_genes")
  )

  # function
  p <- ggplot2::ggplot(
    samples,
    aes(
      x = .data[[group_col]],
      y = nb_detected_genes,
      fill = .data[[group_col]]
    )
  ) +
    ggplot2::geom_boxplot(alpha = 0.5) +
    ggplot2::scale_x_discrete(
      labels = function(x) stringr::str_wrap(x, width = 35)
    ) +
    ggplot2::labs(
      x = "Groups",
      y = "Number of genes detected",
      title = "Number of genes by cohort"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "none"
    )

  return(p)
}

#' @title Helper plot function for identification of outliers
#'
#' @param samples data.table with sample information with perc_detected_genes
#' and a column specifying the cohort.
#' @param group_col String specifying the column with cohort information
#' @param min_perc Numeric. Lower cutoff to identify outliers.
#' @param max_perc Numeric. Upper cutoff to identify outliers
#'
#' @return ggplot object, i.e., beeswarm plot with outlier indication
plot_preprocessing_outliers <- function(
  samples,
  group_col,
  min_perc,
  max_perc
) {
  # checks
  checkmate::assertDataFrame(samples)
  checkmate::assertNames(
    colnames(samples),
    must.include = c(group_col, "perc_detected_genes")
  )
  checkmate::qassert(min_perc, "N1")
  checkmate::qassert(max_perc, "N1")

  # function
  p <- ggplot2::ggplot(
    samples,
    aes(
      x = .data[[group_col]],
      y = perc_detected_genes,
      color = .data[[group_col]]
    )
  ) +
    ggplot2::geom_point(
      position = position_jitter(width = 0.2, height = 0, seed = 123),
      size = 3,
      alpha = 0.7
    ) +
    ggplot2::geom_hline(
      yintercept = min_perc,
      color = "red",
      linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = max_perc,
      color = "red",
      linetype = "dashed"
    ) +
    ggplot2::scale_x_discrete(
      labels = function(x) stringr::str_wrap(x, width = 35)
    ) +
    ggplot2::labs(
      x = "Groups",
      y = "Percentage of genes detected",
      title = "Outlier samples based on %genes detected"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    )

  return(p)
}

#' Helper plot function for Voom normalisation
#'
#' @param voom_object `EList`. Voom object with normalised counts.
#'
#' @return ggplot object, i.e., voom normalisation plot.
plot_voom_normalization <- function(voom_object) {
  # checks
  checkmate::assertClass(voom_object, "EList")

  # function
  voom_plot_data <- data.table::data.table(
    avg_exp = rowMeans(voom_object$E),
    residual_sd = sqrt(matrixStats::rowSds(voom_object$E))
  )

  ggplot2::ggplot(
    data = voom_plot_data,
    mapping = aes(x = avg_exp, residual_sd)
  ) +
    ggplot2::geom_point(size = 0.1) +
    ggplot2::geom_smooth(
      method = "loess",
      span = 0.5,
      se = FALSE,
      color = "red"
    ) +
    ggplot2::labs(
      x = "sqrt(count size)",
      y = "Residual standard deviation",
      title = "Voom: Mean-variance trend"
    ) +
    ggplot2::theme_classic()
}

#' Helper plot function for boxplot of normalized data
#'
#' @param samples data.table with sample information with perc_detected_genes
#' and a column specifying the cohort.
#' @param voom_object `EList`. Voom object with normalised counts.
#' @param group_col String. The grouping column.
#'
#' @return ggplot object, i.e., box plot with expression per sample.
plot_boxplot_normalization <- function(samples, voom_object, group_col) {
  # checks
  checkmate::assertClass(voom_object, "EList")
  checkmate::assertDataTable(samples)
  checkmate::qassert(group_col, "S1")
  checkmate::assertTRUE(group_col %in% names(samples))

  # function
  boxplot_data <- data.table::data.table(
    sample = rep(colnames(voom_object$E), each = nrow(voom_object$E)),
    expression = as.vector(voom_object$E),
    group = rep(as.factor(samples[[group_col]]), each = nrow(voom_object$E))
  )

  p <- ggplot2::ggplot(
    boxplot_data,
    aes(x = sample, y = expression, fill = group)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(
      x = "Samples",
      y = "Normalized log2 expression",
      title = "Boxplot of normalized gene expression across samples"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )

  return(p)
}

#' Return QC plots
#'
#' @description
#' Getter function to extract the QC plots from the [bixverse::bulk_dge()]
#' class. These are added when you run for example
#' [bixverse::qc_bulk_dge()] and [bixverse::normalise_bulk_dge()]. You can
#' either leave the plot choice as `NULL` and provide input when prompted, or
#' you provide the name. The possible plots that might be in the class
#' \itemize{
#'  \item p1_nb_genes_cohort Proportion of non-zero genes for the samples in the
#'  respective cohorts (added after using [bixverse::qc_bulk_dge()]).
#'  \item p2_outliers An outlier plot based on the data from p1, added after
#'  using [bixverse::qc_bulk_dge()].
#'  \item p3_voom_normalization Initial Voom normalisation plot after filtering
#'  lowly expressed genes. Added after using [bixverse::normalise_bulk_dge()].
#'  \item p4_boxplot_normalization Expression levels after normalisation. Added
#'  after using [bixverse::normalise_bulk_dge()].
#'  \item p5_pca_case_control A PCA plot with the chosen case control category.
#'  Added if [bixverse::calculate_pca_bulk_dge()] is run.
#'  \item p6_batch_correction_plot A PCA plot pre and post batch correction
#'  with the case-control category overlayed. Added if
#'  [bixverse::batch_correction_bulk_dge()] is run.
#' }
#'
#' @param object `bulk_dge` class.
#' @param plot_choice Optional string or integer. Index or name of the plate.
#'
#' @return Returns the DGEList stored in the class.
#'
#' @export
get_dge_qc_plot <- S7::new_generic(
  name = "get_dge_qc_plot",
  dispatch_args = "object",
  fun = function(object, plot_choice = NULL) {
    S7::S7_dispatch()
  }
)

#' @method get_dge_qc_plot bulk_dge
#'
#' @export
S7::method(get_dge_qc_plot, bulk_dge) <-
  function(object, plot_choice = NULL) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::bulk_dge"
    )
    checkmate::qassert(plot_choice, c("0", "I1", "S1"))

    plots <- S7::prop(object, "plots")
    plot_names <- names(plots)

    user_choice <- ifelse(
      is.null(plot_choice),
      select_user_option(plot_names),
      plot_choice
    )

    res <- tryCatch(
      {
        plots[[user_choice]]
      },
      error = function(e) {
        warning(paste(
          "You are trying to return a not yet generated or not existing plot.",
          "Returning NULL"
        ))
        return(NULL)
      }
    )

    # Return
    return(res)
  }

### bulk cor -------------------------------------------------------------------

#' @title Plot the highly variable genes
#'
#' @description
#' Plots the median-absolute deviation of the genes and applied thresholds.
#' Expects that [bixverse::preprocess_bulk_coexp()] was run and will throw an
#' error otherwise.
#'
#' @param object The underlying class, see [bixverse::bulk_coexp()].
#' @param bins Integer. Number of bins to plot.
#'
#' @export
plot_hvgs <- S7::new_generic(
  "plot_hvgs",
  "object",
  fun = function(object, bins = 50L) {
    S7::S7_dispatch()
  }
)

#' @method plot_hvgs bulk_coexp
#'
#' @import ggplot2
#'
#' @export
S7::method(plot_hvgs, bulk_coexp) <- function(object, bins = 50L) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(bins, "I1")
  # Early return
  if (is.null(S7::prop(object, "params")[["preprocessing"]])) {
    warning("No pre-processing data found. Returning NULL.")
    return(NULL)
  }
  plot_df <- S7::prop(object, "processed_data")[["feature_meta"]]

  p <- ggplot2::ggplot(data = plot_df, mapping = aes(x = MAD)) +
    ggplot2::geom_histogram(
      mapping = aes(fill = hvg),
      bins = 50L,
      color = "black",
      alpha = 0.7
    ) +
    ggplot2::scale_fill_manual(
      values = setNames(c("orange", "grey"), c(TRUE, FALSE))
    ) +
    ggplot2::ggtitle(
      "Distribution of MAD across the genes",
      subtitle = "And included genes"
    ) +
    ggplot2::theme_minimal()

  mad_threshold <- S7::prop(object, "params")[["preprocessing"]][[
    "mad_threshold"
  ]]

  if (mad_threshold != "not applicable") {
    p <- p +
      ggplot2::geom_vline(
        xintercept = mad_threshold,
        linetype = "dashed",
        color = "darkred"
      )
  }

  return(p)
}


#' Helper plot function for pca with contrasts
#'
#' @param pca_dt data.table. Datatable with PCA and contrast information
#' @param grps Factor or character vector. The group vector.
#'
#' @return ggplot object for the pca
plot_pca <- function(pca_dt, grps) {
  # checks
  checkmate::assertDataTable(pca_dt)
  checkmate::assertTRUE(all(c(grps, "PC_1", "PC_2") %in% colnames(pca_dt)))

  p <- ggplot2::ggplot(
    data = pca_dt,
    mapping = aes(x = PC_1, y = PC_2)
  ) +
    ggplot2::geom_point(
      mapping = aes(col = .data[[grps]]),
      size = 3,
      alpha = 0.7
    ) +
    ggplot2::xlab("PC1") +
    ggplot2::ylab("PC2") +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("PCA with key columns") +
    ggplot2::labs(colour = "Groups:")

  return(p)
}


### plots ----------------------------------------------------------------------

#' Plot the PCA data
#'
#' @param object The underlying class, see [bixverse::bulk_dge()].
#' @param pcs_to_plot String vector of length 2.
#' Will default to `c("PC_1", "PC_2")`.
#' @param cols_to_plot String vector. The columns within the meta-data to plot.
#' Defaults to `c('contrast_info', 'sample_source')`
#' @param ... additional parameters
#'
#' @return A plot if the PCA information was found. `NULL` if no PCA was found.
#'
#' @export
plot_pca_res <- S7::new_generic(
  "plot_pca_res",
  "object",
  fun = function(
    object,
    cols_to_plot = c('contrast_info', 'sample_source'),
    pcs_to_plot = c("PC_1", "PC_2"),
    ...
  ) {
    S7::S7_dispatch()
  }
)

#' @method plot_pca_res bulk_dge
#'
#' @export
#'
#' @import data.table
#' @import ggplot2
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
S7::method(plot_pca_res, bulk_dge) <- function(
  object,
  cols_to_plot = c('contrast_info', 'sample_source'),
  pcs_to_plot = c("PC_1", "PC_2"),
  ...
) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_dge")
  checkmate::qassert(cols_to_plot, "S+")
  checkmate::qassert(pcs_to_plot, "S2")

  pca_dt <- S7::prop(object, "outputs")[['pca']]
  meta_data <- S7::prop(object, "meta_data")

  if (is.null(pca_dt)) {
    warning(paste(
      "No PCA results found in the object.",
      "Did you run calculate_pca_bulk_dge()?",
      "Returning NULL."
    ))
    return(NULL)
  }

  checkmate::assertNames(
    names(pca_dt),
    must.include = pcs_to_plot
  )
  checkmate::assertNames(
    names(meta_data),
    must.include = cols_to_plot
  )

  plot_df <- data.table::merge.data.table(
    pca_dt[, c('sample_id', pcs_to_plot), with = FALSE],
    meta_data[, c('sample_id', cols_to_plot), with = FALSE],
    by.x = 'sample_id',
    by.y = 'sample_id'
  ) %>%
    data.table::melt(id.vars = c('sample_id', pcs_to_plot))

  p <- ggplot2::ggplot(
    data = plot_df,
    mapping = aes(x = .data[[pcs_to_plot[1]]], y = .data[[pcs_to_plot[2]]])
  ) +
    ggplot2::geom_point(mapping = aes(col = value), size = 3, alpha = 0.7) +
    ggplot2::facet_wrap(facets = ~variable, ncol = 3L) +
    ggplot2::xlab("PC1") +
    ggplot2::ylab("PC2") +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("PCA with key columns") +
    ggplot2::labs(colour = "Groups:")

  return(p)
}
