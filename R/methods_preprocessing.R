# pre-processing ---------------------------------------------------------------

## bulk dge --------------------------------------------------------------------

#' QC on the bulk dge data
#'
#' @description
#' This function generates general QC metrics for bulk DGE experiments.
#'
#' @param object The underlying class, see [bixverse::bulk_dge()].
#' @param group_col String. The column in the metadata that will contain the
#' contrast groups. Needs to be part of the metadata stored in the class.
#' @param norm_method String. One of `c("TMM", "TMMwsp", "RLE", "upperquartile",`
#' ` "none")`. Please refer to [edgeR::normLibSizes()].
#' @param outlier_threshold Float. Number of standard deviations in terms of
#' percentage genes detected you allow before removing a sample. Defaults to `2`.
#' @param min_prop Float. Minimum proportion of samples in which the gene has
#' to be identified in.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return Returns the class with the `processed_data` data slot populated and
#' applied parameters added to the `params` slot.
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
preprocess_bulk_dge <- S7::new_generic(
  "preprocess_bulk_dge",
  "object",
  fun = function(
    object,
    group_col,
    norm_method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
    outlier_threshold = 2,
    min_prop = 0.2,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)


#' @method preprocess_bulk_dge bulk_dge
#'
#' @import ggplot2
#'
#' @export
S7::method(preprocess_bulk_dge, bulk_dge) <- function(
  object,
  group_col,
  norm_method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
  outlier_threshold = 2,
  min_prop = 0.2,
  .verbose = TRUE
) {
  # Scope checks
  . <- `:=` <- NULL
  norm_method <- match.arg(norm_method)
  # Checks
  checkmate::assertClass(
    object,
    "bixverse::bulk_dge"
  )
  checkmate::qassert(group_col, "S1")
  checkmate::assertChoice(
    norm_method,
    c("TMM", "TMMwsp", "RLE", "upperquartile", "none")
  )
  checkmate::qassert(outlier_threshold, "N1")
  checkmate::qassert(min_prop, "R1[0,1]")
  checkmate::qassert(.verbose, "B1")

  raw_counts <- S7::prop(object, "raw_counts")
  meta_data <- S7::prop(object, "meta_data")

  checkmate::assertTRUE(group_col %in% colnames(meta_data))

  # Sample outlier removal
  if (.verbose) {
    message("Detecting sample outliers.")
  }
  detected_genes_nb <- data.table::data.table(
    sample_id = colnames(raw_counts),
    nb_detected_genes = matrixStats::colSums2(raw_counts > 0)
  )

  samples <- merge(meta_data, detected_genes_nb, by = "sample_id") %>%
    .[, `:=`(perc_detected_genes = nb_detected_genes / nrow(raw_counts))]

  p1_nb_genes_cohort <- ggplot2::ggplot(
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

  sd_samples = sd(samples$perc_detected_genes, na.rm = TRUE)
  min_perc = mean(samples$perc_detected_genes, na.rm = TRUE) -
    outlier_threshold * sd_samples
  max_perc = mean(samples$perc_detected_genes, na.rm = TRUE) +
    outlier_threshold * sd_samples

  p2_outliers <- ggplot2::ggplot(
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

  outliers <- samples$perc_detected_genes <= min_perc
  if (.verbose) {
    message(sprintf(
      "A total of %i samples are detected as outlier.",
      sum(outliers)
    ))
  }

  samples_red <- samples[!(outliers), ]
  raw_counts <- raw_counts[, unique(samples_red$sample_id)]

  ## Voom normalization
  if (.verbose) {
    message("Removing lowly expressed genes and normalising counts.")
  }
  dge_list <- edgeR::DGEList(raw_counts)

  # Filter lowly expressed genes
  to_keep <- edgeR::filterByExpr(
    y = dge_list,
    min.prop = min_prop,
    group = samples[[group_col]]
  )
  if (.verbose) {
    message(sprintf("A total of %i genes are kept.", sum(to_keep)))
  }

  dge_list <- edgeR::calcNormFactors(
    dge_list[to_keep, ],
    method = norm_method
  )

  samples <- samples_red[[group_col]]

  design <- model.matrix(~ 0 + samples)

  voom_obj <- limma::voom(
    counts = dge_list,
    design = design,
    normalize.method = "quantile",
    plot = TRUE
  )

  p3_voom_normalization <- recordPlot()
  dev.off()

  boxplot(
    voom_obj$E,
    col = as.factor(samples),
    main = "Boxplot of normalized gene expression across samples"
  )
  p4_boxplot_normalization <- recordPlot()
  dev.off()

  S7::prop(object, "plots") <- list(
    p1_nb_genes_cohort = p1_nb_genes_cohort,
    p2_outliers = p2_outliers,
    p3_voom_normalization = p3_voom_normalization,
    p4_boxplot_normalization = p4_boxplot_normalization
  )

  S7::prop(object, "outputs") <- list(
    dge_list = dge_list,
    sample_info = samples_red,
    normalised_counts = voom_obj$E,
    group_col = group_col
  )

  S7::prop(object, "params")[["QC_params"]] <- list(
    normalisation_method = norm_method,
    min_prop = min_prop,
    outlier_threshold = outlier_threshold
  )

  return(object)
}

## bulk_coexp ------------------------------------------------------------------

#' Process the raw data
#'
#' @description
#' Function to do general pre-processing on top of the [bixverse::bulk_coexp()].
#' Options to do scaling, HVG selection, etc.
#'
#' @param object The underlying class, see [bixverse::bulk_coexp()].
#' @param hvg Integer or float. If an integer, the top `hvg` genes will be
#' included; if float, the float has to be between 0 and 1, representing the
#' percentage of genes to include.
#' @param mad_threshold Float. Instead of of selecting number or proportion of
#' genes, you can also provide a mad_threshold.
#' @param scaling Boolean. Shall the data be scaled.
#' @param scaling_type String. You have the option to use normal scaling or
#' robust scaling.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return Returns the class with the `processed_data` data slot populated and
#' applied parameters added to the `params` slot.
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
preprocess_bulk_coexp <- S7::new_generic(
  "preprocess_bulk_coexp",
  "object",
  fun = function(
    object,
    hvg = NULL,
    mad_threshold = NULL,
    scaling = FALSE,
    scaling_type = c("normal", "robust"),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method preprocess_bulk_coexp bulk_coexp
#' @export
S7::method(preprocess_bulk_coexp, bulk_coexp) <- function(
  object,
  hvg = NULL,
  mad_threshold = NULL,
  scaling = FALSE,
  scaling_type = c("normal", "robust"),
  .verbose = TRUE
) {
  # Scope checks
  feature_name <- MAD <- NULL
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(mad_threshold, c("R1", "0"))
  nfeatures <- S7::prop(object, "params")[["original_dim"]][2]
  checkmate::qassert(hvg, c("R1[0,1]", sprintf("I1[0,%i]", nfeatures), "0"))
  checkmate::qassert(scaling, "B1")
  if (scaling) {
    checkmate::assertChoice(scaling_type, c("normal", "robust"))
  }

  mat <- S7::prop(object, "raw_data")

  feature_meta <- data.table::data.table(
    feature_name = colnames(mat),
    mean_exp = colMeans(mat),
    MAD = matrixStats::colMads(mat),
    var_exp = matrixStats::colVars(mat)
  ) %>%
    data.table::setorder(-MAD)

  if (!is.null(hvg) & !is.null(mad_threshold)) {
    choice <- menu(
      choices = c("MAD threshold", "HVG threshold"),
      title = paste(
        "You have provided both a MAD and HVG",
        "threshold. Which one do you want to use?"
      )
    )
    if (choice == 1) {
      hvg <- NULL
    } else {
      mad_threshold <- NULL
    }
  }

  if (is.null(hvg) & is.null(mad_threshold)) {
    hvg <- 1
  }

  if (!is.null(hvg)) {
    no_genes_to_take <-
      ifelse(is.integer(hvg), hvg, ceiling(hvg * ncol(mat)))
    hvg_genes <- feature_meta[1:no_genes_to_take, feature_name]
  } else {
    hvg_genes <- feature_meta[MAD >= mad_threshold, feature_name]
  }

  if (.verbose) {
    message(sprintf("A total of %i genes will be included.", length(hvg_genes)))
  }

  feature_meta[, hvg := feature_name %in% hvg_genes]

  # Process the matrix
  matrix_processed <- mat[, hvg_genes]

  if (scaling) {
    fun <-
      ifelse(scaling_type == "normal", "scale", "bixverse::robust_scale")
    matrix_processed <- rlang::eval_tidy(rlang::quo(apply(
      matrix_processed,
      1,
      !!!rlang::parse_exprs(fun)
    )))
  }

  processing_params <- list(
    mad_threshold = if (is.null(mad_threshold)) {
      "not applicable"
    } else {
      mad_threshold
    },
    hvg = if (is.null(hvg)) {
      "not applicable"
    } else {
      hvg
    },
    scaling = scaling,
    scaling_type = if (!scaling) {
      "not applicable"
    } else {
      scaling_type
    }
  )

  S7::prop(object, "params")[["preprocessing"]] <-
    processing_params
  S7::prop(object, "processed_data")[["processed_data"]] <-
    matrix_processed
  S7::prop(object, "processed_data")[["feature_meta"]] <-
    feature_meta

  object
}

# plots ------------------------------------------------------------------------

## bulk_coexp ------------------------------------------------------------------

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

## bulk dge --------------------------------------------------------------------

#' Return QC plots
#'
#' @description
#' Getter function to extract the QC plots from the [bixverse::bulk_dge()]
#' class. These are added when you run for example
#' [bixverse::preprocess_bulk_dge()]. You can either leave the plot choice as
#' `NULL` and provide input when prompted, or you provide the name. The possible
#' plots that might be in the class
#' \itemize{
#'  \item p1_nb_genes_cohort Proportion of non-zero genes for the samples in the
#'  respective cohorts
#'  \item p2_outliers An outlier plot based on the data from p1
#'  \item p3_voom_normalization Initial Voom normalisation plot after filtering
#'  lowly expressed genes
#'  \item p4_boxplot_normalization Expression levels after normalisation
#'  \item p5_pca_case_control A PCA plot with the chosen case control category.
#'  Added if [bixverse::calculate_pca_bulk_dge()] is run.
#'  \item p6_batch_correction_plot A PCA plot pre and post batch correction
#'  with the case-control category overlayed. Added if
#' [bixverse::batch_correction_bulk_dge()] is run.
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

    # Return
    return(plots[[user_choice]])
  }
