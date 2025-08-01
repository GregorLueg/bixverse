# methods ----------------------------------------------------------------------

## pre-processing --------------------------------------------------------------

### bulk dge -------------------------------------------------------------------

#' QC on the bulk dge data
#'
#' @description
#' This function generates general QC metrics for bulk DGE experiments.
#'
#' @param object The underlying class, see [bixverse::bulk_dge()].
#' @param group_col String. The column in the metadata that will contain the
#' contrast groups. Needs to be part of the metadata stored in the class.
#' @param norm_method String. One of
#' `c("TMM", "TMMwsp", "RLE", "upperquartile", "none")`. Please refer to
#' [edgeR::normLibSizes()].
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


### bulk cor modules -----------------------------------------------------------

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

### plots ----------------------------------------------------------------------

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

    # Return
    return(plots[[user_choice]])
  }

## pca -------------------------------------------------------------------------

#' Calculate PCA on the expression.
#'
#' @description
#' Calculates the principal component on top of the filtered count matrix and
#' adds the information of the first 10 principal components to the outputs.
#'
#' @param object The underlying class, see [bixverse::bulk_dge()].
#' @param scale_genes Boolean. Shall the log(cpm) counts be scaled prior the PCA
#' calculation. Defaults to `FALSE`.
#' @param pcs Integer. Number of PCs to return and add to the outputs slot.
#' @param no_hvg_genes Integer. Number of highly variable genes to include.
#' Defaults to 2500.
#'
#' @return Returns the class with additional data added to the outputs.
#'
#' @export
calculate_pca_bulk_dge <- S7::new_generic(
  "calculate_pca_bulk_dge",
  "object",
  fun = function(
    object,
    scale_genes = FALSE,
    pcs = 10L,
    no_hvg_genes = 2500L
  ) {
    S7::S7_dispatch()
  }
)

#' @method calculate_pca_bulk_dge bulk_dge
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
S7::method(calculate_pca_bulk_dge, bulk_dge) <- function(
  object,
  scale_genes = FALSE,
  pcs = 10L,
  no_hvg_genes = 2500L
) {
  gene_id <- hvg <- var_id <- sample_id <- . <-
    NULL

  # Checks
  checkmate::assertClass(object, "bixverse::bulk_dge")
  checkmate::qassert(scale_genes, "B1")
  checkmate::qassert(pcs, "I1")
  checkmate::qassert(no_hvg_genes, "I1")

  # Early return
  normalised_counts_present <- checkmate::testNames(
    names(S7::prop(object, "outputs")),
    must.include = "normalised_counts"
  )

  if (!normalised_counts_present) {
    warning(paste(
      "No normalised counts found. Did you run preprocess_bulk_dge()?",
      "Returning object as is"
    ))
    return(object)
  }

  normalised_counts <- S7::prop(object, "outputs")[['normalised_counts']]
  contrast_info <- S7::prop(object, "outputs")[['group_col']]
  sample_info <- S7::prop(object, "outputs")[["sample_info"]]

  hvg_data <- list(
    gene_id = rownames(normalised_counts),
    mad = matrixStats::rowMads(normalised_counts)
  ) %>%
    data.table::as.data.table() %>%
    data.table::setorder(-mad)

  hvg_genes <- hvg_data[1:no_hvg_genes, gene_id]

  if (!is.null(S7::prop(object, "variable_info"))) {
    S7::prop(object, "variable_info")[, hvg := var_id %in% hvg_genes]
  }

  input_genes <- t(normalised_counts[hvg_genes, ])
  pca_results <- rs_prcomp(input_genes, scale = scale_genes)
  pcs_to_take <- min(pcs, ncol(pca_results$scores))
  pca_dt <- pca_results$scores[, 1:pcs_to_take] %>%
    `colnames<-`(sprintf("PC_%i", 1:pcs_to_take)) %>%
    data.table::as.data.table() %>%
    .[, sample_id := rownames(input_genes)] %>%
    .[, c("sample_id", sprintf("PC_%i", 1:pcs_to_take)), with = FALSE]

  plot_df <- data.table::merge.data.table(
    pca_dt,
    sample_info[, c('sample_id', contrast_info), with = FALSE],
    by.x = 'sample_id',
    by.y = 'sample_id'
  )

  p5_pca_case_control <- ggplot2::ggplot(
    data = plot_df,
    mapping = aes(x = PC_1, y = PC_2)
  ) +
    ggplot2::geom_point(
      mapping = aes(col = .data[[contrast_info]]),
      size = 3,
      alpha = 0.7
    ) +
    ggplot2::xlab("PC1") +
    ggplot2::ylab("PC2") +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("PCA with key columns") +
    ggplot2::labs(colour = "Groups:")

  pca_anova <- check_pca_grp_differences(
    pc1 = plot_df$PC_1,
    pc2 = plot_df$PC_2,
    grps = plot_df[[contrast_info]]
  )

  ## params
  pca_params <- list(
    hvg_genes = hvg_genes,
    scale = scale,
    pcs_taken = pcs_to_take
  )

  S7::prop(object, "params")[["pca"]] <- pca_params
  S7::prop(object, "outputs")[['pca']] <- pca_dt
  S7::prop(object, "outputs")[['pca_anova']] <- pca_anova
  S7::prop(object, "plots")[['p5_pca_case_control']] <- p5_pca_case_control

  return(object)
}

### helpers --------------------------------------------------------------------

#' Helper function to check if PC1 and/or PC2 distinguish groups
#'
#' @description
#' Runs an ANOVA for the group variable vs PC1 and PC2 and checks if either
#' principal component can separate the groups.
#'
#' @param pc1 Numeric vector. The values of PC1.
#' @param pc2 Numeric vector. The values of PC2.
#' @param grps Factor or character vector. The group vector.
#'
#' @returns A data.table with two columns
#'
#' @import data.table
check_pca_grp_differences <- function(pc1, pc2, grps) {
  # Globals
  lm <- anova <- NULL

  # Checks
  checkmate::qassert(pc1, "N+")
  checkmate::qassert(pc2, "N+")
  checkmate::qassert(grps, c("F+", "S+"))
  # Function body
  pca_anova_dt <- if (length(unique(grps)) > 1) {
    pc1_anova_pval <- anova(lm(pc1 ~ grps))$`Pr(>F)`[1]
    pc2_anova_pval <- anova(lm(pc2 ~ grps))$`Pr(>F)`[1]
    data.table::data.table(
      pc = c("PC1", "PC2"),
      pvalue = c(pc1_anova_pval, pc2_anova_pval)
    )
  } else {
    data.table::data.table(
      pc = c("PC1", "PC2"),
      pvalue = c(NA, NA)
    )
  }

  return(pca_anova_dt)
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

## batch correction ------------------------------------------------------------

#' Run a linear batch correction
#'
#' @description
#' Runs a linear batch correction over the data regressing out batch effects
#' and adding `normalised_counts_corrected` to the object. Should these counts
#' be found by [bixverse::calculate_all_dges()], they will be used for
#' calculations of effect sizes based on Hedge's G.
#'
#' @param object The underlying class, see [bixverse::bulk_dge()].
#' @param contrast_column String. The contrast column in which the groupings
#' are stored. Needs to be found in the meta_data within the properties.
#' @param batch_col String. The column in which the batch effect groups can
#' be found.
#' @param scale_genes Boolean. Shall the log(cpm) counts be scaled prior the PCA
#' calculation. Defaults to `FALSE`.
#' @param no_hvg_genes Integer. Number of highly variable genes to include.
#' Defaults to 2500.
#'
#' @return Returns the class with additional data added to the outputs.
#'
#' @export
batch_correction_bulk_dge <- S7::new_generic(
  "batch_correction_bulk_dge",
  "object",
  fun = function(
    object,
    contrast_column,
    batch_col,
    scale_genes = FALSE,
    no_hvg_genes = 2500L
  ) {
    S7::S7_dispatch()
  }
)

#' @method batch_correction_bulk_dge bulk_dge
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
#' @import patchwork
#' @import ggplot2
S7::method(batch_correction_bulk_dge, bulk_dge) <- function(
  object,
  contrast_column,
  batch_col,
  scale_genes = FALSE,
  no_hvg_genes = 2500L
) {
  gene_id <- . <- sample_id <- NULL

  # Checks
  checkmate::assertClass(object, "bixverse::bulk_dge")
  checkmate::qassert(contrast_column, "S1")
  checkmate::qassert(batch_col, "S1")
  checkmate::qassert(scale_genes, "B1")
  checkmate::qassert(no_hvg_genes, "I1")

  # Body
  # Early return
  normalised_counts_present <- checkmate::testNames(
    names(S7::prop(object, "outputs")),
    must.include = "normalised_counts"
  )

  if (!normalised_counts_present) {
    warning(paste(
      "No normalised counts found. Did you run preprocess_bulk_dge()?",
      "Returning object as is"
    ))
    return(object)
  }

  pca_present <- !is.null(S7::prop(object, "outputs")[['pca']])

  if (!pca_present) {
    message("No PCA data found. Running PCA with default parameters now.")
    object <- calculate_pca_bulk_dge(object)
  }

  normalised_counts <- S7::prop(object, "outputs")[['normalised_counts']]
  sample_info <- S7::prop(object, "outputs")[["sample_info"]]
  pca_dt <- S7::prop(object, "outputs")[['pca']]

  batch_data <- factor(sample_info[[batch_col]])

  design <- model.matrix(~ 0 + factor(sample_info[[contrast_column]]))
  colnames(design) = unique(sample_info[[contrast_column]])

  normalised_counts_corrected <- limma::removeBatchEffect(
    x = normalised_counts,
    batch = batch_data,
    design = design
  )

  # TODO abstract out the HVG genes and PCA calculations into one single
  # function

  hvg_data <- list(
    gene_id = rownames(normalised_counts_corrected),
    mad = matrixStats::rowMads(normalised_counts_corrected)
  ) %>%
    data.table::as.data.table() %>%
    data.table::setorder(-mad)

  hvg_genes <- hvg_data[1:no_hvg_genes, gene_id]

  input_genes <- t(normalised_counts_corrected[hvg_genes, ])
  pca_results <- rs_prcomp(input_genes, scale = scale_genes)

  pca_dt_cor <- pca_results$scores[, 1:2] %>%
    `colnames<-`(sprintf("PC_%i", 1:2)) %>%
    data.table::as.data.table() %>%
    .[, sample_id := rownames(input_genes)] %>%
    merge(
      .,
      sample_info[, c("sample_id", contrast_column), with = FALSE],
      by = 'sample_id'
    )

  pca_dt_uncor <- pca_dt[, c("PC_1", "PC_2", "sample_id"), with = FALSE] %>%
    merge(
      .,
      sample_info[, c("sample_id", contrast_column), with = FALSE],
      by = 'sample_id'
    )

  # Otherwise it continues bugging...
  library(patchwork)

  plot_uncor <- ggplot2::ggplot(
    data = pca_dt_uncor,
    mapping = aes(x = PC_1, y = PC_2)
  ) +
    ggplot2::geom_point(
      mapping = aes(col = .data[[contrast_column]]),
      size = 3,
      alpha = 0.7
    ) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Pre batch correction") +
    ggplot2::xlab("PC1") +
    ggplot2::ylab("PC2") +
    ggplot2::theme(legend.position = "none")

  plot_cor <- ggplot2::ggplot(
    data = pca_dt_cor,
    mapping = aes(x = PC_1, y = PC_2)
  ) +
    ggplot2::geom_point(
      mapping = aes(col = .data[[contrast_column]]),
      size = 3,
      alpha = 0.7
    ) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Post batch correction") +
    ggplot2::xlab("PC1") +
    ggplot2::ylab("PC2") +
    ggplot2::theme(legend.position = "bottom")

  p6_batch_correction_plot <- plot_uncor +
    plot_cor +
    patchwork::plot_annotation(
      title = 'PCA plots pre and post batch effect correction',
      subtitle = 'Batch effect correction via limma::removeBatchEffect()'
    )

  ## params
  batch_cor_params <- list(
    batch_effect_col = batch_col,
    batches = batch_data
  )

  S7::prop(object, "outputs")[[
    'normalised_counts_corrected'
  ]] <- normalised_counts_corrected
  S7::prop(object, "plots")[[
    'p6_batch_correction_plot'
  ]] <- p6_batch_correction_plot
  S7::prop(object, "params")[['batch_cor_params']] <- batch_cor_params

  return(object)
}

## differential gene expression ------------------------------------------------

#' Calculate all possible DGE variants
#'
#' @description
#' This class is a wrapper class around applying various differential gene
#' expression on the data. At a minimum you will need to provide a
#' `contrast_column` that can be found in the meta-data. If you do not provide
#' a vector of contrasts that you wish to test for, every permutation of
#' groups represented in that column will be tested against each other. Two
#' variants of differential gene expression analyses will be applied:
#' \enumerate{
#'   \item A standard Limma Voom approach will be applied. For this specific
#'   approach you can also provide co-variates that will be used in the model
#'   fitting. The function will return the Limma Voom results for every
#'   combination of contrasts found. For more details, please refer to
#'   [bixverse::run_limma_voom()].
#'   \item Secondly, the Hedge's effect size for each gene between all
#'   combinations of the groups will be calculated. The effect sizes can
#'   subsequently be used for e.g., meta-analyses across various studies of
#'   interest.
#' }
#'
#' @param object The underlying class, see [bixverse::bulk_dge()].
#' @param contrast_column String. The contrast column in which the groupings
#' are stored. Needs to be found in the meta_data within the properties.
#' @param filter_column Optional String. If there is a column you wish to use as
#' sub groupings, this can be provided here. An example could be different
#' sampled tissues and you wish to run the DGE analyses within each tissue
#' separately.
#' @param co_variates Optional string vector. Any co-variates you wish to
#' consider during the Limma Voom modelling.
#' @param contrast_list Optional string vector. A vectors that contains the
#' contrast formatted as `"contrast1-contrast2"`. Default `NULL` will create
#' all possible contrast automatically.
#' @param ... Additional parameters to forward to [limma::eBayes()] or
#' [limma::voom()].
#' @param small_sample_correction Can be NULL (automatic determination if a
#' small sample size correction should be applied) or a Boolean.
#' @param .verbose Controls verbosity of the function.
#'
#' @return Returns the class with additional data added to the outputs.
#'
#' @export
calculate_all_dges <- S7::new_generic(
  "calculate_all_dges",
  "object",
  fun = function(
    object,
    contrast_column,
    contrast_list = NULL,
    filter_column = NULL,
    co_variates = NULL,
    small_sample_correction = NULL,
    ...,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calculate_all_dges bulk_dge
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
S7::method(calculate_all_dges, bulk_dge) <- function(
  object,
  contrast_column,
  contrast_list = NULL,
  filter_column = NULL,
  co_variates = NULL,
  small_sample_correction = NULL,
  ...,
  .verbose = TRUE
) {
  . <- subgroup <- NULL

  # First checks
  checkmate::assertClass(object, "bixverse::bulk_dge")
  checkmate::qassert(contrast_column, "S+")
  checkmate::qassert(co_variates, c("S+", "0"))
  checkmate::qassert(filter_column, c("S+", "0"))

  # Early return
  dge_list_present <- checkmate::testNames(
    names(S7::prop(object, "outputs")),
    must.include = "dge_list"
  )

  if (!dge_list_present) {
    warning(paste(
      "No dge_list found. Did you run preprocess_bulk_dge()?",
      "Returning object as is"
    ))
    return(object)
  }

  ## get objects
  all_specified_columns <- setdiff(
    c(contrast_column, co_variates, filter_column),
    NULL
  )
  sample_info <- S7::prop(object, "outputs")[["sample_info"]]
  if (is.null(S7::prop(object, "outputs")[['dge_list']])) {
    dge_list = NULL
  } else {
    dge_list <- S7::prop(object, "outputs")[['dge_list']]
  }

  checkmate::assertTRUE(all(all_specified_columns %in% colnames(sample_info)))

  counts_batch_cor <- S7::prop(object, "outputs")[[
    'normalised_counts_corrected'
  ]]

  norm_counts <- if (is.null(counts_batch_cor)) {
    S7::prop(object, "outputs")[['normalised_counts']]
  } else {
    if (.verbose) {
      message(paste(
        "Found batch corrected counts.",
        "These will be used for effect size calculations"
      ))
    }
    counts_batch_cor
  }

  if (is.null(filter_column)) {
    if (.verbose) {
      message("Calculating the differential expression with limma results.")
    }

    limma_results_final <- run_limma_voom(
      meta_data = sample_info,
      main_contrast = contrast_column,
      dge_list = dge_list,
      normalised_counts = norm_counts,
      contrast_list = contrast_list,
      co_variates = co_variates,
      ...,
      .verbose = .verbose
    ) %>%
      .[, subgroup := NA]

    if (.verbose) {
      message("Calculating the Hedge's G effect sizes.")
    }
    hedges_g_results_final <- hedges_g_dge(
      meta_data = sample_info,
      main_contrast = contrast_column,
      contrast_list = contrast_list,
      normalised_counts = norm_counts,
      small_sample_correction = small_sample_correction,
      .verbose = .verbose
    ) %>%
      .[, subgroup := NA]
  } else {
    if (.verbose) {
      message(paste(
        "Filtering column provided.",
        "Method will run Limma Voom and Hedge's G on the individual data sets."
      ))
    }

    # Iterate through the grps
    groups <- unique(sample_info[[filter_column]])
    results <- purrr::map(groups, \(group) {
      # Filter the meta data and dge list
      sample_info_red <- sample_info[
        eval(parse(
          text = paste0(filter_column, " == '", group, "'")
        ))
      ]
      if (is.null(dge_list)) {
        dge_list_red = NULL
        norm_counts_red <- norm_counts[, sample_info_red$sample_id]
      } else {
        dge_list_red <- dge_list[, sample_info_red$sample_id]
        norm_counts_red <- norm_counts[, sample_info_red$sample_id]
      }

      # Limma Voom
      limma_results <- run_limma_voom(
        meta_data = data.table::copy(sample_info_red),
        main_contrast = contrast_column,
        dge_list = dge_list_red,
        normalised_counts = norm_counts_red,
        co_variates = co_variates,
        contrast_list = contrast_list,
        ...,
        .verbose = .verbose
      ) %>%
        .[, subgroup := group]

      hedges_g_results <- hedges_g_dge(
        meta_data = data.table::copy(sample_info_red),
        main_contrast = contrast_column,
        normalised_counts = norm_counts_red,
        contrast_list = contrast_list,
        small_sample_correction = small_sample_correction,
        .verbose = .verbose
      ) %>%
        .[, subgroup := group]

      return(list(
        limma_results = limma_results,
        hedges_g_results = hedges_g_results
      ))
    })

    # Rbind the data
    limma_results_final <- purrr::map(
      results,
      ~ {
        .[['limma_results']]
      }
    ) %>%
      data.table::rbindlist()

    hedges_g_results_final <- purrr::map(
      results,
      ~ {
        .[['hedges_g_results']]
      }
    ) %>%
      data.table::rbindlist()
  }

  dge_params <- list(
    contrast_column = contrast_column,
    filter_column = filter_column,
    co_variates = co_variates
  )

  S7::prop(object, "outputs")[['limma_voom_res']] <- limma_results_final
  S7::prop(object, "outputs")[['hedges_g_res']] <- hedges_g_results_final
  S7::prop(object, "params")[["dge"]] <- dge_params

  return(object)
}
