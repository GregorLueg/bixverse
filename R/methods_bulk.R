# methods ----------------------------------------------------------------------

## pre-processing --------------------------------------------------------------

### bulk dge -------------------------------------------------------------------

#' QC on the bulk dge data
#'
#' @description
#' This function will do QC on the bulk data and remove outlier samples that
#' show substantially lower expression of genes compared to the rest of the
#' data and remove lowly expressed genes.
#'
#' @param object The underlying class, see [bixverse::bulk_dge()].
#' @param group_col String. The column in the metadata that will contain the
#' contrast groups. Needs to be part of the metadata stored in the class.
#' @param outlier_threshold Float. Number of standard deviations in terms of
#' percentage genes detected you allow before removing a sample. Defaults to `2`.
#' @param min_prop Float. Minimum proportion of samples in which the gene has
#' to be identified in.
#' @param min_count Float. Minimum number of counts (cpm) to be detected in
#' min_prop of the samples (in cohorts defined by groups_coll)
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
qc_bulk_dge <- S7::new_generic(
  "qc_bulk_dge",
  "object",
  fun = function(
    object,
    group_col,
    outlier_threshold = 2,
    min_prop = 0.2,
    min_count = 10,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)


#' @method qc_bulk_dge bulk_dge
#'
#' @import ggplot2
#'
#' @export
S7::method(qc_bulk_dge, bulk_dge) <- function(
  object,
  group_col,
  outlier_threshold = 2,
  min_prop = 0.2,
  min_count = 10,
  .verbose = TRUE
) {
  # Scope checks
  . <- `:=` <- NULL
  # Checks
  checkmate::assertClass(
    object,
    "bixverse::bulk_dge"
  )
  checkmate::qassert(group_col, "S1")
  checkmate::qassert(outlier_threshold, "N1")
  checkmate::qassert(min_prop, "R1[0,1]")
  checkmate::qassert(.verbose, "B1")

  raw_counts <- S7::prop(object, "raw_counts")
  meta_data <- S7::prop(object, "meta_data")

  checkmate::assertTRUE(group_col %in% colnames(meta_data))

  if (.verbose) {
    message("Detecting sample outliers.")
  }
  detected_genes_nb <- data.table::data.table(
    sample_id = colnames(raw_counts),
    nb_detected_genes = matrixStats::colSums2(raw_counts > 0)
  )

  samples <- merge(meta_data, detected_genes_nb, by = "sample_id") %>%
    .[, `:=`(perc_detected_genes = nb_detected_genes / nrow(raw_counts))]

  # plot 1 - number of genes per cohort
  p1_nb_genes_cohort <- plot_preprocessing_genes(samples, group_col)

  ## outlier detection
  sd_samples = sd(samples$perc_detected_genes, na.rm = TRUE)
  min_perc = mean(samples$perc_detected_genes, na.rm = TRUE) -
    outlier_threshold * sd_samples
  max_perc = mean(samples$perc_detected_genes, na.rm = TRUE) +
    outlier_threshold * sd_samples
  outliers <- samples$perc_detected_genes <= min_perc

  # plot 2 - outlier plot based on standard deviations
  p2_outliers <- plot_preprocessing_outliers(
    samples,
    group_col,
    min_perc,
    max_perc
  )

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
    message("Removing lowly expressed genes.")
  }
  dge_list <- edgeR::DGEList(raw_counts)

  # Filter lowly expressed genes
  to_keep <- edgeR::filterByExpr(
    y = dge_list,
    min.prop = min_prop,
    min.count = min_count,
    group = samples[[group_col]]
  )
  dge_list <- dge_list[to_keep, ]
  count_matrix <- dge_list$counts

  if (.verbose) {
    message(sprintf("A total of %i genes are kept.", sum(to_keep)))
  }

  S7::prop(object, "plots") <- list(
    p1_nb_genes_cohort = p1_nb_genes_cohort,
    p2_outliers = p2_outliers
  )

  S7::prop(object, "outputs") <- list(
    dge_list = dge_list,
    sample_info = samples_red,
    group_col = group_col,
    raw_counts_filtered = raw_counts
  )

  S7::prop(object, "params")[["QC_params"]] <- list(
    min_prop = min_prop,
    outlier_threshold = outlier_threshold
  )

  return(object)
}

#' Normalise the count data for DGE.
#'
#' @description
#' This function will apply the CPM + Voom normalisation and can additionally
#' calculate TPM and FPKM values for plotting purposes.
#'
#' @param object The underlying class, see [bixverse::bulk_dge()].
#' @param group_col String. The column in the metadata that will contain the
#' contrast groups. Needs to be part of the metadata stored in the class.
#' @param norm_method String. One of
#' `c("TMM", "TMMwsp", "RLE", "upperquartile", "none")`. Please refer to
#' [edgeR::normLibSizes()].
#' @param calc_tpm Boolean. Output TPM calculation (default = FALSE).
#' @param calc_fpkm Boolean. Output FPKM calculation (default = FALSE).
#' @param gene_lengths Optional named numeric. If you want to calculate TPM
#' or FPKM you need to provide this one. The names need to be the same
#' identifier as used in the counts.
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
normalise_bulk_dge <- S7::new_generic(
  "normalise_bulk_dge",
  "object",
  fun = function(
    object,
    group_col,
    norm_method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
    calc_tpm = FALSE,
    calc_fpkm = FALSE,
    gene_lengths = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)


#' @method normalise_bulk_dge bulk_dge
#'
#' @import ggplot2
#'
#' @export
S7::method(normalise_bulk_dge, bulk_dge) <- function(
  object,
  group_col,
  norm_method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
  calc_tpm = FALSE,
  calc_fpkm = FALSE,
  gene_lengths = NULL,
  .verbose = TRUE
) {
  norm_method <- match.arg(norm_method)

  # scope checks
  . <- `:=` <- NULL

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
  checkmate::qassert(calc_tpm, "B1")
  checkmate::qassert(calc_fpkm, "B1")
  checkmate::qassert(.verbose, "B1")

  # early return
  if (is.null(S7::prop(object, "outputs")[["dge_list"]])) {
    warning(paste(
      "Could not find the DGE list in the object. Did you run qc_bulk_dge()?",
      "Returning class as is."
    ))
    return(object)
  }

  dge_list <- S7::prop(object, "outputs")[["dge_list"]]
  raw_counts_filtered <- S7::prop(object, "outputs")[["raw_counts_filtered"]]
  meta_data_filtered <- S7::prop(object, "outputs")[["sample_info"]]

  # assertion at this point for calculations here
  if (calc_tpm | calc_fpkm) {
    checkmate::assertNumeric(gene_lengths, names = "named", any.missing = FALSE)
    checkmate::assertTRUE(all(
      rownames(raw_counts_filtered) == names(gene_lengths)
    ))
  }

  ## CPM
  dge_list <- edgeR::calcNormFactors(
    dge_list,
    method = norm_method
  )
  groups <- meta_data_filtered[[group_col]]
  design <- model.matrix(~ 0 + groups)
  voom_obj <- limma::voom(
    counts = dge_list,
    design = design,
    normalize.method = "quantile",
    plot = FALSE
  )

  # plot 3 - mean variance trend from voom
  p3_voom_normalization <- plot_voom_normalization(voom_object = voom_obj)

  # plot 4 - boxplots
  p4_boxplot_normalization <- plot_boxplot_normalization(
    samples = meta_data_filtered,
    voom_object = voom_obj,
    group_col = group_col
  )

  ## TPM and FPKM

  tpm_counts <- if (calc_tpm) {
    calculate_tpm(counts = raw_counts_filtered, gene_lengths = gene_lengths)
  } else {
    NULL
  }

  fpkm_counts <- if (calc_fpkm) {
    calculate_rpkm(counts = raw_counts_filtered, gene_lengths = gene_lengths)
  } else {
    NULL
  }

  S7::prop(object, "plots")[["p3_voom_normalization"]] <- p3_voom_normalization
  S7::prop(object, "plots")[[
    "p4_boxplot_normalization"
  ]] <- p4_boxplot_normalization

  S7::prop(object, "outputs")[["normalised_counts"]] <- voom_obj$E
  S7::prop(object, "outputs")[["tpm_counts"]] <- tpm_counts
  S7::prop(object, "outputs")[["fpkm_counts"]] <- fpkm_counts

  S7::prop(object, "params")[["norm"]] <- list(
    normalisation_method = norm_method
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

  p5_pca_case_control <- plot_pca(
    pca_dt = plot_df,
    grps = contrast_info
  )

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

## batch correction ------------------------------------------------------------

#' Run a linear batch correction
#'
#' @description
#' Runs a linear batch correction over the data regressing out batch effects
#' and adding `normalised_counts_corrected` to the object. Should these counts
#' be found by [bixverse::calculate_dge_hedges()], they will be used for
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

  plot_uncor <- plot_pca(
    pca_dt = pca_dt_uncor,
    grps = contrast_column
  ) +
    ggplot2::ggtitle("Pre batch correction")
  plot_cor <- plot_pca(
    pca_dt = pca_dt_cor,
    grps = contrast_column
  ) +
    ggplot2::ggtitle("Post batch correction")

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

#' Calculates the Limma Voom DGE
#'
#' @description
#' This function will apply the Limma Voom DGE workflow. At a minimum you will
#' need to provide `contrast_column` that can be found in the meta-data. If you
#' do not provide a vector of contrasts that you wish to test for, every
#' permutation of groups represented in that column will be tested against each
#' other.
#'
#' @param object The underlying class, see [bixverse::bulk_dge()].
#' @param contrast_column String. The contrast column in which the groupings
#' are stored. Needs to be found in the meta_data within the properties.
#' @param contrast_list Optional string vector. A vectors that contains the
#' contrast formatted as `"contrast1-contrast2"`. Default `NULL` will create
#' all possible contrast automatically.
#' @param filter_column Optional String. If there is a column you wish to use as
#' sub groupings, this can be provided here. An example could be different
#' sampled tissues and you wish to run the DGE analyses within each tissue
#' separately in the data.
#' @param co_variates Optional string vector. Any co-variates you wish to
#' consider during the Limma Voom modelling.
#' @param quantile_norm Boolean. Shall the data also be quantile normalised.
#' Defaults to `FALSE`.
#' @param ... Additional parameters to forward to [limma::eBayes()] or
#' [limma::voom()].
#' @param .verbose Controls verbosity of the function.
#'
#' @return Returns the class with additional data added to the outputs.
#'
#' @export
calculate_dge_limma <- S7::new_generic(
  "calculate_dge_limma",
  "object",
  fun = function(
    object,
    contrast_column,
    contrast_list = NULL,
    filter_column = NULL,
    co_variates = NULL,
    quantile_norm = FALSE,
    ...,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calculate_dge_limma bulk_dge
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
S7::method(calculate_dge_limma, bulk_dge) <- function(
  object,
  contrast_column,
  contrast_list = NULL,
  filter_column = NULL,
  co_variates = NULL,
  quantile_norm = FALSE,
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
      "No dge_list found. Did you run qc_bulk_dge()?",
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
  dge_list <- S7::prop(object, "outputs")[["dge_list"]]

  checkmate::assertTRUE(all(all_specified_columns %in% colnames(sample_info)))

  if (is.null(filter_column)) {
    if (.verbose) {
      message("Calculating the differential expression with Limma Voom.")
    }

    limma_results_final <- run_limma_voom(
      meta_data = sample_info,
      main_contrast = contrast_column,
      dge_list = dge_list,
      contrast_list = contrast_list,
      co_variates = co_variates,
      quantile_norm = quantile_norm,
      ...,
      .verbose = .verbose
    ) %>%
      .[, subgroup := NA]
  } else {
    if (.verbose) {
      message(paste(
        "Filtering column provided. Calculating the differential expression",
        "with Limma Voom across the groups."
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
      dge_list_red <- dge_list[, sample_info_red$sample_id]

      # Limma Voom
      limma_results <- run_limma_voom(
        meta_data = sample_info,
        main_contrast = contrast_column,
        dge_list = dge_list,
        contrast_list = contrast_list,
        co_variates = co_variates,
        quantile_norm = quantile_norm,
        ...,
        .verbose = .verbose
      ) %>%
        .[, subgroup := group]

      return(limma_results)
    })

    # rbind the data
    limma_results_final <- data.table::rbindlist(limma_results)
  }

  dge_params <- list(
    contrast_column = contrast_column,
    filter_column = filter_column,
    co_variates = co_variates
  )

  S7::prop(object, "outputs")[['limma_voom_res']] <- limma_results_final
  S7::prop(object, "params")[["limma_dge"]] <- dge_params

  return(object)
}

#' Calculates the Hedge's G effect size
#'
#' @description
#' This function will calculate the Hedge's G effect size on the normalised
#' counts. Should batch-corrected counts be found, these will be used. At a
#' minimum you will need to provide `contrast_column` that can be found in the
#' meta-data. If you do not provide a vector of contrasts that you wish to test
#' for, every permutation of groups represented in that column will be tested
#' against each other.
#'
#' @param object The underlying class, see [bixverse::bulk_dge()].
#' @param contrast_column String. The contrast column in which the groupings
#' are stored. Needs to be found in the meta_data within the properties.
#' @param contrast_list Optional string vector. A vectors that contains the
#' contrast formatted as `"contrast1-contrast2"`. Default `NULL` will create
#' all possible contrast automatically.
#' @param filter_column Optional String. If there is a column you wish to use as
#' sub groupings, this can be provided here. An example could be different
#' sampled tissues and you wish to run the DGE analyses within each tissue
#' separately in the data.
#' @param .verbose Controls verbosity of the function.
#'
#' @return Returns the class with additional data added to the outputs.
#'
#' @export
calculate_dge_hedges <- S7::new_generic(
  "calculate_dge_hedges",
  "object",
  fun = function(
    object,
    contrast_column,
    contrast_list = NULL,
    filter_column = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calculate_dge_hedges bulk_dge
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
S7::method(calculate_dge_hedges, bulk_dge) <- function(
  object,
  contrast_column,
  contrast_list = NULL,
  filter_column = NULL,
  .verbose = TRUE
) {
  . <- subgroup <- NULL

  # First checks
  checkmate::assertClass(object, "bixverse::bulk_dge")
  checkmate::qassert(contrast_column, "S+")
  checkmate::qassert(filter_column, c("S+", "0"))

  # Early return
  norm_count_present <- checkmate::testNames(
    names(S7::prop(object, "outputs")),
    must.include = "normalised_counts"
  )

  if (!norm_count_present) {
    warning(paste(
      "No normalised found. Did you run preprocess_bulk_dge()",
      "and normalise_bulk_dge()? Returning object as is."
    ))
    return(object)
  }

  ## get objects
  all_specified_columns <- setdiff(
    c(contrast_column, filter_column),
    NULL
  )
  sample_info <- S7::prop(object, "outputs")[["sample_info"]]

  checkmate::assertTRUE(all(all_specified_columns %in% names(sample_info)))

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
      message("Calculating the differential expression based on Hedge's G.")
    }

    hedges_g_results_final <- hedges_g_dge(
      meta_data = sample_info,
      main_contrast = contrast_column,
      contrast_list = contrast_list,
      normalised_counts = norm_counts,
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
      norm_counts_red <- norm_counts[, sample_info_red$sample_id]

      hedges_g_results <- hedges_g_dge(
        meta_data = data.table::copy(sample_info_red),
        main_contrast = contrast_column,
        normalised_counts = norm_counts_red,
        contrast_list = contrast_list,
        .verbose = .verbose
      ) %>%
        .[, subgroup := group]

      return(hedges_g_results)
    })

    hedges_g_results_final <- data.table::rbindlist(results)
  }

  dge_params <- list(
    contrast_column = contrast_column,
    filter_column = filter_column
  )

  S7::prop(object, "outputs")[['hedges_g_res']] <- hedges_g_results_final
  S7::prop(object, "params")[["effect_size_dge"]] <- dge_params

  return(object)
}

# deprecated methods -----------------------------------------------------------

## qc --------------------------------------------------------------------------

#' QC on the bulk dge data (DEPRECATED!)
#'
#' @description
#' This is a deprecated method and will raise an error. Please use
#' [bixverse::qc_bulk_dge()] and [bixverse::normalise_bulk_dge()] instead.
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
#' @return Throws an error, as it is deprecated.
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
  stop(paste(
    "This function has been depracated.",
    "Please use qc_bulk_dge() and normalise_bulk_dge() instead."
  ))
}

## dge -------------------------------------------------------------------------

#' Calculate all possible DGE variants (DEPRECATED!)
#'
#' @description
#' This is a deprecated method and will raise an error. Please use
#' [bixverse::calculate_dge_limma()] and [bixverse::calculate_dge_hedges()].
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
#' @return Throws an error, as it is deprecated.
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
  stop(paste(
    "This function has been depracated.",
    "Please use calculate_dge_limma() and calculate_dge_hedges() instead."
  ))
}
