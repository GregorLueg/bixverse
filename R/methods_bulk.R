# methods ----------------------------------------------------------------------

## helpers ---------------------------------------------------------------------

#' Helper function to check if PC1/2 distinguish groups
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
.check_pca_grp_differences <- function(pc1, pc2, grps) {
  # Checks
  checkmate::qassert(pc1, "N+")
  checkmate::qassert(pc2, "N+")
  checkmate::qassert(grps, c("F+", "S+"))
  # Function body
  pca_anova_dt <- if (length(unique(grps)) > 1) {
    pc1_anova_pval <- anova(lm(pc1 ~ grps))$`Pr(>F)`[1]
    pc2_anova_pval <- anova(lm(pc2 ~ grps))$`Pr(>F)`[1]
    data.table(
      pc = c("PC1", "PC2"),
      pvalue = c(pc1_anova_pval, pc2_anova_pval)
    )
  } else {
    data.table(
      pc = c("PC1", "PC2"),
      pvalue = c(NA, NA)
    )
  }

  return(pca_anova_dt)
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

  p5_pca_case_control <- ggplot(
    data = plot_df,
    mapping = aes(x = PC_1, y = PC_2)
  ) +
    geom_point(mapping = aes(col = .data[[contrast_info]])) +
    xlab("PC1") +
    ylab("PC2") +
    theme_bw() +
    ggtitle("PCA with key columns") +
    labs(colour = "Groups:")

  pca_anova <- .check_pca_grp_differences(
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
S7::method(batch_correction_bulk_dge, bulk_dge) <- function(
  object,
  contrast_column,
  batch_col,
  scale_genes = FALSE,
  no_hvg_genes = 2500L
) {
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

  normalised_counts_corrected <- limma::removeBatchEffect(
    x = normalised_counts,
    batch = batch_data
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

  plot_uncor <- ggplot(data = pca_dt_uncor, mapping = aes(x = PC_1, y = PC_2)) +
    geom_point(mapping = aes(col = .data[[contrast_column]])) +
    theme_bw() +
    ggtitle("Pre batch correction") +
    xlab("PC1") +
    ylab("PC2")

  plot_cor <- ggplot(data = pca_dt_cor, mapping = aes(x = PC_1, y = PC_2)) +
    geom_point(mapping = aes(col = .data[[contrast_column]])) +
    theme_bw() +
    ggtitle("Post batch correction") +
    xlab("PC1") +
    ylab("PC2")

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
#' `contrast_column` that can be found in the meta-data. Every permutation of
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
#' @param co_variates Optional String vector. Any co-variates you wish to
#' consider during the Limma Voom modelling.
#' @param ... Additional parameters to forward to [limma::eBayes()] or
#' [limma::voom()].
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
    filter_column = NULL,
    co_variates = NULL,
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
  filter_column = NULL,
  co_variates = NULL,
  ...,
  .verbose = TRUE
) {
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
  all_specified_columns <- c(contrast_column, co_variates, filter_column)
  sample_info <- S7::prop(object, "outputs")[["sample_info"]]
  dge_list <- S7::prop(object, "outputs")[['dge_list']]
  checkmate::assertTRUE(all(all_specified_columns %in% colnames(sample_info)))

  counts_batch_cor <- S7::prop(object, "outputs")[[
    'normalised_counts_corrected'
  ]]

  norm_counts <- if (is.null(counts_batch_cor)) {
    S7::prop(object, "outputs")[['normalised_counts']]
  } else {
    if (.verbose)
      message(paste(
        "Found batch corrected counts.",
        "These will be used for effect size calculations"
      ))
    counts_batch_cor
  }

  if (is.null(filter_column)) {
    if (.verbose)
      message("Calculating the differential expression with limma results.")
    limma_results_final <- run_limma_voom(
      meta_data = sample_info,
      main_contrast = contrast_column,
      dge_list = dge_list,
      co_variates = co_variates,
      ...,
      .verbose = .verbose
    ) %>%
      .[, subgroup := NA]

    if (.verbose) message("Calculating the Hedge's G effect sizes.")
    hedges_g_results_final <- hedges_g_dge(
      meta_data = sample_info,
      main_contrast = contrast_column,
      normalised_counts = norm_counts,
      .verbose = .verbose
    ) %>%
      .[, subgroup := NA]
  } else {
    if (.verbose)
      message(paste(
        "Filtering column provided.",
        "Method will run Limma Voom and Hedge's G on the individual data sets."
      ))
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
      norm_counts_red <- norm_counts[, sample_info_red$sample_id]

      # Limma Voom
      limma_results <- run_limma_voom(
        meta_data = data.table::copy(sample_info_red),
        main_contrast = contrast_column,
        dge_list = dge_list_red,
        co_variates = co_variates,
        ...,
        .verbose = .verbose
      ) %>%
        .[, subgroup := group]

      hedges_g_results <- hedges_g_dge(
        meta_data = data.table::copy(sample_info_red),
        main_contrast = contrast_column,
        normalised_counts = norm_counts_red,
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
      rbindlist()

    hedges_g_results_final <- purrr::map(
      results,
      ~ {
        .[['hedges_g_results']]
      }
    ) %>%
      rbindlist()
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


# plotting ---------------------------------------------------------------------

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

  p <- ggplot(
    data = plot_df,
    mapping = aes(x = .data[[pcs_to_plot[1]]], y = .data[[pcs_to_plot[2]]])
  ) +
    geom_point(mapping = aes(col = value)) +
    facet_wrap(facets = ~variable, ncol = 3L) +
    xlab("PC1") +
    ylab("PC2") +
    theme_minimal() +
    ggtitle("PCA with key columns") +
    labs(colour = "Groups:")

  return(p)
}
