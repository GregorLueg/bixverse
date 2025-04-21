# methods ----------------------------------------------------------------------

#' Calculate PCA on the expression.
#'
#' @description
#' Calculates the principal component on top of the filtered count matrix and
#' adds the information of the first 10 principal components to the outputs.
#'
#' @param object The underlying class, see [bixverse::bulk_dge()].
#' @param scale Boolean. Shall the log(cpm) counts be scaled prior the PCA
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
    scale = FALSE,
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
  scale = FALSE,
  pcs = 10L,
  no_hvg_genes = 2500L
) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_dge")
  checkmate::qassert(scale, "B1")
  checkmate::qassert(pcs, "I1")
  checkmate::qassert(no_hvg_genes, "I1")

  dge_list <- S7::prop(object, "outputs")[['dge_list']]
  cpm <- edgeR::cpm(dge_list, log = TRUE)

  hvg_data <- data.table::setDT(
    list(
      gene_id = rownames(cpm),
      mad = matrixStats::rowMads(cpm)
    )
  ) %>%
    data.table::setorder(-mad)

  hvg_genes <- hvg_data[1:no_hvg_genes, gene_id]

  # Add the information to the table
  if (!is.null(S7::prop(object, "variable_info"))) {
    S7::prop(object, "variable_info")[, hvg := var_id %in% hvg_genes]
  }

  input_genes <- t(cpm[hvg_genes, ])
  pca_results <- rs_prcomp(input_genes, scale = scale)
  pcs_to_take <- min(pcs, ncol(pca_results$scores))
  pca_dt <- pca_results$scores[, 1:pcs_to_take] %>%
    `colnames<-`(sprintf("PC_%i", 1:pcs_to_take)) %>%
    data.table::as.data.table() %>%
    .[, sample_id := rownames(input_genes)] %>%
    .[, c("sample_id", sprintf("PC_%i", 1:pcs_to_take)), with = FALSE]

  pca_params <- list(
    hvg_genes = hvg_genes,
    scale = scale,
    pcs_taken = pcs_to_take
  )

  S7::prop(object, "params")[["pca"]] <-
    pca_params
  S7::prop(object, "outputs")[['pca_dt']] <- pca_dt

  return(object)
}


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
#' @param gene_filter Optional String vector. You can provide a vector of
#' strings for genes you wish to filter for, e.g., protein-coding genes.
#' @param co_variates Optional String vector. Any co-variates you wish to
#' consider during the Limma Voom modelling.
#' @param ... Additional parameters to forward to [edgeR::filterByExpr()].
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
    gene_filter = NULL,
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
  gene_filter = NULL,
  ...,
  .verbose = TRUE
) {
  # First checks
  checkmate::assertClass(object, "bixverse::bulk_dge")
  meta_data <- S7::prop(object, "meta_data")
  checkmate::qassert(contrast_column, "S+")
  checkmate::qassert(co_variates, c("S+", "0"))
  checkmate::qassert(filter_column, c("S+", "0"))
  checkmate::qassert(gene_filter, c("S+", "F+", "0"))
  all_specified_columns <- c(contrast_column, co_variates, filter_column)
  checkmate::assertTRUE(all(all_specified_columns %in% colnames(meta_data)))
  dge_list <- S7::prop(object, "outputs")[['dge_list']]

  if (!is.null(gene_filter)) {
    if (.verbose)
      message("Using the provided gene filter on the DGEList object")
    dge_list <- dge_list[gene_filter, ]
  }

  if (is.null(gene_filter)) {
    if (.verbose) message("Calculating the Limma Voom DGE results.")
    limma_results_final <- run_limma_voom(
      meta_info = meta_data,
      main_contrast = contrast_column,
      dge_list = dge_list,
      co_variates = co_variates,
      ...,
      .verbose = .verbose
    ) %>%
      .[, subgroup := NA]

    if (.verbose) message("Calculating the Hedge's G effect sizes.")
    hedges_g_results_final <- hedges_g_dge_list(
      meta_info = meta_data,
      main_contrast = contrast_column,
      dge_list = dge_list,
      ...,
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
    groups <- unique(meta_data[[filter_column]])
    results <- purrr::map(groups, \(group) {
      # Filter the meta data and dge list
      meta_data_red <- meta_data[
        eval(parse(
          text = paste0(filter_column, " == '", group, "'")
        ))
      ]
      dge_list_red <- dge_list[, meta_data_red$sample_id]

      # Limma Voom
      limma_results <- run_limma_voom(
        meta_info = data.table::copy(meta_data_red),
        main_contrast = contrast_column,
        dge_list = dge_list_red,
        co_variates = co_variates,
        ...,
        .verbose = .verbose
      ) %>%
        .[, subgroup := group]

      hedges_g_results <- hedges_g_dge_list(
        meta_info = data.table::copy(meta_data_red),
        main_contrast = contrast_column,
        dge_list = dge_list_red,
        ...,
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
    co_variates = co_variates,
    gene_filter = gene_filter
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
    pcs_to_plot = c("PC_1", "PC_2")
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
  pcs_to_plot = c("PC_1", "PC_2")
) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_dge")
  checkmate::qassert(cols_to_plot, "S2")
  checkmate::qassert(pcs_to_plot, "S+")

  pca_dt <- S7::prop(object, "outputs")[['pca_dt']]
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
