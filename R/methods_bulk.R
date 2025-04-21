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
#' @param hvg_genes Integer. Number of highly variable genes to include.
#'
#' @return Returns the class with additional data added to the outputs.
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
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
#' @export
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
#'
#' @import data.table
#' @import ggplot2
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
plot_pca_res <- S7::new_generic(
  "plot_pca_res",
  "object",
  fun = function(
    object,
    cols_to_plot = c('contrast_info', 'sample_source'),
    PCs_to_plot = c("PC_1", "PC_2")
  ) {
    S7::S7_dispatch()
  }
)

#' @method plot_pca_res bulk_dge
#' @export
S7::method(plot_pca_res, bulk_dge) <- function(
  object,
  cols_to_plot = c('contrast_info', 'sample_source'),
  PCs_to_plot = c("PC_1", "PC_2")
) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_dge")
  checkmate::qassert(cols_to_plot, "S2")
  checkmate::qassert(PCs_to_plot, "S+")

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
    must.include = PCs_to_plot
  )
  checkmate::assertNames(
    names(meta_data),
    must.include = cols_to_plot
  )

  plot_df <- data.table::merge.data.table(
    pca_dt[, c('sample_id', PCs_to_plot), with = FALSE],
    meta_data[, c('sample_id', cols_to_plot), with = FALSE],
    by.x = 'sample_id',
    by.y = 'sample_id'
  ) %>%
    data.table::melt(id.vars = c('sample_id', PCs_to_plot))

  p <- ggplot(
    data = plot_df,
    mapping = aes(x = .data[[PCs_to_plot[1]]], y = .data[[PCs_to_plot[2]]])
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

# helpers ----------------------------------------------------------------------

## general helpers -------------------------------------------------------------

#' Fixes contrast names for DGEs
#'
#' @param x Vector of strings or factors.
#'
#' @returns Vector with fixed naming based on R conventions.
#'
#' @export
fix_contrast_names <- function(x) {
  checkmate::qassert(x, c("S+", "F+"))
  as.factor(gsub("\\.", "_", make.names(gsub("[[:punct:]]", "", x))))
}


#' Create all limma contrasts
#'
#' @param limma_fit The fitted limma model, i.e., output of [limma::lmFit()].
#'
#' @returns The Limma contrasts for further usage.
#'
#' @export
all_contrasts <- function(limma_fit) {
  checkmate::assertClass(limma_fit, "MArrayLM")
  coefs_fit <- colnames(coef(limma_fit))
  cb <- combn(coefs_fit, 2, FUN = function(x) {
    paste0(x[[1]], "-", x[[2]])
  })
  contrasts <- limma::makeContrasts(contrasts = cb, levels = coefs_fit)
  return(contrasts)
}

## dge helpers -----------------------------------------------------------------

#' Run Limma Voom
#'
#' @description
#' Wrapper function to run Limma Voom workflows.
#'
#'
#' @param meta_info data.table. The meta information.
#' @param main_contrast String. Which column contains the main groups you want
#' to test dge for.
#' @param dge_list DGEList, see [edgeR::DGEList()].
#' @param co_variates String or NULL. Optional co-variates you wish to consider
#' during model fitting.
#' @param ... Additional parameters to forward to [edgeR::filterByExpr()].
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A data.table with all the DGE results by the identified contrasts.
#'
#' @export
run_limma_voom <- function(
  meta_info,
  main_contrast,
  dge_list,
  co_variates = NULL,
  ...,
  .verbose = TRUE
) {
  variabales <- c(main_contrast, co_variates)
  # Checks
  checkmate::assertDataFrame(meta_info)
  checkmate::qassert(main_contrast, "S1")
  checkmate::assertClass(dge_list, "DGEList")
  checkmate::qassert(co_variates, c("S+", "0"))
  checkmate::assertNames(
    names(meta_info),
    must.include = variabales
  )
  checkmate::qassert(.verbose, "B1")

  # Fix any names
  meta_info[,
    (variables) := lapply(.SD, fix_contrast_names),
    .SDcols = variables
  ]
  if (.verbose)
    message(paste(
      "Fixing any naming issues for the selected main contrast",
      "and any co-variates."
    ))

  model_matrix <- model.matrix(as.formula(model_formula), data = meta_data_red)
  colnames(model_matrix) <- gsub(main_contrast, "", colnames(model_matrix))

  # Filter lowly expressed genes
  to_keep <- edgeR::filterByExpr(
    y = dge_list,
    design = model_matrix,
    ...
  )
  if (.verbose) message(sprintf("A total of %i genes are kept.", sum(to_keep)))

  voom_obj <- limma::voom(
    counts = dge_list[to_keep, ],
    design = model_matrix,
    normalize.method = "quantile",
    plot = TRUE
  )
  limma_fit <- limma::lmFit(voom_obj, model_matrix)

  contrasts <- all_contrasts(limma_fit)

  final_fit <- limma::contrasts.fit(limma_fit, contrasts)
  final_fit <- limma::eBayes(final_fit)

  tested_contrasts <- attributes(contrasts)$dimnames$Contrasts

  all_dge_res <- purrr::map(tested_contrasts, \(coef) {
    top.table <- as.data.table(
      limma::topTable(fit = fit2, coef = coef, sort.by = "P", n = Inf),
      keep.rownames = TRUE
    ) %>%
      setnames(old = 'rn', new = 'gene_id') %>%
      .[, contrast := gsub("-", "_vs_", coef)]
  }) %>%
    rbindlist()

  return(all_dge_res)
}
