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
#' @param .verbose Boolean.
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
    no_hvg_genes = 2500L,
    .verbose = TRUE
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
  no_hvg_genes = 2500L,
  .verbose = TRUE
) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_dge")
  checkmate::qassert(scale, "B1")
  checkmate::qassert(pcs, "I1")
  checkmate::qassert(no_hvg_genes, "I1")
  checkmate::qassert(.verbose, "B1")

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

# helpers ----------------------------------------------------------------------
