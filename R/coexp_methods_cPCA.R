# processing ----

#' Prepare class for contrastive PCA
#'
#' @description
#' This is the generic function for doing the preprocessing for contrastive PCA.
#'
#' @export
contrastive_pca_processing <- S7::new_generic(
  "contrastive_pca_processing",
  "bulk_coexp"
)


#' @name contrastive_pca_processing
#'
#' @description
#' This function will prepare the `bulk_coexp` for subsequent usage of the
#' contrastive PCA function.
#'
#' @usage contrastive_pca_processing(
#'  bulk_coexp = bulk_coexp,
#'  background_matrix = background_matrix,
#'  scale = FALSE
#'  .verbose = TRUE
#' )
#'
#' @param bulk_coexp `bulk_coexp` class.
#' @param background_matrix Numeric matrix.
#' @param scale Boolean. Shall the data be scaled. Defaults to FALSE.
#' @param verbose Boolean. Controls verbosity of the function.
#'
#' @return `bulk_coexp` with additional data in the slots
#'
#' @export
#'
#' @references Abid, et al., Nature Communications, 2018
#'
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
#' @import data.table
#'
#' @method contrastive_pca_processing bulk_coexp
S7::method(contrastive_pca_processing, bulk_coexp) <-
  function(bulk_coexp,
           background_matrix,
           scale = FALSE,
           .verbose = TRUE) {
    # Checks
    checkmate::assertClass(bulk_coexp, "bixverse::bulk_coexp")
    checkmate::assertMatrix(background_matrix, mode = "numeric")
    checkmate::qassert(scale, "B1")
    checkmate::qassert(.verbose, "B1")

    # Function body
    target_mat <- S7::prop(bulk_coexp, "raw_data")
    background_mat <- background_mat

    intersecting_features <- intersect(
      colnames(target_mat),
      colnames(background_mat)
    )
    if (.verbose)
      message(sprintf(
        "A total of %i features/genes were identified",
        length(intersecting_features)
      )
      )

    target_mat <- target_mat[, intersecting_features]
    background_mat <- background_mat[, intersecting_features]

    if (scale) {
      target_mat <- scale(target_mat, scale = scale)
      background_mat <- scale(background_mat, scale = scale)
    }

    target_covar <- rs_covariance(target_mat)
    background_covar <- rs_covariance(background_mat)

    internal_params <- list(
      "intersecting_features" = intersecting_features,
      'scaled_data' = scale
    )

    # Data
    S7::prop(bulk_coexp, "processed_data")[["target_mat"]] <-
      target_mat
    S7::prop(bulk_coexp, "processed_data")[["background_mat"]] <-
      background_mat
    # Covariance matrices
    S7::prop(bulk_coexp, "processed_data")[["target_covar"]] <-
      target_covar
    S7::prop(bulk_coexp, "processed_data")[["background_covar"]] <-
      background_covar


    # Set the object to a cPCA analysis
    S7::prop(bulk_coexp, "params")["detection_method"] <- "cPCA"
    S7::prop(bulk_coexp, "params")[["cPCA_params"]] <- internal_params


    # Return
    bulk_coexp
  }


# methods ----



# plotting ----

#' Plot various alphas for the contrastive PCA
#'
#' @description
#' This function will plot
#'
#' @export
c_pca_plot_alphas <- S7::new_generic("c_pca_plot_alphas", "bulk_coexp")

