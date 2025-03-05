# class ----

bulk_coexp <- S7::new_class(
  # Names, parents
  name = "bulk_coexp",
  parent = bixverse_base_class,

  # Properties
  properties = list(
    raw_data = S7::class_double,
    meta_data = S7::class_data.frame,
    processed_data = S7::class_list,
    outputs = S7::class_list,
    params = S7::class_list,
    final_results = S7::class_data.frame
  ),

  #' Bulk RNAseq co-expression modules
  #'
  #' @description
  #' Class for applying various co-expression module detection methods on top of
  #' bulk RNAseq data.
  #'
  #' @param raw_data The raw count matrix. Rows = samples, columns = features.
  #' @param meta_data Metadata information on the samples. It expects to have a
  #' column sample_id and case_control column.
  #'
  #' @return Returns the `bulk_coexp` class for further operations.
  #'
  #' @export
  constructor = function(raw_data, meta_data) {
    # Checks
    checkmate::assertMatrix(raw_data, mode = "numeric")
    checkmate::assertDataFrame(meta_data)
    checkmate::assertNames(names(meta_data),
                           must.include = c("sample_id", "case_control"))
    checkmate::assertTRUE(all(rownames(raw_data) %in% meta_data$sample_id))

    params <- list(
      "original_dim" = dim(raw_data)
    )

    S7::new_object(
      S7::S7_object(),
      raw_data = raw_data,
      meta_data = data.table::as.data.table(meta_data),
      processed_data = list(),
      outputs = list(),
      params = params,
      final_results = data.table::data.table()
    )
  }
)

# getters ----

#' Return the outputs from bulk_coexp
#'
#' @description
#' Getter function to extract the outputs from the [bixverse::bulk_coexp()]
#' class.
#'
#' @export
get_outputs <- S7::new_generic(
  "get_outputs",
  "bulk_coexp"
)


#' @name get_outputs
#'
#' @description Get the outputs that are stored in the class.
#'
#' @param bulk_coexp The underlying `bulk_coexp` class, see
#' [bixverse::bulk_coexp()].
#'
#' @return Returns the outputs stored in the class.
#'
#' @method get_outputs bulk_coexp
S7::method(get_outputs, bulk_coexp) <-
  function(bulk_coexp) {
    # Checks
    checkmate::assertClass(
      bulk_coexp, "bixverse::bulk_coexp"
    )

    # Return
    return(S7::prop(bulk_coexp, "outputs"))
  }


# prints ----

# S7::method(format, bulk_coexp) <- function(x, ...) {
#   pre_processed <- !purrr::is_empty(S7::prop(bulk_coexp, "processed_data"))
#   co_exp_method <- if (purrr::is_empty(S7::prop(bulk_coexp, "params")[["detection_method"]])) {
#     'not defined yet'
#   } else {
#     S7::prop(bulk_coexp, "params")[["detection_method"]]
#   }
#   outputs_available = !purrr::is_empty(S7::prop(bulk_coexp, "params")[["outputs"]])
#
#   sprintf(
#     "Bulk co-expression class object:\n  Pre-processsed: %b\n  Method: %s\n  Outputs available: %b",
#     pre_processed,
#     co_exp_method,
#     outputs_available
#   )
# }

# general methods ----

#' Process the raw data
#'
#' @description
#' Function to do general pre-processing on top of the [bixverse::bulk_coexp()].
#' Options to do scaling, HVG selection, etc.
#'
#' @export
preprocess_bulk_coexp <- S7::new_generic(
  "preprocess_bulk_coexp",
  "bulk_coexp"
)


#' @name preprocess_bulk_coexp
#'
#' @description Do pre-processing on top of the data. This can include selection
#' of highly variable genes and/or scaling.
#'
#' @param bulk_coexp The underlying `bixverse_base_class` class, see
#' [bixverse::bulk_coexp()].
#' @param hvg ...
#' @param mad_threshold ...
#' @param scaling ...
#' @param scaling_type ...
#'
#' @return Returns the class with the processed_data data slot populated.
#'
#' @method preprocess_bulk_coexp bulk_coexp
S7::method(preprocess_bulk_coexp, bulk_coexp) <-
  function(bulk_coexp,
           hvg = NULL,
           mad_threshold = NULL,
           scaling = F,
           scaling_type = c("normal", "robust")) {
    # Checks
    checkmate::assertClass(bulk_coexp, "bixverse::bulk_coexp")
    checkmate::qassert(mad_threshold, c("R1", "0"))
    nfeatures <- S7::prop(bulk_coexp, "params")[['original_dim']][2]
    checkmate::qassert(hvg, c("R1[0,1]", sprintf("I1[0,%i]", nfeatures), "0"))
    checkmate::qassert(scaling, "B1")
    if (scaling) {
      checkmate::assertChoice(scaling_type, c("normal", "robust"))
    }

    mat <- S7::prop(bulk_coexp, "raw_data")

    feature_meta <- data.table::data.table(
      feature_name = colnames(mat),
      mean_exp = colMeans(mat),
      MAD = matrixStats::colMads(mat),
      var_exp = matrixStats::colVars(mat)
    ) %>%
      data.table::setorder(-MAD)

    if(is.null(hvg) & is.null(mad_threshold)) {
      hvg = 1
    }

    if (!is.null(hvg)) {
      no_genes_to_take <-
        ifelse(is.integer(hvg), hvg, ceiling(hvg * ncol(mat)))
      hvg_genes <- feature_meta[1:no_genes_to_take, feature_name]
    } else {
      hvg_genes <- feature_meta[MAD >= mad_threshold, feature_name]
    }

    feature_meta[, hvg := feature_name %in% hvg_genes]

    # Process the matrix
    matrix_processed <- mat[, hvg_genes]

    if (scaling) {
      fun <-
        ifelse(scaling_type == "normal",
               "scale",
               "bixverse::robust_scale"
        )
      matrix_processed <- rlang::eval_tidy(rlang::quo(apply(
        matrix_processed, 1, !!!rlang::parse_exprs(fun)
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

    S7::prop(bulk_coexp, "params")[['preprocessing']] <- processing_params
    S7::prop(bulk_coexp, "processed_data")[['processed_data']] <- matrix_processed
    S7::prop(bulk_coexp, "processed_data")[['feature_meta']] <- feature_meta

    # Return
    bulk_coexp
  }

