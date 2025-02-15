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
#' @description Get the final results from the
#'
#' @param bulk_coexp The underlying `bixverse_base_class` class, see
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
