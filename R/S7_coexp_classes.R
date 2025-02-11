# class ----

bulk_coexp <- S7::new_class(
  # Names, parents
  name = "bulk_coexp",
  parent = bixverse_generic_class,

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
  #' ...
  #'
  #' @param raw_data ...
  #' @param meta_data ...
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
    checkmate::assertTRUE(all(colnames(raw_data) %in% meta_data$sampleID))

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
