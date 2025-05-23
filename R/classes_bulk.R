# classes ---------------------------------------------------------------------

#' @title Bulk RNAseq co-expression modules
#'
#' @description
#' Class for applying various co-expression module detection methods on top of
#' bulk RNAseq data.
#'
#' @param raw_data The raw count matrix. Rows = samples, columns = features.
#' @param meta_data data.table Metadata information on the samples. Expects to
#' have a `sample_id` column.
#' @param variable_info data.table. Metadata information on the features. This
#' is an optional table.
#'
#' @section Properties:
#' \describe{
#'   \item{raw_data}{A numerical matrix of the provided raw data.}
#'   \item{meta_data}{A data.table with the meta-information about the samples.}
#'   \item{variable_info}{An optional data.table containing the variable info.}
#'   \item{processed_data}{A list in which various types of processed data will
#'   be stored.}
#'   \item{outputs}{A list in which key outputs will be stored.}
#'   \item{params}{A (nested) list that will store all the parameters of the
#'   applied function.}
#'   \item{final_results}{A data.table that will contain the final results.}
#' }
#'
#' @return Returns the `bulk_coexp` class for further operations.
#'
#' @export
bulk_coexp <- S7::new_class(
  # Names, parents
  name = "bulk_coexp",
  parent = bixverse_base_class,

  # Properties
  properties = list(
    raw_data = S7::class_double,
    meta_data = S7::class_data.frame,
    variable_info = S7::class_any,
    processed_data = S7::class_list,
    outputs = S7::class_list,
    params = S7::class_list,
    final_results = S7::class_any
  ),
  constructor = function(raw_data, meta_data, variable_info = NULL) {
    # Checks
    checkmate::assertMatrix(raw_data, mode = "numeric")
    checkmate::assertDataTable(meta_data)
    checkmate::assertNames(
      names(meta_data),
      must.include = c("sample_id", "case_control")
    )
    checkmate::assert(
      checkmate::checkDataTable(variable_info),
      checkmate::checkNull(variable_info)
    )
    checkmate::assertTRUE(all(rownames(raw_data) %in% meta_data$sample_id))

    params <- list(
      "original_dim" = dim(raw_data)
    )

    S7::new_object(
      S7::S7_object(),
      raw_data = raw_data,
      meta_data = meta_data,
      variable_info = variable_info,
      processed_data = list(),
      outputs = list(),
      params = params,
      final_results = data.table::data.table()
    )
  }
)

#' @title Bulk RNAseq differential gene expression class
#'
#' @description
#' Class for coordinating differential gene expression analyses with subsequent
#' GSE in a structured format. Additionally, the class will store the counts in
#' [edgeR::DGEList()] for subsequent processing.
#'
#' @param raw_counts matrix. The raw count matrix. Rows = genes, columns =
#' samples. Note: this is different from the [bixverse::bulk_coexp()] class!
#' @param meta_data data.table. Metadata information on the samples. It expects
#' to have a column sample_id and case_control column.
#' @param variable_info data.table. Metadata information on the features. This
#' is an optional table. Defaults to `NULL`.
#' @param alternative_gene_id String. Optional alternative gene identifier to
#' be used. Must be a column of variable_info!
#'
#' @section Properties:
#' \describe{
#'   \item{raw_counts}{A numerical matrix of the provided raw data.}
#'   \item{meta_data}{A data.table with the meta-information about the samples.}
#'   \item{variable_info}{An optional data.table containing the variable info.}
#'   \item{outputs}{A list in which key outputs will be stored.}
#'   \item{plots}{A list with the plots that are generated during subsequent
#'   QC steps.}
#'   \item{params}{A (nested) list that will store all the parameters of the
#'   applied function.}
#'   \item{final_results}{A list in which final results will be stored.}
#' }
#'
#' @return Returns the `bulk_coexp` class for further operations.
#'
#' @export
bulk_dge <- S7::new_class(
  # Names, parents
  name = "bulk_dge",
  parent = bixverse_base_class,

  # Properties
  properties = list(
    raw_counts = S7::class_numeric,
    meta_data = S7::class_data.frame,
    variable_info = S7::class_any,
    outputs = S7::class_list,
    params = S7::class_list,
    plots = S7::class_list,
    final_results = S7::class_any
  ),
  constructor = function(
    raw_counts,
    meta_data,
    variable_info = NULL,
    alternative_gene_id = NULL
  ) {
    # Checks
    checkmate::assertMatrix(raw_counts, mode = "numeric")
    checkmate::assertDataTable(meta_data)
    checkmate::assertNames(
      names(meta_data),
      must.include = c("sample_id")
    )
    checkmate::assertTRUE(all(colnames(raw_counts) %in% meta_data$sample_id))
    checkmate::assert(
      checkmate::checkDataTable(variable_info),
      checkmate::checkNull(variable_info)
    )
    if (!is.null(variable_info)) {
      checkmate::checkNames(colnames(variable_info), must.include = "var_id")
    }
    if (!is.null(alternative_gene_id)) {
      checkmate::qassert(alternative_gene_id, "S")
      checkmate::assertDataTable(variable_info)
      checkmate::assertTRUE(alternative_gene_id %in% colnames(variable_info))
      rownames(raw_counts) <- variable_info[[alternative_gene_id]]
    }

    params <- list(
      original_dim = dim(raw_counts)
    )

    S7::new_object(
      S7::S7_object(),
      raw_counts = raw_counts,
      meta_data = meta_data,
      variable_info = variable_info,
      outputs = list(),
      plots = list(),
      params = params,
      final_results = list()
    )
  }
)

# additional constructors ------------------------------------------------------

## bulk_dge --------------------------------------------------------------------

#' Wrapper function to generate bulk_dge object from h5ad
#'
#' @description
#' This is a helper function that can be used to create a `bulk_dge` object
#' (see [bixverse::bulk_dge()]) directly from h5ad objects.
#'
#' @param h5_path String. Path to the h5ad object.
#' @param .verbose Controls verbosity of the function
#'
#' @returns `bulk_dge` object.
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
bulk_dge_from_h5ad <- function(
  h5_path,
  .verbose = TRUE
) {
  # Checks
  checkmate::qassert(h5_path, "S1")
  checkmate::assertFileExists(h5_path)

  h5_obj <- anndata_parser$new(h5_path)
  if (.verbose) message("Loading data from the h5ad object")
  c(meta_data, var_info, counts) %<-% h5_obj$get_key_data()
  bulk_dge_obj <- bulk_dge(
    raw_counts = counts,
    meta_data = meta_data,
    variable_info = var_info
  )
  return(bulk_dge_obj)
}

# utils ------------------------------------------------------------------------

## object manipulation ---------------------------------------------------------

#' Remove samples from object
#'
#' @description
#' This function allows to remove certain samples from the object
#'
#' @param object The underlying object, either `bixverse::bulk_coexp` or
#' `bixverse::bulk_dge`.
#' @param samples_to_remove Character vector. The sample identifiers to remove.
#' @param ... Additional arguments to parse to the functions.
#'
#' @return Returns the object with the samples removed. This will regenerated
#' the object from the start and remove any data in it.
#'
#' @export
remove_samples <- S7::new_generic(
  name = "remove_samples",
  dispatch_args = "object",
  fun = function(object, samples_to_remove, ...) {
    S7::S7_dispatch()
  }
)

#' @method remove_samples bulk_dge
#'
#' @export
S7::method(remove_samples, bulk_dge) <-
  function(object, samples_to_remove, ...) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::bulk_dge"
    )
    # Data
    meta_data <- S7::prop(object, "meta_data")
    raw_counts <- S7::prop(object, "raw_counts")
    variable_info <- S7::prop(object, "variable_info")

    meta_data_new <- meta_data[!sample_id %in% samples_to_remove]
    raw_counts_new <- raw_counts[, meta_data_new$sample_id]

    object_new <- bulk_dge(
      raw_counts = raw_counts_new,
      meta_data = meta_data_new,
      variable_info = variable_info
    )

    # Return
    return(object_new)
  }

## common getters --------------------------------------------------------------

#' Return the metadata
#'
#' @description
#' Getter function to extract the metadata from the [bixverse::bulk_coexp()] or
#' [bixverse::bulk_dge()].
#'
#' @param object The underlying object, either `bixverse::bulk_coexp` or
#' `bixverse::bulk_dge`.
#' @param ... Additional arguments to parse to the functions.
#'
#' @return Returns the metadata stored in the class.
#'
#' @export
get_metadata <- S7::new_generic(
  name = "get_metadata",
  dispatch_args = "object",
  fun = function(object, ...) {
    S7::S7_dispatch()
  }
)

#' @method get_metadata bulk_coexp
#'
#' @export
S7::method(get_metadata, bulk_coexp) <-
  function(object, ...) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::bulk_coexp"
    )

    # Return
    return(S7::prop(object, "meta_data"))
  }


#' @method get_metadata bulk_dge
#'
#' @export
S7::method(get_metadata, bulk_dge) <-
  function(object) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::bulk_dge"
    )

    # Return
    return(S7::prop(object, "meta_data"))
  }


#' Return the outputs
#'
#' @description
#' Getter function to extract the outputs from the [bixverse::bulk_coexp()] or
#' [bixverse::bulk_dge()].
#'
#' @param object The underlying object, either `bixverse::bulk_coexp` or
#' `bixverse::bulk_dge`.
#' @param ... Additional arguments to parse to the functions.
#'
#' @return Returns the outputs stored in the class.
#'
#' @export
get_outputs <- S7::new_generic(
  name = "get_outputs",
  dispatch_args = "object",
  fun = function(object, ...) {
    S7::S7_dispatch()
  }
)


#' @method get_outputs bulk_coexp
#'
#' @export
S7::method(get_outputs, bulk_coexp) <-
  function(object, ...) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::bulk_coexp"
    )

    # Return
    return(S7::prop(object, "outputs"))
  }


#' @method get_outputs bulk_dge
#'
#' @export
S7::method(get_outputs, bulk_dge) <-
  function(object, ...) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::bulk_dge"
    )

    # Return
    return(S7::prop(object, "outputs"))
  }

## individual getters ----------------------------------------------------------

### bulk dge class -------------------------------------------------------------

#' Return the DGEList
#'
#' @description
#' Getter function to extract the DGEList from the [bixverse::bulk_dge()] class.
#'
#' @param object `bulk_dge` class.
#'
#' @return Returns the DGEList stored in the class.
#'
#' @export
get_dge_list <- S7::new_generic(
  name = "get_dge_list",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)


#' @method get_dge_list bulk_dge
#'
#' @export
S7::method(get_dge_list, bulk_dge) <-
  function(object) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::bulk_dge"
    )

    # Return
    return(S7::prop(object, "outputs")[['dge_list']])
  }


#' Return the Limma Voom results
#'
#' @description
#' Getter function to extract the Limma Voom results from the
#' [bixverse::bulk_dge()] class.
#'
#' @param object `bulk_dge` class.
#'
#' @return Returns the Limma Voom results. (If found.)
#'
#' @export
get_dge_limma_voom <- S7::new_generic(
  name = "get_dge_limma_voom",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)


#' @method get_dge_limma_voom bulk_dge
#'
#' @export
S7::method(get_dge_limma_voom, bulk_dge) <-
  function(object) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::bulk_dge"
    )

    # Return
    return(S7::prop(object, "outputs")[['limma_voom_res']])
  }


#' Return the effect size results
#'
#' @description
#' Getter function to extract the Effect size results from the
#' [bixverse::bulk_dge()] class.
#'
#' @param object `bulk_dge` class.
#'
#' @return Returns the effect size results.  (If found.)
#'
#' @export
get_dge_effect_sizes <- S7::new_generic(
  name = "get_dge_effect_sizes",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)


#' @method get_dge_effect_sizes bulk_dge
#'
#' @export
S7::method(get_dge_effect_sizes, bulk_dge) <-
  function(object) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::bulk_dge"
    )

    # Return
    return(S7::prop(object, "outputs")[['hedges_g_res']])
  }

## individual setters ----------------------------------------------------------

### bulk dge class -------------------------------------------------------------

#' Change the primary gene identifier of bulk_dge
#'
#' @description
#' Changes the primary gene identifier in the bulk_dge class. To do so, you need
#' to either provide a `variable_info` data.table with the alternative gene
#' identifier you wish to use or it exists already in the object itself. If it
#' exists in the object, that variable_info will be used.
#'
#' @param object `bulk_dge` class, see [bixverse::bulk_dge].
#' @param alternative_gene_id String. The column containing the alternative gene
#' identifier. Must be present in the provided `variable_info` data.table or
#' within the class attributes.
#' @param variable_info Optional data.table with variable information. If
#' `variable_info` is in an attribute of the class, that one will be used.
#'
#' @return The class with modified primary gene identifier.
#'
#' @export
change_gene_identifier <- S7::new_generic(
  name = "change_gene_identifier",
  dispatch_args = "object",
  fun = function(object, alternative_gene_id, variable_info = NULL) {
    S7::S7_dispatch()
  }
)


#' @method change_gene_identifier bulk_dge
#'
#' @export
S7::method(change_gene_identifier, bulk_dge) <-
  function(object, alternative_gene_id, variable_info = NULL) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::bulk_dge"
    )

    if (is.null(S7::prop(object, "variable_info"))) {
      expected_nrows <- S7::prop(object, "params")[['original_dim']][1]
      checkmate::assertDataTable(variable_info, nrows = expected_nrows)
      S7::prop(object, "variable_info") <- variable_info
    }

    variable_info <- S7::prop(object, "variable_info")
    checkmate::assertTRUE(alternative_gene_id %in% colnames(variable_info))

    rownames(S7::prop(object, "raw_counts")) <- variable_info[[
      alternative_gene_id
    ]]

    # Return
    return(object)
  }


#' @method add_new_metadata bulk_dge
#'
#' @export
S7::method(add_new_metadata, bulk_dge) <-
  function(object, new_metadata, ...) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::bulk_dge"
    )
    checkmate::assertNames(
      names(new_metadata),
      must.include = c("sample_id")
    )

    raw_counts <- S7::prop(object, "raw_counts")

    checkmate::assertTRUE(all(colnames(raw_counts) %in% meta_data$sample_id))

    S7::prop(object, "meta_data") <- new_metadata

    return(object)
  }

## prints ----------------------------------------------------------------------

#' @name print.bulk_coexp
#' @title print Method for bulk_coexp object
#'
#' @description
#' Print a bulk_coexp object.
#'
#' @param x An object of class `bulk_coexp`.
#' @param ... Additional arguments (currently not used).
#'
#' @returns Invisibly returns `x`.
#'
#' @method print bulk_coexp
S7::method(print, bulk_coexp) <- function(x, ...) {
  # Pre-processing
  preprocessed <- !is.null(S7::prop(x, "processed_data")[["processed_data"]])
  features <- if (preprocessed) {
    no_features <- S7::prop(x, "processed_data")[["feature_meta"]] %>%
      .[(hvg)] %>%
      nrow()
    sprintf("  Number of HVG: %i.\n", no_features)
  } else {
    ""
  }
  # Prints for different methods
  method_info <- if (is.null(S7::prop(x, "params")[["detection_method"]])) {
    # Case of nothing has been applied
    ""
  } else if (S7::prop(x, "params")[["detection_method"]] == "cPCA") {
    # Contrastive PCA
    no_intersecting_features <- length(S7::prop(x, "params")[["cPCA_params"]][[
      "intersecting_features"
    ]])
    paste0(
      " Detection method: cPCA.\n",
      sprintf("  No of intersecting features: %i.\n", no_intersecting_features)
    )
  } else if (
    S7::prop(x, "params")[["detection_method"]] == "correlation-based"
  ) {
    # For simple correlations
    non_parametric <- S7::prop(x, "params")[["correlation_params"]][[
      "spearman"
    ]]
    graph_generated <- !is.null(S7::prop(x, "params")[["correlation_graph"]][[
      "no_nodes"
    ]])
    paste0(
      " Detection method: correlation based.\n",
      sprintf("  Non-parametric correlation applied: %s.\n", non_parametric),
      sprintf("  Graph generated: %s.\n", graph_generated)
    )
  }

  cat(
    "Bulk co-expression module class (bulk_coexp).\n",
    " Pre-processing done: ",
    preprocessed,
    ".\n",
    features,
    method_info,
    sep = ""
  )

  invisible(x)
}

# TODO write print for bulk_dge
