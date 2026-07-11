# classes ----------------------------------------------------------------------

## BulkCoExp -------------------------------------------------------------------

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
#' @return Returns the `BulkCoExp` class for further operations.
#'
#' @export
BulkCoExp <- S7::new_class(
  # Names, parents
  name = "BulkCoExp",
  parent = BixverseBaseClass,

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
    checkmate::assertMatrix(
      raw_data,
      mode = "numeric",
      row.names = "named",
      col.names = "named"
    )
    checkmate::assertDataTable(meta_data)
    checkmate::assertNames(
      names(meta_data),
      must.include = c("sample_id")
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
#' samples. Note: this is different from the [bixverse::BulkCoExp()] class!
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
#' @return Returns the `BulkDge` class for further operations.
#'
#' @export
BulkDge <- S7::new_class(
  # Names, parents
  name = "BulkDge",
  parent = BixverseBaseClass,

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
    checkmate::assertMatrix(
      raw_counts,
      mode = "numeric",
      row.names = "named",
      col.names = "named"
    )
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

## BulkDge ---------------------------------------------------------------------

#' Wrapper function to generate BulkDge object from h5ad
#'
#' @description
#' This is a helper function that can be used to create a `BulkDge` object
#' (see [bixverse::BulkDge()]) directly from h5ad objects.
#'
#' @param h5_path String. Path to the h5ad object.
#' @param .verbose Controls verbosity of the function
#'
#' @returns `BulkDge` object.
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
bulk_dge_from_h5ad <- function(
  h5_path,
  .verbose = TRUE
) {
  # Globals
  meta_data <- var_info <- counts <- NULL

  # Checks
  checkmate::qassert(h5_path, "S1")
  checkmate::assertFileExists(h5_path)

  h5_obj <- anndata_parser$new(h5_path)
  if (.verbose) {
    message("Loading data from the h5ad object")
  }
  c(meta_data, var_info, counts) %<-% h5_obj$get_key_data()
  bulk_dge_obj <- BulkDge(
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
#' @param object The underlying object, either [bixverse::BulkCoExp()] or
#' [bixverse::BulkDge()].
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

#' @method remove_samples BulkDge
#'
#' @export
S7::method(remove_samples, BulkDge) <-
  function(object, samples_to_remove, ...) {
    sample_id <- NULL
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::BulkDge"
    )
    # Data
    meta_data <- S7::prop(object, "meta_data")
    raw_counts <- S7::prop(object, "raw_counts")
    variable_info <- S7::prop(object, "variable_info")

    meta_data_new <- meta_data[!sample_id %in% samples_to_remove]
    raw_counts_new <- raw_counts[, meta_data_new$sample_id]

    object_new <- BulkDge(
      raw_counts = raw_counts_new,
      meta_data = meta_data_new,
      variable_info = variable_info
    )

    # Return
    return(object_new)
  }


#' Helper to fix meta-data columns to be R conform
#'
#' @description
#' This function will update the specified columns in the metadata of an
#' [bixverse::BulkDge()] or [bixverse::BulkCoExp()] to be conform with R standard
#' naming conventions. This is useful to do before running DGE methods as they
#' expect standardised names.
#'
#' @param object The underlying object, either [bixverse::BulkCoExp()] or
#' [bixverse::BulkDge()].
#' @param col_names Character vector. The columns to fix.
#' @param ... Additional arguments to parse to the functions.
#'
#' @return Returns the object with the respective metadata columns updated.
#'
#' @export
fix_meta_data_column <- S7::new_generic(
  name = "fix_meta_data_column",
  dispatch_args = "object",
  fun = function(object, col_names, ...) {
    S7::S7_dispatch()
  }
)

#' @method fix_meta_data_column BulkDge
#'
#' @export
S7::method(fix_meta_data_column, BulkDge) <-
  function(object, col_names) {
    checkmate::assertClass(
      object,
      "bixverse::BulkDge"
    )
    # Data
    S7::prop(object, "meta_data") <- S7::prop(object, "meta_data")[,
      (col_names) := lapply(.SD, fix_contrast_names),
      .SDcols = col_names
    ]

    # Return
    return(object)
  }


#' Replace values in a metadata column
#'
#' @description
#' This function will update the values in a given metadata column based on
#' what you are providing in terms of replacement.
#'
#' @param object The underlying object, either [bixverse::BulkCoExp()] or
#' [bixverse::BulkDge()].
#' @param column Character vector. The columns for which to replace the
#' values.
#' @param replacement Named character vector. The values with which to replace
#' the data.
#' @param ... Additional arguments to parse to the functions.
#'
#' @return Returns the object with the respective metadata updated.
#'
#' @export
update_metadata_values <- S7::new_generic(
  name = "update_metadata_values",
  dispatch_args = "object",
  fun = function(object, column, replacement, ...) {
    S7::S7_dispatch()
  }
)

#' @method update_metadata_values BulkDge
#'
#' @export
S7::method(update_metadata_values, BulkDge) <-
  function(object, column, replacement) {
    meta_data <- S7::prop(object, "meta_data")

    checkmate::assertClass(
      object,
      "bixverse::BulkDge"
    )
    checkmate::assertTRUE(
      all(meta_data[[column]] %in% names(replacement))
    )

    meta_data[[column]] <- replacement[meta_data[[column]]]

    S7::prop(object, "meta_data") <- meta_data

    # Return
    return(object)
  }

## common getters --------------------------------------------------------------

#' Return the metadata
#'
#' @description
#' Getter function to extract the metadata from the [bixverse::BulkCoExp()] or
#' [bixverse::BulkDge()].
#'
#' @param object The underlying object, either [bixverse::BulkCoExp()] or
#' [bixverse::BulkDge()].
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

#' @method get_metadata BulkCoExp
#'
#' @export
S7::method(get_metadata, BulkCoExp) <-
  function(object, ...) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::BulkCoExp"
    )

    # Return
    return(S7::prop(object, "meta_data"))
  }


#' @method get_metadata BulkDge
#'
#' @export
S7::method(get_metadata, BulkDge) <-
  function(object) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::BulkDge"
    )

    # Return
    return(S7::prop(object, "meta_data"))
  }


#' Return the outputs
#'
#' @description
#' Getter function to extract the outputs from the [bixverse::BulkCoExp()] or
#' [bixverse::BulkDge()].
#'
#' @param object The underlying object, either [bixverse::BulkCoExp()] or
#' [bixverse::BulkDge()].
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


#' @method get_outputs BulkCoExp
#'
#' @export
S7::method(get_outputs, BulkCoExp) <-
  function(object, ...) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::BulkCoExp"
    )

    # Return
    return(S7::prop(object, "outputs"))
  }


#' @method get_outputs BulkDge
#'
#' @export
S7::method(get_outputs, BulkDge) <-
  function(object, ...) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::BulkDge"
    )

    # Return
    return(S7::prop(object, "outputs"))
  }

## individual getters ----------------------------------------------------------

### bulk dge class -------------------------------------------------------------

#' Return the DGEList
#'
#' @description
#' Getter function to extract the DGEList from the [bixverse::BulkDge()] class.
#'
#' @param object `BulkDge` class.
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


#' @method get_dge_list BulkDge
#'
#' @export
S7::method(get_dge_list, BulkDge) <-
  function(object) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::BulkDge"
    )

    # Return
    return(S7::prop(object, "outputs")[['dge_list']])
  }


#' Return the Limma Voom results
#'
#' @description
#' Getter function to extract the Limma Voom results from the
#' [bixverse::BulkDge()] class.
#'
#' @param object `BulkDge` class.
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


#' @method get_dge_limma_voom BulkDge
#'
#' @export
S7::method(get_dge_limma_voom, BulkDge) <-
  function(object) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::BulkDge"
    )

    limma_res <- S7::prop(object, "outputs")[['limma_voom_res']]
    if (is.null(limma_res)) {
      warning(paste(
        "No results found. Did you run get_dge_limma_voom()?",
        "Returning NULL"
      ))
    }

    # Return
    return(limma_res)
  }


#' Return the effect size results
#'
#' @description
#' Getter function to extract the Effect size results from the
#' [bixverse::BulkDge()] class.
#'
#' @param object `BulkDge` class.
#'
#' @return Returns the effect size results. (If found.)
#'
#' @export
get_dge_effect_sizes <- S7::new_generic(
  name = "get_dge_effect_sizes",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @method get_dge_effect_sizes BulkDge
#'
#' @export
S7::method(get_dge_effect_sizes, BulkDge) <-
  function(object) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::BulkDge"
    )

    hedges_g_res <- S7::prop(object, "outputs")[['hedges_g_res']]
    if (is.null(hedges_g_res)) {
      warning(paste(
        "No results found. Did you run calculate_dge_hedges()?",
        "Returning NULL"
      ))
    }

    # Return
    return(hedges_g_res)
  }

#' Return the TPM-normalised counts
#'
#' @description
#' Getter function to extract the TPM-normalised counts from the
#' [bixverse::BulkDge()] class.
#'
#' @param object `BulkDge` class.
#'
#' @return Returns the TPM-normalised counts. (If found.)
#'
#' @export
get_tpm_counts <- S7::new_generic(
  name = "get_tpm_counts",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @method get_tpm_counts BulkDge
#'
#' @export
S7::method(get_tpm_counts, BulkDge) <-
  function(object) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::BulkDge"
    )

    tpm_counts <- S7::prop(object, "outputs")[['tpm_counts']]
    if (is.null(tpm_counts)) {
      warning(paste(
        "No TPM counts found. Did you run normalise_bulk_dge()",
        "with calc_tpm = TRUE? Returning NULL"
      ))
    }

    # Return
    return(tpm_counts)
  }

#' Return the FPKM-normalised counts
#'
#' @description
#' Getter function to extract the FPKM-normalised counts from the
#' [bixverse::BulkDge()] class.
#'
#' @param object `BulkDge` class.
#'
#' @return Returns the FPKM-normalised counts. (If found.)
#'
#' @export
get_fpkm_counts <- S7::new_generic(
  name = "get_fpkm_counts",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @method get_fpkm_counts BulkDge
#'
#' @export
S7::method(get_fpkm_counts, BulkDge) <-
  function(object) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::BulkDge"
    )

    fpkm_counts <- S7::prop(object, "outputs")[['fpkm_counts']]
    if (is.null(fpkm_counts)) {
      warning(paste(
        "No FPKM counts found. Did you run normalise_bulk_dge()",
        "with calc_fpkm = TRUE? Returning NULL"
      ))
    }

    # Return
    return(fpkm_counts)
  }

### bulk coexp class -----------------------------------------------------------

#' Return the epsilon data
#'
#' @description
#' Getter function to extract the `epsilon param ~ power law goodness of fit`
#' data from the [bixverse::BulkCoExp()] class.
#'
#' @param object `BulkCoExp` class.
#'
#' @return Returns the epsilon data. (If found. Otherwise `NULL`).
#'
#' @export
get_epsilon_res <- S7::new_generic(
  name = "get_epsilon_res",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @method get_epsilon_res BulkCoExp
#'
#' @export
S7::method(get_epsilon_res, BulkCoExp) <-
  function(object) {
    # checks
    checkmate::assertClass(
      object,
      "bixverse::BulkCoExp"
    )
    # body
    epsilon_results <- S7::prop(object, "outputs")[['epsilon_data']]
    if (is.null(epsilon_results)) {
      warning(
        paste(
          "No epsilon results found.",
          "Did you run cor_module_check_epsilon()? Returning NULL."
        )
      )
    }

    # return
    return(epsilon_results)
  }

#' @title Return the resolution results
#'
#' @description
#' Getter function to get the resolution results (if available).
#'
#' @param object The class, see [bixverse::BulkCoExp()].
#'
#' @return If resolution results were found, returns the data.table. Otherwise,
#' throws a warning and returns NULL.
#'
#' @export
get_resolution_res <- S7::new_generic(
  name = "get_resolution_res",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @export
#' @method get_resolution_res BulkCoExp
S7::method(get_resolution_res, BulkCoExp) <- function(object) {
  # Checks
  checkmate::assertClass(object, "bixverse::BulkCoExp")
  # Body
  resolution_results <- S7::prop(object, "outputs")[["resolution_results"]]
  if (is.null(resolution_results)) {
    warning(
      paste(
        "No resolution results found.",
        "Did you run cor_module_graph_check_res()? Returning NULL."
      )
    )
  }

  resolution_results
}

## individual setters ----------------------------------------------------------

### bulk dge class -------------------------------------------------------------

#' Change the primary gene identifier of BulkDge
#'
#' @description
#' Changes the primary gene identifier in the BulkDge class. To do so, you need
#' to either provide a `variable_info` data.table with the alternative gene
#' identifier you wish to use or it exists already in the object itself. If it
#' exists in the object, that variable_info will be used.
#'
#' @param object `BulkDge` class, see [bixverse::BulkDge()].
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


#' @method change_gene_identifier BulkDge
#'
#' @export
S7::method(change_gene_identifier, BulkDge) <- function(
  object,
  alternative_gene_id,
  variable_info = NULL
) {
  # Checks
  checkmate::assertClass(
    object,
    "bixverse::BulkDge"
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


#' @method add_new_metadata BulkDge
#'
#' @export
S7::method(add_new_metadata, BulkDge) <-
  function(object, new_metadata, ...) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::BulkDge"
    )
    checkmate::assertNames(
      names(new_metadata),
      must.include = c("sample_id")
    )

    raw_counts <- S7::prop(object, "raw_counts")

    checkmate::assertTRUE(all(colnames(raw_counts) %in% new_metadata$sample_id))

    S7::prop(object, "meta_data") <- new_metadata

    return(object)
  }

## prints ----------------------------------------------------------------------

#' @noRd
S7::method(print, BulkCoExp) <- function(x, ...) {
  . <- hvg <- NULL

  params <- S7::prop(x, "params")
  outputs <- S7::prop(x, "outputs")
  final_results <- S7::prop(x, "final_results")

  n_modules_from_result <- if (
    inherits(final_results, "BulkModuleResult")
  ) {
    length(unique(get_modules(final_results)$module_id))
  } else {
    NA_integer_
  }

  # Pre-processing summary
  preprocessed <- !is.null(S7::prop(x, "processed_data")[["processed_data"]])
  hvg_line <- if (preprocessed) {
    no_features <- S7::prop(x, "processed_data")[["feature_meta"]] %>%
      .[(hvg)] %>%
      nrow()
    sprintf("  Number of HVG: %i.\n", no_features)
  } else {
    ""
  }

  detection_method <- params[["detection_method"]]

  method_info <- if (is.null(detection_method)) {
    ""
  } else if (detection_method == "c_pca") {
    # Contrastive PCA
    no_intersecting_features <- length(params[["c_pca_params"]][[
      "intersecting_features"
    ]])
    paste0(
      " Detection method: contrastive PCA.\n",
      sprintf("  No of intersecting features: %i.\n", no_intersecting_features)
    )
  } else if (detection_method == "correlation-based") {
    # Simple correlations + optional TOM + optional graph + optional CoReMo
    non_parametric <- params[["correlation_params"]][["spearman"]]
    tom_applied <- isTRUE(params[["correlation_params"]][["TOM"]])
    graph_generated <- !is.null(params[["correlation_graph"]][["no_nodes"]])
    coremo_run <- !is.null(outputs[["final_modules"]])
    leiden_run <- !is.null(params[["module_final_gen"]])

    lines <- c(
      " Detection method: correlation-based.\n",
      sprintf("  Spearman correlation: %s.\n", non_parametric),
      sprintf("  TOM applied: %s.\n", tom_applied),
      sprintf("  Correlation graph built: %s.\n", graph_generated)
    )
    if (coremo_run) {
      no_modules <- length(unique(outputs[["final_modules"]]$cluster_id))
      lines <- c(
        lines,
        sprintf("  CoReMo modules identified: %i.\n", no_modules)
      )
    }
    if (leiden_run) {
      lines <- c(
        lines,
        sprintf(
          "  Leiden modules identified: %i.\n",
          n_modules_from_result
        )
      )
    }
    paste(lines, collapse = "")
  } else if (detection_method == "differential correlation-based") {
    non_parametric <- params[["correlation_params"]][["spearman"]]
    n_shared <- params[["correlation_params"]][["no_intersecting_features"]]
    graph_generated <- !is.null(params[["correlation_graph"]][["no_nodes"]])
    modules_found <- !is.na(n_modules_from_result)
    lines <- c(
      " Detection method: differential correlation-based.\n",
      sprintf("  Spearman correlation: %s.\n", non_parametric),
      sprintf("  Shared features between target and background: %s.\n", n_shared),
      sprintf("  Differential correlation graph built: %s.\n", graph_generated)
    )
    if (modules_found) {
      lines <- c(
        lines,
        sprintf(
          "  Differential modules identified: %i.\n",
          n_modules_from_result
        )
      )
    }
    paste(lines, collapse = "")
  } else if (detection_method == "ICA-based") {
    ica_assessment <- params[["ica_stability_assessment"]]
    ica_final <- params[["ica_final_gen"]]
    tested_ncomps <- ica_assessment[["tested_ncomps"]]
    optimal_ncomp <- ica_assessment[["optimal_ncomp"]]
    lines <- c(
      " Detection method: ICA-based.\n"
    )
    if (!is.null(tested_ncomps)) {
      lines <- c(
        lines,
        sprintf("  Tested ncomps: %i.\n", length(tested_ncomps))
      )
    }
    if (!is.null(optimal_ncomp) && !is.na(optimal_ncomp)) {
      lines <- c(
        lines,
        sprintf("  Optimal ncomp: %i.\n", optimal_ncomp)
      )
    }
    if (!is.null(ica_final)) {
      lines <- c(
        lines,
        sprintf("  Stabilised ICA fitted with %i components.\n",
                ica_final[["no_comp"]])
      )
    }
    if (!is.na(n_modules_from_result)) {
      lines <- c(
        lines,
        sprintf("  Modules identified: %i.\n", n_modules_from_result)
      )
    }
    paste(lines, collapse = "")
  } else if (detection_method == "dgrdl-based") {
    fit_params <- params[["fit_params"]]
    grid <- params[["grid_search_params"]]
    lines <- c(
      " Detection method: DGRDL (dual graph-regularised dictionary learning).\n"
    )
    if (!is.null(grid)) {
      total <- length(grid[["tested_dict_sizes"]]) *
        length(grid[["tested_k_neighbours"]]) *
        length(grid[["tested_seeds"]])
      lines <- c(
        lines,
        sprintf("  Grid search: %i parameter combinations tested.\n", total)
      )
    }
    if (!is.null(fit_params)) {
      lines <- c(
        lines,
        sprintf(
          "  Fitted dictionary: size = %i, k_neighbours = %i, sparsity = %i.\n",
          fit_params[["dict_size"]],
          fit_params[["k_neighbours"]],
          fit_params[["sparsity"]]
        )
      )
    }
    if (!is.na(n_modules_from_result)) {
      lines <- c(
        lines,
        sprintf("  Modules identified: %i.\n", n_modules_from_result)
      )
    }
    paste(lines, collapse = "")
  } else if (detection_method == "nmf-based") {
    nmf_params <- params[["nmf_fit"]]
    stabilised <- isTRUE(nmf_params[["stabilised"]])
    lines <- c(
      sprintf(
        " Detection method: %sNMF (HALS).\n",
        if (stabilised) "stabilised " else ""
      )
    )
    if (!is.null(nmf_params)) {
      lines <- c(
        lines,
        sprintf(
          "  k = %i, preprocessing = %s.\n",
          nmf_params[["k"]],
          nmf_params[["preprocessing"]]
        )
      )
      if (stabilised) {
        lines <- c(
          lines,
          sprintf("  n_runs = %i.\n", nmf_params[["n_runs"]])
        )
      } else if (!is.null(nmf_params[["converged"]])) {
        lines <- c(
          lines,
          sprintf(
            "  Converged: %s, iterations: %i.\n",
            nmf_params[["converged"]],
            nmf_params[["n_iter"]]
          )
        )
      }
    }
    if (!is.na(n_modules_from_result)) {
      lines <- c(
        lines,
        sprintf("  Modules identified: %i.\n", n_modules_from_result)
      )
    }
    paste(lines, collapse = "")
  } else {
    sprintf(" Detection method: %s.\n", detection_method)
  }

  cat(
    "Bulk co-expression module class (BulkCoExp).\n",
    " Pre-processing done: ",
    preprocessed,
    ".\n",
    hvg_line,
    method_info,
    sep = ""
  )

  invisible(x)
}

#' @noRd
S7::method(print, BulkDge) <- function(x, ...) {
  params <- S7::prop(x, "params")
  outputs <- S7::prop(x, "outputs")

  original_dim <- params[["original_dim"]]
  n_samples <- nrow(S7::prop(x, "meta_data"))
  variable_info <- S7::prop(x, "variable_info")
  has_variable_info <- !is.null(variable_info)

  # Downstream steps and their marker outputs
  step_flags <- list(
    "qc_bulk_dge()" = !is.null(outputs[["dge_list"]]),
    "normalise_bulk_dge()" = !is.null(outputs[["normalised_counts"]]),
    "batch_correction_bulk_dge()" =
      !is.null(outputs[["normalised_counts_corrected"]]),
    "calculate_pca_bulk_dge()" = !is.null(outputs[["pca"]]),
    "calculate_dge_limma()" = !is.null(outputs[["limma_voom_res"]]),
    "calculate_dge_hedges()" = !is.null(outputs[["hedges_g_res"]]),
    "TPM normalisation" = !is.null(outputs[["tpm_counts"]]),
    "FPKM normalisation" = !is.null(outputs[["fpkm_counts"]])
  )

  step_lines <- purrr::imap_chr(step_flags, \(flag, name) {
    sprintf("  %s: %s.\n", name, flag)
  })

  cat(
    "Bulk differential gene expression class (BulkDge).\n",
    sprintf(
      " Raw counts: %i genes x %i samples.\n",
      original_dim[1],
      original_dim[2]
    ),
    sprintf(" Meta-data rows: %i.\n", n_samples),
    sprintf(" Variable info provided: %s.\n", has_variable_info),
    " Applied steps:\n",
    paste(step_lines, collapse = ""),
    sep = ""
  )

  invisible(x)
}
