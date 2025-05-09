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

# general methods --------------------------------------------------------------

## bulk_coexp ------------------------------------------------------------------

#' Process the raw data
#'
#' @description
#' Function to do general pre-processing on top of the [bixverse::bulk_coexp()].
#' Options to do scaling, HVG selection, etc.
#'
#' @param object The underlying class, see [bixverse::bulk_coexp()].
#' @param hvg Integer or float. If an integer, the top `hvg` genes will be
#' included; if float, the float has to be between 0 and 1, representing the
#' percentage of genes to include.
#' @param mad_threshold Float. Instead of of selecting number or proportion of
#' genes, you can also provide a mad_threshold.
#' @param scaling Boolean. Shall the data be scaled.
#' @param scaling_type String. You have the option to use normal scaling or
#' robust scaling.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return Returns the class with the `processed_data` data slot populated and
#' applied parameters added to the `params` slot.
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
preprocess_bulk_coexp <- S7::new_generic(
  "preprocess_bulk_coexp",
  "object",
  fun = function(
    object,
    hvg = NULL,
    mad_threshold = NULL,
    scaling = FALSE,
    scaling_type = c("normal", "robust"),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method preprocess_bulk_coexp bulk_coexp
#' @export
S7::method(preprocess_bulk_coexp, bulk_coexp) <- function(
  object,
  hvg = NULL,
  mad_threshold = NULL,
  scaling = FALSE,
  scaling_type = c("normal", "robust"),
  .verbose = TRUE
) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(mad_threshold, c("R1", "0"))
  nfeatures <- S7::prop(object, "params")[["original_dim"]][2]
  checkmate::qassert(hvg, c("R1[0,1]", sprintf("I1[0,%i]", nfeatures), "0"))
  checkmate::qassert(scaling, "B1")
  if (scaling) {
    checkmate::assertChoice(scaling_type, c("normal", "robust"))
  }

  mat <- S7::prop(object, "raw_data")

  feature_meta <- data.table::data.table(
    feature_name = colnames(mat),
    mean_exp = colMeans(mat),
    MAD = matrixStats::colMads(mat),
    var_exp = matrixStats::colVars(mat)
  ) %>%
    data.table::setorder(-MAD)

  if (!is.null(hvg) & !is.null(mad_threshold)) {
    choice <- menu(
      choices = c("MAD threshold", "HVG threshold"),
      title = paste(
        "You have provided both a MAD and HVG",
        "threshold. Which one do you want to use?"
      )
    )
    if (choice == 1) {
      hvg <- NULL
    } else {
      mad_threshold <- NULL
    }
  }

  if (is.null(hvg) & is.null(mad_threshold)) {
    hvg <- 1
  }

  if (!is.null(hvg)) {
    no_genes_to_take <-
      ifelse(is.integer(hvg), hvg, ceiling(hvg * ncol(mat)))
    hvg_genes <- feature_meta[1:no_genes_to_take, feature_name]
  } else {
    hvg_genes <- feature_meta[MAD >= mad_threshold, feature_name]
  }

  if (.verbose) {
    message(sprintf("A total of %i genes will be included.", length(hvg_genes)))
  }

  feature_meta[, hvg := feature_name %in% hvg_genes]

  # Process the matrix
  matrix_processed <- mat[, hvg_genes]

  if (scaling) {
    fun <-
      ifelse(scaling_type == "normal", "scale", "bixverse::robust_scale")
    matrix_processed <- rlang::eval_tidy(rlang::quo(apply(
      matrix_processed,
      1,
      !!!rlang::parse_exprs(fun)
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

  S7::prop(object, "params")[["preprocessing"]] <-
    processing_params
  S7::prop(object, "processed_data")[["processed_data"]] <-
    matrix_processed
  S7::prop(object, "processed_data")[["feature_meta"]] <-
    feature_meta

  # Return
  object
}

## bulk_dge --------------------------------------------------------------------

#' QC on the bulk dge data
#'
#' @description
#' This function generates general QC metrics for bulk DGE experiments.
#'
#' @param object The underlying class, see [bixverse::bulk_dge()].
#' @param group_col String. The column in the metadata that will contain the
#' contrast groups. Needs to be part of the metadata stored in the class.
#' @param norm_method String. One of `c("TMM", "TMMwsp", "RLE", "upperquartile",`
#' ` "none")`. Please refer to [edgeR::normLibSizes()].
#' @param outlier_threshold Float. Number of standard deviations in terms of
#' percentage genes detected you allow before removing a sample. Defaults to `2`.
#' @param min_prop Float. Minimum proportion of samples in which the gene has
#' to be identified in.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return ...
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
preprocess_bulk_dge <- S7::new_generic(
  "preprocess_bulk_dge",
  "object",
  fun = function(
    object,
    group_col,
    norm_method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
    outlier_threshold = 2,
    min_prop = 0.2,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)


#' @method preprocess_bulk_dge bulk_dge
#'
#' @import ggplot2
#'
#' @export
S7::method(preprocess_bulk_dge, bulk_dge) <- function(
  object,
  group_col,
  norm_method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
  outlier_threshold = 2,
  min_prop = 0.2,
  .verbose = TRUE
) {
  norm_method <- match.arg(norm_method)
  # Checks
  checkmate::assertClass(
    object,
    "bixverse::bulk_dge"
  )
  checkmate::qassert(group_col, "S1")
  checkmate::assertChoice(
    norm_method,
    c("TMM", "TMMwsp", "RLE", "upperquartile", "none")
  )
  checkmate::qassert(outlier_threshold, "N1")
  checkmate::qassert(min_prop, "R1[0,1]")
  checkmate::qassert(.verbose, "B1")

  raw_counts <- S7::prop(object, "raw_counts")
  meta_data <- S7::prop(object, "meta_data")

  checkmate::assertTRUE(group_col %in% colnames(meta_data))

  # Sample outlier removal
  if (.verbose) message("Detecting sample outliers.")
  detected_genes_nb <- data.table::data.table(
    sample_id = colnames(raw_counts),
    nb_detected_genes = matrixStats::colSums2(raw_counts > 0)
  )

  samples <- merge(meta_data, detected_genes_nb, by = "sample_id") %>%
    .[, `:=`(perc_detected_genes = nb_detected_genes / nrow(raw_counts))]

  p1_nb_genes_cohort <- ggplot(
    samples,
    aes(
      x = .data[[group_col]],
      y = nb_detected_genes,
      fill = .data[[group_col]]
    )
  ) +
    geom_boxplot(alpha = 0.5) +
    labs(
      x = "Groups",
      y = "Number of genes detected",
      title = "Number of genes by cohort"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "none"
    )

  sd_samples = sd(samples$perc_detected_genes, na.rm = TRUE)
  min_perc = mean(samples$perc_detected_genes, na.rm = TRUE) -
    outlier_threshold * sd_samples
  max_perc = mean(samples$perc_detected_genes, na.rm = TRUE) +
    outlier_threshold * sd_samples

  p2_outliers <- ggplot(
    samples,
    aes(x = 1, y = perc_detected_genes, color = .data[[group_col]])
  ) +
    geom_point(
      position = position_jitter(width = 0.2, height = 0, seed = 123),
      size = 3,
      alpha = 0.7
    ) +
    geom_hline(yintercept = min_perc, color = "red", linetype = "dashed") +
    geom_hline(yintercept = max_perc, color = "red", linetype = "dashed") +
    labs(
      x = "",
      y = "Percentage of genes detected",
      title = "Outlier samples based on %genes detected"
    ) +
    theme_classic() +
    theme(legend.position = "bottom", axis.text.x = element_blank())

  outliers <- samples$perc_detected_genes <= min_perc
  if (.verbose)
    message(sprintf(
      "A total of %i samples are detected as outlier.",
      sum(outliers)
    ))

  samples_red <- samples[!(outliers), ]
  raw_counts <- raw_counts[, unique(samples_red$sample_id)]

  ## Voom normalization
  if (.verbose)
    message("Removing lowly expressed genes and normalising counts.")
  dge_list <- edgeR::DGEList(raw_counts)

  # Filter lowly expressed genes
  to_keep <- edgeR::filterByExpr(
    y = dge_list,
    min.prop = min_prop,
    group = samples[[group_col]]
  )
  if (.verbose) message(sprintf("A total of %i genes are kept.", sum(to_keep)))

  dge_list <- edgeR::calcNormFactors(
    dge_list[to_keep, ],
    method = norm_method
  )

  samples <- samples_red[[group_col]]

  design <- model.matrix(~ 0 + samples)

  voom_obj <- limma::voom(
    counts = dge_list,
    design = design,
    normalize.method = "quantile",
    plot = TRUE
  )

  p3_voom_normalization <- recordPlot()
  dev.off()

  boxplot(
    voom_obj$E,
    col = as.factor(samples),
    main = "Boxplot of normalized gene expression across samples"
  )
  p4_boxplot_normalization <- recordPlot()
  dev.off()

  S7::prop(object, "plots") <- list(
    p1_nb_genes_cohort = p1_nb_genes_cohort,
    p2_outliers = p2_outliers,
    p3_voom_normalization = p3_voom_normalization,
    p4_boxplot_normalization = p4_boxplot_normalization
  )

  S7::prop(object, "outputs") <- list(
    dge_list = dge_list,
    sample_info = samples_red,
    normalised_counts = voom_obj$E,
    group_col = group_col
  )

  S7::prop(object, "params")[["QC_params"]] <- list(
    normalisation_method = norm_method,
    min_prop = min_prop,
    outlier_threshold = outlier_threshold
  )

  return(object)
}


# plots ------------------------------------------------------------------------

## bulk_coexp ------------------------------------------------------------------

#' @title Plot the highly variable genes
#'
#' @description
#' Plots the median-absolute deviation of the genes and applied thresholds.
#' Expects that [bixverse::preprocess_bulk_coexp()] was run and will throw an
#' error otherwise.
#'
#' @param object The underlying class, see [bixverse::bulk_coexp()].
#' @param bins Integer. Number of bins to plot.
#'
plot_hvgs <- S7::new_generic(
  "plot_hvgs",
  "object",
  fun = function(object, bins = 50L) {
    S7::S7_dispatch()
  }
)

#' @method plot_hvgs bulk_coexp
#'
#' @import ggplot2
#'
#' @export
S7::method(plot_hvgs, bulk_coexp) <- function(object, bins = 50L) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(bins, "I1")
  # Early return
  if (is.null(S7::prop(object, "params")[["preprocessing"]])) {
    warning("No pre-processing data found. Returning NULL.")
    return(NULL)
  }
  plot_df <- S7::prop(object, "processed_data")[["feature_meta"]]

  p <- ggplot(data = plot_df, mapping = aes(x = MAD)) +
    geom_histogram(
      mapping = aes(fill = hvg),
      bins = 50L,
      color = "black",
      alpha = 0.7
    ) +
    scale_fill_manual(values = setNames(c("orange", "grey"), c(TRUE, FALSE))) +
    ggtitle(
      "Distribution of MAD across the genes",
      subtitle = "And included genes"
    ) +
    theme_minimal()

  mad_threshold <- S7::prop(object, "params")[["preprocessing"]][[
    "mad_threshold"
  ]]

  if (mad_threshold != "not applicable") {
    p <- p +
      geom_vline(
        xintercept = mad_threshold,
        linetype = "dashed",
        color = "darkred"
      )
  }

  return(p)
}

## bulk dge --------------------------------------------------------------------

#' Return QC plots
#'
#' @description
#' Getter function to extract the QC plots from the [bixverse::bulk_dge()]
#' class. These are added when you run for example
#' [bixverse::preprocess_bulk_dge()]. You can either leave the plot choice as
#' `NULL` and provide input when prompted, or you provide the name. The possible
#' plots that might be in the class
#' \itemize{
#'  \item p1_nb_genes_cohort Proportion of non-zero genes for the samples in the
#'  respective cohorts
#'  \item p2_outliers An outlier plot based on the data from p1
#'  \item p3_voom_normalization Initial Voom normalisation plot after filtering
#'  lowly expressed genes
#'  \item p4_boxplot_normalization Expression levels after normalisation
#'  \item p5_pca_case_control A PCA plot with the chosen case control category.
#'  Added if [bixverse::calculate_pca_bulk_dge()] is run.
#'  \item p6_batch_correction_plot A PCA plot pre and post batch correction
#'  with the case-control category overlayed. Added if [bixverse::batch_correction_bulk_dge()]
#'  is run.
#' }
#'
#' @param object `bulk_dge` class.
#' @param plot_choice Optional string or integer. Index or name of the plate.
#'
#' @return Returns the DGEList stored in the class.
#'
#' @export
get_dge_qc_plot <- S7::new_generic(
  name = "get_dge_qc_plot",
  dispatch_args = "object",
  fun = function(object, plot_choice = NULL) {
    S7::S7_dispatch()
  }
)

#' @method get_dge_qc_plot bulk_dge
#'
#' @export
S7::method(get_dge_qc_plot, bulk_dge) <-
  function(object, plot_choice = NULL) {
    # Checks
    checkmate::assertClass(
      object,
      "bixverse::bulk_dge"
    )
    checkmate::qassert(plot_choice, c("0", "I1", "S1"))

    plots <- S7::prop(object, "plots")
    plot_names <- names(plots)

    user_choice <- ifelse(
      is.null(plot_choice),
      select_user_option(plot_names),
      plot_choice
    )

    # Return
    return(plots[[user_choice]])
  }
