# deprecated calls -------------------------------------------------------------

## bulk coexp ------------------------------------------------------------------

#' @title Bulk RNAseq co-expression modules (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This constructor has been renamed to [bixverse::BulkCoExp()].
#'
#' @param raw_data The raw count matrix. Rows = samples, columns = features.
#' @param meta_data data.table. Metadata information on the samples. Expects to
#' have a `sample_id` column.
#' @param variable_info data.table. Metadata information on the features. This
#' is an optional table.
#'
#' @return Returns a [bixverse::BulkCoExp()] object.
#'
#' @keywords internal
#' @importFrom lifecycle deprecate_warn
#' @export
bulk_coexp <- function(raw_data, meta_data, variable_info = NULL) {
  lifecycle::deprecate_warn(
    when = "0.3.0",
    what = "bulk_coexp()",
    with = "BulkCoExp()"
  )
  BulkCoExp(
    raw_data = raw_data,
    meta_data = meta_data,
    variable_info = variable_info
  )
}

## bulk dge --------------------------------------------------------------------

#' @title Bulk RNAseq differential gene expression class (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This constructor has been renamed to [bixverse::BulkDge()].
#'
#' @param raw_counts matrix. The raw count matrix. Rows = genes, columns =
#' samples.
#' @param meta_data data.table. Metadata information on the samples. Expects to
#' have a `sample_id` column.
#' @param variable_info data.table. Metadata information on the features.
#' Defaults to `NULL`.
#' @param alternative_gene_id String. Optional alternative gene identifier to
#' be used. Must be a column of `variable_info`.
#'
#' @return Returns a [bixverse::BulkDge()] object.
#'
#' @keywords internal
#' @importFrom lifecycle deprecate_warn
#' @export
bulk_dge <- function(
  raw_counts,
  meta_data,
  variable_info = NULL,
  alternative_gene_id = NULL
) {
  lifecycle::deprecate_warn(
    when = "0.3.0",
    what = "bulk_dge()",
    with = "BulkDge()"
  )
  BulkDge(
    raw_counts = raw_counts,
    meta_data = meta_data,
    variable_info = variable_info,
    alternative_gene_id = alternative_gene_id
  )
}

## network diffusion -----------------------------------------------------------

#' @title Network diffusion class (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This constructor has been renamed to [bixverse::NetworkDiffusions()].
#'
#' @param edge_data_frame data.table that contains the edge information. It is
#' expected to have the columns 'from' and 'to'.
#' @param weighted Boolean. Is the graph weighted. If set to TRUE, the
#' `edge_data_frame` needs to have a weight column.
#' @param directed Boolean. Shall the graph be stored as directed.
#'
#' @return Returns a [bixverse::NetworkDiffusions()] object.
#'
#' @keywords internal
#' @importFrom lifecycle deprecate_warn
#' @export
network_diffusions <- function(edge_data_frame, weighted, directed) {
  lifecycle::deprecate_warn(
    when = "0.3.0",
    what = "network_diffusions()",
    with = "NetworkDiffusions()"
  )
  NetworkDiffusions(
    edge_data_frame = edge_data_frame,
    weighted = weighted,
    directed = directed
  )
}

## rbh graph -------------------------------------------------------------------

#' @title Reciprocal best hit graph (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This constructor has been renamed to [bixverse::RbhGraph()].
#'
#' @param module_results data.table or list containing the gene modules.
#' @param rbh_type String. One of `c("cor", "set")`.
#' @param dataset_col The column indicating the data set/method origin.
#' @param module_col The column storing the module names.
#' @param value_col The column storing the genes within each module.
#'
#' @return Returns a [bixverse::RbhGraph()] object.
#'
#' @keywords internal
#' @importFrom lifecycle deprecate_warn
#' @export
rbh_graph <- function(
  module_results,
  rbh_type = c("cor", "set"),
  dataset_col = NULL,
  module_col = NULL,
  value_col = NULL
) {
  lifecycle::deprecate_warn(
    when = "0.3.0",
    what = "rbh_graph()",
    with = "RbhGraph()"
  )
  RbhGraph(
    module_results = module_results,
    rbh_type = rbh_type,
    dataset_col = dataset_col,
    module_col = module_col,
    value_col = value_col
  )
}

## snf -------------------------------------------------------------------------

#' @title Similarity network fusion (deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This constructor has been renamed to [bixverse::Snf()].
#'
#' @param data Optional data to transform into adjacency data. Can be a
#' data.table (categorical/mixed) or a matrix (continuous).
#' @param data_name Optional string. Name of the data modality.
#' @param snf_params List. The SNF parameters, see [bixverse::params_snf()].
#'
#' @return Returns a [bixverse::Snf()] object.
#'
#' @keywords internal
#' @importFrom lifecycle deprecate_warn
#' @export
snf <- function(data = NULL, data_name = NULL, snf_params = params_snf()) {
  lifecycle::deprecate_warn(
    when = "0.3.0",
    what = "snf()",
    with = "Snf()"
  )
  Snf(data = data, data_name = data_name, snf_params = snf_params)
}
