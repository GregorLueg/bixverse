# classes related to reference generation for single cell

# symphony ---------------------------------------------------------------------

## class -----------------------------------------------------------------------

#' @title bixverse SymphonyReference class
#'
#' @description
#' Holds a Symphony reference: PCA loadings, per-HVG scaling stats, soft
#' cluster centroids and the cached compression terms (Nr, C) needed to map
#' queries without the reference cells. The Harmony-corrected reference
#' embedding (`z_corr`) is retained for downstream label transfer. Reference
#' cell labels (e.g. cell type annotations) can be snapshotted at build time
#' via `label_columns` in [build_symphony_ref()] or attached post-hoc via
#' [add_symphony_labels()]. For details on the method, refer to Kang et al.
#'
#' @section Properties:
#' \describe{
#'   \item{hvg_gene_names}{Character vector of HVG gene names in reference
#'   loading order.}
#'   \item{gene_means}{Per-HVG mean of the normalised reference data.}
#'   \item{gene_sds}{Per-HVG standard deviation of the normalised reference
#'   data.}
#'   \item{loadings}{PCA gene loadings matrix (n_hvgs x d).}
#'   \item{z_orig}{Pre-Harmony PCA scores (N x d). `NULL` in slim
#'   references.}
#'   \item{z_corr}{Post-Harmony corrected embedding (N x d). Always kept,
#'   even in slim references, since it backs kNN label transfer.}
#'   \item{r}{Soft cluster assignments (K x N). `NULL` in slim references.}
#'   \item{centroids}{Cosine-normalised reference centroids (K x d). Used for
#'   query soft clustering.}
#'   \item{nr}{Reference cluster sizes; row-sums of `r` (length K).}
#'   \item{c_cache}{Cached `R * Z_corr` compression term (K x d).}
#'   \item{no_pcs}{Number of principal components.}
#'   \item{harmony_backend}{Which Harmony variant was used to build the
#'   reference (`"v1"` or `"v2"`).}
#'   \item{batch_vars}{Names of the batch variables used during reference
#'   construction.}
#'   \item{slim}{Logical; if `TRUE`, `z_orig` and `r` are dropped to reduce
#'   memory footprint.}
#'   \item{labels}{Optional `data.table` of reference cell labels with
#'   `nrow(z_corr)` rows, one column per label. `NULL` if no labels were
#'   stored.}
#' }
#'
#' @return Returns the `SymphonyReference` class for further operations.
#'
#' @references Kang et al., Nat. Commun., 2021
#'
#' @export
SymphonyReference <- S7::new_class(
  name = "SymphonyReference",
  properties = list(
    hvg_gene_names = S7::class_character,
    gene_means = S7::class_numeric,
    gene_sds = S7::class_numeric,
    loadings = S7::class_any,
    z_orig = S7::class_any,
    z_corr = S7::class_any,
    r = S7::class_any,
    centroids = S7::class_any,
    nr = S7::class_numeric,
    c_cache = S7::class_any,
    no_pcs = S7::class_integer,
    harmony_backend = S7::class_character,
    batch_vars = S7::class_character,
    slim = S7::class_logical,
    labels = S7::class_any
  )
)

## methods ---------------------------------------------------------------------

### getters --------------------------------------------------------------------

#' Getter for the PCA loadings of a Symphony reference
#'
#' @param object `SymphonyReference` class.
#'
#' @return The PCA gene loadings matrix (n_hvgs x d).
#'
#' @export
get_symphony_loadings <- S7::new_generic("get_symphony_loadings", "object")

S7::method(get_symphony_loadings, SymphonyReference) <- function(object) {
  S7::prop(object, "loadings")
}

#' Getter for the corrected embedding of a Symphony reference
#'
#' @param object `SymphonyReference` class.
#'
#' @return The post-Harmony corrected embedding (N x d).
#'
#' @export
get_symphony_z_corr <- S7::new_generic("get_symphony_z_corr", "object")

S7::method(get_symphony_z_corr, SymphonyReference) <- function(object) {
  S7::prop(object, "z_corr")
}

#' Getter for the HVG gene names of a Symphony reference
#'
#' @param object `SymphonyReference` class.
#'
#' @return Character vector of HVG gene names in reference loading order.
#'
#' @export
get_symphony_hvg_names <- S7::new_generic("get_symphony_hvg_names", "object")

S7::method(get_symphony_hvg_names, SymphonyReference) <- function(object) {
  S7::prop(object, "hvg_gene_names")
}

#' Getter for the stored labels of a Symphony reference
#'
#' @param object `SymphonyReference` class.
#'
#' @return A `data.table` of reference cell labels in `z_corr` row order, or
#' `NULL` if no labels are stored.
#'
#' @export
get_symphony_labels <- S7::new_generic("get_symphony_labels", "object")

S7::method(get_symphony_labels, SymphonyReference) <- function(object) {
  S7::prop(object, "labels")
}

## primitives ------------------------------------------------------------------

### print ----------------------------------------------------------------------

#' @name print.SymphonyReference
#'
#' @title print Method for SymphonyReference object
#'
#' @description
#' Print a SymphonyReference object.
#'
#' @param x An object of class `SymphonyReference`.
#' @param ... Additional arguments (currently not used).
#'
#' @returns Invisibly returns `x`.
#'
#' @method print SymphonyReference
#'
#' @keywords internal
S7::method(print, SymphonyReference) <- function(x, ...) {
  loadings <- S7::prop(x, "loadings")
  centroids <- S7::prop(x, "centroids")
  labels <- S7::prop(x, "labels")
  label_str <- if (is.null(labels)) {
    "none"
  } else {
    paste(names(labels), collapse = ", ")
  }
  cat(
    "Symphony reference\n",
    sprintf("  Harmony backend: %s\n", S7::prop(x, "harmony_backend")),
    sprintf("  No HVGs: %i\n", nrow(loadings)),
    sprintf("  No PCs: %i\n", ncol(loadings)),
    sprintf("  No clusters: %i\n", nrow(centroids)),
    sprintf(
      "  Batch variables: %s\n",
      paste(S7::prop(x, "batch_vars"), collapse = ", ")
    ),
    sprintf("  Slim: %s\n", S7::prop(x, "slim")),
    sprintf("  Labels: %s\n", label_str),
    sep = ""
  )
  invisible(x)
}
