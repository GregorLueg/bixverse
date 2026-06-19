# symphony class ---------------------------------------------------------------

## class -----------------------------------------------------------------------

#' @title bixverse SymphonyReference class
#'
#' @description
#' Holds a Symphony reference: PCA loadings, per-HVG scaling stats, soft
#' cluster centroids and the cached compression terms (Nr, C) needed to map
#' queries without the reference cells. The Harmony-corrected reference
#' embedding (`z_corr`) is retained for downstream label transfer.
#'
#' @section Properties:
#' \describe{
#'   \item{hvg_gene_names}{Character vector of HVG gene names in reference
#'   loading order.}
#'   \item{gene_means}{Per-HVG mean of the normalised reference data.}
#'   \item{gene_sds}{Per-HVG standard deviation of the normalised reference
#'   data.}
#'   \item{loadings}{PCA gene loadings matrix (n_hvgs x d).}
#'   \item{z_orig}{Pre-Harmony PCA scores (N x d).}
#'   \item{z_corr}{Post-Harmony corrected embedding (N x d). `NULL` in slim
#'   references.}
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
#'   \item{slim}{Logical; if `TRUE`, `z_orig`, `z_corr` and `r` are dropped
#'   to reduce memory footprint.}
#' }
#'
#' @return Returns the `SymphonyReference` class for further operations.
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
    slim = S7::class_logical
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
#' @return The post-Harmony corrected embedding (N x d), or `NULL` for slim
#' references.
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
    sep = ""
  )
  invisible(x)
}
