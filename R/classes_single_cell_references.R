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
#' @param hvg_gene_names Character vector of HVG gene names in reference
#' loading order.
#' @param gene_means Per-HVG mean of the normalised reference data.
#' @param gene_sds Per-HVG standard deviation of the normalised reference data.
#' @param loadings PCA gene loadings matrix (n_hvgs x d).
#' @param z_corr Post-Harmony corrected embedding (N x d).
#' @param centroids Cosine-normalised reference centroids (K x d).
#' @param nr Reference cluster sizes; row-sums of `r` (length K).
#' @param c_cache Cached `R * Z_corr` compression term (K x d).
#' @param no_pcs Number of principal components.
#' @param harmony_backend Which Harmony variant was used (`"v1"` or `"v2"`).
#' @param batch_vars Names of the batch variables used during reference
#' construction.
#' @param slim Logical; if `TRUE`, `z_orig` and `r` are dropped. Default
#' `FALSE`.
#' @param z_orig Pre-Harmony PCA scores (N x d). `NULL` in slim references.
#' @param r Soft cluster assignments (K x N). `NULL` in slim references.
#' @param labels Optional `data.table` of reference cell labels with
#' `nrow(z_corr)` rows, one column per label. `NULL` if no labels stored.
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
  ),
  constructor = function(
    hvg_gene_names,
    gene_means,
    gene_sds,
    loadings,
    z_corr,
    centroids,
    nr,
    c_cache,
    no_pcs,
    harmony_backend,
    batch_vars,
    slim = FALSE,
    z_orig = NULL,
    r = NULL,
    labels = NULL
  ) {
    S7::new_object(
      S7::S7_object(),
      hvg_gene_names = hvg_gene_names,
      gene_means = as.numeric(gene_means),
      gene_sds = as.numeric(gene_sds),
      loadings = loadings,
      z_orig = z_orig,
      z_corr = z_corr,
      r = r,
      centroids = centroids,
      nr = as.numeric(nr),
      c_cache = c_cache,
      no_pcs = as.integer(no_pcs),
      harmony_backend = harmony_backend,
      batch_vars = batch_vars,
      slim = slim,
      labels = labels
    )
  }
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
get_symphony_loadings <- S7::new_generic(
  name = "get_symphony_loadings",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

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
get_symphony_z_corr <- S7::new_generic(
  name = "get_symphony_z_corr",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

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
get_symphony_hvg_names <- S7::new_generic(
  name = "get_symphony_hvg_names",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

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
get_symphony_labels <- S7::new_generic(
  name = "get_symphony_labels",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

S7::method(get_symphony_labels, SymphonyReference) <- function(object) {
  S7::prop(object, "labels")
}

## primitives ------------------------------------------------------------------

### print ----------------------------------------------------------------------

#' @noRd
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
