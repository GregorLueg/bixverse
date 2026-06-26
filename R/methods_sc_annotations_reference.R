# cell type reference and annotation methods -----------------------------------

## sc type ---------------------------------------------------------------------

#' Calculate ScType scores per cell
#'
#' @description
#' Implements the approach from
#'
#' @param object `SingleCells` or `SingleCellsMultiModal`
#' @param cell_marker A list, see [prepare_cell_markers()].
#' @param sensitivity Boolean. Shall shared marker genes be downweighted (like
#' in the original reference). Defaults to `TRUE`.
#' @param weight_floor Optional numeric. A value between 0 to 1 and sets the
#' floor for the weights if `sensitivity = TRUE`. If not provided, defaults to
#' `0` as in the original implementation.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns An `ScTypeResults` results class
#'
#' @export
calc_sc_type_scores <- S7::new_generic(
  name = "calc_sc_type_scores",
  dispatch_args = "object",
  fun = function(
    object,
    cell_marker_list,
    sensitivity = TRUE,
    weight_floor = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calc_sc_type_scores SingleCells
S7::method(calc_sc_type_scores, SingleCells) <- function(
  object,
  cell_marker_list,
  sensitivity = TRUE,
  weight_floor = NULL,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertCellMarkerList(cell_marker_list)
  checkmate::qassert(.verbose, c("I1", "B1"))

  res <- rs_sc_type(
    f_path = get_rust_count_gene_f_path(object),
    cell_indices = get_cells_to_keep(object),
    cell_markers = cell_marker_list,
    sensitivity = sensitivity,
    weight_floor = weight_floor,
    verbose = parse_verbosity(.verbose)
  )

  class(res) <- "ScTypeResults"

  return(res)
}

## symphony --------------------------------------------------------------------

### build symphony reference ---------------------------------------------------

#' Build a Symphony reference from a SingleCells object
#'
#' @description
#' Runs sparse PCA over the provided HVGs, runs Harmony for batch correction,
#' and compresses the result into the cached terms used at query time. The
#' Harmony version is auto-detected from the class of `harmony_params`.
#'
#' @param object `SingleCells` (the reference).
#' @param batch_column String. Primary batch column in the obs table.
#' @param additional_batch_columns Optional character vector.
#' @param hvg Integer vector. R-style 1-based HVG indices into the gene
#' universe of `object`. Must be provided explicitly.
#' @param harmony_params List. Output of [bixverse::params_sc_harmony()] or
#' [bixverse::params_sc_harmony_v2()].
#' @param pca_params List. Output of [bixverse::params_sc_pca()].
#' @param no_pcs Integer. Number of principal components.
#' @param slim Boolean. If `TRUE`, drops `z_orig` and `r` from the returned
#' reference. `z_corr` is always kept (needed for kNN label transfer).
#' @param seed Integer.
#' @param .verbose Boolean or integer.
#'
#' @return A [SymphonyReference] object.
#'
#' @export
build_symphony_ref <- S7::new_generic(
  name = "build_symphony_ref",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    additional_batch_columns = NULL,
    hvg,
    harmony_params = params_sc_harmony(),
    pca_params = params_sc_pca(),
    no_pcs = 30L,
    slim = FALSE,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method build_symphony_ref SingleCells
#'
#' @export
S7::method(build_symphony_ref, SingleCells) <- function(
  object,
  batch_column,
  additional_batch_columns = NULL,
  hvg,
  harmony_params = params_sc_harmony(),
  pca_params = params_sc_pca(),
  no_pcs = 30L,
  slim = FALSE,
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(additional_batch_columns, c("S+", "0"))
  checkmate::qassert(hvg, "I+")
  assertScPca(pca_params)
  checkmate::qassert(no_pcs, "I1[1,)")
  checkmate::qassert(slim, "B1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # detect harmony backend from class attribute
  harmony_version <- if (inherits(harmony_params, "params_sc_harmony_v2")) {
    assertScHarmonyParamsV2(harmony_params)
    "v2"
  } else if (inherits(harmony_params, "params_sc_harmony")) {
    assertScHarmonyParams(harmony_params)
    "v1"
  } else {
    stop(
      paste(
        "harmony_params must come from params_sc_harmony()",
        "or params_sc_harmony_v2()"
      )
    )
  }

  # resolve HVG names from the reference's gene universe
  gene_names <- get_gene_names(object)
  if (any(hvg < 1L) || any(hvg > length(gene_names))) {
    stop("hvg contains out-of-range indices.")
  }
  hvg_gene_names <- gene_names[hvg]
  hvg_rust <- as.integer(hvg) - 1L

  # batch labels (same pattern as harmony_sc)
  batch_columns_all <- c(batch_column, additional_batch_columns)
  batch_index_ls <- lapply(batch_columns_all, function(col) {
    vals <- unlist(object[[col]])
    as.integer(factor(vals)) - 1L
  })

  cells_to_keep <- get_cells_to_keep(object)
  checkmate::assertTRUE(all(
    purrr::map_dbl(batch_index_ls, length) == length(cells_to_keep)
  ))

  # auto-determine harmony k
  if (is.null(harmony_params$k)) {
    harmony_params$k <- as.integer(min(round(length(cells_to_keep) / 30), 100L))
    if (.verbose) {
      message(sprintf(
        " Auto-determined number of Harmony clusters: %d",
        harmony_params$k
      ))
    }
  }

  ref_rs <- rs_build_symphony_ref(
    f_path = get_rust_count_gene_f_path(object),
    f_path_cell = get_rust_count_cell_f_path(object),
    cell_indices = cells_to_keep,
    hvg_indices = hvg_rust,
    batch_labels = batch_index_ls,
    pca_params = pca_params,
    no_pcs = as.integer(no_pcs),
    harmony_params = harmony_params,
    harmony_version = harmony_version,
    seed = as.integer(seed),
    verbose = parse_verbosity(.verbose)
  )

  z_orig <- if (slim) NULL else ref_rs$z_orig
  r_field <- if (slim) NULL else ref_rs$r

  SymphonyReference(
    hvg_gene_names = hvg_gene_names,
    gene_means = as.numeric(ref_rs$gene_means),
    gene_sds = as.numeric(ref_rs$gene_sds),
    loadings = ref_rs$loadings,
    z_orig = z_orig,
    z_corr = ref_rs$z_corr,
    r = r_field,
    centroids = ref_rs$centroids,
    nr = as.numeric(ref_rs$nr),
    c_cache = ref_rs$c,
    no_pcs = as.integer(no_pcs),
    harmony_backend = harmony_version,
    batch_vars = batch_columns_all,
    slim = slim
  )
}

### map symphony query ---------------------------------------------------------

#' Map a SingleCells query onto a Symphony reference
#'
#' @description
#' Projects query cells through the reference's PCA loadings and applies the
#' cached MoE batch correction. Gene matching is by name against the
#' reference's stored HVG names; unmatched HVGs become zero columns. Results
#' are attached to the query as embeddings.
#'
#' @param reference `SymphonyReference`.
#' @param query `SingleCells` query.
#' @param batch_column Optional string. If `NULL`, no batch correction is
#' applied (z_corr = z_pca).
#' @param additional_batch_columns Optional character vector.
#' @param params List. Output of [bixverse::params_symphony_map()].
#' @param .verbose Boolean or integer.
#'
#' @return The `query` object with embeddings `"symphony"` (z_corr),
#' `"symphony_pca"` (z_pca) and `"symphony_r"` (soft cluster assignments,
#' transposed to N_q x K).
#'
#' @export
map_symphony_query <- S7::new_generic(
  name = "map_symphony_query",
  dispatch_args = "reference",
  fun = function(
    reference,
    query,
    batch_column = NULL,
    additional_batch_columns = NULL,
    params = params_symphony_map(),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method map_symphony_query SymphonyReference
#'
#' @export
S7::method(map_symphony_query, SymphonyReference) <- function(
  reference,
  query,
  batch_column = NULL,
  additional_batch_columns = NULL,
  params = params_symphony_map(),
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(query, SingleCells))
  checkmate::qassert(batch_column, c("S1", "0"))
  checkmate::qassert(additional_batch_columns, c("S+", "0"))
  assertSymphonyMap(params)
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # gene matching: reference HVG names -> query gene index (0-based, NA if absent)
  query_gene_names <- get_gene_names(query)
  hvg_names <- S7::prop(reference, "hvg_gene_names")
  query_idx_1based <- match(hvg_names, query_gene_names)
  n_missing <- sum(is.na(query_idx_1based))
  if (n_missing > 0L && .verbose) {
    message(sprintf(
      " %d / %d reference HVGs not found in query; treated as zero columns.",
      n_missing,
      length(hvg_names)
    ))
  }
  ref_to_query_gene_map <- as.integer(query_idx_1based - 1L) # NA preserved

  # query batch labels
  batch_index_ls <- if (is.null(batch_column)) {
    list()
  } else {
    batch_columns_all <- c(batch_column, additional_batch_columns)
    lapply(batch_columns_all, function(col) {
      vals <- unlist(query[[col]])
      as.integer(factor(vals)) - 1L
    })
  }

  cells_to_keep <- get_cells_to_keep(query)
  if (length(batch_index_ls) > 0L) {
    checkmate::assertTRUE(all(
      purrr::map_dbl(batch_index_ls, length) == length(cells_to_keep)
    ))
  }

  res <- rs_symphony_map_query(
    f_path_query = get_rust_count_gene_f_path(query),
    cell_indices_query = cells_to_keep,
    gene_means = S7::prop(reference, "gene_means"),
    gene_sds = S7::prop(reference, "gene_sds"),
    loadings = S7::prop(reference, "loadings"),
    centroids = S7::prop(reference, "centroids"),
    nr = S7::prop(reference, "nr"),
    c_cache = S7::prop(reference, "c_cache"),
    ref_to_query_gene_map = ref_to_query_gene_map,
    batch_labels_query = batch_index_ls,
    params_symphony = params,
    verbose = parse_verbosity(.verbose)
  )

  z_corr <- res$z_corr
  colnames(z_corr) <- sprintf("symphony_%s", seq_len(ncol(z_corr)))

  z_pca <- res$z_pca
  colnames(z_pca) <- sprintf("symphony_pca_%s", seq_len(ncol(z_pca)))

  r_embd <- t(res$r) # Rust returns K x N_q; transpose to N_q x K
  colnames(r_embd) <- sprintf("symphony_r_%s", seq_len(ncol(r_embd)))

  query <- set_embedding(query, embd = z_corr, name = "symphony")
  query <- set_embedding(query, embd = z_pca, name = "symphony_pca")
  query <- set_embedding(query, embd = r_embd, name = "symphony_r")

  return(query)
}

### label transfer -------------------------------------------------------------

#' Transfer labels from a Symphony reference to a query via kNN majority vote
#'
#' @description
#' kNN majority vote on the Harmony-corrected embeddings. The reference's
#' `z_corr` is the searchable index; the query's `"symphony"` embedding is
#' what gets queried. Reference labels are read from `reference_object`'s
#' obs in `cells_to_keep` order — that must match the cells used at
#' reference-build time.
#'
#' Distances on `z_corr` are typically Euclidean — set `knn_params$ann_dist`
#' to `"euclidean"` unless you have a reason to use cosine.
#'
#' @param reference `SymphonyReference` (must have `z_corr`).
#' @param reference_object `SingleCells` used to build the reference.
#' @param query `SingleCells` with a `"symphony"` embedding attached.
#' @param label_column String. Column in the reference obs table.
#' @param knn_params List. Output of [params_sc_knn()].
#' @param seed Integer.
#' @param .verbose Boolean or integer.
#'
#' @return A data.table with columns `predicted_<label_column>` and
#' `confidence_<label_column>`, in `get_cells_to_keep(query)` order.
#'
#' @export
transfer_labels_symphony <- S7::new_generic(
  name = "transfer_labels_symphony",
  dispatch_args = "reference",
  fun = function(
    reference,
    reference_object,
    query,
    label_column,
    knn_params = params_sc_knn(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method transfer_labels_symphony SymphonyReference
#'
#' @export
S7::method(transfer_labels_symphony, SymphonyReference) <- function(
  reference,
  reference_object,
  query,
  label_column,
  knn_params = params_sc_knn(),
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(reference_object, SingleCells))
  checkmate::assertTRUE(S7::S7_inherits(query, SingleCells))
  checkmate::qassert(label_column, "S1")
  assertScKnn(knn_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  ref_z_corr <- S7::prop(reference, "z_corr")
  if (is.null(ref_z_corr)) {
    stop("SymphonyReference has no z_corr. Rebuild with slim = FALSE.")
  }

  query_embd <- get_embedding(query, "symphony")
  if (is.null(query_embd)) {
    stop("Query has no 'symphony' embedding. Run map_symphony_query first.")
  }

  ref_obs <- get_sc_obs(reference_object, cols = label_column, filtered = TRUE)
  ref_labels_raw <- unlist(ref_obs[[label_column]])

  if (nrow(ref_z_corr) != length(ref_labels_raw)) {
    stop(sprintf(
      paste(
        "z_corr has %d rows but reference labels have %d entries;",
        "cells_to_keep on the reference object may have changed",
        "since the reference was built."
      ),
      nrow(ref_z_corr),
      length(ref_labels_raw)
    ))
  }

  ref_label_factor <- factor(ref_labels_raw)
  ref_label_idx <- as.integer(ref_label_factor) - 1L

  res <- rs_transfer_labels_symphony(
    reference_z_corr = ref_z_corr,
    query_z_corr = query_embd,
    reference_labels = ref_label_idx,
    n_labels = as.integer(nlevels(ref_label_factor)),
    knn_params = knn_params,
    seed = as.integer(seed),
    verbose = parse_verbosity(.verbose)
  )

  predicted <- levels(ref_label_factor)[res$predicted + 1L]

  out <- data.table::data.table(
    predicted = predicted,
    confidence = res$confidence
  )
  data.table::setnames(
    out,
    c("predicted", "confidence"),
    c(paste0("predicted_", label_column), paste0("confidence_", label_column))
  )
  out
}
