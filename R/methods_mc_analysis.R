# meta cell analysis methods ---------------------------------------------------

## aucell ----------------------------------------------------------------------

#' @method aucell_sc MetaCells
#'
#' @export
S7::method(aucell_sc, MetaCells) <- function(
  object,
  gs_list,
  auc_type = c("wilcox", "auroc"),
  streaming = FALSE,
  .verbose = TRUE
) {
  auc_type <- match.arg(auc_type)

  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::assertList(gs_list, types = "character", names = "named")
  checkmate::assertChoice(auc_type, c("wilcox", "auroc"))
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # 0-indexed gene indices for Rust
  gs_indices <- purrr::map(gs_list, \(e) {
    get_gene_indices(x = object, gene_ids = e, rust_index = FALSE) - 1L
  })

  sparse_data <- mc_counts_to_list(object, assay = "raw")

  auc_res <- rs_mc_aucell(
    sparse_data = sparse_data,
    gs_list = gs_indices,
    auc_type = auc_type,
    verbose = parse_verbosity(.verbose)
  )

  colnames(auc_res) <- names(gs_list)
  rownames(auc_res) <- S7::prop(object, "obs_table")$meta_cell_id

  return(auc_res)
}

## scenic ----------------------------------------------------------------------

### filtering ------------------------------------------------------------------

#' @method scenic_gene_filter_sc MetaCells
#'
#' @export
S7::method(scenic_gene_filter_sc, MetaCells) <- function(
  object,
  scenic_params = params_scenic(),
  cells_to_take = NULL,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  assertScenicParams(scenic_params)
  checkmate::qassert(cells_to_take, c("S+", "0"))
  checkmate::qassert(.verbose, "B1")

  cell_indices <- if (is.null(cells_to_take)) {
    NULL
  } else {
    get_cell_indices(object, cell_ids = cells_to_take, rust_index = FALSE)
  }

  counts <- get_sc_counts(
    object,
    assay = "raw",
    cell_indices = cell_indices
  )

  total_counts <- Matrix::colSums(counts)
  expressed_frac <- Matrix::colSums(counts != 0) / nrow(counts)

  passing <- which(
    total_counts >= scenic_params$min_counts &
      expressed_frac >= scenic_params$min_cells
  )

  if (.verbose) {
    message(sprintf(
      "SCENIC gene filter: %d / %d genes pass.",
      length(passing),
      ncol(counts)
    ))
  }

  passing_gene_ids <- S7::prop(object, "var_table")$gene_id[passing]

  return(passing_gene_ids)
}

### scenic grn -----------------------------------------------------------------

#' @method scenic_grn_sc MetaCells
#'
#' @export
S7::method(scenic_grn_sc, MetaCells) <- function(
  object,
  tf_ids,
  scenic_params = params_scenic(),
  genes_to_take = NULL,
  cells_to_take = NULL,
  streaming = FALSE,
  random_seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::qassert(tf_ids, "S+")
  assertScenicParams(scenic_params)
  checkmate::qassert(genes_to_take, c("S+", "0"))
  checkmate::qassert(cells_to_take, c("S+", "0"))
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(random_seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # resolve cells
  cell_indices <- if (is.null(cells_to_take)) {
    NULL
  } else {
    get_cell_indices(object, cell_ids = cells_to_take, rust_index = FALSE)
  }

  # resolve target genes
  if (is.null(genes_to_take)) {
    if (.verbose) {
      message("No target genes supplied, running gene filter...")
    }
    genes_to_take <- scenic_gene_filter_sc(
      object,
      scenic_params = scenic_params,
      cells_to_take = cells_to_take,
      .verbose = .verbose
    )
  }

  gene_idx <- get_gene_indices(
    object,
    gene_ids = genes_to_take,
    rust_index = FALSE
  )

  # resolve TFs
  all_gene_names <- S7::prop(object, "var_table")$gene_id
  tf_found <- tf_ids[tf_ids %in% all_gene_names]
  n_dropped <- length(tf_ids) - length(tf_found)
  if (n_dropped > 0 && .verbose) {
    warning(sprintf(
      "%d TF identifier(s) not found in the object and dropped.",
      n_dropped
    ))
  }

  if (length(tf_found) == 0) {
    stop("No provided TF identifiers match genes in the object.")
  }

  # TFs that survived the gene filter
  tf_in_targets <- intersect(tf_found, genes_to_take)

  if (length(tf_in_targets) == 0) {
    stop(
      "No TFs remain after intersecting with target gene indices. ",
      "Consider relaxing min_counts / min_cells thresholds."
    )
  }

  # subset matrix to filter-passing genes (and selected cells if any)
  sparse_data <- mc_counts_to_list(
    object,
    cell_indices = cell_indices,
    gene_indices = gene_idx,
    assay = "raw"
  )

  # 0-indexed positions within the filtered matrix
  gene_indices_rust <- seq_along(genes_to_take) - 1L
  tf_indices_rust <- match(tf_in_targets, genes_to_take) - 1L

  if (.verbose) {
    n_cells <- if (is.null(cell_indices)) {
      S7::prop(object, "dims")[1]
    } else {
      length(cell_indices)
    }
    message(sprintf(
      "SCENIC: %d target genes, %d TFs, %d cells",
      length(gene_indices_rust),
      length(tf_indices_rust),
      n_cells
    ))
  }

  importance_matrix <- rs_mc_scenic(
    sparse_data = sparse_data,
    tf_indices = as.integer(tf_indices_rust),
    scenic_params = scenic_params,
    seed = random_seed,
    verbose = parse_verbosity(.verbose)
  )

  rownames(importance_matrix) <- genes_to_take
  colnames(importance_matrix) <- tf_in_targets

  result <- new_scenic_grn(
    importance_matrix = importance_matrix,
    gene_ids = genes_to_take,
    tf_ids = tf_in_targets,
    params = scenic_params
  )

  return(result)
}
