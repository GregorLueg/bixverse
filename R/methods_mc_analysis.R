# meta cell analysis methods ---------------------------------------------------

## aucell ----------------------------------------------------------------------

#' @method aucell_sc MetaCells
#'
#' @export
S7::method(aucell_sc, MetaCells) <- function(
  object,
  gs_list,
  aucell_params = params_sc_aucell(),
  streaming = NULL,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::assertList(gs_list, types = "character", names = "named")
  assertScAucell(aucell_params)
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # 0-indexed gene indices for Rust
  gs_indices <- purrr::map(gs_list, \(e) {
    get_gene_indices(x = object, gene_ids = e, rust_index = TRUE)
  })

  sparse_data <- mc_counts_to_list(object, assay = "raw")

  auc_res <- rs_mc_aucell(
    sparse_data = sparse_data,
    gs_list = gs_indices,
    aucell_params = aucell_params,
    verbose = parse_verbosity(.verbose)
  )

  colnames(auc_res) <- names(gs_list)
  rownames(auc_res) <- S7::prop(object, "obs_table")$meta_cell_id

  return(auc_res)
}

## vision ----------------------------------------------------------------------

### calculate the scores -------------------------------------------------------

# generics found in base_generics_sc.R

#' @method vision_sc MetaCells
#'
#' @export
S7::method(vision_sc, MetaCells) <- function(
  object,
  gs_list,
  streaming = NULL,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::assertList(gs_list, types = "list", names = "named")
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  vision_gs_clean <- purrr::map(gs_list, \(ls) {
    lapply(ls, FUN = get_gene_indices, x = object, rust_index = TRUE)
  })

  # VISION scores the normalised counts; `"raw"` here would score the raw ones
  # without complaining
  vision_res <- rs_mc_vision(
    sparse_data = mc_counts_to_list(object, assay = "norm"),
    gs_list = vision_gs_clean,
    verbose = parse_verbosity(.verbose)
  )

  colnames(vision_res) <- names(gs_list)
  rownames(vision_res) <- S7::prop(object, "obs_table")$meta_cell_id

  return(vision_res)
}

### vision with auto-correlation -----------------------------------------------

# generics found in base_generics_sc.R

#' @method vision_w_autocor_sc MetaCells
#'
#' @export
S7::method(vision_w_autocor_sc, MetaCells) <- function(
  object,
  gs_list,
  embd_to_use,
  no_embd_to_use = NULL,
  use_knn = TRUE,
  vision_params = params_sc_vision(),
  streaming = NULL,
  random_seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::assertList(gs_list, types = "list", names = "named")
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  checkmate::qassert(use_knn, "B1")
  assertScVision(vision_params)
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(random_seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # scores cover every meta cell and every gene, so only the embedding and the
  # kNN are taken off the prep
  prepped <- .prep_mc_embd_knn(
    object = object,
    embd_to_use = embd_to_use,
    use_knn = use_knn,
    no_embd_to_use = no_embd_to_use,
    cells_to_take = NULL,
    genes_to_take = NULL
  )

  vision_gs_clean <- purrr::map(gs_list, \(ls) {
    lapply(ls, FUN = get_gene_indices, x = object, rust_index = TRUE)
  })

  if (.verbose) {
    message(sprintf(
      "Generating %i random gene set clusters with a total of %s permutations.",
      vision_params$n_cluster,
      vision_params$n_perm
    ))
  }

  c(random_gs, cluster_membership) %<-%
    with(
      vision_params,
      .generate_null_perm_gs(
        gs_list = gs_list,
        expr_genes = S7::prop(object, "var_table")$gene_id,
        n_perm = n_perm,
        n_comp = n_cluster,
        random_seed = random_seed
      )
    )

  random_gs_clean <- purrr::map(random_gs, \(rs) {
    lapply(
      rs,
      FUN = function(ls) {
        lapply(ls, FUN = get_gene_indices, x = object, rust_index = TRUE)
      }
    )
  })

  vision_res <- rs_mc_vision_with_autocorrelation(
    sparse_data = mc_counts_to_list(object, assay = "norm"),
    embd = prepped$embd,
    knn_data = prepped$knn_data,
    gs_list = vision_gs_clean,
    random_gs_list = random_gs_clean,
    vision_params = vision_params,
    cluster_membership = as.integer(cluster_membership),
    verbose = parse_verbosity(.verbose),
    seed = random_seed
  )

  auto_cor_dt <- data.table::as.data.table(vision_res$autocor_res)[,
    gene_set_name := names(vision_gs_clean)
  ][, c("gene_set_name", "auto_cor", "p_val", "fdr"), with = FALSE]

  vision_matrix <- vision_res$vision_mat

  colnames(vision_matrix) <- names(gs_list)
  rownames(vision_matrix) <- S7::prop(object, "obs_table")$meta_cell_id

  result <- list(vision_matrix = vision_matrix, auto_cor_dt = auto_cor_dt)

  return(result)
}

## hotspot ---------------------------------------------------------------------

### auto-correlation -----------------------------------------------------------

# generics found in methods_sc_pathways.R

#' @method hotspot_autocor_sc MetaCells
#'
#' @export
S7::method(hotspot_autocor_sc, MetaCells) <- function(
  object,
  embd_to_use = "pca",
  use_knn = TRUE,
  hotspot_params = params_sc_hotspot(),
  no_embd_to_use = NULL,
  cells_to_take = NULL,
  genes_to_take = NULL,
  streaming = NULL,
  random_seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(use_knn, "B1")
  assertScHotspot(hotspot_params)
  checkmate::qassert(no_embd_to_use, c("0", "I1"))
  checkmate::qassert(cells_to_take, c("S+", "0"))
  checkmate::qassert(genes_to_take, c("S+", "0"))
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(random_seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  c(embd, cells_to_take, genes_to_take, knn_data) %<-%
    .prep_mc_embd_knn(
      object = object,
      embd_to_use = embd_to_use,
      use_knn = use_knn,
      no_embd_to_use = no_embd_to_use,
      cells_to_take = cells_to_take,
      genes_to_take = genes_to_take
    )

  hotspot_auto_cor <- rs_mc_hotspot_autocor(
    sparse_data = mc_counts_to_list(object, assay = "raw"),
    embd = embd,
    knn_data = knn_data,
    hotspot_params = hotspot_params,
    cells_to_keep = get_cell_indices(
      object,
      cell_ids = cells_to_take,
      rust_index = TRUE
    ),
    genes_to_use = get_gene_indices(
      object,
      gene_ids = genes_to_take,
      rust_index = TRUE
    ),
    verbose = parse_verbosity(.verbose),
    seed = random_seed
  )

  gene_ids <- S7::prop(object, "var_table")$gene_id[
    as.integer(hotspot_auto_cor$gene_idx) + 1L
  ]

  hotspot_auto_cor_dt <- data.table::as.data.table(hotspot_auto_cor)[,
    gene_id := gene_ids
  ][, c("gene_id", "gaerys_c", "z_score", "pval", "fdr")]

  return(hotspot_auto_cor_dt)
}

### gene <> gene local correlations --------------------------------------------

#' @method hotspot_gene_cor_sc MetaCells
#'
#' @export
S7::method(hotspot_gene_cor_sc, MetaCells) <- function(
  object,
  embd_to_use = "pca",
  use_knn = TRUE,
  hotspot_params = params_sc_hotspot(),
  no_embd_to_use = NULL,
  cells_to_take = NULL,
  genes_to_take = NULL,
  streaming = NULL,
  working_mem_gb = 4,
  random_seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(use_knn, "B1")
  assertScHotspot(hotspot_params)
  checkmate::qassert(no_embd_to_use, c("0", "I1"))
  checkmate::qassert(cells_to_take, c("S+", "0"))
  checkmate::qassert(genes_to_take, c("S+", "0"))
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(random_seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  c(embd, cells_to_take, genes_to_take, knn_data) %<-%
    .prep_mc_embd_knn(
      object = object,
      embd_to_use = embd_to_use,
      use_knn = use_knn,
      no_embd_to_use = no_embd_to_use,
      cells_to_take = cells_to_take,
      genes_to_take = genes_to_take
    )

  hotspot_gene_cor <- rs_mc_hotspot_gene_cor(
    sparse_data = mc_counts_to_list(object, assay = "raw"),
    embd = embd,
    knn_data = knn_data,
    hotspot_params = hotspot_params,
    cells_to_keep = get_cell_indices(
      object,
      cell_ids = cells_to_take,
      rust_index = TRUE
    ),
    genes_to_use = get_gene_indices(
      object,
      gene_ids = genes_to_take,
      rust_index = TRUE
    ),
    verbose = parse_verbosity(.verbose),
    seed = random_seed
  )

  hotspot_obj <- new_sc_hotspot_res(
    hotspot_res = hotspot_gene_cor,
    used_genes = genes_to_take,
    used_cells = cells_to_take
  )

  return(hotspot_obj)
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
  streaming = NULL,
  random_seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::qassert(tf_ids, "S+")
  assertScenicParams(scenic_params)
  checkmate::qassert(genes_to_take, c("S+", "0"))
  checkmate::qassert(cells_to_take, c("S+", "0"))
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(random_seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # resolve cells
  cell_indices <- if (is.null(cells_to_take)) {
    NULL
  } else {
    get_cell_indices(object, cell_ids = cells_to_take, rust_index = FALSE)
  }

  # reduce this here for meta cells...
  if (scenic_params$min_samples_leaf >= 20) {
    if (.verbose) {
      message(paste(
        "The mean leafs per sample is set quite high for meta cells.",
        "Reducing to 10L."
      ))
    }
    scenic_params$min_samples_leaf <- 10L
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

  all_gene_names <- S7::prop(object, "var_table")$gene_id

  # a supplied `genes_to_take` can name genes the object does not hold.
  # `get_gene_indices()` drops those silently, which would leave the row labels
  # of the importance matrix below one longer than the matrix itself
  genes_to_take <- .match_features(genes_to_take, all_gene_names)

  gene_idx <- get_gene_indices(
    object,
    gene_ids = genes_to_take,
    rust_index = FALSE
  )

  # resolve TFs
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

## nmf -------------------------------------------------------------------------

### helpers --------------------------------------------------------------------

#' Resolve cell and gene selection for NMF on MetaCells
#'
#' @param object `MetaCells` class.
#' @param cell_ids Optional string vector. The cells to include.
#' @param gene_ids Optional string vector. The genes to include.
#'
#' @returns A list with the following items:
#' \itemize{
#'  \item cell_indices - The Rust-based cell indices.
#'  \item gene_indices - The Rust-based gene indices.
#'  \item cell_ids - The cell identifiers (1-based)
#'  \item gene_ids - The gene identifiers (1-based).
#' }
#'
#' @keywords internal
#'
#' @keywords internal
.resolve_mc_nmf_selection <- function(object, cell_ids, gene_ids) {
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::qassert(cell_ids, c("S+", "0"))
  checkmate::qassert(gene_ids, c("S+", "0"))

  obs_table <- S7::prop(object, "obs_table")
  var_table <- S7::prop(object, "var_table")

  cell_indices_1b <- if (is.null(cell_ids)) {
    seq_len(nrow(obs_table))
  } else {
    get_cell_indices(object, cell_ids = cell_ids, rust_index = FALSE)
  }

  gene_indices_1b <- if (is.null(gene_ids)) {
    get_hvg(object)
  } else {
    get_gene_indices(object, gene_ids = gene_ids, rust_index = FALSE)
  }

  resolved_cell_ids <- obs_table$meta_cell_id[cell_indices_1b]
  resolved_gene_ids <- var_table$gene_id[gene_indices_1b]

  list(
    cell_indices_1b = as.integer(cell_indices_1b),
    gene_indices_1b = as.integer(gene_indices_1b),
    # Rust-side 0-based for storage on the result class
    cell_indices_rust = as.integer(cell_indices_1b - 1L),
    cell_ids = resolved_cell_ids,
    gene_ids = resolved_gene_ids
  )
}

### methods --------------------------------------------------------------------

#' @method nmf_sc MetaCells
#'
#' @export
S7::method(nmf_sc, MetaCells) <- function(
  object,
  k,
  cell_ids = NULL,
  gene_ids = NULL,
  preprocessing = "none",
  use_second_layer = TRUE,
  nmf_hals_params = params_nmf_hals(),
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::qassert(k, "I1[1,)")
  checkmate::qassert(cell_ids, c("0", "S+"))
  checkmate::qassert(gene_ids, c("0", "S+"))
  checkmate::assertChoice(preprocessing, c("none", "sd", "sqrt_sd"))
  checkmate::qassert(use_second_layer, "B1")
  assertNmfHals(nmf_hals_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  sel <- .resolve_mc_nmf_selection(object, cell_ids, gene_ids)

  assay <- if (use_second_layer) "norm" else "raw"

  count_list <- mc_counts_to_list(
    object = object,
    cell_indices = sel$cell_indices_1b,
    gene_indices = sel$gene_indices_1b,
    assay = assay
  )

  nmf_res <- rs_nmf_single_mc(
    sparse_data = count_list,
    k = k,
    preprocessing = preprocessing,
    use_second_layer = use_second_layer,
    nmf_hals_params = nmf_hals_params,
    seed = seed,
    verbose = parse_verbosity(.verbose)
  )

  params <- c(
    nmf_hals_params,
    list(
      k = k,
      preprocessing = preprocessing,
      use_second_layer = use_second_layer,
      seed = seed
    )
  )

  new_nmf_result(
    nmf_res = nmf_res,
    gene_ids = sel$gene_ids,
    cell_ids = sel$cell_ids,
    cell_indices = sel$cell_indices_rust,
    source_class = "MetaCells",
    params = params
  )
}


#' @method stabilised_nmf_sc MetaCells
#'
#' @export
S7::method(stabilised_nmf_sc, MetaCells) <- function(
  object,
  k,
  cell_ids = NULL,
  gene_ids = NULL,
  preprocessing = "none",
  use_second_layer = TRUE,
  nmf_hals_params = params_nmf_hals(),
  n_runs = 30L,
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::qassert(k, "I1[1,)")
  checkmate::qassert(cell_ids, c("0", "S+"))
  checkmate::qassert(gene_ids, c("0", "S+"))
  checkmate::assertChoice(preprocessing, c("none", "sd", "sqrt_sd"))
  checkmate::qassert(use_second_layer, "B1")
  assertNmfHals(nmf_hals_params)
  checkmate::qassert(n_runs, "I1[1,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  sel <- .resolve_mc_nmf_selection(object, cell_ids, gene_ids)

  assay <- if (use_second_layer) "norm" else "raw"

  count_list <- mc_counts_to_list(
    object = object,
    cell_indices = sel$cell_indices_1b,
    gene_indices = sel$gene_indices_1b,
    assay = assay
  )

  nmf_res <- rs_nmf_multi_mc(
    sparse_data = count_list,
    k = k,
    preprocessing = preprocessing,
    use_second_layer = use_second_layer,
    nmf_hals_params = nmf_hals_params,
    n_runs = n_runs,
    seed = seed,
    verbose = parse_verbosity(.verbose)
  )

  params <- c(
    nmf_hals_params,
    list(
      k = k,
      preprocessing = preprocessing,
      use_second_layer = use_second_layer,
      n_runs = n_runs,
      seed = seed
    )
  )

  new_stabilised_nmf_result(
    nmf_res = nmf_res,
    gene_ids = sel$gene_ids,
    cell_ids = sel$cell_ids,
    cell_indices = sel$cell_indices_rust,
    source_class = "MetaCells",
    params = params
  )
}

#' @method consensus_nmf_sc MetaCells
#'
#' @export
S7::method(consensus_nmf_sc, MetaCells) <- function(
  object,
  k,
  cell_ids = NULL,
  gene_ids = NULL,
  preprocessing = "none",
  use_second_layer = TRUE,
  nmf_hals_params = params_nmf_hals(),
  nmf_consensus_params = params_nmf_consensus(),
  n_runs = 30L,
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::qassert(k, "I1[2,)")
  checkmate::qassert(cell_ids, c("0", "S+"))
  checkmate::qassert(gene_ids, c("0", "S+"))
  checkmate::assertChoice(preprocessing, c("none", "sd", "sqrt_sd"))
  checkmate::qassert(use_second_layer, "B1")
  assertNmfHals(nmf_hals_params)
  assertNmfConsensus(nmf_consensus_params)
  checkmate::qassert(n_runs, "I1[2,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  sel <- .resolve_mc_nmf_selection(object, cell_ids, gene_ids)

  .warn_consensus_target_w(
    nmf_consensus_params,
    n_samples = length(sel$cell_indices_1b),
    k = k,
    n_runs = n_runs
  )

  assay <- if (use_second_layer) "norm" else "raw"

  count_list <- mc_counts_to_list(
    object = object,
    cell_indices = sel$cell_indices_1b,
    gene_indices = sel$gene_indices_1b,
    assay = assay
  )

  nmf_res <- .run_consensus_nmf(
    .rs_call = rs_nmf_consensus_mc,
    nmf_consensus_params = nmf_consensus_params,
    seed = seed,
    sparse_data = count_list,
    k = k,
    preprocessing = preprocessing,
    use_second_layer = use_second_layer,
    nmf_hals_params = nmf_hals_params,
    n_runs = n_runs,
    verbose = parse_verbosity(.verbose)
  )

  params <- c(
    nmf_hals_params,
    list(
      k = k,
      preprocessing = preprocessing,
      use_second_layer = use_second_layer,
      nmf_consensus_params = nmf_consensus_params,
      n_runs = n_runs,
      seed = seed
    )
  )

  new_consensus_nmf_result(
    nmf_res = nmf_res,
    gene_ids = sel$gene_ids,
    cell_ids = sel$cell_ids,
    cell_indices = sel$cell_indices_rust,
    source_class = "MetaCells",
    params = params
  )
}

#' @method nmf_k_sweep_sc MetaCells
#'
#' @export
S7::method(nmf_k_sweep_sc, MetaCells) <- function(
  object,
  k_range,
  cell_ids = NULL,
  gene_ids = NULL,
  preprocessing = "none",
  use_second_layer = TRUE,
  nmf_hals_params = params_nmf_hals(),
  nmf_consensus_params = params_nmf_consensus(),
  n_runs = 30L,
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  k_range <- .assert_nmf_k_range(k_range)
  checkmate::qassert(cell_ids, c("0", "S+"))
  checkmate::qassert(gene_ids, c("0", "S+"))
  checkmate::assertChoice(preprocessing, c("none", "sd", "sqrt_sd"))
  checkmate::qassert(use_second_layer, "B1")
  assertNmfHals(nmf_hals_params)
  assertNmfConsensus(nmf_consensus_params)
  checkmate::qassert(n_runs, "I1[2,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  sel <- .resolve_mc_nmf_selection(object, cell_ids, gene_ids)

  .warn_consensus_target_w(
    nmf_consensus_params,
    n_samples = length(sel$cell_indices_1b),
    k = max(k_range),
    n_runs = n_runs
  )

  assay <- if (use_second_layer) "norm" else "raw"

  count_list <- mc_counts_to_list(
    object = object,
    cell_indices = sel$cell_indices_1b,
    gene_indices = sel$gene_indices_1b,
    assay = assay
  )

  sweep_res <- rs_nmf_k_sweep_mc(
    sparse_data = count_list,
    k_range = k_range,
    preprocessing = preprocessing,
    use_second_layer = use_second_layer,
    nmf_hals_params = nmf_hals_params,
    nmf_consensus_params = .inject_consensus_seed(nmf_consensus_params, seed),
    n_runs = n_runs,
    seed = seed,
    verbose = parse_verbosity(.verbose)
  )

  new_nmf_k_sweep_result(
    sweep_res = sweep_res,
    source_class = "MetaCells",
    params = c(
      nmf_hals_params,
      list(
        k_range = k_range,
        preprocessing = preprocessing,
        use_second_layer = use_second_layer,
        nmf_consensus_params = nmf_consensus_params,
        n_runs = n_runs,
        seed = seed
      )
    )
  )
}
