# methods that work with signatures and/or generate them -----------------------

## module scoring --------------------------------------------------------------

#' Calculate module activity scores
#'
#' @description
#' Implements the approach from Tirosh et al into bixverse, i.e., the
#' AddModuleScore functionality from Seurat. For each module (gene set),
#' computes the average expression of genes in the set minus the average
#' expression of randomly selected control genes from the same expression bins.
#' Genes are binned based on their average expression across cells to ensure
#' controls are expression-matched.
#'
#' @param object `SingleCells` class.
#' @param gs_list Named list. The elements have the gene identifiers of the
#' respective gene sets.
#' @param n_bins Integer. Number of bins to use. Defaults to `24L`.
#' @param n_ctrl Integer. Number of control genes to use per gene for each gene
#' set. Defaults to `100L`.
#' @param seed Integer. The random seed.
#' @param streaming Optional Boolean. Shall the data be streamed in. Useful for
#' larger data sets where you wish to avoid loading in the whole data. If
#' `NULL`, will automatically detect.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns Returns a `ScMatrixRes` with the module scores.
#'
#' @references
#' Tirosh et al, Science (2016)
#'
#' @export
module_scores_sc <- S7::new_generic(
  name = "aucell_sc",
  dispatch_args = "object",
  fun = function(
    object,
    gs_list,
    n_bins = 24L,
    n_ctrl = 100L,
    seed = 42L,
    streaming = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method module_scores_sc SingleCells
#'
#' @export
S7::method(module_scores_sc, SingleCells) <- function(
  object,
  gs_list,
  n_bins = 24L,
  n_ctrl = 100L,
  seed = 42L,
  streaming = NULL,
  .verbose = TRUE
) {
  # checks
  checkmate::checkTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::assertList(gs_list, types = "character", names = "named")
  checkmate::qassert(n_bins, "I1")
  checkmate::qassert(n_ctrl, "I1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(.verbose, c("B1", "I1[0,1]"))

  streaming <- auto_streaming(
    n_cells = nrow(object),
    streaming = streaming,
    .verbose = .verbose
  )

  # get the gene indices
  gs_list <- purrr::map(gs_list, \(e) {
    get_gene_indices(x = object, gene_ids = e, rust_index = TRUE)
  })

  # let's get rusty...
  module_res <- rs_module_scoring(
    f_path_cells = get_rust_count_cell_f_path(object),
    f_path_genes = get_rust_count_gene_f_path(object),
    gs_list = gs_list,
    cells_to_keep = get_cells_to_keep(object),
    nbin = n_bins,
    ctrl = n_ctrl,
    seed = seed,
    streaming = streaming,
    verbose = parse_verbosity(.verbose)
  )

  colnames(module_res) <- names(gs_list)
  rownames(module_res) <- get_cell_names(object, filtered = TRUE)

  module_res <- new_sc_matrix(
    res = module_res,
    cell_indices = get_cells_to_keep(object)
  )

  return(module_res)
}

## aucell ----------------------------------------------------------------------

#' @method aucell_sc SingleCells
#'
#' @export
S7::method(aucell_sc, SingleCells) <- function(
  object,
  gs_list,
  auc_type = c("wilcox", "auroc"),
  streaming = NULL,
  .verbose = TRUE
) {
  auc_type <- match.arg(auc_type)

  streaming <- auto_streaming(
    n_cells = nrow(object),
    streaming = streaming,
    .verbose = .verbose
  )

  # checks
  checkmate::checkTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::assertList(gs_list, types = "character", names = "named")
  checkmate::assertChoice(auc_type, c("wilcox", "auroc"))
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # get the gene indices
  gs_list <- purrr::map(gs_list, \(e) {
    get_gene_indices(x = object, gene_ids = e, rust_index = TRUE)
  })

  auc_res <- rs_aucell(
    f_path = get_rust_count_cell_f_path(object),
    gs_list = gs_list,
    cells_to_keep = get_cells_to_keep(object),
    auc_type = auc_type,
    streaming = streaming,
    verbose = parse_verbosity(.verbose)
  )

  colnames(auc_res) <- names(gs_list)
  rownames(auc_res) <- get_cell_names(object, filtered = TRUE)

  auc_res <- new_sc_matrix(
    res = auc_res,
    cell_indices = get_cells_to_keep(object)
  )

  return(auc_res)
}

## vision ----------------------------------------------------------------------

### helper ---------------------------------------------------------------------

#' Generate random gene sets for VISION p-value calculations
#'
#' @description
#' This function will generate random gene sets given the provided gs_list.
#' Under the hood, it uses the same approach as in DeTomaso, et al. and does
#' not generate a random signature per given signature, but `n_comp`
#' representative ones based on size and balance (positive and negative genes)
#' via k-means clustering. The original gene sets are than matched to the
#' closest cluster. The authors observed that this sufficed to estimate
#' significance, see DeTomaso, et al.
#'
#' @param gs_list Named nested list for which to calculate the local
#' auto-correlations. The elements have the gene identifiers of the respective
#' gene sets and have the option to have a `"pos"` and `"neg"` gene sets. The
#' names need to be part of the variables of the `SingleCells` class.
#' @param expr_genes Character vector. Represents the genes expressed in the
#' experiment.
#' @param n_perm Integer. Number of random permutations to generate.
#' @param n_comp Integer. Number of k-means cluster to identify.
#' @param random_seed Integer. For reproducibility purposes.
#' @param no_cores Optional integer. Number of sessions to use for the
#' [mirai::mirai_map()] approach during generation of the random gene sets.
#' If not provided, will default to half of the available cores with a maximum
#' of `8L`.
#'
#' @returns A list with the following elements:
#' \itemize{
#'   \item random_signatures - Nested list representing the random permutations.
#'   \item clusters - Association of original gene set to random permutation
#'   set.
#' }
#'
#' @keywords internal
generate_null_perm_gs <- function(
  gs_list,
  expr_genes,
  n_perm = 500L,
  n_comp = 5L,
  random_seed = 42L,
  no_cores = NULL
) {
  # checks
  checkmate::assertList(gs_list, types = "list", names = "named")
  checkmate::qassert(expr_genes, "S+")
  checkmate::qassert(n_perm, "I1")
  checkmate::qassert(n_comp, "I1")
  checkmate::qassert(random_seed, "I1")
  checkmate::qassert(no_cores, c("0", "I1"))

  if (is.null(no_cores)) {
    no_cores <- get_cores()
  }

  # function
  sig_data_signed <- purrr::map(gs_list, \(gs) {
    all_genes <- c(gs$pos, gs$neg)
    signs <- c(rep(1, length(gs$pos)), rep(-1, length(gs$neg)))
    signs
  })

  sig_sizes <- purrr::map_dbl(sig_data_signed, length)
  sig_sizes <- log10(sig_sizes)

  sig_balance <- purrr::map_dbl(sig_data_signed, \(sig) {
    sum(sig == 1) / length(sig)
  })

  sig_vars <- cbind(sig_sizes, sig_balance)

  if (nrow(sig_vars) <= n_comp) {
    n_comp <- nrow(sig_vars)
    centers <- sig_vars
    clusters <- as.factor(seq_len(nrow(sig_vars)))
    names(clusters) <- rownames(sig_vars)
  } else {
    if (nrow(unique(sig_vars)) <= n_comp) {
      n_comp <- nrow(unique(sig_vars))
    }

    km <- kmeans(sig_vars, n_comp)
    centers <- km$centers

    levels <- as.character(seq(n_comp))
    clusters <- factor(km$cluster, levels = levels)
  }

  centers[, "sig_sizes"] <- round(10**centers[, "sig_sizes"])

  mirai::daemons(no_cores)

  random_sigs <- mirai::mirai_map(
    seq_len(nrow(centers)),
    .f = function(i, n_perm, centers, random_seed, gene_names) {
      size <- centers[i, "sig_sizes"]
      balance <- centers[i, "sig_balance"]
      n_pos_genes <- ceiling(balance * size)
      n_neg_genes <- size - n_pos_genes

      lapply(seq_len(n_perm), function(iter) {
        set.seed(random_seed + (i - 1) * 1e6 + iter)
        genes <- sample(gene_names, size)

        if (n_pos_genes == size) {
          list(pos = genes)
        } else if (n_neg_genes == size) {
          list(neg = genes)
        } else {
          list(
            pos = genes[1:n_pos_genes],
            neg = genes[(n_pos_genes + 1):size]
          )
        }
      })
    },
    .args = list(
      n_perm = n_perm,
      centers = centers,
      random_seed = random_seed,
      gene_names = expr_genes
    )
  )[]

  mirai::daemons(0)

  res <- list(
    random_signatures = random_sigs,
    clusters = clusters
  )

  return(res)
}

### calculate the scores -------------------------------------------------------

#' Calculate VISION scores
#'
#' @description
#' Calculates an VISION-type scores for pathways based on DeTomaso, et al.
#' Compared to other score types, you can also calculate delta-type scores
#' between positive and negative gene indices, think epithelial vs mesenchymal
#' gene signature, etc.
#'
#' @param object `SingleCells` class.
#' @param gs_list Named nested list. The elements have the gene identifiers of
#' the respective gene sets and have the option to have a `"pos"` and `"neg"`
#' gene sets. The names need to be part of the variables of the
#' `SingleCells` class.
#' @param streaming Optional Boolean. Shall the data be streamed in. Useful for
#' larger data sets where you wish to avoid loading in the whole data. If
#' `NULL`, will automatically detect.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns Returns a `ScMatrixRes` with the VISION scores.
#'
#' @references DeTomaso, et al., Nat. Commun., 2019
#'
#' @export
vision_sc <- S7::new_generic(
  name = "vision_sc",
  dispatch_args = "object",
  fun = function(
    object,
    gs_list,
    streaming = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method vision_sc SingleCells
#'
#' @export
S7::method(vision_sc, SingleCells) <- function(
  object,
  gs_list,
  streaming = NULL,
  .verbose = TRUE
) {
  # checks
  checkmate::checkTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::assertList(gs_list, types = "list", names = "named")
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  streaming <- auto_streaming(
    n_cells = nrow(object),
    streaming = streaming,
    .verbose = .verbose
  )

  # no one sees this...
  vision_gs_clean <- purrr::map(gs_list, \(ls) {
    lapply(ls, FUN = get_gene_indices, x = object, rust_index = TRUE)
  })

  vision_res <- rs_vision(
    f_path = get_rust_count_cell_f_path(object),
    gs_list = vision_gs_clean,
    cells_to_keep = get_cells_to_keep(object),
    streaming = streaming,
    verbose = parse_verbosity(.verbose)
  )

  colnames(vision_res) <- names(gs_list)
  rownames(vision_res) <- get_cell_names(object, filtered = TRUE)

  vision_res <- new_sc_matrix(
    res = vision_res,
    cell_indices = get_cells_to_keep(object)
  )

  return(vision_res)
}

### vision with auto-correlation -----------------------------------------------

#' Calculate VISION scores (with auto-correlation scores)
#'
#' @description
#' Calculates VISION-type scores for pathways based on DeTomaso, et al.
#' Compared to other score types, you can also calculate delta-type scores
#' between positive and negative gene indices, think epithelial vs mesenchymal
#' gene signature, etc. Additionally, this function also calculates the auto-
#' correlation values, answering the question if a given signature shows non-
#' random enrichment on the kNN graph. The kNN graph (and distance measures)
#' will be generated on-the-fly based on the embedding you wish to use.
#'
#' @param object `SingleCells` class.
#' @param gs_list Named nested list. The elements have the gene identifiers of
#' the respective gene sets and have the option to have a `"pos"` and `"neg"`
#' gene sets. The names need to be part of the variables of the
#' `SingleCells` class.
#' @param vision_params List with vision parameters, see
#' [bixverse::params_sc_vision()] with the following elements:
#' \itemize{
#'   \item n_perm - Integer. Number of random permutations
#'   \item n_cluster - Integer. Number of random clusters to generate to
#'   associate each set with.
#'   \item knn - List of kNN parameters. See [bixverse::params_knn_defaults()]
#'   for available parameters and their defaults.
#' }
#' @param embd_to_use String. The embedding to use. Whichever you chose, it
#' needs to be part of the object.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param use_knn Boolean. Shall the internal kNN be used. If set to yes, you
#' need to ensure consistency.
#' @param random_seed Integer. The random seed.
#' @param streaming Optional Boolean. Shall the data be streamed in. Useful for
#' larger data sets where you wish to avoid loading in the whole data. If
#' `NULL`, will automatically detect.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @return Matrix of cells x signatures with the VISION pathway scores as
#' values.
#'
#' @references DeTomaso, et al., Nat. Commun., 2019
#'
#' @export
vision_w_autocor_sc <- S7::new_generic(
  name = "vision_w_autocor_sc",
  dispatch_args = "object",
  fun = function(
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
    S7::S7_dispatch()
  }
)

#' @method vision_w_autocor_sc SingleCells
#'
#' @export
S7::method(vision_w_autocor_sc, SingleCells) <- function(
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
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  assertScVision(vision_params)
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  streaming <- auto_streaming(
    n_cells = nrow(object),
    streaming = streaming,
    .verbose = .verbose
  )

  # get the embedding
  checkmate::assertTRUE(embd_to_use %in% get_available_embeddings(object))
  embd <- get_embedding(x = object, embd_name = embd_to_use)

  if (!is.null(no_embd_to_use)) {
    to_take <- min(c(no_embd_to_use, ncol(embd)))
    embd <- embd[, 1:to_take]
  }

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
      generate_null_perm_gs(
        gs_list = gs_list,
        expr_genes = get_gene_names(object),
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

  knn_data <- if (use_knn) {
    get_knn_obj(object)
  } else {
    NULL
  }

  vision_res <- rs_vision_with_autocorrelation(
    f_path = get_rust_count_cell_f_path(object),
    embd = embd,
    knn_data = knn_data,
    gs_list = vision_gs_clean,
    random_gs_list = random_gs_clean,
    vision_params = vision_params,
    cells_to_keep = get_cells_to_keep(object),
    cluster_membership = as.integer(cluster_membership),
    streaming = streaming,
    verbose = parse_verbosity(.verbose),
    seed = random_seed
  )

  auto_cor_dt <- data.table::as.data.table(vision_res$autocor_res)[,
    gene_set_name := names(vision_gs_clean)
  ][, c("gene_set_name", "auto_cor", "p_val", "fdr"), with = FALSE]

  vision_matrix <- vision_res$vision_mat

  colnames(vision_matrix) <- names(gs_list)
  rownames(vision_matrix) <- get_cell_names(object, filtered = TRUE)

  result <- list(vision_matrix = vision_matrix, auto_cor_dt = auto_cor_dt)

  return(result)
}

## hotspot ---------------------------------------------------------------------

### gene auto-correlations -----------------------------------------------------

#' Calculate the local auto-correlation of a gene
#'
#' @description
#' This method implements the HotSpot approach (see DeTomaso, et al.) to
#' calculate the auto-correlation of a given gene in the kNN graph based on
#' the chosen embedding. This can be used to identify genes that have strong
#' local correlations and vary across the kNN graph.
#'
#' @param object `SingleCells` class.
#' @param embd_to_use String. The embedding to use. Defaults to `"pca"`.
#' @param use_knn Boolean. Shall the internal kNN be used. If set to yes, you
#' need to ensure consistency. If you provide `cells_to_take`, the function
#' will regenerate the kNN graph with these cells.
#' @param hotspot_params List with hotspot parameters, see
#' [bixverse::params_sc_hotspot()] with the following elements:
#' \itemize{
#'   \item model - String. Which of the available models to use for the
#'   gene expression. Choices are one of `c("danb", "normal", "bernoulli")`.
#'   \item normalise - Boolean. Shall the data be normalised.
#'   \item knn - List of kNN parameters. See [bixverse::params_knn_defaults()]
#'   for available parameters and their defaults.
#' }
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param cells_to_take Optional string vector. If you want to only use
#' selected cells. If `NULL` will default to all cells_to_keep in the class.
#' @param genes_to_take Optional string vector. If you wish to limit the
#' search to a subset of genes. If `NULL` will default to all genes in the
#' class.
#' @param streaming Optional Boolean. Shall the data be streamed in. Useful for
#' larger data sets where you wish to avoid loading in the whole data. If
#' `NULL`, will automatically detect.
#' @param random_seed Integer. Used for reproducibility.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns A data.table with the auto-correlations on a per gene basis and
#' various statistics.
#'
#' @details
#' Should a gene not be found in sufficient cells, the gene will be
#' automatically filtered out from the results. This can occur for example
#' if you have filtered out the cells that contain a given gene. The underlying
#' genes are still available, but the cells that might contain them are not
#' included.
#'
#' @export
hotspot_autocor_sc <- S7::new_generic(
  name = "hotspot_autocor_sc",
  dispatch_args = "object",
  fun = function(
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
    S7::S7_dispatch()
  }
)

#' @method hotspot_autocor_sc SingleCells
#'
#' @export
S7::method(hotspot_autocor_sc, SingleCells) <- function(
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
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(embd_to_use, "S1")
  assertScHotspot(hotspot_params)
  checkmate::qassert(no_embd_to_use, c("0", "I1"))
  checkmate::qassert(cells_to_take, c("S+", "0"))
  checkmate::qassert(genes_to_take, c("S+", "0"))
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(random_seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # get the embedding
  checkmate::assertTRUE(embd_to_use %in% get_available_embeddings(object))
  embd <- get_embedding(x = object, embd_name = embd_to_use)

  if (!is.null(no_embd_to_use)) {
    to_take <- min(c(no_embd_to_use, ncol(embd)))
    embd <- embd[, 1:to_take]
  }

  if (is.null(cells_to_take)) {
    cells_to_take <- get_cell_names(object, filtered = TRUE)
  } else {
    # if the user overwrites this specifically, recreate the kNN data internally
    use_knn <- FALSE
  }

  streaming <- auto_streaming(
    n_cells = length(cells_to_take),
    streaming = streaming,
    .verbose = .verbose
  )

  if (is.null(genes_to_take)) {
    genes_to_take <- get_gene_names(object)
  }

  if (length(intersect(rownames(embd), cells_to_take)) == 0) {
    stop(
      paste(
        "There is no intersection in the provided cell names and embedding",
        "space. Please double check your parameters!"
      )
    )
  }

  knn_data <- if (use_knn) {
    get_knn_obj(object)
  } else {
    NULL
  }

  hotspot_auto_cor <- rs_hotspot_autocor(
    f_path_genes = get_rust_count_gene_f_path(object),
    f_path_cells = get_rust_count_cell_f_path(object),
    embd = embd[cells_to_take, ],
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
    streaming = streaming,
    verbose = parse_verbosity(.verbose),
    seed = random_seed
  )

  hotspot_auto_cor_dt <- data.table::as.data.table(hotspot_auto_cor)[,
    gene_id := get_gene_names_from_idx(
      object,
      gene_idx = as.integer(hotspot_auto_cor$gene_idx),
      rust_based = TRUE
    )
  ][, c("gene_id", "gaerys_c", "z_score", "pval", "fdr")]

  return(hotspot_auto_cor_dt)
}

### gene <> gene local correlations --------------------------------------------

#' Calculate the local pairwise gene-gene correlation
#'
#' @description
#' This method implements the HotSpot approach (see DeTomaso, et al.) to
#' calculate the local gene-gene correlations and their Z-scores.
#'
#' @param object `SingleCells` class.
#' @param embd_to_use String. The embedding to use. Defaults to `"pca"`.
#' @param use_knn Boolean. Shall the internal kNN be used. If set to yes, you
#' need to ensure consistency. If you provide `cells_to_take`, the function
#' will regenerate the kNN graph with these cells.
#' @param hotspot_params List with hotspot parameters, see
#' [bixverse::params_sc_hotspot()] with the following elements:
#' \itemize{
#'   \item model - String. Which of the available models to use for the
#'   gene expression. Choices are one of `c("danb", "normal", "bernoulli")`.
#'   \item normalise - Boolean. Shall the data be normalised.
#'   \item knn - List of kNN parameters. See [bixverse::params_knn_defaults()]
#'   for available parameters and their defaults.
#' }
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param cells_to_take Optional string vector. If you want to only use
#' selected cells. If `NULL` will default to all cells_to_keep in the class.
#' @param genes_to_take Optional string vector. If you wish to limit the
#' search to a subset of genes. If `NULL` will default to all genes in the
#' class.
#' @param streaming Optional Boolean. Shall the data be streamed in. Useful for
#' larger data sets where you wish to avoid loading in the whole data. If
#' `NULL`, will automatically detect.
#' @param working_mem_gb Numeric. Approximate working memory (GB) the streaming
#'  pair path may use for resident gene panels. Ignored when `streaming` is
#'  `FALSE`. Larger values mean fewer disk re-reads. Note this excludes the two
#' dense N_genes x N_genes output matrices, which scale with `genes_to_use`.
#' Defaults to `4` (4 GB of memory allocated).
#' @param random_seed Integer. Used for reproducibility.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns A `sc_hotspot` class that can be used for subsequent analysis.
#'
#' @details
#' Should a gene not be found in sufficient cells, the pairs with this gene
#' will be set to 0. Please ensure prior to running the function that you
#' are only calculating gene-gene auto-correlations that occur in sufficient
#' cells.
#'
#' @export
hotspot_gene_cor_sc <- S7::new_generic(
  name = "hotspot_gene_cor_sc",
  dispatch_args = "object",
  fun = function(
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
    S7::S7_dispatch()
  }
)

#' @method hotspot_gene_cor_sc SingleCells
#'
#' @export
S7::method(hotspot_gene_cor_sc, SingleCells) <- function(
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
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(embd_to_use, "S1")
  assertScHotspot(hotspot_params)
  checkmate::qassert(no_embd_to_use, c("0", "I1"))
  checkmate::qassert(cells_to_take, c("S+", "0"))
  checkmate::qassert(genes_to_take, c("S+", "0"))
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(random_seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # get the embedding
  checkmate::assertTRUE(embd_to_use %in% get_available_embeddings(object))
  embd <- get_embedding(x = object, embd_name = embd_to_use)

  if (!is.null(no_embd_to_use)) {
    to_take <- min(c(no_embd_to_use, ncol(embd)))
    embd <- embd[, 1:to_take]
  }

  if (is.null(cells_to_take)) {
    cells_to_take <- get_cell_names(object, filtered = TRUE)
  } else {
    # if the user overwrites this specifically, recreate the kNN data internally
    use_knn <- FALSE
  }

  streaming <- auto_streaming(
    n_cells = length(cells_to_take),
    streaming = streaming,
    .verbose = .verbose
  )

  if (is.null(genes_to_take)) {
    genes_to_take <- get_gene_names(object)
  }

  if (length(intersect(rownames(embd), cells_to_take)) == 0) {
    stop(
      paste(
        "There is no intersection in the provided cell names and embedding",
        "space. Please double check your parameters!"
      )
    )
  }

  knn_data <- if (use_knn) {
    get_knn_obj(object)
  } else {
    NULL
  }

  hotspot_gene_cor <- rs_hotspot_gene_cor(
    f_path_genes = get_rust_count_gene_f_path(object),
    f_path_cells = get_rust_count_cell_f_path(object),
    embd = embd[cells_to_take, ],
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
    streaming = streaming,
    working_mem_gb = working_mem_gb,
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

### scenic gene filter ---------------------------------------------------------

# generics found in base_generics_sc.R

S7::method(scenic_gene_filter_sc, SingleCells) <- function(
  object,
  scenic_params = params_scenic(),
  cells_to_take = NULL,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScenicParams(scenic_params)
  checkmate::qassert(cells_to_take, c("S+", "0"))
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  if (is.null(cells_to_take)) {
    cells_to_take <- get_cell_names(object, filtered = TRUE)
  }

  cell_indices <- get_cell_indices(
    object,
    cell_ids = cells_to_take,
    rust_index = TRUE
  )

  passing_indices <- rs_scenic_gene_filter(
    f_path_genes = get_rust_count_gene_f_path(object),
    cell_indices = cell_indices,
    scenic_params = scenic_params,
    verbose = parse_verbosity(.verbose)
  )

  passing_gene_ids <- get_gene_names_from_idx(
    object,
    gene_idx = as.integer(passing_indices),
    rust_based = TRUE
  )

  return(passing_gene_ids)
}

### scenic GRN inference -------------------------------------------------------

#' @method scenic_grn_sc SingleCells
#'
#' @export
S7::method(scenic_grn_sc, SingleCells) <- function(
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
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(tf_ids, "S+")
  assertScenicParams(scenic_params)
  checkmate::qassert(genes_to_take, c("S+", "0"))
  checkmate::qassert(cells_to_take, c("S+", "0"))
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(random_seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # resolve cells
  if (is.null(cells_to_take)) {
    cells_to_take <- get_cell_names(object, filtered = TRUE)
  }

  streaming <- auto_streaming(
    n_cells = length(cells_to_take),
    streaming = streaming,
    .verbose = .verbose
  )

  cell_indices <- get_cell_indices(
    object,
    cell_ids = cells_to_take,
    rust_index = TRUE
  )

  # resolve target genes
  if (is.null(genes_to_take)) {
    if (.verbose) {
      message("No target genes supplied, running gene filter...")
    }
    genes_to_take <- scenic_gene_filter_sc(
      object,
      scenic_params = scenic_params,
      cells_to_take = cells_to_take,
      .verbose = parse_verbosity(.verbose)
    )
  }

  gene_indices <- get_gene_indices(
    object,
    gene_ids = genes_to_take,
    rust_index = TRUE
  )

  # resolve TF indices, dropping any not found in the object
  all_gene_names <- get_gene_names(object)
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

  tf_all_indices <- get_gene_indices(
    object,
    gene_ids = tf_found,
    rust_index = TRUE
  )

  # intersect TFs with target genes so only TFs passing the gene filter
  # are used as predictors
  tf_indices_red <- intersect(tf_all_indices, gene_indices)

  if (length(tf_indices_red) == 0) {
    stop(
      "No TFs remain after intersecting with target gene indices. ",
      "Consider relaxing min_counts / min_cells thresholds.",
      "Returning NULL"
    )
    return(NULL)
  }

  if (.verbose) {
    message(sprintf(
      "SCENIC: %d target genes, %d TFs, %d cells",
      length(gene_indices),
      length(tf_indices_red),
      length(cell_indices)
    ))
  }

  # run inference
  scenic_fn <- if (streaming) rs_scenic_grn_streaming else rs_scenic_grn
  importance_matrix <- scenic_fn(
    f_path_genes = get_rust_count_gene_f_path(object),
    cell_indices = cell_indices,
    gene_indices = gene_indices,
    tf_indices = as.integer(tf_indices_red),
    scenic_params = scenic_params,
    seed = random_seed,
    verbose = parse_verbosity(.verbose)
  )

  # label the matrix
  tf_names <- get_gene_names_from_idx(
    object,
    gene_idx = as.integer(tf_indices_red),
    rust_based = TRUE
  )
  rownames(importance_matrix) <- genes_to_take
  colnames(importance_matrix) <- tf_names

  # wrap in results class
  result <- new_scenic_grn(
    importance_matrix = importance_matrix,
    gene_ids = genes_to_take,
    tf_ids = tf_names,
    params = scenic_params
  )

  return(result)
}
