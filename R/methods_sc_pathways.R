# methods that work with signatures and/or generate them -----------------------

## aucell ----------------------------------------------------------------------

#' Calculate AUC scores (akin to AUCell)
#'
#' @description
#' Calculates an AUC-type score akin to AUCell across the gene sets. You have
#' the options to calculate the AUC. Two options here: calculate this
#' with proper AUROC calculations (useful for marker gene expression, use the
#' `"auroc"` version) or based on the Mann-Whitney statistic (useful for pathway
#' activity measurs, use the `"wilcox"`). Data can be streamed in chunks of 50k
#' cells per or loaded in in one go.
#'
#' @param object `single_cell_exp` class.
#' @param gs_list Named list. The elements have the gene identifiers of the
#' respective gene sets.
#' @param auc_type String. Which type of AUC to calculate. Choice of
#' `c("wilcox", "auroc")`.
#' @param streaming Boolean. Shall the cell data be streamed in. Useful for
#' larger data sets.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return data.table with the DGE results from the test.
#'
#' @export
aucell_sc <- S7::new_generic(
  name = "aucell_sc",
  dispatch_args = "object",
  fun = function(
    object,
    gs_list,
    auc_type = c("wilcox", "auroc"),
    streaming = FALSE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method aucell_sc single_cell_exp
#'
#' @export
S7::method(aucell_sc, single_cell_exp) <- function(
  object,
  gs_list,
  auc_type = c("wilcox", "auroc"),
  streaming = FALSE,
  .verbose = TRUE
) {
  auc_type <- match.arg(auc_type)

  # checks
  checkmate::checkTRUE(S7::S7_inherits(object, single_cell_exp))
  checkmate::assertList(gs_list, types = "character", names = "named")
  checkmate::assertChoice(auc_type, c("wilcox", "auroc"))
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(.verbose, "B1")

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
    verbose = .verbose
  )

  colnames(auc_res) <- names(gs_list)
  rownames(auc_res) <- get_cell_names(object, filtered = TRUE)

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
#' names need to be part of the variables of the `single_cell_exp` class.
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
#' @param object `single_cell_exp` class.
#' @param gs_list Named nested list. The elements have the gene identifiers of
#' the respective gene sets and have the option to have a `"pos"` and `"neg"`
#' gene sets. The names need to be part of the variables of the
#' `single_cell_exp` class.
#' @param streaming Boolean. Shall the cell data be streamed in. Useful for
#' larger data sets.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return Matrix of cells x signatures with the VISION pathway scores as
#' values.
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
    streaming = FALSE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method vision_sc single_cell_exp
#'
#' @export
S7::method(vision_sc, single_cell_exp) <- function(
  object,
  gs_list,
  streaming = FALSE,
  .verbose = TRUE
) {
  # checks
  checkmate::checkTRUE(S7::S7_inherits(object, single_cell_exp))
  checkmate::assertList(gs_list, types = "list", names = "named")
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(.verbose, "B1")

  # no one sees this...
  vision_gs_clean <- purrr::map(gs_list, \(ls) {
    lapply(ls, FUN = get_gene_indices, x = object, rust_index = TRUE)
  })

  vision_res <- rs_vision(
    f_path = get_rust_count_cell_f_path(object),
    gs_list = vision_gs_clean,
    cells_to_keep = get_cells_to_keep(object),
    streaming = streaming,
    verbose = .verbose
  )

  colnames(vision_res) <- names(gs_list)
  rownames(vision_res) <- get_cell_names(object, filtered = TRUE)

  return(vision_res)
}

### vision with auto-correlation -----------------------------------------------

#' Calculate VISION scores (with auto-correlation scores)
#'
#' @description
#' Calculates an VISION-type scores for pathways based on DeTomaso, et al.
#' Compared to other score types, you can also calculate delta-type scores
#' between positive and negative gene indices, think epithelial vs mesenchymal
#' gene signature, etc. Additionally, this function also calculates the auto-
#' correlation values, answering the question if a given signature shows non-
#' random enrichment on the kNN graph. The kNN graph (and distance measures)
#' will be generated on-the-fly based on the embedding you wish to use.
#'
#' @param object `single_cell_exp` class.
#' @param gs_list Named nested list. The elements have the gene identifiers of
#' the respective gene sets and have the option to have a `"pos"` and `"neg"`
#' gene sets. The names need to be part of the variables of the
#' `single_cell_exp` class.
#' @param vision_params List with vision parameters, see
#' [bixverse::params_sc_vision()] with the following elements:
#' \itemize{
#'   \item n_perm - Integer. Number of random permutations
#'   \item n_cluster - Integer. Number of random clusters to generate to
#'   associate each set with.
#'   \item k - Number of neighbours for the kNN search. Only relevant if you
#'   set regenerate_knn to `TRUE`.
#'   \item knn_method - String. Which kNN algorithm to use. One of
#'   `c("annoy", "hnsw", "nndescent")`. Defaults to `"annoy"`. Only relevant if
#'   you set regenerate_knn to `TRUE`.
#'   \item ann_dist - String. Distance metric for the approximate neighbour
#'   search. One of `c("cosine", "euclidean")`. Defaults to `"cosine"`. Only
#'   relevant if you set regenerate_knn to `TRUE`.
#'   \item n_trees - Integer. Number of trees to use for the annoy algorithm.
#'   Only relevant if you set regenerate_knn to `TRUE`.
#'   \item search_budget - Integer. Search budget per tree for the annoy
#'   algorithm. Only relevant if you set regenerate_knn to `TRUE`.
#'   \item nn_max_iter - Integer. Maximum iterations for NN Descent. Only
#'   relevant if you set regenerate_knn to `TRUE` and use
#'   `"nndescent"`.
#'   \item rho - Numeric. Sampling rate for NN Descent. Only relevant if you
#'   set regenerate_knn to `TRUE` and use `"nndescent"`.
#'   \item delta - Numeric. Early termination criterion for NN Descent. Only
#'   relevant if you set regenerate_knn to `TRUE` and use `"nndescent"`.
#' }
#' @param embd_to_use String. The embedding to use. Whichever you chose, it
#' needs to be part of the object.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param random_seed Integer. The random seed.
#' @param streaming Boolean. Shall the cell data be streamed in. Useful for
#' larger data sets.
#' @param .verbose Boolean. Controls the verbosity of the function.
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
    vision_params = params_sc_vision(),
    streaming = FALSE,
    random_seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method vision_w_autocor_sc single_cell_exp
#'
#' @export
S7::method(vision_w_autocor_sc, single_cell_exp) <- function(
  object,
  gs_list,
  embd_to_use,
  no_embd_to_use = NULL,
  vision_params = params_sc_vision(),
  streaming = FALSE,
  random_seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  assertScVision(vision_params)
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(.verbose, "B1")

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

  vision_res <- rs_vision_with_autocorrelation(
    f_path = get_rust_count_cell_f_path(object),
    embd = embd,
    gs_list = vision_gs_clean,
    random_gs_list = random_gs_clean,
    vision_params = vision_params,
    cells_to_keep = get_cells_to_keep(object),
    cluster_membership = as.integer(cluster_membership),
    streaming = streaming,
    verbose = .verbose,
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
#' @param object `single_cell_exp` class.
#' @param embd_to_use String. The embedding to use. Defaults to `"pca"`.
#' @param hotspot_params List with vision parameters, see
#' [bixverse::params_sc_hotspot()] with the following elements:
#' \itemize{
#'   \item model - String. Which of the available models to use for the
#'   gene expression. Choices are one of `c("danb", "normal", "bernoulli")`.
#'   \item normalise - Boolean. Shall the data be normalised.
#'   \item k - Number of neighbours for the kNN search. Only relevant if you
#'   set regenerate_knn to `TRUE`.
#'   \item knn_method - String. Which kNN algorithm to use. One of
#'   `c("annoy", "hnsw", "nndescent")`. Defaults to `"annoy"`. Only relevant if
#'   you set regenerate_knn to `TRUE`.
#'   \item ann_dist - String. Distance metric for the approximate neighbour
#'   search. One of `c("cosine", "euclidean")`. Defaults to `"cosine"`. Only
#'   relevant if you set regenerate_knn to `TRUE`.
#'   \item n_trees - Integer. Number of trees to use for the annoy algorithm.
#'   Only relevant if you set regenerate_knn to `TRUE`.
#'   \item search_budget - Integer. Search budget per tree for the annoy
#'   algorithm. Only relevant if you set regenerate_knn to `TRUE`.
#'   \item nn_max_iter - Integer. Maximum iterations for NN Descent. Only
#'   relevant if you set regenerate_knn to `TRUE` and use
#'   `"nndescent"`.
#'   \item rho - Numeric. Sampling rate for NN Descent. Only relevant if you
#'   set regenerate_knn to `TRUE` and use `"nndescent"`.
#'   \item delta - Numeric. Early termination criterion for NN Descent. Only
#'   relevant if you set regenerate_knn to `TRUE` and use `"nndescent"`.
#' }
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param cells_to_take Optional string vector. If you want to only use
#' selected cells. If `NULL` will default to all cells_to_keep in the class.
#' @param genes_to_take Optional string vector. If you wish to limit the
#' search to a subset of genes. If `NULL` will default to all genes in the
#' class.
#' @param streaming Boolean. Shall the data be streamed in. Useful for larger
#' data sets.
#' @param random_seed Integer. Used for reproducibility.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A data.table with the auto-correlations on a per gene basis and
#' various statistics.
#'
#' @export
hotspot_autocor_sc <- S7::new_generic(
  name = "hotspot_autocor_sc",
  dispatch_args = "object",
  fun = function(
    object,
    embd_to_use = "pca",
    hotspot_params = params_sc_hotspot(),
    no_embd_to_use = NULL,
    cells_to_take = NULL,
    genes_to_take = NULL,
    streaming = FALSE,
    random_seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method hotspot_autocor_sc single_cell_exp
#'
#' @export
S7::method(hotspot_autocor_sc, single_cell_exp) <- function(
  object,
  embd_to_use = "pca",
  hotspot_params = params_sc_hotspot(),
  no_embd_to_use = NULL,
  cells_to_take = NULL,
  genes_to_take = NULL,
  streaming = FALSE,
  random_seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))
  checkmate::qassert(embd_to_use, "S1")
  assertScHotspot(hotspot_params)
  checkmate::qassert(no_embd_to_use, c("0", "I1"))
  checkmate::qassert(cells_to_take, c("S+", "0"))
  checkmate::qassert(genes_to_take, c("S+", "0"))
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(random_seed, "I1")
  checkmate::qassert(.verbose, "B1")

  # get the embedding
  checkmate::assertTRUE(embd_to_use %in% get_available_embeddings(object))
  embd <- get_embedding(x = object, embd_name = embd_to_use)

  if (!is.null(no_embd_to_use)) {
    to_take <- min(c(no_embd_to_use, ncol(embd)))
    embd <- embd[, 1:to_take]
  }

  if (is.null(cells_to_take)) {
    cells_to_take <- get_cell_names(object, filtered = TRUE)
  }

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

  hotspot_auto_cor <- rs_hotspot_autocor(
    f_path_genes = get_rust_count_gene_f_path(object),
    f_path_cells = get_rust_count_cell_f_path(object),
    embd = embd[cells_to_take, ],
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
    verbose = .verbose,
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
#' @param object `single_cell_exp` class.
#' @param embd_to_use String. The embedding to use. Defaults to `"pca"`.
#' @param hotspot_params List with vision parameters, see
#' [bixverse::params_sc_hotspot()] with the following elements:
#' \itemize{
#'   \item model - String. Which of the available models to use for the
#'   gene expression. Choices are one of `c("danb", "normal", "bernoulli")`.
#'   \item normalise - Boolean. Shall the data be normalised.
#'   \item k - Number of neighbours for the kNN search. Only relevant if you
#'   set regenerate_knn to `TRUE`.
#'   \item knn_method - String. Which kNN algorithm to use. One of
#'   `c("annoy", "hnsw", "nndescent")`. Defaults to `"annoy"`. Only relevant if
#'   you set regenerate_knn to `TRUE`.
#'   \item ann_dist - String. Distance metric for the approximate neighbour
#'   search. One of `c("cosine", "euclidean")`. Defaults to `"cosine"`. Only
#'   relevant if you set regenerate_knn to `TRUE`.
#'   \item n_trees - Integer. Number of trees to use for the annoy algorithm.
#'   Only relevant if you set regenerate_knn to `TRUE`.
#'   \item search_budget - Integer. Search budget per tree for the annoy
#'   algorithm. Only relevant if you set regenerate_knn to `TRUE`.
#'   \item nn_max_iter - Integer. Maximum iterations for NN Descent. Only
#'   relevant if you set regenerate_knn to `TRUE` and use
#'   `"nndescent"`.
#'   \item rho - Numeric. Sampling rate for NN Descent. Only relevant if you
#'   set regenerate_knn to `TRUE` and use `"nndescent"`.
#'   \item delta - Numeric. Early termination criterion for NN Descent. Only
#'   relevant if you set regenerate_knn to `TRUE` and use `"nndescent"`.
#' }
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param cells_to_take Optional string vector. If you want to only use
#' selected cells. If `NULL` will default to all cells_to_keep in the class.
#' @param genes_to_take Optional string vector. If you wish to limit the
#' search to a subset of genes. If `NULL` will default to all genes in the
#' class.
#' @param streaming Boolean. Shall the data be streamed in. Useful for larger
#' data sets.
#' @param random_seed Integer. Used for reproducibility.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A `sc_hotspot` class that can be used for subsequent analysis.
#'
#' @export
hotspot_gene_cor_sc <- S7::new_generic(
  name = "hotspot_gene_cor_sc",
  dispatch_args = "object",
  fun = function(
    object,
    embd_to_use = "pca",
    hotspot_params = params_sc_hotspot(),
    no_embd_to_use = NULL,
    cells_to_take = NULL,
    genes_to_take = NULL,
    streaming = FALSE,
    random_seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method hotspot_gene_cor_sc single_cell_exp
#'
#' @export
S7::method(hotspot_gene_cor_sc, single_cell_exp) <- function(
  object,
  embd_to_use = "pca",
  hotspot_params = params_sc_hotspot(),
  no_embd_to_use = NULL,
  cells_to_take = NULL,
  genes_to_take = NULL,
  streaming = FALSE,
  random_seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))
  checkmate::qassert(embd_to_use, "S1")
  assertScHotspot(hotspot_params)
  checkmate::qassert(no_embd_to_use, c("0", "I1"))
  checkmate::qassert(cells_to_take, c("S+", "0"))
  checkmate::qassert(genes_to_take, c("S+", "0"))
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(random_seed, "I1")
  checkmate::qassert(.verbose, "B1")

  # get the embedding
  checkmate::assertTRUE(embd_to_use %in% get_available_embeddings(object))
  embd <- get_embedding(x = object, embd_name = embd_to_use)

  if (!is.null(no_embd_to_use)) {
    to_take <- min(c(no_embd_to_use, ncol(embd)))
    embd <- embd[, 1:to_take]
  }

  if (is.null(cells_to_take)) {
    cells_to_take <- get_cell_names(object, filtered = TRUE)
  }

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

  hotspot_gene_cor <- rs_hotspot_gene_cor(
    f_path_genes = get_rust_count_gene_f_path(object),
    f_path_cells = get_rust_count_cell_f_path(object),
    embd = embd[cells_to_take, ],
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
    verbose = .verbose,
    seed = random_seed
  )

  hotspot_obj <- new_sc_hotspot_res(
    hotspot_res = hotspot_gene_cor,
    used_genes = genes_to_take,
    used_cells = cells_to_take
  )

  return(hotspot_obj)
}
