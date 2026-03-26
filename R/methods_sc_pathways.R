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
#' @param streaming Boolean. Shall the cell and gene data be streamed in.
#' Useful for larger data sets.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @returns description
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
    streaming = FALSE,
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
  streaming = FALSE,
  .verbose = TRUE
) {
  # checks
  checkmate::checkTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::assertList(gs_list, types = "character", names = "named")
  checkmate::qassert(n_bins, "I1")
  checkmate::qassert(n_ctrl, "I1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(.verbose, "B1")

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
    verbose = .verbose
  )

  colnames(module_res) <- names(gs_list)
  rownames(module_res) <- get_cell_names(object, filtered = TRUE)

  module_res
}

## aucell ----------------------------------------------------------------------

#' Calculate AUC scores (akin to AUCell)
#'
#' @description
#' Calculates an AUC-type score akin to AUCell across the gene sets, see Aibar
#' et al. You have the options to calculate the AUC. Two options here: calculate
#' this with proper AUROC calculations (useful for marker gene expression, use
#' the `"auroc"` version) or based on the Mann-Whitney statistic (useful for
#' pathway activity measurs, use the `"wilcox"`). Data can be streamed in chunks
#' of 50k cells per or loaded in in one go.
#'
#' @param object `SingleCells` class.
#' @param gs_list Named list. The elements have the gene identifiers of the
#' respective gene sets.
#' @param auc_type String. Which type of AUC to calculate. Choice of
#' `c("wilcox", "auroc")`.
#' @param streaming Boolean. Shall the cell data be streamed in. Useful for
#' larger data sets.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return AUCell results in form of a matrix that is cells x gene sets.
#'
#' @export
#'
#' @references Aibar, et al., Nat Methods, 2017
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

#' @method aucell_sc SingleCells
#'
#' @export
S7::method(aucell_sc, SingleCells) <- function(
  object,
  gs_list,
  auc_type = c("wilcox", "auroc"),
  streaming = FALSE,
  .verbose = TRUE
) {
  auc_type <- match.arg(auc_type)

  # checks
  checkmate::checkTRUE(S7::S7_inherits(object, SingleCells))
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

#' @method vision_sc SingleCells
#'
#' @export
S7::method(vision_sc, SingleCells) <- function(
  object,
  gs_list,
  streaming = FALSE,
  .verbose = TRUE
) {
  # checks
  checkmate::checkTRUE(S7::S7_inherits(object, SingleCells))
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

#' @method vision_w_autocor_sc SingleCells
#'
#' @export
S7::method(vision_w_autocor_sc, SingleCells) <- function(
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
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
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
#' @param object `SingleCells` class.
#' @param embd_to_use String. The embedding to use. Defaults to `"pca"`.
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
#' @param streaming Boolean. Shall the data be streamed in. Useful for larger
#' data sets.
#' @param random_seed Integer. Used for reproducibility.
#' @param .verbose Boolean. Controls verbosity of the function.
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

#' @method hotspot_autocor_sc SingleCells
#'
#' @export
S7::method(hotspot_autocor_sc, SingleCells) <- function(
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
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
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
#' @param object `SingleCells` class.
#' @param embd_to_use String. The embedding to use. Defaults to `"pca"`.
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
#' @param streaming Boolean. Shall the data be streamed in. Useful for larger
#' data sets.
#' @param random_seed Integer. Used for reproducibility.
#' @param .verbose Boolean. Controls verbosity of the function.
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

#' @method hotspot_gene_cor_sc SingleCells
#'
#' @export
S7::method(hotspot_gene_cor_sc, SingleCells) <- function(
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
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
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

## scenic ----------------------------------------------------------------------

### scenic gene filter ---------------------------------------------------------

#' Filter genes for SCENIC GRN inference
#'
#' @description
#' Filters genes by minimum total counts and minimum expressed-cell fraction
#' using the SCENIC inclusion criteria. Returns a character vector of gene
#' identifiers passing both filters.
#'
#' @param object `SingleCells` class.
#' @param scenic_params List. SCENIC parameters, see
#' [bixverse::params_scenic()]. Only `min_counts` and `min_cells` are used
#' by this function.
#' @param cells_to_take Optional string vector. Cell identifiers to restrict
#' to. If `NULL`, defaults to all filtered cells in the class.
#' @param .verbose Boolean. Controls verbosity. Defaults to `TRUE`.
#'
#' @returns A character vector of gene identifiers passing the SCENIC
#' inclusion criteria.
#'
#' @export
scenic_gene_filter_sc <- S7::new_generic(
  name = "scenic_gene_filter_sc",
  dispatch_args = "object",
  fun = function(
    object,
    scenic_params = params_scenic(),
    cells_to_take = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method scenic_gene_filter_sc SingleCells
#'
#' @export
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
  checkmate::qassert(.verbose, "B1")

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
    verbose = .verbose
  )

  passing_gene_ids <- get_gene_names_from_idx(
    object,
    gene_idx = as.integer(passing_indices),
    rust_based = TRUE
  )

  return(passing_gene_ids)
}

### scenic GRN inference -------------------------------------------------------

#' Run SCENIC GRN inference
#'
#' @description
#' Runs SCENIC GRN inference on the provided genes using the specified
#' transcription factors as predictors. Returns a `ScenicGrn` object
#' containing the TF-gene importance matrix for further processing.
#'
#' @param object `SingleCells` class.
#' @param tf_ids Character vector. Transcription factor gene identifiers to
#' use as predictors. Must be a subset of gene identifiers present in the
#' object.
#' @param scenic_params List. SCENIC parameters, see
#' [bixverse::params_scenic()].
#' @param genes_to_take Optional character vector. Target gene identifiers.
#' If `NULL`, genes are selected automatically via
#' [bixverse::scenic_gene_filter_sc()] using the `min_counts` and `min_cells`
#' thresholds in `scenic_params`.
#' @param cells_to_take Optional string vector. Cell identifiers to restrict
#' to. If `NULL`, defaults to all filtered cells in the class.
#' @param streaming Boolean. Whether to use the streaming implementation to
#' bound memory usage. Useful for large datasets. Defaults to `FALSE`.
#' @param random_seed Integer. Used for reproducibility. Defaults to `42L`.
#' @param .verbose Boolean. Controls verbosity. Defaults to `TRUE`.
#'
#' @returns A `ScenicGrn` object.
#'
#' @details
#' TF identifiers that are not found in the object's gene list are silently
#' dropped with a warning indicating how many were removed. TF indices are
#' intersected with the target gene indices so that TFs not passing the gene
#' filter are excluded from the predictor set but remain as potential targets
#' if present in `genes_to_take`. You have the option to generate the TF-gene
#' importance values with three distinct methods. For the `random_forest` and
#' the `extratrees` version, a batching strategy is applied in the default
#' settings. Correlated genes are identified and clustered together via
#' k-means clustering on the feature loadings of the PCA. These are then
#' divided into batches of `gene_batch_size` and the regression learners
#' are leveraging multi-target regression to fit all genes in the batch in one
#' go. This massively accelerates the algorithm and the importance values per
#' gene-TF pair are calculated then individually. Due to the batching by
#' similar gene, the signal dilution is limited. If you wish to run the
#' traditional approach, you can set gene_batch_size to `1L` or use the
#' `grnboost2` learner that can only fit one gene at a given time.
#'
#' @export
scenic_grn_sc <- S7::new_generic(
  name = "scenic_grn_sc",
  dispatch_args = "object",
  fun = function(
    object,
    tf_ids,
    scenic_params = params_scenic(),
    genes_to_take = NULL,
    cells_to_take = NULL,
    streaming = FALSE,
    random_seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method scenic_grn_sc SingleCells
#'
#' @export
S7::method(scenic_grn_sc, SingleCells) <- function(
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
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(tf_ids, "S+")
  assertScenicParams(scenic_params)
  checkmate::qassert(genes_to_take, c("S+", "0"))
  checkmate::qassert(cells_to_take, c("S+", "0"))
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(random_seed, "I1")
  checkmate::qassert(.verbose, "B1")

  # resolve cells
  if (is.null(cells_to_take)) {
    cells_to_take <- get_cell_names(object, filtered = TRUE)
  }

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
      .verbose = .verbose
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
    verbose = .verbose
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
