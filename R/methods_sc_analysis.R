# single cell analysis methods -------------------------------------------------

## dges ------------------------------------------------------------------------

### find markers ---------------------------------------------------------------

#' Calculate DGE between two cell groups
#'
#' @description
#' This function can be used to calculate differentially expressed genes
#' between two groups in the single cell data. At the moment, it has only
#' an implementation for the Wilcox-based rank statistic.
#'
#' @param object `single_cell_exp` class.
#' @param cells_1 String. The names of the cells in group 1. Need to be part
#' of the cell names in the object, see [bixverse::get_cell_names()].
#' @param cells_2 String. The names of the cells in group 2. Need to be part
#' of the cell names in the object, see [bixverse::get_cell_names()].
#' @param method String. Which method to use for the calculations of the DGE.
#' At the moment the only option is `"wilcox"`, but the parameter is reserved
#' for future features.
#' @param alternative String. Test alternative. One of
#' `c("twosided", "greater", "less")`. Function will default to `"twosided"`.
#' @param min_prop Numeric. The minimum proportion of cells that need to express
#' the gene to be tested in any of the two groups.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return data.table with the DGE results from the test.
#'
#' @export
find_markers_sc <- S7::new_generic(
  name = "find_markers_sc",
  dispatch_args = "object",
  fun = function(
    object,
    cells_1,
    cells_2,
    method = c("wilcox"),
    alternative = c("twosided", "greater", "less"),
    min_prop = 0.05,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method find_markers_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(find_markers_sc, single_cell_exp) <- function(
  object,
  cells_1,
  cells_2,
  method = c("wilcox"),
  alternative = c("twosided", "greater", "less"),
  min_prop = 0.05,
  .verbose = TRUE
) {
  alternative <- match.arg(alternative)

  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  assertCellsExist(object, cells_1)
  assertCellsExist(object, cells_2)
  checkmate::assertChoice(method, c("wilcox"))
  checkmate::assertChoice(alternative, c("twosided", "greater", "less"))
  checkmate::qassert(min_prop, "N1[0, 1]")
  checkmate::qassert(.verbose, "B1")

  dge_results <- switch(
    method,
    "wilcox" = rs_calculate_dge_mann_whitney(
      f_path = get_rust_count_cell_f_path(object),
      cell_indices_1 = get_cell_indices(
        x = object,
        cell_ids = cells_1,
        rust_index = TRUE
      ),
      cell_indices_2 = get_cell_indices(
        x = object,
        cell_ids = cells_2,
        rust_index = TRUE
      ),
      min_prop = min_prop,
      alternative = alternative,
      verbose = .verbose
    )
  )

  dge_dt <- data.table::as.data.table(dge_results[c(
    "lfc",
    "prop1",
    "prop2",
    "z_scores",
    "p_values",
    "fdr"
  )])[,
    gene_id := get_gene_names(object)[dge_results$genes_to_keep]
  ]

  data.table::setcolorder(
    dge_dt,
    c(
      "gene_id",
      "lfc",
      "prop1",
      "prop2",
      "z_scores",
      "p_values",
      "fdr"
    )
  )

  return(dge_dt)
}

### find all markers -----------------------------------------------------------

#' Find all markers
#'
#' @description
#' This function can be used to run differential gene expression for every
#' group of an unsupervised clustering method for example. You specify a column
#' and the function will start calculating differential gene expression of the
#' first cluster vs. everything else, second cluster vs. everything else, etc.
#' The function will automatically downsample everything else to a random set
#' of 100,000 cells if it should exceed that. This automatic downsampling can
#' be turned off however.
#'
#' @param object `single_cell_exp` class.
#' @param column_of_interest String. The column you wish to use to identify
#' the markers between all combination. Needs to be in the obs table
#' @param method String. Which method to use for the calculations of the DGE.
#' At the moment the only option is `"wilcox"`, but the parameter is reserved
#' for future features.
#' @param alternative String. Test alternative. One of
#' `c("twosided", "greater", "less")`. This function will default to
#' `"greater"`, i.e., genes upregulated in the group.
#' @param min_prop Numeric. The minimum proportion of cells that need to express
#' the gene to be tested in any of the two groups.
#' @param downsampling Boolean. If the other group exceeds 100,000 cells, a
#' random subsample of 100,000 cells will be used.
#' @param seed Integer. Seed that is used for the downsampling.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return data.table with the DGE results from the test.
#'
#' @export
find_all_markers_sc <- S7::new_generic(
  name = "find_markers_sc",
  dispatch_args = "object",
  fun = function(
    object,
    column_of_interest,
    method = "wilcox",
    alternative = c("greater", "less", "twosided"),
    min_prop = 0.05,
    downsampling = TRUE,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method find_all_markers_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(find_all_markers_sc, single_cell_exp) <- function(
  object,
  column_of_interest,
  method = "wilcox",
  alternative = c("greater", "less", "twosided"),
  min_prop = 0.05,
  downsampling = TRUE,
  seed = 42L,
  .verbose = TRUE
) {
  alternative <- match.arg(alternative)

  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::qassert(column_of_interest, "S1")
  checkmate::assertChoice(method, c("wilcox"))
  checkmate::assertChoice(alternative, c("twosided", "greater", "less"))
  checkmate::qassert(min_prop, "N1[0, 1]")
  checkmate::qassert(downsampling, "B1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  obs_data <- object[[c("cell_id", column_of_interest)]][
    !is.na(get(column_of_interest))
  ]

  unique_groups <- unique(obs_data[[column_of_interest]])

  dge_dts <- vector(mode = "list", length = length(unique_groups))

  # loops are my friend...
  for (i in seq_along(unique_groups)) {
    if (.verbose) {
      message(sprintf(
        "Processing group %i out of %i.",
        i,
        length(unique_groups)
      ))
    }

    cell_ids_i <- obs_data[get(column_of_interest) == unique_groups[i], cell_id]
    cell_ids_not_i <- obs_data[
      get(column_of_interest) != unique_groups[i],
      cell_id
    ]
    if (downsampling && length(cell_ids_not_i) > 100000) {
      if (.verbose) {
        message(
          paste(
            " Large number of cells in reference group found.",
            "Downsampling to 100,000 cells."
          )
        )
      }
      set.seed(seed + i)
      cell_ids_not_i <- sample(cell_ids_not_i, 100000)
    }

    dge_results_i <- switch(
      method,
      "wilcox" = rs_calculate_dge_mann_whitney(
        f_path = get_rust_count_cell_f_path(object),
        cell_indices_1 = get_cell_indices(
          x = object,
          cell_ids = cell_ids_i,
          rust_index = TRUE
        ),
        cell_indices_2 = get_cell_indices(
          x = object,
          cell_ids = cell_ids_not_i,
          rust_index = TRUE
        ),
        min_prop = min_prop,
        alternative = alternative,
        verbose = FALSE
      )
    )

    dge_dt_i <- data.table::as.data.table(dge_results_i[c(
      "lfc",
      "prop1",
      "prop2",
      "z_scores",
      "p_values",
      "fdr"
    )])[,
      `:=`(
        gene_id = get_gene_names(object)[dge_results_i$genes_to_keep],
        grp = unique_groups[i]
      )
    ]

    data.table::setcolorder(
      dge_dt_i,
      c(
        "grp",
        "gene_id",
        "lfc",
        "prop1",
        "prop2",
        "z_scores",
        "p_values",
        "fdr"
      )
    )

    dge_dts[[i]] <- dge_dt_i
  }

  dge_dt_final <- data.table::rbindlist(dge_dts)

  return(dge_dt_final)
}

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

## VISION ----------------------------------------------------------------------

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

### VISION with auto-correlation -----------------------------------------------

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
#' ...
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
