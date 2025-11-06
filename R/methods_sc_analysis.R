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

### background distribution ----------------------------------------------------

#### helper --------------------------------------------------------------------

#' Generate random signatures for a null distribution by permuting the data
#'
#' @importFrom stats runif
#' @param eData the data to use for the permutations
#' @param sigData list of signature objects
#' random signature sizes
#' @param num the number of signatures to generate
#' @return A list with two items:
#'
#'   randomSigs: a list of lists of Signature objects.  Each sub-list represents
#'     permutation signatures generated for a specific size/balance
#'
#'   sigAssignments: named factor vector assigning signatures to random background
#'     groups
generatePermutationNull <- function(eData, sigData, num) {
  exp_genes <- rownames(eData)

  sigSize <- vapply(
    sigData,
    function(s) {
      return(length(s@sigDict))
    },
    1
  )

  sigSize <- log10(sigSize)
  sigBalance <- vapply(
    sigData,
    function(s) {
      positive <- sum(s@sigDict >= 0)
      balance <- positive / length(s@sigDict)
      return(balance)
    },
    1
  )

  sigBalance[sigBalance < 0.5] <- 1 - sigBalance[sigBalance < 0.5]

  sigVars <- cbind(sigSize, sigBalance)

  n_components <- 5 # TODO: choose number of components better

  if (nrow(sigVars) <= n_components) {
    n_components <- nrow(sigVars)
    centers <- sigVars
    clusters <- as.factor(seq_len(nrow(sigVars)))
    names(clusters) <- rownames(sigVars)
  } else {
    if (nrow(unique(sigVars)) <= n_components) {
      n_components <- nrow(unique(sigVars))
    }

    km <- kmeans(sigVars, n_components)
    centers <- km$centers

    levels <- as.character(seq(n_components))
    clusters <- factor(km$cluster, levels = levels)
  }

  # Re-order the centers
  row_i <- order(centers[, "sigSize"], centers[, "sigBalance"])

  centers <- centers[row_i, , drop = FALSE]
  levels(clusters) <- as.character(order(row_i))
  rownames(centers) <- as.character(seq_len(n_components))

  # undo the log scaling
  centers[, "sigSize"] <- round(10**centers[, "sigSize"])

  message(
    "Creating ",
    nrow(centers),
    " background signature groups with the following parameters:"
  )
  print(centers) # How do I do this with 'message'??
  message("  signatures per group: ", num)

  randomSigs <- list()
  randomSigAssignments <- character()

  for (cluster_i in rownames(centers)) {
    size <- centers[cluster_i, "sigSize"]
    balance <- centers[cluster_i, "sigBalance"]

    for (j in 1:num) {
      newSigGenes <- sample(exp_genes, min(size, length(exp_genes)))

      upGenes <- floor(balance * size)
      remainder <- (balance * size) %% 1
      if (runif(1, 0, 1) < remainder) {
        upGenes <- upGenes + 1
      }
      newSigSigns <- c(rep(1, upGenes), rep(-1, size - upGenes))

      names(newSigSigns) <- newSigGenes
      newSig <- Signature(
        newSigSigns,
        paste0("RANDOM_BG_", cluster_i, "_", j),
        "x"
      )
      randomSigs[[newSig@name]] <- newSig
      randomSigAssignments[[newSig@name]] <- cluster_i
    }
  }

  randomSigAssignments <- as.factor(randomSigAssignments)

  return(
    list(
      randomSigs = randomSigs,
      sigAssignments = clusters,
      randomSigAssignments = randomSigAssignments
    )
  )
}

#### main ----------------------------------------------------------------------

#' Generate null distribution for VISION scores
#'
#' @description
#' For the calculation of auto-correlation p-values, one needs a background
#' score distribution. This one can be generate with this function. This is
#' based on the VISION method from DeTomaso, et al.
#'
#' @param object `single_cell_exp` class.
#' @param gs_list Named nested list. Same as used for `vision_sc`.
#' @param n_permutations Integer. Number of random signatures per group.
#' @param streaming Boolean. Shall the cell data be streamed in.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @return List with three elements:
#' \itemize{
#'   \item random_scores - Matrix of random signature scores (cells x
#'   random_sigs)
#'   \item sig_assignments - Factor vector mapping real signatures to background
#'   groups
#'   \item random_sig_assignments - Factor vector mapping random sigs to
#'   background groups
#' }
#'
#' @references DeTomaso, et al., Nat. Commun., 2019
#'
#' @export
vision_null_distr_sc <- S7::new_generic(
  name = "vision_null_distr_sc",
  dispatch_args = "object"
)

S7::method(vision_null_distr_sc, single_cell_exp) <- function(
  object,
  gs_list,
  n_permutations = 100,
  streaming = FALSE,
  .verbose = TRUE
) {
  expr_genes <- get_gene_names(object)

  sig_data <- purrr::map(gs_list, \(gs) {
    all_genes <- c(gs$pos, gs$neg)
    signs <- c(rep(1, length(gs$pos)), rep(-1, length(gs$neg)))
    names(signs) <- expr_genes[all_genes + 1]
    signs
  })

  # Generate random signatures using existing function
  null_data <- generatePermutationNull(
    eData = expr_genes,
    sigData = sig_data,
    num = n_permutations
  )

  # Convert random signatures back to index format
  random_gs_list <- purrr::map(null_data$randomSigs, \(sig) {
    pos_genes <- names(sig@sigDict)[sig@sigDict > 0]
    neg_genes <- names(sig@sigDict)[sig@sigDict < 0]
    list(
      pos = match(pos_genes, expr_genes) - 1, # Convert to 0-indexed
      neg = match(neg_genes, expr_genes) - 1
    )
  })

  # Calculate VISION scores for random signatures
  random_scores <- rs_vision(
    f_path = get_rust_count_cell_f_path(object),
    gs_list = random_gs_list,
    cells_to_keep = get_cells_to_keep(object),
    streaming = streaming,
    verbose = .verbose
  )

  colnames(random_scores) <- names(random_gs_list)
  rownames(random_scores) <- get_cell_names(object, filtered = TRUE)

  list(
    random_scores = random_scores,
    sig_assignments = null_data$sigAssignments,
    random_sig_assignments = null_data$randomSigAssignments
  )
}
