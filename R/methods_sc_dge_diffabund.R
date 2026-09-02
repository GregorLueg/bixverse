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
#' @param object `SingleCells` or `SingleCellsSubset` class.
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
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
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

#' @method find_markers_sc ScOrScSubset
#'
#' @importFrom zeallot %<-%
#' @importFrom magrittr %>%
S7::method(find_markers_sc, ScOrScSubset) <- function(
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
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) ||
      S7::S7_inherits(object, SingleCellsSubset)
  )
  assertCellsExist(object, cells_1)
  assertCellsExist(object, cells_2)
  checkmate::assertChoice(method, c("wilcox"))
  checkmate::assertChoice(alternative, c("twosided", "greater", "less"))
  checkmate::qassert(min_prop, "N1[0, 1]")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

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
      verbose = parse_verbosity(.verbose)
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
#' @param object `SingleCells` or `SingleCellsSubset` class.
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
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @return data.table with the DGE results from the test.
#'
#' @export
find_all_markers_sc <- S7::new_generic(
  name = "find_all_markers_sc",
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

#' @method find_all_markers_sc ScOrScSubset
#'
#' @importFrom zeallot %<-%
#' @importFrom magrittr %>%
S7::method(find_all_markers_sc, ScOrScSubset) <- function(
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
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) ||
      S7::S7_inherits(object, SingleCellsSubset)
  )
  checkmate::qassert(column_of_interest, "S1")
  checkmate::assertChoice(method, c("wilcox"))
  checkmate::assertChoice(alternative, c("twosided", "greater", "less"))
  checkmate::qassert(min_prop, "N1[0, 1]")
  checkmate::qassert(downsampling, "B1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

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
        verbose = 0L
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

### find specific markers ------------------------------------------------------

#' Find markers that are specific to a cell group
#'
#' @description
#' This function scores one group of cells against every other group of a column
#' separately, instead of against all of them pooled together. That difference
#' matters for markers. A pooled test is dominated by whichever rival
#' contributes the most cells, so a gene that is high in the reference and just
#' as high in one small rival still comes out looking like a clean marker. Here
#' the gene has to hold up against every rival, and the per-gene summaries
#' (`min_auroc`, `median_auroc`, `min_rank`) tell you whether it does.
#'
#' Leave `reference_group` as `NULL` to run every group of the column as the
#' reference in turn, or name one group to only get that arm.
#'
#' The summaries rank on AUROC rather than the p-value on purpose. Group sizes
#' vary a lot in practice and p-values scale with the group sizes, so a large
#' rival would otherwise crowd out a small one regardless of effect size.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param column_of_interest String. The column in the obs table holding the
#' groups, e.g. a cell type annotation or a clustering.
#' @param reference_group Optional string. The group to use as the reference. If
#' `NULL`, every group of `column_of_interest` is used as the reference in turn.
#' @param method String. Which method to use for the calculations of the DGE.
#' At the moment the only option is `"wilcox"`, but the parameter is reserved
#' for future features.
#' @param alternative String. Test alternative. One of
#' `c("twosided", "greater", "less")`. This function will default to
#' `"greater"`, i.e., genes upregulated in the reference group.
#' @param min_prop Numeric. The minimum proportion of cells that need to express
#' the gene in at least one of the groups. Applied once, globally, so every
#' comparison's FDR is calculated over the same gene set.
#' @param downsampling Boolean. If any group exceeds 100,000 cells, a random
#' subsample of 100,000 cells is used for it. The subsample is drawn once per
#' group, so a group is represented by the same cells in every arm.
#' @param seed Integer. Seed that is used for the downsampling.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @return A `ScSpecificMarkers` class with the following elements
#' \itemize{
#'   \item summary - data.table. Per gene and reference group, the summaries
#'   across all rivals: `prop_ref`, `median_auroc`, `min_auroc`, `mean_auroc`,
#'   `max_auroc`, `worst_rival` (the rival achieving `min_auroc`), `min_rank`
#'   (best AUROC rank the gene reaches against any single rival), the
#'   Simes-combined p-value with its FDR, and the maximum p-value with its FDR.
#'   \item per_comparison - data.table. Per gene, reference group and rival, the
#'   underlying `auroc`, `lfc`, `prop_ref`, `prop_rival`, `z_scores`, `p_values`
#'   and `fdr`.
#'   \item params - List. The parameters the run used.
#' }
#'
#' @references Soneson and Robinson, Nat Methods, 2018; Lun, et al.,
#' F1000Research, 2016
#'
#' @export
find_specific_markers_sc <- S7::new_generic(
  name = "find_specific_markers_sc",
  dispatch_args = "object",
  fun = function(
    object,
    column_of_interest,
    reference_group = NULL,
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

#' @method find_specific_markers_sc ScOrScSubset
S7::method(find_specific_markers_sc, ScOrScSubset) <- function(
  object,
  column_of_interest,
  reference_group = NULL,
  method = "wilcox",
  alternative = c("greater", "less", "twosided"),
  min_prop = 0.05,
  downsampling = TRUE,
  seed = 42L,
  .verbose = TRUE
) {
  alternative <- match.arg(alternative)

  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) ||
      S7::S7_inherits(object, SingleCellsSubset)
  )
  checkmate::qassert(column_of_interest, "S1")
  checkmate::qassert(reference_group, c("S1", "0"))
  checkmate::assertChoice(method, c("wilcox"))
  checkmate::assertChoice(alternative, c("twosided", "greater", "less"))
  checkmate::qassert(min_prop, "N1[0, 1]")
  checkmate::qassert(downsampling, "B1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  .assert_group_by(object, column_of_interest)

  obs_data <- get_sc_obs(
    object,
    cols = c("cell_idx", column_of_interest),
    filtered = TRUE
  )

  grp_vec <- as.character(obs_data[[column_of_interest]])
  annotated <- !is.na(grp_vec)

  if (!all(annotated)) {
    warning(sprintf(
      "Dropping %i cell(s) with no '%s' annotation.",
      sum(!annotated),
      column_of_interest
    ))
  }

  # obs cell_idx is 1-based, Rust wants 0-based. split() orders the groups by
  # name, which keeps the downsampling seed stable across the reference arms
  cell_groups <- split(
    as.integer(obs_data$cell_idx[annotated] - 1L),
    grp_vec[annotated]
  )

  if (length(cell_groups) < 2L) {
    stop(sprintf(
      "Column '%s' has fewer than two groups. Nothing to compare against.",
      column_of_interest
    ))
  }

  if (length(cell_groups) == 2L) {
    warning(
      paste(
        "Only two groups found. With a single rival this collapses to the",
        "pairwise case, see find_markers_sc()."
      )
    )
  }

  if (!is.null(reference_group)) {
    checkmate::assertChoice(reference_group, names(cell_groups))
  }

  grp_names <- names(cell_groups)

  if (downsampling) {
    cell_groups <- purrr::imap(cell_groups, \(indices, grp) {
      if (length(indices) <= 100000) {
        return(indices)
      }
      if (.verbose) {
        message(sprintf(
          " Large number of cells in group '%s'. Downsampling to 100,000.",
          grp
        ))
      }
      set.seed(seed + which(grp_names == grp))
      sample(indices, 100000)
    })
  }

  ref_groups <- if (is.null(reference_group)) {
    grp_names
  } else {
    reference_group
  }

  # the per-arm progress messages cover normal verbosity, so only hand the Rust
  # side the detailed level where its per-comparison timings are wanted
  rust_verbosity <- if (parse_verbosity(.verbose) >= 2L) 2L else 0L

  summary_dts <- vector(mode = "list", length = length(ref_groups))
  per_comparison_dts <- vector(mode = "list", length = length(ref_groups))

  for (i in seq_along(ref_groups)) {
    if (.verbose) {
      message(sprintf(
        "Processing reference group %i out of %i.",
        i,
        length(ref_groups)
      ))
    }

    rival_grps <- setdiff(grp_names, ref_groups[i])

    res_i <- switch(
      method,
      "wilcox" = rs_calculate_dge_one_vs_many(
        f_path = get_rust_count_cell_f_path(object),
        cell_indices_ref = cell_groups[[ref_groups[i]]],
        cell_indices_other = cell_groups[rival_grps],
        min_prop = min_prop,
        alternative = alternative,
        verbose = rust_verbosity
      )
    )

    if (!any(res_i$genes_to_keep)) {
      warning(sprintf(
        "No gene passed the proportion filter for group '%s'. Skipping it.",
        ref_groups[i]
      ))
      next
    }

    melted_i <- .melt_one_vs_many_res(
      rs_res = res_i,
      gene_names = get_gene_names(object),
      ref_grp = ref_groups[i],
      rival_grps = rival_grps
    )

    summary_dts[[i]] <- melted_i$summary
    per_comparison_dts[[i]] <- melted_i$per_comparison
  }

  summary_dt <- data.table::rbindlist(summary_dts)

  # early return if every arm was skipped by the proportion filter
  if (nrow(summary_dt) == 0L) {
    warning(
      paste(
        "No gene passed the proportion filter for any group.",
        "Consider lowering min_prop. Returning NULL."
      )
    )
    return(NULL)
  }

  params <- list(
    method = method,
    alternative = alternative,
    min_prop = min_prop,
    column_of_interest = column_of_interest,
    reference_group = reference_group,
    downsampling = downsampling,
    seed = seed
  )

  new_sc_specific_markers(
    summary = summary_dt,
    per_comparison = data.table::rbindlist(per_comparison_dts),
    params = params
  )
}

## differential abundance ------------------------------------------------------

### miloR ----------------------------------------------------------------------

#' Generate an miloR abundance object for differential abundance testing
#'
#' @description
#' This function implements the miloR differential abundance testing on top
#' of the kNN graph. The general idea of the approach is to use the kNN graph
#' generated from the single cell data, generate representative neighbourhoods
#' and calculate differential abundances within these neighbourhoods. For
#' further details on the method, please refer to Dann, et al. This function
#' will take a `SingleCells` class, run the neighbourhood detection,
#' count the occurrences of a sample and return a `sc_miloR` class for
#' subsequent differential abundance testing and further annotations.
#'
#' @param object `SingleCells` class.
#' @param sample_id_col Character. The column in the obs table representing
#' the sample identifier to count.
#' @param embd_to_use Character. The embedding to use for the refinement
#' procedure. Please use the same here as you used to generate the neighbours!
#' Defaults to `"pca"`.
#' @param no_embd_to_use Optional integer. If you only want to use a subset of
#' the embedding.
#' @param miloR_params A list, please see [bixverse::params_sc_miloR()]. The
#' list has the following parameters:
#' \itemize{
#'   \item prop - Numeric. Proportion of cells to sample as neighbourhood
#'   indices. Must be in (0,1).
#'   \item k_refine - Integer. Number of neighbours to use for refinement.
#'   \item refinement_strategy - String. Strategy for refining sampled indices.
#'   One of `c("approximate", "bruteforce", "index")`.
#'   \item index_type - String. Type of kNN index to use. One of
#'   `c("nndescent", "ivf", "hnsw", "annoy", "exhaustive")`.
#'   \item knn - List of kNN parameters. See [bixverse::params_knn_defaults()]
#'   for available parameters and their defaults.
#' }
#' @param seed Integer. Seed for reproducibility
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @references Dann, et al., Nat Biotechnol, 2022
#'
#' @export
get_miloR_abundances_sc <- S7::new_generic(
  name = "get_miloR_abundances_sc",
  dispatch_args = "object",
  fun = function(
    object,
    sample_id_col,
    embd_to_use = "pca",
    no_embd_to_use = NULL,
    miloR_params = params_sc_miloR(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_miloR_abundances_sc SingleCells
#'
#' @export
S7::method(get_miloR_abundances_sc, SingleCells) <- function(
  object,
  sample_id_col,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  miloR_params = params_sc_miloR(),
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(sample_id_col, "S1")
  assertScMiloR(miloR_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # hard tier: the neighbourhood counts are aligned to the kept cells
  assert_sc_state(object, artefacts = c(embd_to_use, "knn"))

  samples <- unlist(object[[sample_id_col]], use.names = FALSE)

  # a cell with no sample label cannot be counted into any sample, and dropping
  # it silently would deflate the neighbourhood counts unevenly
  if (anyNA(samples)) {
    stop(sprintf(
      "Column '%s' holds %i missing value(s). Every cell needs a sample.",
      sample_id_col,
      sum(is.na(samples))
    ))
  }

  samples <- factor(samples)

  embd <- get_embedding(x = object, embd_name = embd_to_use)

  if (!is.null(no_embd_to_use)) {
    to_take <- min(c(no_embd_to_use, ncol(embd)))
    embd <- embd[, 1:to_take]
  }

  knn_data <- get_knn_mat(object)

  # early return if no kNN graph was generated
  if (is.null(knn_data)) {
    warning(
      paste(
        "No kNN data could be found on the object.",
        "Please run find_neighbours_sc(). Returning NULL."
      )
    )
    return(NULL)
  }

  # check k_refine and dimensionality
  if (miloR_params$k_refine > ncol(embd)) {
    warning(
      paste(
        "You chose a k refinement larger than the number of features",
        "in the embedding. Will change to ncol(embd)."
      )
    )
    miloR_params[["k_refine"]] <- ncol(embd)
  }

  milor_res <- rs_make_milor_nhoods(
    embd = embd,
    knn_indices = knn_data,
    sample_ids = as.integer(samples) - 1L,
    n_samples = nlevels(samples),
    milor_params = miloR_params,
    seed = seed,
    verbose = parse_verbosity(.verbose)
  )

  nhoods <- Matrix::sparseMatrix(
    i = milor_res$nhoods_i + 1,
    j = milor_res$nhoods_j + 1,
    x = milor_res$nhoods_x,
    dims = c(milor_res$nrows, milor_res$ncols)
  )

  # Rust returns the full neighbourhoods x samples grid with explicit zeros,
  # so unlike the old table() this cannot come back ragged
  sample_counts <- milor_res$sample_counts
  dimnames(sample_counts) <- list(
    seq_len(nrow(sample_counts)),
    levels(samples)
  )

  params <- miloR_params
  params[["used_emb"]] <- embd_to_use
  params[["no_embd_to_use"]] <- no_embd_to_use
  params[["index_cell"]] <- milor_res$index_cell

  miloR_obj <- new_sc_miloR_res(
    nhoods = nhoods,
    sample_counts = sample_counts,
    spatial_dist = milor_res$kth_distances,
    nhood_overlap = milor_res$nhood_overlap,
    params = params
  )

  miloR_obj
}

### meld -----------------------------------------------------------------------

#' Run MELD signal smoothing for differential abundance estimation
#'
#' @description
#' This function implements MELD  for estimating sample-associated density on a
#' cell manifold. The general idea is to smooth a binary sample indicator matrix
#' over the kNN graph via spectral filtering, yielding per-cell likelihood
#' estimates for each sample condition. For further details on the method,
#' please refer to Burkhardt, et al. This function will take a `SingleCells`
#' class and return the smoothed density estimates per cell and condition.
#'
#' @param object `SingleCells` class.
#' @param sample_id_col Character. The column in the obs table representing
#' the sample identifier.
#' @param embd_to_use Character. The embedding to use for kNN graph
#' construction. Please use the same here as you used to generate the
#' neighbours. Defaults to `"pca"`.
#' @param no_embd_to_use Optional integer. If you only want to use a subset of
#' the embedding dimensions.
#' @param meld_params A list, please see [bixverse::params_meld()]. The list
#' has the following parameters:
#' \itemize{
#'   \item beta - Numeric. Smoothing strength; larger values produce smoother
#'   densities.
#'   \item offset - Numeric. Shift of the filter centre in the rescaled
#'   spectrum. Must be in `[0, 1]`.
#'   \item order - Numeric. Filter falloff sharpness; larger values approach a
#'   square low-pass.
#'   \item filter - Character. Filter family. One of `c("heat", "laplacian")`.
#'   \item chebyshev_order - Integer. Number of Chebyshev coefficients. Must
#'   be >= 2.
#'   \item lap_type - Character. Type of Laplacian. One of
#'   `c("combinatorial", "normalised")`.
#'   \item normalise_indicators - Logical. If `TRUE`, indicator columns are
#'   divided by their column sum before filtering.
#'   \item landmark - Logical. Whether to use landmark approximation. Useful
#'   for large data sets.
#'   \item n_landmarks - Integer. Number of landmarks to use.
#'   \item knn - List of kNN parameters. See [bixverse::params_knn_defaults()]
#'   for available parameters and their defaults.
#' }
#' @param seed Integer. Seed for reproducibility.
#' @param .verbose Boolean or integer. Controls verbosity and returns run
#' times. `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` ->
#' detailed verbosity.
#'
#' @return A list with:
#' \itemize{
#'  \item raw_scores - The raw MELD scores
#'  \item norm_scores - Negative values were clamped to 0 and the rows L1
#'  normalised. This yields probability-like values.
#' }
#'
#' @references Burkhardt, et al. Nat. Biotechnol., 2021.
#'
#' @export
meld_sc <- S7::new_generic(
  name = "meld_sc",
  dispatch_args = "object",
  fun = function(
    object,
    sample_id_col,
    embd_to_use = "pca",
    no_embd_to_use = NULL,
    meld_params = params_meld(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method meld_sc SingleCells
#'
#' @export
#'
#' @importFrom zeallot %<-%
#' @importFrom magrittr %>%
S7::method(meld_sc, SingleCells) <- function(
  object,
  sample_id_col,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  meld_params = params_meld(),
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(sample_id_col, "S1")
  assertMeldParams(meld_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # hard tier: the densities are returned per kept cell
  assert_sc_state(object, artefacts = c(embd_to_use, "knn"))

  samples <- as.factor(unlist(object[[sample_id_col]], use.names = FALSE))

  embd <- get_embedding(x = object, embd_name = embd_to_use)

  if (!is.null(no_embd_to_use)) {
    to_take <- min(c(no_embd_to_use, ncol(embd)))
    embd <- embd[, 1:to_take]
  }

  knn_data <- get_knn_obj(object)

  if (is.null(knn_data)) {
    warning(
      paste(
        "No kNN data could be found on the object.",
        "The kNN graph will be regenerated with the specified embedding."
      )
    )
    return(NULL)
  }

  meld_res <- rs_meld_sc(
    embd = embd,
    knn_data = knn_data,
    meld_params = meld_params,
    labels = as.integer(samples),
    n_labels = as.integer(length(levels(samples))),
    seed = seed,
    verbose = parse_verbosity(.verbose)
  )

  colnames(meld_res$raw_scores) <- colnames(meld_res$norm_scores) <- levels(
    samples
  )
  rownames(meld_res$raw_scores) <- rownames(
    meld_res$norm_scores
  ) <- get_cell_names(object, filtered = TRUE)

  meld_res
}

### NEBULA ---------------------------------------------------------------------

#' Run NEBULA on single cells
#'
#' @description
#' Fits NEBULA's negative binomial gamma mixed model, which is the test to
#' reach for when the cells are not independent because they came from a
#' handful of donors. The variance is split into a subject-level random effect
#' and a cell-level overdispersion, so the donor structure is modelled rather
#' than ignored, and the nominal false discovery rate holds.
#'
#' The design is a formula evaluated against the obs table, so anything in
#' there can go into the model. Cells with a missing design or subject value
#' are dropped with a warning.
#'
#' Genes are streamed and fitted in batches, which is exact: NEBULA is
#' gene-independent. A fit is milliseconds to seconds per gene though, so
#' running it over the full gene axis is rarely what you want. Restrict
#' `genes_to_use` to the highly variable genes or a candidate set.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param subject_col String. The column in the obs table holding the subject
#' (donor) identifier. This is what the random effect is over. Not the same
#' thing as a sample or a batch, unless they happen to coincide.
#' @param design Formula. The experimental design, evaluated against the obs
#' table, e.g. `~ condition` or `~ condition + age`. Include the intercept.
#' @param coef Optional integer or character. Which coefficient of the design
#' the Wald test reports, as a 1-based column position or a column name.
#' Defaults to the last column.
#' @param contrast Optional numeric vector. One weight per design column.
#' Mutually exclusive with `coef`.
#' @param genes_to_use Optional character vector. The genes to fit. Defaults to
#' every gene in the object, which is usually too many.
#' @param offset Optional numeric vector. Strictly positive scaling factor per
#' cell, aligned to the cells that survive the design. Defaults to `NULL`,
#' which uses the library sizes.
#' @param nebula_params A list, see [bixverse::params_nebula()]. The list has
#' the following parameters:
#' \itemize{
#'   \item nebula_method - String. One of `c("ln", "hl")`.
#'   \item min_sigma, max_sigma - Numeric. Bounds on the subject-level
#'   overdispersion.
#'   \item min_phi, max_phi - Numeric. Bounds on the cell-level overdispersion.
#'   \item cutoff_cell - Numeric. When to refit both overdispersions.
#'   \item kappa - Numeric. When to trust the stage-one subject overdispersion.
#'   \item cpc - Numeric. Minimum mean count per cell for a gene to be tested.
#'   \item mincp - Integer. Minimum number of cells expressing a gene.
#'   \item reml - Boolean. Restricted maximum likelihood.
#'   \item eps - Numeric. Optimiser stopping tolerance.
#'   \item gene_batch_size - Integer. Genes read and fitted per batch.
#'   \item shrink_dispersion - Boolean. Empirical Bayes shrinkage of the
#'   cell-level overdispersions.
#' }
#' @param .verbose Boolean or integer. Controls verbosity and returns run
#' times. `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` ->
#' detailed verbosity.
#'
#' @returns A `ScNebula` class, see [bixverse::new_sc_nebula_res()], with
#' \itemize{
#'   \item results - data.table. One row per gene that survived NEBULA's
#'   expression filter, with the Wald test and both overdispersions.
#'   \item coefficients - Numeric matrix of genes x coefficients.
#'   \item se - Numeric matrix of genes x coefficients.
#'   \item params - List. The parameters the run used.
#' }
#'
#' @references He, et al., Commun Biol, 2021
#'
#' @export
nebula_sc <- S7::new_generic(
  name = "nebula_sc",
  dispatch_args = "object",
  fun = function(
    object,
    subject_col,
    design,
    coef = NULL,
    contrast = NULL,
    genes_to_use = NULL,
    offset = NULL,
    nebula_params = params_nebula(),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method nebula_sc ScOrScSubset
S7::method(nebula_sc, ScOrScSubset) <- function(
  object,
  subject_col,
  design,
  coef = NULL,
  contrast = NULL,
  genes_to_use = NULL,
  offset = NULL,
  nebula_params = params_nebula(),
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) ||
      S7::S7_inherits(object, SingleCellsSubset)
  )
  checkmate::qassert(subject_col, "S1")
  checkmate::assertFormula(design)
  checkmate::qassert(genes_to_use, c("0", "S+"))
  checkmate::qassert(offset, c("0", "N+"))
  assertNebulaParams(nebula_params)
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  obs <- get_sc_obs(
    object,
    cols = unique(c("cell_idx", subject_col, all.vars(design))),
    filtered = TRUE
  )

  inputs <- .nebula_design(
    obs = obs,
    design = design,
    subject_col = subject_col
  )
  tested <- .resolve_tested_sc(
    design = inputs$design_mat,
    coef = coef,
    contrast = contrast
  )

  # obs cell_idx is 1-based, Rust wants 0-based global positions
  cells_to_keep <- as.integer(inputs$obs$cell_idx - 1L)

  if (!is.null(offset)) {
    checkmate::assertNumeric(
      offset,
      len = length(cells_to_keep),
      lower = .Machine$double.eps,
      any.missing = FALSE
    )
  }

  # Rust wants one subject label per cell in the store, not per selected cell.
  # Unselected positions are never read, so they can stay at zero.
  n_cells_total <- nrow(get_sc_obs(object, cols = "cell_idx", filtered = FALSE))
  subject_ids <- integer(n_cells_total)
  subject_ids[inputs$obs$cell_idx] <- as.integer(inputs$subject_fct) - 1L

  gene_indices <- if (is.null(genes_to_use)) {
    seq_along(get_gene_names(object)) - 1L
  } else {
    get_gene_indices(x = object, gene_ids = genes_to_use, rust_index = TRUE)
  }

  res <- rs_nebula_sc(
    f_path_genes = get_rust_count_gene_f_path(object),
    f_path_cells = get_rust_count_cell_f_path(object),
    cells_to_keep = cells_to_keep,
    gene_indices = as.integer(gene_indices),
    subject_ids = subject_ids,
    design = inputs$design_mat,
    offset = offset,
    nebula_params = c(nebula_params, tested),
    verbose = parse_verbosity(.verbose)
  )

  params <- list(
    subject_col = subject_col,
    design = deparse1(design),
    tested = if (is.null(tested$coef)) {
      "contrast"
    } else {
      colnames(inputs$design_mat)[tested$coef + 1L]
    },
    n_cells = length(cells_to_keep),
    n_subjects = nlevels(inputs$subject_fct),
    n_genes_requested = length(gene_indices),
    nebula_params = nebula_params
  )

  .nebula_res_to_class(
    res = res,
    gene_names = get_gene_names(object),
    design_mat = inputs$design_mat,
    params = params
  )
}
