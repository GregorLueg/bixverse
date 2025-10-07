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

  obs_data <- object[[c("cell_id", column_of_interest)]]

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

## gene set calculations -------------------------------------------------------
