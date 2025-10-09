# methods related to batch corrections -----------------------------------------

#' Calculate kBET scores
#'
#' @description
#' This function calculates the k-nearest neighbou batch-effect test (kBET).
#' Briefly, the function leverages a ChiSquare statistic to calculate the
#' differences in batch proportions observed in the neighbourhood of a given
#' cell with the overall batch proportions. If the test is significant for that
#' cell it indicates poor mixing there. Large number of positive tests indicate
#' bad mixing overall. For more details, please see Büttner et al.
#'
#' @param object `single_cell_exp` class.
#' @param batch_column String. The column with the batch information in the
#' obs data of the class.
#' @param threshold Numeric. Number between 0 and 1. Below this threshold, the
#' test is considered significant. Defaults to `0.05`.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @returns A list with the following elements
#' \itemize{
#'   \item kbet_score - Number of significant tests over all cells. 0 indicates
#'   perfect mixing, 1 indicates basically zero mixing between batches.
#'   \item significant_tests - Boolean indicating for which cells the statistic
#'   was below the threshold
#'   \item chisquare_pvals - The p-values of the ChiSquare test.
#' }
#'
#' @export
#'
#' @references Büttner, et al., Nat. Methods, 2019
calculate_kbet_sc <- S7::new_generic(
  name = "calculate_kbet_sc",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    threshold = 0.05,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)


#' @method calculate_kbet_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(calculate_kbet_sc, single_cell_exp) <- function(
  object,
  batch_column,
  threshold = 0.05,
  .verbose = TRUE
) {
  # check
  checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))
  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(threshold, "N1[0, 1]")
  checkmate::qassert(.verbose, "B1")

  # function
  batch_index <- unlist(object[[batch_column]])

  if (!length(levels(factor(batch_index))) > 1) {
    warning("The batch column only has one batch. Returning NULL")
    return(NULL)
  }

  knn_mat <- get_knn_mat(object)

  if (is.null(knn_mat)) {
    warning("No kNN matrix found to calculate kBET. Returning NULL")
    return(NULL)
  }

  rs_res <- rs_kbet(
    knn_mat = knn_mat,
    batch_vector = as.integer(factor(batch_index))
  )

  res <- list(
    kbet_score = sum(rs_res <= threshold) / length(rs_res),
    significant_tests = rs_res <= threshold,
    chisquare_pvals = rs_res
  )

  return(res)
}
