# anndata parser ---------------------------------------------------------------

# TODO Add option to extract the var table (if exists) and interact with other
# data. Future problem.

#' Class for Anndata
#'
#' @description
#' This class helps dealing with h5ad objects from Python AnnData. You have
#' several options will allow for retrieval of the underlying data.
#'
#' @export
#'
#' @import data.table
#'
#' @examples
#' # parse counts and obs back out of a small h5ad file
#' set.seed(42)
#' counts <- matrix(rpois(50, 5), nrow = 10, ncol = 5)
#' obs <- data.table::data.table(sample_id = sprintf("cell_%i", 1:10))
#' var <- data.table::data.table(var_id = sprintf("gene_%i", 1:5))
#' h5_path <- tempfile(fileext = ".h5ad")
#' write_h5ad_sc_dense(h5_path, counts, obs, var, .verbose = FALSE)
#' parser <- AnnDataParser$new(h5_path)
#' dim(parser$get_raw_counts())
#' head(parser$get_obs_table())
#' unlink(h5_path)
AnnDataParser <- R6::R6Class(
  # Class name
  classname = "AnnDataParser",
  public = list(
    #' @description Initialises the Anndata Parser.
    #'
    #' @param h5_path String. Path to the h5 file.
    #'
    #' @returns Returns the initialised class.
    initialize = function(h5_path) {
      # Checks
      checkmate::qassert(h5_path, "S1")
      checkmate::assertFileExists(h5_path)
      private$h5_path <- h5_path
    },
    #' @description Returns the observation table with all the data from the
    #' h5ad file.
    #'
    #' @returns data.table. The found observations are returned. The pandas index
    #' will be named `sample_id`. Categorical columns are returned as factors.
    get_obs_table = function() {
      private$get_frame("obs", "sample_id")
    },
    #' @description Returns the variable table with all the data from the
    #' h5ad file.
    #'
    #' @returns data.table. The found variables are returned. The pandas index
    #' will be named `var_id`. Categorical columns are returned as factors.
    get_var_info = function() {
      private$get_frame("var", "var_id")
    },
    #' @description Returns the counts that are stored in `X` slot of the
    #' anndata object.
    #'
    #' @returns Returns the count matrix with samples = columns and rows =
    #' features.
    get_raw_counts = function() {
      # TODO Needs a conditional to also deal with the sparse versions that
      # for sure are saved into h5
      raw_counts <- rs_h5_read(private$h5_path, "X")
      dimnames(raw_counts) <- list(
        rs_h5ad_frame(private$h5_path, "var")$idx,
        rs_h5ad_frame(private$h5_path, "obs")$idx
      )
      raw_counts
    },
    #' @description Wrapper function that returns a list of the stored count
    #' data and the metadata found in the h5ad file.
    #'
    #' @returns List with following elements:
    #' \itemize{
    #'  \item metadata - metadata from the respective h5ad file
    #'  \item var_info - metadata on the variables from the respective h5ad
    #'  file.
    #'  \item counts - counts that were found in the h5ad file.
    #' }
    get_key_data = function() {
      counts <- self$get_raw_counts()
      meta_data <- self$get_obs_table()
      var_info <- self$get_var_info()
      return(
        list(
          metadata = meta_data,
          var_info = var_info,
          counts = counts
        )
      )
    }
  ),

  private = list(
    h5_path = NULL,
    # Reads a dataframe group with its index as the first column `id_name`
    get_frame = function(group, id_name) {
      frame <- rs_h5ad_frame(private$h5_path, group)
      dt <- data.table::data.table(frame$idx)
      data.table::setnames(dt, id_name)
      if (length(frame$cols) > 0L) {
        dt <- cbind(dt, data.table::as.data.table(frame$cols))
      }
      dt
    }
  )
)
