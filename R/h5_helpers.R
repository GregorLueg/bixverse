# anndata parser ---------------------------------------------------------------

# TODO Add option to extract the var table (if exists) and interact with other
# data. Future problem.

#' @title Class for Anndata
#'
#' @description
#' This class helps dealing with h5ad objects from Python AnnData. You have
#' several options will allow for retrieval of the underlying data.
#'
#' @export
#'
#' @import data.table
anndata_parser <- R6::R6Class(
  # Class name
  classname = "anndata_parser",
  public = list(
    #' @description Initialises the Anndata Parser.
    #'
    #' @param h5_path String. Path to the h5 file.
    #'
    #' @return Returns the initialised class.
    initialize = function(h5_path) {
      # Checks
      checkmate::qassert(h5_path, "S1")
      checkmate::assertFileExists(h5_path)
      h5_content <- rhdf5::h5ls(
        h5_path
      ) %>%
        data.table::setDT()
      rhdf5::H5close()
      # Populate the slots
      private$h5_content <- h5_content
      private$h5_path <- h5_path
    },
    #' @description Returns the observation table with all the data from the
    #' h5ad file.
    #'
    #' @return data.table. The found observations are returned. The pandas index
    #' will be named `sample_id`. Remaining columns (if found) will be returned
    #' as factors due to the way the data is stored in h5.
    get_obs_table = function() {
      obs_index <- private$get_obs_index()
      obs_grps <- private$h5_content[group %like% "/obs/"] %>%
        .[, full_path := paste(group, name, sep = "/")]

      categories_paths <- obs_grps[name == "categories", full_path]
      codes_paths <- obs_grps[name == "codes", full_path]

      categories <- purrr::map(categories_paths, \(path) {
        as.character(rhdf5::h5read(file = private$h5_path, name = path))
      })
      codes <- purrr::map(codes_paths, \(path) {
        as.integer(rhdf5::h5read(file = private$h5_path, name = path)) + 1
      })

      colnames <- unique(gsub("/obs/", "", obs_grps[["group"]]))

      if (length(colnames) > 0) {
        obs <- purrr::map2(categories, codes, \(cat, code) {
          factor(cat[code], levels = cat)
        }) %>%
          data.table::setDT() %>%
          `colnames<-`(colnames) %>%
          .[, sample_id := obs_index] %>%
          .[, c("sample_id", colnames), with = FALSE]
      } else {
        obs <- data.table(
          sample_id = obs_index
        )
      }

      rhdf5::H5close()

      return(obs)
    },
    #' @description Returns the variable table with all the data from the
    #' h5ad file.
    #'
    #' @return data.table. The found observations are returned. The pandas index
    #' will be named `var_id`. Remaining columns (if found) will be returned as
    #' factors due to the way the data is stored in h5.
    get_var_info = function() {
      var_index <- private$get_var_index()
      var_grps <- private$h5_content[group %like% "/var/"] %>%
        .[, full_path := paste(group, name, sep = "/")]

      categories_paths <- var_grps[name == "categories", full_path]
      codes_paths <- var_grps[name == "codes", full_path]

      categories <- purrr::map(categories_paths, \(path) {
        as.character(rhdf5::h5read(file = private$h5_path, name = path))
      })
      codes <- purrr::map(codes_paths, \(path) {
        as.integer(rhdf5::h5read(file = private$h5_path, name = path)) + 1
      })

      colnames <- unique(gsub("/var/", "", var_grps[["group"]]))

      if (length(colnames) > 0) {
        vars <- purrr::map2(categories, codes, \(cat, code) {
          code[code == 0] <- NA
          factor(cat[code], levels = cat)
        }) %>%
          data.table::setDT() %>%
          `colnames<-`(colnames) %>%
          .[, var_id := var_index] %>%
          .[, c("var_id", colnames), with = FALSE]
      } else {
        vars <- data.table(
          var_id = var_index
        )
      }

      rhdf5::H5close()

      return(vars)
    },
    #' @description Returns the counts that are stored in `X` slot of the
    #' anndata object.
    #'
    #' @return Returns the count matrix with samples = columns and rows =
    #' features.
    get_raw_counts = function() {
      # TODO Needs a conditional to also deal with the sparse versions that
      # for sure are saved into h5
      var_names <- private$get_var_index()
      obs_names <- private$get_obs_index()
      raw_counts <- rhdf5::h5read(
        file = private$h5_path,
        name = "/X"
      ) %>%
        `rownames<-`(var_names) %>%
        `colnames<-`(obs_names)

      rhdf5::H5close()

      return(raw_counts)
    },
    #' @description Wrapper function that returns a list of the stored count
    #' data and the metadata found in the h5ad file.
    #'
    #' @return List with following elements:
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
    h5_content = NULL,
    h5_path = NULL,
    # Returns the index stored for the samples
    get_obs_index = function() {
      obs_index = private$h5_content[
        group == "/obs" & otype == "H5I_DATASET"
      ] %>%
        .[, path := paste(group, name, sep = "/")] %>%
        .[, path]
      sample_ids <- c(rhdf5::h5read(
        file = private$h5_path,
        name = obs_index
      ))

      rhdf5::H5close()

      return(sample_ids)
    },
    # Returns the index stored for the variables
    get_var_index = function() {
      var_index = private$h5_content[
        group == "/var" & otype == "H5I_DATASET"
      ] %>%
        .[, path := paste(group, name, sep = "/")] %>%
        .[, path]

      var_ids <- c(rhdf5::h5read(
        file = private$h5_path,
        name = var_index
      ))

      rhdf5::H5close()

      return(var_ids)
    }
  )
)
