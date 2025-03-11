# class - symmetric cor matrix ----

#' @title Class for symmetric cor matrices
#'
#' @description
#' The class allows to store the upper triangular matrix of a symmetric
#' correlation matrix in an memory-efficient form and return data.table if need
#' be.
upper_triangular_cor_mat <- R6::R6Class(
  # Class name
  classname = "upper_triangular_cor_mat",

  # Public
  public = list(
    #' @description Initialises the R6 class.
    #'
    #' @param cor_coef Numerical vector. The correlation coefficients of the
    #' upper triangular correlation matrix stored as a row-major vector
    #' @param features String vector. The features of the correlation matrix.
    #' @param shift Integer. Was a shift applied during the calculation of the
    #' upper triangular matrix. Typically 0 (diagonal included) or 1 (diagonal
    #' not included)
    #'
    #' @return Returns the initialised class.
    initialize = function(cor_coef, features, shift) {
      # Checks
      checkmate::qassert(cor_coef, "N+")
      checkmate::qassert(features, "S+")
      checkmate::qassert(shift, "I1")
      # Gauss trick to calculate what would be expected in terms of number of
      # features given the length of the correlation coefficients and shift
      no_features <- length(features) - shift
      expected_features <- floor(no_features / 2) * (no_features + 1)
      if (no_features %% 2 == 1) {
        expected_features <- expected_features + ceiling(no_features / 2)
      }
      checkmate::assertTRUE(expected_features == length(cor_coef))
      # Populate the slots
      private$correlations <- cor_coef
      private$features <- features
      private$shift <- shift
    },


    #' @description Print the class
    #'
    #' @return Returns the initialised class
    print = function() {
      cat("R6-based class for storing upper triangular, symmetric matrices",
          sep = "")
      cat(" Size of the numerical vector: ",
          length(private$correlations),
          ".\n",
          sep = "")
      cat(" Applied shift: ", private$shift, ".", sep = "")
    },

    #' @description Returns the data in form of a data.table.
    #'
    #' @param factor Boolean. Shall the string columns be transformed into
    #' factors. Reduces size of the object; however, takes longer to generate.
    #' @param .verbose Boolean. Controls verbosity.
    #'
    #' @return A data.table with three columns:
    #' \itemize{
    #' \item feature_a: The name of the first feature in the correlation matrix.
    #' \item feature_b: The name of the second feature in the correlation
    #' matrix.
    #' \item cor: The correlation coefficients between these two features.
    #' }
    get_data_table = function(factor = FALSE, .verbose = TRUE) {
      # Checks
      checkmate::qassert(factor, "B1")
      checkmate::qassert(.verbose, "B1")
      if (.verbose)
        message("Generating data.table format of the correlation matrix.")

      data <- list(
        feature_a = private$get_feature_a(factor = factor),
        feature_b = private$get_feature_b(factor = factor),
        cor = private$correlations
      )

      data.table::setDT(data)

      return(data)
    }
  ),
  # Private
  private = list(
    # Slots
    correlations = NULL,
    features = NULL,
    shift = NULL,
    # Private functions
    get_feature_a = function(factor = FALSE) {
      # Checks
      checkmate::qassert(factor, "B1")
      # Get data
      total_len <- length(private$features)
      shift <- private$shift

      feature_a <- purrr::map(1:total_len, \(idx) {
        rep(private$features[idx], total_len - idx + 1 - shift)
      })
      feature_a <- do.call(c, feature_a)
      if (factor)
        feature_a <- factor(feature_a)

      return(feature_a)
    },

    get_feature_b = function(factor = FALSE) {
      # Checks
      checkmate::qassert(factor, "B1")
      # Get data
      total_len <- length(private$features)
      shift <- private$shift
      feature_b <- purrr::map(1:total_len, \(idx) {
        start_point <- idx + shift
        if (start_point <= total_len)
          private$features[start_point:total_len]
        else
          character(0)
      })
      feature_b <- do.call(c, feature_b)
      if (factor)
        feature_b <- factor(do.call(c, feature_b))

      return(feature_b)
    }
  )
)
