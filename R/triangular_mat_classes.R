# class - symmetric cor matrix ----

#' @title Class for symmetric correlation matrices
#'
#' @description
#' The class allows to store the upper triangular matrix of a symmetric
#' correlation matrix in an memory-efficient form and return a data.table or
#' dense R matrix if need be.
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
    },

    #' @description Return the full correlation matrix.
    #'
    #' @param .verbose Boolean. Controls verbosity.
    #'
    #' @return Returns the correlation matrix as a dense R matrix.
    get_cor_matrix = function(.verbose = TRUE) {
      checkmate::qassert(.verbose, "B1")

      if (.verbose)
        message("Generating the full matrix format of the correlation matrix.")

      n <- length(private$features)
      shift <- private$shift

      mat <- rs_upper_triangle_to_dense(private$correlations, shift = shift, n = n)
      colnames(mat) <- rownames(mat) <- private$features

      return(mat)
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


# class - symmetric differential cor matrix ----


#' @title Class for symmetric differential correlation matrices
#'
#' @description
#' The class allows to store the upper triangular matrices of of differential
#' correlation results in an memory-efficient form and return a data.table or
#' dense R matrix for a given parameter if need be.
upper_triangle_diffcor_mat <- R6::R6Class(
  # Class name
  inherit = upper_triangular_cor_mat,
  classname = "upper_triangular_diffcor_mat",

  public = list(
    #' @description Initialises the R6 class.
    #'
    #' @param diff_cor_res A list of differential correlation results.
    #' @param features String vector. The features of the correlation matrices.
    #'
    #' @return Returns the initialised class.
    initialize = function(diff_cor_res, features) {
      # Checks
      checkmate::assertList(diff_cor_res, types = 'double')
      checkmate::assertNames(names(diff_cor_res),
                             must.include = c("r_a", "r_b", "z_score", "p_val"))
      checkmate::assertTRUE(length(unique(purrr::map_dbl(
        diff_cor_res, length
      ))) == 1)
      # New
      private$correlation_a <- diff_cor_res[['r_a']]
      private$correlation_b <- diff_cor_res[['r_b']]
      private$z_scores <- diff_cor_res[['z_score']]
      private$p_val <- diff_cor_res[['p_val']]
      # From the parent class
      private$features <- features
      private$shift <- 1L
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
    #' \item cor_a: The correlation coefficients for data set a between the
    #' features.
    #' \item cor_b: The correlation coefficients for data set b between the
    #' features.
    #' \item z_score: The differential correlation z-score.
    #' \item p_val: The p-value of the differential correlation.
    #' }
    get_data_table = function(factor = FALSE, .verbose = TRUE) {
      # Checks
      checkmate::qassert(factor, "B1")
      checkmate::qassert(.verbose, "B1")
      if (.verbose)
        message("Generating data.table format of the differential correlation data")

      data <- list(
        feature_a = super$get_feature_a(factor = factor),
        feature_b = super$get_feature_b(factor = factor),
        cor_a = private$correlation_a,
        cor_b = private$correlation_b,
        z_score = private$z_scores,
        p_val = private$p_val
      )

      data.table::setDT(data)

      # A cor of 1 causes arctanh to become infinity in the Fisher transformation
      # which causes crazy p-values... Not useful to return
      data <- data[p_val > 0]

      return(data)
    },

    #' @description Return the full correlation matrix of either sample.
    #'
    #' @param to_ret String. Option of `("cor_a", "cor_b")`, pending on which
    #' correlation matrix you want to retrieve.
    #' @param .verbose Boolean. Controls verbosity.
    #'
    #' @return Returns the specified correlation matrix as a dense R matrix.
    get_cor_matrix = function(to_ret = c("cor_a", "cor_b"),
                              .verbose = TRUE) {
      checkmate::assertChoice(to_ret, choices = c("cor_a", "cor_b"))
      checkmate::qassert(.verbose, "B1")

      if (.verbose)
        message("Generating the full matrix format of the desired data.")

      n <- length(private$features)
      shift <- private$shift

      data <- switch(
        to_ret,
        "cor_a" = private$correlation_a,
        "cor_b" = private$correlation_b
      )

      mat <- rs_upper_triangle_to_dense(data, shift = shift, n = n)
      colnames(mat) <- rownames(mat) <- private$features

      return(mat)
    }
  ),

  private = list(
    correlation_a = NULL,
    correlation_b = NULL,
    z_scores = NULL,
    p_val = NULL
  )
)

