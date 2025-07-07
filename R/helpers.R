# parallelisation stuff --------------------------------------------------------

#' Helper function for core detection.
#'
#' @description
#' Identifies the number of cores/workers to use for parallel tasks. Will
#' default to half of the available cores, to a maximum of `abs_max_workers`.
#'
#' @param abs_max_workers Integer. Absolute maximum number of cores to use.
#' Defaults to `8L`.
#'
#' @returns Integer. Number of cores to use.
get_cores <- function(abs_max_workers = 8L) {
  checkmate::qassert(abs_max_workers, "I1")
  cores <- as.integer(parallel::detectCores() / 2)
  cores <- ifelse(cores >= abs_max_workers, abs_max_workers, cores)
  return(cores)
}

# user options -----------------------------------------------------------------

#' Helper function for user option selection
#'
#' @param options Character vector. The options the user has to choose from.
#'
#' @returns The chosen option
select_user_option <- function(options) {
  cat("Please select an option:\n")
  for (i in seq_along(options)) {
    cat(i, ": ", options[i], "\n")
  }

  while (TRUE) {
    cat("\nEnter the number of your choice: ")
    user_choice <- as.integer(readLines(n = 1))

    if (
      !is.na(user_choice) && user_choice > 0 && user_choice <= length(options)
    ) {
      return(options[user_choice])
    } else {
      cat("Invalid selection. Please try again.\n")
    }
  }
}

# gene set libraries -----------------------------------------------------------

#' Get the Gene Ontology data human
#'
#' @returns data.table with the data ready for [bixverse::gene_ontology_data()].
#'
#' @export
get_go_human_data <- function() {
  path <- system.file("extdata", "go_data_hs.parquet", package = "bixverse")
  if (path != "") {
    go_data_dt <- arrow::read_parquet(path) %>%
      data.table::setDT()

    return(go_data_dt)
  } else {
    stop("The expected .parquet file was not found.")
  }
}

# sparse matrices --------------------------------------------------------------

#' Transform an upper triangle-stored matrix to a sparse one
#'
#' @param upper_triangle_vals Numerical vector. The values of the upper triangle
#' stored in a row major format.
#' @param shift Integer. Did you shift the diagonal up.
#' @param n Integer. Number of columns and rows of the symmetric matrix.
#'
#' @returns The sparse matrix.
#'
#' @export
upper_triangle_to_sparse <- function(upper_triangle_vals, shift, n) {
  # checks
  checkmate::qassert(upper_triangle_vals, "N+")
  checkmate::qassert(shift, "I1")
  checkmate::qassert(n, "I1")

  # functions
  data <- rs_upper_triangle_to_sparse(upper_triangle_vals, as.integer(shift), n)

  matrix = Matrix::sparseMatrix(
    i = data$row_indices + 1,
    p = data$col_ptr,
    x = unlist(data$data),
    dims = c(n, n)
  )

  return(matrix)
}

# inflection points ------------------------------------------------------------

#' Identify the inflection point for elbow-like data
#'
#' @description
#' This function will identify the index of the inflection point of an arbitrary
#' series `x ~ y` via the biggest increase in the first derivative. Useful to
#' identify key points in elbow plots.
#'
#' @param x,y The x and y values.
#' @param span The span parameter for the loess function.
#'
#' @return A list containing:
#' \itemize{
#'   \item inflection_idx - Index of the inflection point
#'   \item gradient_change - Absolute change in the first derivative
#' }
get_inflection_point <- function(x, y, span = 0.5) {
  # Checks
  checkmate::assertNumeric(x, len = length(y))
  checkmate::assertNumeric(y, len = length(x))
  checkmate::qassert(span, "R+[0,1]")
  # Function body
  span <- max(0.1, min(1.0, span))
  fit <- tryCatch(
    {
      fit <- loess(
        y ~ x,
        span = span,
        family = "gaussian",
        degree = 2L,
        normalize = TRUE
      )
      list(fit = fit, warning = FALSE)
    },
    warning = function(w) {
      return(list(fit = NULL, warning = TRUE))
    }
  )

  # early return if the inflection point function is unhappy
  if (fit$warning) {
    return(list(
      inflection_idx = NULL,
      gradient_change = NULL
    ))
  }

  py <- predict(fit$fit, x)

  n <- length(x)
  gradient <- numeric(n)
  gradient[1] <- (py[2] - py[1]) / (x[2] - x[1])
  gradient[n] <- (py[n] - py[n - 1]) / (x[n] - x[n - 1])

  for (i in 2:(n - 1)) {
    gradient[i] <- (py[i + 1] - py[i - 1]) / (x[i + 1] - x[i - 1])
  }

  gradient_change <- diff(gradient)
  # One after the point with the biggest delta
  inflection_idx <- which.max(gradient_change) + 1

  return(
    list(inflection_idx = inflection_idx, gradient_change = gradient_change)
  )
}
