# sparse matrices --------------------------------------------------------------

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
