# write h5ad type formats ------------------------------------------------------

#' Helper function to write data to h5ad format
#'
#' @param f_path String. The filepath to which to save the data
#' @param counts Sparse matrix. Needs to be of class `dgRMatrix` or
#' `dgCMatrix`.
#' @param obs data.table. The observations. Needs to have
#' `nrow(obs) == nrow(counts)`.
#' @param var data.table. The variable data. Needs to have
#' `ncol(var) == ncol(counts)`.
#' @param overwrite Boolean. Shall any found h5ad file be overwritten.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return Returns invisible
#'
#' @export
write_h5ad_sc <- function(
  f_path,
  counts,
  obs,
  var,
  overwrite = TRUE,
  .verbose = TRUE
) {
  # checks
  checkmate::assertPathForOutput(f_path, overwrite = TRUE)
  checkmate::assert(
    checkmate::testClass(counts, "dgRMatrix"),
    checkmate::testClass(counts, "dgCMatrix")
  )
  checkmate::assertDataTable(
    obs,
    min.rows = nrow(counts),
    max.rows = nrow(counts)
  )
  checkmate::assertDataTable(
    var,
    min.rows = ncol(counts),
    max.rows = ncol(counts)
  )
  checkmate::qassert(overwrite, "B1")
  checkmate::qassert(.verbose, "B1")

  if (file.exists(f_path) & !overwrite) {
    stop("The h5ad file already exists and overwrite = FALSE.")
  } else if (file.exists(f_path)) {
    file.remove(f_path)
  }

  rhdf5::h5createFile(f_path)

  if (.verbose) {
    message("Writing the counts to h5ad.")
  }
  # Write X (sparse matrix)
  rhdf5::h5createGroup(f_path, "X")

  if (inherits(counts, "dgCMatrix")) {
    # CSC format
    rhdf5::h5write(counts@x, f_path, "X/data")
    rhdf5::h5write(counts@i, f_path, "X/indices")
    rhdf5::h5write(counts@p, f_path, "X/indptr")
  } else if (inherits(counts, "dgRMatrix")) {
    # CSR format
    rhdf5::h5write(counts@x, f_path, "X/data")
    rhdf5::h5write(counts@j, f_path, "X/indices")
    rhdf5::h5write(counts@p, f_path, "X/indptr")
  }

  if (.verbose) {
    message("Writing the obs to h5ad.")
  }

  rhdf5::h5createGroup(f_path, "obs")
  rhdf5::h5write(obs[[1]], f_path, "obs/_index")
  for (col in names(obs)[-1]) {
    rhdf5::h5write(obs[[col]], f_path, paste0("obs/", col))
  }

  if (.verbose) {
    message("Writing the var to h5ad.")
  }

  rhdf5::h5createGroup(f_path, "var")
  rhdf5::h5write(var[[1]], f_path, "var/_index")
  for (col in names(var)[-1]) {
    rhdf5::h5write(var[[col]], f_path, paste0("var/", col))
  }

  rhdf5::h5closeAll()

  invisible()
}
