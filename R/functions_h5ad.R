# h5ad functions and helpers ---------------------------------------------------

## scan multiple files for ingestion -------------------------------------------

#' Pre-scan multiple h5ad files for multi-sample loading
#'
#' @param h5_paths Character vector of file paths to h5ad files. If names are
#' provided, these will be used as experimental identifiers.
#' @param gene_universe One of "intersection" or "union".
#' @param var_index String. The name within the h5ad var part in which the
#' variable names are stored. Defaults to `"_index"`.
#'
#' @return A list with:
#' \itemize{
#'  \item universe - Character vector of gene names in the universe
#'  \item file_tasks - List of per-file task structures, each containing:
#'  exp_id, h5_path, cs_type, no_cells, no_genes, gene_local_to_universe
#'  (integer vector, NA for genes not in universe, 0-indexed into universe)
#' }
#'
#' @export
prescan_h5ad_files <- function(
  h5_paths,
  gene_universe = c("intersection", "union"),
  var_index = "_index"
) {
  # checks
  gene_universe <- match.arg(gene_universe)
  checkmate::assertCharacter(h5_paths, min.len = 2L)
  invisible(lapply(h5_paths, checkmate::assertFileExists))
  checkmate::qassert(var_index, "S1")

  exp_ids <- if (is.null(names(h5_paths))) {
    tools::file_path_sans_ext(basename(h5_paths))
  } else {
    names(h5_paths)
  }
  checkmate::assertCharacter(exp_ids, len = length(h5_paths), unique = TRUE)

  h5_paths <- path.expand(h5_paths)

  # collect per-file metadata
  file_meta <- vector("list", length(h5_paths))

  for (i in seq_along(h5_paths)) {
    meta <- get_h5ad_dimensions(h5_paths[[i]])
    gene_names <- as.vector(
      rhdf5::h5read(h5_paths[[i]], sprintf("var/%s", var_index))
    )
    rhdf5::h5closeAll()

    file_meta[[i]] <- list(
      exp_id = exp_ids[[i]],
      h5_path = h5_paths[[i]],
      cs_type = meta$type,
      no_cells = meta$dims[["obs"]],
      no_genes = meta$dims[["var"]],
      gene_names = gene_names
    )
  }

  # compute gene universe
  all_gene_sets <- lapply(file_meta, `[[`, "gene_names")

  universe <- if (gene_universe == "intersection") {
    Reduce(intersect, all_gene_sets)
  } else {
    Reduce(union, all_gene_sets)
  }

  # preserve a stable ordering (sorted)
  universe <- sort(universe)
  universe_lookup <- setNames(seq_along(universe) - 1L, universe)

  # build per-file mappings
  file_tasks <- lapply(file_meta, function(fm) {
    mapping <- universe_lookup[fm$gene_names]
    # genes not in universe become NA
    mapping <- as.integer(unname(mapping))

    list(
      exp_id = fm$exp_id,
      h5_path = fm$h5_path,
      cs_type = fm$cs_type,
      no_cells = fm$no_cells,
      no_genes = fm$no_genes,
      gene_local_to_universe = mapping
    )
  })

  names(file_tasks) <- exp_ids

  list(
    universe = universe,
    universe_size = length(universe),
    file_tasks = file_tasks
  )
}

## single h5ad files -----------------------------------------------------------

### dimensions -----------------------------------------------------------------

#' Helper function to get the dimensions and compressed sparse format
#'
#' @param f_path File path to the `.h5ad` file.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item dims - Dimensions of the stored data in the h5ad file.
#'   \item type - Was the data stored in CSR (indptr = cells) or CSC (indptr =
#'   genes).
#' }
get_h5ad_dimensions <- function(f_path) {
  # checks
  checkmate::assertFileExists(f_path)

  # function
  h5_content <- rhdf5::h5ls(
    f_path
  ) %>%
    data.table::setDT()

  on.exit(tryCatch(rhdf5::h5closeAll(), error = function(e) invisible()))

  no_obs <- h5_content[
    group == "/obs" & otype == "H5I_DATASET"
  ][1, as.numeric(dim)]

  no_var <- h5_content[
    group == "/var" & otype == "H5I_DATASET"
  ][1, as.numeric(dim)]

  indptr <- h5_content[
    group == "/X" & name == "indptr",
    as.numeric(dim)
  ]

  cs_format <- ifelse(no_var + 1 == indptr, "CSC", "CSR")

  return(list(
    dims = setNames(c(as.integer(no_obs), as.integer(no_var)), c("obs", "var")),
    type = cs_format
  ))
}

### eda ------------------------------------------------------------------------

#' Read obs and var tables and metadata from an h5ad file
#'
#' @description
#' Useful for exploring the data stored in an h5ad file.
#'
#' @param f_path File path to the `.h5ad` file.
#'
#' @return A list with:
#' \itemize{
#'   \item obs - data.table of cell-level metadata
#'   \item var - data.table of gene-level metadata
#'   \item dims - named integer vector c(obs, var)
#'   \item type - "CSR" or "CSC"
#' }
#'
#' @export
read_h5ad_metadata <- function(f_path) {
  checkmate::assertFileExists(f_path)
  on.exit(tryCatch(rhdf5::h5closeAll(), error = function(e) invisible()))

  h5_content <- rhdf5::h5ls(f_path) |> data.table::setDT()
  meta <- get_h5ad_dimensions(f_path)

  .read_group_as_dt <- function(group_path) {
    entries <- h5_content[group == group_path]
    groups <- h5_content[group == group_path & otype == "H5I_GROUP", name]

    cols <- list()

    # categorical columns (stored as categories + codes)
    for (g in groups) {
      sub_path <- paste0(group_path, "/", g)
      cats <- rhdf5::h5read(f_path, paste0(sub_path, "/categories"))
      codes <- rhdf5::h5read(f_path, paste0(sub_path, "/codes")) + 1L
      cols[[g]] <- factor(cats[codes], levels = cats)
    }

    # direct datasets
    direct <- entries[otype == "H5I_DATASET" & name != "_index", name]
    for (d in direct) {
      cols[[d]] <- as.vector(rhdf5::h5read(f_path, paste0(group_path, "/", d)))
    }

    idx <- as.vector(rhdf5::h5read(f_path, paste0(group_path, "/_index")))
    dt <- data.table::as.data.table(cols)
    dt[, .id := idx]
    data.table::setcolorder(dt, c(".id", setdiff(names(dt), ".id")))
    dt
  }

  list(
    obs = .read_group_as_dt("/obs"),
    var = .read_group_as_dt("/var"),
    dims = meta$dims,
    type = meta$type
  )
}

#' Read summary statistics from the X slot of an h5ad file
#'
#' @param f_path File path to the `.h5ad` file.
#' @param n_sample Number of non-zero values to sample for the preview. NULL
#' reads all.
#'
#' @return A list with:
#' \itemize{
#'   \item stats - named vector: min, max, mean, median, and fraction of values
#'   that are whole numbers
#'   \item is_integer_valued - logical; TRUE if >99% of sampled non-zero values
#'   are whole numbers
#'   \item type - "CSR" or "CSC"
#'   \item dims - named integer vector c(obs, var)
#'   \item sample - numeric vector of sampled non-zero values
#' }
#'
#' @export
read_h5ad_x_summary <- function(f_path, n_sample = 10000L) {
  checkmate::assertFileExists(f_path)
  on.exit(tryCatch(rhdf5::h5closeAll(), error = function(e) invisible()))

  meta <- get_h5ad_dimensions(f_path)

  vals <- rhdf5::h5read(f_path, "/X/data")

  if (!is.null(n_sample) && length(vals) > n_sample) {
    vals_sample <- sample(vals, n_sample)
  } else {
    vals_sample <- vals
  }

  whole_number_frac <- mean(vals_sample == floor(vals_sample))

  list(
    stats = c(
      min = min(vals_sample),
      max = max(vals_sample),
      mean = mean(vals_sample),
      median = median(vals_sample),
      whole_number_frac = whole_number_frac
    ),
    is_integer_valued = whole_number_frac > 0.99,
    type = meta$type,
    dims = meta$dims,
    sample = vals_sample
  )
}
