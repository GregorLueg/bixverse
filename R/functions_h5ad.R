# h5ad functions and helpers ---------------------------------------------------

#' Detect which slot holds raw integer counts in an h5ad file
#'
#' Samples the non-zero values of each candidate slot and returns the first one
#' that is integer-valued. Slots are checked in the given order, so dedicated
#' count slots take precedence over `X`.
#'
#' @param f_path File path to the `.h5ad` file.
#' @param candidates Character vector of slots to test, any of "layers.counts",
#' "raw.X", "X". Order defines priority.
#' @param n_sample Number of values to sample per slot.
#' @param threshold Minimum fraction of non-zero values that must be whole
#' numbers for a slot to count as raw.
#'
#' @returns The detected slot name, or NULL if none qualifies.
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#' # the synthetic writer only fills /X, and it holds raw counts
#' data <- generate_single_cell_test_data(
#'   syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
#' )
#' f_path <- tempfile(fileext = ".h5ad")
#' write_h5ad_sc(f_path, data$counts, data$obs, data$var, .verbose = FALSE)
#' detect_raw_count_slot(f_path)
#'
#' unlink(f_path)
detect_raw_count_slot <- function(
  f_path,
  candidates = c("layers.counts", "raw.X", "X"),
  n_sample = 10000L,
  threshold = 0.99
) {
  checkmate::assertFileExists(f_path)
  candidates <- match.arg(
    candidates,
    c("layers.counts", "raw.X", "X"),
    several.ok = TRUE
  )
  checkmate::qassert(n_sample, "X1[1,)")
  checkmate::qassert(threshold, "N1[0,1]")

  rs_h5ad_raw_slot(
    f_path = f_path,
    candidates = candidates,
    n_sample = as.integer(n_sample),
    threshold = threshold
  )
}

## scan multiple files for ingestion -------------------------------------------

#' Pre-scan multiple h5ad files for multi-sample loading
#'
#' @param h5_paths Character vector of file paths to h5ad files. If names are
#' provided, these will be used as experimental identifiers.
#' @param gene_universe One of "intersection" or "union".
#' @param var_index String. The name within the h5ad var part in which the
#' variable names are stored. Defaults to `"_index"`.
#' @param raw_count_slot Where raw counts live. `"auto"` detects per file via
#' [detect_raw_count_slot()]; otherwise one of `"X"`, `"raw.X"`,
#' `"layers.counts"`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A list with:
#' \itemize{
#'  \item universe - Character vector of gene names in the universe
#'  \item file_tasks - List of per-file task structures, each containing:
#'  exp_id, h5_path, cs_type, no_cells, no_genes, gene_local_to_universe
#'  (integer vector, NA for genes not in universe, 0-indexed into universe)
#' }
#'
#' @export
#'
#' @examples
#' # build the gene universe across two files before a multi-file load
#' data <- generate_single_cell_test_data(
#'   syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
#' )
#' files <- c(a = tempfile(fileext = ".h5ad"), b = tempfile(fileext = ".h5ad"))
#' for (f in files) {
#'   write_h5ad_sc(f, data$counts, data$obs, data$var, .verbose = FALSE)
#' }
#' tasks <- prescan_h5ad_files(h5_paths = files, .verbose = FALSE)
#' tasks$universe_size
#'
#' unlink(files)
prescan_h5ad_files <- function(
  h5_paths,
  gene_universe = c("intersection", "union"),
  var_index = "_index",
  raw_count_slot = c("auto", "X", "raw.X", "layers.counts"),
  .verbose = TRUE
) {
  gene_universe <- match.arg(gene_universe)
  raw_count_slot <- match.arg(raw_count_slot)
  checkmate::assertCharacter(h5_paths, min.len = 2L)
  invisible(lapply(h5_paths, checkmate::assertFileExists))
  checkmate::qassert(var_index, "S1")
  checkmate::qassert(.verbose, "B1")
  checkmate::assertChoice(
    raw_count_slot,
    c("auto", "X", "raw.X", "layers.counts")
  )

  exp_ids <- if (is.null(names(h5_paths))) {
    tools::file_path_sans_ext(basename(h5_paths))
  } else {
    names(h5_paths)
  }
  checkmate::assertCharacter(exp_ids, len = length(h5_paths), unique = TRUE)

  h5_paths <- path.expand(h5_paths)

  file_meta <- vector("list", length(h5_paths))

  if (.verbose) {
    cli::cli_progress_bar("Scanning h5ad files", total = length(h5_paths))
  }

  for (i in seq_along(h5_paths)) {
    if (.verbose) {
      cli::cli_progress_update(status = exp_ids[[i]])
    }

    meta <- get_h5ad_dimensions(h5_paths[[i]])
    gene_names <- as.vector(
      rs_h5_read(h5_paths[[i]], sprintf("var/%s", var_index))
    )

    resolved_slot <- if (raw_count_slot == "auto") {
      detect_raw_count_slot(h5_paths[[i]])
    } else {
      raw_count_slot
    }

    if (is.null(resolved_slot)) {
      warning(sprintf(
        "Could not detect a raw count slot for %s",
        h5_paths[[i]]
      ))
      resolved_slot <- NA_character_
    }

    file_meta[[i]] <- list(
      exp_id = exp_ids[[i]],
      h5_path = h5_paths[[i]],
      cs_type = meta$type,
      no_cells = meta$dims[["obs"]],
      no_genes = meta$dims[["var"]],
      gene_names = gene_names,
      raw_slot = resolved_slot
    )
  }

  if (.verbose) {
    cli::cli_progress_done()
  }

  all_gene_sets <- lapply(file_meta, `[[`, "gene_names")

  universe <- if (gene_universe == "intersection") {
    Reduce(intersect, all_gene_sets)
  } else {
    Reduce(union, all_gene_sets)
  }

  universe <- sort(universe, method = "radix")
  universe_lookup <- setNames(seq_along(universe) - 1L, universe)

  file_tasks <- lapply(file_meta, function(fm) {
    mapping <- universe_lookup[fm$gene_names]
    mapping <- as.integer(unname(mapping))

    list(
      exp_id = fm$exp_id,
      h5_path = fm$h5_path,
      cs_type = fm$cs_type,
      no_cells = fm$no_cells,
      no_genes = fm$no_genes,
      raw_slot = fm$raw_slot,
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

#' Helper function to get the dimensions and storage format
#'
#' Distinguishes sparse (`CSR`/`CSC`) from dense storage. For dense `X` the
#' orientation is inferred by comparing the matrix dims against the obs and
#' var lengths (`DENSE_ROW` = cells x genes, `DENSE_COL` = genes x cells).
#' Ties (`no_obs == no_var`) fall back to the AnnData convention (`DENSE_ROW`).
#'
#' @param f_path File path to the `.h5ad` file.
#'
#' @returns A list with `dims` (named integer `c(obs, var)`) and `type` (one of
#'   `"CSR"`, `"CSC"`, `"DENSE_ROW"`, `"DENSE_COL"`).
#'
#' @export
#'
#' @examples
#' # dimensions and storage layout without reading the counts
#' data <- generate_single_cell_test_data(
#'   syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
#' )
#' f_path <- tempfile(fileext = ".h5ad")
#' write_h5ad_sc(f_path, data$counts, data$obs, data$var, .verbose = FALSE)
#' get_h5ad_dimensions(f_path)
#'
#' unlink(f_path)
get_h5ad_dimensions <- function(f_path) {
  checkmate::assertFileExists(f_path)

  res <- rs_h5ad_dimensions(f_path = f_path)

  list(
    dims = c(obs = as.integer(res$n_obs), var = as.integer(res$n_var)),
    type = res$format
  )
}

### eda ------------------------------------------------------------------------

#' Read obs and var tables and metadata from an h5ad file
#'
#' @description
#' Useful for exploring the data stored in an h5ad file.
#'
#' @param f_path File path to the `.h5ad` file.
#'
#' @returns A list with:
#' \itemize{
#'   \item obs - data.table of cell-level metadata
#'   \item var - data.table of gene-level metadata
#'   \item dims - named integer vector c(obs, var)
#'   \item type - "CSR" or "CSC"
#' }
#'
#' @export
#'
#' @examples
#' # obs and var tables straight out of the file
#' data <- generate_single_cell_test_data(
#'   syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
#' )
#' f_path <- tempfile(fileext = ".h5ad")
#' write_h5ad_sc(f_path, data$counts, data$obs, data$var, .verbose = FALSE)
#' meta <- read_h5ad_metadata(f_path)
#' head(meta$var, 3)
#'
#' unlink(f_path)
read_h5ad_metadata <- function(f_path) {
  checkmate::assertFileExists(f_path)

  meta <- get_h5ad_dimensions(f_path)

  .read_group_as_dt <- function(group_path) {
    frame <- rs_h5ad_frame(f_path, group_path)
    dt <- data.table::as.data.table(frame$cols)

    if (!is.null(frame$idx)) {
      dt[, .id := frame$idx]
    } else {
      n <- if (length(frame$cols) > 0L) length(frame$cols[[1L]]) else 0L
      dt[, .id := seq_len(n)]
    }

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
#' @returns A list with:
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
#'
#' @examples
#' # check whether /X holds raw counts before loading it
#' data <- generate_single_cell_test_data(
#'   syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
#' )
#' f_path <- tempfile(fileext = ".h5ad")
#' write_h5ad_sc(f_path, data$counts, data$obs, data$var, .verbose = FALSE)
#' summary_x <- read_h5ad_x_summary(f_path)
#' summary_x$stats
#' summary_x$is_integer_valued
#'
#' unlink(f_path)
read_h5ad_x_summary <- function(f_path, n_sample = 10000L) {
  checkmate::assertFileExists(f_path)

  meta <- get_h5ad_dimensions(f_path)

  vals <- rs_h5_read(f_path, "X/data")

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

# h5 10x input functions and helpers -------------------------------------------

## helpers ---------------------------------------------------------------------

#' Detect 10x h5 version and read dimensions
#'
#' @param f_path Path to the 10x CellRanger h5 file.
#'
#' @returns A list with `version` (`"v2"`/`"v3"`), `n_cells` and `n_genes`.
#'
#' @keywords internal
get_tenx_h5_metadata <- function(f_path) {
  checkmate::assertFileExists(f_path)

  rs_tenx_h5_metadata(f_path)
}

#' Detect likely isotype-control features by name pattern
#'
#' @param feature_names Character. ADT feature names (matrix colnames).
#' @param pattern String. Case-insensitive regex. Defaults to `"isotype"`.
#'
#' @returns Character vector of matching names, for inspection before passing
#' to [add_adt_counts_sc()] as `isotype_names`.
#'
#' @export
#'
#' @examples
#' # find the isotype controls among the ADT features
#' features <- c("CD3", "CD19", "IgG1_isotype", "IgG2a_isotype")
#' detect_adt_isotypes(features)
detect_adt_isotypes <- function(feature_names, pattern = "isotype") {
  checkmate::qassert(feature_names, "S+")
  checkmate::qassert(pattern, "S1")

  grep(pattern, feature_names, ignore.case = TRUE, value = TRUE)
}

#' Return the ADT feature names removing the isotypes
#'
#' @param feature_names Character. ADT feature names (matrix colnames or
#' from [get_adt_names()]).
#' @param pattern String. Case-insensitive regex. Defaults to `"isotype"`.
#'
#' @returns Character vector of ADT features, but anything with `"isotype"`.
#'
#' @export
#'
#' @examples
#' # the same features with the isotype controls dropped
#' features <- c("CD3", "CD19", "IgG1_isotype", "IgG2a_isotype")
#' remove_adt_isotypes(features)
remove_adt_isotypes <- function(feature_names, pattern = "isotype") {
  checkmate::qassert(feature_names, "S+")
  checkmate::qassert(pattern, "S1")

  feature_names[!grepl(pattern, feature_names, ignore.case = TRUE)]
}

## multi h5 files --------------------------------------------------------------

#' Pre-scan multiple 10x CellRanger h5 files for multi-sample loading
#'
#' @description
#' Walks each input file, reads the per-file feature ids, restricts to the
#' chosen `feature_type` (V3 only; ignored for V2), builds the intersection
#' or union universe of gene ids, and returns the file tasks expected by
#' [bixverse::load_multi_tenx_h5()].
#'
#' V2 and V3 files can be mixed in a single batch. Gene matching is done on
#' the feature `id` (ensembl-style) since gene symbols collide.
#'
#' @param h5_paths Character vector of file paths to 10x h5 files. If names
#' are provided, these will be used as experimental identifiers.
#' @param feature_type String. Modality to keep across all files (V3 only;
#' ignored for V2). Defaults to `"Gene Expression"`.
#' @param gene_universe One of `"intersection"` or `"union"`.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @returns A list with:
#' \itemize{
#'   \item universe - Character vector of gene ids in the universe.
#'   \item universe_size - Length of the universe.
#'   \item file_tasks - Named list of per-file task structures, each
#'   containing `exp_id`, `h5_path`, `version`, `no_cells`, `no_genes`,
#'   `feature_type` and `gene_local_to_universe` (integer vector, `NA`
#'   for features outside the universe / non-target modality, 0-indexed).
#' }
#'
#' @export
#'
#' @examples
#' # gene universe across two 10x h5 files
#' data <- generate_single_cell_test_data(
#'   syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
#' )
#' features <- data.table::data.table(
#'   id = data$var$gene_id,
#'   name = data$var$ensembl_id,
#'   feature_type = "Gene Expression"
#' )
#' files <- c(a = tempfile(fileext = ".h5"), b = tempfile(fileext = ".h5"))
#' for (f in files) {
#'   write_tenx_h5_sc(f, data$counts, data$obs$cell_id, features)
#' }
#' scan_res <- prescan_tenx_h5_files(h5_paths = files, .verbose = FALSE)
#' scan_res$universe_size
#'
#' unlink(files)
prescan_tenx_h5_files <- function(
  h5_paths,
  feature_type = "Gene Expression",
  gene_universe = c("intersection", "union"),
  .verbose = TRUE
) {
  gene_universe <- match.arg(gene_universe)
  checkmate::assertCharacter(h5_paths, min.len = 2L)
  invisible(lapply(h5_paths, checkmate::assertFileExists))
  checkmate::qassert(feature_type, "S1")
  checkmate::qassert(.verbose, "B1")

  exp_ids <- if (is.null(names(h5_paths))) {
    tools::file_path_sans_ext(basename(h5_paths))
  } else {
    names(h5_paths)
  }
  checkmate::assertCharacter(exp_ids, len = length(h5_paths), unique = TRUE)

  h5_paths <- path.expand(h5_paths)

  file_meta <- vector("list", length(h5_paths))

  if (.verbose) {
    cli::cli_progress_bar("Scanning 10x h5 files", total = length(h5_paths))
  }

  for (i in seq_along(h5_paths)) {
    if (.verbose) {
      cli::cli_progress_update(status = exp_ids[[i]])
    }

    meta <- get_tenx_h5_metadata(h5_paths[[i]])

    if (meta$version == "v3") {
      all_ids <- rs_h5_read(h5_paths[[i]], "matrix/features/id")
      all_types <- rs_h5_read(h5_paths[[i]], "matrix/features/feature_type")
      keep_mask <- all_types == feature_type

      if (!any(keep_mask)) {
        stop(sprintf(
          "No features of type '%s' found in %s. Available: %s",
          feature_type,
          h5_paths[[i]],
          paste(unique(all_types), collapse = ", ")
        ))
      }
    } else {
      all_ids <- rs_h5_read(h5_paths[[i]], "genes")
      keep_mask <- rep(TRUE, length(all_ids))
    }

    file_meta[[i]] <- list(
      exp_id = exp_ids[[i]],
      h5_path = h5_paths[[i]],
      version = meta$version,
      no_cells = meta$n_cells,
      no_genes = meta$n_genes,
      all_ids = all_ids,
      keep_mask = keep_mask
    )
  }

  if (.verbose) {
    cli::cli_progress_done()
  }

  per_file_target_ids <- lapply(
    file_meta,
    function(fm) fm$all_ids[fm$keep_mask]
  )

  universe <- if (gene_universe == "intersection") {
    Reduce(intersect, per_file_target_ids)
  } else {
    Reduce(union, per_file_target_ids)
  }

  if (length(universe) == 0L) {
    stop("Gene universe across inputs is empty.")
  }

  universe <- sort(universe, method = "radix")
  universe_lookup <- setNames(seq_along(universe) - 1L, universe)

  file_tasks <- lapply(file_meta, function(fm) {
    mapping <- rep(NA_integer_, fm$no_genes)
    in_target <- fm$keep_mask
    mapping[in_target] <- as.integer(unname(
      universe_lookup[fm$all_ids[in_target]]
    ))

    list(
      exp_id = fm$exp_id,
      h5_path = fm$h5_path,
      version = fm$version,
      no_cells = fm$no_cells,
      no_genes = fm$no_genes,
      feature_type = feature_type,
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

## eda -------------------------------------------------------------------------

#' Read barcode and feature tables and metadata from a 10x h5 file
#'
#' @description
#' Useful for exploring the data stored in a 10x CellRanger v2/v3 h5 file,
#' including the breakdown of feature types in a multi-modal file.
#'
#' @param f_path File path to the 10x `.h5` file.
#'
#' @returns A list with:
#' \itemize{
#'   \item obs - data.table of barcodes
#'   \item var - data.table of features (id, name, and feature_type for v3)
#'   \item dims - named integer vector c(obs, var)
#'   \item version - "v2" or "v3"
#'   \item feature_types - named integer vector of feature_type counts (v3
#'   only, otherwise NULL)
#' }
#'
#' @export
#'
#' @examples
#' # barcodes, features and the feature type breakdown
#' data <- generate_single_cell_test_data(
#'   syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
#' )
#' f_path <- tempfile(fileext = ".h5")
#' write_tenx_h5_sc(
#'   f_path = f_path,
#'   counts = data$counts,
#'   barcodes = data$obs$cell_id,
#'   features = data.table::data.table(
#'     id = data$var$gene_id,
#'     name = data$var$ensembl_id,
#'     feature_type = "Gene Expression"
#'   )
#' )
#' meta <- read_tenx_h5_metadata(f_path)
#' meta$dims
#' meta$feature_types
#'
#' unlink(f_path)
read_tenx_h5_metadata <- function(f_path) {
  checkmate::assertFileExists(f_path)

  meta <- get_tenx_h5_metadata(f_path)
  version <- meta$version

  barcodes_path <- if (version == "v3") "matrix/barcodes" else "barcodes"
  obs <- data.table::data.table(.id = rs_h5_read(f_path, barcodes_path))

  var <- if (version == "v3") {
    data.table::data.table(
      .id = rs_h5_read(f_path, "matrix/features/id"),
      name = rs_h5_read(f_path, "matrix/features/name"),
      feature_type = rs_h5_read(f_path, "matrix/features/feature_type")
    )
  } else {
    data.table::data.table(
      .id = rs_h5_read(f_path, "genes"),
      name = rs_h5_read(f_path, "gene_names")
    )
  }

  feature_types <- if (version == "v3") {
    table(var$feature_type) |> (\(x) setNames(as.integer(x), names(x)))()
  } else {
    NULL
  }

  list(
    obs = obs,
    var = var,
    dims = c(obs = meta$n_cells, var = meta$n_genes),
    version = version,
    feature_types = feature_types
  )
}

## i/o -------------------------------------------------------------------------

### single file ----------------------------------------------------------------

#' Read in 10x h5 ADT data
#'
#' @description
#' Helper function to load in ADT counts from h5. Leverages Rust under the
#' hood for faster filtering and reading.
#'
#' @param f_path String. File to the h5 file from which to read the ADT counts.
#' @param feature_type String. The feature type to return. Defaults here to
#' `"Antibody Capture"`
#'
#' @returns A dense matrix of cells x features
#'
#' @export
#'
#' @examples
#' # pull the antibody capture layer out of a multi-modal 10x file
#' data <- generate_single_cell_test_data(
#'   syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
#' )
#' adt <- generate_single_cell_test_data_adt(
#'   params_sc_synthetic_data_adt(n_cells = 200L)
#' )
#' f_path <- tempfile(fileext = ".h5")
#' write_tenx_h5_sc(
#'   f_path = f_path,
#'   counts = cbind(data$counts, as(adt$counts, "RsparseMatrix")),
#'   barcodes = data$obs$cell_id,
#'   features = data.table::data.table(
#'     id = c(data$var$gene_id, colnames(adt$counts)),
#'     name = c(data$var$ensembl_id, colnames(adt$counts)),
#'     feature_type = rep(
#'       c("Gene Expression", "Antibody Capture"),
#'       c(ncol(data$counts), ncol(adt$counts))
#'     )
#'   )
#' )
#' adt_counts <- read_tenx_h5_adt(f_path)
#' dim(adt_counts)
#'
#' unlink(f_path)
read_tenx_h5_adt <- function(f_path, feature_type = "Antibody Capture") {
  # checks
  checkmate::assertFileExists(f_path)

  # function
  meta <- get_tenx_h5_metadata(path.expand(f_path))

  res <- rs_read_tenx_h5_modality(
    f_path = path.expand(f_path),
    version = meta$version,
    feature_type = feature_type
  )

  m <- res$counts
  rownames(m) <- res$barcodes
  colnames(m) <- res$features
  m
}

### multiple files -------------------------------------------------------------

#' Read in 10x h5 ADT data from multiple files
#'
#' @description
#' Multi-file counterpart to [bixverse::read_tenx_h5_adt()]. Reads the same
#' modality from each input, stacks the cells (rows), and prefixes each
#' barcode with its `exp_id` so the result matches the cell_id convention
#' used by [bixverse::load_multi_tenx_h5()]. The feature space is either
#' the intersection or union of features across inputs; missing features
#' in the union case are filled with zero.
#'
#' @param h5_paths Character vector of file paths to 10x h5 files. If names
#' are provided, these will be used as `exp_id`s; otherwise the file basename
#' is used.
#' @param feature_type String. The feature type to return. Defaults to
#' `"Antibody Capture"`.
#' @param gene_universe One of `"intersection"` or `"union"`.
#'
#' @returns A dense matrix of cells x features with `exp_id_barcode`
#' rownames and feature names as colnames.
#'
#' @export
#'
#' @examples
#' # stack the ADT layer of two files, barcodes prefixed by exp_id
#' data <- generate_single_cell_test_data(
#'   syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
#' )
#' adt <- generate_single_cell_test_data_adt(
#'   params_sc_synthetic_data_adt(n_cells = 200L)
#' )
#' features <- data.table::data.table(
#'   id = c(data$var$gene_id, colnames(adt$counts)),
#'   name = c(data$var$ensembl_id, colnames(adt$counts)),
#'   feature_type = rep(
#'     c("Gene Expression", "Antibody Capture"),
#'     c(ncol(data$counts), ncol(adt$counts))
#'   )
#' )
#' counts <- cbind(data$counts, as(adt$counts, "RsparseMatrix"))
#' files <- c(a = tempfile(fileext = ".h5"), b = tempfile(fileext = ".h5"))
#' for (f in files) {
#'   write_tenx_h5_sc(f, counts, data$obs$cell_id, features)
#' }
#' adt_counts <- read_multi_tenx_h5_adt(files)
#' dim(adt_counts)
#' head(rownames(adt_counts), 2)
#'
#' unlink(files)
read_multi_tenx_h5_adt <- function(
  h5_paths,
  feature_type = "Antibody Capture",
  gene_universe = c("intersection", "union")
) {
  gene_universe <- match.arg(gene_universe)
  checkmate::assertCharacter(h5_paths, min.len = 2L)
  invisible(lapply(h5_paths, checkmate::assertFileExists))
  checkmate::qassert(feature_type, "S1")

  exp_ids <- if (is.null(names(h5_paths))) {
    tools::file_path_sans_ext(basename(h5_paths))
  } else {
    names(h5_paths)
  }
  checkmate::assertCharacter(exp_ids, len = length(h5_paths), unique = TRUE)

  h5_paths <- path.expand(h5_paths)

  per_file <- vector("list", length(h5_paths))

  for (i in seq_along(h5_paths)) {
    meta <- get_tenx_h5_metadata(h5_paths[[i]])
    res <- rs_read_tenx_h5_modality(
      f_path = h5_paths[[i]],
      version = meta$version,
      feature_type = feature_type
    )
    m <- res$counts
    rownames(m) <- paste(exp_ids[[i]], res$barcodes, sep = "_")
    colnames(m) <- res$features
    per_file[[i]] <- m
  }

  feature_sets <- lapply(per_file, colnames)
  features <- if (gene_universe == "intersection") {
    Reduce(intersect, feature_sets)
  } else {
    Reduce(union, feature_sets)
  }

  if (length(features) == 0L) {
    stop("Feature universe across inputs is empty.")
  }

  total_rows <- sum(vapply(per_file, nrow, integer(1L)))
  combined <- matrix(0, nrow = total_rows, ncol = length(features))
  colnames(combined) <- features
  row_names <- character(total_rows)

  offset <- 0L
  for (m in per_file) {
    n <- nrow(m)
    rng <- (offset + 1L):(offset + n)
    in_common <- intersect(colnames(m), features)
    combined[rng, in_common] <- as.matrix(m[, in_common, drop = FALSE])
    row_names[rng] <- rownames(m)
    offset <- offset + n
  }
  rownames(combined) <- row_names

  combined
}
