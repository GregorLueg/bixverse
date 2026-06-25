# h5ad functions and helpers ---------------------------------------------------

## helpers ---------------------------------------------------------------------

#' Resolve the index column from an h5ad group
#'
#' Checks for the `_index` attribute first (new anndata spec), then falls back
#' to the `_index` dataset (old spec). Returns NULL if neither exists.
#'
#' @param f_path Path to the h5ad file.
#' @param group_path The HDF5 group (e.g. "/obs", "/var").
#' @param h5_content data.table from rhdf5::h5ls().
#'
#' @return A list with `idx` (character vector or NULL) and `idx_col` (the
#'   dataset name used, or NULL).
#'
#' @keywords internal
.resolve_h5_index <- function(f_path, group_path, h5_content) {
  grp_attrs <- tryCatch(
    rhdf5::h5readAttributes(f_path, group_path),
    error = function(e) list()
  )

  if ("_index" %in% names(grp_attrs)) {
    idx_col <- grp_attrs[["_index"]]
    idx <- as.vector(rhdf5::h5read(f_path, paste0(group_path, "/", idx_col)))
    return(list(idx = idx, idx_col = idx_col))
  }

  has_ds <- "_index" %in%
    h5_content[
      group == group_path & otype == "H5I_DATASET",
      name
    ]

  if (has_ds) {
    idx <- as.vector(rhdf5::h5read(f_path, paste0(group_path, "/_index")))
    return(list(idx = idx, idx_col = "_index"))
  }

  list(idx = NULL, idx_col = NULL)
}

#' Read a sample of values from a matrix slot in an h5ad file
#'
#' Handles both sparse (group with a `data` dataset) and dense (direct dataset)
#' storage. Returns NULL if the slot does not exist.
#'
#' @keywords internal
.read_slot_value_sample <- function(
  f_path,
  slot_path,
  h5_content,
  n_sample = 10000L
) {
  full <- data.table::fifelse(
    h5_content$group == "/",
    paste0("/", h5_content$name),
    paste0(h5_content$group, "/", h5_content$name)
  )

  if (!slot_path %in% full) {
    return(NULL)
  }

  otype <- h5_content$otype[full == slot_path][1L]

  if (otype == "H5I_GROUP") {
    data_path <- paste0(slot_path, "/data")
    if (!data_path %in% full) {
      return(NULL)
    }
    len <- as.numeric(
      strsplit(h5_content$dim[full == data_path][1L], " x ", fixed = TRUE)[[
        1L
      ]][1L]
    )
    idx <- seq_len(min(len, n_sample))
    return(as.vector(rhdf5::h5read(f_path, data_path, index = list(idx))))
  }

  # dense dataset
  dims <- as.integer(
    strsplit(h5_content$dim[full == slot_path][1L], " x ", fixed = TRUE)[[1L]]
  )
  side <- max(1L, as.integer(ceiling(sqrt(n_sample))))
  idx <- lapply(dims, function(d) seq_len(min(d, side)))
  as.vector(rhdf5::h5read(f_path, slot_path, index = idx))
}

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
#' @return The detected slot name, or NULL if none qualifies.
#'
#' @export
#'
#' @keywords internal
detect_raw_count_slot <- function(
  f_path,
  candidates = c("layers.counts", "raw.X", "X"),
  n_sample = 10000L,
  threshold = 0.99
) {
  checkmate::assertFileExists(f_path)
  on.exit(tryCatch(rhdf5::h5closeAll(), error = function(e) invisible()))

  slot_paths <- c(X = "/X", raw.X = "/raw/X", layers.counts = "/layers/counts")
  candidates <- match.arg(candidates, names(slot_paths), several.ok = TRUE)

  h5_content <- rhdf5::h5ls(f_path) |> data.table::setDT()

  for (slot in candidates) {
    vals <- .read_slot_value_sample(
      f_path,
      slot_paths[[slot]],
      h5_content,
      n_sample
    )
    if (is.null(vals) || length(vals) == 0L) {
      next
    }
    nz <- vals[vals != 0]
    if (length(nz) == 0L) {
      next
    }
    if (mean(nz == floor(nz)) >= threshold) return(slot)
  }

  NULL
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
      rhdf5::h5read(h5_paths[[i]], sprintf("var/%s", var_index))
    )
    rhdf5::h5closeAll()

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

  universe <- sort(universe)
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
#' @return A list with `dims` (named integer `c(obs, var)`) and `type` (one of
#'   `"CSR"`, `"CSC"`, `"DENSE_ROW"`, `"DENSE_COL"`).
get_h5ad_dimensions <- function(f_path) {
  checkmate::assertFileExists(f_path)

  h5_content <- rhdf5::h5ls(f_path) %>%
    data.table::setDT()

  on.exit(tryCatch(rhdf5::h5closeAll(), error = function(e) invisible()))

  no_obs <- h5_content[
    group == "/obs" & otype == "H5I_DATASET"
  ][1, as.numeric(dim)]

  no_var <- h5_content[
    group == "/var" & otype == "H5I_DATASET"
  ][1, as.numeric(dim)]

  x_row <- h5_content[group == "/" & name == "X"]

  if (nrow(x_row) == 0L) {
    stop("Could not locate /X in the h5ad file.")
  }

  cs_format <- if (x_row$otype == "H5I_GROUP") {
    indptr <- h5_content[
      group == "/X" & name == "indptr",
      as.numeric(dim)
    ]
    if (length(indptr) == 0L) {
      stop("/X is a group but contains no indptr dataset.")
    }
    ifelse(no_var + 1 == indptr, "CSC", "CSR")
  } else {
    # h5ls dim is R-view; read native order directly
    fid <- rhdf5::H5Fopen(f_path)
    did <- rhdf5::H5Dopen(fid, "X")
    sid <- rhdf5::H5Dget_space(did)
    x_dims <- rhdf5::H5Sget_simple_extent_dims(sid)$size
    rhdf5::H5Sclose(sid)
    rhdf5::H5Dclose(did)
    rhdf5::H5Fclose(fid)

    if (length(x_dims) != 2L) {
      stop(sprintf("Unexpected /X rank: %d", length(x_dims)))
    }
    if (x_dims[1] == no_obs && x_dims[2] == no_var) {
      "DENSE_COL"
    } else if (x_dims[1] == no_var && x_dims[2] == no_obs) {
      "DENSE_ROW"
    } else if (no_obs == no_var) {
      "DENSE_ROW"
    } else {
      stop(sprintf(
        "/X dense dims (%d x %d) match neither obs/var (%d / %d) ordering.",
        x_dims[1],
        x_dims[2],
        no_obs,
        no_var
      ))
    }
  }

  list(
    dims = setNames(c(as.integer(no_obs), as.integer(no_var)), c("obs", "var")),
    type = cs_format
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

    new_cat_names <- h5_content[
      group == paste0(group_path, "/__categories") & otype == "H5I_DATASET",
      name
    ]
    has_new_cats <- length(new_cat_names) > 0L

    groups <- entries[
      otype == "H5I_GROUP" & name != "__categories",
      name
    ]

    # resolve index first so we can exclude it from direct datasets
    idx_info <- .resolve_h5_index(f_path, group_path, h5_content)

    cols <- list()

    # old-format categoricals
    for (g in groups) {
      sub_path <- paste0(group_path, "/", g)
      sub_entries <- h5_content[
        group == sub_path & otype == "H5I_DATASET",
        name
      ]
      if (all(c("categories", "codes") %in% sub_entries)) {
        cats <- rhdf5::h5read(f_path, paste0(sub_path, "/categories"))
        codes <- rhdf5::h5read(f_path, paste0(sub_path, "/codes"))
        codes[codes < 0L] <- NA_integer_
        cols[[g]] <- factor(cats[codes + 1L], levels = cats)
      }
    }

    # direct datasets — exclude _index dataset and the attribute-referenced col
    skip <- c("_index", idx_info$idx_col)
    direct <- entries[
      otype == "H5I_DATASET" & !name %in% skip,
      name
    ]

    for (d in direct) {
      raw <- as.vector(rhdf5::h5read(f_path, paste0(group_path, "/", d)))
      if (has_new_cats && d %in% new_cat_names) {
        categories <- as.vector(
          rhdf5::h5read(f_path, paste0(group_path, "/__categories/", d))
        )
        raw[raw < 0L] <- NA_integer_
        cols[[d]] <- factor(categories[raw + 1L], levels = categories)
      } else {
        cols[[d]] <- raw
      }
    }

    dt <- data.table::as.data.table(cols)

    if (!is.null(idx_info$idx)) {
      dt[, .id := idx_info$idx]
    } else {
      n <- if (length(cols) > 0L) length(cols[[1L]]) else 0L
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

# h5 10x input functions and helpers -------------------------------------------

## helpers ---------------------------------------------------------------------

#' Detect 10x h5 version and read dimensions
#'
#' @param f_path Path to the 10x CellRanger h5 file.
#'
#' @return A list with `version` (`"v2"`/`"v3"`), `n_cells` and `n_genes`.
#'
#' @keywords internal
get_tenx_h5_metadata <- function(f_path) {
  checkmate::assertFileExists(f_path)

  content <- rhdf5::h5ls(f_path) |> data.table::setDT()
  on.exit(
    tryCatch(rhdf5::h5closeAll(), error = function(e) invisible()),
    add = TRUE
  )

  is_v3 <- nrow(content[group == "/matrix" & name == "data"]) > 0L
  is_v2 <- nrow(content[group == "/" & name == "data"]) > 0L &&
    nrow(content[group == "/" & name == "genes"]) > 0L

  version <- if (is_v3) {
    "v3"
  } else if (is_v2) {
    "v2"
  } else {
    stop(paste(
      "Could not detect 10x h5 version: neither v3 ('/matrix/data')",
      "nor v2 ('/data' + '/genes') layout found."
    ))
  }

  shape_path <- if (version == "v3") "matrix/shape" else "shape"
  shape <- as.integer(rhdf5::h5read(f_path, shape_path))

  list(version = version, n_genes = shape[1], n_cells = shape[2])
}

#' Detect likely isotype-control features by name pattern
#'
#' @param feature_names Character. ADT feature names (matrix colnames).
#' @param pattern String. Case-insensitive regex. Defaults to `"isotype"`.
#'
#' @return Character vector of matching names, for inspection before passing
#' to [add_adt_counts_sc()] as `isotype_names`.
#'
#' @export
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
#' @return Character vector of ADT features, but anything with `"isotype"`.
#'
#' @export
remove_adt_isotypes <- function(feature_names, pattern = "isotype") {
  checkmate::qassert(feature_names, "S+")
  checkmate::qassert(pattern, "S1")

  feature_names[!grepl(pattern, feature_names)]
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
#' @return A list with:
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
      all_ids <- as.character(
        rhdf5::h5read(h5_paths[[i]], "matrix/features/id")
      )
      all_types <- as.character(
        rhdf5::h5read(h5_paths[[i]], "matrix/features/feature_type")
      )
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
      all_ids <- as.character(rhdf5::h5read(h5_paths[[i]], "genes"))
      keep_mask <- rep(TRUE, length(all_ids))
    }
    rhdf5::h5closeAll()

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

  universe <- sort(universe)
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
#' @return A list with:
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
read_tenx_h5_metadata <- function(f_path) {
  checkmate::assertFileExists(f_path)
  on.exit(tryCatch(rhdf5::h5closeAll(), error = function(e) invisible()))

  meta <- get_tenx_h5_metadata(f_path)
  version <- meta$version

  barcodes_path <- if (version == "v3") "matrix/barcodes" else "barcodes"
  obs <- data.table::data.table(
    .id = as.character(rhdf5::h5read(f_path, barcodes_path))
  )

  var <- if (version == "v3") {
    data.table::data.table(
      .id = as.character(rhdf5::h5read(f_path, "matrix/features/id")),
      name = as.character(rhdf5::h5read(f_path, "matrix/features/name")),
      feature_type = as.character(
        rhdf5::h5read(f_path, "matrix/features/feature_type")
      )
    )
  } else {
    data.table::data.table(
      .id = as.character(rhdf5::h5read(f_path, "genes")),
      name = as.character(rhdf5::h5read(f_path, "gene_names"))
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
