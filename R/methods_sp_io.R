# spatial i/o ------------------------------------------------------------------

## helpers ---------------------------------------------------------------------

### image headers --------------------------------------------------------------

#' Read the pixel dimensions out of a PNG header
#'
#' @description
#' Pulls width and height straight out of the IHDR chunk. Avoids pulling `png`
#' (a Suggests) into the ingest path just to learn how big a file is.
#'
#' @param path String. Path to the PNG.
#'
#' @return Named integer vector with `width` and `height`, or `NULL` if the
#' file does not look like a PNG.
#'
#' @keywords internal
.png_dimensions <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)

  hdr <- readBin(con, "raw", n = 24L)
  if (length(hdr) < 24L) {
    return(NULL)
  }
  png_magic <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  if (!identical(hdr[1:8], png_magic)) {
    return(NULL)
  }
  if (!identical(rawToChar(hdr[13:16]), "IHDR")) {
    return(NULL)
  }

  as_u32 <- \(x) sum(as.integer(x) * c(16777216, 65536, 256, 1))

  c(width = as_u32(hdr[17:20]), height = as_u32(hdr[21:24]))
}

#' Read the pixel dimensions out of a baseline TIFF header
#'
#' @description
#' Walks the first IFD and reads tags 256 (`ImageWidth`) and 257
#' (`ImageLength`). Space Ranger ships `cytassist_image.tiff` and no R TIFF
#' reader is a hard dependency here, but the CytAssist scale factor has to be
#' derived from the dimensions of the file that actually shipped, so the header
#' gets parsed by hand. BigTIFF is not handled and returns `NULL`.
#'
#' @param path String. Path to the TIFF.
#'
#' @return Named integer vector with `width` and `height`, or `NULL` if the
#' header could not be parsed.
#'
#' @keywords internal
.tiff_dimensions <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)

  hdr <- readBin(con, "raw", n = 8L)
  if (length(hdr) < 8L) {
    return(NULL)
  }

  endian <- switch(rawToChar(hdr[1:2]), II = "little", MM = "big", NULL)
  if (is.null(endian)) {
    return(NULL)
  }

  magic <- readBin(
    hdr[3:4],
    "integer",
    size = 2L,
    signed = FALSE,
    endian = endian
  )
  if (magic != 42L) {
    return(NULL)
  }

  ifd_offset <- readBin(hdr[5:8], "integer", size = 4L, endian = endian)
  seek(con, where = ifd_offset, origin = "start")
  n_entries <- readBin(
    con,
    "integer",
    size = 2L,
    signed = FALSE,
    endian = endian
  )
  if (length(n_entries) == 0L || n_entries < 1L) {
    return(NULL)
  }
  entries <- readBin(con, "raw", n = 12L * n_entries)
  if (length(entries) < 12L * n_entries) {
    return(NULL)
  }

  # each entry is tag (2), type (2), count (4), value or offset (4). SHORT
  # values sit left-justified in the value field, so the first two bytes are
  # the value under either byte order.
  dims <- c(width = NA_integer_, height = NA_integer_)
  for (i in seq_len(n_entries)) {
    off <- (i - 1L) * 12L
    tag <- readBin(
      entries[(off + 1L):(off + 2L)],
      "integer",
      size = 2L,
      signed = FALSE,
      endian = endian
    )
    if (!tag %in% c(256L, 257L)) {
      next
    }
    type <- readBin(
      entries[(off + 3L):(off + 4L)],
      "integer",
      size = 2L,
      signed = FALSE,
      endian = endian
    )
    val_raw <- entries[(off + 9L):(off + 12L)]
    value <- if (type == 3L) {
      readBin(
        val_raw[1:2],
        "integer",
        size = 2L,
        signed = FALSE,
        endian = endian
      )
    } else {
      readBin(val_raw, "integer", size = 4L, endian = endian)
    }
    dims[[if (tag == 256L) "width" else "height"]] <- as.integer(value)
  }

  if (anyNA(dims)) NULL else dims
}

#' Read the pixel dimensions of a spatial image
#'
#' @description
#' Dispatches on the file extension. Only PNG and TIFF are supported, which
#' covers everything the Visium scale factor derivation needs.
#'
#' @param path String. Path to the image.
#'
#' @return Named integer vector with `width` and `height`, or `NULL`.
#'
#' @keywords internal
.image_dimensions <- function(path) {
  checkmate::assertFileExists(path)

  switch(
    tolower(tools::file_ext(path)),
    png = .png_dimensions(path),
    tif = ,
    tiff = .tiff_dimensions(path),
    NULL
  )
}

### visium directory scan ------------------------------------------------------

#' Resolve the count matrix directory of a Space Ranger output
#'
#' @param visium_dir String. Space Ranger output directory.
#' @param matrix_type String. One of `c("auto", "raw", "filtered")`.
#'
#' @return String. Path to the chosen `*_feature_bc_matrix` directory.
#'
#' @keywords internal
.visium_matrix_dir <- function(visium_dir, matrix_type) {
  candidates <- c(
    filtered = file.path(visium_dir, "filtered_feature_bc_matrix"),
    raw = file.path(visium_dir, "raw_feature_bc_matrix")
  )
  present <- candidates[dir.exists(candidates)]

  if (length(present) == 0L) {
    stop(sprintf(
      paste(
        "Neither 'filtered_feature_bc_matrix/' nor 'raw_feature_bc_matrix/'",
        "found in %s."
      ),
      visium_dir
    ))
  }

  if (matrix_type == "auto") {
    return(unname(present[1L]))
  }

  if (!matrix_type %in% names(present)) {
    stop(sprintf(
      "No '%s_feature_bc_matrix/' directory in %s.",
      matrix_type,
      visium_dir
    ))
  }

  return(unname(present[[matrix_type]]))
}

#' Find the image files Space Ranger writes into `spatial/`
#'
#' @param spatial_dir String. The `spatial/` directory.
#'
#' @return Named list of absolute paths, keyed by resolution. Keys with no
#' matching file on disk are dropped.
#'
#' @keywords internal
.visium_image_sources <- function(spatial_dir) {
  patterns <- list(
    lowres = "^tissue_lowres_image\\.png$",
    hires = "^tissue_hires_image\\.png$",
    cytassist = "^cytassist_image\\.tiff?$"
  )

  found <- lapply(patterns, \(pat) {
    hits <- list.files(
      spatial_dir,
      pattern = pat,
      full.names = TRUE,
      ignore.case = TRUE
    )
    if (length(hits) == 1L) hits else NULL
  })

  found[!purrr::map_lgl(found, is.null)]
}

#' Read a Visium tissue positions file
#'
#' @description
#' Space Ranger 2.x writes `tissue_positions.csv` with a header row, 1.x wrote
#' `tissue_positions_list.csv` without one. Both carry the same six columns in
#' the same order, so the header is detected from the first field and the
#' names are imposed either way.
#'
#' @param spatial_dir String. The `spatial/` directory.
#'
#' @return A `data.table` with `barcode`, `in_tissue`, `array_row`,
#' `array_col`, `pxl_row_in_fullres` and `pxl_col_in_fullres`.
#'
#' @keywords internal
.read_visium_positions <- function(spatial_dir) {
  checkmate::assertDirectoryExists(spatial_dir)

  pos_path <- locate(
    spatial_dir,
    "^tissue_positions(_list)?\\.csv(\\.gz)?$"
  )

  target_cols <- c(
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres"
  )

  con <- if (grepl("\\.gz$", pos_path, ignore.case = TRUE)) {
    gzfile(pos_path, open = "rt")
  } else {
    file(pos_path, open = "rt")
  }
  first_line <- readLines(con, n = 1L)
  close(con)

  has_hdr <- grepl("^\\s*\"?barcode\"?\\s*,", first_line, ignore.case = TRUE)

  positions <- data.table::fread(
    file = pos_path,
    header = has_hdr,
    showProgress = FALSE
  )

  if (!has_hdr) {
    if (ncol(positions) != length(target_cols)) {
      stop(sprintf(
        "Headerless positions file %s has %i columns, expected %i.",
        pos_path,
        ncol(positions),
        length(target_cols)
      ))
    }
    data.table::setnames(positions, target_cols)
  } else {
    data.table::setnames(positions, to_snake_case(names(positions)))
    missing_cols <- setdiff(target_cols, names(positions))
    if (length(missing_cols) > 0L) {
      stop(sprintf(
        "Positions file %s is missing column(s): %s",
        pos_path,
        paste(missing_cols, collapse = ", ")
      ))
    }
    positions <- positions[, target_cols, with = FALSE]
  }

  data.table::set(
    positions,
    j = "barcode",
    value = as.character(positions[["barcode"]])
  )

  return(positions[])
}

#' Parse `scalefactors_json.json` and derive the CytAssist scale factor
#'
#' @description
#' Everything in the JSON is passed through untouched. On CytAssist runs one
#' extra key, `cytassist_scalef`, is derived: `regist_target_img_scalef`
#' describes the registration target, not the `cytassist_image.tiff` Space
#' Ranger ships, which is smaller. The ratio is recovered from the pixel
#' dimensions of the files on disk rather than assumed, so a future Space
#' Ranger that ships the registration-sized frame lands on
#' `regist_target_img_scalef` unchanged.
#'
#' @param spatial_dir String. The `spatial/` directory.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return Named list of scale factors.
#'
#' @keywords internal
.visium_scale_factors <- function(spatial_dir, .verbose = TRUE) {
  sf_path <- locate(spatial_dir, "^scalefactors_json\\.json$")
  scale_factors <- as.list(jsonlite::fromJSON(sf_path))

  regist <- scale_factors[["regist_target_img_scalef"]]
  sources <- .visium_image_sources(spatial_dir)
  cyt_path <- sources[["cytassist"]]

  if (is.null(regist) || is.null(cyt_path)) {
    return(scale_factors)
  }

  fullres_dims <- .visium_fullres_frame(sources, scale_factors)
  cyt_dims <- .image_dimensions(cyt_path)

  if (is.null(fullres_dims) || is.null(cyt_dims)) {
    warning(paste(
      "Could not read the image dimensions needed to derive the CytAssist",
      "scale factor. Falling back to regist_target_img_scalef, which is",
      "correct only if Space Ranger shipped the registration-sized frame."
    ))
    scale_factors[["cytassist_scalef"]] <- regist
    return(scale_factors)
  }

  # the registration target is the full-res frame scaled by
  # regist_target_img_scalef; the shipped CytAssist file is that frame at some
  # other resolution. Both are square-cropped to the longest side.
  target_side <- max(fullres_dims) * regist
  ratio <- max(cyt_dims) / target_side

  if (!is.finite(ratio) || ratio <= 0 || ratio > 1.5) {
    warning(sprintf(
      paste(
        "Derived CytAssist/registration-target size ratio of %.4f is out of",
        "range. Falling back to regist_target_img_scalef."
      ),
      ratio
    ))
    scale_factors[["cytassist_scalef"]] <- regist
    return(scale_factors)
  }

  scale_factors[["cytassist_scalef"]] <- regist * ratio

  if (.verbose && abs(ratio - 1) > 1e-3) {
    message(sprintf(
      paste(
        " CytAssist image is %.0f x %.0f, the registration target is %.0f px",
        "on its longest side. Scale factor derived as %.6g, not the",
        "regist_target_img_scalef of %.6g."
      ),
      as.numeric(cyt_dims[["width"]]),
      as.numeric(cyt_dims[["height"]]),
      as.numeric(target_side),
      scale_factors[["cytassist_scalef"]],
      regist
    ))
  }

  return(scale_factors)
}

#' Recover the full-resolution frame size from a downsampled Visium image
#'
#' @description
#' Space Ranger never ships the full-resolution image, but it ships two known
#' downsamples of it. Dividing either one by its scale factor recovers the
#' frame the `pxl_*_in_fullres` coordinates live in.
#'
#' @param sources Named list of image paths from [.visium_image_sources()].
#' @param scale_factors Named list of parsed scale factors.
#'
#' @return Named numeric vector with `width` and `height`, or `NULL`.
#'
#' @keywords internal
.visium_fullres_frame <- function(sources, scale_factors) {
  candidates <- list(
    c(key = "hires", scalef = "tissue_hires_scalef"),
    c(key = "lowres", scalef = "tissue_lowres_scalef")
  )

  for (cand in candidates) {
    path <- sources[[cand[["key"]]]]
    scalef <- scale_factors[[cand[["scalef"]]]]
    if (is.null(path) || is.null(scalef) || !is.finite(scalef) || scalef <= 0) {
      next
    }
    dims <- .image_dimensions(path)
    if (is.null(dims)) {
      next
    }
    return(c(
      width = dims[["width"]] / scalef,
      height = dims[["height"]] / scalef
    ))
  }

  return(NULL)
}

#' Join a Visium positions table onto the barcode list
#'
#' @description
#' This is the load-bearing join of the whole reader. `barcodes.tsv.gz` is
#' lexicographically sorted, `tissue_positions.csv` is in array raster order.
#' They hold the same barcodes in different orders, so anything that zips the
#' two by row index scrambles every coordinate in the dataset and nothing
#' downstream notices. Match on the barcode string.
#'
#' @param barcodes Character vector of barcodes, in matrix column order.
#' @param positions `data.table` from [.read_visium_positions()].
#' @param exp_id String. The experiment identifier.
#'
#' @return A `data.table` with one row per barcode, in barcode order, carrying
#' the spatial columns and `exp_id`.
#'
#' @keywords internal
.visium_obs_dt <- function(barcodes, positions, exp_id) {
  # checks
  checkmate::qassert(barcodes, "S+")
  checkmate::assertDataTable(positions)
  checkmate::qassert(exp_id, "S1")

  idx <- match(barcodes, positions[["barcode"]])
  if (anyNA(idx)) {
    stop(sprintf(
      "%i barcode(s) of experiment '%s' have no row in the positions file.",
      sum(is.na(idx)),
      exp_id
    ))
  }

  data.table::data.table(
    barcode = barcodes,
    in_tissue = as.integer(positions[["in_tissue"]][idx]),
    array_row = as.integer(positions[["array_row"]][idx]),
    array_col = as.integer(positions[["array_col"]][idx]),
    pxl_row_in_fullres = as.integer(positions[["pxl_row_in_fullres"]][idx]),
    pxl_col_in_fullres = as.integer(positions[["pxl_col_in_fullres"]][idx]),
    exp_id = exp_id
  )
}

#' Read a 10x features file into a var table
#'
#' @param features_path String. Path to `features.tsv(.gz)` or `genes.tsv(.gz)`.
#'
#' @return A `data.table`. The first column is always the gene identifier;
#' three-column files get `gene_id`, `symbol`, `feature_type`.
#'
#' @keywords internal
.read_tenx_features <- function(features_path) {
  checkmate::assertFileExists(features_path)

  delim <- if (grepl("\\.tsv(\\.gz)?$", features_path, ignore.case = TRUE)) {
    "\t"
  } else {
    ","
  }
  features <- data.table::fread(
    file = features_path,
    sep = delim,
    header = FALSE,
    showProgress = FALSE
  )

  known <- c("gene_id", "symbol", "feature_type")
  if (ncol(features) <= length(known)) {
    data.table::setnames(features, known[seq_len(ncol(features))])
  } else {
    data.table::setnames(features, seq_along(known), known)
  }

  return(features[])
}

#' Copy the Visium images into the data directory
#'
#' @description
#' The layout contract puts images under `<dir_data>/images/<exp_id>/` and
#' stores paths relative to `dir_data` so the directory stays portable. With
#' `copy_images = FALSE` nothing is copied and no image is registered, because
#' a path outside `dir_data` cannot be expressed under that contract.
#'
#' @param dir_data String. The `SpatialSpot` data directory.
#' @param exp_id String. The experiment identifier.
#' @param sources Named list of source image paths.
#' @param slide_file Optional string. User-supplied full-resolution slide.
#' @param copy_images Boolean. Whether to copy at all.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return Named list of paths relative to `dir_data`.
#'
#' @keywords internal
.register_visium_images <- function(
  dir_data,
  exp_id,
  sources,
  slide_file,
  copy_images,
  .verbose = TRUE
) {
  if (!is.null(slide_file)) {
    sources[["fullres"]] <- slide_file
  }

  if (length(sources) == 0L) {
    return(list())
  }

  if (!copy_images) {
    if (.verbose) {
      message(paste(
        " copy_images = FALSE, so no images are registered.",
        "Image paths have to sit under the data directory."
      ))
    }
    return(list())
  }

  target_dir <- file.path(dir_data, "images", exp_id)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

  rel_paths <- purrr::imap(sources, \(src, key) {
    rel <- file.path(
      "images",
      exp_id,
      sprintf("%s.%s", key, tolower(tools::file_ext(src)))
    )
    ok <- file.copy(src, file.path(dir_data, rel), overwrite = TRUE)
    if (!ok) {
      stop(sprintf("Failed to copy image '%s' into %s.", src, target_dir))
    }
    rel
  })

  return(rel_paths)
}

#' Scan a single Visium directory
#'
#' @description
#' Resolves the file layout, decompresses the mtx where needed and reads
#' everything that is cheap to read. Shared by the single and multi loaders.
#'
#' @param visium_dir String. Space Ranger output directory.
#' @param exp_id String. The experiment identifier.
#' @param matrix_type String. One of `c("auto", "raw", "filtered")`.
#' @param temp_dir String. Directory decompressed files land in.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return A named list describing the input. `mtx_temp` is non-`NULL` when a
#' file was decompressed and has to be cleaned up by the caller.
#'
#' @keywords internal
.visium_scan_dir <- function(
  visium_dir,
  exp_id,
  matrix_type,
  temp_dir,
  .verbose = TRUE
) {
  # checks
  checkmate::assertDirectoryExists(visium_dir)
  checkmate::qassert(exp_id, "S1")
  checkmate::assertChoice(matrix_type, c("auto", "raw", "filtered"))
  checkmate::assertDirectoryExists(temp_dir)

  visium_dir <- path.expand(visium_dir)
  matrix_dir <- .visium_matrix_dir(visium_dir, matrix_type)
  spatial_dir <- file.path(visium_dir, "spatial")
  if (!dir.exists(spatial_dir)) {
    stop(sprintf("No 'spatial/' directory in %s.", visium_dir))
  }

  mtx_path <- locate(matrix_dir, "\\.mtx(\\.gz)?$")
  features_path <- locate(matrix_dir, "(features|genes)\\.(tsv|csv)(\\.gz)?$")
  barcodes_path <- locate(matrix_dir, "barcodes\\.(tsv|csv)(\\.gz)?$")

  mtx_temp <- NULL
  if (grepl("\\.gz$", mtx_path, ignore.case = TRUE)) {
    if (.verbose) {
      message(sprintf(" Decompressing the mtx file for '%s'.", exp_id))
    }
    mtx_temp <- gunzip_to_temp(mtx_path, temp_dir)
    mtx_path <- mtx_temp
  }

  barcodes <- as.character(
    data.table::fread(
      file = barcodes_path,
      header = FALSE,
      showProgress = FALSE
    )[[1L]]
  )
  features <- .read_tenx_features(features_path)

  list(
    exp_id = exp_id,
    visium_dir = visium_dir,
    spatial_dir = spatial_dir,
    mtx_path = mtx_path,
    mtx_temp = mtx_temp,
    features_path = features_path,
    barcodes_path = barcodes_path,
    cells_as_rows = FALSE,
    barcodes = barcodes,
    features = features,
    gene_ids = as.character(features[[1L]]),
    positions = .read_visium_positions(spatial_dir),
    scale_factors = .visium_scale_factors(spatial_dir, .verbose = .verbose),
    image_sources = .visium_image_sources(spatial_dir)
  )
}

#' Resolve the in-tissue spots to mtx cell indices
#'
#' @description
#' The allowlist the Rust mtx reader takes indexes the mtx file's own cell
#' axis, which for 10x output is the order of `barcodes.tsv.gz` and nothing
#' else. `tissue_positions.csv` is in array raster order, so this joins on the
#' barcode string and returns positions in the barcode vector, never positions
#' in the positions table.
#'
#' @param barcodes Character vector of barcodes, in matrix column order.
#' @param positions `data.table` from [.read_visium_positions()].
#' @param exp_id String. The experiment identifier, for the error message.
#'
#' @return Integer vector of 0-based mtx cell indices with `in_tissue == 1`.
#'
#' @keywords internal
.visium_in_tissue_allowlist <- function(barcodes, positions, exp_id) {
  # checks
  checkmate::qassert(barcodes, "S+")
  checkmate::assertDataTable(positions)
  checkmate::qassert(exp_id, "S1")

  idx <- match(barcodes, positions[["barcode"]])
  if (anyNA(idx)) {
    stop(sprintf(
      "%i barcode(s) of experiment '%s' have no row in the positions file.",
      sum(is.na(idx)),
      exp_id
    ))
  }

  in_tissue <- as.integer(positions[["in_tissue"]][idx])
  allowlist <- which(in_tissue == 1L) - 1L

  if (length(allowlist) == 0L) {
    stop(sprintf("No spot of experiment '%s' sits on tissue.", exp_id))
  }

  return(as.integer(allowlist))
}

#' Apply the in-tissue filter to a loaded spatial object
#'
#' @description
#' Sets the `to_keep` flag from the `in_tissue` column of the obs table. Off
#' tissue spots stay on disk and stay addressable, they simply drop out of
#' every filtered accessor.
#'
#' @param object A [bixverse::SpatialSpot()].
#' @param in_tissue_only Boolean. Whether to restrict to `in_tissue == 1`.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return The updated object.
#'
#' @keywords internal
.apply_in_tissue_filter <- function(object, in_tissue_only, .verbose = TRUE) {
  duckdb_con <- get_sc_duckdb(object)
  obs_tbl <- duckdb_con$get_obs_table(cols = c("cell_idx", "in_tissue"))

  keep <- if (in_tissue_only) {
    as.integer(obs_tbl[["cell_idx"]][obs_tbl[["in_tissue"]] == 1L])
  } else {
    as.integer(obs_tbl[["cell_idx"]])
  }

  if (length(keep) == 0L) {
    stop("No spots left after the in-tissue filter.")
  }

  # silent when nothing was dropped: `load_visium()` already excludes
  # off-tissue spots at ingest, so this is a no-op there and saying so twice
  # reads like the filter ran twice.
  if (.verbose && in_tissue_only && length(keep) < nrow(obs_tbl)) {
    message(sprintf(
      " Keeping %i of %i spots that sit on tissue.",
      length(keep),
      nrow(obs_tbl)
    ))
  }

  set_cells_to_keep(object, sort(keep))
}

## methods ---------------------------------------------------------------------

### single visium --------------------------------------------------------------

#' Load a Space Ranger Visium output into `SpatialSpot`
#'
#' @description
#' Reads one Space Ranger output directory into a [bixverse::SpatialSpot()].
#' The counts go through the streaming mtx reader into the Rust binary format,
#' the barcodes, tissue positions and `exp_id` go into the DuckDB obs table,
#' the images are copied into the data directory and the whole lot is
#' registered as a [bixverse::new_spatial_sample()].
#'
#' Spot coordinates are joined onto the barcodes by barcode string, never by
#' position: `barcodes.tsv.gz` is sorted lexicographically while
#' `tissue_positions.csv` is in array raster order.
#'
#' @param object `SpatialSpot` class.
#' @param visium_dir String. Space Ranger output directory. Must contain a
#' `spatial/` directory and at least one `*_feature_bc_matrix/` directory.
#' @param sp_io_param List. Output of [bixverse::params_sp_visium_io()]. A
#' list with the following elements:
#' \itemize{
#'   \item in_tissue_only - Boolean. Keep only spots with `in_tissue == 1`.
#'   \item matrix_type - String. Which count matrix directory to read.
#'   \item slide_file - Optional string. Full-resolution scanner slide.
#' }
#' @param sc_qc_param List. Output of [bixverse::params_sc_min_quality()]. A
#' list with the following elements:
#' \itemize{
#'   \item min_unique_genes - Integer. Minimum number of genes to be detected
#'   in the spot to be included.
#'   \item min_lib_size - Integer. Minimum library size in the spot to be
#'   included.
#'   \item min_cells - Integer. Minimum number of spots a gene needs to be
#'   detected in to be included.
#'   \item target_size - Float. Target size to normalise to.
#' }
#' @param exp_id String. The experiment identifier written into the obs table
#' and used as the key of the resulting `SpatialSample`.
#' @param copy_images Boolean. Copy the `spatial/` images into
#' `<dir_data>/images/<exp_id>/`. With `FALSE` no image is registered, since
#' the class stores image paths relative to the data directory. Defaults to
#' `TRUE`.
#' @param streaming Integer. CSR-to-CSC conversion mode. `0L` -> in-memory,
#' `1L` -> light streaming, `2L` -> heavy streaming with memory upper
#' boundaries. Defaults to `1L`.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @section In-tissue filtering:
#' `in_tissue_only` is applied at ingest, not afterwards. Off-tissue spots are
#' resolved to mtx cell indices by barcode and handed to the Rust reader as an
#' allowlist, so they never reach the gene statistics and never get written.
#' Gene-level QC in `sc_qc_param`, `min_cells` in particular, therefore sees
#' in-tissue spots only, and `dims[1]` is the in-tissue count rather than the
#' full capture area.
#'
#' With `in_tissue_only = FALSE` every spot is written, exactly as before.
#'
#' @section One sample per object:
#' This loader rewrites the whole store: the binary count files, the obs table
#' and the var table are all replaced. Registration is the only additive part,
#' so calling it twice with two different `exp_id`s would leave an object that
#' claims two samples while holding one sample's counts, with the first
#' sample's cached results still keyed against a gene axis that no longer
#' exists. It therefore errors when the object already carries a different
#' sample. Use [bixverse::prescan_visium_dirs()] plus
#' [bixverse::load_multi_visium()] for several sections.
#'
#' Re-running with the **same** `exp_id` is fine and is the intended way to
#' redo an ingest with different QC.
#'
#' @section Image resolutions:
#' `lowres` and `hires` map onto `tissue_lowres_scalef` and
#' `tissue_hires_scalef`. `cytassist` maps onto the derived `cytassist_scalef`,
#' **not** onto `regist_target_img_scalef`, which describes the larger
#' registration target rather than the file Space Ranger ships. `fullres` is
#' the user-supplied slide at scale `1.0`.
#'
#' @return The class with updated shape information, a populated DuckDB and one
#' registered `SpatialSample`.
#'
#' @export
load_visium <- S7::new_generic(
  name = "load_visium",
  dispatch_args = "object",
  fun = function(
    object,
    visium_dir,
    sp_io_param = params_sp_visium_io(),
    sc_qc_param = params_sc_min_quality(),
    exp_id,
    copy_images = TRUE,
    streaming = 1L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_visium SpatialSpot
#'
#' @export
S7::method(load_visium, SpatialSpot) <- function(
  object,
  visium_dir,
  sp_io_param = params_sp_visium_io(),
  sc_qc_param = params_sc_min_quality(),
  exp_id,
  copy_images = TRUE,
  streaming = 1L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SpatialSpot))
  checkmate::assertDirectoryExists(visium_dir)
  assertSpVisiumIo(sp_io_param)
  assertScMinQC(sc_qc_param)
  checkmate::qassert(exp_id, "S1")
  checkmate::qassert(copy_images, "B1")
  checkmate::qassert(streaming, "I1")
  checkmate::assertTRUE(streaming %in% c(0L, 1L, 2L))
  checkmate::qassert(.verbose, "B1")

  # The ingest below is destructive (counts, obs and var are all rewritten)
  # while `add_spatial_sample()` merges. A second `exp_id` would therefore
  # register a sample whose data has just been overwritten by the incoming one.
  # Guard before anything is written; re-running the same `exp_id` is fine.
  registered <- get_sample_ids(object)
  if (length(registered) > 0L && !identical(registered, exp_id)) {
    stop(sprintf(
      paste(
        "The object already holds sample(s) %s. `load_visium()` rewrites the",
        "whole store, so it cannot add '%s' on top. Use",
        "`prescan_visium_dirs()` and `load_multi_visium()` for several",
        "sections."
      ),
      toString(sprintf("'%s'", registered)),
      exp_id
    ))
  }

  temp_dir <- tempfile(pattern = "bixverse_visium_")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  if (.verbose) {
    message(sprintf("Scanning the Visium directory for '%s'.", exp_id))
  }

  scan <- .visium_scan_dir(
    visium_dir = visium_dir,
    exp_id = exp_id,
    matrix_type = sp_io_param$matrix_type,
    temp_dir = temp_dir,
    .verbose = .verbose
  )

  rust_con <- get_sc_rust_ptr(object)

  # off-tissue spots are excluded before the counts are written, so gene QC
  # never sees them. Anything after the fact would still compute `min_cells`
  # against the full capture area.
  file_res <- if (sp_io_param$in_tissue_only) {
    allowlist <- .visium_in_tissue_allowlist(
      barcodes = scan$barcodes,
      positions = scan$positions,
      exp_id = exp_id
    )
    if (.verbose) {
      message(sprintf(
        "Ingesting %i of %i spots that sit on tissue.",
        length(allowlist),
        length(scan$barcodes)
      ))
    }
    rust_con$mtx_to_file_streaming_subset(
      mtx_path = scan$mtx_path,
      qc_params = sc_qc_param,
      cells_as_rows = scan$cells_as_rows,
      cell_allowlist = allowlist,
      verbose = .verbose
    )
  } else {
    rust_con$mtx_to_file_streaming(
      mtx_path = scan$mtx_path,
      qc_params = sc_qc_param,
      cells_as_rows = scan$cells_as_rows,
      verbose = .verbose
    )
  }

  .dispatch_gene_based_data(
    rust_con = rust_con,
    streaming = streaming,
    batch_size = 1000L,
    max_genes_in_memory = 2000L,
    cell_batch_size = 100000L,
    .verbose = .verbose
  )

  gene_nnz <- rust_con$get_nnz_genes(gene_indices = NULL)
  gene_nnz_dt <- data.table::data.table(no_cells_exp = gene_nnz)

  duckdb_con <- get_sc_duckdb(object)

  if (.verbose) {
    message("Joining the tissue positions onto the barcodes.")
  }

  obs_dt <- .visium_obs_dt(
    barcodes = scan$barcodes,
    positions = scan$positions,
    exp_id = exp_id
  )

  duckdb_con$populate_obs_from_data.table(
    obs_dt = obs_dt,
    filter = as.integer(file_res$cell_indices + 1L)
  )
  duckdb_con$populate_var_from_data.table(
    var_dt = scan$features,
    filter = as.integer(file_res$gene_indices + 1L)
  )

  cell_res_dt <- data.table::setDT(file_res[c("nnz", "lib_size")])

  duckdb_con$add_data_obs(new_data = cell_res_dt)
  duckdb_con$add_data_var(new_data = gene_nnz_dt)
  duckdb_con$set_to_keep_column()

  S7::prop(object, "dims") <- as.integer(rust_con$get_shape())
  object <- set_cell_mapping(
    x = object,
    cell_map = duckdb_con$get_obs_index_map()
  )
  object <- set_gene_mapping(
    x = object,
    gene_map = duckdb_con$get_var_index_map()
  )

  object <- .apply_in_tissue_filter(
    object = object,
    in_tissue_only = sp_io_param$in_tissue_only,
    .verbose = .verbose
  )

  if (.verbose) {
    message("Registering the spatial sample.")
  }

  image_paths <- .register_visium_images(
    dir_data = S7::prop(object, "dir_data"),
    exp_id = exp_id,
    sources = scan$image_sources,
    slide_file = sp_io_param$slide_file,
    copy_images = copy_images,
    .verbose = .verbose
  )

  # per sample, not object-wide, matching `load_multi_visium()`
  kept_exp_ids <- duckdb_con$get_obs_table(
    cols = "exp_id",
    filtered = TRUE
  )[["exp_id"]]

  sample <- new_spatial_sample(
    exp_id = exp_id,
    technology = "visium",
    n_spots = sum(kept_exp_ids == exp_id),
    scale_factors = scan$scale_factors,
    image_paths = image_paths
  )

  object <- add_spatial_sample(object, sample)

  return(object)
}

### multiple visium ------------------------------------------------------------

#' Prescan multiple Visium directories for a multi-load
#'
#' @description
#' Walks each Space Ranger output directory, builds the **intersection** of
#' gene identifiers across inputs (matched on the first column of the features
#' file), decompresses the mtx files into a temporary directory and assembles
#' the per-input tasks [bixverse::load_multi_visium()] expects.
#'
#' The count matrix directory is resolved per input, preferring
#' `filtered_feature_bc_matrix/` over `raw_feature_bc_matrix/`.
#'
#' @param dirs Character vector of Space Ranger output directories. Length >= 2.
#' @param exp_ids Character vector of experiment identifiers, one per
#' directory. Must be unique.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return A list with:
#' \itemize{
#'   \item universe - Character vector of gene IDs in the intersection.
#'   \item universe_size - Length of the universe.
#'   \item file_tasks - List of per-input task lists.
#'   \item temp_files - Character vector of temp files created during
#'   decompression; the caller should `unlink()` these after use.
#' }
#'
#' @export
prescan_visium_dirs <- function(dirs, exp_ids, .verbose = TRUE) {
  # checks
  checkmate::assertCharacter(dirs, min.len = 2L)
  for (d in dirs) {
    checkmate::assertDirectoryExists(d)
  }
  checkmate::assertCharacter(exp_ids, len = length(dirs), unique = TRUE)
  checkmate::qassert(.verbose, "B1")

  temp_dir <- tempfile(pattern = "bixverse_visium_prescan_")
  dir.create(temp_dir)

  if (.verbose) {
    cli::cli_progress_bar(
      "Prescanning Visium directories",
      total = length(dirs)
    )
  }

  file_tasks <- vector("list", length(dirs))
  for (i in seq_along(dirs)) {
    if (.verbose) {
      cli::cli_progress_update(status = exp_ids[i])
    }
    file_tasks[[i]] <- .visium_scan_dir(
      visium_dir = dirs[i],
      exp_id = exp_ids[i],
      matrix_type = "auto",
      temp_dir = temp_dir,
      .verbose = FALSE
    )
  }

  if (.verbose) {
    cli::cli_progress_done()
  }

  universe <- Reduce(intersect, lapply(file_tasks, \(t) t$gene_ids))
  if (length(universe) == 0L) {
    stop("Gene intersection across inputs is empty.")
  }

  for (i in seq_along(file_tasks)) {
    match_idx <- match(file_tasks[[i]]$gene_ids, universe)
    file_tasks[[i]]$gene_local_to_universe <- as.integer(match_idx - 1L)
  }

  names(file_tasks) <- exp_ids

  list(
    universe = universe,
    universe_size = length(universe),
    file_tasks = file_tasks,
    temp_files = unlist(lapply(file_tasks, \(t) t$mtx_temp), use.names = FALSE)
  )
}

#' Load multiple Visium directories into a single `SpatialSpot`
#'
#' @description
#' Takes the result of [bixverse::prescan_visium_dirs()] and loads all inputs
#' into one experiment with global gene QC and sequential spot indexing. The
#' feature space is the **intersection** of input gene IDs. Each input is
#' registered as its own `SpatialSample`.
#'
#' @param object `SpatialSpot` class.
#' @param prescan_result Output of [bixverse::prescan_visium_dirs()].
#' @param sp_io_param List. Output of [bixverse::params_sp_visium_io()].
#' `matrix_type` was already resolved during the prescan and is ignored here.
#' `slide_file` has no meaning across several inputs and must be `NULL`.
#' @param sc_qc_param List. Output of [bixverse::params_sc_min_quality()].
#' @param copy_images Boolean. Copy the images of every input into
#' `<dir_data>/images/<exp_id>/`. Defaults to `TRUE`.
#' @param streaming Integer. CSR-to-CSC conversion mode. `0L` -> in-memory,
#' `1L` -> light streaming, `2L` -> heavy streaming. Defaults to `1L`.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @section In-tissue filtering:
#' Unlike [bixverse::load_visium()], this path still applies `in_tissue_only`
#' through the `to_keep` flag after the counts are written. The multi-file
#' reader has no cell allowlist, so global gene QC across inputs sees the
#' off-tissue counts.
#'
#' @return The class with updated shape, a populated DuckDB and one
#' `SpatialSample` per input.
#'
#' @export
load_multi_visium <- S7::new_generic(
  name = "load_multi_visium",
  dispatch_args = "object",
  fun = function(
    object,
    prescan_result,
    sp_io_param = params_sp_visium_io(),
    sc_qc_param = params_sc_min_quality(),
    copy_images = TRUE,
    streaming = 1L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_multi_visium SpatialSpot
#'
#' @export
S7::method(load_multi_visium, SpatialSpot) <- function(
  object,
  prescan_result,
  sp_io_param = params_sp_visium_io(),
  sc_qc_param = params_sc_min_quality(),
  copy_images = TRUE,
  streaming = 1L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SpatialSpot))
  assertSpVisiumIo(sp_io_param)
  assertScMinQC(sc_qc_param)
  checkmate::assertList(prescan_result)
  checkmate::assertTRUE(all(
    c("universe", "universe_size", "file_tasks", "temp_files") %in%
      names(prescan_result)
  ))
  checkmate::assertNull(
    sp_io_param$slide_file,
    .var.name = "sp_io_param$slide_file"
  )
  checkmate::qassert(copy_images, "B1")
  checkmate::qassert(streaming, "I1")
  checkmate::assertTRUE(streaming %in% c(0L, 1L, 2L))
  checkmate::qassert(.verbose, "B1")

  on.exit(
    if (length(prescan_result$temp_files) > 0L) {
      unlink(prescan_result$temp_files)
    },
    add = TRUE
  )

  rust_con <- get_sc_rust_ptr(object)

  rust_tasks <- lapply(prescan_result$file_tasks, \(t) {
    list(
      exp_id = t$exp_id,
      mtx_path = t$mtx_path,
      cells_as_rows = t$cells_as_rows,
      gene_local_to_universe = t$gene_local_to_universe
    )
  })

  file_res <- rust_con$multi_mtx_to_file(
    file_tasks = rust_tasks,
    universe_size = as.integer(prescan_result$universe_size),
    qc_params = sc_qc_param,
    verbose = .verbose
  )

  .dispatch_gene_based_data(
    rust_con = rust_con,
    streaming = streaming,
    batch_size = 1000L,
    max_genes_in_memory = 2000L,
    cell_batch_size = 100000L,
    .verbose = .verbose
  )

  gene_nnz <- rust_con$get_nnz_genes(gene_indices = NULL)
  gene_nnz_dt <- data.table::data.table(no_cells_exp = gene_nnz)

  duckdb_con <- get_sc_duckdb(object)

  if (.verbose) {
    message("Joining the tissue positions onto the barcodes.")
  }

  obs_dt <- data.table::rbindlist(
    lapply(seq_along(file_res$per_file), \(i) {
      res <- file_res$per_file[[i]]
      task <- prescan_result$file_tasks[[res$exp_id]]
      full_obs <- .visium_obs_dt(
        barcodes = task$barcodes,
        positions = task$positions,
        exp_id = res$exp_id
      )
      full_obs[as.integer(res$cell_indices + 1L)]
    })
  )

  duckdb_con$populate_obs_from_data.table(obs_dt = obs_dt)

  final_gene_names <- prescan_result$universe[file_res$global_gene_indices + 1L]
  duckdb_con$populate_var_minimal(final_gene_names = final_gene_names)

  cell_res_dt <- data.table::rbindlist(lapply(file_res$per_file, \(f) {
    data.table::data.table(nnz = f$nnz, lib_size = f$lib_size)
  }))

  duckdb_con$add_data_obs(new_data = cell_res_dt)
  duckdb_con$add_data_var(new_data = gene_nnz_dt)
  duckdb_con$set_to_keep_column()

  S7::prop(object, "dims") <- as.integer(rust_con$get_shape())
  object <- set_cell_mapping(
    x = object,
    cell_map = duckdb_con$get_obs_index_map()
  )
  object <- set_gene_mapping(
    x = object,
    gene_map = duckdb_con$get_var_index_map()
  )

  object <- .apply_in_tissue_filter(
    object = object,
    in_tissue_only = sp_io_param$in_tissue_only,
    .verbose = .verbose
  )

  if (.verbose) {
    message("Registering the spatial samples.")
  }

  kept_exp_ids <- duckdb_con$get_obs_table(
    cols = "exp_id",
    filtered = TRUE
  )[["exp_id"]]

  for (res in file_res$per_file) {
    task <- prescan_result$file_tasks[[res$exp_id]]
    image_paths <- .register_visium_images(
      dir_data = S7::prop(object, "dir_data"),
      exp_id = res$exp_id,
      sources = task$image_sources,
      slide_file = NULL,
      copy_images = copy_images,
      .verbose = .verbose
    )
    object <- add_spatial_sample(
      object,
      new_spatial_sample(
        exp_id = res$exp_id,
        technology = "visium",
        n_spots = sum(kept_exp_ids == res$exp_id),
        scale_factors = task$scale_factors,
        image_paths = image_paths
      )
    )
  }

  return(object)
}
