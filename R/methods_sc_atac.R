# single cell atac methods -----------------------------------------------------

## helpers ---------------------------------------------------------------------

#' Map cell ids to the barcodes of a fragment file
#'
#' @param ids Character vector. Cell ids in store order.
#' @param barcode_map `NULL`, a named character vector (names are cell ids,
#' values the fragment barcodes) or a function taking and returning a
#' character vector.
#'
#' @returns Character vector of fragment barcodes, one per id.
#'
#' @keywords internal
.fragment_barcodes <- function(ids, barcode_map) {
  ids <- as.character(ids)
  if (is.null(barcode_map)) {
    return(ids)
  }
  if (is.function(barcode_map)) {
    out <- as.character(barcode_map(ids))
    checkmate::assertCharacter(out, len = length(ids), any.missing = FALSE)
    return(out)
  }
  checkmate::assertCharacter(barcode_map, names = "named")
  missing <- setdiff(ids, names(barcode_map))
  if (length(missing)) {
    stop(sprintf(
      "barcode_map lacks %i cell ids, e.g. '%s'.",
      length(missing),
      missing[1]
    ))
  }

  unname(barcode_map[ids])
}

#' Check a fragment file and note what it lacks for a fast read
#'
#' @param fragments String. The fragment file.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @returns Invisible `NULL`. Errors when the file is missing, warns when the
#' tabix index is older than the file.
#'
#' @keywords internal
.check_fragments <- function(fragments, .verbose = TRUE) {
  checkmate::assertFileExists(fragments)

  con <- file(fragments, "rb")
  head <- readBin(con, "raw", n = 16L)
  close(con)
  tbi <- paste0(fragments, ".tbi")
  gzip <- length(head) >= 2 &&
    head[1] == as.raw(0x1f) &&
    head[2] == as.raw(0x8b)
  bgzf <- gzip && length(head) >= 14 && rawToChar(head[13:14]) == "BC"

  if (bgzf && file.exists(tbi)) {
    if (file.mtime(tbi) < file.mtime(fragments)) {
      warning(sprintf(
        paste(
          "%s is older than the fragment file. If the file was rewritten,",
          "regenerate the index (tabix -p bed); a stale one is refused or",
          "reads the wrong chromosomes."
        ),
        tbi
      ))
    }
    return(invisible(NULL))
  }
  if (.verbose) {
    what <- if (bgzf) {
      "no tabix index next to it"
    } else if (gzip) {
      "plain gzip, which inflates on one thread"
    } else {
      "uncompressed"
    }
    message(sprintf(
      "%s: %s, so passes read the whole file at once.",
      basename(fragments),
      what
    ))
  }

  invisible(NULL)
}

#' Warn about cells without fragments and fragments off the targets
#'
#' @param res List. A Rust result with `cells_missing`,
#' `off_target_fraction` and `off_target_chroms`.
#' @param n_cells Integer. Number of requested cells.
#' @param what String. What the targets are, for the message.
#'
#' @returns Invisible `NULL`.
#'
#' @keywords internal
.report_coverage <- function(res, n_cells, what) {
  missing <- as.integer(res$cells_missing)
  if (n_cells > 0 && missing / n_cells > .ATAC_MISSING_CELLS_WARN) {
    warning(sprintf(
      paste(
        "%i of %i cells have no fragment in the file. If the cell ids and",
        "the fragment barcodes differ (a '-1' suffix, a sample prefix), pass",
        "barcode_map."
      ),
      missing,
      n_cells
    ))
  }
  off <- res$off_target_fraction
  if (off > .ATAC_OFF_TARGET_WARN) {
    warning(sprintf(
      paste(
        "%.0f%% of the fragments lie on chromosomes without %s (%s). Check",
        "the chromosome naming (chr1 against 1) and the genome build."
      ),
      100 * off,
      what,
      paste(res$off_target_chroms, collapse = ", ")
    ))
  }

  invisible(NULL)
}

#' `chrom`, `start` and `end` of a BED file or a table of regions
#'
#' @param regions String (a BED file, plain or gzip) or a data.frame with
#' `chrom`, `start` and `end` (0-based, half-open).
#'
#' @returns A list with `chrom`, `start` and `end`.
#'
#' @keywords internal
.regions_from <- function(regions) {
  if (checkmate::testString(regions)) {
    checkmate::assertFileExists(regions)
    return(rs_atac_read_bed(path.expand(regions)))
  }
  checkmate::assertDataFrame(regions)
  missing <- setdiff(c("chrom", "start", "end"), names(regions))
  if (length(missing)) {
    stop(sprintf(
      "The regions lack the column(s) %s.",
      paste(missing, collapse = ", ")
    ))
  }

  list(
    chrom = as.character(regions$chrom),
    start = as.integer(regions$start),
    end = as.integer(regions$end)
  )
}

#' The features to count fragments into
#'
#' @param peaks `NULL`, a BED file or a data.frame of regions.
#' @param tile_width Optional integer. Tile width in bp, instead of `peaks`.
#' @param chrom_sizes Optional named integer vector of chromosome lengths,
#' for the tiles.
#'
#' @returns A list with `chrom`, `start`, `end` and `kind` (`"peak"` or
#' `"tile"`).
#'
#' @keywords internal
.atac_features <- function(peaks, tile_width, chrom_sizes) {
  checkmate::qassert(tile_width, c("0", "X1[1,)"))
  exactly_one <- "Pass exactly one of peaks and tile_width."
  if (!is.null(peaks) && !is.null(tile_width)) {
    stop(exactly_one)
  }
  if (is.null(peaks)) {
    if (is.null(tile_width)) {
      stop(exactly_one)
    }
    if (length(chrom_sizes) == 0) {
      stop("Tiles need chrom_sizes.")
    }
    checkmate::assertIntegerish(chrom_sizes, lower = 1, names = "named")
    tiles <- rs_atac_genome_tiles(
      chroms = names(chrom_sizes),
      sizes = as.integer(chrom_sizes),
      width = as.integer(tile_width)
    )
    return(c(tiles, kind = "tile"))
  }

  c(.regions_from(peaks), kind = "peak")
}

#' `chr:start-end` ids of regions, 0-based half-open
#'
#' @param regions List with `chrom`, `start` and `end`.
#'
#' @returns Character vector.
#'
#' @keywords internal
.region_ids <- function(regions) {
  sprintf("%s:%i-%i", regions$chrom, regions$start, regions$end)
}

#' Count a fragment file into a store and report on the pass
#'
#' @param fragments String. The fragment file.
#' @param store_dir String. Directory of the store to write.
#' @param barcodes Character vector. Cells in store order.
#' @param regions List with `chrom`, `start` and `end` of the features.
#' @param mode String. The count mode.
#' @param kind String. One of `c("peak", "tile", "gene_activity")`.
#' @param target_size Numeric. Target size of the initial log-CPM norm.
#' @param max_mem_gb Optional numeric. Memory budget of the feature-major
#' conversion.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @returns The Rust result: `lib_size`, `nnz` and the coverage keys.
#'
#' @keywords internal
.count_store <- function(
  fragments,
  store_dir,
  barcodes,
  regions,
  mode,
  kind,
  target_size,
  max_mem_gb,
  .verbose
) {
  .check_fragments(fragments, .verbose)
  dir.create(store_dir, showWarnings = FALSE, recursive = TRUE)
  res <- rs_atac_fragments_to_store(
    fragments = fragments,
    f_path_cells = file.path(store_dir, "counts_cells.bin"),
    f_path_genes = file.path(store_dir, "counts_genes.bin"),
    barcodes = barcodes,
    chrom = regions$chrom,
    start = as.integer(regions$start),
    end = as.integer(regions$end),
    mode = mode,
    kind = kind,
    target_size = target_size,
    max_mem_gb = max_mem_gb,
    verbose = as.integer(.verbose)
  )
  .report_coverage(res, length(barcodes), "features")

  res
}

#' Fraction of fragments in peaks from the ATAC library sizes
#'
#' @description
#' In the store's count unit: with insertion counts the share of insertions
#' in peaks (two per fragment), otherwise Signac's fragments in peaks over
#' fragments, which can exceed 1 when a fragment overlaps two peaks.
#'
#' @param x `SingleCellsAtac` or `SingleCellsMultiModal`.
#' @param n_fragments Numeric vector. Fragments per cell, every cell.
#'
#' @returns Numeric vector, or `NULL` for tiles, which cover everything.
#'
#' @keywords internal
.atac_frip <- function(x, n_fragments) {
  info <- S7::prop(x, "atac_info")
  if (is.null(info) || info$kind != "peak") {
    return(NULL)
  }
  lib_col <- .atac_lib_col(x)
  lib <- get_sc_obs(x, cols = lib_col)[[lib_col]]
  per_fragment <- if (info$count_mode == "insertion") 2 else 1

  as.numeric(lib) / (per_fragment * as.numeric(n_fragments))
}

#' Recompute `atac_frip` when a fragment QC left `atac_n_fragments`
#'
#' @param x `SingleCellsAtac` or `SingleCellsMultiModal`.
#'
#' @returns Invisible `NULL`.
#'
#' @keywords internal
.refresh_frip <- function(x) {
  duckdb_con <- get_sc_duckdb(x)
  obs_cols <- duckdb_con$get_obs_cols()
  if (!"atac_n_fragments" %in% obs_cols) {
    return(invisible(NULL))
  }
  frip <- .atac_frip(x, get_sc_obs(x, cols = "atac_n_fragments")[[1]])
  if (is.null(frip)) {
    if ("atac_frip" %in% obs_cols) {
      duckdb_con$drop_columns(table = "obs", cols = "atac_frip")
    }
    return(invisible(NULL))
  }
  duckdb_con$add_data_obs(data.table::data.table(atac_frip = frip))

  invisible(NULL)
}

#' 0-based codes of an obs column over the kept cells
#'
#' @param x The single cell object.
#' @param col String. The obs column.
#'
#' @returns A list with `codes` (0-based integer, one per kept cell, in kept
#' order) and `levels`.
#'
#' @keywords internal
.kept_obs_codes <- function(x, col) {
  checkmate::qassert(col, "S1")
  obs <- get_sc_obs(x, cols = col, filtered = TRUE)
  if (!col %in% names(obs)) {
    stop(sprintf("Column '%s' not found in obs.", col))
  }
  values <- obs[[col]]
  if (anyNA(values)) {
    stop(sprintf("Column '%s' has missing values in the kept cells.", col))
  }
  f <- factor(values)

  list(codes = as.integer(f) - 1L, levels = levels(f))
}

#' TSSs of the GTF genes of the biotypes
#'
#' @param gtf String. The GTF.
#' @param biotypes Optional character vector. `NULL` keeps all.
#'
#' @returns The Rust result: `chrom`, `pos`, `strand`, `gene_id`,
#' `gene_name`.
#'
#' @keywords internal
.atac_tss <- function(gtf, biotypes) {
  checkmate::assertFileExists(gtf)
  checkmate::qassert(biotypes, c("0", "S+"))
  rs_atac_tss_positions(path.expand(gtf), biotypes %||% character())
}

#' Write the ATAC tables of a multi-modal object after a new store
#'
#' @description
#' Writes `var_atac`, the per-cell totals to obs and the count mode, resets
#' the ATAC state and cache, and recomputes a FRiP from an earlier fragment
#' QC against the new counts.
#'
#' @param object `SingleCellsMultiModal` class.
#' @param feature_ids Character vector. One id per feature.
#' @param regions List with `chrom`, `start` and `end`.
#' @param lib_size,nnz Numeric vectors. Per-cell totals, every cell.
#' @param count_mode String. The count mode.
#'
#' @returns The object with the ATAC store attached.
#'
#' @keywords internal
.write_atac_tables <- function(
  object,
  feature_ids,
  regions,
  lib_size,
  nnz,
  count_mode
) {
  duckdb_con <- get_sc_duckdb(object)
  duckdb_con$populate_var_atac_from_data.table(data.table::data.table(
    feature_id = feature_ids,
    chrom = regions$chrom,
    start = as.integer(regions$start),
    end = as.integer(regions$end)
  ))
  duckdb_con$add_data_obs(data.table::data.table(
    atac_lib_size = as.numeric(lib_size),
    atac_nnz = as.numeric(nnz)
  ))

  store_dir <- .atac_store_dir(object)
  .write_count_mode(store_dir, count_mode)
  paths <- .atac_paths(object)
  S7::prop(object, "peak_connection") <- SingleCellCountData$new(
    f_path_cells = paths$cells,
    f_path_genes = paths$genes
  )
  S7::prop(object, "atac_info") <- .read_atac_info(store_dir)
  S7::prop(object, "atac_cache") <- new_sc_cache()
  .refresh_frip(object)

  object
}

#' Refuse to replace an ATAC store unless asked
#'
#' @param object `SingleCellsMultiModal` class.
#' @param overwrite Boolean.
#'
#' @returns Invisible `NULL`.
#'
#' @keywords internal
.check_atac_overwrite <- function(object, overwrite) {
  if (!is.null(S7::prop(object, "atac_info")) && !overwrite) {
    stop(paste(
      "This experiment has an ATAC store. Pass overwrite = TRUE to replace",
      "it."
    ))
  }

  invisible(NULL)
}

## ingest ----------------------------------------------------------------------

### 10x ------------------------------------------------------------------------

#' Add the peaks of a 10x multiome h5 as the ATAC store
#'
#' @description
#' Reads the `Peaks` rows of a Cell Ranger ARC
#' `filtered_feature_bc_matrix.h5` into a peak store under
#' `<dir_data>/atac`, for every cell of the RNA store and in its order, so
#' the RNA cell filter applies unchanged. The features go to the DuckDB table
#' `var_atac`, the per-cell totals to obs as `atac_lib_size` and `atac_nnz`.
#' The 10x peak matrix counts insertions.
#'
#' @param object `SingleCellsMultiModal` class, its RNA loaded from the same
#' h5 (or with the same barcodes).
#' @param h5_path String. The Cell Ranger ARC h5.
#' @param overwrite Boolean. Replace an existing ATAC store.
#' @param max_mem_gb Optional numeric. Memory budget of the peak-major
#' conversion.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @returns The object with the ATAC store attached.
#'
#' @export
add_atac_tenx_sc <- S7::new_generic(
  name = "add_atac_tenx_sc",
  dispatch_args = "object",
  fun = function(
    object,
    h5_path,
    overwrite = FALSE,
    max_mem_gb = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method add_atac_tenx_sc SingleCellsMultiModal
#'
#' @export
S7::method(add_atac_tenx_sc, SingleCellsMultiModal) <- function(
  object,
  h5_path,
  overwrite = FALSE,
  max_mem_gb = NULL,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsMultiModal))
  checkmate::assertFileExists(h5_path)
  checkmate::qassert(overwrite, "B1")
  checkmate::qassert(max_mem_gb, c("0", "N1(0,)"))
  checkmate::qassert(.verbose, "B1")
  .check_atac_overwrite(object, overwrite)

  # function body
  store_dir <- .atac_store_dir(object)
  dir.create(store_dir, showWarnings = FALSE)
  paths <- .atac_paths(object)
  res <- rs_atac_tenx_peaks_to_store(
    h5_path = path.expand(h5_path),
    f_path_cells = paths$cells,
    f_path_genes = paths$genes,
    barcodes = get_cell_names(object, filtered = FALSE),
    max_mem_gb = max_mem_gb,
    verbose = as.integer(.verbose)
  )

  .write_atac_tables(
    object = object,
    feature_ids = res$name,
    regions = res[c("chrom", "start", "end")],
    lib_size = res$lib_size,
    nnz = res$nnz,
    count_mode = "insertion"
  )
}

### fragments ------------------------------------------------------------------

#' Count a fragment file into peaks or genome tiles as the ATAC store
#'
#' @description
#' Counts every cell of the RNA store, in its order, so the RNA cell filter
#' applies unchanged. Pass exactly one of `peaks` and `tile_width`. The
#' features go to `var_atac` with `chr:start-end` ids (0-based, half-open),
#' the per-cell totals to obs as `atac_lib_size` and `atac_nnz`.
#'
#' @param object `SingleCellsMultiModal` class.
#' @param fragments String. The fragment file (BGZF, gzip or plain). With a
#' tabix index next to it, passes read one chromosome at a time.
#' @param peaks Optional. A BED file or a data.frame with `chrom`, `start`,
#' `end` (0-based, half-open), e.g. from [bixverse::call_peaks_sc()].
#' @param tile_width Optional integer. Tile width in bp, instead of `peaks`.
#' @param chrom_sizes Optional named integer vector of chromosome lengths,
#' for the tiles.
#' @param mode String. What one count is. `"insertion"` (each Tn5 insertion,
#' Cell Ranger and ArchR), `"paired-insertion"` (once per fragment and
#' feature, SnapATAC2) or `"fragment"` (each overlapping fragment, Signac).
#' @param barcode_map Optional. Cell id to fragment barcode, where they
#' differ: a named character vector (names are the cell ids) or a function.
#' @param overwrite Boolean. Replace an existing ATAC store.
#' @param max_mem_gb Optional numeric. Memory budget of the peak-major
#' conversion.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @returns The object with the ATAC store attached.
#'
#' @export
add_atac_fragments_sc <- S7::new_generic(
  name = "add_atac_fragments_sc",
  dispatch_args = "object",
  fun = function(
    object,
    fragments,
    peaks = NULL,
    tile_width = NULL,
    chrom_sizes = NULL,
    mode = c("insertion", "paired-insertion", "fragment"),
    barcode_map = NULL,
    overwrite = FALSE,
    max_mem_gb = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method add_atac_fragments_sc SingleCellsMultiModal
#'
#' @export
S7::method(add_atac_fragments_sc, SingleCellsMultiModal) <- function(
  object,
  fragments,
  peaks = NULL,
  tile_width = NULL,
  chrom_sizes = NULL,
  mode = c("insertion", "paired-insertion", "fragment"),
  barcode_map = NULL,
  overwrite = FALSE,
  max_mem_gb = NULL,
  .verbose = TRUE
) {
  mode <- match.arg(mode)
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsMultiModal))
  checkmate::qassert(fragments, "S1")
  checkmate::assertChoice(mode, .ATAC_COUNT_MODES)
  checkmate::qassert(overwrite, "B1")
  checkmate::qassert(max_mem_gb, c("0", "N1(0,)"))
  checkmate::qassert(.verbose, "B1")

  # function body
  features <- .atac_features(peaks, tile_width, chrom_sizes)
  .check_atac_overwrite(object, overwrite)
  res <- .count_store(
    fragments = path.expand(fragments),
    store_dir = .atac_store_dir(object),
    barcodes = .fragment_barcodes(
      get_cell_names(object, filtered = FALSE),
      barcode_map
    ),
    regions = features,
    mode = mode,
    kind = features$kind,
    target_size = .ATAC_TF_IDF_SCALE,
    max_mem_gb = max_mem_gb,
    .verbose = .verbose
  )

  .write_atac_tables(
    object = object,
    feature_ids = .region_ids(features),
    regions = features,
    lib_size = res$lib_size,
    nnz = res$nnz,
    count_mode = mode
  )
}

### cell calling ---------------------------------------------------------------

#' Read the cell barcodes Cell Ranger ATAC called
#'
#' @param path String. `singlecell.csv` (rows with `is__cell_barcode == 1`)
#' or a filtered `*_bc_matrix.h5` (its `matrix/barcodes`).
#'
#' @returns Character vector of barcodes, in file order.
#'
#' @export
read_cellranger_barcodes <- function(path) {
  checkmate::assertFileExists(path)
  path <- path.expand(path)

  ext <- tools::file_ext(path)
  if (ext == "h5") {
    return(as.character(rs_h5_read(path, "matrix/barcodes")))
  }
  if (ext != "csv") {
    stop(sprintf(
      "Expected singlecell.csv or a filtered h5, got %s.",
      basename(path)
    ))
  }
  frame <- data.table::fread(path, colClasses = "character")
  if (!"is__cell_barcode" %in% names(frame)) {
    stop(sprintf("%s has no is__cell_barcode column.", basename(path)))
  }
  calls <- frame$barcode[frame$is__cell_barcode == "1"]

  calls[calls != "NO_BARCODE"]
}

#' Fragments per barcode, every barcode in the file
#'
#' @description
#' For cell calling: plot the counts against their rank (a knee plot) and
#' pick the `min_fragments` of [bixverse::load_atac_fragments()].
#'
#' @param fragments String. The fragment file.
#' @param exclude_chroms Character vector. Chromosomes not counted.
#'
#' @returns A data.table with `barcode` and `n_fragments`, most fragments
#' first.
#'
#' @export
barcode_fragment_counts <- function(
  fragments,
  exclude_chroms = c("chrM", "MT", "chrMT", "M")
) {
  checkmate::assertFileExists(fragments)
  checkmate::qassert(exclude_chroms, "S*")

  res <- rs_atac_barcode_fragment_counts(
    path.expand(fragments),
    exclude_chroms
  )
  out <- data.table::data.table(
    barcode = res$barcode,
    n_fragments = res$n_fragments
  )
  data.table::setorderv(out, c("n_fragments", "barcode"), c(-1L, 1L))

  out
}

#' The cells of an ATAC-only experiment, from exactly one source
#'
#' @param fragments String. The fragment file.
#' @param barcodes,cellranger,min_fragments The sources, exactly one of them
#' not `NULL`.
#' @param exclude_chroms Character vector. Chromosomes the `min_fragments`
#' count skips.
#'
#' @returns Character vector of cell barcodes.
#'
#' @keywords internal
.resolve_atac_barcodes <- function(
  fragments,
  barcodes,
  cellranger,
  min_fragments,
  exclude_chroms
) {
  given <- c(!is.null(barcodes), !is.null(cellranger), !is.null(min_fragments))
  if (sum(given) != 1L) {
    stop("Pass exactly one of barcodes, cellranger and min_fragments.")
  }
  out <- if (!is.null(barcodes)) {
    as.character(barcodes)
  } else if (!is.null(cellranger)) {
    read_cellranger_barcodes(cellranger)
  } else {
    counts <- barcode_fragment_counts(fragments, exclude_chroms)
    sort(counts$barcode[counts$n_fragments >= min_fragments], method = "radix")
  }
  if (length(out) == 0L) {
    stop("No cell barcode to load.")
  }

  out
}

#' Write obs, var and the count mode of a fresh ATAC-only store
#'
#' @param object `SingleCellsAtac` class.
#' @param barcodes Character vector. Cells in store order.
#' @param feature_ids Character vector. One id per feature.
#' @param regions List with `chrom`, `start` and `end`.
#' @param lib_size,nnz Numeric vectors. Per-cell totals.
#' @param count_mode String. The count mode.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @returns The loaded `SingleCellsAtac`.
#'
#' @keywords internal
.finish_atac_root <- function(
  object,
  barcodes,
  feature_ids,
  regions,
  lib_size,
  nnz,
  count_mode,
  .verbose
) {
  rust_con <- get_sc_rust_ptr(object)
  rust_con$set_from_file()

  duckdb_con <- get_sc_duckdb(object)
  duckdb_con$populate_obs_from_data.table(data.table::data.table(
    cell_id = barcodes
  ))
  duckdb_con$populate_var_from_data.table(data.table::data.table(
    gene_id = feature_ids,
    chrom = regions$chrom,
    start = as.integer(regions$start),
    end = as.integer(regions$end),
    no_cells_exp = as.integer(rust_con$get_nnz_genes(gene_indices = NULL))
  ))
  duckdb_con$add_data_obs(data.table::data.table(
    nnz = as.numeric(nnz),
    lib_size = as.numeric(lib_size)
  ))
  duckdb_con$set_to_keep_column()

  dir_data <- S7::prop(object, "dir_data")
  .write_count_mode(dir_data, count_mode)
  .stamp_store_cell_ids(dir_data)

  object <- load_existing(object, .verbose = FALSE)

  empty <- sum(lib_size == 0)
  if (empty > 0 && .verbose) {
    message(sprintf(
      paste(
        "%i cells have no count in any feature; they sit at the origin of",
        "the LSI. Filter them (lib_size > 0) unless they are meant to be",
        "there."
      ),
      empty
    ))
  }

  object
}

### atac only loaders ----------------------------------------------------------

#' Load an ATAC-only experiment from a fragment file
#'
#' @description
#' The cells come from exactly one of `barcodes`, `cellranger` or
#' `min_fragments`; with the last, every barcode with at least that many
#' fragments (outside `exclude_chroms`), sorted, as Cell Ranger orders them.
#' See [bixverse::barcode_fragment_counts()] for the knee plot to pick it
#' from. Pass exactly one of `peaks` and `tile_width`.
#'
#' @param object `SingleCellsAtac` class.
#' @param fragments String. The fragment file (BGZF, gzip or plain).
#' @param peaks Optional. A BED file or a data.frame with `chrom`, `start`,
#' `end`.
#' @param tile_width Optional integer. Tile width in bp, instead of `peaks`.
#' @param chrom_sizes Optional named integer vector of chromosome lengths,
#' for the tiles.
#' @param mode String. What one count is, see
#' [bixverse::add_atac_fragments_sc()].
#' @param barcodes Optional character vector. The cell barcodes.
#' @param cellranger Optional string. Cell Ranger ATAC's `singlecell.csv` or
#' filtered h5.
#' @param min_fragments Optional integer. Call cells with at least this many
#' fragments.
#' @param exclude_chroms Character vector. Chromosomes the `min_fragments`
#' count skips.
#' @param max_mem_gb Optional numeric. Memory budget of the feature-major
#' conversion.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @returns The loaded experiment. obs holds `cell_id`, `lib_size`, `nnz`,
#' `to_keep`; var the features as `gene_id` (`chr:start-end`) with `chrom`,
#' `start`, `end` and `no_cells_exp`.
#'
#' @export
load_atac_fragments <- S7::new_generic(
  name = "load_atac_fragments",
  dispatch_args = "object",
  fun = function(
    object,
    fragments,
    peaks = NULL,
    tile_width = NULL,
    chrom_sizes = NULL,
    mode = c("insertion", "paired-insertion", "fragment"),
    barcodes = NULL,
    cellranger = NULL,
    min_fragments = NULL,
    exclude_chroms = c("chrM", "MT", "chrMT", "M"),
    max_mem_gb = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_atac_fragments SingleCellsAtac
#'
#' @export
S7::method(load_atac_fragments, SingleCellsAtac) <- function(
  object,
  fragments,
  peaks = NULL,
  tile_width = NULL,
  chrom_sizes = NULL,
  mode = c("insertion", "paired-insertion", "fragment"),
  barcodes = NULL,
  cellranger = NULL,
  min_fragments = NULL,
  exclude_chroms = c("chrM", "MT", "chrMT", "M"),
  max_mem_gb = NULL,
  .verbose = TRUE
) {
  mode <- match.arg(mode)
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsAtac))
  checkmate::assertFileExists(fragments)
  checkmate::assertChoice(mode, .ATAC_COUNT_MODES)
  checkmate::qassert(barcodes, c("0", "S+"))
  checkmate::qassert(cellranger, c("0", "S1"))
  checkmate::qassert(min_fragments, c("0", "X1[1,)"))
  checkmate::qassert(exclude_chroms, "S*")
  checkmate::qassert(max_mem_gb, c("0", "N1(0,)"))
  checkmate::qassert(.verbose, "B1")

  # function body
  fragments <- path.expand(fragments)
  cells <- .resolve_atac_barcodes(
    fragments = fragments,
    barcodes = barcodes,
    cellranger = cellranger,
    min_fragments = min_fragments,
    exclude_chroms = exclude_chroms
  )
  features <- .atac_features(peaks, tile_width, chrom_sizes)
  res <- .count_store(
    fragments = fragments,
    store_dir = S7::prop(object, "dir_data"),
    barcodes = cells,
    regions = features,
    mode = mode,
    kind = features$kind,
    target_size = .ATAC_TF_IDF_SCALE,
    max_mem_gb = max_mem_gb,
    .verbose = .verbose
  )

  .finish_atac_root(
    object = object,
    barcodes = cells,
    feature_ids = .region_ids(features),
    regions = features,
    lib_size = res$lib_size,
    nnz = res$nnz,
    count_mode = mode,
    .verbose = .verbose
  )
}

#' Load an ATAC-only experiment from a 10x peak matrix h5
#'
#' @param object `SingleCellsAtac` class.
#' @param h5_path String. Cell Ranger ATAC's `filtered_peak_bc_matrix.h5`, or
#' the Cell Ranger ARC h5 (its `Peaks` rows).
#' @param barcodes Optional character vector. Cells to load; `NULL` loads
#' every barcode of the h5.
#' @param max_mem_gb Optional numeric. Memory budget of the feature-major
#' conversion.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @returns The loaded experiment, as [bixverse::load_atac_fragments()]
#' describes it; the counts are insertions.
#'
#' @export
load_atac_tenx <- S7::new_generic(
  name = "load_atac_tenx",
  dispatch_args = "object",
  fun = function(
    object,
    h5_path,
    barcodes = NULL,
    max_mem_gb = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_atac_tenx SingleCellsAtac
#'
#' @export
S7::method(load_atac_tenx, SingleCellsAtac) <- function(
  object,
  h5_path,
  barcodes = NULL,
  max_mem_gb = NULL,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsAtac))
  checkmate::assertFileExists(h5_path)
  checkmate::qassert(barcodes, c("0", "S+"))
  checkmate::qassert(max_mem_gb, c("0", "N1(0,)"))
  checkmate::qassert(.verbose, "B1")

  # function body
  h5_path <- path.expand(h5_path)
  cells <- if (is.null(barcodes)) {
    as.character(rs_h5_read(h5_path, "matrix/barcodes"))
  } else {
    barcodes
  }
  dir_data <- S7::prop(object, "dir_data")
  res <- rs_atac_tenx_peaks_to_store(
    h5_path = h5_path,
    f_path_cells = file.path(dir_data, "counts_cells.bin"),
    f_path_genes = file.path(dir_data, "counts_genes.bin"),
    barcodes = cells,
    max_mem_gb = max_mem_gb,
    verbose = as.integer(.verbose)
  )

  .finish_atac_root(
    object = object,
    barcodes = cells,
    feature_ids = res$name,
    regions = res[c("chrom", "start", "end")],
    lib_size = res$lib_size,
    nnz = res$nnz,
    count_mode = "insertion",
    .verbose = .verbose
  )
}

### merge ----------------------------------------------------------------------

#' Merge per-sample fragment files into one
#'
#' @description
#' As `cellranger-arc aggr` does: each file's barcodes get its prefix and
#' suffix, so the same barcode in two samples stays two cells; give every
#' file its own. The output is a sorted BGZF file with a tabix index, which
#' every ATAC function reads as one experiment.
#'
#' @param fragments Character vector. The fragment files.
#' @param out String. The merged file; the index goes next to it as `.tbi`.
#' @param prefixes Optional character vector. Barcode prefix per file, e.g.
#' `c("s1_", "s2_")`.
#' @param suffixes Optional character vector. Barcode suffix per file, e.g.
#' `c("-1", "-2")`.
#'
#' @returns A data.table with `fragments`, `prefix`, `suffix` and
#' `n_fragments` per input.
#'
#' @export
merge_fragments <- function(fragments, out, prefixes = NULL, suffixes = NULL) {
  checkmate::assertFileExists(fragments)
  checkmate::qassert(out, "S1")
  checkmate::qassert(prefixes, c("0", "S+"))
  checkmate::qassert(suffixes, c("0", "S+"))

  n <- length(fragments)
  if (is.null(prefixes) && is.null(suffixes)) {
    stop(paste(
      "Pass prefixes or suffixes, so barcodes of different files stay",
      "apart."
    ))
  }
  pre <- prefixes %||% rep("", n)
  suf <- suffixes %||% rep("", n)
  if (length(pre) != n || length(suf) != n) {
    stop("prefixes and suffixes need one entry per fragment file.")
  }
  res <- rs_atac_merge_fragments(
    paths = path.expand(fragments),
    prefixes = pre,
    suffixes = suf,
    out = path.expand(out)
  )

  data.table::data.table(
    fragments = fragments,
    prefix = pre,
    suffix = suf,
    n_fragments = res$n_fragments
  )
}

## gene activity ---------------------------------------------------------------

#' Count fragments over gene body plus promoter into a gene activity store
#'
#' @description
#' Signac's `GeneActivity()`: per gene, the longest span of its biotype,
#' under `max_width`, extended `upstream` of the TSS; each overlapping
#' fragment counts once. The result is a [bixverse::SingleCells()] over every
#' cell of `object`, in its order and with its obs and cell filter, so it
#' serves as the query of [bixverse::transfer_labels_seurat_sc()] and of the
#' RNA tools. Methods that assume RNA counts refuse it or warn. Written to the
#' default directory, later cell filters on `object` reach it; reopen it with
#' [bixverse::get_gene_activity()].
#'
#' @param object `SingleCells` class whose cell ids are the fragment
#' barcodes (a `SingleCellsMultiModal` or a `SingleCellsAtac`).
#' @param fragments String. The fragment file (BGZF, gzip or plain).
#' @param gtf String. Gene annotation (GENCODE, Ensembl, or a
#' transcript-only GTF).
#' @param dir_out Optional string. Output directory. `NULL` uses
#' `gene_activity` inside the experiment directory.
#' @param ga_params List. Output of [bixverse::params_atac_gene_activity()].
#' @param target_size Optional numeric. Target of the log-CPM norm. `NULL`
#' uses the median RNA `lib_size` of the kept cells, as the Signac vignettes
#' do, or 1e4 on an ATAC-only experiment, whose `lib_size` counts peaks.
#' @param barcode_map Optional. Cell id to fragment barcode, where they
#' differ: a named character vector or a function.
#' @param overwrite Boolean. Replace an existing store at `dir_out`.
#' @param max_mem_gb Optional numeric. Memory budget of the gene-major
#' conversion.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @returns The gene activity experiment as a loaded `SingleCells`. var holds
#' `gene_id` (the GTF gene ID without its version, as Cell Ranger keys the
#' RNA), `gene_name` (made unique), `chrom`, `start` and `end`; obs's
#' `lib_size` and `nnz` are the gene activity's.
#'
#' @export
#'
#' @references Stuart, et al., Nat Methods, 2021
gene_activity_sc <- S7::new_generic(
  name = "gene_activity_sc",
  dispatch_args = "object",
  fun = function(
    object,
    fragments,
    gtf,
    dir_out = NULL,
    ga_params = params_atac_gene_activity(),
    target_size = NULL,
    barcode_map = NULL,
    overwrite = FALSE,
    max_mem_gb = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method gene_activity_sc SingleCells
#'
#' @export
S7::method(gene_activity_sc, SingleCells) <- function(
  object,
  fragments,
  gtf,
  dir_out = NULL,
  ga_params = params_atac_gene_activity(),
  target_size = NULL,
  barcode_map = NULL,
  overwrite = FALSE,
  max_mem_gb = NULL,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::assertFileExists(fragments)
  checkmate::assertFileExists(gtf)
  checkmate::qassert(dir_out, c("0", "S1"))
  assertAtacGeneActivityParams(ga_params)
  checkmate::qassert(target_size, c("0", "N1(0,)"))
  checkmate::qassert(overwrite, "B1")
  checkmate::qassert(max_mem_gb, c("0", "N1(0,)"))
  checkmate::qassert(.verbose, "B1")

  # function body
  out <- dir_out %||% file.path(S7::prop(object, "dir_data"), "gene_activity")
  out <- path.expand(out)
  if (file.exists(file.path(out, "counts_genes.bin")) && !overwrite) {
    stop(sprintf(
      "A store already exists at '%s'. Pass overwrite = TRUE to replace it.",
      out
    ))
  }

  obs <- get_sc_obs(object, filtered = FALSE)
  if (is.null(target_size)) {
    if (S7::S7_inherits(object, SingleCellsAtac)) {
      target_size <- .ATAC_TF_IDF_SCALE
    } else {
      if (!"lib_size" %in% names(obs)) {
        stop("obs has no lib_size; pass target_size.")
      }
      kept <- get_cells_to_keep(object) + 1L
      target_size <- stats::median(obs$lib_size[kept])
    }
  }

  regions <- rs_atac_gene_activity_regions(
    path = path.expand(gtf),
    upstream = ga_params$upstream,
    downstream = ga_params$downstream,
    max_width = ga_params$max_width,
    biotypes = ga_params$biotypes %||% character()
  )
  if (dir.exists(out)) {
    unlink(
      file.path(out, c("counts_cells.bin", "counts_genes.bin", "sc_duckdb.db"))
    )
  }
  res <- .count_store(
    fragments = path.expand(fragments),
    store_dir = out,
    barcodes = .fragment_barcodes(obs$cell_id, barcode_map),
    regions = regions,
    mode = "fragment",
    kind = "gene_activity",
    target_size = target_size,
    max_mem_gb = max_mem_gb,
    .verbose = .verbose
  )

  # unversioned, as Cell Ranger keys the RNA, so RNA and gene activity share
  # gene_id; `_PAR_Y` survives and keeps the chrY copies apart
  gene_id <- sub("\\.\\d+", "", regions$gene_id)
  if (anyDuplicated(gene_id)) {
    stop("Gene IDs collide once their versions are stripped.")
  }

  ga <- SingleCells(dir_data = out)
  rust_con <- get_sc_rust_ptr(ga)
  rust_con$set_from_file()
  duckdb_con <- get_sc_duckdb(ga)
  data.table::set(obs, j = "lib_size", value = as.numeric(res$lib_size))
  data.table::set(obs, j = "nnz", value = as.numeric(res$nnz))
  drop <- intersect(c("cell_idx", "to_keep"), names(obs))
  if (length(drop)) {
    obs[, (drop) := NULL]
  }
  duckdb_con$populate_obs_from_data.table(obs)
  duckdb_con$populate_var_from_data.table(data.table::data.table(
    gene_id = gene_id,
    gene_name = regions$name,
    chrom = regions$chrom,
    start = regions$start,
    end = regions$end,
    no_cells_exp = as.integer(rust_con$get_nnz_genes(gene_indices = NULL))
  ))
  duckdb_con$set_to_keep_column()
  .stamp_store_cell_ids(out)

  ga <- load_existing(ga, .verbose = FALSE)
  set_cells_to_keep(ga, get_cell_names(object, filtered = TRUE))
}

#' Open an experiment's gene activity store
#'
#' @description
#' Opens the store [bixverse::gene_activity_sc()] wrote, with the cell filter
#' of `object`.
#'
#' @param object `SingleCells` class the gene activity was counted for.
#' @param dir_data Optional string. Where it was written; `NULL` for the
#' default `<dir_data>/gene_activity`.
#'
#' @returns The gene activity experiment as a `SingleCells`. Errors if there
#' is none, or if it was counted for other cells (a re-ingest since).
#'
#' @export
get_gene_activity <- function(object, dir_data = NULL) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(dir_data, c("0", "S1"))

  out <- dir_data %||%
    file.path(S7::prop(object, "dir_data"), "gene_activity")
  out <- path.expand(out)
  if (!file.exists(file.path(out, "counts_cells.bin"))) {
    stop(sprintf(
      "No gene activity store at %s. Run gene_activity_sc().",
      out
    ))
  }
  theirs <- rs_cell_id_hash(file.path(out, "counts_cells.bin"))
  ours <- rs_cell_id_hash(get_rust_count_cell_f_path(object))
  if (!is.null(theirs) && !is.null(ours) && theirs != ours) {
    stop(sprintf(
      paste(
        "The gene activity at %s was counted for other cells. Run",
        "gene_activity_sc(..., overwrite = TRUE) again."
      ),
      out
    ))
  }

  ga <- load_existing(SingleCells(dir_data = out), .verbose = FALSE)
  set_cells_to_keep(ga, get_cell_names(object, filtered = TRUE))
}

## qc --------------------------------------------------------------------------

#' Per-cell QC from the fragment file
#'
#' @description
#' Adds `atac_n_fragments`, `atac_nucleosome_signal`, `atac_tss_enrichment`
#' (Signac's fast definition) and `atac_mito_fraction` to obs for every cell,
#' plus `atac_frip` on a peak store: its library size over the fragments, in
#' the store's count unit, so with insertion counts the share of insertions in
#' peaks (two per fragment) and otherwise Signac's fragments in peaks over
#' fragments.
#'
#' The nucleosome signal is mononucleosomal over nucleosome-free fragments,
#' as in Signac: infinite for a cell with mononucleosomal but no
#' nucleosome-free fragments (a bad cell, which any upper threshold drops),
#' and 0 rather than Signac's `NaN` for a cell with neither.
#'
#' @param object `SingleCellsAtac` or `SingleCellsMultiModal` class.
#' @param fragments String. The fragment file.
#' @param gtf String. Gene annotation for the TSSs.
#' @param biotypes Optional character vector. Gene biotypes whose TSSs count
#' (Signac: protein coding). `NULL` keeps all.
#' @param mito_chroms Character vector. Names of the mitochondrial genome.
#' @param barcode_map Optional. Cell id to fragment barcode, where they
#' differ.
#'
#' @returns The object with the QC columns added to obs.
#'
#' @export
#'
#' @references Stuart, et al., Nat Methods, 2021
fragment_qc_sc <- S7::new_generic(
  name = "fragment_qc_sc",
  dispatch_args = "object",
  fun = function(
    object,
    fragments,
    gtf,
    biotypes = "protein_coding",
    mito_chroms = c("chrM", "MT", "chrMT", "M"),
    barcode_map = NULL
  ) {
    S7::S7_dispatch()
  }
)

#' @method fragment_qc_sc AtacHost
S7::method(fragment_qc_sc, AtacHost) <- function(
  object,
  fragments,
  gtf,
  biotypes = "protein_coding",
  mito_chroms = c("chrM", "MT", "chrMT", "M"),
  barcode_map = NULL
) {
  # checks
  checkmate::assertFileExists(fragments)
  checkmate::qassert(biotypes, c("0", "S+"))
  checkmate::qassert(mito_chroms, "S*")

  # function body
  fragments <- path.expand(fragments)
  .check_fragments(fragments, TRUE)
  tss <- .atac_tss(gtf, biotypes)
  barcodes <- .fragment_barcodes(
    get_cell_names(object, filtered = FALSE),
    barcode_map
  )
  qc <- rs_atac_fragment_qc(
    fragments = fragments,
    barcodes = barcodes,
    tss_chrom = tss$chrom,
    tss_pos = tss$pos,
    mito_chroms = mito_chroms
  )
  .report_coverage(qc, length(barcodes), "TSSs")

  ns <- qc$nucleosome_signal
  ns[is.nan(ns)] <- 0
  new_obs <- data.table::data.table(
    atac_n_fragments = qc$n_fragments,
    atac_nucleosome_signal = ns,
    atac_tss_enrichment = qc$tss_enrichment,
    atac_mito_fraction = qc$mito_fraction
  )
  frip <- .atac_frip(object, qc$n_fragments)
  if (!is.null(frip)) {
    data.table::set(new_obs, j = "atac_frip", value = frip)
  }
  get_sc_duckdb(object)$add_data_obs(new_obs)

  object
}

#' TSS insertion profile and fragment length histogram of the kept cells
#'
#' @description
#' One pass over the fragment file, per level of `group_by` (SnapATAC2's
#' `TSS_profile` and `frag_size_distr`, which these reproduce): insertions are
#' a fragment's start and `end - 1`, offsets are strand aware, positive
#' downstream.
#'
#' @param object `SingleCellsAtac` or `SingleCellsMultiModal` class.
#' @param fragments String. The fragment file.
#' @param gtf String. Gene annotation for the TSSs.
#' @param group_by Optional string. obs column to split by; `NULL` for one
#' group, `"all"`.
#' @param half_width Integer. Profile half width in bp.
#' @param max_len Integer. Longest fragment length in the histogram.
#' @param flank Integer. Bases at each end of the profile whose mean is the
#' background of the enrichment.
#' @param biotypes Optional character vector. Gene biotypes whose TSSs
#' count. `NULL` keeps all.
#' @param barcode_map Optional. Cell id to fragment barcode, where they
#' differ.
#'
#' @returns A list with
#' \itemize{
#'   \item tss - data.table with `group`, `offset`, `insertions` and
#'   `enrichment` (insertions over the mean of both flanks).
#'   \item lengths - data.table with `group`, `length`, `fragments` and
#'   `fraction` (of the group's fragments up to `max_len`).
#' }
#'
#' @export
fragment_profiles_sc <- S7::new_generic(
  name = "fragment_profiles_sc",
  dispatch_args = "object",
  fun = function(
    object,
    fragments,
    gtf,
    group_by = NULL,
    half_width = 2000L,
    max_len = 1000L,
    flank = 100L,
    biotypes = "protein_coding",
    barcode_map = NULL
  ) {
    S7::S7_dispatch()
  }
)

#' @method fragment_profiles_sc AtacHost
S7::method(fragment_profiles_sc, AtacHost) <- function(
  object,
  fragments,
  gtf,
  group_by = NULL,
  half_width = 2000L,
  max_len = 1000L,
  flank = 100L,
  biotypes = "protein_coding",
  barcode_map = NULL
) {
  # checks
  checkmate::assertFileExists(fragments)
  checkmate::qassert(group_by, c("0", "S1"))
  checkmate::qassert(half_width, "I1[1,)")
  checkmate::qassert(max_len, "I1[1,)")
  checkmate::qassert(flank, "I1[1,)")
  checkmate::assertTRUE(flank <= half_width)
  checkmate::qassert(biotypes, c("0", "S+"))

  # function body
  kept <- get_cells_to_keep(object)
  groups <- if (is.null(group_by)) {
    list(codes = integer(length(kept)), levels = "all")
  } else {
    .kept_obs_codes(object, group_by)
  }
  tss <- .atac_tss(gtf, biotypes)
  barcodes <- .fragment_barcodes(
    get_cell_names(object, filtered = TRUE),
    barcode_map
  )
  res <- rs_atac_fragment_profiles(
    fragments = path.expand(fragments),
    barcodes = barcodes,
    groups = groups$codes,
    tss_chrom = tss$chrom,
    tss_pos = tss$pos,
    tss_strand = tss$strand,
    half_width = half_width,
    max_len = max_len
  )
  .report_coverage(res, length(barcodes), "TSSs")

  prof <- res$tss
  width <- 2L * half_width + 1L
  background <- rowMeans(cbind(
    prof[, seq_len(flank), drop = FALSE],
    prof[, (width - flank + 1L):width, drop = FALSE]
  ))
  enrich <- prof / background
  n_groups <- length(groups$levels)
  tss_dt <- data.table::data.table(
    group = factor(rep(groups$levels, each = width), levels = groups$levels),
    offset = rep(seq(-half_width, half_width), times = n_groups),
    insertions = as.vector(t(prof)),
    enrichment = as.vector(t(enrich))
  )

  lengths <- res$lengths[, -1L, drop = FALSE]
  frac <- lengths / rowSums(lengths)
  len_dt <- data.table::data.table(
    group = factor(rep(groups$levels, each = max_len), levels = groups$levels),
    length = rep(seq_len(max_len), times = n_groups),
    fragments = as.vector(t(lengths)),
    fraction = as.vector(t(frac))
  )

  list(tss = tss_dt, lengths = len_dt)
}

## peak calling ----------------------------------------------------------------

#' Call peaks on the kept cells, pooled or per group
#'
#' @description
#' A native port of MACS3 `callpeak --nomodel` without a control track, one
#' call per group, as Signac's `CallPeaks(group.by = ...)` runs MACS. Against
#' MACS3 3.1.0 the peak coordinates, summits and pileups are exact and the
#' scores agree to MACS's printed precision. Pass the result as `peaks` to
#' [bixverse::add_atac_fragments_sc()] or
#' [bixverse::load_atac_fragments()] to count against it.
#'
#' With a tabix index next to the fragment file (`.tbi`, which Cell Ranger
#' writes), the file is read one chromosome at a time and the sorted tags
#' wait in a scratch file (4 bytes per tag) between passes. Without an index
#' every tag is held at once while reading.
#'
#' @param object `SingleCells` class; its kept cells are pooled.
#' @param fragments String. The fragment file.
#' @param group_by Optional string. obs column to call per level of
#' (clusters, samples); `NULL` pools all kept cells into one group, `"all"`.
#' @param combine Boolean. Merge overlapping peaks across groups into one
#' set, with the groups each was called in. `FALSE` returns every group's
#' peaks.
#' @param peak_params List. Output of [bixverse::params_atac_call_peaks()].
#' @param barcode_map Optional. Cell id to fragment barcode, where they
#' differ.
#' @param max_mem_gb Optional numeric. Tag budget of the reading pass; `NULL`
#' reads one chromosome at a time.
#' @param tmp_dir Optional string. Directory of the scratch file; `NULL` for
#' the system temporary directory.
#'
#' @returns A data.table. With `combine`: `chrom`, `start`, `end` (0-based,
#' half-open) and `peak_called_in`, the comma-separated groups. Without: one
#' row per group and peak, with `group`, `summit` and the summit's `pileup`,
#' `pscore` (`-log10 p`), `fold` and `qscore` (`-log10 q`).
#'
#' @export
#'
#' @references Zhang, et al., Genome Biol, 2008
call_peaks_sc <- S7::new_generic(
  name = "call_peaks_sc",
  dispatch_args = "object",
  fun = function(
    object,
    fragments,
    group_by = NULL,
    combine = TRUE,
    peak_params = params_atac_call_peaks(),
    barcode_map = NULL,
    max_mem_gb = NULL,
    tmp_dir = NULL
  ) {
    S7::S7_dispatch()
  }
)

#' @method call_peaks_sc SingleCells
#'
#' @export
S7::method(call_peaks_sc, SingleCells) <- function(
  object,
  fragments,
  group_by = NULL,
  combine = TRUE,
  peak_params = params_atac_call_peaks(),
  barcode_map = NULL,
  max_mem_gb = NULL,
  tmp_dir = NULL
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::assertFileExists(fragments)
  checkmate::qassert(group_by, c("0", "S1"))
  checkmate::qassert(combine, "B1")
  assertAtacCallPeaksParams(peak_params)
  checkmate::qassert(max_mem_gb, c("0", "N1(0,)"))
  checkmate::qassert(tmp_dir, c("0", "S1"))

  # function body
  fragments <- path.expand(fragments)
  .check_fragments(fragments, TRUE)
  kept <- get_cells_to_keep(object)
  groups <- if (is.null(group_by)) {
    list(codes = integer(length(kept)), levels = "all")
  } else {
    .kept_obs_codes(object, group_by)
  }
  res <- rs_atac_call_peaks(
    fragments = fragments,
    barcodes = .fragment_barcodes(
      get_cell_names(object, filtered = TRUE),
      barcode_map
    ),
    groups = groups$codes,
    tag_source = peak_params$tag_source,
    keep_dup = peak_params$keep_dup,
    shift = peak_params$shift,
    extsize = peak_params$extsize,
    genome_size = peak_params$genome_size,
    llocal = peak_params$llocal,
    q_value = peak_params$q_value,
    min_length = peak_params$min_length,
    max_gap = peak_params$max_gap,
    max_mem_gb = max_mem_gb,
    tmp_dir = tmp_dir
  )
  peaks <- data.table::data.table(
    group = factor(groups$levels[res$group + 1L], levels = groups$levels),
    chrom = res$chrom,
    start = res$start,
    end = res$end,
    summit = res$summit,
    pileup = res$pileup,
    pscore = res$pscore,
    fold = res$fold,
    qscore = res$qscore
  )
  if (nrow(peaks) == 0L) {
    warning(paste(
      "No peak was called. Check the barcodes (barcode_map) and the",
      "q_value."
    ))
  }
  if (!combine) {
    return(peaks)
  }

  merged <- rs_atac_reduce_regions(
    chrom = peaks$chrom,
    start = peaks$start,
    end = peaks$end
  )
  group <- as.character(peaks$group)
  called_in <- purrr::map_chr(merged$members, \(members) {
    present <- unique(group[members + 1L])
    paste(groups$levels[groups$levels %in% present], collapse = ",")
  })

  data.table::data.table(
    chrom = merged$chrom,
    start = merged$start,
    end = merged$end,
    peak_called_in = called_in
  )
}

#' Drop the peaks that overlap a blacklist region
#'
#' @param peaks data.frame. Regions with `chrom`, `start`, `end`, e.g. from
#' [bixverse::call_peaks_sc()].
#' @param blacklist A BED file (plain or gzip) or a data.frame of regions,
#' e.g. the ENCODE blacklist.
#'
#' @returns The peaks that touch no blacklist region, in order.
#'
#' @export
remove_blacklisted <- function(peaks, blacklist) {
  checkmate::assertDataFrame(peaks)

  p <- .regions_from(peaks)
  b <- .regions_from(blacklist)
  hit <- rs_atac_overlaps_any(
    chrom = p$chrom,
    start = p$start,
    end = p$end,
    other_chrom = b$chrom,
    other_start = b$start,
    other_end = b$end
  )

  data.table::as.data.table(peaks)[!hit]
}

#' The closest gene of every ATAC feature
#'
#' @description
#' Signac's `ClosestFeature()`. Written to the ATAC feature table as
#' `closest_gene`, `closest_gene_id`, `distance_to_gene` (bases to the gene
#' body, 0 inside or touching it) and `distance_to_tss` (feature centre to
#' TSS, positive downstream).
#'
#' @param object `SingleCellsAtac` or `SingleCellsMultiModal` class.
#' @param gtf String. Gene annotation.
#' @param biotypes Optional character vector. Biotypes to consider. `NULL`
#' keeps all.
#'
#' @returns The object with the four columns added to the ATAC feature
#' table.
#'
#' @export
annotate_peaks_sc <- S7::new_generic(
  name = "annotate_peaks_sc",
  dispatch_args = "object",
  fun = function(object, gtf, biotypes = "protein_coding") {
    S7::S7_dispatch()
  }
)

#' @method annotate_peaks_sc AtacHost
S7::method(annotate_peaks_sc, AtacHost) <- function(
  object,
  gtf,
  biotypes = "protein_coding"
) {
  # checks
  checkmate::assertFileExists(gtf)
  checkmate::qassert(biotypes, c("0", "S+"))
  .require_atac(object)

  # function body
  var <- .atac_var(object, cols = c("chrom", "start", "end"))
  res <- rs_atac_closest_genes(
    gtf = path.expand(gtf),
    biotypes = biotypes %||% character(),
    chrom = var$chrom,
    start = as.integer(var$start),
    end = as.integer(var$end)
  )
  .add_atac_var(
    object,
    data.table::data.table(
      closest_gene = res$gene_name,
      closest_gene_id = res$gene_id,
      distance_to_gene = res$distance,
      distance_to_tss = res$tss_distance
    )
  )

  object
}

## multiplets ------------------------------------------------------------------

#' Flag multiplets from fragment overlaps (AMULET)
#'
#' @description
#' A diploid nucleus yields at most two overlapping fragments at a locus. Per
#' cell, the loci with more are counted (autosomes, restriction regions
#' dropped), sites shared by too many cells are dropped as repeats, and the
#' count is tested against the mean (Poisson, Benjamini-Hochberg).
#' Implemented from the paper; on the test data it reproduces the reference
#' implementation's counts and p-values. Needs no clustering and no RNA, so
#' it works for ATAC-only experiments. Run it per library: multiplets are a
#' property of one capture.
#'
#' Adds `amulet_n_overlaps`, `amulet_n_sites`, `amulet_p`, `amulet_q` and
#' `amulet_multiplet` to obs for every cell of the store.
#'
#' @param object `SingleCells` class; every cell of its store is tested.
#' @param fragments String. The fragment file.
#' @param restrict_beds Character vector. BED files of regions to ignore:
#' repeats, segmental duplications, blacklists. The paper uses all three.
#' @param amulet_params List. Output of [bixverse::params_atac_amulet()].
#' @param barcode_map Optional. Cell id to fragment barcode, where they
#' differ.
#' @param max_mem_gb Optional numeric. Read budget; `NULL` reads one
#' chromosome at a time (with a tabix index).
#'
#' @returns The object with the AMULET columns added to obs.
#'
#' @export
#'
#' @references Thibodeau, et al., Genome Biol, 2021
amulet_sc <- S7::new_generic(
  name = "amulet_sc",
  dispatch_args = "object",
  fun = function(
    object,
    fragments,
    restrict_beds = character(),
    amulet_params = params_atac_amulet(),
    barcode_map = NULL,
    max_mem_gb = NULL
  ) {
    S7::S7_dispatch()
  }
)

#' @method amulet_sc SingleCells
#'
#' @export
S7::method(amulet_sc, SingleCells) <- function(
  object,
  fragments,
  restrict_beds = character(),
  amulet_params = params_atac_amulet(),
  barcode_map = NULL,
  max_mem_gb = NULL
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::assertFileExists(fragments)
  checkmate::qassert(restrict_beds, "S*")
  if (length(restrict_beds)) {
    checkmate::assertFileExists(restrict_beds)
  }
  assertAtacAmuletParams(amulet_params)
  checkmate::qassert(max_mem_gb, c("0", "N1(0,)"))

  # function body
  fragments <- path.expand(fragments)
  .check_fragments(fragments, TRUE)
  res <- rs_atac_amulet(
    fragments = fragments,
    barcodes = .fragment_barcodes(
      get_cell_names(object, filtered = FALSE),
      barcode_map
    ),
    restrict_beds = path.expand(restrict_beds),
    max_overlap = amulet_params$max_overlap,
    max_fragment_len = amulet_params$max_fragment_len,
    repeat_q = amulet_params$repeat_q,
    multiplet_q = amulet_params$multiplet_q,
    chroms = amulet_params$chroms,
    max_mem_gb = max_mem_gb,
    count_trailing_loci = amulet_params$count_trailing_loci
  )
  get_sc_duckdb(object)$add_data_obs(data.table::data.table(
    amulet_n_overlaps = res$n_overlaps,
    amulet_n_sites = res$n_sites,
    amulet_p = res$p_value,
    amulet_q = res$q_value,
    amulet_multiplet = res$multiplet
  ))

  object
}

## tf-idf and features ---------------------------------------------------------

#' Switch the ATAC store to TF-IDF over the kept cells
#'
#' @description
#' Signac's method 1, `log1p(x / lib_size * scale * idf)`, with the IDF over
#' the kept cells. Header metadata only; no count is rewritten. Run again
#' after the kept cells change. Everything in the ATAC cache (LSI,
#' embeddings, kNN) was built on the old norm and is dropped; rerun
#' [bixverse::calculate_lsi_sc()].
#'
#' @param object `SingleCellsAtac` or `SingleCellsMultiModal` class with
#' ATAC.
#' @param scale Numeric. Scale factor of the term frequencies.
#'
#' @returns The object with the store switched to TF-IDF.
#'
#' @export
#'
#' @references Stuart, et al., Nat Methods, 2021
run_tfidf_sc <- S7::new_generic(
  name = "run_tfidf_sc",
  dispatch_args = "object",
  fun = function(object, scale = 1e4) {
    S7::S7_dispatch()
  }
)

#' @method run_tfidf_sc AtacHost
S7::method(run_tfidf_sc, AtacHost) <- function(object, scale = 1e4) {
  # checks
  checkmate::qassert(scale, "N1(0,)")
  info <- .require_atac(object)

  # function body
  paths <- .atac_paths(object)
  rs_atac_set_tf_idf(
    f_path_cells = paths$cells,
    f_path_genes = paths$genes,
    cells = get_cells_to_keep(object),
    scale = scale
  )
  info$tf_idf <- TRUE
  S7::prop(object, "atac_info") <- info
  object <- .sc_cache_put(object, new_sc_cache(), .atac_modality(object))

  object
}

#' Find the most accessible features
#'
#' @description
#' Signac's `FindTopFeatures()`. Totals are summed over the kept cells; the
#' ATAC feature table gains `count`, `percentile` and `top_feature`, which
#' [bixverse::calculate_lsi_sc()] then uses.
#'
#' @param object `SingleCellsAtac` or `SingleCellsMultiModal` class with
#' ATAC.
#' @param min_cutoff `"qN"` keeps features above the N-th percentile of the
#' totals, a number those with more counts, `NULL` or `"all"` all.
#'
#' @returns The object with the columns added to the ATAC feature table.
#'
#' @export
find_top_features_sc <- S7::new_generic(
  name = "find_top_features_sc",
  dispatch_args = "object",
  fun = function(object, min_cutoff = "q5") {
    S7::S7_dispatch()
  }
)

#' @method find_top_features_sc AtacHost
S7::method(find_top_features_sc, AtacHost) <- function(
  object,
  min_cutoff = "q5"
) {
  # checks
  checkmate::qassert(min_cutoff, c("0", "S1", "N1[0,)"))
  .require_atac(object)

  # function body
  res <- rs_atac_feature_totals(
    f_path_genes = .atac_paths(object)$genes,
    cells = get_cells_to_keep(object),
    min_cutoff = if (is.null(min_cutoff)) NULL else as.character(min_cutoff)
  )
  .add_atac_var(
    object,
    data.table::data.table(
      count = res$count,
      percentile = res$percentile,
      top_feature = res$top_feature
    )
  )

  object
}

## lsi -------------------------------------------------------------------------

#' LSI over the TF-IDF of the kept cells
#'
#' @description
#' Signac's `RunSVD()`. The embedding (components z-scored over the cells)
#' goes to the ATAC cache as `"lsi"`, for [bixverse::find_neighbours_sc()],
#' [bixverse::harmony_sc()], [bixverse::umap_sc()] and
#' [bixverse::generate_wnn_graph_sc()] with `modality = "atac"`,
#' `embd_to_use = "lsi"`. The loadings, singular values and features are kept
#' for [bixverse::project_lsi_sc()].
#'
#' @param object `SingleCellsAtac` or `SingleCellsMultiModal` class, after
#' [bixverse::run_tfidf_sc()].
#' @param features Optional integer vector. 1-based features. `NULL` takes
#' `top_feature` from the feature table if [bixverse::find_top_features_sc()]
#' ran, else all.
#' @param lsi_params List. Output of [bixverse::params_atac_lsi()].
#' @param seed Integer. Seed of the solver.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @returns The object with the `"lsi"` embedding and the LSI fit.
#'
#' @export
#'
#' @references Stuart, et al., Nat Methods, 2021
calculate_lsi_sc <- S7::new_generic(
  name = "calculate_lsi_sc",
  dispatch_args = "object",
  fun = function(
    object,
    features = NULL,
    lsi_params = params_atac_lsi(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calculate_lsi_sc AtacHost
S7::method(calculate_lsi_sc, AtacHost) <- function(
  object,
  features = NULL,
  lsi_params = params_atac_lsi(),
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::qassert(features, c("0", "X+[1,)"))
  assertAtacLsiParams(lsi_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")
  info <- .require_atac(object)
  if (!info$tf_idf) {
    stop("Run run_tfidf_sc() before calculate_lsi_sc().")
  }

  # function body
  features <- if (is.null(features)) {
    var <- .atac_var(object)
    if ("top_feature" %in% names(var)) {
      which(var$top_feature) - 1L
    } else {
      seq_len(nrow(var)) - 1L
    }
  } else {
    as.integer(features) - 1L
  }
  paths <- .atac_paths(object)
  res <- rs_atac_lsi(
    f_path_cells = paths$cells,
    f_path_genes = paths$genes,
    cells = get_cells_to_keep(object),
    features = features,
    n_comp = lsi_params$n_comp,
    solver = lsi_params$solver,
    scale_embeddings = TRUE,
    scale_max = lsi_params$scale_max,
    seed = seed,
    verbose = as.integer(.verbose)
  )

  embd <- res$embedding
  dropped <- integer()
  if (lsi_params$drop_first) {
    embd <- embd[, -1L, drop = FALSE]
    dropped <- 0L
  }
  colnames(embd) <- sprintf(
    "LSI_%i",
    setdiff(seq_len(ncol(res$embedding)), dropped + 1L)
  )
  object <- set_embedding(
    object,
    embd = embd,
    name = "lsi",
    modality = .atac_modality(object)
  )

  info$depth_cor <- res$depth_cor
  info$lsi_dropped <- dropped
  info$lsi_features <- features
  info$lsi_loadings <- res$loadings
  info$lsi_singular_values <- res$singular_values
  S7::prop(object, "atac_info") <- info

  object
}

#' Project the query's kept cells onto the reference LSI
#'
#' @description
#' Seurat's `ProjectSVD(mode = "lsi")` as
#' `FindTransferAnchors(reduction = "lsiproject")` runs it: the query's own
#' TF-IDF times the reference loadings, divided by the singular values, with
#' no centring or scaling. The reference `"lsi"` embedding is z-scored, the
#' projection is not; L2 normalise both per cell before comparing them, as
#' Seurat does. Components the reference dropped are dropped here too. Either
#' side may be multi-modal or ATAC-only.
#'
#' @param reference `SingleCellsAtac` or `SingleCellsMultiModal` class, after
#' [bixverse::calculate_lsi_sc()].
#' @param query `SingleCellsAtac` or `SingleCellsMultiModal` class to
#' project, after [bixverse::run_tfidf_sc()], counted over the reference
#' features.
#'
#' @returns The query with the projection in its ATAC cache as
#' `"lsi_projected"`.
#'
#' @export
project_lsi_sc <- S7::new_generic(
  name = "project_lsi_sc",
  dispatch_args = "reference",
  fun = function(reference, query) {
    S7::S7_dispatch()
  }
)

#' @method project_lsi_sc AtacHost
S7::method(project_lsi_sc, AtacHost) <- function(reference, query) {
  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(query, SingleCellsAtac) ||
      S7::S7_inherits(query, SingleCellsMultiModal)
  )
  ref <- .require_atac(reference)
  if (is.null(ref$lsi_loadings)) {
    stop("Run calculate_lsi_sc() on the reference first.")
  }
  if (!.require_atac(query)$tf_idf) {
    stop("Run run_tfidf_sc() on the query first.")
  }

  # function body
  key <- \(v) paste(v$chrom, v$start, v$end, sep = "\r")
  ref_keys <- key(.atac_var(reference, cols = c("chrom", "start", "end")))
  query_keys <- key(.atac_var(query, cols = c("chrom", "start", "end")))
  wanted <- ref_keys[ref$lsi_features + 1L]
  lookup <- match(wanted, query_keys)
  if (anyNA(lookup)) {
    first <- strsplit(wanted[is.na(lookup)][1], "\r", fixed = TRUE)[[1]]
    stop(sprintf(
      paste(
        "%i reference LSI features are not in the query, e.g. %s:%s-%s.",
        "Count the query over the reference features."
      ),
      sum(is.na(lookup)),
      first[1],
      first[2],
      first[3]
    ))
  }

  paths <- .atac_paths(query)
  embd <- rs_atac_project_lsi(
    f_path_cells = paths$cells,
    f_path_genes = paths$genes,
    cells = get_cells_to_keep(query),
    features = lookup - 1L,
    loadings = ref$lsi_loadings,
    singular_values = ref$lsi_singular_values
  )
  keep <- setdiff(seq_len(ncol(embd)), ref$lsi_dropped + 1L)
  embd <- embd[, keep, drop = FALSE]
  colnames(embd) <- sprintf("LSI_%i", keep)

  set_embedding(
    query,
    embd = embd,
    name = "lsi_projected",
    modality = .atac_modality(query)
  )
}

#' rLSI integration
#'
#' @description
#' Signac's `FindIntegrationAnchors(reduction = "rlsi")` followed by
#' `IntegrateEmbeddings()`. Per batch: IDF over its own cells, TF-IDF and LSI
#' of `features`; anchors between every batch pair in the reciprocal LSI
#' spaces; then `embedding` is corrected along the anchor tree. Written to
#' the ATAC cache as `"rlsi"`.
#'
#' Signac integrates all LSI components and leaves component 1 out
#' afterwards: `calculate_lsi_sc(x, lsi_params = params_atac_lsi(drop_first =
#' FALSE))`, `run_rlsi_sc(x, ...)`, then columns 2 onwards of `"rlsi"`. With
#' the default `drop_first = TRUE` the depth component is left out of the
#' integration too, which changes the correction weights.
#'
#' @param object `SingleCellsAtac` or `SingleCellsMultiModal` class, after
#' [bixverse::calculate_lsi_sc()].
#' @param batch_column String. obs column with the batch of each cell.
#' @param embedding String. The ATAC embedding to integrate.
#' @param features Optional integer vector. 1-based ATAC features of the
#' per-batch LSIs. `NULL` uses those of [bixverse::calculate_lsi_sc()].
#' @param rlsi_params List. Output of [bixverse::params_atac_rlsi()].
#' @param seed Integer. Seed of the solvers and the kNN backend.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @returns The object with the `"rlsi"` embedding.
#'
#' @export
#'
#' @references Stuart, et al., Nat Methods, 2021
run_rlsi_sc <- S7::new_generic(
  name = "run_rlsi_sc",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    embedding = "lsi",
    features = NULL,
    rlsi_params = params_atac_rlsi(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method run_rlsi_sc AtacHost
S7::method(run_rlsi_sc, AtacHost) <- function(
  object,
  batch_column,
  embedding = "lsi",
  features = NULL,
  rlsi_params = params_atac_rlsi(),
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(embedding, "S1")
  checkmate::qassert(features, c("0", "X+[1,)"))
  assertAtacRlsiParams(rlsi_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")
  info <- .require_atac(object)

  # function body
  modality <- .atac_modality(object)
  embd <- get_embedding(object, embd_name = embedding, modality = modality)
  features <- if (is.null(features)) {
    if (is.null(info$lsi_features)) {
      stop("Run calculate_lsi_sc() first or pass features.")
    }
    info$lsi_features
  } else {
    as.integer(features) - 1L
  }
  batch <- .kept_obs_codes(object, batch_column)
  if (length(batch$levels) < 2L) {
    stop(sprintf("Column '%s' has fewer than two batches.", batch_column))
  }
  paths <- .atac_paths(object)
  knn_params <- rlsi_params[
    setdiff(
      names(rlsi_params),
      c(
        "dims",
        "n_comp",
        "k_anchor",
        "k_score",
        "k_weight",
        "sd",
        "l2_norm",
        "solver"
      )
    )
  ]
  out <- rs_atac_rlsi(
    f_path_cells = paths$cells,
    f_path_genes = paths$genes,
    cells = get_cells_to_keep(object),
    features = features,
    batch = batch$codes,
    embedding = unname(embd),
    dims = rlsi_params$dims - 1L,
    n_comp = rlsi_params$n_comp,
    solver = rlsi_params$solver,
    k_anchor = rlsi_params$k_anchor,
    k_score = rlsi_params$k_score,
    k_weight = rlsi_params$k_weight,
    sd = rlsi_params$sd,
    l2_norm = rlsi_params$l2_norm,
    knn_params = knn_params,
    seed = seed,
    verbose = as.integer(.verbose)
  )
  colnames(out) <- colnames(embd)

  set_embedding(
    object,
    embd = out,
    name = "rlsi",
    modality = modality,
    from = embedding
  )
}

#' Correlation of each LSI component with the library size
#'
#' @param object `SingleCellsAtac` or `SingleCellsMultiModal` class, after
#' [bixverse::calculate_lsi_sc()].
#'
#' @returns A data.table with `component` (1-based), `depth_cor` and
#' `dropped`.
#'
#' @export
get_lsi_depth_cor <- function(object) {
  info <- .require_atac(object)
  if (is.null(info$depth_cor)) {
    stop("Run calculate_lsi_sc() first.")
  }
  n <- length(info$depth_cor)

  data.table::data.table(
    component = seq_len(n),
    depth_cor = info$depth_cor,
    dropped = (seq_len(n) - 1L) %in% info$lsi_dropped
  )
}

## peak to gene ----------------------------------------------------------------

#' Link peaks to nearby genes by correlation over groups of cells
#'
#' @description
#' ArchR's `addPeak2GeneLinks()`. Peak and gene counts are summed per level of
#' `group_by` (fine clusters, or meta cell assignments, of the kept cells),
#' log-normalised per group, `log2(x / total * 1e4 + 1)`, and correlated over
#' the groups for every peak within `max_distance` of a gene's TSS. p-values
#' are the two-sided t-test of each correlation, FDR Benjamini-Hochberg over
#' all pairs. Unlike Signac's `LinkPeaks()` there is no GC-matched
#' background, which needs the genome sequence.
#'
#' @param object `SingleCellsMultiModal` class with RNA and ATAC.
#' @param gtf String. Gene annotation for the TSSs, matched to the RNA on
#' `gene_id` without its version.
#' @param group_by String. obs column whose levels are the groups; ArchR uses
#' a few hundred groups of ~100 cells.
#' @param max_distance Integer. Largest distance from a peak to a TSS.
#' @param biotypes Optional character vector. Gene biotypes to link. `NULL`
#' keeps all.
#' @param genes Optional character vector. Gene ids to restrict to; `NULL`
#' for every gene with a TSS.
#' @param min_cells Integer. Groups with fewer kept cells are left out.
#' @param fdr Optional numeric. Keep links below this FDR; `NULL` keeps every
#' tested pair.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @returns A data.table with `peak_id`, `gene_id`, `gene_name`, `distance`
#' (peak centre minus TSS), `r`, `p` and `fdr`.
#'
#' @export
#'
#' @references Granja, et al., Nat Genet, 2021
peak_to_gene_sc <- S7::new_generic(
  name = "peak_to_gene_sc",
  dispatch_args = "object",
  fun = function(
    object,
    gtf,
    group_by,
    max_distance = 250000L,
    biotypes = "protein_coding",
    genes = NULL,
    min_cells = 10L,
    fdr = 0.05,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method peak_to_gene_sc SingleCellsMultiModal
#'
#' @export
S7::method(peak_to_gene_sc, SingleCellsMultiModal) <- function(
  object,
  gtf,
  group_by,
  max_distance = 250000L,
  biotypes = "protein_coding",
  genes = NULL,
  min_cells = 10L,
  fdr = 0.05,
  .verbose = TRUE
) {
  # checks
  checkmate::assertFileExists(gtf)
  checkmate::qassert(group_by, "S1")
  checkmate::qassert(max_distance, "I1[1,)")
  checkmate::qassert(biotypes, c("0", "S+"))
  checkmate::qassert(genes, c("0", "S+"))
  checkmate::qassert(min_cells, "I1[1,)")
  checkmate::qassert(fdr, c("0", "N1(0,1]"))
  checkmate::qassert(.verbose, "B1")
  .require_atac(object)

  # function body
  codes <- .kept_obs_codes(object, group_by)
  kept <- get_cells_to_keep(object)
  groups <- purrr::map(seq_along(codes$levels) - 1L, \(g) {
    kept[codes$codes == g]
  })
  groups <- groups[lengths(groups) >= min_cells]
  if (length(groups) < 3L) {
    stop(sprintf(
      "%i groups of at least %i cells; need 3.",
      length(groups),
      min_cells
    ))
  }

  tss <- .atac_tss(gtf, biotypes)
  tss_dt <- data.table::data.table(
    gene_id = sub("\\.\\d+", "", tss$gene_id),
    gene_name = tss$gene_name,
    chrom = tss$chrom,
    pos = tss$pos
  )
  tss_dt <- tss_dt[!duplicated(tss_dt$gene_id)]
  rna_ids <- get_gene_names(object)
  tss_dt <- tss_dt[tss_dt$gene_id %in% rna_ids]
  if (!is.null(genes)) {
    tss_dt <- tss_dt[tss_dt$gene_id %in% genes]
  }

  rna <- rs_pseudobulk_cells_dense(
    f_path = get_rust_count_cell_f_path(object),
    cell_indices_ls = groups,
    assay = "raw",
    verbose = as.integer(.verbose)
  )[, match(tss_dt$gene_id, rna_ids), drop = FALSE]
  atac <- rs_pseudobulk_cells_dense(
    f_path = .atac_paths(object)$cells,
    cell_indices_ls = groups,
    assay = "raw",
    verbose = as.integer(.verbose)
  )
  var <- .atac_var(object, cols = c("feature_id", "chrom", "start", "end"))
  res <- rs_atac_peak_to_gene(
    peak_counts = atac,
    gene_counts = rna,
    chrom = var$chrom,
    start = as.integer(var$start),
    end = as.integer(var$end),
    tss_chrom = tss_dt$chrom,
    tss_pos = tss_dt$pos,
    max_distance = max_distance,
    scale = 1e4
  )

  out <- data.table::data.table(
    peak_id = var$feature_id[res$peak + 1L],
    gene_id = tss_dt$gene_id[res$gene + 1L],
    gene_name = tss_dt$gene_name[res$gene + 1L],
    distance = res$distance,
    r = res$r,
    p = res$p,
    fdr = res$fdr
  )
  if (!is.null(fdr)) {
    # computed outside `[`, where data.table would read `fdr` as the column
    keep <- out[["fdr"]] < fdr
    out <- out[keep]
  }

  out
}
