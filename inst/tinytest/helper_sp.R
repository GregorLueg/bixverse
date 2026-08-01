# shared fixtures for the spatial tests ----------------------------------------

# tinytest has no helper mechanism of its own, but it does run every test file
# with the working directory set to the directory the file sits in. So the
# spatial test files pull this in with `source("helper_sp.R", local = TRUE)` and
# everything below lands in that file's own evaluation environment.
#
# Everything here is a function. Nothing is built at source time, so a file that
# only needs the capability probes does not pay for a Visium fixture it never
# looks at.

## capability probes -----------------------------------------------------------

#' Is the `png` package around to write fixture images?
sp_has_png <- function() {
  requireNamespace("png", quietly = TRUE)
}

#' Can the Rust image backend actually read an image?
#'
#' The image bindings need the `spatial-image` feature of `bixverse-rs`, which
#' in turn needs OpenSlide at build and run time. A build without it either has
#' no `rs_sp_image_*` at all or fails on the first call, so the probe is a real
#' call on a throwaway PNG rather than a version check.
sp_image_backend_ok <- function() {
  if (!sp_has_png()) {
    return(FALSE)
  }
  probe <- file.path(tempdir(), "bixverse_sp_image_probe.png")
  ok <- tryCatch(
    {
      if (!file.exists(probe)) {
        png::writePNG(array(0.5, dim = c(8L, 8L, 3L)), target = probe)
      }
      meta <- rs_sp_image_metadata(probe)
      is.list(meta) && isTRUE(meta$width == 8L)
    },
    error = function(e) FALSE
  )

  isTRUE(ok)
}

#' Path to the local Space Ranger example data, or `NULL`
#'
#' Not shipped with the package. The tests that use it are the ones pinning
#' against Space Ranger's own numbers, and they skip when it is not there.
sp_example_visium_dir <- function() {
  path <- path.expand("~/Desktop/visium_example_data")
  if (
    dir.exists(path) &&
      file.exists(file.path(path, "spatial", "tissue_positions.csv"))
  ) {
    path
  } else {
    NULL
  }
}

## lattice geometry ------------------------------------------------------------

#' Lay out a Visium hex or a Visium HD square lattice
#'
#' Raster order, i.e. row by row. On the honeycomb `array_col` steps by two
#' inside a row and odd rows are offset by one, so `array_row %% 2 ==
#' array_col %% 2` holds everywhere and an interior spot has exactly six
#' neighbours. Getting that parity wrong drops it to four.
#'
#' @param layout String. `"hex"` or `"square"`.
#' @param n_rows Integer. Number of lattice rows.
#' @param n_cols Integer. Number of spots per row.
#'
#' @return A `data.table` with `array_row`, `array_col`, `pxl_row_in_fullres`
#' and `pxl_col_in_fullres`.
sp_lattice_positions <- function(layout = c("hex", "square"), n_rows, n_cols) {
  layout <- match.arg(layout)

  grid <- data.table::rbindlist(lapply(seq_len(n_rows) - 1L, \(r) {
    array_col <- if (layout == "hex") {
      seq(r %% 2L, by = 2L, length.out = n_cols)
    } else {
      seq_len(n_cols) - 1L
    }
    data.table::data.table(array_row = r, array_col = array_col)
  }))

  # pixel pitch is deliberately anisotropic on the honeycomb so that a distance
  # based layout cannot accidentally reproduce the lattice answer
  grid[, pxl_col_in_fullres := as.integer(200L + array_col * 50L)]
  grid[,
    pxl_row_in_fullres := as.integer(
      200L + array_row * if (layout == "hex") 87L else 50L
    )
  ]

  return(grid)
}

#' Reference neighbour lists for a lattice, in plain R
#'
#' Independent of the Rust builder: an exact table lookup over the array
#' coordinates, which is what the lattice layouts are supposed to be.
#'
#' @param array_row,array_col Integer vectors. Array coordinates, one per spot,
#' in the order the graph is indexed.
#' @param layout String. `"hex"` or `"square"`.
#' @param connectivity Integer. `4L` or `8L`. Square lattice only.
#'
#' @return List of integer vectors, 0-based local neighbour indices, ascending.
sp_lattice_neighbours <- function(
  array_row,
  array_col,
  layout = c("hex", "square"),
  connectivity = 4L
) {
  layout <- match.arg(layout)

  offsets <- if (layout == "hex") {
    cbind(
      c(-1L, -1L, 0L, 0L, 1L, 1L),
      c(-1L, 1L, -2L, 2L, -1L, 1L)
    )
  } else if (connectivity == 4L) {
    cbind(c(-1L, 0L, 0L, 1L), c(0L, -1L, 1L, 0L))
  } else {
    cbind(
      c(-1L, -1L, -1L, 0L, 0L, 1L, 1L, 1L),
      c(-1L, 0L, 1L, -1L, 1L, -1L, 0L, 1L)
    )
  }

  key <- sprintf("%i_%i", array_row, array_col)
  lookup <- stats::setNames(seq_along(key) - 1L, key)

  lapply(seq_along(array_row), \(i) {
    wanted <- sprintf(
      "%i_%i",
      array_row[i] + offsets[, 1L],
      array_col[i] + offsets[, 2L]
    )
    hits <- lookup[wanted]
    sort(as.integer(hits[!is.na(hits)]))
  })
}

## the visium fixture ----------------------------------------------------------

#' Default synthetic counts for the fixture
#'
#' Poisson noise everywhere, a north/south gradient planted on gene 1 and a
#' flat but high-mean gene 2 as a negative control.
#'
#' Both planted genes are deliberately modest. Moran's I runs on normalised
#' counts by default, so a gene big enough to move the library size drags its
#' own spatial structure into every other gene through the denominator: WP6's
#' first fixture had a pure-noise gene scoring 0.19 for exactly that reason.
#'
#' @param positions `data.table`. Lattice rows in the order the counts columns
#' are in.
#' @param n_genes Integer. Number of genes.
#' @param n_rows Integer. Number of lattice rows, used to split north from
#' south.
#'
#' @return An integer matrix, genes x spots.
sp_default_counts <- function(positions, n_genes, n_rows) {
  n_spots <- nrow(positions)
  north <- positions$array_row < n_rows / 2L

  counts <- matrix(
    stats::rpois(n_genes * n_spots, lambda = 4),
    nrow = n_genes,
    ncol = n_spots
  )
  counts[1L, ] <- ifelse(north, 30L, 3L)
  counts[2L, ] <- stats::rpois(n_spots, lambda = 20)

  return(counts)
}

#' Write a synthetic Space Ranger output directory
#'
#' Enough of the layout for [bixverse::load_visium()] to swallow it: an mtx
#' triplet directory, `tissue_positions.csv`, `scalefactors_json.json` and a
#' hires PNG when `png` is installed.
#'
#' The barcode ordering trap is baked in on purpose. `barcodes.tsv.gz` is
#' lexicographic while `tissue_positions.csv` is in array raster order, so
#' anything that joins the two by position rather than by barcode scrambles the
#' coordinates.
#'
#' @param visium_dir String. Directory to write into. Created if missing.
#' @param layout String. `"hex"` or `"square"`.
#' @param n_rows,n_cols Integer. Lattice shape.
#' @param n_genes Integer. Number of genes.
#' @param seed Integer. Seed for the barcodes and the counts.
#' @param counts_fun Function or `NULL`. Takes `(positions, n_genes, n_rows)`
#' in barcode order and returns a genes x spots count matrix.
#' @param image_size Integer. Side of the square hires PNG.
#'
#' @return A list describing the fixture.
sp_make_visium <- function(
  visium_dir,
  layout = c("hex", "square"),
  n_rows = 12L,
  n_cols = 12L,
  n_genes = 60L,
  seed = 11L,
  counts_fun = NULL,
  image_size = 200L
) {
  layout <- match.arg(layout)
  if (is.null(counts_fun)) {
    counts_fun <- sp_default_counts
  }

  matrix_dir <- file.path(visium_dir, "raw_feature_bc_matrix")
  spatial_dir <- file.path(visium_dir, "spatial")
  dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)

  grid <- sp_lattice_positions(layout, n_rows, n_cols)
  n_spots <- nrow(grid)

  set.seed(seed)
  barcodes_raster <- vapply(
    seq_len(n_spots),
    \(i) {
      sprintf(
        "%s-1",
        paste(sample(c("A", "C", "G", "T"), 16L, TRUE), collapse = "")
      )
    },
    character(1)
  )
  stopifnot("Duplicated fixture barcodes" = !anyDuplicated(barcodes_raster))
  barcodes_matrix <- sort(barcodes_raster)

  positions <- data.table::data.table(
    barcode = barcodes_raster,
    in_tissue = 1L,
    array_row = grid$array_row,
    array_col = grid$array_col,
    pxl_row_in_fullres = grid$pxl_row_in_fullres,
    pxl_col_in_fullres = grid$pxl_col_in_fullres
  )
  data.table::setkey(positions, barcode)

  pos_matrix_order <- positions[barcodes_matrix]

  set.seed(seed + 11L)
  counts <- counts_fun(pos_matrix_order, n_genes, n_rows)
  gene_ids <- sprintf("ENSGSP%05i", seq_len(n_genes))

  triplets <- which(counts > 0, arr.ind = TRUE)
  con <- gzfile(file.path(matrix_dir, "matrix.mtx.gz"), "w")
  writeLines("%%MatrixMarket matrix coordinate integer general", con)
  writeLines(
    sprintf("%d %d %d", nrow(counts), ncol(counts), nrow(triplets)),
    con
  )
  writeLines(
    sprintf(
      "%d %d %d",
      triplets[, "row"],
      triplets[, "col"],
      as.integer(counts[triplets])
    ),
    con
  )
  close(con)

  data.table::fwrite(
    data.table::data.table(barcodes_matrix),
    file.path(matrix_dir, "barcodes.tsv.gz"),
    col.names = FALSE
  )
  data.table::fwrite(
    data.table::data.table(
      gene_ids,
      sprintf("SYM%03i", seq_along(gene_ids)),
      "Gene Expression"
    ),
    file.path(matrix_dir, "features.tsv.gz"),
    sep = "\t",
    col.names = FALSE
  )
  data.table::fwrite(
    positions,
    file.path(spatial_dir, "tissue_positions.csv"),
    col.names = TRUE
  )

  scale_factors <- list(
    regist_target_img_scalef = 0.2,
    tissue_hires_scalef = 0.1,
    tissue_lowres_scalef = 0.05,
    fiducial_diameter_fullres = 100,
    spot_diameter_fullres = 40
  )
  writeLines(
    jsonlite::toJSON(scale_factors, auto_unbox = TRUE),
    file.path(spatial_dir, "scalefactors_json.json")
  )

  # hires and lowres both, at the ratio the scale factors imply, so a test can
  # ask whether the resolution actually changed the measurement
  png_written <- FALSE
  if (sp_has_png()) {
    set.seed(seed + 22L)
    lowres_size <- max(8L, as.integer(image_size / 2L))
    png::writePNG(
      array(
        stats::runif(image_size * image_size * 3L),
        dim = c(image_size, image_size, 3L)
      ),
      target = file.path(spatial_dir, "tissue_hires_image.png")
    )
    png::writePNG(
      array(
        stats::runif(lowres_size * lowres_size * 3L),
        dim = c(lowres_size, lowres_size, 3L)
      ),
      target = file.path(spatial_dir, "tissue_lowres_image.png")
    )
    png_written <- TRUE
  }

  list(
    visium_dir = visium_dir,
    layout = layout,
    n_rows = n_rows,
    n_cols = n_cols,
    n_spots = n_spots,
    n_genes = n_genes,
    gene_ids = gene_ids,
    barcodes = barcodes_matrix,
    positions = positions,
    pos_matrix_order = pos_matrix_order,
    counts = counts,
    scale_factors = scale_factors,
    image_size = image_size,
    png_written = png_written
  )
}

#' Ingest a fixture written by [sp_make_visium()]
#'
#' QC is switched off so the gene and spot counts stay exactly what went in,
#' which every downstream shape assertion relies on.
#'
#' @param fixture List. Output of [sp_make_visium()].
#' @param data_dir String. Where the on-disk store goes. Wiped first.
#' @param exp_id String. Experiment identifier.
#'
#' @return A [bixverse::SpatialSpot()].
sp_load_visium_fixture <- function(fixture, data_dir, exp_id = "sample_a") {
  unlink(data_dir, recursive = TRUE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  object <- SpatialSpot(dir_data = data_dir)
  load_visium(
    object,
    visium_dir = fixture$visium_dir,
    sp_io_param = params_sp_visium_io(),
    sc_qc_param = params_sc_min_quality(
      min_unique_genes = 0L,
      min_lib_size = 0L,
      min_cells = 0L
    ),
    exp_id = exp_id,
    .verbose = FALSE
  )
}

#' Obs rows of one experiment in graph order
#'
#' Ascending `cell_idx`, which is the order
#' [bixverse::get_spatial_coords()] and the neighbour lists are in.
#'
#' @param object A [bixverse::SpatialSpot()] or subset.
#' @param exp_id String. Experiment identifier.
#'
#' @return A `data.table`.
sp_obs_in_graph_order <- function(object, exp_id) {
  obs <- get_sc_obs(object, filtered = TRUE)
  # mask outside the `[` call: inside `i` the `exp_id` column would shadow the
  # `exp_id` argument and the filter would keep every sample
  keep <- obs[["exp_id"]] == exp_id
  obs <- obs[keep, ]
  data.table::setorderv(obs, "cell_idx")

  return(obs)
}
