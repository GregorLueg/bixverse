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

## the h5ad fixture ------------------------------------------------------------

#' Write a synthetic spatial h5ad
#'
#' Enough of the AnnData layout for [bixverse::load_spatial_h5ad()] to swallow
#' it: a CSR `X`, an `obs` and `var` with `_index`, `obsm/spatial` and
#' optionally a `uns/spatial` library with scale factors and a `lowres` image.
#'
#' Two traps are baked in on purpose.
#'
#' The lattice is anisotropic, so `x` spans several times what `y` does and a
#' transposed read moves every spot rather than mirroring it. The image, when
#' asked for, gets a dark block exactly where the spots are under the correct
#' order, so a transposed read puts them on blank slide.
#'
#' `rhdf5::h5write()` reverses the dimensions of an R array on the way out
#' unless `native = TRUE`, which would hand the reader a 2 x N `obsm/spatial`.
#' Everything with a shape below therefore goes out native.
#'
#' @param h5_path String. File to write. Overwritten.
#' @param layout String. `"hex"` or `"square"`.
#' @param n_rows,n_cols Integer. Lattice shape.
#' @param n_genes Integer. Number of genes.
#' @param seed Integer. Seed for the barcodes and the counts.
#' @param orientation String. `"xy"` or `"yx"`. Column order written into
#' `obsm/spatial`.
#' @param obs_pixel_cols String. `"none"` for no `pxl_*_in_fullres` columns,
#' `"spaceranger"` for the 10x meaning of those names, `"scanpy"` for the
#' swapped meaning `scanpy.read_visium` gives them.
#' @param with_array_cols Boolean. Write `array_row` / `array_col`.
#' @param with_uns Boolean. Write a `uns/spatial` library at all.
#' @param with_image Boolean. Write a `lowres` image under that library.
#' @param image_scalef Float. `tissue_lowres_scalef` for that image.
#' @param in_tissue Integer vector or `NULL`. Per-spot `in_tissue`, in lattice
#' raster order. `NULL` writes all ones.
#' @param low_count_spots Integer vector or `NULL`. 1-based lattice positions
#' to strip down to a single count, so they fall out under any `min_lib_size`
#' above one. What makes the obs/`obsm` row mapping actually get tested.
#'
#' @return A list describing the fixture.
sp_make_spatial_h5ad <- function(
  h5_path,
  layout = c("hex", "square"),
  n_rows = 12L,
  n_cols = 12L,
  n_genes = 40L,
  seed = 11L,
  orientation = c("xy", "yx"),
  obs_pixel_cols = c("none", "spaceranger", "scanpy"),
  with_array_cols = TRUE,
  with_uns = TRUE,
  with_image = FALSE,
  image_scalef = 0.05,
  in_tissue = NULL,
  low_count_spots = NULL
) {
  layout <- match.arg(layout)
  orientation <- match.arg(orientation)
  obs_pixel_cols <- match.arg(obs_pixel_cols)

  grid <- sp_lattice_positions(layout, n_rows, n_cols)
  n_spots <- nrow(grid)

  # x runs from 200 to 200 + 50 * (last array_col); y from 200 in steps of 87.
  # Stretch x so the two ranges cannot be confused for one another.
  x <- as.numeric(grid$pxl_col_in_fullres) * 6
  y <- as.numeric(grid$pxl_row_in_fullres)

  set.seed(seed)
  barcodes <- vapply(
    seq_len(n_spots),
    \(i) {
      sprintf(
        "%s-1",
        paste(sample(c("A", "C", "G", "T"), 16L, TRUE), collapse = "")
      )
    },
    character(1)
  )
  stopifnot("Duplicated fixture barcodes" = !anyDuplicated(barcodes))

  set.seed(seed + 11L)
  counts <- sp_default_counts(grid, n_genes, n_rows)
  if (!is.null(low_count_spots)) {
    counts[, low_count_spots] <- 0L
    counts[3L, low_count_spots] <- 1L
  }
  gene_ids <- sprintf("ENSGSH%05i", seq_len(n_genes))

  # CSR over cells: row `i` is spot `i`, columns are genes.
  cells_by_genes <- t(counts)
  indptr <- c(0L, cumsum(apply(cells_by_genes, 1L, \(r) sum(r > 0))))
  nz <- lapply(seq_len(n_spots), \(i) which(cells_by_genes[i, ] > 0))
  indices <- as.integer(unlist(nz, use.names = FALSE) - 1L)
  data_vals <- as.numeric(unlist(
    lapply(seq_len(n_spots), \(i) cells_by_genes[i, nz[[i]]]),
    use.names = FALSE
  ))

  unlink(h5_path)
  rhdf5::h5createFile(h5_path)
  on.exit(tryCatch(rhdf5::h5closeAll(), error = function(e) invisible()))

  rhdf5::h5createGroup(h5_path, "X")
  rhdf5::h5write(data_vals, h5_path, "X/data")
  rhdf5::h5write(indices, h5_path, "X/indices")
  rhdf5::h5write(as.integer(indptr), h5_path, "X/indptr")

  rhdf5::h5createGroup(h5_path, "obs")
  rhdf5::h5write(barcodes, h5_path, "obs/_index")
  rhdf5::h5write(
    if (is.null(in_tissue)) rep(1L, n_spots) else as.integer(in_tissue),
    h5_path,
    "obs/in_tissue"
  )
  if (with_array_cols) {
    rhdf5::h5write(as.integer(grid$array_row), h5_path, "obs/array_row")
    rhdf5::h5write(as.integer(grid$array_col), h5_path, "obs/array_col")
  }
  if (obs_pixel_cols != "none") {
    # Space Ranger means row = y, column = x. `scanpy.read_visium` renames the
    # positions columns positionally and ends up meaning the opposite.
    under_row <- if (obs_pixel_cols == "spaceranger") y else x
    under_col <- if (obs_pixel_cols == "spaceranger") x else y
    rhdf5::h5write(under_row, h5_path, "obs/pxl_row_in_fullres")
    rhdf5::h5write(under_col, h5_path, "obs/pxl_col_in_fullres")
  }

  rhdf5::h5createGroup(h5_path, "var")
  rhdf5::h5write(gene_ids, h5_path, "var/_index")

  coords <- if (orientation == "xy") cbind(x, y) else cbind(y, x)
  rhdf5::h5createGroup(h5_path, "obsm")
  rhdf5::h5write(coords, h5_path, "obsm/spatial", native = TRUE)

  if (with_uns) {
    rhdf5::h5createGroup(h5_path, "uns")
    rhdf5::h5createGroup(h5_path, "uns/spatial")
    rhdf5::h5createGroup(h5_path, "uns/spatial/fixture_lib")
    rhdf5::h5createGroup(h5_path, "uns/spatial/fixture_lib/scalefactors")
    sf <- list(
      spot_diameter_fullres = 40,
      tissue_hires_scalef = 0.1,
      tissue_lowres_scalef = image_scalef,
      fiducial_diameter_fullres = 100
    )
    for (k in names(sf)) {
      rhdf5::h5write(
        sf[[k]],
        h5_path,
        sprintf("uns/spatial/fixture_lib/scalefactors/%s", k)
      )
    }
    if (with_image) {
      # Big enough for every spot under the right order and only just: the
      # transpose puts x where the height is and runs off the bottom.
      img_h <- as.integer(ceiling(max(y) * image_scalef)) + 6L
      img_w <- as.integer(ceiling(max(x) * image_scalef)) + 6L
      img <- array(1.0, dim = c(img_h, img_w, 3L))
      rows <- pmin(pmax(round(y * image_scalef), 0L), img_h - 1L) + 1L
      cols <- pmin(pmax(round(x * image_scalef), 0L), img_w - 1L) + 1L
      for (i in seq_len(n_spots)) {
        rr <- max(1L, rows[i] - 1L):min(img_h, rows[i] + 1L)
        cc <- max(1L, cols[i] - 1L):min(img_w, cols[i] + 1L)
        img[rr, cc, ] <- 0.2
      }
      rhdf5::h5createGroup(h5_path, "uns/spatial/fixture_lib/images")
      rhdf5::h5write(
        img,
        h5_path,
        "uns/spatial/fixture_lib/images/lowres",
        native = TRUE
      )
    }
  }

  rhdf5::h5closeAll()

  list(
    h5_path = h5_path,
    layout = layout,
    n_rows = n_rows,
    n_cols = n_cols,
    n_spots = n_spots,
    n_genes = n_genes,
    gene_ids = gene_ids,
    barcodes = barcodes,
    grid = grid,
    x = x,
    y = y,
    counts = counts,
    orientation = orientation,
    in_tissue = if (is.null(in_tissue)) rep(1L, n_spots) else in_tissue,
    low_count_spots = low_count_spots
  )
}

#' Ingest a fixture written by [sp_make_spatial_h5ad()]
#'
#' QC is switched off so the spot and gene counts stay what went in.
#'
#' @param fixture List. Output of [sp_make_spatial_h5ad()].
#' @param data_dir String. Where the on-disk store goes. Wiped first.
#' @param exp_id String. Experiment identifier.
#' @param sp_io_param List. Output of [bixverse::params_sp_h5ad_io()].
#' @param min_lib_size Integer. Library size QC floor. Above one this drops the
#' `low_count_spots` of the fixture, which is what exercises the mapping from
#' the loaded obs rows back onto the rows of `obsm/spatial`.
#'
#' @return A [bixverse::SpatialSpot()].
sp_load_h5ad_fixture <- function(
  fixture,
  data_dir,
  exp_id = "h5ad_a",
  sp_io_param = params_sp_h5ad_io(),
  min_lib_size = 0L
) {
  unlink(data_dir, recursive = TRUE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

  object <- SpatialSpot(dir_data = data_dir)
  suppressWarnings(load_spatial_h5ad(
    object,
    h5_path = fixture$h5_path,
    sp_io_param = sp_io_param,
    sc_qc_param = params_sc_min_quality(
      min_unique_genes = 0L,
      min_lib_size = min_lib_size,
      min_cells = 0L
    ),
    exp_id = exp_id,
    .verbose = FALSE
  ))
}
