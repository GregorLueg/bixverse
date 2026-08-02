# spatial h5ad ingest ----------------------------------------------------------

# The orientation assertions are the ones that matter. The fixture lattice is
# anisotropic on purpose, so a transposed read is a different point cloud rather
# than a mirror image, and every check below compares against the planted
# coordinates rather than against a shape.

source("helper_sp.R", local = TRUE)

library(bixverse)

h5_dir <- file.path(tempdir(), "bixverse_sp_h5ad_tests")
dir.create(h5_dir, recursive = TRUE, showWarnings = FALSE)

## the reader --------------------------------------------------------------

### orientation ------------------------------------------------------------

# Same tissue, written both ways round, with the obs pixel columns in the 10x
# meaning. The reader has to land on the same coordinates either way.
fx_xy <- sp_make_spatial_h5ad(
  file.path(h5_dir, "orient_xy.h5ad"),
  orientation = "xy",
  obs_pixel_cols = "spaceranger"
)
fx_yx <- sp_make_spatial_h5ad(
  file.path(h5_dir, "orient_yx.h5ad"),
  orientation = "yx",
  obs_pixel_cols = "spaceranger"
)

read_xy <- rs_sp_read_h5ad_spatial(fx_xy$h5_path, NULL, "xy")
read_yx <- rs_sp_read_h5ad_spatial(fx_yx$h5_path, NULL, "xy")

expect_equal(
  target = "obs_pixel_columns",
  current = read_xy$evidence,
  info = "the 10x-named pixel columns settle the order"
)
expect_equal(
  target = "xy",
  current = read_xy$orientation,
  info = "column 0 is x in the xy fixture"
)
expect_equal(
  target = "yx",
  current = read_yx$orientation,
  info = "column 0 is y in the yx fixture"
)

expect_equal(
  target = read_xy$coords,
  current = read_yx$coords,
  info = "both column orders resolve to the same (x, y) pairs"
)
expect_equal(
  target = fx_xy$x,
  current = as.numeric(read_xy$coords[, 1L]),
  info = "column 1 of the result is the planted x"
)
expect_equal(
  target = fx_xy$y,
  current = as.numeric(read_xy$coords[, 2L]),
  info = "column 2 of the result is the planted y"
)

# The mutant: an inverted swap. If `to_xy` swapped for `Xy` instead of `Yx`,
# the read of the yx fixture would come back with the two columns exchanged,
# so the x range would collapse onto the y range.
expect_true(
  diff(range(read_yx$coords[, 1L])) > 4 * diff(range(read_yx$coords[, 2L])),
  info = "x spans several times what y does, so a swap is detectable"
)

# scanpy renames the tissue positions columns positionally, so a file that kept
# them carries labels meaning the opposite of what 10x means. With no image
# there is nothing better to go on and the reader takes them, which is the
# failure the image evidence exists to prevent.
fx_scanpy <- sp_make_spatial_h5ad(
  file.path(h5_dir, "orient_scanpy.h5ad"),
  orientation = "xy",
  obs_pixel_cols = "scanpy"
)
read_scanpy <- rs_sp_read_h5ad_spatial(fx_scanpy$h5_path, NULL, "xy")
expect_equal(
  target = "yx",
  current = read_scanpy$orientation,
  info = "the swapped labels are believed when nothing else is available"
)

# ... and with the image present the tissue overrules them.
fx_scanpy_img <- sp_make_spatial_h5ad(
  file.path(h5_dir, "orient_scanpy_img.h5ad"),
  orientation = "xy",
  obs_pixel_cols = "scanpy",
  with_image = TRUE
)
read_scanpy_img <- rs_sp_read_h5ad_spatial(fx_scanpy_img$h5_path, NULL, "xy")
expect_true(
  read_scanpy_img$evidence %in% c("image_tissue", "image_frame"),
  info = "the image is consulted before the obs labels"
)
expect_equal(
  target = "xy",
  current = read_scanpy_img$orientation,
  info = "the image gets the order right where the labels do not"
)
expect_equal(
  target = fx_scanpy_img$x,
  current = as.numeric(read_scanpy_img$coords[, 1L]),
  info = "the image-resolved coordinates are the planted ones"
)

# Nothing to go on at all: the parameter is taken and reported as such.
fx_bare <- sp_make_spatial_h5ad(
  file.path(h5_dir, "orient_bare.h5ad"),
  orientation = "xy",
  obs_pixel_cols = "none",
  with_uns = FALSE
)
read_bare <- rs_sp_read_h5ad_spatial(fx_bare$h5_path, NULL, "xy")
expect_equal(
  target = "assumed",
  current = read_bare$evidence,
  info = "the reader says when it had nothing to go on"
)
expect_equal(
  target = "yx",
  current = rs_sp_read_h5ad_spatial(fx_bare$h5_path, NULL, "yx")$orientation,
  info = "the fallback parameter is honoured rather than ignored"
)

### what else is in the file -----------------------------------------------

expect_true(
  read_xy$has_array_indices,
  info = "the fixture carries the lattice indices"
)
expect_equal(
  target = "fixture_lib",
  current = read_xy$library_id,
  info = "the only uns/spatial library is picked without being named"
)
expect_equal(
  target = 4L,
  current = length(read_xy$scale_factor_names),
  info = "every scale factor is passed through"
)
expect_true(
  is.na(read_bare$library_id),
  info = "a file with no uns/spatial reports no library"
)
expect_equal(
  target = 0L,
  current = length(read_bare$scale_factor_names),
  info = "no uns/spatial means no scale factors"
)

fx_no_array <- sp_make_spatial_h5ad(
  file.path(h5_dir, "no_array.h5ad"),
  with_array_cols = FALSE,
  obs_pixel_cols = "spaceranger"
)
expect_false(
  rs_sp_read_h5ad_spatial(fx_no_array$h5_path, NULL, "xy")$has_array_indices,
  info = "a file without array_row/array_col reports it"
)

expect_error(
  rs_sp_read_h5ad_spatial(fx_xy$h5_path, "not_a_library", "xy"),
  info = "a library that is not there is an error, not a silent fallback"
)
expect_error(
  rs_sp_read_h5ad_spatial(fx_xy$h5_path, NULL, "rowcol"),
  info = "a nonsense orientation is rejected"
)

## the layout report -------------------------------------------------------

expect_equal(
  target = c("hex", "square", "knn", "radius"),
  current = bixverse:::.sp_available_layouts(TRUE),
  info = "every layout is available with the array indices"
)
expect_equal(
  target = c("knn", "radius"),
  current = bixverse:::.sp_available_layouts(FALSE),
  info = "the lattice layouts drop out without the array indices"
)

## the join back onto the loaded cells -------------------------------------

# QC drops cells, so the obs table is a subset of the file in file order while
# `obsm/spatial` still holds every row. Zipping by position shifts every
# coordinate from the first dropped cell onwards.
expect_equal(
  target = c(1L, 3L, 4L),
  current = bixverse:::.sp_h5ad_coord_rows(
    obs_ids = c("a", "c", "d"),
    all_ids = c("a", "b", "c", "d"),
    n_coords = 4L
  ),
  info = "the surviving cells map onto their own rows, not onto 1:n"
)
expect_error(
  bixverse:::.sp_h5ad_coord_rows(c("a"), c("a", "a"), 2L),
  info = "duplicated obs identifiers are refused"
)
expect_error(
  bixverse:::.sp_h5ad_coord_rows(c("z"), c("a", "b"), 2L),
  info = "a cell with no row in the index is an error"
)
expect_error(
  bixverse:::.sp_h5ad_coord_rows(c("a"), c("a", "b"), 3L),
  info = "an obs/obsm length mismatch is refused"
)

## end to end --------------------------------------------------------------

data_dir <- file.path(h5_dir, "store")
obj <- sp_load_h5ad_fixture(fx_xy, data_dir, exp_id = "h5ad_a")

expect_true(
  S7::S7_inherits(obj, SpatialSpot),
  info = "the h5ad loads into a SpatialSpot"
)
expect_equal(
  target = fx_xy$n_spots,
  current = as.integer(S7::prop(obj, "dims")[1L]),
  info = "every spot made it in"
)
expect_equal(
  target = "h5ad_a",
  current = get_sample_ids(obj),
  info = "the sample is registered"
)

coords <- get_spatial_coords(obj, exp_id = "h5ad_a", filtered = TRUE)
obs_order <- sp_obs_in_graph_order(obj, "h5ad_a")
planted <- match(obs_order[["cell_id"]], fx_xy$barcodes)

expect_equal(
  target = fx_xy$x[planted],
  current = as.numeric(coords[, "x"]),
  info = "x survives the round trip through the DuckDB in the right order"
)
expect_equal(
  target = fx_xy$y[planted],
  current = as.numeric(coords[, "y"]),
  info = "y survives the round trip through the DuckDB in the right order"
)

# The same tissue written the other way round has to give the same object.
obj_yx <- sp_load_h5ad_fixture(
  fx_yx,
  file.path(h5_dir, "store_yx"),
  exp_id = "h5ad_a"
)
coords_yx <- get_spatial_coords(obj_yx, exp_id = "h5ad_a", filtered = TRUE)
expect_equal(
  target = coords,
  current = coords_yx,
  info = "the ingest is invariant to the on-disk column order"
)

### the QC drop ------------------------------------------------------------

# `load_h5ad()` drops cells that fail QC, so the obs table becomes a subset of
# the file while `obsm/spatial` still has every row. Anything that zips the two
# by position shifts every coordinate from the first dropped spot onwards.
dropped <- c(2L, 5L, 9L, 40L, 77L)
fx_qc <- sp_make_spatial_h5ad(
  file.path(h5_dir, "qc_drop.h5ad"),
  obs_pixel_cols = "spaceranger",
  low_count_spots = dropped
)
obj_qc <- sp_load_h5ad_fixture(
  fx_qc,
  file.path(h5_dir, "store_qc"),
  exp_id = "qc",
  min_lib_size = 10L
)

expect_equal(
  target = fx_qc$n_spots - length(dropped),
  current = as.integer(S7::prop(obj_qc, "dims")[1L]),
  info = "the low-count spots are dropped by QC"
)

coords_qc <- get_spatial_coords(obj_qc, exp_id = "qc", filtered = TRUE)
obs_qc <- sp_obs_in_graph_order(obj_qc, "qc")
planted_qc <- match(obs_qc[["cell_id"]], fx_qc$barcodes)

expect_equal(
  target = fx_qc$x[planted_qc],
  current = as.numeric(coords_qc[, "x"]),
  info = "surviving spots keep their own x, not their neighbour's"
)
expect_equal(
  target = fx_qc$y[planted_qc],
  current = as.numeric(coords_qc[, "y"]),
  info = "surviving spots keep their own y, not their neighbour's"
)
expect_false(
  any(dropped %in% planted_qc),
  info = "the dropped spots really are gone, so the shift is reachable"
)

### the analysis chain -----------------------------------------------------

obj <- build_spatial_graph_sp(
  obj,
  exp_id = "h5ad_a",
  graph_params = params_sp_graph(layout = "knn", k = 6L),
  .verbose = FALSE
)
graph <- get_per_sample_spatial_graph(obj, exp_id = "h5ad_a")
expect_equal(
  target = fx_xy$n_spots,
  current = nrow(graph),
  info = "the kNN graph covers every spot"
)

obj <- build_spatial_graph_sp(
  obj,
  exp_id = "h5ad_a",
  graph_params = params_sp_graph(layout = "hex"),
  .verbose = FALSE
)
hex_graph <- get_per_sample_spatial_graph(obj, exp_id = "h5ad_a")
ref_neighbours <- sp_lattice_neighbours(
  array_row = as.integer(obs_order[["array_row"]]),
  array_col = as.integer(obs_order[["array_col"]]),
  layout = "hex"
)
expect_equal(
  target = sum(lengths(ref_neighbours)),
  current = length(hex_graph@x),
  info = "the hex layout reproduces the reference lattice off the h5ad obs"
)

obj <- morans_i_sp(obj, exp_id = "h5ad_a", .verbose = FALSE)
morans <- get_per_sample_morans_i(obj, exp_id = "h5ad_a")
expect_true(
  data.table::is.data.table(morans),
  info = "Moran's I runs on an h5ad-loaded sample"
)
expect_equal(
  target = fx_xy$n_genes,
  current = nrow(morans),
  info = "Moran's I covers every gene"
)
# Gene 1 carries the planted north/south gradient, gene 2 is flat noise.
gradient <- morans[morans$gene_id == fx_xy$gene_ids[1L], ][["morans_i"]]
flat <- morans[morans$gene_id == fx_xy$gene_ids[2L], ][["morans_i"]]
expect_true(
  gradient > flat,
  info = "the planted spatial gene beats the flat one"
)

obj <- sparkx_sp(obj, exp_id = "h5ad_a", .verbose = FALSE)
sparkx <- get_per_sample_sparkx(obj, exp_id = "h5ad_a")
expect_equal(
  target = fx_xy$n_genes,
  current = nrow(sparkx),
  info = "SPARK-X runs on an h5ad-loaded sample"
)

### the orientation is an isometry -----------------------------------------

# The claim the whole evidence ranking rests on: swapping x and y is a
# reflection, so every pairwise distance survives it and no distance-based
# statistic can tell the two apart. `fx_bare` settles nothing on its own, so
# forcing the two column orders gives a sample and its own transpose.
obj_a <- sp_load_h5ad_fixture(
  fx_bare,
  file.path(h5_dir, "store_iso_a"),
  exp_id = "iso",
  sp_io_param = params_sp_h5ad_io(assume_orientation = "xy")
)
obj_b <- sp_load_h5ad_fixture(
  fx_bare,
  file.path(h5_dir, "store_iso_b"),
  exp_id = "iso",
  sp_io_param = params_sp_h5ad_io(assume_orientation = "yx")
)

coords_a <- get_spatial_coords(obj_a, exp_id = "iso", filtered = TRUE)
coords_b <- get_spatial_coords(obj_b, exp_id = "iso", filtered = TRUE)
expect_equal(
  target = as.numeric(coords_a[, "x"]),
  current = as.numeric(coords_b[, "y"]),
  info = "the forced orders really are each other's transpose"
)

graph_params <- params_sp_graph(layout = "knn", k = 6L)
obj_a <- build_spatial_graph_sp(
  obj_a,
  exp_id = "iso",
  graph_params = graph_params,
  .verbose = FALSE
)
obj_b <- build_spatial_graph_sp(
  obj_b,
  exp_id = "iso",
  graph_params = graph_params,
  .verbose = FALSE
)
expect_equal(
  target = get_per_sample_spatial_graph(obj_a, exp_id = "iso"),
  current = get_per_sample_spatial_graph(obj_b, exp_id = "iso"),
  info = "the kNN graph is blind to a transpose"
)

obj_a <- morans_i_sp(obj_a, exp_id = "iso", .verbose = FALSE)
obj_b <- morans_i_sp(obj_b, exp_id = "iso", .verbose = FALSE)
expect_equal(
  target = get_per_sample_morans_i(obj_a, exp_id = "iso")[["morans_i"]],
  current = get_per_sample_morans_i(obj_b, exp_id = "iso")[["morans_i"]],
  info = "Moran's I is blind to a transpose"
)

obj_a <- sparkx_sp(obj_a, exp_id = "iso", .verbose = FALSE)
obj_b <- sparkx_sp(obj_b, exp_id = "iso", .verbose = FALSE)
expect_equal(
  target = get_per_sample_sparkx(obj_a, exp_id = "iso")[["pval"]],
  current = get_per_sample_sparkx(obj_b, exp_id = "iso")[["pval"]],
  info = "SPARK-X is blind to a transpose"
)

### the coordinates-only sample --------------------------------------------

fx_coords_only <- sp_make_spatial_h5ad(
  file.path(h5_dir, "coords_only.h5ad"),
  obs_pixel_cols = "none",
  with_uns = FALSE
)
obj_bare <- sp_load_h5ad_fixture(
  fx_coords_only,
  file.path(h5_dir, "store_bare"),
  exp_id = "bare"
)

expect_equal(
  target = "bare",
  current = get_sample_ids(obj_bare),
  info = "a sample with no scale factors at all still registers"
)
expect_equal(
  target = list(),
  current = get_sample(obj_bare, "bare")$scale_factors,
  info = "and carries an empty scale factor list"
)

obj_bare <- build_spatial_graph_sp(
  obj_bare,
  exp_id = "bare",
  graph_params = params_sp_graph(layout = "knn", k = 6L),
  .verbose = FALSE
)
obj_bare <- morans_i_sp(obj_bare, exp_id = "bare", .verbose = FALSE)
expect_equal(
  target = fx_coords_only$n_genes,
  current = nrow(get_per_sample_morans_i(obj_bare, exp_id = "bare")),
  info = "Moran's I runs on a coordinates-only sample"
)

img_err <- tryCatch(
  image_features_sp(obj_bare, exp_id = "bare", .verbose = FALSE),
  error = function(e) conditionMessage(e)
)
expect_true(
  is.character(img_err) && grepl("image", img_err, ignore.case = TRUE),
  info = "image features fail on a coordinates-only sample, naming the image"
)

### the warnings -----------------------------------------------------------

fx_no_lattice <- sp_make_spatial_h5ad(
  file.path(h5_dir, "no_lattice.h5ad"),
  with_array_cols = FALSE,
  obs_pixel_cols = "spaceranger"
)
no_lattice_dir <- file.path(h5_dir, "store_no_lattice")
dir.create(no_lattice_dir, recursive = TRUE, showWarnings = FALSE)
expect_warning(
  load_spatial_h5ad(
    SpatialSpot(dir_data = no_lattice_dir),
    h5_path = fx_no_lattice$h5_path,
    sc_qc_param = params_sc_min_quality(
      min_unique_genes = 0L,
      min_lib_size = 0L,
      min_cells = 0L
    ),
    exp_id = "no_lattice",
    .verbose = FALSE
  ),
  pattern = "hex",
  info = "the missing lattice layouts are flagged at load time"
)

### the file with no coordinates -------------------------------------------

no_coords <- file.path(h5_dir, "no_coords.h5ad")
file.copy(fx_xy$h5_path, no_coords, overwrite = TRUE)
rhdf5::h5delete(no_coords, "obsm/spatial")
rhdf5::h5closeAll()

no_coords_err <- tryCatch(
  rs_sp_read_h5ad_spatial(no_coords, NULL, "xy"),
  error = function(e) conditionMessage(e)
)
expect_true(
  grepl("load_visium", no_coords_err),
  info = "a file with no coordinates points at the Space Ranger loader"
)

## the relaxed scale factor validation -------------------------------------

expect_true(
  isTRUE(validate_scale_factors(list(), technology = "visium")),
  info = "an empty scale factor list is legal for a visium sample"
)
expect_true(
  isTRUE(validate_scale_factors(
    list(spot_diameter_fullres = 40),
    technology = "visium"
  )),
  info = "a partial scale factor list is legal"
)
expect_error(
  validate_scale_factors(
    list(spot_diameter_fullres = c(1, 2)),
    technology = "visium"
  ),
  info = "but a scale factor that is there still has to be a single number"
)
expect_error(
  validate_scale_factors(
    list(tissue_hires_scalef = NA_real_),
    technology = "visium"
  ),
  info = "and it still has to be finite"
)
expect_true(
  isTRUE(validate_scale_factors(
    list(tissue_hires_trim_scalef = 0.2),
    technology = "visium"
  )),
  info = "an unknown scale factor key is carried rather than rejected"
)

sample_no_sf <- new_spatial_sample(
  exp_id = "image_free",
  technology = "visium",
  n_spots = 10L,
  scale_factors = list()
)
expect_true(
  isTRUE(validate_spatial_sample(sample_no_sf)),
  info = "a visium sample with no scale factors is a valid sample"
)

## the params wrapper ------------------------------------------------------

expect_equal(
  target = list(
    library_id = NULL,
    assume_orientation = "xy",
    in_tissue_only = TRUE,
    technology = "visium"
  ),
  current = params_sp_h5ad_io(),
  info = "the defaults are the scanpy column order and in-tissue only"
)
expect_error(
  params_sp_h5ad_io(assume_orientation = "rowcol"),
  info = "an unknown column order is rejected"
)
expect_true(
  isTRUE(bixverse:::checkSpH5adIo(params_sp_h5ad_io())),
  info = "the checkmate extension accepts what the wrapper builds"
)
expect_true(
  is.character(bixverse:::checkSpH5adIo(list(library_id = NULL))),
  info = "and rejects a list that is missing fields"
)
expect_true(
  is.character(bixverse:::checkSpH5adIo(
    utils::modifyList(params_sp_h5ad_io(), list(assume_orientation = "nope"))
  )),
  info = "and rejects an out-of-set column order"
)

unlink(h5_dir, recursive = TRUE)
