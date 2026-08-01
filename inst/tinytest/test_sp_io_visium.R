# visium io --------------------------------------------------------------------

library(magrittr)

test_temp_dir <- file.path(tempdir(), "io_visium")
dir.create(test_temp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir))

## fixture ---------------------------------------------------------------------

### parameters -----------------------------------------------------------------

n_array_rows <- 6L
n_array_cols <- 8L
n_spots <- n_array_rows * n_array_cols
n_genes <- 30L

# scale factors chosen so the derived CytAssist scale is a round number:
# hires 200 px / 0.05 -> a 4000 px full-res frame, registration target
# 4000 * 0.2 = 800 px, shipped CytAssist file 400 px, so 0.2 * 0.5 = 0.1.
hires_dims <- c(200L, 165L)
lowres_dims <- c(60L, 50L)
cytassist_dims <- c(400L, 400L)
sf_json <- list(
  regist_target_img_scalef = 0.2,
  tissue_hires_scalef = 0.05,
  tissue_lowres_scalef = 0.015,
  fiducial_diameter_fullres = 100,
  spot_diameter_fullres = 60
)
expected_cytassist_scalef <- 0.1

qc_permissive <- params_sc_min_quality(
  min_unique_genes = 0L,
  min_lib_size = 0L,
  min_cells = 0L
)

### helpers --------------------------------------------------------------------

# A header-only baseline TIFF. Enough for the dimension parser, deliberately
# not enough for anything to read it as an image.
write_stub_tiff <- function(path, width, height) {
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw("II"), con)
  writeBin(42L, con, size = 2L, endian = "little")
  writeBin(8L, con, size = 4L, endian = "little")
  writeBin(2L, con, size = 2L, endian = "little")
  for (entry in list(c(256L, width), c(257L, height))) {
    writeBin(entry[1], con, size = 2L, endian = "little")
    writeBin(4L, con, size = 2L, endian = "little")
    writeBin(1L, con, size = 4L, endian = "little")
    writeBin(entry[2], con, size = 4L, endian = "little")
  }
  writeBin(0L, con, size = 4L, endian = "little")
}

# Real PNG if the package is around, otherwise a signature plus IHDR stub that
# still carries the dimensions.
write_test_png <- function(path, width, height) {
  if (requireNamespace("png", quietly = TRUE)) {
    png::writePNG(array(0.5, dim = c(height, width, 3L)), target = path)
    return(invisible(TRUE))
  }
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(137, 80, 78, 71, 13, 10, 26, 10)), con)
  writeBin(13L, con, size = 4L, endian = "big")
  writeBin(charToRaw("IHDR"), con)
  writeBin(width, con, size = 4L, endian = "big")
  writeBin(height, con, size = 4L, endian = "big")
  invisible(FALSE)
}

write_visium_fixture <- function(
  dir_path,
  barcodes_matrix_order,
  positions,
  counts,
  gene_ids,
  positions_style = c("v2", "v1")
) {
  positions_style <- match.arg(positions_style)

  matrix_dir <- file.path(dir_path, "raw_feature_bc_matrix")
  spatial_dir <- file.path(dir_path, "spatial")
  dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)

  # gene-major mtx, second line a metadata comment like Space Ranger writes
  triplets <- which(counts > 0, arr.ind = TRUE)
  con <- gzfile(file.path(matrix_dir, "matrix.mtx.gz"), "w")
  writeLines("%%MatrixMarket matrix coordinate integer general", con)
  writeLines('%metadata_json: {"software_version": "test"}', con)
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
    data.table::data.table(barcodes_matrix_order),
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

  if (positions_style == "v2") {
    data.table::fwrite(
      positions,
      file.path(spatial_dir, "tissue_positions.csv"),
      col.names = TRUE
    )
  } else {
    data.table::fwrite(
      positions,
      file.path(spatial_dir, "tissue_positions_list.csv"),
      col.names = FALSE
    )
  }

  writeLines(
    jsonlite::toJSON(sf_json, auto_unbox = TRUE),
    file.path(spatial_dir, "scalefactors_json.json")
  )

  write_test_png(
    file.path(spatial_dir, "tissue_hires_image.png"),
    hires_dims[1],
    hires_dims[2]
  )
  write_test_png(
    file.path(spatial_dir, "tissue_lowres_image.png"),
    lowres_dims[1],
    lowres_dims[2]
  )
  write_stub_tiff(
    file.path(spatial_dir, "cytassist_image.tiff"),
    cytassist_dims[1],
    cytassist_dims[2]
  )

  invisible(dir_path)
}

new_test_spatial <- function(name) {
  path <- file.path(test_temp_dir, name)
  unlink(path, recursive = TRUE)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  SpatialSpot(dir_data = path)
}

### data -----------------------------------------------------------------------

set.seed(123L)

random_barcodes <- vapply(
  seq_len(n_spots),
  \(i) {
    paste0(paste(sample(c("A", "C", "G", "T"), 16L, TRUE), collapse = ""), "-1")
  },
  character(1)
)
stopifnot(!anyDuplicated(random_barcodes))

# `barcodes.tsv.gz` is lexicographically sorted, the positions file is in array
# raster order. The fixture uses a cyclic shift between the two, which is a
# derangement: no barcode sits at the same index in both files, so a positional
# zip is wrong for every single spot rather than accidentally right for some.
barcodes_matrix_order <- sort(random_barcodes)
barcodes_raster_order <- barcodes_matrix_order[
  c(seq(2L, n_spots), 1L)
]

positions <- data.table::data.table(
  barcode = barcodes_raster_order,
  in_tissue = 0L,
  array_row = rep(seq_len(n_array_rows) - 1L, each = n_array_cols),
  array_col = rep(seq(0L, by = 2L, length.out = n_array_cols), n_array_rows)
)
positions[, in_tissue := as.integer(array_row %in% c(2L, 3L))]
# pixel coordinates unique per spot, so a mis-join cannot land on the right
# value by accident. The small cross-term stands in for slide rotation.
positions[, pxl_row_in_fullres := 1000L + array_row * 137L + array_col]
positions[, pxl_col_in_fullres := 2000L + array_col * 91L + array_row]
n_in_tissue <- sum(positions$in_tissue)

gene_ids <- sprintf("ENSGTEST%05i", seq_len(n_genes))

set.seed(456L)
counts <- matrix(
  rbinom(n_genes * n_spots, size = 20L, prob = 0.4),
  nrow = n_genes,
  ncol = n_spots
)
counts[counts < 6L] <- 0L
# guarantee every spot and every gene carries at least one count
counts[1L, ] <- 3L
counts[, 1L] <- 3L

visium_dir <- write_visium_fixture(
  file.path(test_temp_dir, "sample_a"),
  barcodes_matrix_order = barcodes_matrix_order,
  positions = positions,
  counts = counts,
  gene_ids = gene_ids
)

## fixture sanity --------------------------------------------------------------

expect_true(
  current = !identical(barcodes_matrix_order, positions$barcode),
  info = "visium io - fixture barcode orders actually differ"
)

expect_equal(
  current = sum(barcodes_matrix_order == positions$barcode),
  target = 0L,
  info = "visium io - fixture barcode orders share no index (derangement)"
)

expect_true(
  current = setequal(barcodes_matrix_order, positions$barcode),
  info = "visium io - fixture barcode sets are identical"
)

## parameters ------------------------------------------------------------------

sp_io_default <- params_sp_visium_io()

expect_true(
  current = sp_io_default$in_tissue_only,
  info = "visium io - in_tissue_only defaults to TRUE"
)

expect_equal(
  current = sp_io_default$matrix_type,
  target = "auto",
  info = "visium io - matrix_type defaults to auto"
)

expect_true(
  current = isTRUE(bixverse:::checkSpVisiumIo(sp_io_default)),
  info = "visium io - checkSpVisiumIo accepts the default params"
)

expect_true(
  current = is.character(bixverse:::checkSpVisiumIo(list(
    in_tissue_only = TRUE
  ))),
  info = "visium io - checkSpVisiumIo rejects an incomplete list"
)

expect_true(
  current = is.character(
    bixverse:::checkSpVisiumIo(
      list(in_tissue_only = "yes", matrix_type = "auto", slide_file = NULL)
    )
  ),
  info = "visium io - checkSpVisiumIo rejects a non-boolean in_tissue_only"
)

expect_true(
  current = is.character(
    bixverse:::checkSpVisiumIo(
      list(in_tissue_only = TRUE, matrix_type = "nope", slide_file = NULL)
    )
  ),
  info = "visium io - checkSpVisiumIo rejects an unknown matrix_type"
)

expect_error(
  current = params_sp_visium_io(in_tissue_only = 1L),
  info = "visium io - params_sp_visium_io rejects a non-boolean flag"
)

## image header parsing --------------------------------------------------------

expect_equal(
  current = unname(bixverse:::.tiff_dimensions(
    file.path(visium_dir, "spatial", "cytassist_image.tiff")
  )),
  target = as.integer(cytassist_dims),
  info = "visium io - TIFF header dimensions parsed"
)

expect_equal(
  current = unname(bixverse:::.png_dimensions(
    file.path(visium_dir, "spatial", "tissue_hires_image.png")
  )),
  target = as.integer(hires_dims),
  info = "visium io - PNG header dimensions parsed"
)

## matrix directory resolution -------------------------------------------------

both_dir <- file.path(test_temp_dir, "both_matrices")
dir.create(
  file.path(both_dir, "filtered_feature_bc_matrix"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  file.path(both_dir, "raw_feature_bc_matrix"),
  recursive = TRUE,
  showWarnings = FALSE
)

expect_equal(
  current = basename(bixverse:::.visium_matrix_dir(both_dir, "auto")),
  target = "filtered_feature_bc_matrix",
  info = "visium io - auto prefers the filtered matrix"
)

expect_equal(
  current = basename(bixverse:::.visium_matrix_dir(both_dir, "raw")),
  target = "raw_feature_bc_matrix",
  info = "visium io - the raw matrix can be requested explicitly"
)

expect_equal(
  current = basename(bixverse:::.visium_matrix_dir(visium_dir, "auto")),
  target = "raw_feature_bc_matrix",
  info = "visium io - auto falls back to the raw matrix"
)

expect_error(
  current = bixverse:::.visium_matrix_dir(visium_dir, "filtered"),
  info = "visium io - a missing filtered matrix errors"
)

## scale factors ---------------------------------------------------------------

scale_factors <- bixverse:::.visium_scale_factors(
  file.path(visium_dir, "spatial"),
  .verbose = FALSE
)

expect_equal(
  current = scale_factors$tissue_hires_scalef,
  target = sf_json$tissue_hires_scalef,
  info = "visium io - hires scale factor passed through untouched"
)

expect_equal(
  current = scale_factors$cytassist_scalef,
  target = expected_cytassist_scalef,
  info = "visium io - CytAssist scale derived from the shipped image size"
)

expect_true(
  current = scale_factors$cytassist_scalef <
    scale_factors$regist_target_img_scalef,
  info = "visium io - CytAssist scale is not regist_target_img_scalef"
)

## loading ---------------------------------------------------------------------

### all spots ------------------------------------------------------------------

object_all <- new_test_spatial("out_all")
object_all <- load_visium(
  object_all,
  visium_dir = visium_dir,
  sp_io_param = params_sp_visium_io(in_tissue_only = FALSE),
  sc_qc_param = qc_permissive,
  exp_id = "sample_a",
  .verbose = FALSE
)

expect_equal(
  current = S7::prop(object_all, "dims"),
  target = c(n_spots, n_genes),
  info = "visium io - all spots and genes land on disk"
)

expect_equal(
  current = length(get_cells_to_keep(object_all)),
  target = n_spots,
  info = "visium io - in_tissue_only = FALSE keeps every spot"
)

expect_equal(
  current = get_sample(object_all, "sample_a")$n_spots,
  target = n_spots,
  info = "visium io - sample records the full spot count"
)

### in tissue only -------------------------------------------------------------

object <- new_test_spatial("out_tissue")
object <- load_visium(
  object,
  visium_dir = visium_dir,
  sp_io_param = params_sp_visium_io(),
  sc_qc_param = qc_permissive,
  exp_id = "sample_a",
  .verbose = FALSE
)

expect_equal(
  current = length(get_cells_to_keep(object)),
  target = n_in_tissue,
  info = "visium io - in_tissue_only drops the off-tissue spots"
)

expect_equal(
  current = length(get_spot_indices_for_exp(object, "sample_a")),
  target = n_in_tissue,
  info = "visium io - spot indices honour the in-tissue filter"
)

expect_equal(
  current = get_sample_ids(object),
  target = "sample_a",
  info = "visium io - the sample is registered under its exp_id"
)

expect_equal(
  current = get_sample(object, "sample_a")$technology,
  target = "visium",
  info = "visium io - technology recorded as visium"
)

expect_equal(
  current = get_sample(object, "sample_a")$n_spots,
  target = n_in_tissue,
  info = "visium io - the sample records its own kept spot count"
)

## the barcode join ------------------------------------------------------------

# The single most important behaviour in the reader. `barcodes.tsv.gz` is
# lexicographic, the positions file is in raster order, and the fixture makes
# the two a derangement of each other. Everything below is asserted per barcode.

obs <- get_sc_duckdb(object_all)$get_obs_table(
  cols = c(
    "cell_idx",
    "cell_id",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres",
    "exp_id"
  )
)
data.table::setorder(obs, cell_idx)

pos_by_barcode <- data.table::copy(positions)
data.table::setkey(pos_by_barcode, barcode)
truth <- pos_by_barcode[obs$cell_id]

expect_equal(
  current = obs$pxl_row_in_fullres,
  target = truth$pxl_row_in_fullres,
  info = "visium io - pxl_row_in_fullres matches the positions file per barcode"
)

expect_equal(
  current = obs$pxl_col_in_fullres,
  target = truth$pxl_col_in_fullres,
  info = "visium io - pxl_col_in_fullres matches the positions file per barcode"
)

expect_equal(
  current = obs$array_row,
  target = truth$array_row,
  info = "visium io - array_row matches the positions file per barcode"
)

expect_equal(
  current = obs$array_col,
  target = truth$array_col,
  info = "visium io - array_col matches the positions file per barcode"
)

expect_equal(
  current = obs$in_tissue,
  target = truth$in_tissue,
  info = "visium io - in_tissue matches the positions file per barcode"
)

expect_equal(
  current = unique(obs$exp_id),
  target = "sample_a",
  info = "visium io - exp_id written into the obs table"
)

# What a positional zip would have produced. If the join were replaced by a
# cbind of barcodes and positions rows, the loaded values would equal this
# instead, and every one of these spots would be wrong.
zipped <- positions[match(obs$cell_id, barcodes_matrix_order)]

expect_equal(
  current = sum(zipped$pxl_row_in_fullres == obs$pxl_row_in_fullres),
  target = 0L,
  info = "visium io - a positional zip disagrees with the join on every spot"
)

expect_true(
  current = !identical(zipped$in_tissue, obs$in_tissue),
  info = "visium io - a positional zip disagrees on in_tissue as well"
)

## coordinate orientation ------------------------------------------------------

# row is the vertical axis, column the horizontal one. Naming is row-first,
# plotting is x-first, so these two are deliberately crossed.
coords <- get_spatial_coords(object, "sample_a")
obs_kept <- get_sc_duckdb(object)$get_obs_table(
  cols = c("cell_idx", "cell_id"),
  filtered = TRUE
)
data.table::setorder(obs_kept, cell_idx)
truth_kept <- pos_by_barcode[obs_kept$cell_id]

expect_equal(
  current = nrow(coords),
  target = n_in_tissue,
  info = "visium io - one coordinate row per kept spot"
)

expect_equal(
  current = as.numeric(coords[, "x"]),
  target = as.numeric(truth_kept$pxl_col_in_fullres),
  info = "visium io - x is pxl_col_in_fullres"
)

expect_equal(
  current = as.numeric(coords[, "y"]),
  target = as.numeric(truth_kept$pxl_row_in_fullres),
  info = "visium io - y is pxl_row_in_fullres"
)

## space ranger 1.x positions --------------------------------------------------

visium_dir_v1 <- write_visium_fixture(
  file.path(test_temp_dir, "sample_v1"),
  barcodes_matrix_order = barcodes_matrix_order,
  positions = positions,
  counts = counts,
  gene_ids = gene_ids,
  positions_style = "v1"
)

expect_true(
  current = file.exists(
    file.path(visium_dir_v1, "spatial", "tissue_positions_list.csv")
  ),
  info = "visium io - 1.x fixture writes the headerless positions file"
)

object_v1 <- new_test_spatial("out_v1")
object_v1 <- load_visium(
  object_v1,
  visium_dir = visium_dir_v1,
  sp_io_param = params_sp_visium_io(in_tissue_only = FALSE),
  sc_qc_param = qc_permissive,
  exp_id = "sample_a",
  .verbose = FALSE
)

obs_v1 <- get_sc_duckdb(object_v1)$get_obs_table(
  cols = c("cell_idx", "cell_id", "pxl_row_in_fullres", "pxl_col_in_fullres")
)
data.table::setorder(obs_v1, cell_idx)

expect_equal(
  current = obs_v1$pxl_row_in_fullres,
  target = obs$pxl_row_in_fullres,
  info = "visium io - headerless 1.x positions parse to the same coordinates"
)

expect_equal(
  current = obs_v1$cell_id,
  target = obs$cell_id,
  info = "visium io - headerless 1.x positions keep the same barcodes"
)

## images ----------------------------------------------------------------------

image_paths <- get_sample(object, "sample_a")$image_paths

expect_equal(
  current = sort(names(image_paths)),
  target = c("cytassist", "hires", "lowres"),
  info = "visium io - the three shipped images are registered"
)

expect_true(
  current = all(
    grepl("^images/sample_a/", unlist(image_paths))
  ),
  info = "visium io - image paths are relative to the data directory"
)

expect_true(
  current = all(file.exists(
    file.path(S7::prop(object, "dir_data"), unlist(image_paths))
  )),
  info = "visium io - the images were copied into the data directory"
)

expect_error(
  current = get_image(object, "sample_a", "cytassist"),
  info = "visium io - get_image refuses the CytAssist TIFF"
)

expect_error(
  current = get_image(object, "sample_a", "fullres"),
  info = "visium io - get_image errors on an unregistered resolution"
)

if (requireNamespace("png", quietly = TRUE)) {
  lowres <- get_image(object, "sample_a", "lowres")
  expect_equal(
    current = dim(lowres)[1:2],
    target = c(lowres_dims[2], lowres_dims[1]),
    info = "visium io - the lowres PNG loads at the expected size"
  )
}

### no copy --------------------------------------------------------------------

object_nocopy <- new_test_spatial("out_nocopy")
object_nocopy <- load_visium(
  object_nocopy,
  visium_dir = visium_dir,
  sp_io_param = params_sp_visium_io(),
  sc_qc_param = qc_permissive,
  exp_id = "sample_a",
  copy_images = FALSE,
  .verbose = FALSE
)

expect_equal(
  current = length(get_sample(object_nocopy, "sample_a")$image_paths),
  target = 0L,
  info = "visium io - copy_images = FALSE registers no images"
)

expect_true(
  current = !dir.exists(
    file.path(S7::prop(object_nocopy, "dir_data"), "images")
  ),
  info = "visium io - copy_images = FALSE touches nothing on disk"
)

### user supplied slide --------------------------------------------------------

slide_path <- file.path(test_temp_dir, "slide.png")
write_test_png(slide_path, 40L, 30L)

object_slide <- new_test_spatial("out_slide")
object_slide <- load_visium(
  object_slide,
  visium_dir = visium_dir,
  sp_io_param = params_sp_visium_io(slide_file = slide_path),
  sc_qc_param = qc_permissive,
  exp_id = "sample_a",
  .verbose = FALSE
)

expect_true(
  current = "fullres" %in%
    names(get_sample(object_slide, "sample_a")$image_paths),
  info = "visium io - a user-supplied slide registers under fullres"
)

## multiple samples ------------------------------------------------------------

visium_dir_b <- write_visium_fixture(
  file.path(test_temp_dir, "sample_b"),
  barcodes_matrix_order = barcodes_matrix_order,
  positions = positions,
  counts = counts,
  gene_ids = gene_ids
)

### the wrong way ------------------------------------------------------------

# `load_visium()` rewrites the counts, the obs table and the var table of the
# whole store, while `add_spatial_sample()` merges. A second `exp_id` therefore
# used to leave an object claiming two samples over one sample's data, with the
# first sample's cached results still keyed against a gene axis that had just
# been replaced. It has to refuse instead.
expect_error(
  current = load_visium(
    object,
    visium_dir = visium_dir_b,
    sp_io_param = params_sp_visium_io(),
    sc_qc_param = qc_permissive,
    exp_id = "sample_b",
    .verbose = FALSE
  ),
  info = "visium io - a second exp_id through load_visium errors"
)

expect_equal(
  current = get_sample_ids(object),
  target = "sample_a",
  info = "visium io - and the refused load leaves the registration alone"
)

# re-running the same exp_id is the QC tweak loop and has to keep working
expect_equal(
  current = get_sample_ids(load_visium(
    object,
    visium_dir = visium_dir,
    sp_io_param = params_sp_visium_io(),
    sc_qc_param = qc_permissive,
    exp_id = "sample_a",
    .verbose = FALSE
  )),
  target = "sample_a",
  info = "visium io - re-running the same exp_id is still allowed"
)

prescan <- prescan_visium_dirs(
  dirs = c(visium_dir, visium_dir_b),
  exp_ids = c("sample_a", "sample_b"),
  .verbose = FALSE
)

expect_equal(
  current = prescan$universe_size,
  target = n_genes,
  info = "visium io - the prescan universe is the gene intersection"
)

expect_equal(
  current = names(prescan$file_tasks),
  target = c("sample_a", "sample_b"),
  info = "visium io - prescan tasks are keyed by exp_id"
)

object_multi <- new_test_spatial("out_multi")
object_multi <- load_multi_visium(
  object_multi,
  prescan_result = prescan,
  sp_io_param = params_sp_visium_io(),
  sc_qc_param = qc_permissive,
  .verbose = FALSE
)

expect_equal(
  current = S7::prop(object_multi, "dims"),
  target = c(2L * n_spots, n_genes),
  info = "visium io - the multi load stacks both inputs"
)

expect_equal(
  current = sort(get_sample_ids(object_multi)),
  target = c("sample_a", "sample_b"),
  info = "visium io - both samples are registered"
)

# the invariant the load_visium guard exists to keep: what the object claims to
# hold and what the store actually holds
expect_equal(
  current = sort(get_sample_ids(object_multi)),
  target = sort(unique(get_sc_obs(object_multi, filtered = TRUE)$exp_id)),
  info = "visium io - the registered samples match the exp_ids in obs"
)

expect_equal(
  current = purrr::map_int(
    get_samples(object_multi),
    \(s) s$n_spots
  )[["sample_b"]],
  target = n_in_tissue,
  info = "visium io - n_spots is counted per sample, not object-wide"
)

expect_equal(
  current = length(get_spot_indices_for_exp(object_multi, "sample_b")),
  target = n_in_tissue,
  info = "visium io - the in-tissue filter applies per sample in a multi load"
)

obs_multi <- get_sc_duckdb(object_multi)$get_obs_table(
  cols = c("cell_idx", "cell_id", "exp_id", "pxl_row_in_fullres")
)
data.table::setorder(obs_multi, cell_idx)
truth_multi <- pos_by_barcode[obs_multi$cell_id]

expect_equal(
  current = obs_multi$pxl_row_in_fullres,
  target = truth_multi$pxl_row_in_fullres,
  info = "visium io - the barcode join holds across a multi load"
)

expect_error(
  current = load_multi_visium(
    new_test_spatial("out_multi_err"),
    prescan_result = prescan,
    sp_io_param = params_sp_visium_io(slide_file = slide_path),
    sc_qc_param = qc_permissive,
    .verbose = FALSE
  ),
  info = "visium io - a slide file is refused in the multi load"
)

## failure modes ---------------------------------------------------------------

empty_dir <- file.path(test_temp_dir, "empty")
dir.create(empty_dir, showWarnings = FALSE)

expect_error(
  current = load_visium(
    new_test_spatial("out_empty"),
    visium_dir = empty_dir,
    sc_qc_param = qc_permissive,
    exp_id = "nope",
    .verbose = FALSE
  ),
  info = "visium io - a directory without a count matrix errors"
)

mismatch_dir <- file.path(test_temp_dir, "mismatch")
mismatch_positions <- data.table::copy(positions)
mismatch_positions[1L, barcode := "NOTAREALBARCODE-1"]
write_visium_fixture(
  mismatch_dir,
  barcodes_matrix_order = barcodes_matrix_order,
  positions = mismatch_positions,
  counts = counts,
  gene_ids = gene_ids
)

expect_error(
  current = load_visium(
    new_test_spatial("out_mismatch"),
    visium_dir = mismatch_dir,
    sc_qc_param = qc_permissive,
    exp_id = "sample_a",
    .verbose = FALSE
  ),
  info = "visium io - a barcode missing from the positions file errors"
)
