# histology image features -----------------------------------------------------

# The image path is the one part of the spatial stack with a hard system
# dependency: the bindings need the `spatial-image` feature of `bixverse-rs`,
# which needs OpenSlide at build and run time. A build without it has no
# working `rs_sp_image_*`, so everything here sits behind a probe that actually
# calls into the backend rather than behind a version check.

source("helper_sp.R", local = TRUE)

test_temp_dir <- file.path(tempdir(), "sp_image")
dir.create(test_temp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir))

image_backend <- sp_image_backend_ok()

if (!image_backend) {
  exit_file("sp image - no Rust image backend (OpenSlide) in this build")
}

## metadata --------------------------------------------------------------------

# A plain PNG reports one pyramid level, a downsample of one and no vendor,
# so the shape of the answer is the same as for a whole-slide image.
flat_png <- file.path(test_temp_dir, "flat.png")
png::writePNG(
  array(stats::runif(40L * 24L * 3L), dim = c(24L, 40L, 3L)),
  target = flat_png
)
meta <- rs_sp_image_metadata(flat_png)

expect_equal(
  current = names(meta),
  target = c(
    "width",
    "height",
    "n_levels",
    "level_dims",
    "downsamples",
    "vendor"
  ),
  info = "sp image - the metadata shape is fixed"
)

expect_equal(
  current = c(meta$width, meta$height),
  target = c(40L, 24L),
  info = "sp image - width is the horizontal axis, height the vertical one"
)

expect_equal(
  current = meta$n_levels,
  target = 1L,
  info = "sp image - a plain PNG is a single pyramid level"
)

expect_equal(
  current = as.numeric(meta$downsamples),
  target = 1,
  info = "sp image - and is not downsampled against itself"
)

expect_true(
  current = is.na(meta$vendor),
  info = "sp image - a plain PNG has no OpenSlide vendor"
)

expect_equal(
  current = dim(meta$level_dims),
  target = c(1L, 2L),
  info = "sp image - level_dims is levels x 2"
)

expect_error(
  current = rs_sp_image_metadata(
    file.path(test_temp_dir, "does_not_exist.png")
  ),
  info = "sp image - a missing file errors rather than returning nothing"
)

## features off a fixture ------------------------------------------------------

fixture <- sp_make_visium(
  file.path(test_temp_dir, "sample_a"),
  layout = "hex",
  n_rows = 8L,
  n_cols = 8L,
  n_genes = 20L,
  seed = 303L,
  image_size = 200L
)
object <- sp_load_visium_fixture(fixture, file.path(test_temp_dir, "out"))

object <- image_features_sp(object, resolution = "hires", .verbose = FALSE)
features <- get_per_sample_image_features(object, "sample_a", "hires")

expected_features <- c(
  "r_mean",
  "r_sd",
  "r_entropy",
  "g_mean",
  "g_sd",
  "g_entropy",
  "b_mean",
  "b_sd",
  "b_entropy",
  "grey_mean",
  "grey_sd",
  "grey_entropy",
  "haem_mean",
  "haem_sd",
  "haem_entropy",
  "eosin_mean",
  "eosin_sd",
  "eosin_entropy",
  "mean_od",
  "glcm_grey_contrast",
  "glcm_grey_homogeneity",
  "glcm_grey_correlation",
  "glcm_grey_energy",
  "glcm_grey_entropy",
  "glcm_haem_contrast",
  "glcm_haem_homogeneity",
  "glcm_haem_correlation",
  "glcm_haem_energy",
  "glcm_haem_entropy"
)

expect_equal(
  current = features$feature_names,
  target = expected_features,
  info = "sp image - the feature block is 29 named columns"
)

expect_equal(
  current = dim(features$values),
  target = c(fixture$n_spots, length(expected_features)),
  info = "sp image - one row per spot, one column per feature"
)

expect_equal(
  current = features$params$n_dropped,
  target = 0L,
  info = "sp image - a 200 px image covers the whole fixture frame"
)

expect_equal(
  current = colnames(features$values),
  target = features$feature_names,
  info = "sp image - the value matrix carries the feature names"
)

expect_equal(
  current = as.integer(rownames(features$values)),
  target = features$cell_idx,
  info = "sp image - and the cell_idx as row names"
)

expect_true(
  current = all(is.finite(features$values)),
  info = "sp image - no feature comes back as NA or infinite"
)

# 8-bit scale, not the unit interval. Worth pinning: everything downstream
# that standardises these features silently assumes one or the other.
channel_means <- features$values[, c("r_mean", "g_mean", "b_mean")]

expect_true(
  current = all(channel_means >= 0 & channel_means <= 255),
  info = "sp image - channel means are on the 8-bit scale"
)

expect_true(
  current = max(channel_means) > 1,
  info = "sp image - and are not rescaled to the unit interval"
)

# optical density is -log10(I / I_max), so a mid-grey tile lands near 0.3
expect_true(
  current = all(features$values[, "mean_od"] > 0),
  info = "sp image - mean optical density is positive"
)

### the join back onto obs -----------------------------------------------------

feature_dt <- get_data(features)

expect_true(
  current = isTRUE(attr(feature_dt, "is_obs")),
  info = "sp image - get_data marks the table for add_sc_new_obs"
)

expect_equal(
  current = names(feature_dt)[1L],
  target = "cell_idx",
  info = "sp image - cell_idx comes first"
)

joined <- add_sc_new_obs(object, feature_dt)
obs <- get_sc_obs(joined, filtered = TRUE)

expect_true(
  current = all(expected_features %in% names(obs)),
  info = "sp image - the features land in the obs table"
)

data.table::setorderv(obs, "cell_idx")

expect_equal(
  current = as.numeric(obs$grey_mean),
  target = as.numeric(
    features$values[as.character(obs$cell_idx), "grey_mean"]
  ),
  tolerance = 1e-6,
  info = "sp image - and land against the right spot"
)

### a selection ----------------------------------------------------------------

expect_equal(
  current = names(get_data(features, columns = c("mean_od", "r_sd"))),
  target = c("cell_idx", "mean_od", "r_sd"),
  info = "sp image - get_data honours a column selection"
)

expect_error(
  current = get_data(features, columns = "not_a_feature"),
  info = "sp image - an unknown feature errors"
)

## parameters bite -------------------------------------------------------------

wide_tile <- get_per_sample_image_features(
  image_features_sp(
    object,
    resolution = "hires",
    image_params = params_sp_image(tile_scale = 3.0),
    .verbose = FALSE
  ),
  "sample_a",
  "hires"
)

expect_true(
  current = !isTRUE(all.equal(
    wide_tile$values[, "grey_sd"],
    features$values[, "grey_sd"]
  )),
  info = "sp image - tile_scale changes what gets measured"
)

coarse_glcm <- get_per_sample_image_features(
  image_features_sp(
    object,
    resolution = "hires",
    image_params = params_sp_image(glcm_levels = 4L),
    .verbose = FALSE
  ),
  "sample_a",
  "hires"
)

expect_true(
  current = !isTRUE(all.equal(
    coarse_glcm$values[, "glcm_grey_entropy"],
    features$values[, "glcm_grey_entropy"]
  )),
  info = "sp image - glcm_levels changes the texture features"
)

expect_equal(
  current = as.numeric(coarse_glcm$values[, "grey_mean"]),
  target = as.numeric(features$values[, "grey_mean"]),
  tolerance = 1e-6,
  info = "sp image - but leaves the first-order statistics alone"
)

swapped_stain <- get_per_sample_image_features(
  image_features_sp(
    object,
    resolution = "hires",
    image_params = params_sp_image(
      stain_haem = c(0.2, 0.8, 0.5),
      stain_eosin = c(0.6, 0.3, 0.7)
    ),
    .verbose = FALSE
  ),
  "sample_a",
  "hires"
)

expect_true(
  current = !isTRUE(all.equal(
    swapped_stain$values[, "haem_mean"],
    features$values[, "haem_mean"]
  )),
  info = "sp image - a custom stain basis changes the deconvolved channels"
)

expect_equal(
  current = as.numeric(swapped_stain$values[, "r_mean"]),
  target = as.numeric(features$values[, "r_mean"]),
  tolerance = 1e-6,
  info = "sp image - and leaves the raw colour channels alone"
)

## resolutions -----------------------------------------------------------------

# The scale factor selects both the pixel space and the effective tile size,
# so the same slide at two resolutions is not the same measurement.
object <- image_features_sp(object, resolution = "lowres", .verbose = FALSE)
cache <- get_sp_cache(object)

expect_equal(
  current = sort(names(cache$per_sample_image_features$sample_a)),
  target = c("hires", "lowres"),
  info = "sp image - features are cached per resolution, not overwritten"
)

expect_error(
  current = image_features_sp(
    object,
    resolution = "fullres",
    .verbose = FALSE
  ),
  info = "sp image - an unregistered resolution errors"
)

## tiles falling off the image -------------------------------------------------

# Spots whose tile leaves the image are dropped inside Rust rather than
# failing the run, so `cell_idx` gets shorter than the sample.
small_fixture <- sp_make_visium(
  file.path(test_temp_dir, "small_image"),
  layout = "hex",
  n_rows = 8L,
  n_cols = 8L,
  n_genes = 20L,
  seed = 404L,
  image_size = 60L
)
small_object <- sp_load_visium_fixture(
  small_fixture,
  file.path(test_temp_dir, "out_small")
)
small_object <- image_features_sp(
  small_object,
  resolution = "hires",
  .verbose = FALSE
)
small_features <- get_per_sample_image_features(
  small_object,
  "sample_a",
  "hires"
)

expect_true(
  current = small_features$params$n_dropped > 0L,
  info = "sp image - spots off a too-small image are dropped"
)

expect_equal(
  current = length(small_features$cell_idx) +
    small_features$params$n_dropped,
  target = small_fixture$n_spots,
  info = "sp image - dropped plus kept accounts for every spot"
)

expect_equal(
  current = nrow(small_features$values),
  target = length(small_features$cell_idx),
  info = "sp image - the value matrix shrinks with the spots"
)

expect_true(
  current = !is.unsorted(small_features$cell_idx),
  info = "sp image - the surviving cell_idx stay in ascending order"
)

## the local example data ------------------------------------------------------

# The two honest facts about the image path, both only checkable on a real
# Space Ranger output:
#
#   1. A spot is 14.4 px across in the hires PNG and 21.6 px in the CytAssist
#      frame. Haralick texture on a 14 px tile is close to noise.
#   2. `image_features_sp(resolution = "cytassist")` works on data where
#      `get_image(..., "cytassist")` errors, because the Rust side reads TIFF
#      and the R side does not.
example_dir <- sp_example_visium_dir()

if (at_home() && !is.null(example_dir)) {
  ex_out <- file.path(test_temp_dir, "example")
  unlink(ex_out, recursive = TRUE)
  dir.create(ex_out, recursive = TRUE, showWarnings = FALSE)

  ex_object <- load_visium(
    SpatialSpot(dir_data = ex_out),
    visium_dir = example_dir,
    sc_qc_param = params_sc_min_quality(
      min_unique_genes = 0L,
      min_lib_size = 0L,
      min_cells = 0L
    ),
    exp_id = "example",
    copy_images = TRUE,
    .verbose = FALSE
  )
  ex_sample <- get_sample(ex_object, "example")
  spot_fullres <- ex_sample$scale_factors[["spot_diameter_fullres"]]

  expect_equal(
    current = spot_fullres *
      bixverse:::.sp_image_scalef(ex_sample, "hires"),
    target = 14.4,
    tolerance = 0.05,
    info = "sp image - a spot is about 14.4 px across in the hires PNG"
  )

  expect_equal(
    current = spot_fullres *
      bixverse:::.sp_image_scalef(ex_sample, "cytassist"),
    target = 21.6,
    tolerance = 0.05,
    info = "sp image - and about 21.6 px in the CytAssist frame"
  )

  expect_error(
    current = get_image(ex_object, "example", "cytassist"),
    info = "sp image - the R side cannot load the CytAssist TIFF"
  )

  cytassist_meta <- rs_sp_image_metadata(
    file.path(
      S7::prop(ex_object, "dir_data"),
      ex_sample$image_paths[["cytassist"]]
    )
  )

  expect_true(
    current = cytassist_meta$width > 0L && cytassist_meta$height > 0L,
    info = "sp image - but the Rust side reads its header"
  )

  ex_object <- image_features_sp(
    ex_object,
    resolution = "cytassist",
    .verbose = FALSE
  )
  cyto_features <- get_per_sample_image_features(
    ex_object,
    "example",
    "cytassist"
  )

  expect_equal(
    current = length(cyto_features$feature_names),
    target = 29L,
    info = "sp image - and cuts tiles out of it all the same"
  )

  expect_true(
    current = all(is.finite(cyto_features$values)),
    info = "sp image - the CytAssist features are usable numbers"
  )
}
