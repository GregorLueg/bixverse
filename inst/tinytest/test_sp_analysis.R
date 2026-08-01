# spatial subset and analysis methods ------------------------------------------

source("helper_sp.R", local = TRUE)

test_temp_dir <- file.path(tempdir(), "sp_analysis")
dir.create(test_temp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir))

## fixture ---------------------------------------------------------------------

# Standard Visium geometry: `array_col` steps by two inside a row and the rows
# are offset by one, so an interior spot has exactly six neighbours. Anything
# that gets the parity wrong drops to four.
#
# The barcode ordering trap comes with the fixture: `barcodes.tsv.gz` is
# lexicographic while `tissue_positions.csv` is in array raster order, so
# anything that joins the two by position scrambles the coordinates.
n_array_rows <- 12L
n_genes <- 60L

fixture <- sp_make_visium(
  file.path(test_temp_dir, "sample_a"),
  layout = "hex",
  n_rows = n_array_rows,
  n_cols = 12L,
  n_genes = n_genes,
  seed = 11L,
  image_size = 200L
)

n_spots <- fixture$n_spots
gene_ids <- fixture$gene_ids

# the image assertions need a readable PNG *and* the Rust image backend behind
# it. Gating on `png` alone hard-errors on a build without OpenSlide.
image_ok <- fixture$png_written && sp_image_backend_ok()

object <- sp_load_visium_fixture(
  fixture,
  file.path(test_temp_dir, "out"),
  exp_id = "sample_a"
)

## parameters ------------------------------------------------------------------

expect_true(
  current = isTRUE(bixverse:::checkSpGraph(params_sp_graph())),
  info = "sp analysis - checkSpGraph accepts the defaults"
)

expect_true(
  current = is.character(bixverse:::checkSpGraph(list(layout = "hex"))),
  info = "sp analysis - checkSpGraph rejects an incomplete list"
)

expect_true(
  current = is.character(bixverse:::checkSpGraph(
    utils::modifyList(params_sp_graph(), list(layout = "nope"))
  )),
  info = "sp analysis - checkSpGraph rejects an unknown layout"
)

expect_true(
  current = is.character(bixverse:::checkSpGraph(
    utils::modifyList(params_sp_graph(), list(connectivity = 5L))
  )),
  info = "sp analysis - checkSpGraph rejects a connectivity other than 4 or 8"
)

expect_error(
  current = params_sp_graph(layout = "radius"),
  info = "sp analysis - a radius layout without a radius errors"
)

expect_true(
  current = isTRUE(bixverse:::checkSpSvg(params_sp_svg())),
  info = "sp analysis - checkSpSvg accepts the defaults"
)

expect_true(
  current = is.character(bixverse:::checkSpSvg(list(assay = "counts"))),
  info = "sp analysis - checkSpSvg rejects an unknown assay"
)

expect_true(
  current = isTRUE(bixverse:::checkSpSparkx(params_sp_sparkx())),
  info = "sp analysis - checkSpSparkx accepts the defaults"
)

expect_true(
  current = isTRUE(bixverse:::checkSpSparkx(params_sp_sparkx(
    kernels = list(list(kernel = "gaussian", bandwidth = 10))
  ))),
  info = "sp analysis - checkSpSparkx accepts an explicit kernel bank"
)

expect_true(
  current = is.character(bixverse:::checkSpSparkx(
    utils::modifyList(
      params_sp_sparkx(),
      list(kernels = list(list(kernel = "triangle", bandwidth = 10)))
    )
  )),
  info = "sp analysis - checkSpSparkx rejects an unknown kernel"
)

expect_true(
  current = isTRUE(bixverse:::checkSpNhood(params_sp_nhood())),
  info = "sp analysis - checkSpNhood accepts the defaults"
)

expect_true(
  current = isTRUE(bixverse:::checkSpImage(params_sp_image())),
  info = "sp analysis - checkSpImage accepts the defaults"
)

expect_error(
  current = params_sp_image(glcm_offsets_dy = 1L),
  info = "sp analysis - half a GLCM offset pair errors"
)

expect_true(
  current = is.character(bixverse:::checkSpImage(
    utils::modifyList(params_sp_image(), list(stain_haem = c(1, 0, 0)))
  )),
  info = "sp analysis - checkSpImage rejects half a stain basis"
)

## the subset class ------------------------------------------------------------

subset_obj <- SpatialSpotSubset(object, exp_id = "sample_a")

expect_true(
  current = S7::S7_inherits(subset_obj, SingleCellsSubset),
  info = "sp analysis - SpatialSpotSubset inherits from SingleCellsSubset"
)

expect_equal(
  current = S7::prop(subset_obj, "dims"),
  target = c(n_spots, n_genes),
  info = "sp analysis - the subset carries every spot of the experiment"
)

expect_equal(
  current = S7::prop(subset_obj, "grouping_column"),
  target = "exp_id",
  info = "sp analysis - the grouping column is exp_id"
)

expect_true(
  current = identical(
    S7::prop(subset_obj, "count_connection"),
    S7::prop(object, "count_connection")
  ),
  info = "sp analysis - the Rust count connection is shared, not copied"
)

subset_to_original <- S7::prop(subset_obj, "subset_to_original")

expect_false(
  current = is.unsorted(subset_to_original),
  info = "sp analysis - subset_to_original is ascending"
)

expect_equal(
  current = as.integer(get_spot_indices_for_exp(object, "sample_a")),
  target = subset_to_original - 1L,
  info = "sp analysis - subset_to_original is the parent 1-indexed position"
)

expect_equal(
  current = as.integer(get_spot_indices_for_exp(subset_obj, "sample_a")),
  target = as.integer(get_spot_indices_for_exp(object, "sample_a")),
  info = "sp analysis - the subset reports the same global 0-based indices"
)

# Anchored on the fixture's own positions file, matched by barcode. Comparing
# the subset's cached coords with `get_spatial_coords(object, ...)` says
# nothing: the subset property was assigned from that very call, so both sides
# move together whatever the parent does.
obs_parent <- sp_obs_in_graph_order(object, "sample_a")
coords_parent <- get_spatial_coords(object, "sample_a")
truth <- fixture$positions[match(
  obs_parent[["cell_id"]],
  fixture$positions[["barcode"]]
)]

expect_equal(
  current = as.numeric(coords_parent[, "x"]),
  target = as.numeric(truth[["pxl_col_in_fullres"]]),
  info = "sp analysis - the coords x come off the positions file by barcode"
)

expect_equal(
  current = as.numeric(coords_parent[, "y"]),
  target = as.numeric(truth[["pxl_row_in_fullres"]]),
  info = "sp analysis - and so do the coords y"
)

expect_equal(
  current = get_spatial_coords(subset_obj, "sample_a"),
  target = coords_parent,
  info = "sp analysis - the subset carries the parent's coordinates"
)

# coordinates have to line up with the obs rows, per spot, not just in count
obs_sub <- get_sc_obs(subset_obj)
coords_sub <- S7::prop(subset_obj, "coords")

expect_equal(
  current = as.numeric(coords_sub[, "x"]),
  target = as.numeric(obs_sub$pxl_col_in_fullres),
  info = "sp analysis - coords x lines up with the subset obs rows"
)

expect_equal(
  current = as.numeric(coords_sub[, "y"]),
  target = as.numeric(obs_sub$pxl_row_in_fullres),
  info = "sp analysis - coords y lines up with the subset obs rows"
)

expect_equal(
  current = get_sample_ids(subset_obj),
  target = "sample_a",
  info = "sp analysis - the subset knows its one experiment"
)

expect_error(
  current = get_sample(subset_obj, "sample_b"),
  info = "sp analysis - asking the subset for another experiment errors"
)

expect_error(
  current = SpatialSpotSubset(object, exp_id = "sample_b"),
  info = "sp analysis - an unknown exp_id errors at construction"
)

## adjacency round trip --------------------------------------------------------

indices_in <- list(c(1L, 2L), c(0L, 2L), c(0L, 1L))
weights_in <- list(c(1, 2), c(3, 4), c(5, 6))
sparse <- bixverse:::.sp_adjacency_to_sparse(indices_in, weights_in, 3L)
round_trip <- bixverse:::.sp_sparse_to_adjacency(sparse)

expect_equal(
  current = round_trip$indices,
  target = indices_in,
  info = "sp analysis - the adjacency indices survive the sparse round trip"
)

expect_equal(
  current = round_trip$weights,
  target = weights_in,
  info = "sp analysis - the adjacency weights survive the sparse round trip"
)

## the spatial graph -----------------------------------------------------------

subset_obj <- build_spatial_graph_sp(subset_obj, .verbose = FALSE)
graph <- get_per_sample_spatial_graph(subset_obj, "sample_a")

expect_true(
  current = inherits(graph, "CsparseMatrix"),
  info = "sp analysis - the graph is cached as a sparse matrix"
)

expect_equal(
  current = dim(graph),
  target = c(n_spots, n_spots),
  info = "sp analysis - the graph is spots x spots"
)

adjacency <- bixverse:::.sp_sparse_to_adjacency(graph)
neighbour_counts <- lengths(adjacency$indices)

expect_equal(
  current = max(neighbour_counts),
  target = 6L,
  info = "sp analysis - an interior hex spot has six neighbours"
)

expect_true(
  current = all(graph@x == 1),
  info = "sp analysis - binary weighting gives weights of one"
)

expect_true(
  current = isTRUE(all.equal(
    as.matrix(graph),
    as.matrix(Matrix::t(graph))
  )),
  info = "sp analysis - the lattice adjacency is reciprocal"
)

# `degree` out of Rust is the row sum of W + t(W), not the neighbour count. On
# a reciprocal binary lattice that is exactly twice it.
sym_degree <- Matrix::rowSums(graph + Matrix::t(graph))

expect_equal(
  current = as.numeric(sym_degree),
  target = as.numeric(2L * neighbour_counts),
  info = "sp analysis - symmetrised degree is twice the neighbour count"
)

expect_error(
  current = morans_i_sp(
    SpatialSpotSubset(object, exp_id = "sample_a"),
    .verbose = FALSE
  ),
  info = "sp analysis - Moran's I without a graph errors"
)

## moran's i -------------------------------------------------------------------

subset_obj <- morans_i_sp(subset_obj, .verbose = FALSE)
morans <- get_per_sample_morans_i(subset_obj, "sample_a")

expect_true(
  current = data.table::is.data.table(morans),
  info = "sp analysis - the Moran's I result is still a data.table"
)

expect_true(
  current = inherits(morans, "SpMoransRes"),
  info = "sp analysis - the Moran's I result carries its own class"
)

expect_equal(
  current = nrow(morans),
  target = n_genes,
  info = "sp analysis - every gene was tested"
)

# `get_gene_names()` is where `gene_id` is taken from in the first place, so
# comparing the two is a tautology. The features file is the external
# reference: a permuted var read breaks this and nothing else in the suite.
expect_equal(
  current = morans$gene_id,
  target = fixture$gene_ids,
  info = "sp analysis - the gene labels are the features file in file order"
)

expect_equal(
  current = morans$gene_idx,
  target = seq_len(n_genes) - 1L,
  info = "sp analysis - paired with the on-disk index they were computed at"
)

expect_true(
  current = attr(morans, "e_i") < 0,
  info = "sp analysis - the null expectation is negative"
)

expect_true(
  current = morans[gene_id == gene_ids[1L]]$morans_i > 0.5,
  info = "sp analysis - the planted north/south gradient scores high"
)

expect_true(
  current = morans[gene_id == gene_ids[1L]]$fdr < 0.01,
  info = "sp analysis - the planted gradient is significant"
)

expect_true(
  current = abs(morans[gene_id == gene_ids[2L]]$morans_i) <
    0.25 * morans[gene_id == gene_ids[1L]]$morans_i,
  info = "sp analysis - the noise gene scores far below the planted one"
)

## sparkx ----------------------------------------------------------------------

subset_obj <- sparkx_sp(subset_obj, .verbose = FALSE)
sparkx <- get_per_sample_sparkx(subset_obj, "sample_a")

expect_true(
  current = inherits(sparkx, "SpSparkxRes"),
  info = "sp analysis - the SPARK-X result carries its own class"
)

expect_equal(
  current = nrow(sparkx),
  target = n_genes,
  info = "sp analysis - SPARK-X tested every gene"
)

expect_equal(
  current = dim(attr(sparkx, "stat_per_kernel")),
  target = c(n_genes, length(attr(sparkx, "kernels"))),
  info = "sp analysis - the per-kernel matrix is genes x kernels"
)

expect_true(
  current = sparkx[gene_id == gene_ids[1L]]$fdr < 0.01,
  info = "sp analysis - SPARK-X picks up the planted gradient"
)

## neighbourhood enrichment ----------------------------------------------------

# The labels come off the fixture's positions file by barcode, not out of
# `.sp_obs_for_exp()`. Building them with the same helper `nhood_enrichment_sp`
# reads them back with means any row-ordering change permutes the write and the
# read identically and the assertions below never notice.
obs_labels <- get_sc_obs(subset_obj)
data.table::setorderv(obs_labels, "cell_idx")
label_rows <- fixture$positions[match(
  obs_labels[["cell_id"]],
  fixture$positions[["barcode"]]
)]
subset_obj[["band"]] <- ifelse(
  label_rows[["array_row"]] < n_array_rows / 2L,
  "north",
  "south"
)

subset_obj <- nhood_enrichment_sp(
  subset_obj,
  label_col = "band",
  nhood_params = params_sp_nhood(n_perm = 200L),
  seed = 7L,
  .verbose = FALSE
)
nhood <- get_per_sample_nhood_enrichment(subset_obj, "sample_a", "band")

expect_true(
  current = inherits(nhood, "SpNhoodRes"),
  info = "sp analysis - the enrichment result carries its own class"
)

expect_equal(
  current = nhood$label_levels,
  target = c("north", "south"),
  info = "sp analysis - the label levels come back in encoding order"
)

expect_equal(
  current = dim(nhood$z_scores),
  target = c(2L, 2L),
  info = "sp analysis - the Z matrix is K x K"
)

expect_true(
  current = all(diag(nhood$z_scores) > 0),
  info = "sp analysis - two spatially separated bands self-associate"
)

expect_true(
  current = all(nhood$z_scores[1L, 2L] < 0),
  info = "sp analysis - the two bands avoid each other"
)

expect_error(
  current = nhood_enrichment_sp(
    subset_obj,
    label_col = "not_a_column",
    .verbose = FALSE
  ),
  info = "sp analysis - an unknown label column errors"
)

## image features --------------------------------------------------------------

if (image_ok) {
  subset_obj <- image_features_sp(
    subset_obj,
    resolution = "hires",
    .verbose = FALSE
  )
  features <- get_per_sample_image_features(subset_obj, "sample_a", "hires")

  expect_true(
    current = inherits(features, "SpImageFeatures"),
    info = "sp analysis - the image feature result carries its own class"
  )

  expect_equal(
    current = nrow(features$values),
    target = length(features$cell_idx),
    info = "sp analysis - one feature row per surviving spot"
  )

  expect_true(
    current = all(features$cell_idx %in% subset_to_original),
    info = "sp analysis - the feature cell_idx sit in the parent index space"
  )

  feature_dt <- get_data(features)

  expect_true(
    current = isTRUE(attr(feature_dt, "is_obs")),
    info = "sp analysis - get_data marks the table for add_sc_new_obs"
  )

  expect_equal(
    current = names(feature_dt)[1L],
    target = "cell_idx",
    info = "sp analysis - get_data puts cell_idx first"
  )

  expect_error(
    current = image_features_sp(
      subset_obj,
      resolution = "fullres",
      .verbose = FALSE
    ),
    info = "sp analysis - an unregistered resolution errors"
  )
}

## the cache -------------------------------------------------------------------

sp_cache <- get_sp_cache(subset_obj)

expect_true(
  current = "per_sample_image_features" %in% names(sp_cache),
  info = "sp analysis - the SpCache has the image feature slot"
)

expect_equal(
  current = names(sp_cache$per_sample_nhood_enrichment$sample_a),
  target = "band",
  info = "sp analysis - enrichment nests by exp_id then label column"
)

if (image_ok) {
  expect_equal(
    current = names(sp_cache$per_sample_image_features$sample_a),
    target = "hires",
    info = "sp analysis - image features nest by exp_id then resolution"
  )
}

# The sp_cache setters and getters are S3 generics targeting S7 classes, so
# they only dispatch when registered through S7::method(). A plain
# `set_per_sample_spatial_graph.SpatialSpot` silently never fires.
expect_null(
  current = get_per_sample_spatial_graph(object, "sample_a"),
  info = "sp analysis - the sp_cache getter dispatches on a SpatialSpot"
)

object <- build_spatial_graph_sp(object, .verbose = FALSE)

expect_true(
  current = inherits(
    get_per_sample_spatial_graph(object, "sample_a"),
    "CsparseMatrix"
  ),
  info = "sp analysis - a SpatialSpot caches its own graph"
)

## pipeline --------------------------------------------------------------------

pipeline <- step_spatial_graph_sp(.verbose = FALSE) %>>%
  step_morans_i_sp(.verbose = FALSE)

expect_equal(
  current = length(pipeline),
  target = 2L,
  info = "sp analysis - the spatial steps chain"
)

expect_error(
  current = validate_pipeline(pipeline, "SingleCells"),
  info = "sp analysis - the spatial steps refuse a plain SingleCells"
)

expect_silent(
  current = validate_pipeline(pipeline, "SpatialSpotSubset"),
  info = "sp analysis - the spatial steps accept a SpatialSpotSubset"
)

piped <- apply_pipeline(pipeline, SpatialSpotSubset(object, "sample_a"))

expect_equal(
  current = get_per_sample_morans_i(piped, "sample_a")$morans_i,
  target = morans$morans_i,
  info = "sp analysis - the pipeline reproduces the direct calls"
)
