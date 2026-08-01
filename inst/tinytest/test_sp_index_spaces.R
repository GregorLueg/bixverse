# global vs local spot index spaces --------------------------------------------

# Every other spatial fixture runs one sample with every spot kept, so
# `subset_to_original` is `1:n` and the global 0-based spot index, the local
# positional index into the graph and the row of the coordinate matrix are all
# the same numbers. The index-space contract degenerates and nothing tests it.
#
# This file is the opposite case: two samples of different shape, interior
# spots dropped from both so `to_keep` is non-contiguous, and the second
# sample's spots starting well past zero in the global index space. Every
# method that mixes the two spaces runs over it.
#
# The lever is the planted structure. Gene 3 is spatially graded in sample A
# and flat noise in sample B, gene 5 the other way round. Anything that hands
# Rust a local index where a global one belongs reads sample A's counts for
# sample B and the two swap over.

source("helper_sp.R", local = TRUE)

test_temp_dir <- file.path(tempdir(), "sp_index_spaces")
dir.create(test_temp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir))

n_genes <- 40L
graded_a <- 3L
graded_b <- 5L

## fixtures --------------------------------------------------------------------

# north/south on gene 3, flat on gene 5
counts_a <- function(positions, n_genes, n_rows) {
  n_spots <- nrow(positions)
  m <- matrix(
    stats::rpois(n_genes * n_spots, lambda = 4),
    nrow = n_genes,
    ncol = n_spots
  )
  m[graded_a, ] <- ifelse(positions$array_row < n_rows / 2L, 20L, 3L)
  m[graded_b, ] <- stats::rpois(n_spots, lambda = 10)
  m
}

# east/west on gene 5, flat on gene 3
counts_b <- function(positions, n_genes, n_rows) {
  n_spots <- nrow(positions)
  m <- matrix(
    stats::rpois(n_genes * n_spots, lambda = 4),
    nrow = n_genes,
    ncol = n_spots
  )
  m[graded_a, ] <- stats::rpois(n_spots, lambda = 10)
  m[graded_b, ] <- ifelse(positions$array_col < 8L, 20L, 3L)
  m
}

fixture_a <- sp_make_visium(
  file.path(test_temp_dir, "sample_a"),
  layout = "hex",
  n_rows = 10L,
  n_cols = 10L,
  n_genes = n_genes,
  seed = 505L,
  counts_fun = counts_a
)

fixture_b <- sp_make_visium(
  file.path(test_temp_dir, "sample_b"),
  layout = "hex",
  n_rows = 6L,
  n_cols = 8L,
  n_genes = n_genes,
  seed = 606L,
  counts_fun = counts_b
)

gene_a <- fixture_a$gene_ids[graded_a]
gene_b <- fixture_a$gene_ids[graded_b]

# A two-tone hires image for sample B: everything left of the boundary close to
# black, everything right of it close to white, so the optical density of a
# tile says which side of the slide its spot sits on. Every other image fixture
# in the suite is `runif` noise, where an x/y swap in the coordinate-to-tile
# mapping gives different-but-finite values that no assertion can tell apart.
hires_scalef <- fixture_b$scale_factors$tissue_hires_scalef
spot_radius_hires <- fixture_b$scale_factors$spot_diameter_fullres *
  hires_scalef /
  2
fixture_b_x <- fixture_b$positions$pxl_col_in_fullres * hires_scalef
dark_boundary <- (min(fixture_b_x) + max(fixture_b_x)) / 2

if (sp_has_png()) {
  tone <- array(0.95, dim = c(fixture_b$image_size, fixture_b$image_size, 3L))
  tone[, seq_len(as.integer(dark_boundary)), ] <- 0.05
  png::writePNG(
    tone,
    target = file.path(
      fixture_b$visium_dir,
      "spatial",
      "tissue_hires_image.png"
    )
  )
}

## loading ---------------------------------------------------------------------

data_dir <- file.path(test_temp_dir, "out")
unlink(data_dir, recursive = TRUE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

prescan <- prescan_visium_dirs(
  dirs = c(fixture_a$visium_dir, fixture_b$visium_dir),
  exp_ids = c("sample_a", "sample_b"),
  .verbose = FALSE
)

object <- load_multi_visium(
  SpatialSpot(dir_data = data_dir),
  prescan_result = prescan,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 0L,
    min_lib_size = 0L,
    min_cells = 0L
  ),
  .verbose = FALSE
)

## a non-contiguous to_keep ----------------------------------------------------

obs_all <- get_sc_obs(object, filtered = TRUE)
data.table::setorderv(obs_all, "cell_idx")
idx_a <- obs_all[obs_all[["exp_id"]] == "sample_a", ][["cell_idx"]]
idx_b <- obs_all[obs_all[["exp_id"]] == "sample_b", ][["cell_idx"]]

expect_true(
  current = min(idx_b) > max(idx_a),
  info = "sp index - the second sample starts past the first in cell_idx"
)

# interior spots, so both samples end up with holes in the middle of their
# global index range rather than a truncated tail
dropped <- c(idx_a[c(17L, 41L, 68L)], idx_b[c(11L, 30L)])
kept <- setdiff(obs_all[["cell_idx"]], dropped)
object <- set_cells_to_keep(object, kept)

n_kept_a <- length(idx_a) - 3L
n_kept_b <- length(idx_b) - 2L

expect_equal(
  current = length(get_spot_indices_for_exp(object, "sample_b")),
  target = n_kept_b,
  info = "sp index - the dropped spots leave the filtered spot set"
)

subset_a <- SpatialSpotSubset(object, exp_id = "sample_a")
subset_b <- SpatialSpotSubset(object, exp_id = "sample_b")
sto_a <- S7::prop(subset_a, "subset_to_original")
sto_b <- S7::prop(subset_b, "subset_to_original")

expect_true(
  current = any(diff(sto_a) > 1L) && any(diff(sto_b) > 1L),
  info = "sp index - subset_to_original has holes in both samples"
)

expect_false(
  current = identical(sto_b, seq_along(sto_b)),
  info = "sp index - and the second sample's is not 1:n"
)

expect_true(
  current = min(sto_b) > length(idx_a),
  info = "sp index - the second sample sits at a non-zero global offset"
)

expect_equal(
  current = as.integer(get_spot_indices_for_exp(subset_b, "sample_b")),
  target = sto_b - 1L,
  info = "sp index - the subset reports the global 0-based indices"
)

# coordinates checked against the fixture's own positions file by barcode,
# rather than against another read of the table they came from
obs_b <- get_sc_obs(subset_b)
data.table::setorderv(obs_b, "cell_idx")
truth_b <- fixture_b$positions[match(
  obs_b[["cell_id"]],
  fixture_b$positions[["barcode"]]
)]
coords_b <- S7::prop(subset_b, "coords")

expect_equal(
  current = as.numeric(coords_b[, "x"]),
  target = as.numeric(truth_b[["pxl_col_in_fullres"]]),
  info = "sp index - the subset coords match the positions file by barcode"
)

expect_equal(
  current = as.numeric(coords_b[, "y"]),
  target = as.numeric(truth_b[["pxl_row_in_fullres"]]),
  info = "sp index - on both axes"
)

## the graph -------------------------------------------------------------------

object <- build_spatial_graph_sp(object, .verbose = FALSE)

expect_equal(
  current = dim(get_per_sample_spatial_graph(object, "sample_a")),
  target = rep(n_kept_a, 2L),
  info = "sp index - the graph is built over the kept spots, not every spot"
)

expect_equal(
  current = dim(get_per_sample_spatial_graph(object, "sample_b")),
  target = rep(n_kept_b, 2L),
  info = "sp index - and per sample, not object-wide"
)

# dropping an interior spot costs each of its neighbours an edge, so the graph
# is not the full lattice with rows deleted
degree_b <- lengths(
  bixverse:::.sp_sparse_to_adjacency(
    get_per_sample_spatial_graph(object, "sample_b")
  )$indices
)
full_degree_b <- lengths(sp_lattice_neighbours(
  fixture_b$positions[["array_row"]],
  fixture_b$positions[["array_col"]],
  "hex"
))

expect_true(
  current = max(degree_b) == 6L &&
    sum(degree_b == 6L) < sum(full_degree_b == 6L),
  info = "sp index - dropping interior spots costs their neighbours an edge"
)

## a graph built over a different spot set -------------------------------------

# N5. The graph lives in the local index space and the matrix records nothing
# about which spots that was. Swapping `to_keep` for a different set of the
# same size used to pair expression, labels and adjacency against three
# different spots and return finite, plausible numbers.
swapped_keep <- sort(c(setdiff(kept, idx_a[20L]), dropped[1L]))
stopifnot(
  "The swap changed the spot count" = length(swapped_keep) == length(kept),
  "The swap changed nothing" = !identical(swapped_keep, kept)
)
# the DuckDB behind the object is an R6 connection, so this rewrites `to_keep`
# for every handle on the same store, `object` included
swapped <- set_cells_to_keep(object, swapped_keep)

expect_error(
  current = morans_i_sp(swapped, exp_id = "sample_a", .verbose = FALSE),
  info = "sp index - Moran's I refuses a graph built over other spots"
)

expect_error(
  current = nhood_enrichment_sp(
    swapped,
    label_col = "array_row",
    exp_id = "sample_a",
    .verbose = FALSE
  ),
  info = "sp index - and so does the neighbourhood enrichment"
)

# sample B's spot set never moved, so its graph is still the right one
expect_silent(
  current = morans_i_sp(swapped, exp_id = "sample_b", .verbose = FALSE),
  info = "sp index - the untouched sample is still fine"
)

object <- set_cells_to_keep(object, kept)

## moran's i -------------------------------------------------------------------

object <- morans_i_sp(object, .verbose = FALSE)
morans_a <- get_per_sample_morans_i(object, "sample_a")
morans_b <- get_per_sample_morans_i(object, "sample_b")

expect_equal(
  current = morans_a$gene_id[which.max(morans_a$morans_i)],
  target = gene_a,
  info = "sp index - sample A's top gene is the one graded in sample A"
)

expect_equal(
  current = morans_b$gene_id[which.max(morans_b$morans_i)],
  target = gene_b,
  info = "sp index - and sample B's is the one graded in sample B"
)

# a local index where a global one belongs reads the other sample's counts, and
# those rows are in barcode order, so nothing scores at all
expect_true(
  current = morans_a[gene_id == gene_a]$morans_i >
    3 * morans_a[gene_id == gene_b]$morans_i &&
    morans_b[gene_id == gene_b]$morans_i >
      3 * morans_b[gene_id == gene_a]$morans_i,
  info = "sp index - each graded gene clears the other sample's by a margin"
)

expect_true(
  current = morans_a[gene_id == gene_a]$fdr < 0.01 &&
    morans_b[gene_id == gene_b]$fdr < 0.01,
  info = "sp index - both planted gradients are significant in their own sample"
)

# `.sp_resolve_genes` with a character selection in reverse var-table order.
# The indices go to Rust as on-disk gene positions and the labels come back by
# position, so a permuted var read would separate the two.
gene_names <- get_gene_names(object)
wanted <- rev(gene_names[c(graded_a, graded_b, 12L)])
selected <- get_per_sample_morans_i(
  morans_i_sp(object, exp_id = "sample_a", genes = wanted, .verbose = FALSE),
  "sample_a"
)

expect_equal(
  current = selected$gene_id,
  target = wanted,
  info = "sp index - a gene selection comes back in the order it was asked for"
)

expect_equal(
  current = selected$gene_idx,
  target = as.integer(match(wanted, gene_names) - 1L),
  info = "sp index - and the gene_idx out of Rust matches the label beside it"
)

expect_equal(
  current = selected$morans_i,
  target = morans_a[match(wanted, morans_a$gene_id)]$morans_i,
  info = "sp index - the selected genes score what the full run gave them"
)

expect_warning(
  current = morans_i_sp(
    object,
    exp_id = "sample_a",
    genes = c(gene_a, "NOT_A_GENE"),
    .verbose = FALSE
  ),
  info = "sp index - an unmatched gene name warns rather than going quiet"
)

## sparkx ----------------------------------------------------------------------

object <- sparkx_sp(object, .verbose = FALSE)
sparkx_a <- get_per_sample_sparkx(object, "sample_a")
sparkx_b <- get_per_sample_sparkx(object, "sample_b")

expect_true(
  current = sparkx_a[gene_id == gene_a]$fdr < 0.01 &&
    sparkx_a[gene_id == gene_b]$fdr > 0.01,
  info = "sp index - SPARK-X picks sample A's graded gene, not sample B's"
)

expect_true(
  current = sparkx_b[gene_id == gene_b]$fdr < 0.01 &&
    sparkx_b[gene_id == gene_a]$fdr > 0.01,
  info = "sp index - and the other way round in sample B"
)

# an explicit kernel bank is validated by checkSpSparkx and never otherwise
# passed to the method, so nothing checks it reaches Rust
one_kernel <- get_per_sample_sparkx(
  sparkx_sp(
    object,
    exp_id = "sample_b",
    sparkx_params = params_sp_sparkx(
      kernels = list(list(kernel = "gaussian", bandwidth = 200))
    ),
    .verbose = FALSE
  ),
  "sample_b"
)

expect_equal(
  current = dim(attr(one_kernel, "stat_per_kernel")),
  target = c(n_genes, 1L),
  info = "sp index - an explicit kernel bank reaches Rust"
)

## neighbourhood enrichment ----------------------------------------------------

# Three levels in sample A, two in sample B, in one obs column. The column is
# written against ascending global `cell_idx` across both samples, so a method
# that pulled its labels globally rather than per sample gets the wrong K as
# well as the wrong labels.
obs_kept <- get_sc_obs(object, filtered = TRUE)
data.table::setorderv(obs_kept, "cell_idx")
band <- ifelse(
  obs_kept[["exp_id"]] == "sample_a",
  c("low", "mid", "high")[cut(
    obs_kept[["array_row"]],
    breaks = c(-1L, 2L, 6L, 100L),
    labels = FALSE
  )],
  ifelse(obs_kept[["array_col"]] < 8L, "west", "east")
)
object[["band"]] <- band

object <- nhood_enrichment_sp(
  object,
  label_col = "band",
  nhood_params = params_sp_nhood(n_perm = 200L),
  seed = 7L,
  .verbose = FALSE
)
nhood_a <- get_per_sample_nhood_enrichment(object, "sample_a", "band")
nhood_b <- get_per_sample_nhood_enrichment(object, "sample_b", "band")

expect_equal(
  current = dim(nhood_a$z_scores),
  target = c(3L, 3L),
  info = "sp index - sample A gets its own three levels"
)

expect_equal(
  current = nhood_b$label_levels,
  target = c("east", "west"),
  info = "sp index - sample B gets only the two levels it carries"
)

expect_true(
  current = all(diag(nhood_a$z_scores) > 0),
  info = "sp index - each of the three bands self-associates"
)

# the outer two bands never touch, the middle one touches both, so a K x K
# indexing slip shows up here and nowhere in a two-level test
expect_true(
  current = nhood_a$z_scores["low", "high"] < nhood_a$z_scores["low", "mid"] &&
    nhood_a$z_scores["low", "high"] < nhood_a$z_scores["mid", "high"],
  info = "sp index - the two outer bands avoid each other most"
)

expect_true(
  current = isTRUE(all.equal(nhood_a$z_scores, t(nhood_a$z_scores))),
  info = "sp index - symmetrise = TRUE gives a symmetric Z matrix"
)

# `symmetrise` is a no-op by construction. The crate credits every edge to both
# off-diagonal cells, so the count matrix is symmetric whatever the graph is,
# and averaging Z with its transpose changes nothing. Pinned here rather than
# left as a parameter nobody has ever passed: it goes red the day the count
# rule becomes directed and the documentation does not follow.
knn_object <- build_spatial_graph_sp(
  object,
  exp_id = "sample_a",
  graph_params = params_sp_graph("knn", k = 3L),
  .verbose = FALSE
)
knn_nhood <- function(symmetrise) {
  get_per_sample_nhood_enrichment(
    nhood_enrichment_sp(
      knn_object,
      label_col = "band",
      exp_id = "sample_a",
      nhood_params = params_sp_nhood(n_perm = 200L, symmetrise = symmetrise),
      seed = 7L,
      .verbose = FALSE
    ),
    "sample_a",
    "band"
  )
}

expect_true(
  current = isTRUE(all.equal(
    knn_nhood(TRUE)$z_scores,
    knn_nhood(FALSE)$z_scores
  )),
  info = "sp index - symmetrise is a no-op on a symmetric count matrix"
)

expect_true(
  current = isTRUE(all.equal(
    knn_nhood(FALSE)$observed,
    t(knn_nhood(FALSE)$observed)
  )),
  info = "sp index - a directed kNN graph still gives symmetric edge counts"
)

## image features --------------------------------------------------------------

if (sp_image_backend_ok()) {
  object <- image_features_sp(object, resolution = "hires", .verbose = FALSE)
  feat_a <- get_per_sample_image_features(object, "sample_a", "hires")
  feat_b <- get_per_sample_image_features(object, "sample_b", "hires")

  expect_true(
    current = all(feat_b$cell_idx %in% sto_b),
    info = "sp index - the feature cell_idx sit in the parent index space"
  )

  # the assertion a `1:n` fixture cannot make: swap the `spot_cell_idx` lookup
  # in new_sp_image_features_res() for the positional index and this goes red
  expect_true(
    current = min(feat_b$cell_idx) > length(idx_a),
    info = "sp index - and carry the second sample's global offset"
  )

  expect_equal(
    current = length(feat_a$cell_idx) + feat_a$params$n_dropped,
    target = n_kept_a,
    info = "sp index - the tile run covers the kept spots of its own sample"
  )

  # tie each feature row back to the pixels of its own spot
  spot_barcode <- obs_b[["cell_id"]][
    match(feat_b$cell_idx, obs_b[["cell_idx"]])
  ]
  spot_x <- fixture_b$positions[["pxl_col_in_fullres"]][
    match(spot_barcode, fixture_b$positions[["barcode"]])
  ] *
    hires_scalef
  od <- feat_b$values[, "mean_od"]
  dark_side <- spot_x < dark_boundary - spot_radius_hires - 1
  light_side <- spot_x > dark_boundary + spot_radius_hires + 1

  expect_true(
    current = sum(dark_side) > 5L && sum(light_side) > 5L,
    info = "sp index - the two-tone image splits the spots into two groups"
  )

  expect_true(
    current = min(od[dark_side]) > max(od[light_side]),
    info = "sp index - every dark-side spot is darker than every light one"
  )
}
