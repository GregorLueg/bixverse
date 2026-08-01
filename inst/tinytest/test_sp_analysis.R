# spatial subset and analysis methods ------------------------------------------

test_temp_dir <- file.path(tempdir(), "sp_analysis")
dir.create(test_temp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir))

## fixture ---------------------------------------------------------------------

### a visium hex lattice -------------------------------------------------------

# Standard Visium geometry: `array_col` steps by two inside a row and the rows
# are offset by one, so an interior spot has exactly six neighbours. Anything
# that gets the parity wrong drops to four.
n_array_rows <- 12L
n_per_row <- 12L
n_spots <- n_array_rows * n_per_row
n_genes <- 60L

sf_json <- list(
  regist_target_img_scalef = 0.2,
  tissue_hires_scalef = 0.1,
  tissue_lowres_scalef = 0.05,
  fiducial_diameter_fullres = 100,
  spot_diameter_fullres = 40
)

grid <- data.table::rbindlist(lapply(seq_len(n_array_rows) - 1L, \(r) {
  data.table::data.table(
    array_row = r,
    array_col = seq(r %% 2L, by = 2L, length.out = n_per_row)
  )
}))
grid[, pxl_col_in_fullres := as.integer(200L + array_col * 50L)]
grid[, pxl_row_in_fullres := as.integer(200L + array_row * 87L)]

set.seed(11L)
barcodes_raster_order <- vapply(
  seq_len(n_spots),
  \(i)
    sprintf(
      "%s-1",
      paste(sample(c("A", "C", "G", "T"), 16L, TRUE), collapse = "")
    ),
  character(1)
)
stopifnot(!anyDuplicated(barcodes_raster_order))
barcodes_matrix_order <- sort(barcodes_raster_order)

positions <- data.table::data.table(
  barcode = barcodes_raster_order,
  in_tissue = 1L,
  array_row = grid$array_row,
  array_col = grid$array_col,
  pxl_row_in_fullres = grid$pxl_row_in_fullres,
  pxl_col_in_fullres = grid$pxl_col_in_fullres
)
data.table::setkey(positions, barcode)

# Counts in matrix (barcode) order, with a planted north/south gradient on gene
# 1 and pure noise on gene 2. The gradient is kept modest on purpose: Moran's I
# runs on normalised counts by default, so a gene big enough to move the library
# size drags spatial structure into every other gene through the denominator.
pos_in_matrix_order <- positions[barcodes_matrix_order]
north <- pos_in_matrix_order$array_row < n_array_rows / 2L

set.seed(22L)
counts <- matrix(
  rpois(n_genes * n_spots, lambda = 4),
  nrow = n_genes,
  ncol = n_spots
)
counts[1L, ] <- ifelse(north, 30L, 3L)
counts[2L, ] <- rpois(n_spots, lambda = 20)
gene_ids <- sprintf("ENSGSP%05i", seq_len(n_genes))

### writing it out -------------------------------------------------------------

write_png_or_stub <- function(path, width, height) {
  if (requireNamespace("png", quietly = TRUE)) {
    set.seed(33L)
    png::writePNG(
      array(runif(height * width * 3L), dim = c(height, width, 3L)),
      target = path
    )
    return(TRUE)
  }
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(137, 80, 78, 71, 13, 10, 26, 10)), con)
  writeBin(13L, con, size = 4L, endian = "big")
  writeBin(charToRaw("IHDR"), con)
  writeBin(width, con, size = 4L, endian = "big")
  writeBin(height, con, size = 4L, endian = "big")
  FALSE
}

visium_dir <- file.path(test_temp_dir, "sample_a")
matrix_dir <- file.path(visium_dir, "raw_feature_bc_matrix")
spatial_dir <- file.path(visium_dir, "spatial")
dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)

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
data.table::fwrite(
  positions,
  file.path(spatial_dir, "tissue_positions.csv"),
  col.names = TRUE
)
writeLines(
  jsonlite::toJSON(sf_json, auto_unbox = TRUE),
  file.path(spatial_dir, "scalefactors_json.json")
)
# hires covers the full-res frame (max x 1350, max y 1157) at 0.1
png_written <- write_png_or_stub(
  file.path(spatial_dir, "tissue_hires_image.png"),
  200L,
  200L
)

### loading --------------------------------------------------------------------

data_dir <- file.path(test_temp_dir, "out")
unlink(data_dir, recursive = TRUE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

object <- SpatialSpot(dir_data = data_dir)
object <- load_visium(
  object,
  visium_dir = visium_dir,
  sp_io_param = params_sp_visium_io(),
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 0L,
    min_lib_size = 0L,
    min_cells = 0L
  ),
  exp_id = "sample_a",
  .verbose = FALSE
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

expect_equal(
  current = get_spatial_coords(subset_obj, "sample_a"),
  target = get_spatial_coords(object, "sample_a"),
  info = "sp analysis - the cached coords match the parent"
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

expect_equal(
  current = morans$gene_id,
  target = get_gene_names(subset_obj),
  info = "sp analysis - the gene ids come back in var table order"
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

obs_labels <- bixverse:::.sp_obs_for_exp(subset_obj, "sample_a")
subset_obj[["band"]] <- ifelse(
  obs_labels$array_row < n_array_rows / 2L,
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

if (png_written) {
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
  info = "sp analysis - the SpCache has the sixth slot"
)

expect_equal(
  current = names(sp_cache$per_sample_nhood_enrichment$sample_a),
  target = "band",
  info = "sp analysis - enrichment nests by exp_id then label column"
)

if (png_written) {
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
