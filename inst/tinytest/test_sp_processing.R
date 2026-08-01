# spatially variable genes -----------------------------------------------------

# Moran's I and SPARK-X. test_sp_analysis.R already checks that both run end to
# end and that the planted gradient comes out on top; this file is about the
# numbers themselves.
#
# Moran's I is pinned against two independent references:
#
#   1. `spdep::moran.test(randomisation = FALSE)`. The Rust implements the
#      Cliff-Ord normality variance specifically to match it. Live when spdep
#      is installed, hard-coded values otherwise.
#   2. Space Ranger's own `spatial_enrichment.csv` on the local example data.
#      Rank correlation only, because Space Ranger builds its own neighbour
#      graph.

source("helper_sp.R", local = TRUE)

test_temp_dir <- file.path(tempdir(), "sp_processing")
dir.create(test_temp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir))

## fixture ---------------------------------------------------------------------

fixture <- sp_make_visium(
  file.path(test_temp_dir, "sample_a"),
  layout = "hex",
  n_rows = 10L,
  n_cols = 10L,
  n_genes = 40L,
  seed = 101L
)

object <- sp_load_visium_fixture(fixture, file.path(test_temp_dir, "out"))
object <- build_spatial_graph_sp(object, .verbose = FALSE)

n_spots <- fixture$n_spots
n_genes <- fixture$n_genes
gene_ids <- fixture$gene_ids

graph <- get_per_sample_spatial_graph(object, "sample_a")
# what the statistics actually see: Rust folds (i, j) and (j, i) together and
# scatters the sum back into both rows, so the graph is always W + t(W)
sym_graph <- graph + Matrix::t(graph)

raw_counts <- get_sc_counts(
  object,
  assay = "raw",
  return_format = "gene",
  .verbose = FALSE
)

expect_equal(
  current = dim(raw_counts),
  target = c(n_spots, n_genes),
  info = "sp processing - the fixture round trips at the expected shape"
)

## the numerical pin -----------------------------------------------------------

# Raw counts, not normalised ones. Moran's I runs on normalised counts by
# default, which couples every gene through the library size and makes an exact
# external comparison impossible to set up honestly.
object <- morans_i_sp(
  object,
  svg_params = params_sp_svg(assay = "raw"),
  .verbose = FALSE
)
morans_raw <- get_per_sample_morans_i(object, "sample_a")

### closed form Cliff-Ord ------------------------------------------------------

# An R re-derivation straight from the definition. Independent of spdep, so the
# variance stays pinned even on a machine without it.

sp_moran_reference <- function(x, w) {
  n <- length(x)
  s0 <- sum(w)
  s1 <- 0.5 * sum((w + Matrix::t(w))^2)
  s2 <- sum((Matrix::rowSums(w) + Matrix::colSums(w))^2)

  z <- x - mean(x)
  i <- (n / s0) * sum(z * as.numeric(w %*% z)) / sum(z^2)
  e_i <- -1 / (n - 1)
  var_i <- (n^2 * s1 - n * s2 + 3 * s0^2) / (s0^2 * (n^2 - 1)) - e_i^2

  list(i = i, e_i = e_i, var_i = var_i, z = (i - e_i) / sqrt(var_i))
}

ref_all <- lapply(seq_len(n_genes), \(g) {
  sp_moran_reference(as.numeric(raw_counts[, g]), sym_graph)
})

expect_equal(
  current = attr(morans_raw, "e_i"),
  target = -1 / (n_spots - 1),
  info = "sp processing - E[I] is -1 / (N - 1)"
)

expect_equal(
  current = attr(morans_raw, "var_i"),
  target = ref_all[[1L]]$var_i,
  tolerance = 1e-10,
  info = "sp processing - Var[I] matches the Cliff-Ord normality closed form"
)

expect_equal(
  current = morans_raw$morans_i,
  target = vapply(ref_all, \(r) r$i, numeric(1)),
  tolerance = 1e-6,
  info = "sp processing - every gene's I matches the closed form"
)

expect_equal(
  current = morans_raw$z,
  target = vapply(ref_all, \(r) r$z, numeric(1)),
  tolerance = 1e-5,
  info = "sp processing - the Z scores match the closed form"
)

expect_equal(
  current = morans_raw$pval,
  target = 2 * stats::pnorm(-abs(vapply(ref_all, \(r) r$z, numeric(1)))),
  tolerance = 1e-6,
  info = "sp processing - the p-value is two-sided normal"
)

### spdep ----------------------------------------------------------------------

# Generated with spdep 1.4.2 on this exact fixture (hex lattice, 10 x 10, seed
# 101, raw assay, binary weights) via
#   spdep::moran.test(x, spdep::mat2listw(W + t(W), style = "B"),
#                     randomisation = FALSE, alternative = "two.sided")
# Kept as literals so the pin survives on a machine without spdep, and checked
# live against spdep whenever it is installed.
pin_gene_ids <- c(
  "ENSGSP00001",
  "ENSGSP00003",
  "ENSGSP00007",
  "ENSGSP00011",
  "ENSGSP00019",
  "ENSGSP00028"
)
spdep_i <- c(
  0.85440613026819923,
  -0.07141824226521157,
  0.01072819836346023,
  0.01753277354708851,
  -0.00207180360437064,
  -0.01232881908880543
)
spdep_z <- c(
  14.3880619704169721,
  -1.0205076338134056,
  0.3466621941995088,
  0.4599112870671428,
  0.1336307304501677,
  -0.0370776044266815
)
spdep_e_i <- -0.0101010101010101
spdep_var_i <- 0.00361021087382269

pinned <- morans_raw[match(pin_gene_ids, morans_raw$gene_id), ]

expect_equal(
  current = pinned$morans_i,
  target = spdep_i,
  tolerance = 1e-6,
  info = "sp processing - Moran's I matches the stored spdep values"
)

expect_equal(
  current = pinned$z,
  target = spdep_z,
  tolerance = 1e-5,
  info = "sp processing - the Z scores match the stored spdep values"
)

expect_equal(
  current = attr(morans_raw, "e_i"),
  target = spdep_e_i,
  tolerance = 1e-10,
  info = "sp processing - E[I] matches the stored spdep value"
)

expect_equal(
  current = attr(morans_raw, "var_i"),
  target = spdep_var_i,
  tolerance = 1e-10,
  info = "sp processing - Var[I] matches the stored spdep value"
)

if (requireNamespace("spdep", quietly = TRUE)) {
  listw <- spdep::mat2listw(sym_graph, style = "B", zero.policy = TRUE)

  live <- lapply(pin_gene_ids, \(gene) {
    spdep::moran.test(
      as.numeric(raw_counts[, gene]),
      listw,
      randomisation = FALSE,
      alternative = "two.sided",
      zero.policy = TRUE
    )
  })

  expect_equal(
    current = pinned$morans_i,
    target = vapply(live, \(m) unname(m$estimate[1L]), numeric(1)),
    tolerance = 1e-6,
    info = "sp processing - Moran's I matches a live spdep::moran.test"
  )

  expect_equal(
    current = attr(morans_raw, "var_i"),
    target = unname(live[[1L]]$estimate[3L]),
    tolerance = 1e-10,
    info = "sp processing - Var[I] matches spdep with randomisation = FALSE"
  )

  expect_equal(
    current = pinned$z,
    target = vapply(live, \(m) unname(m$statistic), numeric(1)),
    tolerance = 1e-5,
    info = "sp processing - the Z scores match spdep"
  )

  expect_equal(
    current = pinned$pval,
    target = vapply(live, \(m) m$p.value, numeric(1)),
    tolerance = 1e-5,
    info = "sp processing - the two-sided p-values match spdep"
  )
}

## moran's i behaviour ---------------------------------------------------------

expect_equal(
  current = morans_raw$gene_id,
  target = get_gene_names(object),
  info = "sp processing - results come back in var table order"
)

expect_true(
  current = all(morans_raw$morans_i <= 1.0001),
  info = "sp processing - I stays at or below one"
)

expect_true(
  current = morans_raw[gene_id == gene_ids[1L]]$morans_i > 0.5,
  info = "sp processing - the planted north/south gradient scores high"
)

### streaming ------------------------------------------------------------------

streamed <- get_per_sample_morans_i(
  morans_i_sp(
    object,
    svg_params = params_sp_svg(assay = "raw"),
    streaming = FALSE,
    .verbose = FALSE
  ),
  "sample_a"
)

expect_equal(
  current = streamed$morans_i,
  target = morans_raw$morans_i,
  info = "sp processing - streaming and in-memory give the same I"
)

### assay choice ---------------------------------------------------------------

morans_norm <- get_per_sample_morans_i(
  morans_i_sp(object, .verbose = FALSE),
  "sample_a"
)

expect_true(
  current = !isTRUE(all.equal(morans_norm$morans_i, morans_raw$morans_i)),
  info = "sp processing - the norm assay is not the raw one"
)

expect_equal(
  current = attr(morans_norm, "var_i"),
  target = attr(morans_raw, "var_i"),
  info = "sp processing - the null depends on the graph, not on the assay"
)

### gene selection -------------------------------------------------------------

subset_by_name <- get_per_sample_morans_i(
  morans_i_sp(
    object,
    genes = gene_ids[c(5L, 1L, 9L)],
    svg_params = params_sp_svg(assay = "raw"),
    .verbose = FALSE
  ),
  "sample_a"
)

expect_equal(
  current = nrow(subset_by_name),
  target = 3L,
  info = "sp processing - a gene selection tests only those genes"
)

expect_equal(
  current = subset_by_name[order(gene_id)]$morans_i,
  target = morans_raw[gene_id %in% gene_ids[c(1L, 5L, 9L)]][
    order(gene_id)
  ]$morans_i,
  info = "sp processing - a subset reproduces the full-run values"
)

subset_by_index <- get_per_sample_morans_i(
  morans_i_sp(
    object,
    genes = c(5L, 1L, 9L),
    svg_params = params_sp_svg(assay = "raw"),
    .verbose = FALSE
  ),
  "sample_a"
)

expect_equal(
  current = sort(subset_by_index$gene_id),
  target = sort(gene_ids[c(1L, 5L, 9L)]),
  info = "sp processing - integer gene indices are 1-based into the var table"
)

expect_error(
  current = morans_i_sp(
    object,
    genes = "not_a_gene",
    .verbose = FALSE
  ),
  info = "sp processing - an unknown gene identifier errors"
)

### constant genes -------------------------------------------------------------

# A gene with zero variance has no defined I. Rust emits NaN and the Z and
# p-value follow.
#
# Known defect, not asserted here: the BH step maps those NaN p-values to an FDR
# of 0, so constant genes come out "significant" and lead the `print` method's
# top-gene table. Filter on `is.finite(morans_i)` before reading the FDR column.
flat_fixture <- sp_make_visium(
  file.path(test_temp_dir, "flat"),
  layout = "hex",
  n_rows = 6L,
  n_cols = 6L,
  n_genes = 5L,
  seed = 55L,
  counts_fun = function(positions, n_genes, n_rows) {
    counts <- matrix(
      3L,
      nrow = n_genes,
      ncol = nrow(positions)
    )
    counts[1L, ] <- ifelse(positions$array_row < n_rows / 2L, 12L, 2L)
    counts
  }
)
flat_object <- sp_load_visium_fixture(
  flat_fixture,
  file.path(test_temp_dir, "out_flat")
)
flat_object <- build_spatial_graph_sp(flat_object, .verbose = FALSE)
flat_object <- morans_i_sp(
  flat_object,
  svg_params = params_sp_svg(assay = "raw"),
  .verbose = FALSE
)
flat_res <- get_per_sample_morans_i(flat_object, "sample_a")

expect_true(
  current = all(is.nan(
    flat_res[gene_id != flat_fixture$gene_ids[1L]]$morans_i
  )),
  info = "sp processing - a constant gene has no defined I"
)

expect_true(
  current = all(is.nan(flat_res[gene_id != flat_fixture$gene_ids[1L]]$z)),
  info = "sp processing - and no Z score either"
)

expect_true(
  current = all(
    is.nan(flat_res[gene_id != flat_fixture$gene_ids[1L]]$pval)
  ),
  info = "sp processing - and no p-value"
)

expect_true(
  current = is.finite(flat_res[gene_id == flat_fixture$gene_ids[1L]]$morans_i),
  info = "sp processing - the one varying gene still gets a number"
)

### the result class on derived tables -----------------------------------------

# `SpMoransRes` sits in front of a data.table, and `[` and merge() carry the
# class onto whatever comes out. A summary built off a result is then still an
# SpMoransRes without the columns the print method wants, so the print method
# has to survive that rather than error in somebody's face.

derived <- merge(
  morans_raw[, c("gene_id", "fdr"), with = FALSE],
  data.table::data.table(
    gene_id = get_per_sample_sparkx(
      sparkx_sp(object, .verbose = FALSE),
      "sample_a"
    )$gene_id,
    sparkx_fdr = 1
  ),
  by = "gene_id"
)

expect_true(
  current = inherits(derived, "SpMoransRes"),
  info = "sp processing - merge carries the result class onto the output"
)

expect_silent(
  current = print(derived[, list(n = .N)]),
  info = "sp processing - printing a derived table without the columns is fine"
)

expect_equal(
  current = class(bixverse:::.sp_strip_result_class(morans_raw)),
  target = c("data.table", "data.frame"),
  info = "sp processing - the result class can be dropped back to a data.table"
)

expect_equal(
  current = class(morans_raw)[1L],
  target = "SpMoransRes",
  info = "sp processing - and dropping it does not touch the original"
)

## space ranger ----------------------------------------------------------------

# Correlation rather than equality: Space Ranger builds its own neighbour graph,
# so this pins the statistic only once the graph roughly agrees. Anything below
# about 0.9 means the two have diverged for a real reason.
#
# Needs the local example data, so it only runs at home.
example_dir <- sp_example_visium_dir()

if (at_home() && !is.null(example_dir)) {
  sr_path <- file.path(example_dir, "spatial", "spatial_enrichment.csv")

  if (file.exists(sr_path)) {
    sr_out <- file.path(test_temp_dir, "space_ranger")
    unlink(sr_out, recursive = TRUE)
    dir.create(sr_out, recursive = TRUE, showWarnings = FALSE)

    sr_object <- load_visium(
      SpatialSpot(dir_data = sr_out),
      visium_dir = example_dir,
      sc_qc_param = params_sc_min_quality(
        min_unique_genes = 0L,
        min_lib_size = 0L,
        min_cells = 0L
      ),
      exp_id = "example",
      copy_images = FALSE,
      .verbose = FALSE
    )
    sr_object <- build_spatial_graph_sp(sr_object, .verbose = FALSE)
    sr_object <- morans_i_sp(sr_object, .verbose = FALSE)
    ours <- get_per_sample_morans_i(sr_object, "example")

    theirs <- data.table::fread(sr_path)
    data.table::setnames(theirs, c("Feature ID", "I"), c("gene_id", "sr_i"))

    joined <- merge(
      data.table::data.table(
        gene_id = ours$gene_id,
        bx_i = ours$morans_i
      ),
      theirs[, c("gene_id", "sr_i"), with = FALSE],
      by = "gene_id"
    )
    joined <- joined[is.finite(bx_i) & is.finite(sr_i)]

    expect_true(
      current = nrow(joined) > 15000L,
      info = "sp processing - the Space Ranger gene list is covered"
    )

    expect_true(
      current = stats::cor(
        joined$bx_i,
        joined$sr_i,
        method = "spearman"
      ) >
        0.9,
      info = "sp processing - Moran's I tracks Space Ranger's own numbers"
    )
  }
}

## sparkx ----------------------------------------------------------------------

object <- sparkx_sp(object, .verbose = FALSE)
sparkx <- get_per_sample_sparkx(object, "sample_a")

expect_equal(
  current = length(attr(sparkx, "kernels")),
  target = 10L,
  info = "sp processing - the default kernel bank is five gaussian, five cosine"
)

expect_equal(
  current = rownames(attr(sparkx, "stat_per_kernel")),
  target = sparkx$gene_id,
  info = "sp processing - the per-kernel matrix is named by gene"
)

expect_equal(
  current = colnames(attr(sparkx, "pval_per_kernel")),
  target = attr(sparkx, "kernels"),
  info = "sp processing - and by kernel"
)

per_kernel_p <- attr(sparkx, "pval_per_kernel")

expect_true(
  current = all(per_kernel_p >= 0 & per_kernel_p <= 1, na.rm = TRUE),
  info = "sp processing - per-kernel p-values sit in the unit interval"
)

expect_true(
  current = all(sparkx$pval >= 0 & sparkx$pval <= 1, na.rm = TRUE),
  info = "sp processing - the Cauchy-combined p-value does too"
)

expect_true(
  current = sparkx[gene_id == gene_ids[1L]]$pval <
    min(sparkx[gene_id != gene_ids[1L]]$pval),
  info = "sp processing - the planted gradient is the strongest SPARK-X hit"
)

### determinism ----------------------------------------------------------------

expect_equal(
  current = get_per_sample_sparkx(
    sparkx_sp(object, seed = 42L, .verbose = FALSE),
    "sample_a"
  )$pval,
  target = sparkx$pval,
  info = "sp processing - SPARK-X is reproducible at a fixed seed"
)

expect_equal(
  current = get_per_sample_sparkx(
    sparkx_sp(object, streaming = FALSE, .verbose = FALSE),
    "sample_a"
  )$pval,
  target = sparkx$pval,
  info = "sp processing - streaming does not change the SPARK-X result"
)

### an explicit kernel bank ----------------------------------------------------

one_kernel <- get_per_sample_sparkx(
  sparkx_sp(
    object,
    sparkx_params = params_sp_sparkx(
      kernels = list(list(kernel = "gaussian", bandwidth = 100))
    ),
    .verbose = FALSE
  ),
  "sample_a"
)

expect_equal(
  current = length(attr(one_kernel, "kernels")),
  target = 1L,
  info = "sp processing - an explicit kernel bank replaces the default one"
)

expect_equal(
  current = dim(attr(one_kernel, "stat_per_kernel")),
  target = c(n_genes, 1L),
  info = "sp processing - the per-kernel matrix follows the bank size"
)

expect_true(
  current = one_kernel[gene_id == gene_ids[1L]]$fdr < 0.01,
  info = "sp processing - a single gaussian kernel still finds the gradient"
)

### sparkx needs no graph ------------------------------------------------------

# Unlike Moran's I, SPARK-X works straight off the coordinates.
graphless <- SpatialSpotSubset(
  sp_load_visium_fixture(fixture, file.path(test_temp_dir, "out_nograph")),
  exp_id = "sample_a"
)

expect_silent(
  current = sparkx_sp(graphless, .verbose = FALSE),
  info = "sp processing - SPARK-X runs without a spatial graph"
)

expect_error(
  current = morans_i_sp(graphless, .verbose = FALSE),
  info = "sp processing - Moran's I does not"
)
