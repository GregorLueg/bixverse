# regression tests ------------------------------------------------------------
#
# Bugs found during the documentation sweep. Every one of these was reachable
# from a plain call with documented arguments and none had coverage, which is
# exactly why they survived. The pattern throughout: the suite passed arguments
# explicitly, so the documented defaults were never exercised.

source("helper_sc.R", local = TRUE)

set.seed(123L)

# rs_range_norm ----------------------------------------------------------------

## inverted output -------------------------------------------------------------

# array_max_min() returns (max, min); rs_range_norm() destructured it as
# (min, max), so the scale came out negative and the normalisation ran
# backwards: the largest input got the smallest output.

norm_res <- rs_range_norm(c(1, 5, 10, 20), max_val = 1, min_val = 0.05)

expect_equal(
  current = norm_res[which.max(c(1, 5, 10, 20))],
  target = 1,
  info = "rs_range_norm maps the largest input to max_val"
)

expect_equal(
  current = norm_res[which.min(c(1, 5, 10, 20))],
  target = 0.05,
  info = "rs_range_norm maps the smallest input to min_val"
)

expect_true(
  current = !is.unsorted(norm_res),
  info = "rs_range_norm is monotonically increasing in its input"
)

## degenerate inputs -----------------------------------------------------------

# array_max_min() indexes arr[0] with no empty guard, so an empty vector
# panicked out of Rust with "index out of bounds: the len is 0 but the index
# is 0".

expect_equal(
  current = length(rs_range_norm(numeric(0), max_val = 1, min_val = 0.05)),
  target = 0L,
  info = "rs_range_norm returns empty for empty input instead of panicking"
)

expect_equal(
  current = rs_range_norm(c(3, 3, 3), max_val = 1, min_val = 0.05),
  target = c(1, 1, 1),
  info = "rs_range_norm handles a constant vector without dividing by zero"
)

# get_diffcor_graph ------------------------------------------------------------

# The Rust panic above surfaced here whenever no gene pair survived the FDR and
# correlation thresholds, because the empty delta_cor vector reached
# rs_range_norm().

diffcor_obj <- local({
  sig <- synthetic_signal_matrix()
  mat <- t(sig$mat)
  target <- mat[sig$group %in% c("group1", "group2"), ]
  background <- mat[sig$group == "group3", ]
  obj <- BulkCoExp(
    target,
    data.table::data.table(sample_id = rownames(target))
  )
  obj <- preprocess_bulk_coexp(obj, hvg = 0.1, .verbose = FALSE)
  diffcor_module_processing(obj, background, .verbose = FALSE)
})

expect_silent(
  current = get_diffcor_graph(
    diffcor_obj,
    fdr_threshold = 0,
    min_cor = 1,
    .verbose = FALSE
  ),
  info = "get_diffcor_graph survives thresholds that admit no pair"
)

# documented defaults ----------------------------------------------------------

# These arguments declare a c(...) default but the method never called
# match.arg(), so the assertChoice fired on the whole vector. Calling with the
# documented default is the entire test.

bulk_obj <- local({
  mat <- t(synthetic_signal_matrix()$mat)
  obj <- BulkCoExp(mat, data.table::data.table(sample_id = rownames(mat)))
  preprocess_bulk_coexp(obj, hvg = 0.3, .verbose = FALSE)
})

expect_silent(
  current = cor_module_processing(bulk_obj, .verbose = FALSE),
  info = "cor_module_processing accepts its documented cor_method default"
)

expect_silent(
  current = cor_module_check_epsilon(
    cor_module_processing(bulk_obj, .verbose = FALSE),
    .verbose = FALSE
  ),
  info = "cor_module_check_epsilon accepts its documented rbf_func default"
)

# cor_module_tom ---------------------------------------------------------------

# Passed shift = 1L into upper_triangular_sym_mat$new(), which asserts logical.
# Every call errored.

expect_silent(
  current = cor_module_tom(
    cor_module_processing(bulk_obj, .verbose = FALSE),
    .verbose = FALSE
  ),
  info = "cor_module_tom passes a logical shift to the triangular matrix"
)

# contrastive PCA --------------------------------------------------------------

cpca_data <- synthetic_c_pca_data()

cpca_obj <- BulkCoExp(
  cpca_data$target,
  data.table::data.table(sample_id = rownames(cpca_data$target))
)

## background matrix scoping ---------------------------------------------------

# The body did `background_mat <- background_mat`. The right hand side is not a
# local (the argument is background_matrix), so R climbed to globalenv(). It
# only worked in the cpca vignette because that vignette happens to define a
# global of that name. Run it in a clean environment with a decoy global set to
# something wrong: if the function reaches out, the result changes.

background_mat <- matrix(0, nrow = 2L, ncol = 2L)

cpca_processed <- contrastive_pca_processing(
  cpca_obj,
  cpca_data$background,
  .verbose = FALSE
)

expect_equal(
  current = dim(cpca_processed@processed_data[["background_mat"]]),
  target = dim(cpca_data$background),
  info = "contrastive_pca_processing uses its argument, not a global"
)

rm(background_mat)

## no_pcs is honoured ----------------------------------------------------------

# colnames came from sprintf("cPC_%i", seq(1:10)), hardcoding ten components,
# so any other no_pcs died in colnames<-.

for (n in c(3L, 5L, 12L)) {
  expect_silent(
    current = contrastive_pca(cpca_processed, alpha = 1.0, no_pcs = n),
    info = sprintf("contrastive_pca works with no_pcs = %d", n)
  )
}

# generate_personalisation_vec -------------------------------------------------

# Two bugs: `diffusion_vec < diffusion_vec / sum(...)` was a comparison, not an
# assignment, so the vector was never normalised; and qassert(node_weights,
# "N1") capped the input at a single node despite the docs and the vectorised
# body both expecting several.

pers_graph <- igraph::graph_from_data_frame(
  data.frame(from = c("a", "b", "c"), to = c("b", "c", "d")),
  directed = TRUE
)

pers_vec <- generate_personalisation_vec(
  pers_graph,
  node_weights = c(a = 3, c = 1)
)

expect_equal(
  current = sum(pers_vec),
  target = 1,
  info = "generate_personalisation_vec normalises to sum one"
)

expect_equal(
  current = unname(pers_vec[["a"]] / pers_vec[["c"]]),
  target = 3,
  info = "generate_personalisation_vec preserves the weight ratio"
)

# find_rbh_communities ---------------------------------------------------------

# Read S7::prop(object, "RbhGraph"), the class name, where the property is
# "rbh_graph". Errored unconditionally, and took plot_resolution_res() with it.

rbh_obj <- local({
  modules <- data.table::data.table(
    origin = rep(c("set_a", "set_b"), each = 20),
    module = rep(c("m1", "m2", "m3", "m4"), each = 10),
    gene = unlist(replicate(4, sample(letters, 10), simplify = FALSE))
  )
  obj <- RbhGraph(
    modules,
    rbh_type = "set",
    dataset_col = "origin",
    module_col = "module",
    value_col = "gene"
  )
  generate_rbh_graph(obj, minimum_similarity = 0)
})

expect_silent(
  current = find_rbh_communities(rbh_obj, parallel = FALSE, .verbose = FALSE),
  info = "find_rbh_communities reads the rbh_graph property"
)

# deprecated snf ---------------------------------------------------------------

# The wrapper constructed Snf(), which does not exist. The class is
# SimilarityNetworkFusion.

snf_res <- local({
  snf_mat <- matrix(rnorm(200), nrow = 20)
  rownames(snf_mat) <- sprintf("sample_%02d", 1:20)
  colnames(snf_mat) <- sprintf("feature_%02d", 1:10)
  suppressWarnings(
    snf(
      data = snf_mat,
      data_name = "continuous",
      snf_params = params_snf(k = 5L)
    )
  )
})

expect_true(
  current = S7::S7_inherits(snf_res, SimilarityNetworkFusion),
  info = "the deprecated snf() constructs a SimilarityNetworkFusion"
)

# bulk PCA gene clamping -------------------------------------------------------

# hvg_data[1:no_hvg_genes, ] produced NAs whenever the object held fewer genes
# than the 2500 default, and the next line died with subscript out of bounds.
# Every synthetic generator in the package produces fewer than 2500.

bulk_dge_obj <- local({
  syn <- synthetic_bulk_cor_matrix()
  meta <- data.table::data.table(
    sample_id = colnames(syn$counts),
    case_control = rep(c("case", "control"), each = 50),
    batch = rep(c("b1", "b2"), length.out = ncol(syn$counts))
  )
  obj <- BulkDge(raw_counts = syn$counts, meta_data = meta)
  obj <- qc_bulk_dge(obj, group_col = "case_control", .verbose = FALSE)
  normalise_bulk_dge(obj, group_col = "case_control", .verbose = FALSE)
})

expect_true(
  current = nrow(bulk_dge_obj@outputs$normalised_counts) < 2500L,
  info = "the fixture has fewer genes than the no_hvg_genes default"
)

expect_silent(
  current = calculate_pca_bulk_dge(bulk_dge_obj),
  info = "calculate_pca_bulk_dge clamps no_hvg_genes to the gene count"
)

expect_silent(
  current = suppressWarnings(batch_correction_bulk_dge(
    bulk_dge_obj,
    contrast_column = "case_control",
    batch_col = "batch"
  )),
  info = "batch_correction_bulk_dge passes no_hvg_genes into its PCA fallback"
)

# synthetic data generator -----------------------------------------------------

## visible return --------------------------------------------------------------

# The body ended on `res <- list(...)` with no trailing `res`, so the result
# came back invisibly, unlike both sibling generators.

expect_true(
  current = withVisible(
    generate_single_cell_test_data(
      syn_data_params = params_sc_synthetic_data(
        n_cells = 50L,
        n_genes = 40L
      )
    )
  )$visible,
  info = "generate_single_cell_test_data returns visibly"
)

## sample arguments are coupled ------------------------------------------------

# n_samples without sample_bias silently produced no sample_id column and only
# surfaced much later as a DuckDB binder error.

expect_error(
  current = params_sc_synthetic_data(n_samples = 6L),
  info = "n_samples without sample_bias is rejected up front"
)

expect_error(
  current = params_sc_synthetic_data(sample_bias = "even"),
  info = "sample_bias without n_samples is rejected up front"
)

expect_true(
  current = "sample_id" %in%
    names(
      generate_single_cell_test_data(
        syn_data_params = params_sc_synthetic_data(
          n_cells = 60L,
          n_genes = 40L,
          n_samples = 3L,
          sample_bias = "even"
        )
      )$obs
    ),
  info = "n_samples plus sample_bias does produce a sample_id column"
)

# ADT isotypes -----------------------------------------------------------------

# detect_adt_isotypes() matched case-insensitively, remove_adt_isotypes() did
# not, so the detector found features the remover then kept.

isotype_feats <- c("CD3", "IgG2a_Isotype", "CD19", "igg1_isotype")

expect_equal(
  current = sort(setdiff(isotype_feats, remove_adt_isotypes(isotype_feats))),
  target = sort(detect_adt_isotypes(isotype_feats)),
  info = "detect_adt_isotypes and remove_adt_isotypes agree on case"
)

# Gene Ontology ----------------------------------------------------------------

# get_go_levels() called igraph::dfs(father = ), deprecated in igraph 2.2.0, so
# every get_go_data_human() call told the user to report an issue.

go_warnings <- NULL

withCallingHandlers(
  invisible(get_go_data_human()),
  warning = function(w) {
    go_warnings <<- c(go_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

expect_false(
  current = any(grepl("father", go_warnings)),
  info = "get_go_data_human does not trip the igraph dfs deprecation"
)

# Cell Ranger round trip -------------------------------------------------------

# write_cellranger_output() wrote mat.mtx while get_cell_ranger_params()
# asserted matrix.mtx, so the package could not read back what it had written.

cellranger_dir <- sc_test_dir("cellranger_roundtrip")

local({
  data <- generate_single_cell_test_data(
    syn_data_params = params_sc_synthetic_data(n_cells = 60L, n_genes = 40L)
  )
  write_cellranger_output(
    f_path = cellranger_dir,
    counts = data$counts,
    obs = data$obs,
    var = data$var,
    rows = "genes",
    format_type = "tsv",
    .verbose = FALSE
  )
})

expect_true(
  current = file.exists(file.path(cellranger_dir, "matrix.mtx")),
  info = "write_cellranger_output writes matrix.mtx, as Cell Ranger does"
)

expect_silent(
  current = get_cell_ranger_params(cellranger_dir, has_hdr = TRUE),
  info = "get_cell_ranger_params reads back what the writer produced"
)

expect_equal(
  current = basename(
    get_cell_ranger_params(cellranger_dir, has_hdr = TRUE)$path_mtx
  ),
  target = "matrix.mtx",
  info = "the round trip resolves the matrix file"
)

## v2 naming still resolves ----------------------------------------------------

# Real Cell Ranger v2 directories call the feature file genes.tsv.

cellranger_v2 <- sc_test_dir("cellranger_v2")

invisible(file.create(
  file.path(cellranger_v2, c("barcodes.tsv", "genes.tsv", "matrix.mtx"))
))

expect_equal(
  current = basename(get_cell_ranger_params(cellranger_v2)$path_var),
  target = "genes.tsv",
  info = "get_cell_ranger_params still resolves the v2 genes.tsv naming"
)

expect_error(
  current = get_cell_ranger_params(sc_test_dir("cellranger_empty")),
  info = "an empty directory is rejected with a message, not a stale assertion"
)

sc_test_cleanup(cellranger_dir, cellranger_v2)

# synthetic_signal_matrix argument types ---------------------------------------

# no_grps and per_group asserted "R1" (double only) while their siblings
# asserted "I1", and all five are documented as Integer. So the documented
# 3L was rejected.

expect_silent(
  current = synthetic_signal_matrix(no_grps = 2L, per_group = 15L),
  info = "synthetic_signal_matrix accepts integers, as documented"
)

expect_silent(
  current = synthetic_signal_matrix(no_grps = 2, per_group = 15),
  info = "synthetic_signal_matrix still accepts doubles"
)

# duplicated stats helpers -----------------------------------------------------

# ot_harmonic_score, robust_scale and calculate_effect_size were each defined
# twice, and the functions_stats.R copy carried an invalid sprintf("%b").
# Alphabetical sourcing meant stats_helpers.R shadowed it.

expect_equal(
  current = length(robust_scale(c(1, 2, 3, 4, 100))),
  target = 5L,
  info = "robust_scale survives the de-duplication"
)

expect_true(
  current = is.numeric(ot_harmonic_score(c(0.9, 0.5, 0.1))),
  info = "ot_harmonic_score survives the de-duplication"
)

expect_silent(
  current = calculate_effect_size(
    mat_a = matrix(rnorm(100), nrow = 10),
    mat_b = matrix(rnorm(100), nrow = 10),
    small_sample_correction = TRUE
  ),
  info = "calculate_effect_size does not hit the invalid sprintf format"
)
