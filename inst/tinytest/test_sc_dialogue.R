# dialogue ---------------------------------------------------------------------

source("helper_sc.R", local = TRUE)

set.seed(123L)

test_temp_dir <- sc_test_dir("dialogue")

## testing parameters ----------------------------------------------------------

# a small permutation null: the point is the pipeline, not the tail of the
# p-value. 20 puts the floor at 0.05.
n_permutations <- 20L
k <- 2L

test_pmd_params <- params_dialogue_pmd(
  k = k,
  n_permutations = n_permutations,
  abn_c = 10L,
  n_genes = 12L,
  seed = 7L
)

## synthetic test data ---------------------------------------------------------

fixture <- dialogue_test_fixture()

n_cell_types <- length(fixture$features)
cell_types <- names(fixture$features)

sc_object <- dialogue_test_object(test_temp_dir, fixture)
sc_obs <- get_sc_obs(sc_object, filtered = TRUE)

## planted programme -----------------------------------------------------------

dialogue_res <- dialogue_sc(
  object = sc_object,
  cell_type_col = "cell_grp",
  sample_col = "sample_id",
  features = fixture$features,
  quality_col = "cell_quality",
  gene_ids = fixture$var$gene_id,
  pmd_params = test_pmd_params,
  .verbose = FALSE
)

expect_inherits(
  dialogue_res,
  "DialogueResult",
  info = "dialogue_sc() returns a DialogueResult"
)

# ground truth decides which programme is the planted one; everything below is
# then a statement about what DIALOGUE said of it
agreement_per_programme <- purrr::map(
  seq_len(k),
  \(p) {
    dialogue_latent_agreement(
      dialogue_res$cca_scores,
      sc_obs,
      fixture$latent,
      p
    )
  }
)
worst_agreement <- purrr::map_dbl(agreement_per_programme, min)
planted_programme <- which.max(worst_agreement)

expect_true(
  max(worst_agreement) > 0.7,
  info = paste(
    "some programme's sample-averaged canonical scores track the planted",
    "latent in every cell type"
  )
)

expect_true(
  all(
    dialogue_latent_agreement(
      dialogue_res$scores,
      sc_obs,
      fixture$latent,
      planted_programme
    ) >
      0.7
  ),
  info = "the final scores, not just the canonical ones, track the latent"
)

planted_rows <- dialogue_res$programmes[
  dialogue_res$programmes$programme == planted_programme,
]

expect_true(
  all(planted_rows$emp_p < 0.1),
  info = "every cell type pair of the planted programme is called significant"
)

expect_equal(
  length(dialogue_res$mcp_cell_types[[planted_programme]]),
  n_cell_types,
  info = "the planted programme spans every cell type"
)

# looser than the Rust fixture's 0.5, and deliberately so: there the planted
# layer goes in untouched, here the refit works off a handful of signature genes
# on log-normalised counts. Under the null this would sit near 0.05, so 0.3 is
# still a real statement.
expect_true(
  all(abs(dialogue_res$refit_fidelity[, planted_programme]) > 0.3),
  info = "the stage three refit stays anchored to the canonical score"
)

## signatures ------------------------------------------------------------------

permissive <- dialogue_res$signatures[
  dialogue_res$signatures$list == "permissive",
]

# every cell type must produce a signature and it must be dominated by that
# cell type's planted genes. Chance is n_planted / n_genes, about 2%. Asserting
# per cell type rather than "somewhere" is what makes this notice two of three
# signatures turning to noise.
for (ct in cell_types) {
  sig <- permissive[
    permissive$programme == planted_programme & permissive$cell_type == ct,
  ]

  expect_true(
    nrow(sig) > 0L,
    info = sprintf("%s produced a signature for the planted programme", ct)
  )
  expect_true(
    mean(sig$gene_id %in% fixture$planted[[ct]]) > 0.8,
    info = sprintf("%s signature is dominated by its planted genes", ct)
  )
}

expect_true(
  nrow(dialogue_res$verdicts) > 0L,
  info = "stage two produced gene associations"
)

expect_true(
  any(dialogue_res$verdicts$coefficient > 0),
  info = "the staged refit kept at least one gene"
)

## shapes and names ------------------------------------------------------------

expect_equal(
  nrow(dialogue_res$programmes),
  k * ncol(utils::combn(n_cell_types, 2L)),
  info = "one row per programme and cell type pair"
)

expect_equal(
  dim(dialogue_res$refit_fidelity),
  c(n_cell_types, k),
  info = "refit fidelity is cell types x programmes"
)

expect_equal(
  names(dialogue_res$scores),
  cell_types,
  info = "the score list is named by cell type"
)

expect_equal(
  purrr::map_int(dialogue_res$scores, nrow),
  purrr::map_int(fixture$features, nrow),
  info = "one score row per cell of that cell type"
)

expect_true(
  all(purrr::map_lgl(dialogue_res$scores, \(m) ncol(m) == k)),
  info = "one score column per programme"
)

expect_true(
  all(
    purrr::imap_lgl(dialogue_res$scores, \(m, ct) {
      setequal(rownames(m), sc_obs$cell_id[sc_obs$cell_grp == ct])
    })
  ),
  info = "score row names are that cell type's cells"
)

expect_true(
  all(dialogue_res$verdicts$gene_id %in% fixture$var$gene_id),
  info = "verdict gene ids resolve against the var table"
)

expect_true(
  setequal(dialogue_res$shared_samples, unique(sc_obs$sample_id)),
  info = "every sample is shared, as the fixture was built"
)

expect_true(
  all(purrr::map_lgl(
    dialogue_res$ws,
    \(m) nrow(m) <= ncol(fixture$features[[1]])
  )),
  info = "the weights cover at most the feature columns handed in"
)

## reproducibility -------------------------------------------------------------

dialogue_rerun <- dialogue_sc(
  object = sc_object,
  cell_type_col = "cell_grp",
  sample_col = "sample_id",
  features = fixture$features,
  quality_col = "cell_quality",
  gene_ids = fixture$var$gene_id,
  pmd_params = test_pmd_params,
  .verbose = FALSE
)

expect_equal(
  dialogue_res$programmes,
  dialogue_rerun$programmes,
  info = "same seed, same empirical p-values"
)

expect_equal(
  dialogue_res$scores,
  dialogue_rerun$scores,
  info = "same seed, same scores"
)

## input validation ------------------------------------------------------------

expect_error(
  dialogue_sc(
    object = sc_object,
    cell_type_col = "cell_grp",
    sample_col = "sample_id",
    features = fixture$features[1],
    quality_col = "cell_quality",
    pmd_params = test_pmd_params,
    .verbose = FALSE
  ),
  pattern = "at least two cell types",
  info = "one cell type is rejected before Rust sees it"
)

expect_error(
  dialogue_sc(
    object = sc_object,
    cell_type_col = "cell_grp",
    sample_col = "sample_id",
    features = stats::setNames(fixture$features, c("a", "b", "c")),
    quality_col = "cell_quality",
    pmd_params = test_pmd_params,
    .verbose = FALSE
  ),
  pattern = "absent from",
  info = "feature names that are not cell types are rejected"
)

expect_error(
  dialogue_sc(
    object = sc_object,
    cell_type_col = "cell_grp",
    sample_col = "sample_id",
    features = purrr::map(fixture$features, \(m) m[, 1L, drop = FALSE]),
    quality_col = "cell_quality",
    pmd_params = test_pmd_params,
    .verbose = FALSE
  ),
  pattern = "at least 2",
  info = "a single feature column is rejected"
)

expect_error(
  dialogue_sc(
    object = sc_object,
    cell_type_col = "cell_grp",
    sample_col = "sample_id",
    features = purrr::map(fixture$features, \(m) unname(m)),
    quality_col = "cell_quality",
    pmd_params = test_pmd_params,
    .verbose = FALSE
  ),
  pattern = "no row names",
  info = "features without row names are rejected, they cannot be aligned"
)

expect_error(
  dialogue_sc(
    object = sc_object,
    cell_type_col = "cell_grp",
    sample_col = "sample_id",
    features = purrr::map(fixture$features, \(m) m[-1L, , drop = FALSE]),
    quality_col = "cell_quality",
    pmd_params = test_pmd_params,
    .verbose = FALSE
  ),
  pattern = "missing",
  info = "features that do not cover every cell are rejected"
)

# collapse the samples so only three survive in every cell type
few_samples_obs <- data.table::copy(fixture$obs)
few_samples_obs[, sample_id := sprintf("sample_%02d", seq_len(.N) %% 3L + 1L)]
few_samples_dir <- sc_test_dir("dialogue_few_samples")
few_samples_object <- dialogue_test_object(
  few_samples_dir,
  fixture,
  obs = few_samples_obs
)

expect_error(
  dialogue_sc(
    object = few_samples_object,
    cell_type_col = "cell_grp",
    sample_col = "sample_id",
    features = fixture$features,
    quality_col = "cell_quality",
    pmd_params = test_pmd_params,
    .verbose = FALSE
  ),
  pattern = "present in every cell type",
  info = "too few shared samples is caught with a message about the design"
)

## parameter wrappers ----------------------------------------------------------

expect_error(params_dialogue_pmd(k = 0L))
expect_error(params_dialogue_pmd(n_permutations = 1L))
expect_error(params_dialogue_pmd(cap = 0.5))
expect_error(params_dialogue_pmd(averaging = "geometric"))
expect_error(params_dialogue_refine(support_p = 0))
expect_error(params_dialogue_refine(min_support_fraction = 1.5))

expect_error(
  bixverse:::assertDialoguePmd(within(params_dialogue_pmd(), rm(cap))),
  info = "a parameter list missing a field is rejected"
)

expect_error(
  bixverse:::assertDialogueHlm(utils::modifyList(
    params_dialogue_hlm(),
    list(satterthwaite = 1.5)
  )),
  info = "a parameter list with the wrong type is rejected"
)

expect_error(
  bixverse:::assertDialogueRefine(utils::modifyList(
    params_dialogue_refine(),
    list(strict_p = 2)
  )),
  info = "a parameter list out of range is rejected"
)

expect_true(
  bixverse:::checkDialoguePmd(params_dialogue_pmd()),
  info = "the defaults pass their own check"
)

## subset ----------------------------------------------------------------------

subset_object <- SingleCellsSubset(
  sc_object,
  grouping_column = "sample_id",
  group = "sample_01"
)

expect_error(
  dialogue_sc(
    object = subset_object,
    cell_type_col = "cell_grp",
    sample_col = "sample_id",
    features = fixture$features,
    quality_col = "cell_quality",
    pmd_params = test_pmd_params,
    .verbose = FALSE
  ),
  pattern = "present in every cell type",
  info = paste(
    "the subset dispatches and reaches the same design check;",
    "one sample cannot support DIALOGUE"
  )
)

## meta cells ------------------------------------------------------------------

# meta cells have to be built *within* samples, otherwise the random intercept
# in stage two has no well-defined level. This is the workflow the docs point
# at, so the test builds it the same way.
mc_dir <- sc_test_dir("dialogue_mc")
mc_sc_object <- dialogue_test_object(mc_dir, fixture)

per_sample_mc <- purrr::map(sort(unique(fixture$obs$sample_id)), \(s) {
  sub <- SingleCellsSubset(
    mc_sc_object,
    grouping_column = "sample_id",
    group = s
  )
  sub <- find_hvg_sc(sub, hvg_no = 50L, .verbose = FALSE)
  sub <- calculate_pca_sc(sub, no_pcs = 10L, .verbose = FALSE)
  sub <- find_neighbours_sc(
    sub,
    neighbours_params = params_sc_neighbours(knn = list(k = 5L)),
    .verbose = FALSE
  )
  generate_bt_meta_cells_sc(
    sub,
    sc_meta_cell_params = params_sc_bt_metacells(target_no_metacells = 15L),
    .verbose = FALSE
  )
})
names(per_sample_mc) <- sort(unique(fixture$obs$sample_id))

mc_object <- merge_meta_cells(per_sample_mc, .verbose = FALSE)
mc_obs <- S7::prop(mc_object, "obs_table")
membership <- mc_obs$original_cell_idx

# dominant cell type per meta cell, as the other meta cell tests do. Meta cells
# built on a within-sample kNN are not pure, so this is an annotation, not a
# partition of the planted design.
mc_object[["cell_grp"]] <- purrr::map_chr(membership, \(i) {
  names(which.max(table(fixture$obs$cell_grp[i])))
})
mc_object[["sample_id"]] <- mc_obs$source_id
mc_object[["cell_quality"]] <- purrr::map_dbl(membership, \(i) {
  mean(fixture$obs$cell_quality[i])
})

mc_obs <- S7::prop(mc_object, "obs_table")
mc_cell_types <- sort(unique(mc_obs$cell_grp))

# meta cell features are the planted features averaged over the cells that went
# into each meta cell
feature_all <- do.call(rbind, fixture$features)[fixture$obs$cell_id, ]
mc_feature_all <- t(vapply(
  membership,
  \(i) colMeans(feature_all[i, , drop = FALSE]),
  numeric(ncol(feature_all))
))
rownames(mc_feature_all) <- mc_obs$meta_cell_id

mc_features <- purrr::map(mc_cell_types, \(ct) {
  mc_feature_all[mc_obs$meta_cell_id[mc_obs$cell_grp == ct], , drop = FALSE]
})
names(mc_features) <- mc_cell_types

mc_res <- dialogue_sc(
  object = mc_object,
  cell_type_col = "cell_grp",
  sample_col = "sample_id",
  features = mc_features,
  quality_col = "cell_quality",
  gene_ids = fixture$var$gene_id,
  pmd_params = params_dialogue_pmd(
    k = k,
    n_permutations = n_permutations,
    abn_c = 2L,
    n_genes = 12L,
    seed = 7L
  ),
  hlm_params = params_dialogue_hlm(min_cells_per_sample = 2L),
  .verbose = FALSE
)

expect_inherits(
  mc_res,
  "DialogueResult",
  info = "mc - dialogue_sc() runs off the in-memory meta cell counts"
)

expect_equal(
  mc_res$source_class,
  "MetaCells",
  info = "mc - the result records where it came from"
)

expect_equal(
  mc_res$cell_types,
  mc_cell_types,
  info = "mc - the cell types come back in the order they went in"
)

# not every sample carries all three cell types once the cells are aggregated
# into 15 meta cells, so the shared set is a subset rather than the whole thing
expect_true(
  all(mc_res$shared_samples %in% mc_obs$sample_id),
  info = "mc - the shared samples resolve back to sample labels"
)

expect_true(
  length(mc_res$shared_samples) >= 5L,
  info = "mc - enough samples survive the intersection to decompose on"
)

expect_equal(
  purrr::map_int(mc_res$scores, nrow),
  purrr::map_int(mc_features, nrow),
  info = "mc - one score row per meta cell of that cell type"
)

expect_true(
  all(purrr::imap_lgl(mc_res$scores, \(m, ct) {
    setequal(rownames(m), mc_obs$meta_cell_id[mc_obs$cell_grp == ct])
  })),
  info = "mc - score row names are meta cell identifiers"
)

expect_true(
  all(mc_res$verdicts$gene_id %in% fixture$var$gene_id),
  info = "mc - verdict gene ids resolve against the var table"
)

# the meta cells mix cell types, so this is a statement about the pipeline
# finding structure, not about recovering the planted genes
expect_true(
  min(mc_res$programmes$emp_p) < 0.1,
  info = "mc - at least one programme is called significant"
)

## clean up --------------------------------------------------------------------

sc_test_cleanup(test_temp_dir, few_samples_dir, mc_dir)
