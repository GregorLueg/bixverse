# cistarget tests --------------------------------------------------------------

## synthetic data --------------------------------------------------------------

# fmt: skip
synth_rankings <- matrix(
  c(
    1L, 5L, 4L, 2L,
    2L, 2L, 1L, 5L,
    4L, 3L, 2L, 1L,
    3L, 1L, 5L, 4L,
    5L, 4L, 3L, 3L
  ),
  nrow = 5,
  ncol = 4,
  byrow = TRUE
)
rownames(synth_rankings) <- sprintf("gene_%i", 1:5)
colnames(synth_rankings) <- sprintf("motif_%i", 1:4)

bad_matrix <- synth_rankings
rownames(bad_matrix) <- NULL

synth_annot <- data.table::data.table(
  motif = rep(sprintf("motif_%i", 1:4), each = 2),
  TF = c(
    "TF1",
    "TF2",
    "TF3",
    "TF4",
    "TF5",
    "TF6",
    "TF7",
    "TF8"
  ),
  direct_annotation = c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE),
  inferred_orthology = c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE),
  inferred_motif_sim = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
  annotationSource = factor(c(
    "directAnnotation",
    "inferredBy_Orthology",
    "directAnnotation",
    "inferredBy_Orthology",
    "inferredBy_MotifSimilarity",
    "inferredBy_MotifSimilarity",
    "directAnnotation",
    "inferredBy_Orthology"
  )),
  description = "test"
)
data.table::setkeyv(synth_annot, c("motif", "TF"))

synth_gs_list <- list(
  test_set_1 = c("gene_1", "gene_2", "gene_3"),
  test_set_2 = c("gene_3", "gene_4", "gene_5"),
  empty_set = c("gene_999") # no overlap
)

## tests -----------------------------------------------------------------------

### helper function ------------------------------------------------------------

mock_rs_result <- list(
  motif_idx = c(1L, 2L),
  nes = c(4.5, 3.2),
  auc = c(0.85, 0.72),
  rank_at_max = c(10L, 15L),
  n_enriched = c(3L, 2L),
  leading_edge = list(c(1L, 2L, 3L), c(3L, 4L))
)

processed <- bixverse:::process_cistarget_res(
  cs_ls = mock_rs_result,
  gs_name = "test_set",
  represented_motifs = colnames(synth_rankings),
  represented_genes = rownames(synth_rankings)
)

expect_true(
  data.table::is.data.table(processed),
  info = "CisTarget - process_cistarget_res returns data.table"
)

expect_equal(
  nrow(processed),
  2L,
  info = "CisTarget - process_cistarget_res correct number of rows"
)

expect_equal(
  processed$gs_name[1],
  "test_set",
  info = "CisTarget - process_cistarget_res preserves gene set name"
)

expect_true(
  all(c("motif", "nes", "auc", "leading_edge_genes") %in% names(processed)),
  info = "CisTarget - process_cistarget_res has expected columns"
)

mock_empty_result <- list(
  motif_idx = integer(0),
  nes = numeric(0),
  auc = numeric(0),
  rank_at_max = integer(0),
  n_enriched = integer(0),
  leading_edge = list()
)

processed_empty <- bixverse:::process_cistarget_res(
  cs_ls = mock_empty_result,
  gs_name = "empty_set",
  represented_motifs = colnames(synth_rankings),
  represented_genes = rownames(synth_rankings)
)

expect_null(
  processed_empty,
  info = "CisTarget - process_cistarget_res returns NULL for empty results"
)

### main function --------------------------------------------------------------

#### errors, warnings ----------------------------------------------------------

expect_error(
  run_cistarget(
    gs_list = list(c("gene_1")),
    rankings = synth_rankings,
    annot_data = synth_annot,
    .verbose = FALSE
  ),
  info = "CisTarget - run_cistarget rejects unnamed gene set list"
)

expect_error(
  run_cistarget(
    gs_list = synth_gs_list,
    rankings = bad_matrix,
    annot_data = synth_annot,
    .verbose = FALSE
  ),
  info = "CisTarget - run_cistarget rejects bad matrix"
)

bad_annot <- synth_annot[, .(motif, TF)]

expect_error(
  run_cistarget(
    gs_list = synth_gs_list,
    rankings = synth_rankings,
    annot_data = bad_annot,
    .verbose = FALSE
  ),
  info = "CisTarget - run_cistarget rejects annotation without required columns"
)

expect_warning(
  run_cistarget(
    gs_list = list(no_overlap = c("fake_gene_999")),
    rankings = synth_rankings,
    annot_data = synth_annot,
    .verbose = FALSE
  ),
  info = "CisTarget - run_cistarget warns about gene sets with zero overlap"
)

#### actual run ----------------------------------------------------------------

result <- suppressWarnings(
  run_cistarget(
    gs_list = synth_gs_list,
    rankings = synth_rankings,
    annot_data = synth_annot,
    cis_target_params = params_cistarget(
      auc_threshold = 1,
      nes_threshold = 0.2
    ),
    .verbose = FALSE
  )
)

expect_true(
  data.table::is.data.table(result),
  info = "CisTarget - run_cistarget returns data.table"
)

expected_cols <- c(
  "gs_name",
  "motif",
  "nes",
  "auc",
  "TF_highConf",
  "TF_lowConf"
)

expect_true(
  all(expected_cols %in% names(result)),
  info = "CisTarget - run_cistarget has expected output columns"
)

expect_true(
  current = !"empty_set" %in% result$gs_name,
  info = "CisTarget - empty set not in results"
)

## cistarget parameters --------------------------------------------------------

cistarget_defaults <- params_cistarget()

expect_equal(
  current = cistarget_defaults$max_rank,
  target = 5000L,
  info = "CisTarget - max_rank defaults to the RcisTarget value"
)

expect_equal(
  current = cistarget_defaults$n_mean,
  target = 100L,
  info = "CisTarget - n_mean defaults to the RcisTarget value"
)

expect_true(
  current = isTRUE(bixverse:::checkCistargetParams(cistarget_defaults)),
  info = "CisTarget - the default params pass their own validator"
)

expect_false(
  current = isTRUE(bixverse:::checkCistargetParams(
    cistarget_defaults[c("auc_threshold", "nes_threshold")]
  )),
  info = "CisTarget - the validator catches a truncated params list"
)

# max_rank is clamped to the ranking depth, so a huge value must not error
deep_max_rank <- suppressWarnings(
  run_cistarget(
    gs_list = synth_gs_list[1:2],
    rankings = synth_rankings,
    annot_data = synth_annot,
    cis_target_params = params_cistarget(
      auc_threshold = 1,
      nes_threshold = 0.2,
      max_rank = 1e6L
    ),
    .verbose = FALSE
  )
)

expect_equal(
  current = deep_max_rank,
  target = result[gs_name != "empty_set"],
  info = "CisTarget - max_rank above the ranking depth is clamped"
)

## regulon building ------------------------------------------------------------

synth_imp <- matrix(
  seq_len(6) / 21,
  nrow = 3,
  ncol = 2,
  dimnames = list(sprintf("gene_%i", 1:3), c("TF1", "TF2"))
)

synth_grn <- bixverse:::new_scenic_grn(
  importance_matrix = synth_imp,
  gene_ids = rownames(synth_imp),
  tf_ids = colnames(synth_imp),
  params = list()
)

synth_grn$tf_to_gene_results <- data.table::data.table(
  tf = c(rep("TF1", 4L), rep("TF2", 2L)),
  gene = c(sprintf("gene_%i", 1:4), sprintf("gene_%i", 5:6)),
  importance = seq(0.6, 0.1, length.out = 6L),
  pairwise_cor = c(0.5, 0.4, 0.3, -0.6, 0.2, 0.1),
  cor_sign = c(1L, 1L, 1L, -1L, 1L, 1L),
  in_leading_edge = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
)

# TF1 keeps three activating leading edge genes plus itself, TF2 is left with
# gene_5 plus itself
regulons <- build_regulons(
  synth_grn,
  min_genes = 3L,
  .verbose = FALSE
)

expect_equal(
  current = names(regulons),
  target = "TF1",
  info = "build_regulons - TF2 drops below min_genes after the filters"
)

expect_equal(
  current = length(regulons$TF1),
  target = 4L,
  info = "build_regulons - three leading edge targets plus the TF"
)

expect_true(
  current = "TF1" %in% regulons$TF1,
  info = "build_regulons - the TF is added to its own regulon"
)

expect_false(
  current = "gene_4" %in% regulons$TF1,
  info = "build_regulons - repressing links are dropped by default"
)

regulons_both <- build_regulons(
  synth_grn,
  min_genes = 1L,
  mode = "both",
  .verbose = FALSE
)

expect_true(
  current = all(c("TF1_pos", "TF1_neg") %in% names(regulons_both)),
  info = "build_regulons - mode 'both' splits the sign into the name"
)

expect_false(
  current = "gene_6" %in% unlist(regulons_both),
  info = "build_regulons - genes outside the leading edge are dropped"
)

regulons_no_le <- build_regulons(
  synth_grn,
  use_leading_edge = FALSE,
  add_tf = FALSE,
  min_genes = 1L,
  .verbose = FALSE
)

expect_true(
  current = "gene_6" %in% regulons_no_le$TF2,
  info = "build_regulons - use_leading_edge = FALSE keeps everything"
)

expect_false(
  current = "TF2" %in% regulons_no_le$TF2,
  info = "build_regulons - add_tf = FALSE leaves the TF out"
)

empty_grn <- bixverse:::new_scenic_grn(
  importance_matrix = synth_imp,
  gene_ids = rownames(synth_imp),
  tf_ids = colnames(synth_imp),
  params = list()
)

expect_warning(
  current = build_regulons(empty_grn, .verbose = FALSE),
  info = "build_regulons - warns when there are no TF to gene pairs"
)

## binarisation ----------------------------------------------------------------

set.seed(42L)
synth_auc <- cbind(
  bimodal = c(rnorm(100L, 0.1, 0.01), rnorm(100L, 0.6, 0.01)),
  flat = rnorm(200L, 0.3, 0.05)
)
rownames(synth_auc) <- sprintf("cell_%03d", 1:200)

binarised <- binarise_regulon_activity(synth_auc, .verbose = FALSE)

expect_equal(
  current = dim(binarised$binary),
  target = dim(synth_auc),
  info = "binarise - the binary matrix keeps the input dimensions"
)

expect_equal(
  current = dimnames(binarised$binary),
  target = dimnames(synth_auc),
  info = "binarise - dimnames are preserved"
)

expect_true(
  current = binarised$thresholds[regulon == "bimodal", bimodal],
  info = "binarise - the two-mode regulon is called bimodal"
)

expect_equal(
  current = binarised$thresholds[regulon == "bimodal", n_cells_on],
  target = 100L,
  info = "binarise - the two-mode regulon switches on the upper mode"
)

expect_false(
  current = binarised$thresholds[regulon == "flat", bimodal],
  info = "binarise - the single-mode regulon falls back to mean + 2 sd"
)

expect_equal(
  current = binarised$thresholds[regulon == "flat", threshold],
  target = mean(synth_auc[, "flat"]) + 2 * sd(synth_auc[, "flat"]),
  tolerance = 1e-8,
  info = "binarise - the unimodal threshold is mean + 2 sd"
)

expect_error(
  current = binarise_regulon_activity(
    synth_auc,
    binarise_params = list(bw_adjust = 1, n_grid = 2)
  ),
  info = "binarise - the params validator rejects a malformed list"
)
