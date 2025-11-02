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
    annot_data = synth_annot
  ),
  info = "CisTarget - run_cistarget rejects unnamed gene set list"
)

expect_error(
  run_cistarget(
    gs_list = synth_gs_list,
    rankings = bad_matrix,
    annot_data = synth_annot
  ),
  info = "CisTarget - run_cistarget rejects bad matrix"
)

bad_annot <- synth_annot[, .(motif, TF)]

expect_error(
  run_cistarget(
    gs_list = synth_gs_list,
    rankings = synth_rankings,
    annot_data = bad_annot
  ),
  info = "CisTarget - run_cistarget rejects annotation without required columns"
)

expect_warning(
  run_cistarget(
    gs_list = list(no_overlap = c("fake_gene_999")),
    rankings = synth_rankings,
    annot_data = synth_annot
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
    )
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
