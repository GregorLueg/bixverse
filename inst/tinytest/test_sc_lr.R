# nichenet end-to-end ----------------------------------------------------------

set.seed(123L)

## tests -----------------------------------------------------------------------

### helpers --------------------------------------------------------------------

#### rank_scale ----------------------------------------------------------------

expect_equal(
  current = bixverse:::.rank_scale(c(1, 2, 3, 4)),
  target = c(0.25, 0.5, 0.75, 1.0),
  info = "rank_scale - monotone vector"
)

expect_equal(
  current = bixverse:::.rank_scale(c(1, 1, 2, 2)),
  target = c(1.5 / 3.5, 1.5 / 3.5, 1.0, 1.0),
  info = "rank_scale - average ties"
)

expect_true(
  current = max(bixverse:::.rank_scale(rnorm(50))) == 1,
  info = "rank_scale - max is one"
)

expect_true(
  current = min(bixverse:::.rank_scale(rnorm(50))) > 0,
  info = "rank_scale - min is positive"
)

#### scaling_zscore ------------------------------------------------------------

expect_equivalent(
  current = mean(bixverse:::.scaling_zscore(rnorm(100))),
  target = 0,
  tolerance = 1e-10,
  info = "scaling_zscore - mean is zero"
)

expect_equivalent(
  current = sd(bixverse:::.scaling_zscore(rnorm(100))),
  target = 1,
  tolerance = 1e-10,
  info = "scaling_zscore - sd is one"
)

expect_equal(
  current = bixverse:::.scaling_zscore(c(1, 1, 1, 1)),
  target = c(0, 0, 0, 0),
  info = "scaling_zscore - zero variance returns zeros"
)

#### scale_quantile_adapted ----------------------------------------------------

x <- 1:100
y <- bixverse:::.scale_quantile_adapted(x, outlier_cutoff = 0.05)

expect_true(
  current = all(y >= 0 & y <= 1),
  info = "scale_quantile_adapted - within [0, 1]"
)

expect_true(
  current = sum(y == 0) > 1 && sum(y == 1) > 1,
  info = "scale_quantile_adapted - clipping flattens the tails"
)

expect_equal(
  current = bixverse:::.scale_quantile_adapted(rep(1, 10)),
  target = rep(0.5, 10),
  info = "scale_quantile_adapted - degenerate input returns 0.5"
)

#### assertions ----------------------------------------------------------------

good_de <- data.table::data.table(
  cluster_id = c("S1", "S1", "S2"),
  gene = c("L1", "L2", "L1"),
  lfc = c(1.0, -0.5, 0.2),
  pval = c(0.01, 0.5, 0.8)
)

expect_true(
  current = data.table::is.data.table(bixverse:::assertDegTable(good_de, "x")),
  info = "assert_de_table - returns invisible data.table on success"
)

bad_de_missing_col <- data.table::copy(good_de)[, lfc := NULL]
expect_error(
  current = bixverse:::assertDegTable(bad_de_missing_col, "x"),
  info = "assert_de_table - missing column"
)

bad_de_bad_pval <- data.table::copy(good_de)
bad_de_bad_pval[, pval := c(0.01, 2, 0.8)]
expect_error(
  current = bixverse:::assertDegTable(bad_de_bad_pval, "x"),
  info = "assert_de_table - pval outside [0, 1]"
)

#### resolve_weights -----------------------------------------------------------

expect_equal(
  current = unname(bixverse:::resolve_weights(NULL, "case_control", TRUE)),
  target = rep(1, 7),
  info = "resolve_weights - case_control defaults"
)

expect_equal(
  current = unname(bixverse:::resolve_weights(NULL, "one_condition", FALSE)),
  target = c(1, 1, 1, 1, 1, 0, 0),
  info = "resolve_weights - one_condition defaults"
)

expect_error(
  current = bixverse:::resolve_weights(NULL, "case_control", FALSE),
  info = "resolve_weights - missing condition_de but non-zero condition weights"
)

### ligand-target influence ----------------------------------------------------

#### synthetic network ---------------------------------------------------------

# Two disjoint signalling components:
#
#   L1 -> SIG1 -> TF1 -> {G1, G2, G3, G4}     (component A)
#   L2 -> SIG2 -> TF2 -> {G5, G6, G7, G8}     (component B)
#   L3 -> SIG3                                 (dead end, no TF)
#   L4 -> SIG4                                 (dead end, no TF)
#
# Expected signal:
# - L1's regulatory potential is positive on G1..G4 and zero on G5..G8.
# - L2's regulatory potential is positive on G5..G8 and zero on G1..G4.
# - L3 and L4 have zero regulatory potential on all target genes.

ppi_network <- data.table::data.table(
  from = c("L1", "SIG1", "L2", "SIG2", "L3", "L4"),
  to = c("SIG1", "TF1", "SIG2", "TF2", "SIG3", "SIG4"),
  weight = 1.0
)

grn_network <- data.table::data.table(
  from = c(rep("TF1", 4), rep("TF2", 4)),
  to = c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8"),
  weight = 1.0
)

ligand_seeds <- list(L1 = "L1", L2 = "L2", L3 = "L3", L4 = "L4")

#### construction --------------------------------------------------------------

# `ltf_cutoff = 0` skips the per-seed quantile threshold so the full PPR mass
# reaches the matmul. With only 18 nodes in the universe, the default cutoff
# of 0.99 would prune everything except the seed itself.

inf_obj <- generate_ligand_target_influence(
  ligand_seeds = ligand_seeds,
  ppi_network = ppi_network,
  grn_network = grn_network,
  params = params_ligand_target(ltf_cutoff = 0)
)

expect_true(
  current = checkmate::testClass(inf_obj, "LigandTargetInfluence"),
  info = "ltf - returns a LigandTargetInfluence object"
)

inf_mat <- get_influence(inf_obj)

# Universe: 10 PPI nodes (L1..L4, SIG1..SIG4, TF1, TF2) + 8 GRN-only genes
expect_equal(
  current = dim(inf_mat),
  target = c(4L, 18L),
  info = "ltf - correct dimensions (4 ligands x 18 nodes)"
)

expect_true(
  current = all(c("L1", "L2", "L3", "L4") %in% rownames(inf_mat)),
  info = "ltf - row names match ligand seeds"
)

expect_true(
  current = all(
    c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8") %in% colnames(inf_mat)
  ),
  info = "ltf - column names include target genes"
)

#### signal recovery -----------------------------------------------------------

tf1_targets <- c("G1", "G2", "G3", "G4")
tf2_targets <- c("G5", "G6", "G7", "G8")
zero_tol <- 1e-9

expect_true(
  current = all(inf_mat["L1", tf1_targets] > 0),
  info = "ltf - L1 has positive regulatory potential on TF1 targets"
)

expect_true(
  current = all(inf_mat["L1", tf2_targets] < zero_tol),
  info = "ltf - L1 has zero regulatory potential on TF2 targets (disjoint)"
)

expect_true(
  current = all(inf_mat["L2", tf2_targets] > 0),
  info = "ltf - L2 has positive regulatory potential on TF2 targets"
)

expect_true(
  current = all(inf_mat["L2", tf1_targets] < zero_tol),
  info = "ltf - L2 has zero regulatory potential on TF1 targets (disjoint)"
)

expect_true(
  current = all(inf_mat["L3", c(tf1_targets, tf2_targets)] < zero_tol),
  info = "ltf - L3 has zero regulatory potential (no TF in its component)"
)

expect_true(
  current = all(inf_mat["L4", c(tf1_targets, tf2_targets)] < zero_tol),
  info = "ltf - L4 has zero regulatory potential (no TF in its component)"
)

# The signal should be symmetric: the magnitude reaching each target from its
# matched ligand should be the same regardless of which component (graphs are
# structurally identical).
expect_equivalent(
  current = inf_mat["L1", "G1"],
  target = inf_mat["L2", "G5"],
  tolerance = 1e-9,
  info = "ltf - symmetric components produce equal magnitudes"
)

#### within-row uniformity -----------------------------------------------------

# All TF1 targets receive the same contribution from L1 because GRN weights
# are equal and they share the same upstream TF.
expect_equivalent(
  current = sd(inf_mat["L1", tf1_targets]),
  target = 0,
  tolerance = 1e-9,
  info = "ltf - equal-weight GRN edges produce uniform target values"
)

### ligand activity scoring ----------------------------------------------------

#### construction --------------------------------------------------------------

activity <- ligand_activity_scores(
  ligand_influence = inf_obj,
  gene_sets = list(set_A = tf1_targets, set_B = tf2_targets)
)

expect_true(
  current = checkmate::testDataTable(activity, nrows = 8L),
  info = "activity - one row per (gene_set, ligand)"
)

expected_cols <- c(
  "gene_set",
  "ligand",
  "auroc",
  "aupr",
  "aupr_corrected",
  "pearson",
  "spearman"
)
expect_true(
  current = all(expected_cols %in% names(activity)),
  info = "activity - expected columns"
)

#### top hits per gene set -----------------------------------------------------

top_A <- activity[gene_set == "set_A"][order(-aupr_corrected)][1]
top_B <- activity[gene_set == "set_B"][order(-aupr_corrected)][1]

expect_true(
  current = top_A$ligand == "L1",
  info = "activity - L1 ranked highest for set_A (TF1 targets)"
)

expect_true(
  current = top_B$ligand == "L2",
  info = "activity - L2 ranked highest for set_B (TF2 targets)"
)

#### perfect separation --------------------------------------------------------

# Positive predictions live exclusively on the matched targets, so positives
# strictly outrank negatives and AUROC/AUPR are exactly 1.
expect_equivalent(
  current = activity[gene_set == "set_A" & ligand == "L1", auroc],
  target = 1.0,
  tolerance = 1e-9,
  info = "activity - L1 achieves perfect AUROC on set_A"
)

expect_equivalent(
  current = activity[gene_set == "set_A" & ligand == "L1", aupr],
  target = 1.0,
  tolerance = 1e-9,
  info = "activity - L1 achieves perfect AUPR on set_A"
)

expect_equivalent(
  current = activity[gene_set == "set_B" & ligand == "L2", auroc],
  target = 1.0,
  tolerance = 1e-9,
  info = "activity - L2 achieves perfect AUROC on set_B"
)

# Pearson/Spearman: prediction vector is proportional to response indicator
# (constant on positives, zero on negatives) -> correlation is exactly 1.
expect_equivalent(
  current = activity[gene_set == "set_A" & ligand == "L1", pearson],
  target = 1.0,
  tolerance = 1e-9,
  info = "activity - L1 achieves perfect Pearson on set_A"
)

expect_equivalent(
  current = activity[gene_set == "set_A" & ligand == "L1", spearman],
  target = 1.0,
  tolerance = 1e-9,
  info = "activity - L1 achieves perfect Spearman on set_A"
)

#### aupr_corrected against random baseline ------------------------------------

# n_pos = 4, n = 18 -> random AUPR = 4 / 18
expect_equivalent(
  current = activity[gene_set == "set_A" & ligand == "L1", aupr_corrected],
  target = 1.0 - 4 / 18,
  tolerance = 1e-9,
  info = "activity - aupr_corrected on set_A = 1 - n_pos / n"
)

#### mismatched ligand is below random -----------------------------------------

expect_true(
  current = activity[gene_set == "set_A" & ligand == "L2", auroc] < 0.5,
  info = "activity - L2 scores below random AUROC on set_A"
)

expect_true(
  current = activity[gene_set == "set_A" & ligand == "L2", aupr_corrected] < 0,
  info = "activity - L2 has negative aupr_corrected on set_A"
)

### prioritise_interactions: synthetic data ------------------------------------

set.seed(123L)

#### setup ---------------------------------------------------------------------

senders <- c("S1", "S2", "S3")
receivers <- "R1"

lr_network <- data.table::data.table(
  ligand = c("L1", "L2", "L3", "L4"),
  receptor = c("REC1", "REC2", "REC1", "REC3")
)

# DE: implant signal in (S1, L1) and matching receptor
celltype_de <- data.table::CJ(
  cluster_id = c(senders, receivers),
  gene = c("L1", "L2", "L3", "L4", "REC1", "REC2", "REC3")
)
celltype_de[, lfc := rnorm(.N, 0, 0.3)]
celltype_de[, pval := runif(.N, 0.1, 1)]
celltype_de[cluster_id == "S1" & gene == "L1", `:=`(lfc = 3.0, pval = 1e-10)]
celltype_de[cluster_id == "R1" & gene == "REC1", `:=`(lfc = 2.5, pval = 1e-9)]

# Expression info: L1 enriched in S1, REC1 enriched in R1
expression_info <- data.table::CJ(
  cluster_id = c(senders, receivers),
  gene = c("L1", "L2", "L3", "L4", "REC1", "REC2", "REC3")
)
expression_info[, avg_expr := runif(.N, 0, 0.1)]
expression_info[cluster_id == "S1" & gene == "L1", avg_expr := 5.0]
expression_info[cluster_id == "R1" & gene == "REC1", avg_expr := 5.0]

# Activity: L1 dominant
ligand_activities <- data.table::data.table(
  ligand = c("L1", "L2", "L3", "L4"),
  aupr_corrected = c(0.6, 0.05, 0.02, 0.01)
)

#### basic run -----------------------------------------------------------------

result <- prioritise_interactions(
  celltype_de = celltype_de,
  expression_info = expression_info,
  ligand_activities = ligand_activities,
  lr_network = lr_network,
  senders_oi = senders,
  receivers_oi = receivers,
  scenario = "one_condition"
)

expect_true(
  current = checkmate::testDataTable(result, min.rows = 1L),
  info = "prioritise - returns a non-empty data.table"
)

expected_cols <- c(
  "sender",
  "receiver",
  "ligand",
  "receptor",
  "prioritisation_score",
  "prioritisation_rank"
)
expect_true(
  current = all(expected_cols %in% names(result)),
  info = "prioritise - expected output columns"
)

expect_true(
  current = !is.unsorted(-result$prioritisation_score),
  info = "prioritise - sorted by descending score"
)

expect_equal(
  current = result$prioritisation_rank,
  target = data.table::frank(-result$prioritisation_score, ties.method = "min"),
  info = "prioritise - rank matches score order"
)

#### signal recovery -----------------------------------------------------------

top_hit <- result[1, .(sender, ligand, receptor)]

expect_true(
  current = top_hit$sender == "S1" &&
    top_hit$ligand == "L1" &&
    top_hit$receptor == "REC1",
  info = "prioritise - top hit recovers implanted signal"
)

#### weight zero skips component -----------------------------------------------

weights_no_activity <- c(
  de_ligand = 1,
  de_receptor = 1,
  activity_scaled = 0,
  exprs_ligand = 1,
  exprs_receptor = 1,
  ligand_condition_specificity = 0,
  receptor_condition_specificity = 0
)

result_no_act <- prioritise_interactions(
  celltype_de = celltype_de,
  expression_info = expression_info,
  ligand_activities = ligand_activities,
  lr_network = lr_network,
  senders_oi = senders,
  receivers_oi = receivers,
  weights = weights_no_activity,
  scenario = "one_condition"
)

expect_false(
  current = "scaled_activity" %in% names(result_no_act),
  info = "prioritise - activity column absent when weight is zero"
)

expect_true(
  current = "scaled_p_val_adapted_ligand" %in% names(result_no_act),
  info = "prioritise - DE column present when weight non-zero"
)

#### condition DE branch -------------------------------------------------------

condition_de <- data.table::data.table(
  gene = c("L1", "L2", "L3", "L4", "REC1", "REC2", "REC3"),
  lfc = c(2.0, 0, 0, 0, 1.8, 0, 0),
  pval = c(1e-6, 0.5, 0.6, 0.7, 1e-5, 0.5, 0.5)
)

result_case <- prioritise_interactions(
  celltype_de = celltype_de,
  expression_info = expression_info,
  ligand_activities = ligand_activities,
  lr_network = lr_network,
  senders_oi = senders,
  receivers_oi = receivers,
  condition_de = condition_de,
  scenario = "case_control"
)

expect_true(
  current = "scaled_p_val_adapted_ligand_group" %in%
    names(result_case) &&
    "scaled_p_val_adapted_receptor_group" %in% names(result_case),
  info = "prioritise - condition columns present in case_control"
)

expect_error(
  current = prioritise_interactions(
    celltype_de = celltype_de,
    expression_info = expression_info,
    ligand_activities = ligand_activities,
    lr_network = lr_network,
    senders_oi = senders,
    receivers_oi = receivers,
    condition_de = NULL,
    scenario = "case_control"
  ),
  info = "prioritise - missing condition_de with case_control errors"
)

#### duplicate ligand activity -------------------------------------------------

bad_la <- rbind(ligand_activities, ligand_activities)

expect_error(
  current = prioritise_interactions(
    celltype_de = celltype_de,
    expression_info = expression_info,
    ligand_activities = bad_la,
    lr_network = lr_network,
    senders_oi = senders,
    receivers_oi = receivers,
    scenario = "one_condition"
  ),
  info = "prioritise - duplicate ligand in activity errors"
)

### end-to-end: ltf + activity + prioritisation --------------------------------

set.seed(123L)

# Use the inf_obj from the ltf section. Wire L1, L2, L3 as ligands available
# in two senders (X, Y), and a single receiver Z whose DEGs match TF1 targets.
# Expected: top hit is (X, L1, Z, REC_L1) because L1 has the highest activity
# (matches the receiver gene set), DE signal, and expression signal.

e2e_lr <- data.table::data.table(
  ligand = c("L1", "L2", "L3"),
  receptor = c("REC_L1", "REC_L2", "REC_L3")
)

e2e_activity <- ligand_activity_scores(
  ligand_influence = inf_obj,
  gene_sets = list(receiver_degs = tf1_targets)
)
# subset to just the LR ligands and drop the gene_set column
e2e_activity <- e2e_activity[
  ligand %in% e2e_lr$ligand,
  .(ligand, aupr_corrected)
]

e2e_celltype_de <- data.table::CJ(
  cluster_id = c("X", "Y", "Z"),
  gene = c("L1", "L2", "L3", "REC_L1", "REC_L2", "REC_L3")
)
e2e_celltype_de[, lfc := rnorm(.N, 0, 0.2)]
e2e_celltype_de[, pval := runif(.N, 0.1, 1)]
# Implant: X strongly upregulates L1, Z strongly upregulates REC_L1
e2e_celltype_de[cluster_id == "X" & gene == "L1", `:=`(lfc = 3.0, pval = 1e-10)]
e2e_celltype_de[
  cluster_id == "Z" & gene == "REC_L1",
  `:=`(lfc = 2.5, pval = 1e-9)
]

e2e_expr <- data.table::CJ(
  cluster_id = c("X", "Y", "Z"),
  gene = c("L1", "L2", "L3", "REC_L1", "REC_L2", "REC_L3")
)
e2e_expr[, avg_expr := runif(.N, 0, 0.1)]
e2e_expr[cluster_id == "X" & gene == "L1", avg_expr := 5.0]
e2e_expr[cluster_id == "Z" & gene == "REC_L1", avg_expr := 5.0]

e2e_result <- prioritise_interactions(
  celltype_de = e2e_celltype_de,
  expression_info = e2e_expr,
  ligand_activities = e2e_activity,
  lr_network = e2e_lr,
  senders_oi = c("X", "Y"),
  receivers_oi = "Z",
  scenario = "one_condition"
)

e2e_top <- e2e_result[1, .(sender, ligand, receptor)]

expect_true(
  current = e2e_top$sender == "X" &&
    e2e_top$ligand == "L1" &&
    e2e_top$receptor == "REC_L1",
  info = "e2e - full pipeline recovers implanted ligand-receptor signal"
)
