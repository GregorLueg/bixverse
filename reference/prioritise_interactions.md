# Prioritise sender-ligand-receiver-receptor interactions

For every (sender, ligand, receiver, receptor) tuple drawn from the LR
network, compute a weighted prioritisation score combining:

- ligand DE in the sender (`de_ligand`)

- receptor DE in the receiver (`de_receptor`)

- ligand activity in the receiver (`activity_scaled`)

- ligand expression specificity across senders (`exprs_ligand`)

- receptor expression specificity across receivers (`exprs_receptor`)

- ligand condition specificity (`ligand_condition_specificity`)

- receptor condition specificity (`receptor_condition_specificity`)

Components with weight zero are skipped (their joins are not performed).

DE tables are supplied by the user. `find_markers_sc` and
`find_all_markers_sc` give a quick Wilcox; pseudo-bulk + DESeq2 / edgeR
is generally preferable for statistical inference.

## Usage

``` r
prioritise_interactions(
  celltype_de,
  expression_info,
  ligand_activities,
  lr_network,
  senders_oi,
  receivers_oi,
  condition_de = NULL,
  weights = NULL,
  scenario = c("case_control", "one_condition")
)
```

## Arguments

- celltype_de:

  A `data.table` with columns `cluster_id`, `gene`, `lfc`, `pval`. One
  row per (cluster, gene). Must cover both senders and receivers.

- expression_info:

  A `data.table` with columns `cluster_id`, `gene`, `avg_expr`. As
  returned by `compute_expression_info_sc`.

- ligand_activities:

  A `data.table` with at least the columns `ligand` and
  `aupr_corrected`. One row per ligand. Subset to a single gene set
  before calling.

- lr_network:

  A `data.table` with columns `ligand`, `receptor`.

- senders_oi:

  Character vector of sender cluster IDs.

- receivers_oi:

  Character vector of receiver cluster IDs.

- condition_de:

  Optional `data.table` with columns `gene`, `lfc`, `pval` from a
  condition contrast (e.g. case vs control). Required when
  `scenario == "case_control"` and the corresponding weights are
  non-zero.

- weights:

  Optional named numeric vector. If `NULL`, defaults are chosen by
  `scenario`. Names: `de_ligand`, `de_receptor`, `activity_scaled`,
  `exprs_ligand`, `exprs_receptor`, `ligand_condition_specificity`,
  `receptor_condition_specificity`.

- scenario:

  `"case_control"` (all weights 1) or `"one_condition"`
  (condition-specificity weights 0). Ignored if `weights` is supplied.

## Value

A `data.table` with one row per surviving (sender, ligand, receiver,
receptor) tuple, sorted by descending `prioritisation_score`, with a
`prioritisation_rank` column.
