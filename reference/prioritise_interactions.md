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

## Examples

``` r
# one sender, one receiver, with the signal planted on L1 -> R1
lr_network <- data.table::data.table(
  ligand = c("L1", "L2"),
  receptor = c("R1", "R2")
)
celltype_de <- data.table::CJ(
  cluster_id = c("sender", "receiver"),
  gene = c("L1", "L2", "R1", "R2")
)
celltype_de[, lfc := c(3, 0.1, 0.1, 0.1, 0.1, 0.1, 2.5, 0.1)]
#> Key: <cluster_id, gene>
#>    cluster_id   gene   lfc
#>        <char> <char> <num>
#> 1:   receiver     L1   3.0
#> 2:   receiver     L2   0.1
#> 3:   receiver     R1   0.1
#> 4:   receiver     R2   0.1
#> 5:     sender     L1   0.1
#> 6:     sender     L2   0.1
#> 7:     sender     R1   2.5
#> 8:     sender     R2   0.1
celltype_de[, pval := c(1e-10, 0.5, 0.5, 0.5, 0.5, 0.5, 1e-9, 0.5)]
#> Key: <cluster_id, gene>
#>    cluster_id   gene   lfc  pval
#>        <char> <char> <num> <num>
#> 1:   receiver     L1   3.0 1e-10
#> 2:   receiver     L2   0.1 5e-01
#> 3:   receiver     R1   0.1 5e-01
#> 4:   receiver     R2   0.1 5e-01
#> 5:     sender     L1   0.1 5e-01
#> 6:     sender     L2   0.1 5e-01
#> 7:     sender     R1   2.5 1e-09
#> 8:     sender     R2   0.1 5e-01
expression_info <- celltype_de[, .(cluster_id, gene, avg_expr = 0.05)]
expression_info[cluster_id == "sender" & gene == "L1", avg_expr := 5]
#> Key: <cluster_id, gene>
#>    cluster_id   gene avg_expr
#>        <char> <char>    <num>
#> 1:   receiver     L1     0.05
#> 2:   receiver     L2     0.05
#> 3:   receiver     R1     0.05
#> 4:   receiver     R2     0.05
#> 5:     sender     L1     5.00
#> 6:     sender     L2     0.05
#> 7:     sender     R1     0.05
#> 8:     sender     R2     0.05
expression_info[cluster_id == "receiver" & gene == "R1", avg_expr := 5]
#> Key: <cluster_id, gene>
#>    cluster_id   gene avg_expr
#>        <char> <char>    <num>
#> 1:   receiver     L1     0.05
#> 2:   receiver     L2     0.05
#> 3:   receiver     R1     5.00
#> 4:   receiver     R2     0.05
#> 5:     sender     L1     5.00
#> 6:     sender     L2     0.05
#> 7:     sender     R1     0.05
#> 8:     sender     R2     0.05
res <- prioritise_interactions(
  celltype_de = celltype_de,
  expression_info = expression_info,
  ligand_activities = data.table::data.table(
    ligand = c("L1", "L2"),
    aupr_corrected = c(0.75, -0.125)
  ),
  lr_network = lr_network,
  senders_oi = "sender",
  receivers_oi = "receiver",
  scenario = "one_condition"
)
head(res[, .(sender, ligand, receiver, receptor, prioritisation_score)])
#>    sender ligand receiver receptor prioritisation_score
#>    <char> <char>   <char>   <char>                <num>
#> 1: sender     L1 receiver       R1            0.8333333
#> 2: sender     L2 receiver       R2            0.5000000
```
