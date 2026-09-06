# Compute ligand activity scores against gene sets

For each gene set, ranks ligands by how well their influence vector
aligns with set membership across a background. Returns AUROC, AUPR,
AUPR corrected against the random baseline, Pearson, and Spearman per
ligand and per gene set.

The influence matrix is restricted to `background` columns before
scoring, so AUROC / AUPR reflect ranking within the background (typical
NicheNet practice: background = genes expressed in the receiver cells,
gene set = the DEGs).

## Usage

``` r
ligand_activity_scores(ligand_influence, gene_sets, background = NULL)
```

## Arguments

- ligand_influence:

  A `LigandTargetInfluence` object.

- gene_sets:

  Either a character vector (one gene set) or a list of character
  vectors. Names of the list are propagated to the output.

- background:

  Character vector or `NULL`. Genes against which gene sets are scored.
  Defaults to all genes in `ligand_influence`. Background members not
  present in the influence matrix are silently dropped.

## Value

A `data.table` with one row per (gene set, ligand) pair and columns
`gene_set`, `ligand`, `auroc`, `aupr`, `aupr_corrected`, `pearson`,
`spearman`.

## Examples

``` r
# rank the ligands against the TF1 target set
ppi <- data.table::data.table(
  from = c("L1", "SIG1", "L2", "SIG2"),
  to = c("SIG1", "TF1", "SIG2", "TF2"),
  weight = 1.0
)
grn <- data.table::data.table(
  from = rep(c("TF1", "TF2"), each = 3),
  to = c("G1", "G2", "G3", "G4", "G5", "G6"),
  weight = 1.0
)
inf <- generate_ligand_target_influence(
  ligand_seeds = list(L1 = "L1", L2 = "L2"),
  ppi_network = ppi,
  grn_network = grn,
  params = params_ligand_target(ltf_cutoff = 0)
)
ligand_activity_scores(
  ligand_influence = inf,
  gene_sets = list(set_A = c("G1", "G2", "G3"))
)
#>    gene_set ligand     auroc  aupr aupr_corrected    pearson   spearman
#>      <char> <char>     <num> <num>          <num>      <num>      <num>
#> 1:    set_A     L1 1.0000000 1.000          0.750  1.0000000  1.0000000
#> 2:    set_A     L2 0.3333333 0.125         -0.125 -0.3333333 -0.3333333
```
