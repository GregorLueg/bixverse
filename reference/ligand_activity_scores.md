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
