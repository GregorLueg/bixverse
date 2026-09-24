# Calculate the NicheNet ligand activity scores

**\[experimental\]**

## Usage

``` r
rs_ligand_activity_scores(ligand_influence, in_gene_sets)
```

## Arguments

- ligand_influence:

  A ligand x background genes matrix that measures the ligand to target
  gene influence.

- in_gene_sets:

  List of logical vectors, one per gene set, each of length
  `ncol(ligand_influence)`. Genes of interest are `TRUE`, the background
  genes `FALSE`.

## Value

A list with one element per gene set, each a list of per-ligand vectors
(`NaN` where the metric is undefined):

- `auroc` - The Area Under the Receiver Operating Characteristic.

- `aupr` - The Area Under the Precision-Recall curve.

- `aupr_corrected` - The corrected AUPR.

- `pearson` - The Pearson correlations.

- `spearman` - The Spearman correlations.
