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

  A list of logicals with the genes of interest being set to `TRUE` and
  the background genes set to `FALSE`.

## Value

A list with internal lists with:

- `auroc` - The Area Under the Receiver Operating Characteristic for
  that ligand

- `aupr` - The Area Under the Precision-Recall curve for that ligand.

- `aupr_corrected` - The corrected AUPR

- `pearson` - The Pearson correlations

- `spearman` - The Spearman correlations
