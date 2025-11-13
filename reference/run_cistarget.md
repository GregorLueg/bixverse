# Main function to run CisTarget

The `bixverse` implementation of the RCisTarget workflow, one of the
algorithms used in SCENIC, see Aibar, et al. You will need motif to
target gene rankings, see
[`read_motif_ranking()`](read_motif_ranking.md) and the motif to TF
annotations, see
[`read_motif_annotation_file()`](read_motif_annotation_file.md).

## Usage

``` r
run_cistarget(
  gs_list,
  rankings,
  annot_data,
  cis_target_params = params_cistarget()
)
```

## Arguments

- gs_list:

  Named list of character vectors. Each element is a gene set containing
  gene identifiers that must match row names in `rankings`.

- rankings:

  Integer matrix. Motif rankings for genes. Row names are gene
  identifiers, column names are motif identifiers. Lower values indicate
  higher regulatory potential.

- annot_data:

  data.table. Motif annotation database mapping motifs to transcription
  factors. Must contain columns: `motif`, `TF`, and `annotationSource`.

- cis_target_params:

  List. Output of [`params_cistarget()`](params_cistarget.md):

  - auc_threshold - Numeric. Proportion of genes to use for AUC
    threshold calculation. Default 0.05 means top 5 percent of genes.

  - nes_threshold - Numeric. Normalised Enrichment Score threshold for
    determining significant motifs. Default is 3.0.

  - rcc_method - Character. Recovery curve calculation method: "approx"
    (approximate, faster) or "icistarget" (exact, slower).

  - high_conf_cats - Character vector. Annotation categories considered
    high confidence (e.g., "directAnnotation", "inferredBy_Orthology").

  - low_conf_cats - Character vector. Annotation categories considered
    lower confidence (e.g., "inferredBy_MotifSimilarity").

## Value

data.table with enriched motifs and corresponding statistics and high &
low confidence TFs for each gene set.

## References

Aibar, et al., Nat Methods, 2017
