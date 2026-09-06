# Main function to run CisTarget

The `bixverse` implementation of the RCisTarget workflow, one of the
algorithms used in SCENIC, see Aibar, et al. You will need motif to
target gene rankings, see
[`read_motif_ranking()`](https://gregorlueg.github.io/bixverse/reference/read_motif_ranking.md)
and the motif to TF annotations, see
[`read_motif_annotation_file()`](https://gregorlueg.github.io/bixverse/reference/read_motif_annotation_file.md).

## Usage

``` r
run_cistarget(
  gs_list,
  rankings,
  annot_data,
  cis_target_params = params_cistarget(),
  .verbose = TRUE
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

  List. Output of
  [`params_cistarget()`](https://gregorlueg.github.io/bixverse/reference/params_cistarget.md):

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

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

data.table with enriched motifs and corresponding statistics and high &
low confidence TFs for each gene set.

## References

Aibar, et al., Nat Methods, 2017

## Examples

``` r
# motif enrichment against a tiny synthetic ranking database
rankings <- matrix(
  c(1L, 5L, 4L, 2L, 2L, 2L, 1L, 5L, 4L, 3L, 2L, 1L, 3L, 1L, 5L, 4L,
    5L, 4L, 3L, 3L),
  nrow = 5,
  byrow = TRUE,
  dimnames = list(sprintf("gene_%i", 1:5), sprintf("motif_%i", 1:4))
)
annot <- data.table::data.table(
  motif = sprintf("motif_%i", 1:4),
  TF = sprintf("TF%i", 1:4),
  annotationSource = factor(c(
    "directAnnotation",
    "inferredBy_Orthology",
    "inferredBy_MotifSimilarity",
    "inferredBy_MotifSimilarity_n_Orthology"
  ))
)
res <- run_cistarget(
  gs_list = list(set_a = c("gene_1", "gene_2", "gene_3")),
  rankings = rankings,
  annot_data = annot,
  cis_target_params = params_cistarget(
    auc_threshold = 1,
    nes_threshold = 0.2
  ),
  .verbose = FALSE
)
res[, c("gs_name", "motif", "nes")]
#> Key: <motif>
#>    gs_name   motif       nes
#>     <char>  <char>     <num>
#> 1:   set_a motif_1 0.7071068
#> 2:   set_a motif_3 0.7071068
```
