# Run the SCENIC motif enrichment

This function will run the motif enrichment based on RCisTarget (or in
this case the internal implementation via
[`run_cistarget()`](https://gregorlueg.github.io/bixverse/reference/run_cistarget.md)).
You need to provide the expected rankings and TF annotations (for
details, please see
[`run_cistarget()`](https://gregorlueg.github.io/bixverse/reference/run_cistarget.md)).
Briefly, this function will run CisTarget and add the results to the
`ScenicGrn` object and add additionally a column `"in_motif"`, for a
given TF to gene set to say if it was part of the motifs associated with
this TF (or not). You have the option to limit this to only the the high
confidence TFs (default), or also include the low confidence TFs (i.e.,
links from TF to motif that are less certain).

## Usage

``` r
tf_to_genes_motif_enrichment(
  x,
  motif_rankings,
  annot_data,
  cis_target_params = params_cistarget(),
  gene_id_to_symbol = NULL,
  only_high_conf_tf = TRUE,
  .verbose = TRUE
)

# S3 method for class 'ScenicGrn'
tf_to_genes_motif_enrichment(
  x,
  motif_rankings,
  annot_data,
  cis_target_params = params_cistarget(),
  gene_id_to_symbol = NULL,
  only_high_conf_tf = TRUE,
  .verbose = TRUE
)
```

## Arguments

- x:

  `ScenicGrn` object for which to generate the TF to gene associations.

- motif_rankings:

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

- gene_id_to_symbol:

  Named character vector. Mapping from gene identifiers used internally
  (e.g., Ensembl IDs) to the identifiers used in the motif rankings
  (e.g., HGNC symbols). Names are internal IDs, values are ranking IDs.
  If `NULL` (default), no mapping is applied and the internal IDs are
  assumed to match the ranking rownames.

- only_high_conf_tf:

  Boolean. Shall only the high confidence TF to motif association be
  used. Defaults to `TRUE`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

Adds a data.table with the first tf to gene results to the class.
