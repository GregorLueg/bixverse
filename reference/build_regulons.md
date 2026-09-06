# Build the final regulons

Turns the TF to gene table into the gene sets you hand to
[`aucell_sc()`](https://gregorlueg.github.io/bixverse/reference/aucell_sc.md).
[`tf_to_genes_motif_enrichment()`](https://gregorlueg.github.io/bixverse/reference/tf_to_genes_motif_enrichment.md)
only flags which pairs survived the CisTarget leading edge, it does not
remove anything, so this is the step that applies that filter.

Following SCENIC, the TF is added back to its own target list and
regulons below `min_genes` are dropped. Note the CisTarget leading edge
usually cuts hard: a module of a few hundred candidate targets typically
ends up as a regulon of a few dozen.

## Usage

``` r
build_regulons(
  x,
  use_leading_edge = TRUE,
  add_tf = TRUE,
  min_genes = 10L,
  mode = c("activating", "repressing", "both"),
  .verbose = TRUE
)

# S3 method for class 'ScenicGrn'
build_regulons(
  x,
  use_leading_edge = TRUE,
  add_tf = TRUE,
  min_genes = 10L,
  mode = c("activating", "repressing", "both"),
  .verbose = TRUE
)
```

## Arguments

- x:

  `ScenicGrn` object. You need to have run
  [`tf_to_genes_motif_enrichment()`](https://gregorlueg.github.io/bixverse/reference/tf_to_genes_motif_enrichment.md)
  for the leading edge filter to do anything.

- use_leading_edge:

  Boolean. Shall only the TF to gene pairs inside the CisTarget leading
  edge be kept. Defaults to `TRUE`.

- add_tf:

  Boolean. Shall the TF be added to its own regulon. Defaults to `TRUE`,
  following SCENIC.

- min_genes:

  Integer. Regulons with fewer genes than this are dropped. Defaults to
  `10L`, the SCENIC value.

- mode:

  String. Which links to keep, based on the `cor_sign` column added by
  [`tf_to_genes_correlations()`](https://gregorlueg.github.io/bixverse/reference/tf_to_genes_correlations.md).
  One of `c("activating", "repressing", "both")`. With `"both"` the
  regulon names get a `"_pos"` or `"_neg"` suffix so the two stay
  distinct. Defaults to `"activating"`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

A named list of character vectors, one per regulon.

## References

Aibar, et al., Nat Methods, 2017

## Examples

``` r
# turn the signed TF to gene table into regulon gene sets
sc <- demo_single_cells()
grn <- scenic_grn_sc(
  sc,
  tf_ids = sprintf("gene_%02d", 1:5),
  scenic_params = params_scenic(
    min_counts = 1L,
    learner_params = list(n_trees = 20L)
  ),
  .verbose = FALSE
)
grn <- identify_tf_to_genes(
  grn,
  method = "top_k",
  k_tfs = 3L,
  .verbose = FALSE
)
grn <- tf_to_genes_correlations(grn, object = sc, .verbose = FALSE)
lengths(build_regulons(
  grn,
  use_leading_edge = FALSE,
  min_genes = 3L,
  .verbose = FALSE
))
#> gene_01 gene_02 gene_03 gene_04 gene_05 
#>      10       4      17      13       5 

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
