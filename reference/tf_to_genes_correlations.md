# Generate TF to gene correlations

This function will calculate the correlations between the identified TF
to gene pairs. You need to have run
[`identify_tf_to_genes()`](https://gregorlueg.github.io/bixverse/reference/identify_tf_to_genes.md)!

Following SCENIC, the correlation is turned into a three-level sign at
`rho_threshold`: `+1` for activating links, `-1` for repressing ones and
`0` for everything in the band between. `mode` then decides which of
those you keep. The default keeps the activating links only, which is
what SCENIC does with `onlyPositiveCorr = TRUE`.

## Usage

``` r
tf_to_genes_correlations(
  x,
  object,
  rho_threshold = 0.03,
  mode = c("activating", "repressing", "both"),
  remove_self = TRUE,
  spearman = TRUE,
  cor_filter = NULL,
  .verbose = TRUE
)

# S3 method for class 'ScenicGrn'
tf_to_genes_correlations(
  x,
  object,
  rho_threshold = 0.03,
  mode = c("activating", "repressing", "both"),
  remove_self = TRUE,
  spearman = TRUE,
  cor_filter = NULL,
  .verbose = TRUE
)
```

## Arguments

- x:

  `ScenicGrn` object for which to generate the TF to gene associations.

- object:

  `SingleCells` or `MetaCells` object that was used to generate the
  original GRNs.

- rho_threshold:

  Float. Absolute correlation above which a TF to gene link counts as
  activating or repressing. Defaults to `0.03`, the SCENIC value.

- mode:

  String. Which links to keep. One of
  `c("activating", "repressing", "both")`. Defaults to `"activating"`.

- remove_self:

  Boolean. Shall self loops (where TF controls its own expression) be
  removed. Defaults to `TRUE`. Note that
  [`build_regulons()`](https://gregorlueg.github.io/bixverse/reference/build_regulons.md)
  adds the TF back to its own regulon, which is also what SCENIC does.

- spearman:

  Boolean. Shall Spearman correlation be used. Defaults to `TRUE`.

- cor_filter:

  Deprecated. Use `rho_threshold` and `mode` instead. If given, it is
  treated as a one-sided lower bound on the correlation and
  `rho_threshold` and `mode` are ignored.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

Adds a `pairwise_cor` and a `cor_sign` column to the TF to gene results.

## References

Aibar, et al., Nat Methods, 2017

## Examples

``` r
# sign the TF to gene links and keep the activating ones
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
head(get_tf_to_gene(grn))
#>         tf    gene importance pairwise_cor cor_sign
#>     <char>  <char>      <num>        <num>    <int>
#> 1: gene_03 gene_01  0.3568220    0.5735584        1
#> 2: gene_04 gene_01  0.1553617    0.5115550        1
#> 3: gene_03 gene_02  0.2651168    0.5171022        1
#> 4: gene_04 gene_02  0.1317756    0.4859316        1
#> 5: gene_04 gene_03  0.1788230    0.5314755        1
#> 6: gene_01 gene_03  0.1435614    0.5735584        1

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
