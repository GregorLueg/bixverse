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
