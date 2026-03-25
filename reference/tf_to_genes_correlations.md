# Generate TF to gene correlations

This function will calculate the correlations between the identified TF
to gene pairs. You need to have run
[`identify_tf_to_genes()`](https://gregorlueg.github.io/bixverse/reference/identify_tf_to_genes.md)!

## Usage

``` r
tf_to_genes_correlations(
  x,
  object,
  cor_filter = NULL,
  remove_self = TRUE,
  spearman = TRUE,
  .verbose = TRUE
)

# S3 method for class 'ScenicGrn'
tf_to_genes_correlations(
  x,
  object,
  cor_filter = NULL,
  remove_self = TRUE,
  spearman = TRUE,
  .verbose = TRUE
)
```

## Arguments

- x:

  `ScenicGrn` object for which to generate the TF to gene associations.

- object:

  `SingleCells` object that was used to generate the original GRNs.

- cor_filter:

  Optional float. If you wish to filter out TF genes below a certain
  correlation. If `NULL` all genes will be kept.

- remove_self:

  Boolean. Shall self loops (where TF controls its own expression) be
  removed. Defaults to `TRUE`.

- spearman:

  Boolean. Shall Spearman correlation be used. Defaults to `TRUE`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

Adds the correlations coefficients between to the TF to gene.
