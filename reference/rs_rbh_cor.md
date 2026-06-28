# Generate reciprocal best hits based on correlations

**\[experimental\]** This function takes list of (named) matrices which
represent for example matrix factorisation results you wish to identify
(k-th) reciprocal best hits (RBH) for. The rows need to represent the
features and the columns the parts you wish to calculate the RBH for.

## Usage

``` r
rs_rbh_cor(module_matrices, k_best, spearman, min_similarity)
```

## Arguments

- module_matrices:

  A list of named matrices. Rows represent features and columns the
  samples you wish to calculate the correlations for.

- k_best:

  Integer. Number of best neighbours to consider. If set to `1L`, this
  behaves as the traditional reciprocal best hit. If you set this to
  `3L` you consider edges if the modules is in the top 3 best modules by
  similarity for each other.

- spearman:

  Shall Spearman correlation be used.

- min_similarity:

  Minimum (absolute) correlations that needs to exist between two terms.

## Value

A list containing:

- origin - The name of the origin of the gene modules.

- target - The name of the target of the gene modules.

- comparisons - Integer vector indicating how many RBH hits were
  identified in this comparison

- origin_modules - Names of the gene modules from the origin.

- target_modules - Names of the gene modules from the target.

- similarity - The similarities between the two respective gene modules.
