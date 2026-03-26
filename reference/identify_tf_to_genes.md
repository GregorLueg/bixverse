# Identify the TF to gene regulation

Generates the TF to gene associations from the importance matrix
produced by the tree-based regression models. Two filtering strategies
are available:

- **Threshold** (`method = "threshold"`): For each gene, computes mean +
  `n_sd` \* SD of the importance scores across all TFs and retains only
  pairs exceeding that threshold. This adapts to the per-gene importance
  distribution and is less sensitive to differences between learners.

- **Top-k** (`method = "top_k"`): Selects the top `k_tfs` TFs per gene
  and/or the top `k_genes` genes per TF. Both margins can be combined
  (union). At least one of `k_tfs` or `k_genes` must be provided.

Both methods accept an optional `min_importance` floor.

## Usage

``` r
identify_tf_to_genes(
  x,
  method = c("threshold", "top_k"),
  k_tfs = NULL,
  k_genes = NULL,
  n_sd = 2,
  min_importance = NULL,
  .verbose = TRUE
)

# S3 method for class 'ScenicGrn'
identify_tf_to_genes(
  x,
  method = c("threshold", "top_k"),
  k_tfs = NULL,
  k_genes = NULL,
  n_sd = 2,
  min_importance = NULL,
  .verbose = TRUE
)
```

## Arguments

- x:

  `ScenicGrn` object.

- method:

  Character. Either `"top_k"` or `"threshold"`.

- k_tfs:

  Optional integer. Top TFs per gene (only used when
  `method = "top_k"`).

- k_genes:

  Optional integer. Top genes per TF (only used when
  `method = "top_k"`).

- n_sd:

  Numeric. Number of standard deviations above the per-gene mean to use
  as the threshold (only used when `method = "threshold"`). Default is
  `2`.

- min_importance:

  Optional numeric in \[0, 1\]. Absolute minimum importance score for
  inclusion.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The `ScenicGrn` object with a `tf_to_gene_results` data.table added.
