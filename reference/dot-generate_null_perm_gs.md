# Generate random gene sets for VISION p-value calculations

This function will generate random gene sets given the provided gs_list.
Under the hood, it uses the same approach as in DeTomaso, et al. and
does not generate a random signature per given signature, but `n_comp`
representative ones based on size and balance (positive and negative
genes) via k-means clustering. The original gene sets are than matched
to the closest cluster. The authors observed that this sufficed to
estimate significance, see DeTomaso, et al.

## Usage

``` r
.generate_null_perm_gs(
  gs_list,
  expr_genes,
  n_perm = 500L,
  n_comp = 5L,
  random_seed = 42L
)
```

## Arguments

- gs_list:

  Named nested list for which to calculate the local auto-correlations.
  The elements have the gene identifiers of the respective gene sets and
  have the option to have a `"pos"` and `"neg"` gene sets. The names
  need to be part of the variables of the `SingleCells` class.

- expr_genes:

  Character vector. Represents the genes expressed in the experiment.

- n_perm:

  Integer. Number of random permutations to generate.

- n_comp:

  Integer. Number of k-means cluster to identify.

- random_seed:

  Integer. For reproducibility purposes.

## Value

A list with the following elements:

- random_signatures - Nested list representing the random permutations.

- clusters - Association of original gene set to random permutation set.
