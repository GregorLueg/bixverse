# Generate an miloR abundance object for differential abundance testing

This function implements the miloR differential abundance testing on top
of the kNN graph. The general idea of the approach is to use the kNN
graph generated from the single cell data, generate representative
neighbourhoods and calculate differential abundances within these
neighbourhoods. For further details on the method, please refer to Dann,
et al. This function will take an `single_cell_exp` class, run the
neighbourhood detection, count the occurrences of a sample and return a
`sc_miloR` class for subsequent differential abundance testing and
further annotations.

## Usage

``` r
get_miloR_abundances_sc(
  object,
  sample_id_col,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  miloR_params = params_sc_miloR(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `single_cell_exp` class.

- sample_id_col:

  Character. The column in the obs table representing the sample
  identifier to count.

- embd_to_use:

  Character. The embedding to use for the refinement procedure. Please
  use the same here as you used to generate the neighbours! Defaults to
  `"pca"`.

- no_embd_to_use:

  Optional integer. If you only want to use a subset of the embedding.

- miloR_params:

  A list, please see `params_sc_milor()`. The list has the following
  parameters:

  - prop - Numeric. Proportion of cells to sample as neighbourhood
    indices. Must be in (0,1). Defaults to `0.2`.

  - k_refine - Integer. Number of neighbours to use for refinement.
    Defaults to `20L`.

  - refinement_strategy - String. Strategy for refining sampled indices.
    One of `c("approximate", "bruteforce", "index")`. Defaults to
    `"index"`.

  - index_type - String. Type of kNN index to use. One of
    `c("annoy", "hnsw")`. Defaults to `"annoy"`.

  - k - Integer. Number of neighbours to consider. Defaults to `15L`.

  - knn_method - String. One of `c("annoy", "hnsw")`. The method to use
    for the approximate nearest neighbour search. Defaults to `"annoy"`.
    Note: `"nndescent"` is not supported for MiloR!

  - ann_dist - String. One of `c("cosine", "euclidean")`. The distance
    metric to be used for the approximate neighbour search. Defaults to
    `"euclidean"`.

  - search_budget - Integer. Search budget per tree for Annoy. Defaults
    to `100L`.

  - n_trees - Integer. Number of trees to generate for Annoy. Defaults
    to `100L`.

  - nn_max_iter - Integer. Maximum iterations for NNDescent. Defaults to
    `15L`.

  - rho - Numeric. Sampling rate for NNDescent. Defaults to `1.0`.

  - delta - Numeric. Early termination criterion for NNDescent. Defaults
    to `0.001`.

- seed:

  Integer. Seed for reproducibility

- .verbose:

  Boolean. Controls verbosity of the method.

## References

Dann, et al., Nat Biotechnol, 2022
