# Generate an miloR abundance object for differential abundance testing

This function implements the miloR differential abundance testing on top
of the kNN graph. The general idea of the approach is to use the kNN
graph generated from the single cell data, generate representative
neighbourhoods and calculate differential abundances within these
neighbourhoods. For further details on the method, please refer to Dann,
et al. This function will take a `SingleCells` class, run the
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

  `SingleCells` class.

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

  A list, please see
  [`params_sc_miloR()`](https://gregorlueg.github.io/bixverse/reference/params_sc_miloR.md).
  The list has the following parameters:

  - prop - Numeric. Proportion of cells to sample as neighbourhood
    indices. Must be in (0,1).

  - k_refine - Integer. Number of neighbours to use for refinement.

  - refinement_strategy - String. Strategy for refining sampled indices.
    One of `c("approximate", "bruteforce", "index")`.

  - index_type - String. Type of kNN index to use. One of
    `c("annoy", "hnsw")`.

  - knn - List of kNN parameters. See
    [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
    for available parameters and their defaults. Note: `knn_method`
    cannot be `"exhaustive"` for MiloR as it basically boils down to
    `"bruteforce"`.

- seed:

  Integer. Seed for reproducibility

- .verbose:

  Boolean. Controls verbosity of the method.

## References

Dann, et al., Nat Biotechnol, 2022
