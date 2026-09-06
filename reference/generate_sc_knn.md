# Generate a new SingleCellNearestNeighbour from data

Generate a new SingleCellNearestNeighbour from data

## Usage

``` r
generate_sc_knn(
  data,
  neighbours_params = params_sc_neighbours(),
  seed = 42L,
  .validate_index = FALSE,
  .verbose = TRUE
)
```

## Arguments

- data:

  Numerical matrix. Samplex x features. The embedding matrix from which
  to generate the kNN data.

- neighbours_params:

  List. Output of
  [`params_sc_neighbours()`](https://gregorlueg.github.io/bixverse/reference/params_sc_neighbours.md).
  A list with the following items:

  - full_snn - Boolean. Shall the full shared nearest neighbour graph be
    generated that generates edges between all cells instead of between
    only neighbours. (Not used in this function.)

  - pruning - Numeric. Weights below this threshold will be set to 0 in
    the generation of the sNN graph. (Not used in this function.)

  - snn_similarity - String. One of `c("rank", "jaccard")`. Defines how
    the weight from the SNN graph is calculated. For details, please see
    [`params_sc_neighbours()`](https://gregorlueg.github.io/bixverse/reference/params_sc_neighbours.md).
    (Not used in this function.)

  - knn - List of kNN parameters. See
    [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
    for available parameters and their defaults.

- seed:

  Integer. Random seed for reproducibility.

- .validate_index:

  Boolean. Shall an exhaustive search against a subset of cells be run
  to validate the approximate nearest neighbour index.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

The `SingleCellNearestNeighbour` for downstream usage.

## Examples

``` r
# kNN over a random embedding; the rows have to carry cell names
set.seed(42L)
embd <- matrix(rnorm(500 * 10), nrow = 500)
rownames(embd) <- sprintf("cell_%03d", 1:500)
knn <- generate_sc_knn(embd, .verbose = FALSE)
dim(get_knn_mat(knn))
#> [1] 500  15
```
