# Generates random permutation vectors

Generates random permutation vectors

## Usage

``` r
generate_perm_diffusion_vecs(
  graph,
  diffusion_vec,
  bins = 25L,
  iters = 1000L,
  random_seed = 10101L
)
```

## Arguments

- graph:

  igraph. The graph for which to generate the random diffusion vectors

- diffusion_vec:

  Named numeric. The initial diffusion vector.

- bins:

  Integer. Number of bins to use for node degree aware sampling.

- iters:

  Integer. Number of random permutations to generate.

- random_seed:

  Integer. Random seed.

## Value

List with the permutations.
