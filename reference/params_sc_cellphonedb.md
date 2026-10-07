# Parameters for the CellPhoneDB analysis

Defaults match CellPhoneDB v5: 1000 permutations and a 10% expression
threshold.

## Usage

``` r
params_sc_cellphonedb(
  n_perm = 1000L,
  threshold = 0.1,
  seed = 42L,
  perm_batch = NULL
)
```

## Arguments

- n_perm:

  Integer. Number of cluster label permutations for the p-values.
  Defaults to `1000L`.

- threshold:

  Numeric. Both partners need a fraction of expressing cells strictly
  above this value in their cluster. Defaults to `0.1`.

- seed:

  Integer. Seed for the permutations. Defaults to `42L`.

- perm_batch:

  Integer or `NULL`. Permutations that share one pass over the gene
  data. `NULL` uses the Rust default (16). Defaults to `NULL`.

## Value

A named list with the following elements:

- n_perm - Integer. Number of cluster label permutations for the
  p-values. Defaults to `1000L`.

- threshold - Numeric. Both partners need a fraction of expressing cells
  strictly above this value in their cluster. Defaults to `0.1`.

- seed - Integer. Seed for the permutations. Defaults to `42L`.

- perm_batch - Integer or `NULL`. Permutations that share one pass over
  the gene data. `NULL` uses the Rust default (16). Defaults to `NULL`.
