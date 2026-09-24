# Run NMF (HALS) on MetaCells

**\[experimental\]** Assumes that the sparse data is pre-filtered for
the cells/genes you wish to include. Indices in the sparse data need to
be 0-indexed. Both data layers hold the supplied values, so which assay
NMF runs on is decided by what is passed in, not by `use_second_layer`.

## Usage

``` r
rs_nmf_single_mc(
  sparse_data,
  k,
  preprocessing,
  use_second_layer,
  nmf_hals_params,
  seed,
  verbose
)
```

## Arguments

- sparse_data:

  A named list with `data`, `indptr`, `indices`, `nrow`, `ncol` and
  `cs_type`. Shape is (metacells, genes).

- k:

  Integer. Number of latent factors to return.

- preprocessing:

  String. One of `c("none", "sd", "sqrt_sd")`.

- use_second_layer:

  Boolean. Shall the second data layer be used.

- nmf_hals_params:

  Named list. Contains the NMF parameters, see
  [`params_nmf_hals()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_hals.md).

- seed:

  Integer. Random seed for initialisation.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items

- w - The `W` matrix of shape `n_meta_cells x k`.

- h - The `H` matrix of shape `k x n_genes`.

- final_loss - Final squared Frobenius reconstruction loss.

- n_iter - Number of iterations the algorithm ran for.

- converged - Did the NMF algorithm converge.
