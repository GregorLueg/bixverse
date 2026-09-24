# SCENIC on MetaCells

**\[experimental\]** Assumes that the sparse data is pre-filtered for
the cells and genes you wish to include: every column is a target gene.
The regressors read the second layer, which here is an `f32` cast of the
supplied counts.

## Usage

``` r
rs_mc_scenic(sparse_data, tf_indices, scenic_params, seed, verbose)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `cs_type`. Shape is (metacells, genes).

- tf_indices:

  Integer vector. 0-indexed(!) column positions of the transcription
  factors within `sparse_data`.

- scenic_params:

  Named list. Contains all of the parameters needed for SCENIC, see
  [`params_scenic()`](https://gregorlueg.github.io/bixverse/reference/params_scenic.md).

- seed:

  Integer. Controls reproducibility of the function.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A genes x TFs importance matrix, rows in column order of `sparse_data`,
columns in the order of `tf_indices`.
