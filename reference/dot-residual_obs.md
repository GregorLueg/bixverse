# Observation rows for the selected cells, in the selected order

The group labels and the covariates are matched to cells by position, so
the observation rows have to arrive in the same order as `cell_indices`.
[`get_sc_obs()`](https://gregorlueg.github.io/bixverse/reference/get_sc_obs.md)
reads from DuckDB without an `ORDER BY`, so the order it returns is not
something to rely on. Reordering here is cheap and removes a failure
that would otherwise be silent: a model fitted against shuffled
covariates gives plausible residuals and a wrong embedding.

## Usage

``` r
.residual_obs(object, cell_indices)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- cell_indices:

  Integer. The 0-based cells, in the order Rust gets them.

## Value

The observation table for those cells, row `i` being `cell_indices[i]`.
