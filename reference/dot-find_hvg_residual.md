# Residual-based HVG selection

Shared body of the residual branch in
[`find_hvg_sc()`](https://gregorlueg.github.io/bixverse/reference/find_hvg_sc.md).
The selection itself happens in Rust, so this resolves the fit, calls
over and writes the result back.

## Usage

``` r
.find_hvg_residual(
  object,
  hvg_no,
  write_var,
  gene_batch_size = NULL,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- hvg_no:

  Integer. Variable features per group.

- write_var:

  Boolean. Write the residual variance to the variable table.

- gene_batch_size:

  Integer or `NULL`. Genes held in memory per batch.

- .verbose:

  Boolean or Integer. Controls verbosity.

## Value

The object with the HVGs set.
