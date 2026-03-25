# Filter genes for SCENIC GRN inference

Filters genes by minimum total counts and minimum expressed-cell
fraction using the SCENIC inclusion criteria. Returns a character vector
of gene identifiers passing both filters.

## Usage

``` r
scenic_gene_filter_sc(
  object,
  scenic_params = params_scenic(),
  cells_to_take = NULL,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- scenic_params:

  List. SCENIC parameters, see
  [`params_scenic()`](https://gregorlueg.github.io/bixverse/reference/params_scenic.md).
  Only `min_counts` and `min_cells` are used by this function.

- cells_to_take:

  Optional string vector. Cell identifiers to restrict to. If `NULL`,
  defaults to all filtered cells in the class.

- .verbose:

  Boolean. Controls verbosity. Defaults to `TRUE`.

## Value

A character vector of gene identifiers passing the SCENIC inclusion
criteria.
