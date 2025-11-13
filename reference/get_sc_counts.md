# Getter the counts

Getter the counts

## Usage

``` r
get_sc_counts(
  object,
  assay = c("raw", "norm"),
  return_format = c("cell", "gene"),
  cell_indices = NULL,
  gene_indices = NULL,
  use_cells_to_keep = TRUE,
  .verbose = TRUE
)
```

## Arguments

- object:

  `single_cell_exp`, `meta_cells` class.

- assay:

  String. Which slot to return. One of `c("raw", "norm")`. Defaults to
  `"raw"`.

- return_format:

  String. One of `c("cell", "gene")`. Return data in cell-centric
  compressed format (CSR) or gene-centric compressed format (CSC).
  Defaults to `"cell"`. Not relevant for `meta_cells`.

- cell_indices:

  Optional cell indices.

- gene_indices:

  Optional gene indices.

- use_cells_to_keep:

  Boolean. Shall cells to keep be found in the class, shall the counts
  be reduced to these. Not relevant for `meta_cells`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The counts table
