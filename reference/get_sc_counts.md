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
  modality = c("rna", "adt"),
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells`, `MetaCells`, `SingleCellsMultiModal` class.

- assay:

  String. Which slot to return. One of `c("raw", "norm")`. Defaults to
  `"raw"`.

- return_format:

  String. One of `c("cell", "gene")`. Return data in cell-centric
  compressed format (CSR) or gene-centric compressed format (CSC).
  Defaults to `"cell"`. Not relevant for `MetaCells`.

- cell_indices:

  Optional cell indices.

- gene_indices:

  Optional gene indices.

- use_cells_to_keep:

  Boolean. Shall cells to keep be found in the class, shall the counts
  be reduced to these. Not relevant for `MetaCells`.

- modality:

  String. The modality to return. One of `c("rna", "adt")`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The counts table

## Examples

``` r
# raw counts for the first ten genes, cell-centric (CSR)
sc <- demo_single_cells(prepped = FALSE)
counts <- get_sc_counts(sc, gene_indices = 1:10, .verbose = FALSE)
dim(counts)
#> [1] 500  10

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
