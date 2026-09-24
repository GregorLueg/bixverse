# Turn a corrected gene-major store into a `SingleCells`

The corrected counts come out gene-major only, and a `SingleCells` needs
the cell-major twin plus the database as well. The gene axis is the
model's, so the variable table is rebuilt from the genes that survived
rather than copied across: index `j` of the new store is a different
gene from index `j` of the old one.

## Usage

``` r
.build_corrected_object(
  object,
  dir_out,
  res,
  cell_indices,
  gene_batch_size,
  .verbose
)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- dir_out:

  String. Directory the store lives in.

- res:

  List. The result of
  [`rs_sct_corrected_counts()`](https://gregorlueg.github.io/bixverse/reference/rs_sct_corrected_counts.md).

- cell_indices:

  Integer. The 0-based cells that were written.

- gene_batch_size:

  Integer or `NULL`. Genes read per batch.

- .verbose:

  Boolean or Integer. Controls verbosity.

## Value

The new `SingleCells` over the corrected counts.
