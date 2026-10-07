# Check a ligand-receptor database

Checkmate extension for checking a ligand-receptor database as returned
by
[`get_cellphonedb_db()`](https://gregorlueg.github.io/bixverse/reference/get_cellphonedb_db.md).

## Usage

``` r
checkLrDb(x)
```

## Arguments

- x:

  The data.table to check. Needs `interaction_id`, `partner_a` and
  `partner_b` as character columns and `genes_a` and `genes_b` as list
  columns of non-empty character vectors (the subunits).

## Value

`TRUE` if the check was successful, otherwise an error message.
