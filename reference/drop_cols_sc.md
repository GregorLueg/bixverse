# Drop columns from the obs or var table

Drops the named columns from the obs or var table of single cell-related
classes. Protected identifier and bookkeeping columns (`cell_idx`,
`cell_id`, `to_keep` for obs; `gene_idx`, `gene_id` for var) are refused
with a warning. Columns that do not exist trigger a warning and are
skipped.

## Usage

``` r
drop_cols_sc(object, table = c("obs", "var"), cols)
```

## Arguments

- object:

  `SingleCells` (or other compatible) class.

- table:

  String. One of `c("obs", "var")`.

- cols:

  Character vector. Column names to drop.

## Value

Invisible self.
