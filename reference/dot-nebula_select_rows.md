# Restrict a NEBULA obs table to selected cells

Shared by
[`nebula_sc()`](https://gregorlueg.github.io/bixverse/reference/nebula_sc.md),
[`nebula_mc()`](https://gregorlueg.github.io/bixverse/reference/nebula_mc.md)
and the GPU counterpart. Ids that are not in the table, e.g. cells that
failed quality control, are reported and ignored.

## Usage

``` r
.nebula_select_rows(obs, id_col, ids)
```

## Arguments

- obs:

  data.table. The obs table of the object.

- id_col:

  String. The identifier column to match on, `"cell_id"` or
  `"meta_cell_id"`.

- ids:

  Character vector. The identifiers to keep.

## Value

`obs` restricted to the rows in `ids`, in its original order.
