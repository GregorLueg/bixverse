# Assert that a grouping column exists in the obs table

Assert that a grouping column exists in the obs table

## Usage

``` r
.assert_group_by(object, group_by)
```

## Arguments

- object:

  A `SingleCells` object.

- group_by:

  Character. Column name to validate.

## Value

Invisibly `NULL`. Stops with an informative message if the column is
absent.
