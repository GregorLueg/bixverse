# Align a group mapping to ids as a factor

Align a group mapping to ids as a factor

## Usage

``` r
.as_group_factor(groups, ids)
```

## Arguments

- groups:

  Optional named character vector or factor. `NULL` puts everything into
  a single group `"all"`.

- ids:

  Character vector. Ids to align to.

## Value

Factor of the same length as `ids`, without unused levels.
