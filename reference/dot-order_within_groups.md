# Order indices group by group

Order indices group by group

## Usage

``` r
.order_within_groups(groups, order_fun)
```

## Arguments

- groups:

  Factor. Groups in the original order.

- order_fun:

  Function. Takes the integer indices of one group and returns them
  reordered.

## Value

Integer vector of indices, groups in level order.
