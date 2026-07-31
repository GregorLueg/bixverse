# Row-bind the observation tables of meta cell objects

Row-bind the observation tables of meta cell objects

## Usage

``` r
.merge_mc_obs(inputs, source_ids, prefix_ids)
```

## Arguments

- inputs:

  List of `MetaCells` objects.

- source_ids:

  Character vector with the source identifiers.

- prefix_ids:

  Boolean. Prefix the meta cell identifiers.

## Value

A data.table with the merged observations, re-indexed.
