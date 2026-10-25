# Prepare a source/target network for the ULM / MLM activity scores

Mirrors decoupler's network pruning: targets absent from the features
are dropped, then sources with fewer than `tmin` remaining targets.
Unlike decoupler, duplicated `(source, target)` pairs are not an error;
the first one is kept with a warning.

## Usage

``` r
prep_lm_network(net, features, tmin)
```

## Arguments

- net:

  data.frame. Must have the columns `source` and `target` and optionally
  a numeric `weight` column. Absent weights are set to 1.

- features:

  Character vector. The feature identifiers; the position of a feature
  is its (1-indexed) index in the data.

- tmin:

  Integer. Minimum number of targets per source after the intersection
  with `features`.

## Value

A named list with one element per surviving source, each a list with
`idx` (0-indexed feature positions) and `weight`. `NULL` with a warning
if no source survives.
