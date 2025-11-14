# Identify hotspot gene clusters

Identify hotspot gene clusters

## Usage

``` r
set_hotspot_membership(x, fdr_threshold = 0.05, min_size = 10L)

# S3 method for class 'sc_hotspot'
set_hotspot_membership(x, fdr_threshold = 0.05, min_size = 10L)
```

## Arguments

- x:

  An object to generate the hotspot gene clusters for.

- fdr_threshold:

  Numeric. The maximum FDR for a given gene-gene local correlation to be
  included.

- min_size:

  Integer. Minimum cluster size.
