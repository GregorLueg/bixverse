# `chrom`, `start` and `end` of a BED file or a table of regions

`chrom`, `start` and `end` of a BED file or a table of regions

## Usage

``` r
.regions_from(regions)
```

## Arguments

- regions:

  String (a BED file, plain or gzip) or a data.frame with `chrom`,
  `start` and `end` (0-based, half-open).

## Value

A list with `chrom`, `start` and `end`.
