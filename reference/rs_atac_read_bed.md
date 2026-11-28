# Read a BED file

**\[experimental\]** Reads the first three columns of a BED file, plain
or gzip.

## Usage

``` r
rs_atac_read_bed(path)
```

## Arguments

- path:

  String. The BED file.

## Value

A list with `chrom`, `start` (0-based) and `end` (exclusive).
