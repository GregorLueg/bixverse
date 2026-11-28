# Merge overlapping and abutting regions

**\[experimental\]**

## Usage

``` r
rs_atac_reduce_regions(chrom, start, end)
```

## Arguments

- chrom:

  Character vector. Chromosomes.

- start:

  Integer vector. 0-based starts.

- end:

  Integer vector. Exclusive ends.

## Value

A list with `chrom`, `start` and `end` of the merged regions and
`members`, a list with the input indices (0-based) each one absorbed.
