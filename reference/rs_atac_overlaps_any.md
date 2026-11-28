# Flag regions that overlap any of another set

**\[experimental\]**

## Usage

``` r
rs_atac_overlaps_any(chrom, start, end, other_chrom, other_start, other_end)
```

## Arguments

- chrom, start, end:

  The regions to test.

- other_chrom, other_start, other_end:

  The set to test against.

## Value

Logical vector, one flag per region.
