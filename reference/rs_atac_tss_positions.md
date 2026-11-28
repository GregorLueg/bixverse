# TSS positions from a GTF

**\[experimental\]**

## Usage

``` r
rs_atac_tss_positions(path, biotypes)
```

## Arguments

- path:

  String. The GTF, plain or gzip.

- biotypes:

  Character vector. Biotypes to keep, empty for all.

## Value

A list with `chrom`, `pos` (0-based), `strand` (`-1` minus, `1` plus),
`gene_id` and `gene_name`, one entry per kept gene.
