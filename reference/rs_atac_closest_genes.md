# The closest gene body to every region

**\[experimental\]**

## Usage

``` r
rs_atac_closest_genes(gtf, biotypes, chrom, start, end)
```

## Arguments

- gtf:

  String. The GTF.

- biotypes:

  Character vector. Biotypes to keep, empty for all.

- chrom, start, end:

  The regions.

## Value

A list with `gene_name` and `gene_id` (`NA` on a chromosome without
genes), `distance` (0 when overlapping) and `tss_distance` (region
centre to TSS, positive downstream).
