# Gene body plus promoter regions from a GTF

**\[experimental\]** The regions of Signac's `GeneActivity()`.

## Usage

``` r
rs_atac_gene_activity_regions(path, upstream, downstream, max_width, biotypes)
```

## Arguments

- path:

  String. The GTF.

- upstream:

  Integer. Bases added upstream.

- downstream:

  Integer. Bases added downstream.

- max_width:

  Integer. Genes this long or longer are dropped.

- biotypes:

  Character vector. Biotypes to keep, empty for all.

## Value

A list with `chrom`, `start`, `end`, `name` (unique) and `gene_id` (as
in the GTF).
