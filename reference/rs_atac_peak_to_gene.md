# Peak to gene links over groups of cells

**\[experimental\]** ArchR's `addPeak2GeneLinks()` on pseudobulked
counts.

## Usage

``` r
rs_atac_peak_to_gene(
  peak_counts,
  gene_counts,
  chrom,
  start,
  end,
  tss_chrom,
  tss_pos,
  max_distance,
  scale
)
```

## Arguments

- peak_counts:

  Numeric matrix. Groups x peaks, summed counts.

- gene_counts:

  Numeric matrix. Groups x genes, summed counts.

- chrom, start, end:

  The peaks.

- tss_chrom:

  Character vector. TSS chromosome per gene column.

- tss_pos:

  Integer vector. TSS position per gene column.

- max_distance:

  Integer. Largest peak to TSS distance.

- scale:

  Numeric. Normalisation scale.

## Value

A list with `peak`, `gene` (0-based column indices), `distance`, `r`,
`p` and `fdr`, one entry per tested pair.
