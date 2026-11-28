# Per-cell QC from a fragment file

**\[experimental\]** Signac's fragment-based metrics: fragments per
cell, nucleosome signal, TSS enrichment and the mitochondrial fraction.

## Usage

``` r
rs_atac_fragment_qc(fragments, barcodes, tss_chrom, tss_pos, mito_chroms)
```

## Arguments

- fragments:

  String. The fragment file.

- barcodes:

  Character vector. Cells, output order.

- tss_chrom:

  Character vector. TSS chromosomes.

- tss_pos:

  Integer vector. TSS positions (0-based).

- mito_chroms:

  Character vector. Names of the mitochondrial genome.

## Value

A list with per-cell `n_fragments`, `nucleosome_signal`,
`tss_enrichment` and `mito_fraction`, plus the coverage keys (off
target: on a chromosome without a TSS).
