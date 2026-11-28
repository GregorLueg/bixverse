# TSS profiles and fragment length histograms per group of cells

**\[experimental\]** SnapATAC2's `TSS_profile` and `frag_size_distr`.

## Usage

``` r
rs_atac_fragment_profiles(
  fragments,
  barcodes,
  groups,
  tss_chrom,
  tss_pos,
  tss_strand,
  half_width,
  max_len
)
```

## Arguments

- fragments:

  String. The fragment file.

- barcodes:

  Character vector. Cells.

- groups:

  Integer vector. Group of each cell, `0..n_groups`.

- tss_chrom:

  Character vector. TSS chromosomes.

- tss_pos:

  Integer vector. TSS positions (0-based).

- tss_strand:

  Integer vector. TSS strands, `-1` minus, else plus.

- half_width:

  Integer. Profile half width.

- max_len:

  Integer. Longest fragment with its own bin.

## Value

A list with `tss` (groups x `2 * half_width + 1`), `lengths` (groups x
`max_len + 1`, bin 0 the longer fragments) and the coverage keys.
