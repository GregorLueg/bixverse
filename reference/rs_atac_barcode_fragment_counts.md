# Fragments per barcode

**\[experimental\]** Counts every barcode in the file.

## Usage

``` r
rs_atac_barcode_fragment_counts(fragments, exclude_chroms)
```

## Arguments

- fragments:

  String. The fragment file.

- exclude_chroms:

  Character vector. Chromosomes not counted.

## Value

A list with `barcode` (sorted) and `n_fragments`.
