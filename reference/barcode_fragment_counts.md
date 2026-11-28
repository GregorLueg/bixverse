# Fragments per barcode, every barcode in the file

For cell calling: plot the counts against their rank (a knee plot) and
pick the `min_fragments` of
[`load_atac_fragments()`](https://gregorlueg.github.io/bixverse/reference/load_atac_fragments.md).

## Usage

``` r
barcode_fragment_counts(
  fragments,
  exclude_chroms = c("chrM", "MT", "chrMT", "M")
)
```

## Arguments

- fragments:

  String. The fragment file.

- exclude_chroms:

  Character vector. Chromosomes not counted.

## Value

A data.table with `barcode` and `n_fragments`, most fragments first.
