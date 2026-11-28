# Merge per-sample fragment files into one

As `cellranger-arc aggr` does: each file's barcodes get its prefix and
suffix, so the same barcode in two samples stays two cells; give every
file its own. The output is a sorted BGZF file with a tabix index, which
every ATAC function reads as one experiment.

## Usage

``` r
merge_fragments(fragments, out, prefixes = NULL, suffixes = NULL)
```

## Arguments

- fragments:

  Character vector. The fragment files.

- out:

  String. The merged file; the index goes next to it as `.tbi`.

- prefixes:

  Optional character vector. Barcode prefix per file, e.g.
  `c("s1_", "s2_")`.

- suffixes:

  Optional character vector. Barcode suffix per file, e.g.
  `c("-1", "-2")`.

## Value

A data.table with `fragments`, `prefix`, `suffix` and `n_fragments` per
input.
