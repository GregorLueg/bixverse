# Merge fragment files into one indexed file

**\[experimental\]** Barcodes get a prefix and suffix per input, so
samples stay apart.

## Usage

``` r
rs_atac_merge_fragments(paths, prefixes, suffixes, out)
```

## Arguments

- paths:

  Character vector. The inputs.

- prefixes:

  Character vector. Barcode prefix per input.

- suffixes:

  Character vector. Barcode suffix per input.

- out:

  String. Output path; the index goes to `<out>.tbi`.

## Value

A list with `n_fragments` per input and `chroms` of the output.
