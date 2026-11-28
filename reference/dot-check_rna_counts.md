# Refuse or warn on a method that assumes RNA counts

Peak, tile, SNP and protein stores are refused. A gene activity store is
refused for methods built on an mRNA count model (`refuse = TRUE`) and
warned about for the rest, which gene length and fragments shared by
neighbouring genes bias.

## Usage

``` r
.check_rna_counts(object, method, refuse)
```

## Arguments

- object:

  The single cell object.

- method:

  String. The method name, for the message.

- refuse:

  Boolean. Refuse a gene activity store instead of warning.

## Value

Invisible `NULL`.
