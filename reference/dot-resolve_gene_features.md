# Resolve feature names to gene indices, keeping the two in step

[`get_gene_indices()`](https://gregorlueg.github.io/bixverse/reference/get_gene_indices.md)
warns and silently drops anything it cannot match, so the caller is left
holding a `features` vector longer than the indices it got back.
Anything that labels a Rust result by position has to use the surviving
names rather than the requested ones, otherwise the labels slide off the
values. This is the in-memory
[`.match_features()`](https://gregorlueg.github.io/bixverse/reference/dot-match_features.md)
for the streaming paths.

## Usage

``` r
.resolve_gene_features(object, features)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- features:

  Character vector. The requested gene ids.

## Value

A list with `features`, the surviving gene ids in request order, and
`indices`, their 0-indexed positions for Rust.
