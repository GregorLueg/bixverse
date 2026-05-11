# Transform the counts to a Rust-specific list

Helper function to transform the counts from the `MetaCells` into a Rust
specific list.

## Usage

``` r
mc_counts_to_list(
  object,
  cell_indices = NULL,
  gene_indices = NULL,
  assay = c("raw", "norm")
)
```

## Arguments

- object:

  `MetaCells` class.

- cell_indices:

  Optional integer. Defines the indices of the (meta)cells to extract.

- gene_indices:

  Optional integer. Defines the indices of the genes to extract.

- assay:

  String. One of `c("raw", "norm")`
