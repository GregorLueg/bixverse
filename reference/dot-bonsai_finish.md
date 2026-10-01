# Build the BonsaiTree and time the whole call

Build the BonsaiTree and time the whole call

## Usage

``` r
.bonsai_finish(
  rs_res,
  started,
  cell_idx,
  cell_names,
  leaf_sizes,
  genes_in,
  gene_ids,
  bonsai_params
)
```

## Arguments

- rs_res:

  List. The raw return of the Rust entry point.

- started:

  POSIXct. When the Rust call began.

- cell_idx:

  Integer. The cells or metacells of the leaves (0-indexed).

- cell_names:

  Character vector. Their names.

- leaf_sizes:

  Integer. Cells behind each leaf.

- genes_in:

  Integer. The candidate genes (0-indexed).

- gene_ids:

  Character vector. Identifiers of `genes_in`, same order.

- bonsai_params:

  List. See
  [`params_sc_bonsai()`](https://gregorlueg.github.io/bixverse/reference/params_sc_bonsai.md).

## Value

A `BonsaiTree`, with the `total` row appended to its timings.
