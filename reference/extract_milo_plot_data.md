# Extract the miloR neighbourhood graph positioned on an embedding

Places every tested neighbourhood of a
[`get_miloR_abundances_sc()`](https://gregorlueg.github.io/bixverse/reference/get_miloR_abundances_sc.md)
result at the embedding coordinate of its index cell, as Milo's
`plotNhoodGraphDA()` does, and returns the neighbourhood graph as node
and edge tables. Two neighbourhoods are connected by the number of cells
they share.

The logFC is returned as is. Masking non-significant neighbourhoods is
left to the plot, `is_sig` carries the call.

## Usage

``` r
extract_milo_plot_data(
  object,
  milo_res,
  embedding = "umap",
  alpha = 0.1,
  overlap = 1L,
  ...
)
```

## Arguments

- object:

  A single cell class. The one the miloR result was built on.

- milo_res:

  `miloR` class. Needs to have been through
  [`test_nhoods()`](https://gregorlueg.github.io/bixverse/reference/test_nhoods.md).

- embedding:

  String. Name of the embedding to position the nodes in.

- alpha:

  Numeric. Spatial FDR threshold for `is_sig`. Defaults to `0.1`.

- overlap:

  Integer. Minimum number of shared cells for an edge. Defaults to `1L`.

- ...:

  Additional arguments forwarded to
  [`extract_embedding_data()`](https://gregorlueg.github.io/bixverse/reference/extract_embedding_data.md)
  and onward to
  [`get_embedding()`](https://gregorlueg.github.io/bixverse/reference/get_embedding.md)
  (e.g. `modality`).

## Value

A list with the embedding stored as an `embedding` attribute and

- nodes - data.table with `Nhood`, `dim_1`, `dim_2`, `size` (cells in
  the neighbourhood), every column of the differential abundance results
  and `is_sig`.

- edges - data.table with `from`, `to`, `weight` (shared cells) and the
  segment coordinates `x`, `y`, `xend`, `yend`.

## References

Dann, et al., Nat Biotechnol, 2022
