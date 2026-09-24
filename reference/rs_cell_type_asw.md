# Calculate cell type silhouette width from an embedding

**\[experimental\]** Average silhouette width on cell type labels,
rescaled to `[0, 1]` via `(s + 1) / 2`. Higher values mean cell types
stay separated.

## Usage

``` r
rs_cell_type_asw(embedding, labels, max_cells, verbose, seed)
```

## Arguments

- embedding:

  Numeric matrix. Cells x dimensions.

- labels:

  Integer vector. The cell type per cell. The codes need not be 0-based
  or contiguous.

- max_cells:

  Integer or NULL. If not NULL, subsample to this many cells for
  performance. If NULL, all cells are used.

- verbose:

  Boolean. Controls verbosity of the function.

- seed:

  Integer. Seed for subsampling reproducibility.

## Value

A list with the following items

- per_cell - Per-cell rescaled silhouette scores

- mean_asw - Mean rescaled silhouette width

- median_asw - Median rescaled silhouette width

## References

Luecken, et al., Nat Methods, 2022
