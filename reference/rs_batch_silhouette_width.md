# Calculate batch silhouette width from an embedding

Computes the average silhouette width on batch labels using pairwise
distances in the embedding space. Values near 0 indicate good batch
mixing, values near 1 indicate batch separation.

## Usage

``` r
rs_batch_silhouette_width(embedding, batch_vector, max_cells, seed)
```

## Arguments

- embedding:

  Numeric matrix. The embedding to assess (e.g. PCA or corrected
  embedding). Rows are cells, columns are dimensions.

- batch_vector:

  Integer vector. The integers indicate to which batch a given cell
  belongs.

- max_cells:

  Integer or NULL. If not NULL, subsample to this many cells for
  performance. Defaults to 5000.

- seed:

  Integer. Seed for subsampling reproducibility.

## Value

A list with the following items

- per_cell - Per-cell silhouette scores

- mean_asw - Mean silhouette width

- median_asw - Median silhouette width
