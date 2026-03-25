# Calculate batch average silhouette width

Computes the average silhouette width (ASW) on batch labels using
pairwise Euclidean distances in the embedding space. For each cell, the
function estimates the mean distance to cells of the same batch (a) and
the mean distance to cells of the nearest other batch (b), then computes
s = (b - a) / max(a, b).

Values near 0 indicate good batch mixing, values near 1 indicate batch
separation, and negative values suggest overcorrection. This metric is
best suited for embedding-based correction methods (e.g. Harmony,
fastMNN). For graph-based methods like BBKNN, consider using
[`calculate_batch_lisi_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_batch_lisi_sc.md)
instead.

## Usage

``` r
calculate_batch_asw_sc(
  object,
  batch_column,
  embd_to_use = "pca",
  max_cells = 5000L,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- batch_column:

  String. The column with the batch information in the obs data of the
  class.

- embd_to_use:

  String. Which embedding to compute the ASW on. One of
  `c("pca", "harmony", "mnn")`. Defaults to `"pca"`.

- max_cells:

  Integer or `NULL`. If not `NULL`, subsample to this many cells for
  performance. The pairwise distance computation is O(n^2), so
  subsampling is recommended for large datasets. Defaults to `5000L`.

- seed:

  Integer. Seed for subsampling reproducibility.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

A `BatchSilhouetteScores` object with the following elements

- per_cell - Per-cell silhouette scores in `[-1, 1]`.

- mean_asw - Mean silhouette width across all cells.

- median_asw - Median silhouette width across all cells.

- n_batches - Number of batches in the data.

- embedding_used - Which embedding the ASW was computed on.
