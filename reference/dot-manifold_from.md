# Parents of a manifold embedding (UMAP, t-SNE, PHATE)

These three read their source embedding from `cache_modality` but write
the result under `modality`, which differ for `"wnn"`. Both parents are
therefore spelled out modality qualified rather than left to the write
modality.

## Usage

``` r
.manifold_from(embd_to_use, cache_modality, modality, has_knn)
```

## Arguments

- embd_to_use:

  String. Name of the source embedding.

- cache_modality:

  String. Modality the source embedding was read from.

- modality:

  String. Modality the result is written to.

- has_knn:

  Boolean. Whether a cached kNN fed the manifold.

## Value

Character vector of parent artefact names.
