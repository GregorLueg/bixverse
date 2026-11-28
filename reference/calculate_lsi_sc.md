# LSI over the TF-IDF of the kept cells

Signac's `RunSVD()`. The embedding (components z-scored over the cells)
goes to the ATAC cache as `"lsi"`, for
[`find_neighbours_sc()`](https://gregorlueg.github.io/bixverse/reference/find_neighbours_sc.md),
[`harmony_sc()`](https://gregorlueg.github.io/bixverse/reference/harmony_sc.md),
[`umap_sc()`](https://gregorlueg.github.io/bixverse/reference/umap_sc.md)
and
[`generate_wnn_graph_sc()`](https://gregorlueg.github.io/bixverse/reference/generate_wnn_graph_sc.md)
with `modality = "atac"`, `embd_to_use = "lsi"`. The loadings, singular
values and features are kept for
[`project_lsi_sc()`](https://gregorlueg.github.io/bixverse/reference/project_lsi_sc.md).

## Usage

``` r
calculate_lsi_sc(
  object,
  features = NULL,
  lsi_params = params_atac_lsi(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCellsAtac` or `SingleCellsMultiModal` class, after
  [`run_tfidf_sc()`](https://gregorlueg.github.io/bixverse/reference/run_tfidf_sc.md).

- features:

  Optional integer vector. 1-based features. `NULL` takes `top_feature`
  from the feature table if
  [`find_top_features_sc()`](https://gregorlueg.github.io/bixverse/reference/find_top_features_sc.md)
  ran, else all.

- lsi_params:

  List. Output of
  [`params_atac_lsi()`](https://gregorlueg.github.io/bixverse/reference/params_atac_lsi.md).

- seed:

  Integer. Seed of the solver.

- .verbose:

  Boolean. Controls verbosity.

## Value

The object with the `"lsi"` embedding and the LSI fit.

## References

Stuart, et al., Nat Methods, 2021
