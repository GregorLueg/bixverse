# Collect the cached embeddings and graphs of a `SingleCells` object

PCA factors go to `obsm/X_pca` and the loadings to `varm/PCs`, every
other embedding to `obsm/X_<name>`, and the sNN graph to
`obsp/connectivities` together with the `uns/neighbors` block ScanPy
looks for. Anything that is not cached, or that no longer matches the
cells being written, is skipped.

## Usage

``` r
.h5ad_export_cache(object, .verbose = TRUE)
```

## Arguments

- object:

  `SingleCells` class.

- .verbose:

  Boolean. Controls the verbosity of the function.

## Value

A list with `obsm` and `varm` (named lists of double matrices), `obsp`
(named list of CSR matrices as `indptr`, `indices`, `data`) and
`uns_json` (string).
