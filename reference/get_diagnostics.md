# Get diagnostics from a BulkModuleResult

Returns one diagnostic by key, or the whole named list if `which` is
`NULL`. Contents are method-specific: CoReMo stores `stability` and
`cluster_quality`; Leiden stores `resolution_used` and `modularity`; ICA
stores `ica_meta` and `stability_scores`; DGRDL stores
`feature_laplacian`, `sample_laplacian`, and the grid-search table; NMF
stores `final_loss`, `n_iter`, `converged`, and (for stabilised runs)
`losses`, `best_idx`.

## Usage

``` r
get_diagnostics(object, which = NULL)
```

## Arguments

- object:

  A `BulkModuleResult`.

- which:

  String or `NULL`. Name of the diagnostic to return. If `NULL`, returns
  the whole list.

## Value

The requested diagnostic, the named list, or `NULL` (with warning) if
`which` is not among the stored diagnostic keys.
