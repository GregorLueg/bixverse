# rLSI integration

Signac's `FindIntegrationAnchors(reduction = "rlsi")` followed by
`IntegrateEmbeddings()`. Per batch: IDF over its own cells, TF-IDF and
LSI of `features`; anchors between every batch pair in the reciprocal
LSI spaces; then `embedding` is corrected along the anchor tree. Written
to the ATAC cache as `"rlsi"`.

Signac integrates all LSI components and leaves component 1 out
afterwards:
`calculate_lsi_sc(x, lsi_params = params_atac_lsi(drop_first = FALSE))`,
`run_rlsi_sc(x, ...)`, then columns 2 onwards of `"rlsi"`. With the
default `drop_first = TRUE` the depth component is left out of the
integration too, which changes the correction weights.

## Usage

``` r
run_rlsi_sc(
  object,
  batch_column,
  embedding = "lsi",
  features = NULL,
  rlsi_params = params_atac_rlsi(),
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCellsAtac` or `SingleCellsMultiModal` class, after
  [`calculate_lsi_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_lsi_sc.md).

- batch_column:

  String. obs column with the batch of each cell.

- embedding:

  String. The ATAC embedding to integrate.

- features:

  Optional integer vector. 1-based ATAC features of the per-batch LSIs.
  `NULL` uses those of
  [`calculate_lsi_sc()`](https://gregorlueg.github.io/bixverse/reference/calculate_lsi_sc.md).

- rlsi_params:

  List. Output of
  [`params_atac_rlsi()`](https://gregorlueg.github.io/bixverse/reference/params_atac_rlsi.md).

- seed:

  Integer. Seed of the solvers and the kNN backend.

- .verbose:

  Boolean. Controls verbosity.

## Value

The object with the `"rlsi"` embedding.

## References

Stuart, et al., Nat Methods, 2021
