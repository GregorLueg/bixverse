# Calculate the PCA on top of the normalised ADT counts

This function will run PCA - via (randomised) SVD - on the normalised
counts and add the PCA results to the ScCache for the ADT counts.

## Usage

``` r
calculate_pca_adt_sc(
  object,
  no_pcs,
  features = NULL,
  randomised_svd = FALSE,
  seed = 42L
)
```

## Arguments

- object:

  `SingleCellsMultiModal` class with ADT counts added.

- no_pcs:

  Integer. Number of PCs to calculate.

- features:

  Optional string vector. If you want to subset to a specific set of ADT
  probes (for example to exclude isotypes).

- randomised_svd:

  Boolean. Shall randomised SVD be used. Faster, but less precise.

- seed:

  Integer. Controls reproducibility. Only relevant if
  `randomised_svd = TRUE`.

## Value

The function will add the PCA factors, loadings and singular values for
the ADT data to the object.

## Examples

``` r
# PCA over the CLR-normalised protein counts
rna <- generate_single_cell_test_data()
adt <- generate_single_cell_test_data_adt()
dir <- tempfile("bixverse_mm")
dir.create(dir)
object <- load_r_data(
  SingleCellsMultiModal(dir_data = dir),
  counts = rna$counts,
  obs = rna$obs,
  var = rna$var,
  sc_qc_param = params_sc_min_quality(min_unique_genes = 5L),
  .verbose = FALSE
)
object <- add_adt_counts_sc(object, adt_counts = adt$counts, method = "clr")
object <- calculate_pca_adt_sc(object, no_pcs = 10L)
dim(get_pca_factors(object, modality = "adt"))
#> [1] 1000   10

unlink(dir, recursive = TRUE, force = TRUE)
```
