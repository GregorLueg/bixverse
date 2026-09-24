# Write scTransform-corrected counts to a new store

Reverses the residual transform with every latent variable, the library
size included, held at its median. The depth structure goes, the
per-sample intercept stays.

The result is a new store on disk, not a layer on `object`. Its gene
axis is the model's, so it is narrower than the source and the indices
do not line up, which is why the observation and variable tables are
rebuilt rather than copied.

## Usage

``` r
sct_corrected_counts_sc(
  object,
  dir_out = NULL,
  build_cell_store = TRUE,
  overwrite = FALSE,
  gene_batch_size = NULL,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class with a scTransform fit
  attached. The analytic Pearson model has no corrected-count
  equivalent.

- dir_out:

  String or `NULL`. Directory to write to. `NULL` uses `sct_corrected`
  inside the object's own data directory.

- build_cell_store:

  Boolean. Also write the `counts_cells.bin` companion and the database,
  giving back a `SingleCells` rather than a path. Costs a second pass
  over the data.

- overwrite:

  Boolean. Overwrite an existing store in `dir_out`.

- gene_batch_size:

  Integer or `NULL`. Genes held in memory per batch.

- .verbose:

  Boolean or Integer. Controls verbosity.

## Value

With `build_cell_store = TRUE` a new `SingleCells` over the corrected
counts, otherwise the path of the gene-major file, invisibly.

## Examples

``` r
# corrected counts as a fresh object
sc <- demo_single_cells(prepped = FALSE)
sc <- fit_residuals_sc(sc, .verbose = FALSE)
corrected <- sct_corrected_counts_sc(sc, .verbose = FALSE)
dim(corrected)
#> [1] 500  50

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
