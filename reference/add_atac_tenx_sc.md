# Add the peaks of a 10x multiome h5 as the ATAC store

Reads the `Peaks` rows of a Cell Ranger ARC
`filtered_feature_bc_matrix.h5` into a peak store under
`<dir_data>/atac`, for every cell of the RNA store and in its order, so
the RNA cell filter applies unchanged. The features go to the DuckDB
table `var_atac`, the per-cell totals to obs as `atac_lib_size` and
`atac_nnz`. The 10x peak matrix counts insertions.

## Usage

``` r
add_atac_tenx_sc(
  object,
  h5_path,
  overwrite = FALSE,
  max_mem_gb = NULL,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCellsMultiModal` class, its RNA loaded from the same h5 (or
  with the same barcodes).

- h5_path:

  String. The Cell Ranger ARC h5.

- overwrite:

  Boolean. Replace an existing ATAC store.

- max_mem_gb:

  Optional numeric. Memory budget of the peak-major conversion.

- .verbose:

  Boolean. Controls verbosity.

## Value

The object with the ATAC store attached.
