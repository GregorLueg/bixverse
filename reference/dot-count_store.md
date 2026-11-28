# Count a fragment file into a store and report on the pass

Count a fragment file into a store and report on the pass

## Usage

``` r
.count_store(
  fragments,
  store_dir,
  barcodes,
  regions,
  mode,
  kind,
  target_size,
  max_mem_gb,
  .verbose
)
```

## Arguments

- fragments:

  String. The fragment file.

- store_dir:

  String. Directory of the store to write.

- barcodes:

  Character vector. Cells in store order.

- regions:

  List with `chrom`, `start` and `end` of the features.

- mode:

  String. The count mode.

- kind:

  String. One of `c("peak", "tile", "gene_activity")`.

- target_size:

  Numeric. Target size of the initial log-CPM norm.

- max_mem_gb:

  Optional numeric. Memory budget of the feature-major conversion.

- .verbose:

  Boolean. Controls verbosity.

## Value

The Rust result: `lib_size`, `nnz` and the coverage keys.
