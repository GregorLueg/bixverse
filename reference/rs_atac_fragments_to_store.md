# Count fragments into features and write the store

**\[experimental\]**

## Usage

``` r
rs_atac_fragments_to_store(
  fragments,
  f_path_cells,
  f_path_genes,
  barcodes,
  chrom,
  start,
  end,
  mode,
  kind,
  target_size,
  max_mem_gb,
  verbose
)
```

## Arguments

- fragments:

  String. The fragment file.

- f_path_cells:

  String. Cell-based file to write.

- f_path_genes:

  String. Gene-based file to write.

- barcodes:

  Character vector. Cells in store order.

- chrom, start, end:

  The features.

- mode:

  String. `"insertion"`, `"paired-insertion"` or `"fragment"`.

- kind:

  String. `"peak"`, `"tile"` or `"gene_activity"`.

- target_size:

  Numeric. Target size of the initial log-CPM norm.

- max_mem_gb:

  Optional numeric. Budget of the gene file conversion.

- verbose:

  Integer. `0L` - quiet; `1L` - normal; `2L` - detailed.

## Value

A list with `lib_size` and `nnz` per cell, plus `cells_missing`,
`off_target_fraction` and `off_target_chroms` of the pass.
