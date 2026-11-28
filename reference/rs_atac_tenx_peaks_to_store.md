# The `Peaks` modality of a 10x multiome h5 into a peak store

**\[experimental\]**

## Usage

``` r
rs_atac_tenx_peaks_to_store(
  h5_path,
  f_path_cells,
  f_path_genes,
  barcodes,
  max_mem_gb,
  verbose
)
```

## Arguments

- h5_path:

  String. The Cell Ranger ARC h5.

- f_path_cells:

  String. Cell-based file to write.

- f_path_genes:

  String. Gene-based file to write.

- barcodes:

  Character vector. Cells in store order.

- max_mem_gb:

  Optional numeric. Budget of the gene file conversion.

- verbose:

  Integer. `0L` - quiet; `1L` - normal; `2L` - detailed.

## Value

A list with `name`, `chrom`, `start`, `end` per peak and `lib_size`,
`nnz` per cell.
