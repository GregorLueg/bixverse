# Identify HVGs (batch aware)

This is a helper function to identify highly variable genes in a
batch-aware manner. At the moment the implementation has only the
VST-based version (known as Seurat v3). The other methods will be
implemented in the future. This function will calculate the HVG per
given experimental batch and you can choose the way to combine them. The
choices are union (of Top x HVG per batch), based on the average
variance per batch or only take genes that are amongst the Top X HVG in
all batches. Important. The function returns 0-indices for the genes!

## Usage

``` r
find_hvg_batch_aware_sc(
  object,
  batch_column,
  hvg_no = 2000L,
  gene_comb_method = c("union", "average", "intersection"),
  hvg_params = params_sc_hvg(),
  streaming = NULL,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- batch_column:

  String. The column name of the batch column in the obs table.

- hvg_no:

  Integer. Number of highly variable genes to include. Defaults to
  `2000L`.

- gene_comb_method:

  String. One of `c("union", "average", "intersection")`. The method to
  combine the HVG across the different batches. Defaults to `"union"`.

- hvg_params:

  List, see
  [`params_sc_hvg()`](https://gregorlueg.github.io/bixverse/reference/params_sc_hvg.md).
  This list contains

  - method - Which method to use. One of
    `c("vst", "meanvarbin", "dispersion")`

  - loess_span - The span for the loess function to standardise the
    variance

  - num_bin - Integer. Not yet implemented.

  - bin_method - String. One of `c("equal_width", "equal_freq")`. Not
    implemented yet.

- streaming:

  Optional Boolean. Shall the data be streamed in. Useful for larger
  data sets where you wish to avoid loading in the whole data. If
  `NULL`, will automatically detect.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

This function will return a list with:

- hvg_genes - The gene names of the HVGs.

- hvg_gene_idx - The (0-index) gene features.

- batch_hvg_data - data.table with the detailed information of the
  variance per batch.

## Examples

``` r
# highly variable genes taken as the union over the batches
sc <- demo_single_cells(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 600L, n_genes = 50L, n_batches = 3L
  )
)
hvg <- find_hvg_batch_aware_sc(
  sc,
  hvg_no = 20L,
  batch_column = "batch_index",
  gene_comb_method = "union",
  .verbose = FALSE
)
head(hvg$hvg_genes)
#>        38        37        43        20        36        35 
#> "gene_39" "gene_38" "gene_44" "gene_21" "gene_37" "gene_36" 

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
