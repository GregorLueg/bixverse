# Identify HVGs (batch aware)

This is a helper function to identify highly variable genes in a
batch-aware manner. At the moment the implementation has only the
VST-based version (known as Seurat v3). The other methods will be
implemented in the future. This function will calculate the HVG per
given experimental batch and you can choose the way to combine them. The
choices are union (of Top x HVG per batch), based on the average
variance per batch or only take genes that are amongst the Top X HVG in
all batches.

## Usage

``` r
find_hvg_batch_aware_sc(
  object,
  batch_column,
  hvg_no = 2000L,
  gene_comb_method = c("union", "average", "intersection"),
  hvg_params = params_sc_hvg(),
  streaming = FALSE,
  .verbose = TRUE
)
```

## Arguments

- object:

  `single_cell_exp` class.

- batch_column:

  String. The column name of the batch column in the obs table.

- hvg_no:

  Integer. Number of highly variable genes to include. Defaults to
  `2000L`.

- gene_comb_method:

  String. One of `c("union", "average", "intersection")`. The method to
  combine the HVG across the different batches. Defaults to `"union"`.

- hvg_params:

  List, see [`params_sc_hvg()`](params_sc_hvg.md). This list contains

  - method - Which method to use. One of
    `c("vst", "meanvarbin", "dispersion")`

  - loess_span - The span for the loess function to standardise the
    variance

  - num_bin - Integer. Not yet implemented.

  - bin_method - String. One of `c("equal_width", "equal_freq")`. Not
    implemented yet.

- streaming:

  Boolean. Shall the genes be streamed in. Useful for larger data sets
  where you wish to avoid loading in the whole data. Defaults to
  `FALSE`.

- .verbose:

  Boolean. Controls verbosity and returns run times.

## Value

This function will return a list with:

- hvg_indices - The indices of the batch-aware HVG.

- batch_hvg_data - data.table with the detailed information of the
  variance per batch.
