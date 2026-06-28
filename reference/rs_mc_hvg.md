# Meta cells highly variable genes

**\[experimental\]** Calculates highly variable genes for MetaCells or
more generally speaking sparse data. This is happening in-memory
compared to the (usually much) larger single cell data sets.

## Usage

``` r
rs_mc_hvg(sparse_data, hvg_method, loess_span, binning, n_bins, clip_max)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `format`.

- hvg_method:

  String. Which HVG detection method to use. Options are
  `c("vst", "meanvarbin", "dispersion")`.

- loess_span:

  Numeric. The span parameter for the loess function (only used for
  `"vst"`).

- binning:

  String. The binning strategy for the `meanvarbin` and `dispersion`
  methods. One of `c("equal_width", "equal_frequency")`.

- n_bins:

  Integer. Number of bins for the `meanvarbin` and `dispersion` methods.

- clip_max:

  Optional clipping number. Defaults to `sqrt(no_cells)` if not provided
  (only used for `"vst"`).

## Value

A list with the HVG statistics. If `hvg_method == "vst"`:

- mean - The average expression of the gene.

- var - The variance of the gene.

- var_exp - The expected variance of the gene.

- var_std - The standardised variance of the gene.

For `"meanvarbin"` and `"dispersion"`:

- mean - The average expression of the gene.

- dispersion - The dispersion of the gene.

- dispersion_scaled - The scaled dispersion per bin per gene.

- bin - The bin of the gene.
