# Default parameters for generation of synthetic CellSweep data

Shapes the fixture with a planted ambient profile that
[`generate_cellsweep_test_data()`](https://gregorlueg.github.io/bixverse/reference/generate_cellsweep_test_data.md)
builds. The defaults give 600 real barcodes over 3 cell types plus 2000
empty droplets, on 200 genes.

## Usage

``` r
params_sc_synthetic_cellsweep(
  n_real = 600L,
  n_empty = 2000L,
  n_genes = 200L,
  n_celltypes = 3L,
  n_markers = 20L,
  marker_weight = 25,
  ambient_dominance = 0.6,
  alpha_mean = 0.3,
  alpha_sd = 0.12,
  real_lib_size = 3000L,
  empty_lib_size = 120L
)
```

## Arguments

- n_real:

  Integer. Number of real barcodes. Cell types are assigned round-robin
  over them. Defaults to `600L`.

- n_empty:

  Integer. Number of empty droplets. The ambient profile is estimated
  off these, so at least 30 and preferably a lot more. Defaults to
  `2000L`.

- n_genes:

  Integer. Number of genes. Defaults to `200L`.

- n_celltypes:

  Integer. Number of cell types. Defaults to `3L`.

- n_markers:

  Integer. Width of each cell type's marker block. The blocks are
  contiguous and disjoint, so `n_markers * n_celltypes` has to fit into
  `n_genes`. Defaults to `20L`.

- marker_weight:

  Numeric. Enrichment of a marker gene over background in its own cell
  type's profile. Must exceed 1. Defaults to `25.0`.

- ambient_dominance:

  Numeric. Fraction of the soup coming from the first cell type. The
  remainder is flat background. Defaults to `0.6`.

- alpha_mean:

  Numeric. Mean planted ambient fraction across real barcodes. Defaults
  to `0.3`.

- alpha_sd:

  Numeric. Spread of the planted ambient fraction. Defaults to `0.12`.

- real_lib_size:

  Integer. Expected library size of a real barcode. Defaults to `3000L`.

- empty_lib_size:

  Integer. Expected library size of an empty droplet. Defaults to
  `120L`.

## Value

A named list with the following elements:

- n_real - Integer. Number of real barcodes. Cell types are assigned
  round-robin over them. Defaults to `600L`.

- n_empty - Integer. Number of empty droplets. The ambient profile is
  estimated off these, so at least 30 and preferably a lot more.
  Defaults to `2000L`.

- n_genes - Integer. Number of genes. Defaults to `200L`.

- n_celltypes - Integer. Number of cell types. Defaults to `3L`.

- n_markers - Integer. Width of each cell type's marker block. The
  blocks are contiguous and disjoint, so `n_markers * n_celltypes` has
  to fit into `n_genes`. Defaults to `20L`.

- marker_weight - Numeric. Enrichment of a marker gene over background
  in its own cell type's profile. Must exceed 1. Defaults to `25.0`.

- ambient_dominance - Numeric. Fraction of the soup coming from the
  first cell type. The remainder is flat background. Defaults to `0.6`.

- alpha_mean - Numeric. Mean planted ambient fraction across real
  barcodes. Defaults to `0.3`.

- alpha_sd - Numeric. Spread of the planted ambient fraction. Defaults
  to `0.12`.

- real_lib_size - Integer. Expected library size of a real barcode.
  Defaults to `3000L`.

- empty_lib_size - Integer. Expected library size of an empty droplet.
  Defaults to `120L`.

## Details

The soup is the first cell type plus flat background rather than a
mixture of every cell type profile. A soup sitting in the span of the
cell type profiles makes the contamination fraction unidentifiable, and
the fixture would then be testing the repulsion term rather than the
model.
