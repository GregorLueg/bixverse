# Single cell test data with a planted ambient profile

This function generates synthetic data for CellSweep test purposes.
Every real barcode is a two-component multinomial: a planted fraction
`alpha` of its library comes from the soup, the rest from its own cell
type profile. Empty droplets are pure soup at a much smaller library
size, which is what the ambient profile is estimated off. Real barcodes
come first, the empty droplets after.

## Usage

``` r
generate_cellsweep_test_data(
  syn_data_params = params_sc_synthetic_cellsweep(),
  seed = 42L
)
```

## Arguments

- syn_data_params:

  List. Contains the parameters for the generation of the synthetic
  data, see:
  [`params_sc_synthetic_cellsweep()`](https://gregorlueg.github.io/bixverse/reference/params_sc_synthetic_cellsweep.md).

- seed:

  Integer. The seed for the generation of the synthetic data.

## Value

List with the following items

- counts - dgRMatrix with cells x genes.

- obs - data.table with `cell_id`, `cell_grp` (`NA` for the empty
  droplets), `sample_id`, `is_empty` and `alpha_true` (`NA` for the
  empty droplets).

- var - data.table that contains the gene information.

- ambient_true - Numeric vector. The soup the empty droplets were drawn
  from, named by gene and summing to one.

- celltype_profiles_true - Numeric matrix of cell types x genes. Each
  row sums to one.

## Details

The empty droplets carry no cell type label, which is exactly what
[`cellsweep_sc()`](https://gregorlueg.github.io/bixverse/reference/cellsweep_sc.md)
keys off: barcodes that are neither empty nor annotated are excluded
from the fit. Load the counts with a fully permissive
[`params_sc_min_quality()`](https://gregorlueg.github.io/bixverse/reference/params_sc_min_quality.md),
otherwise the ingest deletes the empty droplets the model trains on.

## Examples

``` r
# a small synthetic experiment with a planted soup
data <- generate_cellsweep_test_data(
  syn_data_params = params_sc_synthetic_cellsweep(
    n_real = 60L,
    n_empty = 200L,
    n_genes = 60L
  )
)
dim(data$counts)
#> [1] 260  60
head(data$obs, 3)
#>     cell_id    cell_grp sample_id is_empty alpha_true
#>      <char>      <char>    <char>   <lgcl>      <num>
#> 1: cell_001 cell_type_1 sample_01    FALSE 0.37124013
#> 2: cell_002 cell_type_2 sample_01    FALSE 0.09390341
#> 3: cell_003 cell_type_3 sample_01    FALSE 0.23890524
```
