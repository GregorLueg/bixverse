# Read obs and var tables and metadata from an h5ad file

Useful for exploring the data stored in an h5ad file.

## Usage

``` r
read_h5ad_metadata(f_path)
```

## Arguments

- f_path:

  File path to the `.h5ad` file.

## Value

A list with:

- obs - data.table of cell-level metadata

- var - data.table of gene-level metadata

- dims - named integer vector c(obs, var)

- type - "CSR" or "CSC"

## Examples

``` r
# obs and var tables straight out of the file
data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
)
f_path <- tempfile(fileext = ".h5ad")
write_h5ad_sc(f_path, data$counts, data$obs, data$var, .verbose = FALSE)
meta <- read_h5ad_metadata(f_path)
head(meta$var, 3)
#>        .id ensembl_id
#>     <char>     <char>
#> 1: gene_01     ens_01
#> 2: gene_02     ens_02
#> 3: gene_03     ens_03

unlink(f_path)
```
