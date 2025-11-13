# Load in h5ad data via Rust

Loads in h5ad data within Rust and automatically converts the data into
CSR with cells x genes.

## Usage

``` r
rs_h5ad_data(f_path, cs_type, nrows, ncols, cell_quality, verbose)
```

## Arguments

- f_path:

  File path. The path to the h5ad file.

- cs_type:

  String. Is the data stored in CSC or CSR.

- nrows:

  Integer. Number of rows in the file.

- ncols:

  Integer. Number of columns in the file.

- cell_quality:

  List. Specifiying the cell quality. Please refer to
  [`params_sc_min_quality()`](params_sc_min_quality.md).

- verbose:

  Boolean. Controls verbosity of the function

## Value

A list with:

- data - The data of the sparse matrix stored on the h5ad file.

- indices - The indices of the sparse matrix stored in the h5ad file.

- indptr - The indptr of the sparse matrix stored in the h5ad file.

- no_genes - No of genes in the sparse matrix (i.e., ncol).

- no_cells - No of cells in the sparse matrix (i.e., nrow).
