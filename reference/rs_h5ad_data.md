# Load in h5ad data via Rust

**\[experimental\]** Loads in h5ad data within Rust and automatically
converts the data into CSR with cells x genes.

## Usage

``` r
rs_h5ad_data(f_path, cs_type, nrows, ncols, cell_quality, slot, verbose)
```

## Arguments

- f_path:

  File path. The path to the h5ad file.

- cs_type:

  String. One of `c("csr", "csc")`. How the data is stored in the file.
  Other values raise an error.

- nrows:

  Integer. Number of rows in the file.

- ncols:

  Integer. Number of columns in the file.

- cell_quality:

  List. Specifying the cell quality. Please refer to
  [`params_sc_min_quality()`](https://gregorlueg.github.io/bixverse/reference/params_sc_min_quality.md).

- slot:

  String. In which slot the raw data can be found. One of
  `c("X", "raw", "layers.counts")`. Unknown strings default to `"X"`.

- verbose:

  Boolean. Controls verbosity of the function

## Value

A list with the CSR data (cells x genes) of the cells and genes passing
`cell_quality`:

- data - The counts of the sparse matrix.

- indices - The 0-based gene indices of the sparse matrix.

- indptr - The index pointers of the sparse matrix.

- no_genes - No of genes in the sparse matrix (i.e., ncol).

- no_cells - No of cells in the sparse matrix (i.e., nrow).
