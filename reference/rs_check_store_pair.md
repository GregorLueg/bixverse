# Error unless a cell file and a gene file belong to the same store

**\[experimental\]** Reads the headers only.

## Usage

``` r
rs_check_store_pair(f_path_cells, f_path_genes)
```

## Arguments

- f_path_cells:

  String. Cell-based binary.

- f_path_genes:

  String. Gene-based binary.

## Value

Invisible `NULL`, or an error naming the field that differs.
