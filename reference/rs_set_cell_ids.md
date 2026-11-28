# Record the cell IDs of a v5 store in both headers

**\[experimental\]** Only their hash is stored. Stores that share a cell
axis (a gene activity store and its source) carry the same one.

## Usage

``` r
rs_set_cell_ids(f_path_cells, f_path_genes, ids)
```

## Arguments

- f_path_cells:

  String. Cell-based binary.

- f_path_genes:

  String. Gene-based binary.

- ids:

  Character vector. Cell IDs in store order.

## Value

Invisible `NULL`. Errors for a store that is not v5 or a wrong number of
IDs.
