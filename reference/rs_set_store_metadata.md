# Set caller metadata on both files of a v5 store

**\[experimental\]** Given keys replace, all others are kept.

## Usage

``` r
rs_set_store_metadata(f_path_cells, f_path_genes, keys, values)
```

## Arguments

- f_path_cells:

  String. Cell-based binary.

- f_path_genes:

  String. Gene-based binary.

- keys:

  Character vector. The keys.

- values:

  Character vector. The values, same length as `keys`.

## Value

Invisible `NULL`. Errors for a multi-line key or value, or a v4 store.
