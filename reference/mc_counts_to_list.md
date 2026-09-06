# Transform the counts to a Rust-specific list

Helper function to transform the counts from the `MetaCells` into a Rust
specific list.

## Usage

``` r
mc_counts_to_list(
  object,
  cell_indices = NULL,
  gene_indices = NULL,
  assay = c("raw", "norm")
)
```

## Arguments

- object:

  `MetaCells` class.

- cell_indices:

  Optional integer. Defines the indices of the (meta)cells to extract.

- gene_indices:

  Optional integer. Defines the indices of the genes to extract.

- assay:

  String. One of `c("raw", "norm")`

## Examples

``` r
# the CSR representation Rust expects
sc <- demo_single_cells()
mc <- generate_bt_meta_cells_sc(
  sc,
  sc_meta_cell_params = params_sc_bt_metacells(target_no_metacells = 50L),
  .verbose = FALSE
)
names(mc_counts_to_list(mc, assay = "raw"))
#> [1] "indptr"  "indices" "data"    "cs_type" "nrow"    "ncol"   

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
