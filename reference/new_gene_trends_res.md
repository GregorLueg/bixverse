# Helper function to generate the gene trend results

Takes the raw Rust output of
[`rs_gene_trends()`](https://gregorlueg.github.io/bixverse/reference/rs_gene_trends.md)
and melts the per-branch grid matrices into one long table, mapping the
0-indexed branch cells back onto cell names.

## Usage

``` r
new_gene_trends_res(rs_res, features, used_cells, branches, params)
```

## Arguments

- rs_res:

  List. The raw return of
  [`rs_gene_trends()`](https://gregorlueg.github.io/bixverse/reference/rs_gene_trends.md).

- features:

  Character vector. The genes that were fitted, in column order of the
  expression matrix.

- used_cells:

  Character vector. The cells the fit ran over, in row order of the
  expression matrix.

- branches:

  Character vector. Branch names, in column order of the fate
  probability matrix.

- params:

  List. The parameters the run used.

## Value

Generates the `GeneTrendsRes` class.
