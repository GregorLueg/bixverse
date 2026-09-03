# Wrap the NEBULA fits into a result class

Names the coefficient matrices, joins the per-gene test onto the gene
identifiers and hands back a
[`new_sc_nebula_res()`](https://gregorlueg.github.io/bixverse/reference/new_sc_nebula_res.md).

## Usage

``` r
.nebula_res_to_class(res, gene_names, design_mat, params)
```

## Arguments

- res:

  List. The raw return of
  [`rs_nebula_sc()`](https://gregorlueg.github.io/bixverse/reference/rs_nebula_sc.md)
  /
  [`rs_nebula_mc()`](https://gregorlueg.github.io/bixverse/reference/rs_nebula_mc.md).

- gene_names:

  Character vector. All gene identifiers of the object, in store order,
  so `res$gene_idx` indexes into it.

- design_mat:

  Numeric matrix. The design that was fitted, for its column names.

- params:

  Named list. The parameters of the run.

## Value

A `ScNebula` class.
