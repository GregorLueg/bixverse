# Helper function to generate HotSpot data

Wrapper class that stores the HotSpot

## Usage

``` r
new_sc_hotspot_res(hotspot_res, used_genes, used_cells)
```

## Arguments

- hotspot_res:

  Named list with the following items:

  - z - Pairwise gene-gene matrix with the Z-scores.

  - cor - Pairwise gene-gene matrix with the local correlation
    coefficients.

- used_genes:

  Character vector. The used genes to generate the matrices.

- used_cells:

  Character vector. The used cells to generate the matrices.

## Value

Generates the `sc_hotspot` class.
