# Helper function to generate the PAGA results

Takes the raw Rust output of
[`rs_paga()`](https://gregorlueg.github.io/bixverse/reference/rs_paga.md)
and transforms the two abstracted graphs into sparse matrices with the
cluster levels as dimnames.

## Usage

``` r
new_paga_res(rs_res, cluster_levels, cluster_col, modality)
```

## Arguments

- rs_res:

  List. The raw return of
  [`rs_paga()`](https://gregorlueg.github.io/bixverse/reference/rs_paga.md).

- cluster_levels:

  Character vector. The cluster levels, in the order the 0-indexed
  partition labels referred to.

- cluster_col:

  String. The obs column the clustering came from.

- modality:

  String. The modality the kNN graph came from.

## Value

Generates the `PagaRes` class.
