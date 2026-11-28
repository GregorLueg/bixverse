# What the features of an object's store are

What the features of an object's store are

## Usage

``` r
.store_feature_kind(object)
```

## Arguments

- object:

  The single cell object.

## Value

String. One of `"gene"`, `"peak"`, `"tile"`, `"snp"`, `"protein"` or
`"gene_activity"`. Meta cells carry no store and read as `"gene"`.
