# Map cell ids to the barcodes of a fragment file

Map cell ids to the barcodes of a fragment file

## Usage

``` r
.fragment_barcodes(ids, barcode_map)
```

## Arguments

- ids:

  Character vector. Cell ids in store order.

- barcode_map:

  `NULL`, a named character vector (names are cell ids, values the
  fragment barcodes) or a function taking and returning a character
  vector.

## Value

Character vector of fragment barcodes, one per id.
