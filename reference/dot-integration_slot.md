# Helper to fetch an integration result slot from other_data

Helper to fetch an integration result slot from other_data

## Usage

``` r
.integration_slot(x, modality)
```

## Arguments

- x:

  `SingleCellsMultiModal` class.

- modality:

  String. The integration slot name, e.g. `"wnn"`.

## Value

The integration result list stored under `other_data[[modality]]`.
