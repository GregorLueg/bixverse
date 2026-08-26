# Validate one cell type's feature matrix and align its rows

Validate one cell type's feature matrix and align its rows

## Usage

``` r
.assert_dialogue_features(feature, cell_ids, cell_type)
```

## Arguments

- feature:

  Numeric matrix. The cell type's features.

- cell_ids:

  Character. The cell type's cell identifiers, in the order the indices
  were resolved in.

- cell_type:

  String. Name of the cell type, for the error message.

## Value

The matrix, rows reordered to `cell_ids`.
