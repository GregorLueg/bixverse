# Convert SingleCellNearestNeighbour to manifoldsR NearestNeighbours

Converts the bixverse single cell kNN class into the manifoldsR
`NearestNeighbours` class for use in UMAP, tSNE, etc. Note that the
bixverse class stores 0-indexed indices whereas manifoldsR expects
1-indexed, so the conversion handles that.

## Usage

``` r
sc_knn_to_nearest_neighbours(x)
```

## Arguments

- x:

  `SingleCellNearestNeighbour` class.

## Value

A `NearestNeighbours` class compatible with manifoldsR.
