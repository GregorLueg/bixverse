# Copy an object's cell filter into its default gene activity store

The gene activity store at `<dir>/gene_activity` shares the cell axis of
the experiment; without this, a filter set after it was written would
not reach it.

## Usage

``` r
.sync_gene_activity_filter(x)
```

## Arguments

- x:

  The single cell object owning the gene activity.

## Value

Invisible `NULL`.
