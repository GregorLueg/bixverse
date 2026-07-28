# Get the module membership from a BulkModuleResult

Returns the data.table of gene to module_id assignments. The exact
columns depend on the method that produced the result (CoReMo adds
`sign` and `stability`; NMF/ICA/DGRDL add `loading`, `sign` and the
thresholding score; Leiden adds only `module_id`).

`gene` is **not** a unique key for the matrix factorisation methods.
ICA, NMF and DGRDL assign membership by keeping the tails of each
component's loading distribution, so a gene loading strongly on three
components appears in three rows, and a gene in no tail appears in none.
That is the point of a factorisation. The partition-based methods
(CoReMo, Leiden) do emit one row per gene. Do not assume uniqueness
without checking `method`.

## Usage

``` r
get_modules(object)
```

## Arguments

- object:

  A `BulkModuleResult`.

## Value

data.table with at minimum `gene` and `module_id` columns. One row per
(gene, module) pair.
