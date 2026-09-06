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

## Examples

``` r
# gene to module assignments from an NMF fit
syn <- synthetic_bulk_cor_matrix()
mat <- log1p(t(syn$counts))
meta <- data.table::data.table(sample_id = rownames(mat))
object <- BulkCoExp(raw_data = mat, meta_data = meta)
object <- preprocess_bulk_coexp(object, hvg = 500L, .verbose = FALSE)
object <- nmf_bulk(object, k = 3L, .verbose = FALSE)
res <- S7::prop(object, "final_results")
head(get_modules(res))
#>        gene module_id  loading   sign        z
#>      <char>    <char>    <num> <char>    <num>
#> 1: gene_107   comp_01 42.79332    pos 3.774586
#> 2: gene_240   comp_02 29.67022    pos 9.529616
#> 3: gene_201   comp_02 27.39015    pos 8.528474
#> 4: gene_269   comp_02 26.93528    pos 8.328749
#> 5: gene_291   comp_02 26.75724    pos 8.250575
#> 6: gene_288   comp_02 26.25955    pos 8.032049
```
