# Calculate PCA on the expression.

Calculates the principal component on top of the filtered count matrix
and adds the information of the first 10 principal components to the
outputs.

## Usage

``` r
calculate_pca_bulk_dge(
  object,
  scale_genes = FALSE,
  pcs = 10L,
  no_hvg_genes = 2500L
)
```

## Arguments

- object:

  The underlying class, see
  [`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md).

- scale_genes:

  Boolean. Shall the log(cpm) counts be scaled prior the PCA
  calculation. Defaults to `FALSE`.

- pcs:

  Integer. Number of PCs to return and add to the outputs slot.

- no_hvg_genes:

  Integer. Number of highly variable genes to include. Defaults to 2500.

## Value

Returns the class with additional data added to the outputs.

## Examples

``` r
# PCA over the most variable genes of the normalised counts
syn <- synthetic_bulk_cor_matrix()
meta <- data.table::data.table(
  sample_id = colnames(syn$counts),
  case_control = rep(c("case", "control"), each = 50)
)
object <- BulkDge(raw_counts = syn$counts, meta_data = meta)
object <- qc_bulk_dge(object, group_col = "case_control", .verbose = FALSE)
object <- normalise_bulk_dge(
  object,
  group_col = "case_control",
  .verbose = FALSE
)
#> calcNormFactors has been renamed to normLibSizes
object <- calculate_pca_bulk_dge(object, no_hvg_genes = 500L)
head(get_outputs(object)$pca)
#>     sample_id      PC_1       PC_2       PC_3       PC_4        PC_5
#>        <char>     <num>      <num>      <num>      <num>       <num>
#> 1:   sample_1  9.948872  -9.247662  4.4376363 -1.1538336  0.07117582
#> 2:  sample_10 -8.356527 -11.333461  0.5720476 -2.9418083 -0.39054751
#> 3: sample_100 15.136687  -3.158151 -1.7093319  3.0947308 -1.37346954
#> 4:  sample_11  6.557195 -14.894399  0.9352057 -2.3963880  1.47430400
#> 5:  sample_12  6.294540  15.042584  6.9301394 -1.6259489 -0.85228072
#> 6:  sample_13  6.874570   1.309428  8.4740366 -0.1374491 -1.74678772
#>           PC_6       PC_7       PC_8        PC_9       PC_10
#>          <num>      <num>      <num>       <num>       <num>
#> 1: -0.83484905  0.5471809 -0.9427867  0.68840292 -0.63813058
#> 2:  0.42850161  4.5141059  0.8725173  0.09963631 -0.16641037
#> 3: -3.24209198 -4.1073248 -1.4664193 -0.47897795 -0.03883322
#> 4: -0.43851199  1.7793772 -0.5956544  0.85235802 -0.12901565
#> 5: -1.77071439  1.4409903 -1.6038953 -0.48971647  2.61043156
#> 6: -0.02339393 -2.2806380 -1.1023069  0.10153030 -3.30639581
```
