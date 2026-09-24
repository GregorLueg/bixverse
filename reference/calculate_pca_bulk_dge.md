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
object <- calculate_pca_bulk_dge(object, no_hvg_genes = 500L)
head(get_outputs(object)$pca)
#>     sample_id      PC_1        PC_2       PC_3       PC_4      PC_5       PC_6
#>        <char>     <num>       <num>      <num>      <num>     <num>      <num>
#> 1:   sample_1 11.414541 -13.0043249  5.0594914 -1.3855692  2.028329 -0.9098379
#> 2:  sample_10 -8.725271 -11.5502655 -1.4669694 -4.3679483 -2.154628 -0.4249633
#> 3: sample_100 14.554498  -2.8328726 -2.1071670  3.4928711  6.069259  1.8988263
#> 4:  sample_11  6.563648 -17.6102384 -0.8532572 -3.2007390 -1.690664  0.3787470
#> 5:  sample_12  6.538132  14.1635016 10.2197952 -1.6011918  1.140860  0.7129805
#> 6:  sample_13  6.941648  -0.4137485  9.0801575  0.3289563  2.110583 -0.6537719
#>          PC_7       PC_8       PC_9      PC_10
#>         <num>      <num>      <num>      <num>
#> 1:  0.1849586 -2.6250983  0.7826299  2.2234994
#> 2:  3.3238704 -0.2113417 -0.3072612 -0.8666567
#> 3: -1.1091071  1.3089212 -1.4279141  2.3160770
#> 4:  2.7972027 -0.2006078 -0.7214862 -0.6428317
#> 5:  2.5048927 -0.5685659 -2.0599492 -2.2567787
#> 6: -0.6846346 -0.6171280 -5.7555942  0.4322864
```
