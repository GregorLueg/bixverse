# Extract the TF to gene data from the ScenicGrn object

Extract the TF to gene data from the ScenicGrn object

## Usage

``` r
get_tf_to_gene(x)

# S3 method for class 'ScenicGrn'
get_tf_to_gene(x)
```

## Arguments

- x:

  `ScenicGrn` object from which to extract the TF to gene data.table.

## Value

data.table with TF to gene information

## Examples

``` r
# the TF to gene links that survived the importance filter
sc <- demo_single_cells()
grn <- scenic_grn_sc(
  sc,
  tf_ids = sprintf("gene_%02d", 1:5),
  scenic_params = params_scenic(
    min_counts = 1L,
    learner_params = list(n_trees = 20L)
  ),
  .verbose = FALSE
)
grn <- identify_tf_to_genes(
  grn,
  method = "top_k",
  k_tfs = 3L,
  .verbose = FALSE
)
head(get_tf_to_gene(grn))
#>         tf    gene importance
#>     <char>  <char>      <num>
#> 1: gene_01 gene_01  0.4643152
#> 2: gene_03 gene_01  0.3568220
#> 3: gene_04 gene_01  0.1553617
#> 4: gene_02 gene_02  0.4926334
#> 5: gene_03 gene_02  0.2651168
#> 6: gene_04 gene_02  0.1317756

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
