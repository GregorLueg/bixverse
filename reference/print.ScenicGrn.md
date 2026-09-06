# Print a ScenicGrn object

Print a ScenicGrn object

## Usage

``` r
# S3 method for class 'ScenicGrn'
print(x, ...)
```

## Arguments

- x:

  A `ScenicGrn` object.

- ...:

  Further arguments passed to or from other methods.

## Value

Invisibly returns `x`.

## Examples

``` r
# what the GRN inference produced and which steps have run
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
print(grn)
#> ScenicGrn (GRN results)
#>   No genes:                 50 
#>   No TFs:                   5 
#>   Applied learner:          randomforest 
#>   TF to gene generated:     FALSE 
#>   CisTarget res generated:  FALSE 

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
