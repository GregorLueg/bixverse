# Helper function to check if PC1 and/or PC2 distinguish groups

Runs an ANOVA for the group variable vs PC1 and PC2 and checks if either
principal component can separate the groups.

## Usage

``` r
check_pca_grp_differences(pc1, pc2, grps)
```

## Arguments

- pc1:

  Numeric vector. The values of PC1.

- pc2:

  Numeric vector. The values of PC2.

- grps:

  Factor or character vector. The group vector.

## Value

A data.table with two columns
