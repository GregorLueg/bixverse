# Pull covariate columns out of the observation table

The design takes numerics only. A factor would need dummy coding, which
changes the rank of the design matrix, so it is refused rather than
guessed at. The library size is never a covariate: it enters the model
as a fixed offset, which is what separates v2 from v1.

## Usage

``` r
.residual_covariates(obs, covariate_columns, method)
```

## Arguments

- obs:

  data.table. The observation table for the selected cells.

- covariate_columns:

  Character or `NULL`. Columns to use.

- method:

  String. The method being fitted.

## Value

A named list of numeric vectors, empty when there are none.
