# Similarity network fusion

This class is designed to generate and store the adjacency matrices for
subsequent usage in similarity network fusion. There are methods to add
different types of adjacency matrices (pending the type of data, it can
be continuous, categorical and/or mixed). Subsequently, there are
methods to generate the final network.

## Usage

``` r
snf(data = NULL, data_name = NULL, snf_params = params_snf())
```

## Arguments

- data:

  Optional data to already transform into adjacency data. Any data
  supplied will be assumed to be samples x features. The provided data
  taype can be a data.table (for categorical and/or mixed types) or a
  matrix (for continous types). If you provide a data.table, the
  function will assume the first column are the sample identifiers.
  Please ensure that setting.

- data_name:

  Optional string. Name of the data modality.

- snf_params:

  List. The SNF parameters, see [`params_snf()`](params_snf.md). The
  list contains the following elements:

  - k - Integer. Number of neighbours to consider.

  - t - Integer. Number of iterations for the SNF algorithm.

  - mu - Float. Normalisation factor for the Gaussian kernel width.

  - alpha - Float. Normalisation parameter controlling the fusion
    strength.

  - normalise - Boolean. Shall continuous values be Z-scored.

  - distance_metric - String. One of
    `c("euclidean", "manhattan", "canberra", "cosine")`. Which distance
    metric to use for the continuous calculations. In case of pure
    categorical, Hamming will be used, for mixed data types Gower
    distance is used.

  The parameters will be internally stored for subsequent usage in other
  functions.

## Value

Returns the `snf` class for further operations.

## Properties

- adj_matrices:

  A list containing the processed adjacency matrices.

- snf_adj:

  Matrix. The final adjacency of the fused network.

- params:

  A (nested) list that will store all the parameters of the applied
  function.

- final_results:

  data.table. Contains final results.
