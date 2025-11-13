# Generate a sparse dictionary with DGRDL

This is the Rust implementation of dual graph regularised dictionary
learning in the implementation of Pan, et al., Cell Systems, 2022.

## Usage

``` r
rs_sparse_dict_dgrdl(x, dgrdl_params, seed, verbose)
```

## Arguments

- x:

  Numerical matrix. Rows = samples, columns = features.

- dgrdl_params:

  A list with the parameters for the algorithm. Expects the following
  items.

  - sparsity - Sparsity constraint (max non-zero coefficients per
    signal).

  - dict_size - Size of the dictionary.

  - alpha - Float. Sample context regularisation weight. The higher the
    stronger the regularisation.

  - beta - Float. Feature context regularisation weight. The higher the
    stronger the regularisation.

  - max_iter - Integer. Maximum iteration for the algorithm.

  - k_neighbours - Integer. Number of k neighbours for the sample and
    feature Laplacian matrix for the regularisation

  - admm_iter Integer. Number of iterations for using alternating
    direction method of multipliers (ADMM).

  - rho Float. ADMM step size.

- seed:

  Integer. Seed for the initialisation of the algorithm.

- verbose:

  Boolean. Controls the verbosity of the function and reports timing of
  individual steps.

## Value

A list with the following elements:

- dictionary - The dictionary of samples x dict_size.

- coefficients - The feature loadings of size dict_size x features.

- feature_laplacian - The KNN graph laplacian of the features in a
  sparse format list.

- sample_laplacian - The KNN graph laplacian of the samples in a sparse
  format list.
