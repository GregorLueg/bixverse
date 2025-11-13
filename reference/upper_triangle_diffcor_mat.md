# Class for symmetric differential correlation matrices

The class allows to store the upper triangular matrices of of
differential correlation results in an memory-efficient form and return
a data.table or dense R matrix for a given parameter if need be.

## Super class

[`bixverse::upper_triangular_sym_mat`](upper_triangular_sym_mat.md) -\>
`upper_triangle_diffcor_mat`

## Methods

### Public methods

- [`upper_triangle_diffcor_mat$new()`](#method-upper_triangle_diffcor_mat-new)

- [`upper_triangle_diffcor_mat$get_data_table()`](#method-upper_triangle_diffcor_mat-get_data_table)

- [`upper_triangle_diffcor_mat$get_cor_matrix()`](#method-upper_triangle_diffcor_mat-get_cor_matrix)

- [`upper_triangle_diffcor_mat$clone()`](#method-upper_triangle_diffcor_mat-clone)

Inherited methods

- [`bixverse::upper_triangular_sym_mat$get_data()`](upper_triangular_sym_mat.html#method-get_data)
- [`bixverse::upper_triangular_sym_mat$get_sparse_matrix()`](upper_triangular_sym_mat.html#method-get_sparse_matrix)
- [`bixverse::upper_triangular_sym_mat$get_sym_matrix()`](upper_triangular_sym_mat.html#method-get_sym_matrix)
- [`bixverse::upper_triangular_sym_mat$print()`](upper_triangular_sym_mat.html#method-print)

------------------------------------------------------------------------

### Method `new()`

Initialises the R6 class.

#### Usage

    upper_triangle_diffcor_mat$new(diff_cor_res, features)

#### Arguments

- `diff_cor_res`:

  A list of differential correlation results.

- `features`:

  String vector. The features of the correlation matrices.

#### Returns

Returns the initialised class.

------------------------------------------------------------------------

### Method `get_data_table()`

Returns the data in form of a data.table.

#### Usage

    upper_triangle_diffcor_mat$get_data_table(factor = FALSE, .verbose = TRUE)

#### Arguments

- `factor`:

  Boolean. Shall the string columns be transformed into factors. Reduces
  size of the object; however, takes longer to generate.

- `.verbose`:

  Boolean. Controls verbosity.

#### Returns

A data.table with three columns:

- feature_a: The name of the first feature in the correlation matrix.

- feature_b: The name of the second feature in the correlation matrix.

- cor_a: The correlation coefficients for data set a between the
  features.

- cor_b: The correlation coefficients for data set b between the
  features.

- z_score: The differential correlation z-score.

- p_val: The p-value of the differential correlation.

------------------------------------------------------------------------

### Method `get_cor_matrix()`

Return the full correlation matrix of either sample.

#### Usage

    upper_triangle_diffcor_mat$get_cor_matrix(
      to_ret = c("cor_a", "cor_b"),
      .verbose = TRUE
    )

#### Arguments

- `to_ret`:

  String. Option of `("cor_a", "cor_b")`, pending on which correlation
  matrix you want to retrieve.

- `.verbose`:

  Boolean. Controls verbosity.

#### Returns

Returns the specified correlation matrix as a dense R matrix.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    upper_triangle_diffcor_mat$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
