# Class for symmetric correlation matrices

The class allows to store the upper triangular matrix of a symmetric
matrix (think correlation matrix, distance matrix, etc.) in an
memory-efficient form and return a data.table or dense (or sparse) R
matrix if need be.

## Methods

### Public methods

- [`upper_triangular_sym_mat$new()`](#method-upper_triangular_sym_mat-initialize)

- [`upper_triangular_sym_mat$print()`](#method-upper_triangular_sym_mat-print)

- [`upper_triangular_sym_mat$get_data_table()`](#method-upper_triangular_sym_mat-get_data_table)

- [`upper_triangular_sym_mat$get_sym_matrix()`](#method-upper_triangular_sym_mat-get_sym_matrix)

- [`upper_triangular_sym_mat$get_sparse_matrix()`](#method-upper_triangular_sym_mat-get_sparse_matrix)

- [`upper_triangular_sym_mat$get_data()`](#method-upper_triangular_sym_mat-get_data)

- [`upper_triangular_sym_mat$clone()`](#method-upper_triangular_sym_mat-clone)

------------------------------------------------------------------------

### `upper_triangular_sym_mat$new()`

Initialises the R6 class.

#### Usage

    upper_triangular_sym_mat$new(values, features, shift)

#### Arguments

- `values`:

  Numerical vector. The correlation coefficients of the upper triangular
  correlation matrix stored as a row-major vector

- `features`:

  String vector. The features of the correlation matrix.

- `shift`:

  Bollean Was a shift applied during the calculation of the upper
  triangular matrix. If `FALSE`, the diagonal was included, if `TRUE`,
  the diagonal was removed.

#### Returns

Returns the initialised class.

------------------------------------------------------------------------

### `upper_triangular_sym_mat$print()`

Print the class

#### Usage

    upper_triangular_sym_mat$print()

#### Returns

Returns the initialised class

------------------------------------------------------------------------

### `upper_triangular_sym_mat$get_data_table()`

Returns the data in form of a data.table.

#### Usage

    upper_triangular_sym_mat$get_data_table(factor = FALSE, .verbose = TRUE)

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

- cor: The correlation coefficients between these two features.

------------------------------------------------------------------------

### `upper_triangular_sym_mat$get_sym_matrix()`

Return the full correlation matrix.

#### Usage

    upper_triangular_sym_mat$get_sym_matrix(.verbose = TRUE)

#### Arguments

- `.verbose`:

  Boolean. Controls verbosity.

#### Returns

Returns the correlation matrix as a dense R matrix.

------------------------------------------------------------------------

### `upper_triangular_sym_mat$get_sparse_matrix()`

Return a sparse version of the correlation matrix

#### Usage

    upper_triangular_sym_mat$get_sparse_matrix(.verbose = TRUE)

#### Arguments

- `.verbose`:

  Boolean. Controls verbosity

#### Returns

The sparse matrix.

------------------------------------------------------------------------

### `upper_triangular_sym_mat$get_data()`

Return the correlation data and shift

#### Usage

    upper_triangular_sym_mat$get_data()

#### Returns

A list with

- cor_data - Numeric vector. The values.

- features - String. The feature names.

- n_features - Integer. Number of initial features.

- shift - Integer. The applied shift.

------------------------------------------------------------------------

### `upper_triangular_sym_mat$clone()`

The objects of this class are cloneable with this method.

#### Usage

    upper_triangular_sym_mat$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# store the off-diagonal of a correlation matrix and get it back
set.seed(42)
cor_mat <- cor(matrix(rnorm(80), nrow = 20, ncol = 4))
object <- upper_triangular_sym_mat$new(
  values = cor_mat[lower.tri(cor_mat)],
  features = sprintf("gene_%i", 1:4),
  shift = TRUE
)
object$get_data_table(.verbose = FALSE)
#>    feature_a feature_b         sim
#>       <char>    <char>       <num>
#> 1:    gene_1    gene_2  0.43694936
#> 2:    gene_1    gene_3  0.04324370
#> 3:    gene_1    gene_4 -0.03960938
#> 4:    gene_2    gene_3  0.03976509
#> 5:    gene_2    gene_4 -0.06555478
#> 6:    gene_3    gene_4  0.33289052
dim(object$get_sym_matrix(.verbose = FALSE))
#> [1] 4 4
```
