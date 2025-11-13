# Class for symmetric correlation matrices

The class allows to store the upper triangular matrix of a symmetric
matrix (think correlation matrix, distance matrix, etc.) in an
memory-efficient form and return a data.table or dense (or sparse) R
matrix if need be.

## Methods

### Public methods

- [`upper_triangular_sym_mat$new()`](#method-upper_triangular_sym_mat-new)

- [`upper_triangular_sym_mat$print()`](#method-upper_triangular_sym_mat-print)

- [`upper_triangular_sym_mat$get_data_table()`](#method-upper_triangular_sym_mat-get_data_table)

- [`upper_triangular_sym_mat$get_sym_matrix()`](#method-upper_triangular_sym_mat-get_sym_matrix)

- [`upper_triangular_sym_mat$get_sparse_matrix()`](#method-upper_triangular_sym_mat-get_sparse_matrix)

- [`upper_triangular_sym_mat$get_data()`](#method-upper_triangular_sym_mat-get_data)

- [`upper_triangular_sym_mat$clone()`](#method-upper_triangular_sym_mat-clone)

------------------------------------------------------------------------

### Method `new()`

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

  Integer. Was a shift applied during the calculation of the upper
  triangular matrix. Typically 0 (diagonal included) or 1 (diagonal not
  included)

#### Returns

Returns the initialised class.

------------------------------------------------------------------------

### Method [`print()`](https://rdrr.io/r/base/print.html)

Print the class

#### Usage

    upper_triangular_sym_mat$print()

#### Returns

Returns the initialised class

------------------------------------------------------------------------

### Method `get_data_table()`

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

### Method `get_sym_matrix()`

Return the full correlation matrix.

#### Usage

    upper_triangular_sym_mat$get_sym_matrix(.verbose = TRUE)

#### Arguments

- `.verbose`:

  Boolean. Controls verbosity.

#### Returns

Returns the correlation matrix as a dense R matrix.

------------------------------------------------------------------------

### Method `get_sparse_matrix()`

Return a sparse version of the correlation matrix

#### Usage

    upper_triangular_sym_mat$get_sparse_matrix(.verbose = TRUE)

#### Arguments

- `.verbose`:

  Boolean. Controls verbosity

#### Returns

The sparse matrix.

------------------------------------------------------------------------

### Method `get_data()`

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

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    upper_triangular_sym_mat$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
