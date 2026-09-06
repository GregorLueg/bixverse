# Class for Anndata

This class helps dealing with h5ad objects from Python AnnData. You have
several options will allow for retrieval of the underlying data.

## Methods

### Public methods

- [`AnnDataParser$new()`](#method-AnnDataParser-initialize)

- [`AnnDataParser$get_obs_table()`](#method-AnnDataParser-get_obs_table)

- [`AnnDataParser$get_var_info()`](#method-AnnDataParser-get_var_info)

- [`AnnDataParser$get_raw_counts()`](#method-AnnDataParser-get_raw_counts)

- [`AnnDataParser$get_key_data()`](#method-AnnDataParser-get_key_data)

- [`AnnDataParser$clone()`](#method-AnnDataParser-clone)

------------------------------------------------------------------------

### `AnnDataParser$new()`

Initialises the Anndata Parser.

#### Usage

    AnnDataParser$new(h5_path)

#### Arguments

- `h5_path`:

  String. Path to the h5 file.

#### Returns

Returns the initialised class.

------------------------------------------------------------------------

### `AnnDataParser$get_obs_table()`

Returns the observation table with all the data from the h5ad file.

#### Usage

    AnnDataParser$get_obs_table()

#### Returns

data.table. The found observations are returned. The pandas index will
be named `sample_id`. Remaining columns (if found) will be returned as
factors due to the way the data is stored in h5.

------------------------------------------------------------------------

### `AnnDataParser$get_var_info()`

Returns the variable table with all the data from the h5ad file.

#### Usage

    AnnDataParser$get_var_info()

#### Returns

data.table. The found observations are returned. The pandas index will
be named `var_id`. Remaining columns (if found) will be returned as
factors due to the way the data is stored in h5.

------------------------------------------------------------------------

### `AnnDataParser$get_raw_counts()`

Returns the counts that are stored in `X` slot of the anndata object.

#### Usage

    AnnDataParser$get_raw_counts()

#### Returns

Returns the count matrix with samples = columns and rows = features.

------------------------------------------------------------------------

### `AnnDataParser$get_key_data()`

Wrapper function that returns a list of the stored count data and the
metadata found in the h5ad file.

#### Usage

    AnnDataParser$get_key_data()

#### Returns

List with following elements:

- metadata - metadata from the respective h5ad file

- var_info - metadata on the variables from the respective h5ad file.

- counts - counts that were found in the h5ad file.

------------------------------------------------------------------------

### `AnnDataParser$clone()`

The objects of this class are cloneable with this method.

#### Usage

    AnnDataParser$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# parse counts and obs back out of a small h5ad file
set.seed(42)
counts <- matrix(rpois(50, 5), nrow = 10, ncol = 5)
obs <- data.table::data.table(sample_id = sprintf("cell_%i", 1:10))
var <- data.table::data.table(var_id = sprintf("gene_%i", 1:5))
h5_path <- tempfile(fileext = ".h5ad")
write_h5ad_sc_dense(h5_path, counts, obs, var, .verbose = FALSE)
parser <- AnnDataParser$new(h5_path)
dim(parser$get_raw_counts())
#> [1]  5 10
head(parser$get_obs_table())
#>    sample_id
#>       <char>
#> 1:    cell_1
#> 2:    cell_2
#> 3:    cell_3
#> 4:    cell_4
#> 5:    cell_5
#> 6:    cell_6
unlink(h5_path)
```
