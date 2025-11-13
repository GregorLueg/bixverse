# Class for Anndata

This class helps dealing with h5ad objects from Python AnnData. You have
several options will allow for retrieval of the underlying data.

## Methods

### Public methods

- [`anndata_parser$new()`](#method-anndata_parser-new)

- [`anndata_parser$get_obs_table()`](#method-anndata_parser-get_obs_table)

- [`anndata_parser$get_var_info()`](#method-anndata_parser-get_var_info)

- [`anndata_parser$get_raw_counts()`](#method-anndata_parser-get_raw_counts)

- [`anndata_parser$get_key_data()`](#method-anndata_parser-get_key_data)

- [`anndata_parser$clone()`](#method-anndata_parser-clone)

------------------------------------------------------------------------

### Method `new()`

Initialises the Anndata Parser.

#### Usage

    anndata_parser$new(h5_path)

#### Arguments

- `h5_path`:

  String. Path to the h5 file.

#### Returns

Returns the initialised class.

------------------------------------------------------------------------

### Method `get_obs_table()`

Returns the observation table with all the data from the h5ad file.

#### Usage

    anndata_parser$get_obs_table()

#### Returns

data.table. The found observations are returned. The pandas index will
be named `sample_id`. Remaining columns (if found) will be returned as
factors due to the way the data is stored in h5.

------------------------------------------------------------------------

### Method `get_var_info()`

Returns the variable table with all the data from the h5ad file.

#### Usage

    anndata_parser$get_var_info()

#### Returns

data.table. The found observations are returned. The pandas index will
be named `var_id`. Remaining columns (if found) will be returned as
factors due to the way the data is stored in h5.

------------------------------------------------------------------------

### Method `get_raw_counts()`

Returns the counts that are stored in `X` slot of the anndata object.

#### Usage

    anndata_parser$get_raw_counts()

#### Returns

Returns the count matrix with samples = columns and rows = features.

------------------------------------------------------------------------

### Method `get_key_data()`

Wrapper function that returns a list of the stored count data and the
metadata found in the h5ad file.

#### Usage

    anndata_parser$get_key_data()

#### Returns

List with following elements:

- metadata - metadata from the respective h5ad file

- var_info - metadata on the variables from the respective h5ad file.

- counts - counts that were found in the h5ad file.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    anndata_parser$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
