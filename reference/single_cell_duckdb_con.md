# Class for storing single cell experimental data in DuckDB (nightly!)

This class wraps up the DB connection and methods to interact with the
observation/cell metadata and the feature/gene metadata.

## Value

Invisible self and populates the internal obs table. Function to
populate the var table from plain text files

Invisible self and populates the internal var table. Function to
populate the obs table from R

Invisible self and populates the internal obs table. Function to
populate the var table from R

Invisible self and populates the internal obs table.

## Super class

[`bixverse::single_cell_duckdb_base`](single_cell_duckdb_base.md) -\>
`single_cell_duckdb_con`

## Methods

### Public methods

- [`single_cell_duckdb_con$populate_obs_from_h5()`](#method-single_cell_duckdb_con-populate_obs_from_h5)

- [`single_cell_duckdb_con$populate_vars_from_h5()`](#method-single_cell_duckdb_con-populate_vars_from_h5)

- [`single_cell_duckdb_con$populate_obs_from_plain_text()`](#method-single_cell_duckdb_con-populate_obs_from_plain_text)

- [`single_cell_duckdb_con$populate_var_from_plain_text()`](#method-single_cell_duckdb_con-populate_var_from_plain_text)

- [`single_cell_duckdb_con$populate_obs_from_data.table()`](#method-single_cell_duckdb_con-populate_obs_from_data.table)

- [`single_cell_duckdb_con$populate_var_from_data.table()`](#method-single_cell_duckdb_con-populate_var_from_data.table)

- [`single_cell_duckdb_con$clone()`](#method-single_cell_duckdb_con-clone)

Inherited methods

- [`bixverse::single_cell_duckdb_base$add_data_obs()`](single_cell_duckdb_base.html#method-add_data_obs)
- [`bixverse::single_cell_duckdb_base$add_data_var()`](single_cell_duckdb_base.html#method-add_data_var)
- [`bixverse::single_cell_duckdb_base$filter_var_table()`](single_cell_duckdb_base.html#method-filter_var_table)
- [`bixverse::single_cell_duckdb_base$get_cells_to_keep()`](single_cell_duckdb_base.html#method-get_cells_to_keep)
- [`bixverse::single_cell_duckdb_base$get_obs_index_map()`](single_cell_duckdb_base.html#method-get_obs_index_map)
- [`bixverse::single_cell_duckdb_base$get_obs_table()`](single_cell_duckdb_base.html#method-get_obs_table)
- [`bixverse::single_cell_duckdb_base$get_var_index_map()`](single_cell_duckdb_base.html#method-get_var_index_map)
- [`bixverse::single_cell_duckdb_base$get_vars_table()`](single_cell_duckdb_base.html#method-get_vars_table)
- [`bixverse::single_cell_duckdb_base$initialize()`](single_cell_duckdb_base.html#method-initialize)
- [`bixverse::single_cell_duckdb_base$join_data_obs()`](single_cell_duckdb_base.html#method-join_data_obs)
- [`bixverse::single_cell_duckdb_base$set_cells_to_keep()`](single_cell_duckdb_base.html#method-set_cells_to_keep)
- [`bixverse::single_cell_duckdb_base$set_to_keep_column()`](single_cell_duckdb_base.html#method-set_to_keep_column)

------------------------------------------------------------------------

### Method `populate_obs_from_h5()`

This function populates the obs table from an h5 file (if found).

#### Usage

    single_cell_duckdb_con$populate_obs_from_h5(h5_path, filter = NULL)

#### Arguments

- `h5_path`:

  String. Path to the h5 file from which to load in the observation
  table.

- `filter`:

  Optional integer. Positions of obs to read in from file.

#### Returns

Returns invisible self. As a side effect, it will load in the obs data
from the h5ad file into the DuckDB.

------------------------------------------------------------------------

### Method `populate_vars_from_h5()`

This populates the vars table from an h5 file (if found.)

#### Usage

    single_cell_duckdb_con$populate_vars_from_h5(h5_path, filter = NULL)

#### Arguments

- `h5_path`:

  String. Path to the h5 file from which to load in the observation
  table.

- `filter`:

  Optional integer. Positions of obs to read in from file.

#### Returns

Returns invisible self. As a side effect, it will load in the obs data
from the h5ad file into the DuckDB. Function to populate the obs table
from plain text files

------------------------------------------------------------------------

### Method `populate_obs_from_plain_text()`

#### Usage

    single_cell_duckdb_con$populate_obs_from_plain_text(
      f_path,
      has_hdr,
      filter = NULL
    )

#### Arguments

- `f_path`:

  File path to the plain text file.

- `has_hdr`:

  Boolean. Does the flat file have a header.

- `filter`:

  Optional integer. Positions of obs to read in from file.

------------------------------------------------------------------------

### Method `populate_var_from_plain_text()`

#### Usage

    single_cell_duckdb_con$populate_var_from_plain_text(
      f_path,
      has_hdr,
      filter = NULL
    )

#### Arguments

- `f_path`:

  File path to the plain text file.

- `has_hdr`:

  Boolean. Does the flat file have a header.

- `filter`:

  Optional integer. Positions of obs to read in from file.

------------------------------------------------------------------------

### Method `populate_obs_from_data.table()`

#### Usage

    single_cell_duckdb_con$populate_obs_from_data.table(obs_dt, filter = NULL)

#### Arguments

- `obs_dt`:

  data.table

- `filter`:

  Optional integer. Row indices to keep

------------------------------------------------------------------------

### Method `populate_var_from_data.table()`

#### Usage

    single_cell_duckdb_con$populate_var_from_data.table(var_dt, filter = NULL)

#### Arguments

- `var_dt`:

  data.table

- `filter`:

  Optional integer. Row indices to keep

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    single_cell_duckdb_con$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
