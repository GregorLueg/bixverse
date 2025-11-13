# Base class for the single cell DuckDB connection

This is the base class for the single cell experiment DuckDB connection,
containing standard functions such as retrievel of tables, connection
checks, etc.

## Methods

### Public methods

- [`single_cell_duckdb_base$new()`](#method-single_cell_duckdb_base-new)

- [`single_cell_duckdb_base$get_obs_table()`](#method-single_cell_duckdb_base-get_obs_table)

- [`single_cell_duckdb_base$get_vars_table()`](#method-single_cell_duckdb_base-get_vars_table)

- [`single_cell_duckdb_base$get_obs_index_map()`](#method-single_cell_duckdb_base-get_obs_index_map)

- [`single_cell_duckdb_base$get_var_index_map()`](#method-single_cell_duckdb_base-get_var_index_map)

- [`single_cell_duckdb_base$get_cells_to_keep()`](#method-single_cell_duckdb_base-get_cells_to_keep)

- [`single_cell_duckdb_base$set_cells_to_keep()`](#method-single_cell_duckdb_base-set_cells_to_keep)

- [`single_cell_duckdb_base$filter_var_table()`](#method-single_cell_duckdb_base-filter_var_table)

- [`single_cell_duckdb_base$add_data_obs()`](#method-single_cell_duckdb_base-add_data_obs)

- [`single_cell_duckdb_base$join_data_obs()`](#method-single_cell_duckdb_base-join_data_obs)

- [`single_cell_duckdb_base$set_to_keep_column()`](#method-single_cell_duckdb_base-set_to_keep_column)

- [`single_cell_duckdb_base$add_data_var()`](#method-single_cell_duckdb_base-add_data_var)

- [`single_cell_duckdb_base$clone()`](#method-single_cell_duckdb_base-clone)

------------------------------------------------------------------------

### Method `new()`

Initialises the Singe Cell DuckDB connection

#### Usage

    single_cell_duckdb_base$new(db_dir, db_name = "sc_duckdb.db")

#### Arguments

- `db_dir`:

  String. Path to where to store the db.

- `db_name`:

  String. The name of the DB. Defaults to `"sc_duckdb.db"`.

#### Returns

Returns the initialised class

------------------------------------------------------------------------

### Method `get_obs_table()`

Returns the full observation table from the DuckDB

#### Usage

    single_cell_duckdb_base$get_obs_table(
      indices = NULL,
      cols = NULL,
      filtered = FALSE
    )

#### Arguments

- `indices`:

  Optional cell/obs indices.

- `cols`:

  Optional column names to return.

- `filtered`:

  Boolean. Whether to return all cells or filtered to to_keep cells.

#### Returns

The observation table (if found) as a data.table with optionally
selected indices and/or columns.

------------------------------------------------------------------------

### Method `get_vars_table()`

Returns the full var table from the DuckDB.

#### Usage

    single_cell_duckdb_base$get_vars_table(indices = NULL, cols = NULL)

#### Arguments

- `indices`:

  Optional gene/var indices.

- `cols`:

  Optional column names to return.

#### Returns

The var table (if found) as a data.table with optionally selected
indices and/or columns.

------------------------------------------------------------------------

### Method `get_obs_index_map()`

Returns a mapping between cell index and cell names/barcodes.

#### Usage

    single_cell_duckdb_base$get_obs_index_map()

#### Returns

A named numeric containing the cell index mapping.

------------------------------------------------------------------------

### Method `get_var_index_map()`

Returns a mapping between variable index and variable names.

#### Usage

    single_cell_duckdb_base$get_var_index_map()

#### Returns

A named numeric containing the gene index mapping.

------------------------------------------------------------------------

### Method [`get_cells_to_keep()`](get_cells_to_keep.md)

Returns the indices of the cells that have to_keep = TRUE in the DB.

#### Usage

    single_cell_duckdb_base$get_cells_to_keep()

#### Returns

The index positions (1-index) of the cells to keep. Filter the obs table
and reset the cell idx

------------------------------------------------------------------------

### Method [`set_cells_to_keep()`](set_cells_to_keep.md)

#### Usage

    single_cell_duckdb_base$set_cells_to_keep(cell_idx_to_keep)

#### Arguments

- `cell_idx_to_keep`:

  Integer vector with the cell indices to keep. Needs to be 1-indexed!

#### Returns

Invisible self after updating the to_keep column in the DuckDB. Filter
the var table and reset the gene idx

------------------------------------------------------------------------

### Method `filter_var_table()`

#### Usage

    single_cell_duckdb_base$filter_var_table(filter_vec)

#### Arguments

- `filter_vec`:

  Boolean vector that will be used to filter the var table.

------------------------------------------------------------------------

### Method `add_data_obs()`

Add new data to the obs table in the DuckDB

#### Usage

    single_cell_duckdb_base$add_data_obs(new_data)

#### Arguments

- `new_data`:

  A data.table with additional new columns. The order needs to be the
  same as the original in the obs table.

#### Returns

Invisible self while adding the new columns to the obs table in the
DuckDB.

------------------------------------------------------------------------

### Method `join_data_obs()`

Left join new data to the obs table in the DuckDB by cell_idx

#### Usage

    single_cell_duckdb_base$join_data_obs(new_data)

#### Arguments

- `new_data`:

  A data.table with a cell_idx column to join on.

#### Returns

Invisible self while left joining the new data to the obs table in the
DuckDB.

------------------------------------------------------------------------

### Method `set_to_keep_column()`

Indepenent of the loader, set the to_keep column to `TRUE` initially

#### Usage

    single_cell_duckdb_base$set_to_keep_column()

------------------------------------------------------------------------

### Method `add_data_var()`

Add the information which genes pass threshold to the DuckDB.

#### Usage

    single_cell_duckdb_base$add_data_var(new_data)

#### Arguments

- `new_data`:

  A data.table with additional new columns. The order needs to be the
  same as the original in the var table.

#### Returns

Invisible self while adding the new columns to the var table in the
DuckDB.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    single_cell_duckdb_base$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
