# Base class for the single cell DuckDB connection

This is the base class for the single cell experiment DuckDB connection,
containing standard functions such as retrievel of tables, connection
checks, etc.

## Methods

### Public methods

- [`SingleCellDuckDBBase$new()`](#method-SingleCellDuckDBBase-new)

- [`SingleCellDuckDBBase$get_obs_table()`](#method-SingleCellDuckDBBase-get_obs_table)

- [`SingleCellDuckDBBase$get_vars_table()`](#method-SingleCellDuckDBBase-get_vars_table)

- [`SingleCellDuckDBBase$get_obs_index_map()`](#method-SingleCellDuckDBBase-get_obs_index_map)

- [`SingleCellDuckDBBase$get_var_index_map()`](#method-SingleCellDuckDBBase-get_var_index_map)

- [`SingleCellDuckDBBase$get_cells_to_keep()`](#method-SingleCellDuckDBBase-get_cells_to_keep)

- [`SingleCellDuckDBBase$get_obs_cols()`](#method-SingleCellDuckDBBase-get_obs_cols)

- [`SingleCellDuckDBBase$set_cells_to_keep()`](#method-SingleCellDuckDBBase-set_cells_to_keep)

- [`SingleCellDuckDBBase$filter_var_table()`](#method-SingleCellDuckDBBase-filter_var_table)

- [`SingleCellDuckDBBase$add_data_obs()`](#method-SingleCellDuckDBBase-add_data_obs)

- [`SingleCellDuckDBBase$join_data_obs()`](#method-SingleCellDuckDBBase-join_data_obs)

- [`SingleCellDuckDBBase$set_to_keep_column()`](#method-SingleCellDuckDBBase-set_to_keep_column)

- [`SingleCellDuckDBBase$add_data_var()`](#method-SingleCellDuckDBBase-add_data_var)

- [`SingleCellDuckDBBase$rename_column()`](#method-SingleCellDuckDBBase-rename_column)

- [`SingleCellDuckDBBase$clone()`](#method-SingleCellDuckDBBase-clone)

------------------------------------------------------------------------

### Method `new()`

Initialises the Singe Cell DuckDB connection

#### Usage

    SingleCellDuckDBBase$new(db_dir, db_name = "sc_duckdb.db")

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

    SingleCellDuckDBBase$get_obs_table(
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

    SingleCellDuckDBBase$get_vars_table(indices = NULL, cols = NULL)

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

    SingleCellDuckDBBase$get_obs_index_map()

#### Returns

A named numeric containing the cell index mapping.

------------------------------------------------------------------------

### Method `get_var_index_map()`

Returns a mapping between variable index and variable names.

#### Usage

    SingleCellDuckDBBase$get_var_index_map()

#### Returns

A named numeric containing the gene index mapping.

------------------------------------------------------------------------

### Method [`get_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/get_cells_to_keep.md)

Returns the indices of the cells that have to_keep = TRUE in the DB.

#### Usage

    SingleCellDuckDBBase$get_cells_to_keep()

#### Returns

The index positions (1-index) of the cells to keep.

------------------------------------------------------------------------

### Method `get_obs_cols()`

Returns the available column names in the obs table.

#### Usage

    SingleCellDuckDBBase$get_obs_cols()

#### Returns

A character vector of column names. Filter the obs table and reset the
cell idx

------------------------------------------------------------------------

### Method [`set_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/set_cells_to_keep.md)

#### Usage

    SingleCellDuckDBBase$set_cells_to_keep(cell_idx_to_keep)

#### Arguments

- `cell_idx_to_keep`:

  Integer vector with the cell indices to keep. Needs to be 1-indexed!

#### Returns

Invisible self after updating the to_keep column in the DuckDB. Filter
the var table and reset the gene idx

------------------------------------------------------------------------

### Method `filter_var_table()`

#### Usage

    SingleCellDuckDBBase$filter_var_table(filter_vec)

#### Arguments

- `filter_vec`:

  Boolean vector that will be used to filter the var table.

------------------------------------------------------------------------

### Method `add_data_obs()`

Add new data to the obs table in the DuckDB

#### Usage

    SingleCellDuckDBBase$add_data_obs(new_data)

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

    SingleCellDuckDBBase$join_data_obs(new_data)

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

    SingleCellDuckDBBase$set_to_keep_column()

------------------------------------------------------------------------

### Method `add_data_var()`

Add the information which genes pass threshold to the DuckDB.

#### Usage

    SingleCellDuckDBBase$add_data_var(new_data)

#### Arguments

- `new_data`:

  A data.table with additional new columns. The order needs to be the
  same as the original in the var table.

#### Returns

Invisible self while adding the new columns to the var table in the
DuckDB.

------------------------------------------------------------------------

### Method `rename_column()`

Rename a column in the obs or var table.

#### Usage

    SingleCellDuckDBBase$rename_column(table = c("obs", "var"), old, new)

#### Arguments

- `table`:

  String. Either "obs" or "var".

- `old`:

  String. The current column name.

- `new`:

  String. The desired new column name.

#### Returns

Invisible self.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SingleCellDuckDBBase$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
