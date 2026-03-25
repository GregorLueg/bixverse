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

[`bixverse::SingleCellDuckDBBase`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.md)
-\> `SingleCellDuckDB`

## Methods

### Public methods

- [`SingleCellDuckDB$populate_obs_from_h5()`](#method-SingleCellDuckDB-populate_obs_from_h5)

- [`SingleCellDuckDB$populate_vars_from_h5()`](#method-SingleCellDuckDB-populate_vars_from_h5)

- [`SingleCellDuckDB$populate_obs_from_multi_h5()`](#method-SingleCellDuckDB-populate_obs_from_multi_h5)

- [`SingleCellDuckDB$populate_vars_from_h5_reordered()`](#method-SingleCellDuckDB-populate_vars_from_h5_reordered)

- [`SingleCellDuckDB$populate_obs_from_plain_text()`](#method-SingleCellDuckDB-populate_obs_from_plain_text)

- [`SingleCellDuckDB$populate_var_from_plain_text()`](#method-SingleCellDuckDB-populate_var_from_plain_text)

- [`SingleCellDuckDB$populate_obs_from_data.table()`](#method-SingleCellDuckDB-populate_obs_from_data.table)

- [`SingleCellDuckDB$populate_var_from_data.table()`](#method-SingleCellDuckDB-populate_var_from_data.table)

- [`SingleCellDuckDB$clone()`](#method-SingleCellDuckDB-clone)

Inherited methods

- [`bixverse::SingleCellDuckDBBase$add_data_obs()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-add_data_obs)
- [`bixverse::SingleCellDuckDBBase$add_data_var()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-add_data_var)
- [`bixverse::SingleCellDuckDBBase$filter_var_table()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-filter_var_table)
- [`bixverse::SingleCellDuckDBBase$get_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-get_cells_to_keep)
- [`bixverse::SingleCellDuckDBBase$get_obs_cols()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-get_obs_cols)
- [`bixverse::SingleCellDuckDBBase$get_obs_index_map()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-get_obs_index_map)
- [`bixverse::SingleCellDuckDBBase$get_obs_table()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-get_obs_table)
- [`bixverse::SingleCellDuckDBBase$get_var_index_map()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-get_var_index_map)
- [`bixverse::SingleCellDuckDBBase$get_vars_table()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-get_vars_table)
- [`bixverse::SingleCellDuckDBBase$initialize()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-initialize)
- [`bixverse::SingleCellDuckDBBase$join_data_obs()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-join_data_obs)
- [`bixverse::SingleCellDuckDBBase$rename_column()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-rename_column)
- [`bixverse::SingleCellDuckDBBase$set_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-set_cells_to_keep)
- [`bixverse::SingleCellDuckDBBase$set_to_keep_column()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-set_to_keep_column)

------------------------------------------------------------------------

### Method `populate_obs_from_h5()`

This function populates the obs table from an h5 file (if found).

#### Usage

    SingleCellDuckDB$populate_obs_from_h5(
      h5_path,
      filter = NULL,
      cell_id_col = NULL
    )

#### Arguments

- `h5_path`:

  String. Path to the h5 file from which to load in the observation
  table.

- `filter`:

  Optional integer. Positions of obs to read in from file.

- `cell_id_col`:

  Optional string. If you can see that the first column in the h5ad obs
  table is NOT the cell identifier, you can provide here the correct
  column name to use as cell_name.

#### Returns

Returns invisible self. As a side effect, it will load in the obs data
from the h5ad file into the DuckDB.

------------------------------------------------------------------------

### Method `populate_vars_from_h5()`

This populates the vars table from an h5 file (if found.)

#### Usage

    SingleCellDuckDB$populate_vars_from_h5(h5_path, filter = NULL)

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

### Method `populate_obs_from_multi_h5()`

Populate the obs table from multiple h5ad files.

#### Usage

    SingleCellDuckDB$populate_obs_from_multi_h5(per_file_info, cell_id_col = NULL)

#### Arguments

- `per_file_info`:

  List of lists, each containing: h5_path, exp_id, cell_filter
  (1-indexed integer vector of cells to keep).

- `cell_id_col`:

  Optional string. Column name for cell identifiers.

#### Returns

Invisible self. Populates the obs table in DuckDB.

------------------------------------------------------------------------

### Method `populate_vars_from_h5_reordered()`

Populate the var table from an h5ad file, filtered and reordered to
match a target gene set.

#### Usage

    SingleCellDuckDB$populate_vars_from_h5_reordered(h5_path, final_gene_names)

#### Arguments

- `h5_path`:

  String. Path to the reference h5ad file.

- `final_gene_names`:

  Character vector. Gene names in the desired final order.

#### Returns

Invisible self. Populates the var table in DuckDB. Function to populate
the obs table from plain text files

------------------------------------------------------------------------

### Method `populate_obs_from_plain_text()`

#### Usage

    SingleCellDuckDB$populate_obs_from_plain_text(f_path, has_hdr, filter = NULL)

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

    SingleCellDuckDB$populate_var_from_plain_text(f_path, has_hdr, filter = NULL)

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

    SingleCellDuckDB$populate_obs_from_data.table(obs_dt, filter = NULL)

#### Arguments

- `obs_dt`:

  data.table

- `filter`:

  Optional integer. Row indices to keep

------------------------------------------------------------------------

### Method `populate_var_from_data.table()`

#### Usage

    SingleCellDuckDB$populate_var_from_data.table(var_dt, filter = NULL)

#### Arguments

- `var_dt`:

  data.table

- `filter`:

  Optional integer. Row indices to keep

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    SingleCellDuckDB$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
