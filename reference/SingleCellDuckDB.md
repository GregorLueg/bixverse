# Class for storing single cell experimental data in DuckDB (nightly!)

This class wraps up the DB connection and methods to interact with the
observation/cell metadata and the feature/gene metadata.

## Super class

[`SingleCellDuckDBBase`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.md)
-\> `SingleCellDuckDB`

## Methods

### Public methods

- [`SingleCellDuckDB$populate_obs_from_h5ad()`](#method-SingleCellDuckDB-populate_obs_from_h5ad)

- [`SingleCellDuckDB$populate_vars_from_h5ad()`](#method-SingleCellDuckDB-populate_vars_from_h5ad)

- [`SingleCellDuckDB$populate_obs_from_multi_h5ad()`](#method-SingleCellDuckDB-populate_obs_from_multi_h5ad)

- [`SingleCellDuckDB$populate_var_minimal()`](#method-SingleCellDuckDB-populate_var_minimal)

- [`SingleCellDuckDB$populate_obs_from_plain_text()`](#method-SingleCellDuckDB-populate_obs_from_plain_text)

- [`SingleCellDuckDB$populate_var_from_plain_text()`](#method-SingleCellDuckDB-populate_var_from_plain_text)

- [`SingleCellDuckDB$populate_obs_from_data.table()`](#method-SingleCellDuckDB-populate_obs_from_data.table)

- [`SingleCellDuckDB$populate_var_from_data.table()`](#method-SingleCellDuckDB-populate_var_from_data.table)

- [`SingleCellDuckDB$populate_var_adt_from_data.table()`](#method-SingleCellDuckDB-populate_var_adt_from_data.table)

- [`SingleCellDuckDB$populate_obs_from_tenx_h5()`](#method-SingleCellDuckDB-populate_obs_from_tenx_h5)

- [`SingleCellDuckDB$populate_vars_from_tenx_h5()`](#method-SingleCellDuckDB-populate_vars_from_tenx_h5)

- [`SingleCellDuckDB$populate_obs_from_multi_tenx_h5()`](#method-SingleCellDuckDB-populate_obs_from_multi_tenx_h5)

- [`SingleCellDuckDB$populate_obs_from_multi_duckdb()`](#method-SingleCellDuckDB-populate_obs_from_multi_duckdb)

- [`SingleCellDuckDB$populate_vars_from_duckdb_reordered()`](#method-SingleCellDuckDB-populate_vars_from_duckdb_reordered)

- [`SingleCellDuckDB$populate_obs_from_multi_plain_text()`](#method-SingleCellDuckDB-populate_obs_from_multi_plain_text)

- [`SingleCellDuckDB$clone()`](#method-SingleCellDuckDB-clone)

Inherited methods

- [`SingleCellDuckDBBase$add_data_obs()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-add_data_obs)
- [`SingleCellDuckDBBase$add_data_var()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-add_data_var)
- [`SingleCellDuckDBBase$drop_columns()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-drop_columns)
- [`SingleCellDuckDBBase$filter_var_table()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-filter_var_table)
- [`SingleCellDuckDBBase$get_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-get_cells_to_keep)
- [`SingleCellDuckDBBase$get_obs_cols()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-get_obs_cols)
- [`SingleCellDuckDBBase$get_obs_index_map()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-get_obs_index_map)
- [`SingleCellDuckDBBase$get_obs_table()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-get_obs_table)
- [`SingleCellDuckDBBase$get_var_index_map()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-get_var_index_map)
- [`SingleCellDuckDBBase$get_vars_adt_table()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-get_vars_adt_table)
- [`SingleCellDuckDBBase$get_vars_table()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-get_vars_table)
- [`SingleCellDuckDBBase$initialize()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-initialize)
- [`SingleCellDuckDBBase$join_data_obs()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-join_data_obs)
- [`SingleCellDuckDBBase$rename_column()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-rename_column)
- [`SingleCellDuckDBBase$set_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-set_cells_to_keep)
- [`SingleCellDuckDBBase$set_to_keep_column()`](https://gregorlueg.github.io/bixverse/reference/SingleCellDuckDBBase.html#method-set_to_keep_column)

------------------------------------------------------------------------

### `SingleCellDuckDB$populate_obs_from_h5ad()`

This function populates the obs table from an h5 file (if found).

#### Usage

    SingleCellDuckDB$populate_obs_from_h5ad(
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

### `SingleCellDuckDB$populate_vars_from_h5ad()`

This populates the vars table from an h5 file (if found.)

#### Usage

    SingleCellDuckDB$populate_vars_from_h5ad(h5_path, filter = NULL)

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

### `SingleCellDuckDB$populate_obs_from_multi_h5ad()`

Populate the obs table from multiple h5ad files.

#### Usage

    SingleCellDuckDB$populate_obs_from_multi_h5ad(
      per_file_info,
      cell_id_col = NULL
    )

#### Arguments

- `per_file_info`:

  List of lists, each containing: h5_path, exp_id, cell_filter
  (1-indexed integer vector of cells to keep).

- `cell_id_col`:

  Optional string. Column name for cell identifiers.

#### Returns

Invisible self. Populates the obs table in DuckDB.

------------------------------------------------------------------------

### `SingleCellDuckDB$populate_var_minimal()`

Populate a minimal var table for multi-file ingestion: gene_idx and
gene_id only. Per-file var annotations are intentionally not merged;
downstream annotation is the caller's responsibility (e.g. via an
external reference and a future enrich_var()).

#### Usage

    SingleCellDuckDB$populate_var_minimal(final_gene_names)

#### Arguments

- `final_gene_names`:

  Character vector. Gene ids in the final order (matches the Rust binary
  gene axis).

#### Returns

Invisible self. Function to populate the obs table from plain text files

------------------------------------------------------------------------

### `SingleCellDuckDB$populate_obs_from_plain_text()`

#### Usage

    SingleCellDuckDB$populate_obs_from_plain_text(f_path, has_hdr, filter = NULL)

#### Arguments

- `f_path`:

  File path to the plain text file.

- `has_hdr`:

  Boolean. Does the flat file have a header.

- `filter`:

  Optional integer. Positions of obs to read in from file.

#### Returns

Invisible self and populates the internal obs table. Function to
populate the var table from plain text files

------------------------------------------------------------------------

### `SingleCellDuckDB$populate_var_from_plain_text()`

#### Usage

    SingleCellDuckDB$populate_var_from_plain_text(f_path, has_hdr, filter = NULL)

#### Arguments

- `f_path`:

  File path to the plain text file.

- `has_hdr`:

  Boolean. Does the flat file have a header.

- `filter`:

  Optional integer. Positions of obs to read in from file.

#### Returns

Invisible self and populates the internal var table. Function to
populate the obs table from R

------------------------------------------------------------------------

### `SingleCellDuckDB$populate_obs_from_data.table()`

#### Usage

    SingleCellDuckDB$populate_obs_from_data.table(obs_dt, filter = NULL)

#### Arguments

- `obs_dt`:

  data.table

- `filter`:

  Optional integer. Row indices to keep

#### Returns

Invisible self and populates the internal obs table. Function to
populate the var table from R

------------------------------------------------------------------------

### `SingleCellDuckDB$populate_var_from_data.table()`

#### Usage

    SingleCellDuckDB$populate_var_from_data.table(var_dt, filter = NULL)

#### Arguments

- `var_dt`:

  data.table

- `filter`:

  Optional integer. Row indices to keep

#### Returns

Invisible self and populates the internal obs table. Function to
populate the var_adt table from R

------------------------------------------------------------------------

### `SingleCellDuckDB$populate_var_adt_from_data.table()`

#### Usage

    SingleCellDuckDB$populate_var_adt_from_data.table(var_dt, filter = NULL)

#### Arguments

- `var_dt`:

  data.table with `feature_idx` and `feature_id` columns.

- `filter`:

  Optional integer. Row indices to keep.

#### Returns

Invisible self and populates the internal var_adt table.

------------------------------------------------------------------------

### `SingleCellDuckDB$populate_obs_from_tenx_h5()`

Populate the obs table from the barcodes in a 10x CellRanger h5 file.

#### Usage

    SingleCellDuckDB$populate_obs_from_tenx_h5(h5_path, version, filter = NULL)

#### Arguments

- `h5_path`:

  String. Path to the 10x h5 file.

- `version`:

  String. One of `"v2"` or `"v3"`.

- `filter`:

  Optional integer. 1-indexed positions of cells to keep.

#### Returns

Invisible self. Populates the obs table in DuckDB.

------------------------------------------------------------------------

### `SingleCellDuckDB$populate_vars_from_tenx_h5()`

Populate the var table from the features in a 10x CellRanger h5 file.

Features are read in full file order (all modalities); `filter` selects
the gene-expression features that survived QC.

#### Usage

    SingleCellDuckDB$populate_vars_from_tenx_h5(h5_path, version, filter = NULL)

#### Arguments

- `h5_path`:

  String. Path to the 10x h5 file.

- `version`:

  String. One of `"v2"` or `"v3"`.

- `filter`:

  Optional integer. 1-indexed positions of features to keep.

#### Returns

Invisible self. Populates the var table in DuckDB.

------------------------------------------------------------------------

### `SingleCellDuckDB$populate_obs_from_multi_tenx_h5()`

Populate the obs table from multiple 10x CellRanger h5 files.

Reads each input's barcodes, applies the cell filter, prefixes `cell_id`
with `exp_id`, and rbindlists.

#### Usage

    SingleCellDuckDB$populate_obs_from_multi_tenx_h5(per_file_info)

#### Arguments

- `per_file_info`:

  List of lists; each must contain `h5_path` (string), `version` (`"v2"`
  or `"v3"`), `exp_id` (string) and `cell_filter` (1-indexed integer
  vector).

#### Returns

Invisible self.

------------------------------------------------------------------------

### `SingleCellDuckDB$populate_obs_from_multi_duckdb()`

Populate the obs table by merging obs from multiple source DuckDBs.

Reads obs rows where `to_keep = TRUE` from each source, ordered by
`cell_idx` ascending (matching the order used in the Rust bin merge).
Prefixes `cell_id` with `exp_id`, intersects column names across all
inputs, and rbindlists.

#### Usage

    SingleCellDuckDB$populate_obs_from_multi_duckdb(per_file_info)

#### Arguments

- `per_file_info`:

  List of lists; each must contain `db_path` (string) and `exp_id`
  (string).

#### Returns

Invisible self.

------------------------------------------------------------------------

### `SingleCellDuckDB$populate_vars_from_duckdb_reordered()`

Populate the var table from a source DuckDB, filtered and reordered to
match a target gene set.

#### Usage

    SingleCellDuckDB$populate_vars_from_duckdb_reordered(
      source_db_path,
      final_gene_names
    )

#### Arguments

- `source_db_path`:

  String. Path to the source DuckDB.

- `final_gene_names`:

  Character vector. Gene names in the desired final order.

#### Returns

Invisible self.

------------------------------------------------------------------------

### `SingleCellDuckDB$populate_obs_from_multi_plain_text()`

Populate obs from multiple plain-text barcode files.

Reads each input's barcodes file, filters to the cells that passed QC,
prefixes `cell_id` with `exp_id`, intersects column names across inputs,
and rbindlists.

#### Usage

    SingleCellDuckDB$populate_obs_from_multi_plain_text(per_file_info)

#### Arguments

- `per_file_info`:

  List of lists; each must contain `f_path` (string), `exp_id` (string),
  `has_hdr` (boolean), `cell_filter` (1-indexed integer vector).

#### Returns

Invisible self.

------------------------------------------------------------------------

### `SingleCellDuckDB$clone()`

The objects of this class are cloneable with this method.

#### Usage

    SingleCellDuckDB$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
